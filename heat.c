#include "heat.h"
#include "settings.h"
#include <stdio.h>
#include <math.h>

/* ===== РЕЖИМ 1: Сплайн из .dat ===== */
#if USE_SPLINE

static int nPoints = 0;
static double Tdat[MAX_POINTS];
static double lamdat[MAX_POINTS];
static double Cpdat[MAX_POINTS];
static double rhodat[MAX_POINTS];

int ReadProperties(const char* filename) {
    FILE* f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Не удалось открыть файл %s\n", filename);
        return -1;
    }

    char line[256];
    nPoints = 0;

    while (fgets(line, sizeof(line), f) && nPoints < MAX_POINTS) {
        if (line[0] == '#' || line[0] == '\n') {
            continue;
        }
        sscanf(line, "%lf %lf %lf %lf",
               &Tdat[nPoints], &rhodat[nPoints],
               &Cpdat[nPoints], &lamdat[nPoints]);
        Cpdat[nPoints] *= 1000.0; /* kJ -> J */
        nPoints++;
    }

    fclose(f);
    return nPoints;
}

double Interpolate(double T, const double Tdat[], const double val[], int n) {
    if (T <= Tdat[0])     return val[0];
    if (T >= Tdat[n - 1]) return val[n - 1];

    for (int k = 0; k < n - 1; k++) {
        if (T >= Tdat[k] && T <= Tdat[k + 1]) {
            return val[k] + (val[k + 1] - val[k]) / (Tdat[k + 1] - Tdat[k]) * (T - Tdat[k]);
        }
    }

    return val[n - 1];
}

double clambda(double T) {
    return Interpolate(T, Tdat, lamdat, nPoints);
}

double cCv(double T) {
    return Interpolate(T, Tdat, Cpdat, nPoints);
}

double crho(double T) {
    return Interpolate(T, Tdat, rhodat, nPoints);
}

double cqv(double T) {
    (void)T;
    return 0.0;
}

double cq(double t) {
    (void)t;
    return Q_LASER;
}

/* ===== РЕЖИМ 0: Полиномы SiO₂ (плавленый кварц, 300–1000 K) ===== */
#else

int ReadProperties(const char* filename) {
    (void)filename;
    return 0;
}

/* λ(T) = 1.38 − 6e-4·(T−300) + 2e-7·(T−300)²  [Вт/(м·К)] */
double clambda(double T) {
    double dT = T - 300.0;
    return 1.38 - 6e-4 * dT + 2e-7 * dT * dT;
}

/* Cv(T) = 745 + 0.5·(T−300) − 2e-4·(T−300)²  [Дж/(кг·К)] */
double cCv(double T) {
    double dT = T - 300.0;
    return 745.0 + 0.5 * dT - 2e-4 * dT * dT;
}

/* ρ(T) ≈ 2200 кг/м³ (слабая зависимость) */
double crho(double T) {
    (void)T;
    return 2200.0;
}

double cqv(double T) {
    (void)T;
    return 0.0;
}

double cq(double t) {
    return LASER_A * (exp(t) - 1.0);
}

/* β(T) ≈ 0.5e-6 1/K (коэф. линейного теплового расширения SiO₂) */
double cbeta(double T) {
    (void)T;
    return 0.5e-6;
}

/* E(T) ≈ 73e9 − 1e7·(T−300)  [Па] (модуль Юнга fused silica) */
double cE(double T) {
    return 73e9 - 1e7 * (T - 300.0);
}

double computeTension(double T_avg) {
    double dT = T_avg - T_INIT;
    return TENSION_0 - cE(T_avg) * H_THICK * cbeta(T_avg) * dT;
}

#endif

/* ===== Функции решателя ===== */

void makeMesh(struct cell_t *cells, struct face_t *faces, double length) {
    const double h = length / NC;

    for (int i = 1; i <= NC; ++i) {
        const double x = (i - 1) * h + h / 2.0;
        cells[i].x = x;
    }

    for (int i = 0; i < NF; ++i) {
        const int cl = i;
        const int cr = i + 1;
        faces[i].cl = cl;
        faces[i].cr = cr;
    }
}

void init(struct cell_t cells[NT], struct face_t faces[NF], double Tinit) {
    for (int i = 1; i <= NC; ++i) {
        const double r  = crho(Tinit);
        const double l  = clambda(Tinit);
        const double cv = cCv(Tinit);
        const double E  = r * cv * Tinit;

        cells[i].rho    = r;
        cells[i].lambda = l;
        cells[i].Cv     = cv;
        cells[i].E      = E;
        cells[i].T      = Tinit;
    }
    for (int f = 0; f < NF; ++f)
        faces[f].q = 0.0;
}

void applyBoundaries(struct cell_t cells[NT], double Tleft, double Tright, double t) {
#ifdef USE_LASER_FLUX
    const double q = cq(t);
    Tleft  = cells[1].T  + H / cells[1].lambda  * q;
#else
    (void)t;
#endif

    cells[0].T      = Tleft;
    cells[0].rho    = crho(Tleft);
    cells[0].lambda = clambda(Tleft);
    cells[0].Cv     = cCv(Tleft);
    cells[0].E      = cells[0].rho * cells[0].Cv * Tleft;

#ifdef USE_LASER_FLUX
    Tright = cells[NC].T - H / cells[NC].lambda * q;
#endif

    cells[NC + 1].T      = Tright;
    cells[NC + 1].rho    = crho(Tright);
    cells[NC + 1].lambda = clambda(Tright);
    cells[NC + 1].Cv     = cCv(Tright);
    cells[NC + 1].E      = cells[NC + 1].rho * cells[NC + 1].Cv * Tright;
}

void computeFluxes(struct cell_t cells[NT], struct face_t faces[NF], double tau) {
    for (int c = 0; c < NT; ++c)
        cells[c].q = 0.0;

    for (int f = 0; f < NF; ++f) {
        const int    cl  = faces[f].cl;
        const int    cr  = faces[f].cr;
        const double Tl  = cells[cl].T;
        const double Tr  = cells[cr].T;
        const double Tf  = (Tl + Tr) / 2.0;

        const double lf   = clambda(Tf);
        const double rf   = crho(Tf);
        const double cvf  = cCv(Tf);
        const double tau_c = 3.0 * lf / (rf * cvf * C_PHONON * C_PHONON);

        /* Каттанео (неявная): q^{n+1}·(1 + τ/τ_c) = q^n − (τ/τ_c)·λ·(Tr−Tl)/h
         * При τ >> τ_c (τ_c ≈ 8e-14 с для SiO₂) → вырождается в Фурье.
         * Явная схема неустойчива при τ/τ_c ≈ 7e9 (NaN на первом шаге). */
        const double ratio = tau / tau_c;
        faces[f].q = (faces[f].q - ratio * lf * (Tr - Tl) / H) / (1.0 + ratio);

        cells[cl].q += faces[f].q;
        cells[cr].q -= faces[f].q;
    }
}

void integrate(struct cell_t cells[NT], double tau) {
    for (int c = 1; c <= NC; ++c) {
        cells[c].E -= (tau / H) * cells[c].q;
    }
}

void computeSources(struct cell_t cells[NT], double tau) {
    for (int c = 1; c <= NC; ++c) {
        cells[c].qv = cqv(cells[c].T);
        cells[c].E += tau * cells[c].qv;
    }
}

void updateTemp(struct cell_t cells[NT]) {
    for (int c = 1; c <= NC; ++c) {
        cells[c].T = cells[c].E / (cells[c].rho * cells[c].Cv);
    }
}

void updateParams(struct cell_t cells[NT]) {
    for (int c = 1; c <= NC; ++c) {
        cells[c].rho    = crho(cells[c].T);
        cells[c].lambda = clambda(cells[c].T);
        cells[c].Cv     = cCv(cells[c].T);
    }
}

/* ===== Тепловой момент и изгибная жёсткость ===== */

double computeThermalMoment(struct cell_t cells[NT]) {
    double M = 0.0;
    for (int i = 1; i <= NC; ++i) {
        const double T  = cells[i].T;
        const double zi = cells[i].x - H_THICK / 2.0;
        M += zi * (-cE(T) * cbeta(T) * (T - T_INIT)) * H;
    }
    return M;
}

double computeBendingStiffness(struct cell_t cells[NT]) {
    double D = 0.0;
    for (int i = 1; i <= NC; ++i) {
        const double T  = cells[i].T;
        const double zi = cells[i].x - H_THICK / 2.0;
        D += zi * zi * cE(T) * H;
    }
    return D;
}

double averageTemperature(struct cell_t cells[NT]) {
    double sum = 0.0;
    for (int i = 1; i <= NC; ++i) sum += cells[i].T;
    return sum / NC;
}

/* ===== Аналитическое решение (ряд Фурье) ===== */

#define K_TERMS 100

void analyticalSolution(double *x, double *Ta, int n, double t) {
    const double lam = clambda(T_INIT);
    const double r   = crho(T_INIT);
    const double cv  = cCv(T_INIT);
    const double alpha = lam / (r * cv);

    for (int i = 0; i < n; ++i) {
        /* Стационарная часть */
        double T = (T2 - T1) / H_THICK * x[i] + T1;

        /* Нестационарная часть (ряд Фурье) */
        for (int k = 1; k <= K_TERMS; ++k) {
            const double lambda_k = (PI * k / H_THICK) * (PI * k / H_THICK);
            const double sign_k = (k % 2 == 0) ? 1.0 : -1.0; /* (-1)^k */
            const double Ak = (2.0 / (PI * k)) *
                ((T_INIT - T1) * (1.0 - sign_k) + (T2 - T1) * sign_k);
            T += Ak * sin(PI * k * x[i] / H_THICK) * exp(-alpha * lambda_k * t);
        }

        Ta[i] = T;
    }
}
