#include "heat.h"
#include "settings.h"
#include <stdio.h>

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

/* ===== РЕЖИМ 0: Готовые формулы ===== */
#else

int ReadProperties(const char* filename) {
    (void)filename;
    return 0;
}

double clambda(double T) {
    (void)T;
    return 401.0;
}

double cCv(double T) {
    (void)T;
    return 386.0;
}

double crho(double T) {
    (void)T;
    return 8930.0;
}

double cqv(double T) {
    (void)T;
    return 0.0;
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

void init(struct cell_t cells[NT], double Tinit) {
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
}

void applyBoundaries(struct cell_t cells[NT], double Tleft, double Tright) {
    cells[0].T      = Tleft;
    cells[0].rho    = crho(Tleft);
    cells[0].lambda = clambda(Tleft);
    cells[0].Cv     = cCv(Tleft);
    cells[0].E      = cells[0].rho * cells[0].Cv * Tleft;

    cells[NC + 1].T      = Tright;
    cells[NC + 1].rho    = crho(Tright);
    cells[NC + 1].lambda = clambda(Tright);
    cells[NC + 1].Cv     = cCv(Tright);
    cells[NC + 1].E      = cells[NC + 1].rho * cells[NC + 1].Cv * Tright;
}

static double computeFlux(double lambda, double Tl, double Tr, double h) {
    return -lambda * (Tr - Tl) / h;
}

void computeFluxes(struct cell_t cells[NT], struct face_t faces[NF]) {
    for (int c = 0; c < NT; ++c) {
        cells[c].q = 0.0;
    }

    for (int f = 0; f < NF; ++f) {
        const double Tl = cells[faces[f].cl].T;
        const double Tr = cells[faces[f].cr].T;
        const double Tf = (Tl + Tr) / 2.0;

        const double lambdaf = clambda(Tf);
        const double q = computeFlux(lambdaf, Tl, Tr, H);

        cells[faces[f].cl].q += q;
        cells[faces[f].cr].q -= q;
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
