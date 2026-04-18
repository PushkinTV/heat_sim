#include "static.h"
#include "heat.h"
#include "settings.h"
#include <math.h>

double lambda(int i) {
    return -4.0 / (H_BEAM * H_BEAM)
           * sin(PI * i / (2.0 * (N_BEAM - 1.)))
           * sin(PI * i / (2.0 * (N_BEAM - 1.)));
}

double phi(int i, int n) {
    const double norm = sqrt(2.0 / (N_BEAM - 1.));
    const double f    = sin(PI * n * i / (N_BEAM - 1.));
    return norm * f;
}

void fourier(const double* input, double* output) {
    for (int i = 0; i < N_BEAM; ++i) {
        double sum = 0.0;
        for (int n = 0; n < N_BEAM; ++n) {
            sum += input[n] * phi(i, n);
        }
        output[i] = sum;
    }
}

void calculateStaticDeflection(double* u, double tau) {
    double wi[N_BEAM];
    double wip1[N_BEAM];
    double cs[N_BEAM];
    double csp1[N_BEAM];
    double pi[N_BEAM];
    double pe[N_BEAM];
    /* eps = E(T_INIT)·h³/(12·L²·T) — изгибная жёсткость (безразмерная) */
    double eps = (cE(T_INIT) * H_THICK * H_THICK * H_THICK)
                 / (12.0 * L_BEAM * L_BEAM * TENSION_0);

    for (int i = 0; i < N_BEAM; ++i)
        wi[i] = 0.0;

    double r_tol = 1.0;
    for (int s = 0; s < 1000 && r_tol > R_TOL; ++s) {
        fourier(wi, cs);

        for (int i = 0; i < N_BEAM; ++i)
            pi[i] = Px(i);

        fourier(pi, pe);

        for (int i = 0; i < N_BEAM; ++i) {
            csp1[i] = (cs[i] + ALPHA_LOAD * pe[i] * tau)
                      / (1.0 - lambda(i) * tau
                         + eps * tau * lambda(i) * lambda(i));
        }

        fourier(csp1, wip1);

        wip1[0]          = 0.0;
        wip1[N_BEAM - 1] = 0.0;

        r_tol = fabs(wip1[M_BEAM] - wi[M_BEAM]);

        for (int i = 0; i < N_BEAM; ++i)
            wi[i] = wip1[i];
    }

    for (int i = 0; i < N_BEAM; ++i)
        u[i] = wi[i];
}

double W(double x) {
    return ALPHA_LOAD / 4.0 * (x * (1.0 - x)
           + 1.0 / PI / PI * sin(PI * x) * sin(PI * x));
}

double Px(double i) {
    double x = i * H_BEAM;
    return sin(PI * x) * sin(PI * x);
}

/* Один шаг схемы Эйлера-Бернулли.
 * eps = D_th/(L²·T) — изгибная жёсткость, пересчитывается снаружи.
 * thermal_load = T·M/D — равномерная нагрузка от теплового момента. */
void performStep(double* wi, double* wim1, double tau, double eps, double thermal_load) {
    double wip1[N_BEAM];
    double cs[N_BEAM],  csp1[N_BEAM],  csm1[N_BEAM];
    double p_arr[N_BEAM], pe[N_BEAM];

    fourier(wim1, csm1);
    fourier(wi,   cs);

    for (int i = 0; i < N_BEAM; ++i)
        p_arr[i] = ALPHA_LOAD * Px(i) + thermal_load;
    fourier(p_arr, pe);

    for (int i = 0; i < N_BEAM; ++i) {
        const double lam = lambda(i);
        csp1[i] = (2.0*cs[i] - csm1[i] + tau*tau*pe[i])
                  / (1.0 - lam*tau*tau + eps*tau*tau*lam*lam);
    }

    fourier(csp1, wip1);
    wip1[0]          = 0.0;
    wip1[N_BEAM - 1] = 0.0;

    for (int i = 0; i < N_BEAM; ++i) {
        wim1[i] = wi[i];
        wi[i]   = wip1[i];
    }
}
