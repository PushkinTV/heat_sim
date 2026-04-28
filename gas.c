#include "gas.h"
#include "settings.h"
#include <math.h>
#include <stdio.h>

#define T_ENV    300.0      /* K, комнатная температура */
#define NU_AIR   1.57e-5    /* м²/с */
#define LAM_AIR  0.026      /* Вт/(м·К) */
#define PR_AIR   0.71

#define NU_LAM   4.36       /* Nu ламинар, постоянный поток */
#define RE_LAM   2300.0
#define RE_TURB  4000.0

void analyzeConvection(double T_wall, double q_wall) {
    printf("\n=== Конвективный теплообмен (газ, п.8) ===\n");
    printf("T_wall = %.1f K   T_env = %.1f K   q_wall = %.3e Вт/м²\n",
           T_wall, T_ENV, q_wall);

    /* Шаг 1-2: из известного потока найти Nu */
    double h_conv = q_wall / (T_wall - T_ENV);
    double Nu     = h_conv * D_GAP / LAM_AIR;

    /* Шаг 3-4: итерационная проверка согласованности режима.
     * Начинаем с предположения о турбулентности.
     * Если Re не попадает в зону предположенного режима — переключаем. */
    int turb = (Nu > NU_LAM);   /* 1 = турбулент, 0 = ламинар */
    double Re = 0.0;
    const char *reg = "";

    for (int iter = 0; iter < 4; ++iter) {
        if (turb) {
            Re  = pow(Nu / (0.023 * pow(PR_AIR, 0.4)), 1.0 / 0.8);
            reg = "turb";
            if (Re < RE_LAM) {
                printf("  итер %d: предположен турбулент → Re=%.0f < 2300, переключаю на ламинар\n",
                       iter + 1, Re);
                turb = 0;
                continue;
            }
        } else {
            /* Ламинарный: Nu фиксирован = NU_LAM, берём Re_max */
            Re  = RE_LAM;
            reg = "lam";
            if (Nu > NU_LAM) {
                printf("  итер %d: предположен ламинар → Nu=%.2f > %.2f, переключаю на турбулент\n",
                       iter + 1, Nu, NU_LAM);
                turb = 1;
                continue;
            }
        }
        break;  /* режим согласован */
    }

    if (Re >= RE_LAM && Re <= RE_TURB)
        printf("  ПРЕДУПРЕЖДЕНИЕ: переходный режим (2300 ≤ Re ≤ 4000)\n");

    double u_e = Re * NU_AIR / D_GAP;

    printf("Nu = %.2f   Re = %.1f   режим: %s   h = %.2f Вт/(м²·К)   u_e_min = %.3f м/с\n",
           Nu, Re, reg, h_conv, u_e);
}
