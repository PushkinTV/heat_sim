#include "gas.h"
#include "settings.h"
#include <math.h>
#include <stdio.h>

#define T_ENV    300.0      /* K, комнатная температура */
#define NU_AIR   1.57e-5    /* м²/с */
#define LAM_AIR  0.026      /* Вт/(м·К) */
#define PR_AIR   0.71

/* Nu для канала с постоянным тепловым потоком (ламинар) */
#define NU_LAM   4.36

void analyzeConvection(double T_wall, double q_wall) {
    printf("\n=== Конвективный теплообмен (газ, п.8) ===\n");
    printf("T_wall = %.1f K   T_env = %.1f K   q_wall = %.3e Вт/м²\n",
           T_wall, T_ENV, q_wall);

    /* Шаг 1: коэффициент теплоотдачи из известного потока и перепада */
    double h_conv = q_wall / (T_wall - T_ENV);

    /* Шаг 2: число Нуссельта */
    double Nu = h_conv * D_GAP / LAM_AIR;

    double Re;
    const char *reg;

    /* Шаг 3: определить режим и найти Re */
    if (Nu <= NU_LAM) {
        /* Ламинарный режим покрывает нужный Nu — Re ≤ 2300, берём максимум */
        Re  = 2300.0;
        reg = "lam";
    } else {
        /* Турбулентный: Nu = 0.023 · Re^0.8 · Pr^0.4  →  Re = ... */
        Re  = pow(Nu / (0.023 * pow(PR_AIR, 0.4)), 1.0 / 0.8);
        reg = (Re > 4000.0) ? "turb" : "trans";
    }

    double u_e = Re * NU_AIR / D_GAP;

    /* Шаг 4: результат */
    printf("Nu = %.2f   Re = %.1f   режим: %s   h = %.2f Вт/(м²·К)   u_e_min = %.3f м/с\n",
           Nu, Re, reg, h_conv, u_e);

    /* Проверки согласованности */
    if (Nu > NU_LAM && Re < 2300.0)
        printf("ПРЕДУПРЕЖДЕНИЕ: Nu > Nu_lam, но Re < 2300 — режим несогласован\n");
    else if (Re >= 2300.0 && Re <= 4000.0)
        printf("ПРЕДУПРЕЖДЕНИЕ: переходный режим (2300 < Re ≤ 4000)\n");
}
