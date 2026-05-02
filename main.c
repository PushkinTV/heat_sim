#include "heat.h"
#include "static.h"
#include "settings.h"
#include "gas.h"
#include <stdio.h>
#include <math.h>
#include <time.h>

static double elapsed_ms(struct timespec start, struct timespec end) {
    return (end.tv_sec - start.tv_sec) * 1000.0
         + (end.tv_nsec - start.tv_nsec) / 1.0e6;
}

int main(void) {
    int n = ReadProperties(DATA_FILE);
    if (n < 0) return 1;

    struct cell_t cells[NT];
    struct face_t faces[NF];

    makeMesh(cells, faces, H_THICK);
    init(cells, faces, T_INIT);
    applyBoundaries(cells, T1, T2, 0.0);

    /* --- Механический прогиб (без температуры) — статика, один раз ---
     * eps_mech = E·h³/(12·L²·T), extra_load = 0 → −w'' = α·sin²(πx). */
    double w_mech[N_BEAM];
    {
        const double eps_mech = (cE(T_INIT) * H_THICK * H_THICK * H_THICK)
                              / (12.0 * L_BEAM * L_BEAM * TENSION_0);
        calculateStaticDeflection(w_mech, TAU_BEAM_MID, eps_mech, 0.0);

        FILE *fm = fopen("deflection_mech.dat", "w");
        fprintf(fm, "# x_norm\tw_mech\n");
        for (int i = 0; i < N_BEAM; ++i)
            fprintf(fm, "%.6e\t%.6e\n", (double)i * H_BEAM, w_mech[i]);
        fclose(fm);
    }

    /* Таймеры */
    struct timespec ts_start, ts_end;
    double ms_fluxes = 0, ms_integrate = 0, ms_sources = 0;
    double ms_temp = 0, ms_params = 0, ms_bc = 0;

    struct timespec ts_total_start;
    clock_gettime(CLOCK_MONOTONIC, &ts_total_start);

    double time = 0.0;
    for (int counter = 0; time <= INTEGRATE_TIME; ++counter) {
        /* Диффузионный CFL: τ = k·h²/(2·α_max) */
        double alpha_max = 0.0;
        double tau_c_min = 1e30;
        for (int c = 1; c <= NC; ++c) {
            double a   = cells[c].lambda / (cells[c].rho * cells[c].Cv);
            double tac = 3.0 * cells[c].lambda /
                         (cells[c].rho * cells[c].Cv * C_PHONON * C_PHONON);
            if (a   > alpha_max) alpha_max = a;
            if (tac < tau_c_min) tau_c_min = tac;
        }
        double tau = CFL_K * H * H / (2.0 * alpha_max);
        (void)tau_c_min;

        /* Тепловой решатель */
        clock_gettime(CLOCK_MONOTONIC, &ts_start);
        computeFluxes(cells, faces, tau);
        clock_gettime(CLOCK_MONOTONIC, &ts_end);
        ms_fluxes += elapsed_ms(ts_start, ts_end);

        clock_gettime(CLOCK_MONOTONIC, &ts_start);
        integrate(cells, tau);
        clock_gettime(CLOCK_MONOTONIC, &ts_end);
        ms_integrate += elapsed_ms(ts_start, ts_end);

        clock_gettime(CLOCK_MONOTONIC, &ts_start);
        computeSources(cells, tau);
        clock_gettime(CLOCK_MONOTONIC, &ts_end);
        ms_sources += elapsed_ms(ts_start, ts_end);

        clock_gettime(CLOCK_MONOTONIC, &ts_start);
        updateTemp(cells);
        clock_gettime(CLOCK_MONOTONIC, &ts_end);
        ms_temp += elapsed_ms(ts_start, ts_end);

        clock_gettime(CLOCK_MONOTONIC, &ts_start);
        updateParams(cells);
        clock_gettime(CLOCK_MONOTONIC, &ts_end);
        ms_params += elapsed_ms(ts_start, ts_end);

        clock_gettime(CLOCK_MONOTONIC, &ts_start);
        applyBoundaries(cells, T1, T2, time);
        clock_gettime(CLOCK_MONOTONIC, &ts_end);
        ms_bc += elapsed_ms(ts_start, ts_end);

        time += tau;

        if (counter % 1000 == 0) {
            double M_th = computeThermalMoment(cells);
            double D_th = computeBendingStiffness(cells);
            printf("t=%.3e  T_top=%.1f  T_bot=%.1f  M=%.3e  D=%.3e\n",
                   time, cells[1].T, cells[NC].T, M_th, D_th);
        }
    }

    /* === Термоупругий прогиб — статика на финальных коэффициентах ===
     * Тепловое поле к концу пробега практически выходит на квазиравновесие;
     * берём итоговые M_th, D_th, T_tens и считаем равновесный прогиб
     * расширенным calculateStaticDeflection (с extra_load = T·M/D). */
    double M_th_fin = computeThermalMoment(cells);
    double D_th_fin = computeBendingStiffness(cells);
    double T_avg_fin = averageTemperature(cells);
    double tens_fin = computeTension(T_avg_fin);
    double eps_th = D_th_fin / (L_BEAM * L_BEAM * tens_fin);
    double thermal_load_fin = TENSION_0 * M_th_fin / D_th_fin;

    double w_thermo[N_BEAM];
    calculateStaticDeflection(w_thermo, TAU_BEAM_MID, eps_th, thermal_load_fin);

    struct timespec ts_total_end;
    clock_gettime(CLOCK_MONOTONIC, &ts_total_end);
    double ms_total = elapsed_ms(ts_total_start, ts_total_end);

    /* === Profiling === */
    printf("\n=== Profiling ===\n");
    printf("Total:           %.3f ms\n", ms_total);
    printf("computeFluxes:   %.3f ms\n", ms_fluxes);
    printf("integrate:       %.3f ms\n", ms_integrate);
    printf("computeSources:  %.3f ms\n", ms_sources);
    printf("updateTemp:      %.3f ms\n", ms_temp);
    printf("updateParams:    %.3f ms\n", ms_params);
    printf("applyBoundaries: %.3f ms\n", ms_bc);

    /* T_profile.dat */
    FILE *fp = fopen("T_profile.dat", "w");
    fprintf(fp, "# x\tT_numerical\n");
    for (int i = 1; i <= NC; ++i)
        fprintf(fp, "%.6e\t%.6e\n", cells[i].x, cells[i].T);
    fclose(fp);

    /* deflection_thermo.dat */
    FILE *ft = fopen("deflection_thermo.dat", "w");
    fprintf(ft, "# x_norm\tw_thermo\n");
    for (int i = 0; i < N_BEAM; ++i)
        fprintf(ft, "%.6e\t%.6e\n", (double)i * H_BEAM, w_thermo[i]);
    fclose(ft);

    /* === Результаты === */
    printf("\n=== Термоупругость (статическое равновесие) ===\n");
    printf("M (тепловой момент)   = %.6e Н\n",    M_th_fin);
    printf("D (изгибная жёсткость)= %.6e Н·м²\n", D_th_fin);
    printf("T_tens (натяжение)    = %.6e Н/м\n",  tens_fin);
    printf("eps_th                = %.6e\n",      eps_th);
    printf("thermal_load          = %.6e\n",      thermal_load_fin);
    printf("w_thermo[M_BEAM]      = %.6e\n",      w_thermo[M_BEAM]);
    printf("w_mech  [M_BEAM]      = %.6e\n",      w_mech[M_BEAM]);
    printf("ratio w_thermo/w_mech = %.3f\n",
           w_thermo[M_BEAM] / w_mech[M_BEAM]);
    printf("T_top=%.2f K  T_bot=%.2f K  ΔT=%.2f K\n",
           cells[1].T, cells[NC].T, cells[1].T - cells[NC].T);

    /* п.8: конвективный теплообмен с газом */
    double q_laser_final = LASER_A * (exp(INTEGRATE_TIME / LASER_T_REF) - 1.0);
    analyzeConvection(cells[NC].T, q_laser_final);

    return 0;
}
