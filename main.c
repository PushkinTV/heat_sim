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

    /* --- Механический прогиб (без температуры) — один раз --- */
    double w_mech[N_BEAM];
    calculateStaticDeflection(w_mech, TAU_BEAM_MID);
    {
        FILE *fm = fopen("deflection_mech.dat", "w");
        fprintf(fm, "# x_norm\tw_mech\n");
        for (int i = 0; i < N_BEAM; ++i)
            fprintf(fm, "%.6e\t%.6e\n", (double)i * H_BEAM, w_mech[i]);
        fclose(fm);
    }

    /* --- Термомеханический прогиб — старт из нуля --- */
    double wi_thermo[N_BEAM], wim1_thermo[N_BEAM];
    for (int i = 0; i < N_BEAM; ++i) {
        wi_thermo[i]   = 0.0;
        wim1_thermo[i] = 0.0;
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

        /* Термомеханическая связь */
        double M_th  = computeThermalMoment(cells);
        double D_th  = computeBendingStiffness(cells);
        double T_avg = averageTemperature(cells);
        double tens  = computeTension(T_avg);
        double eps_t = D_th / (L_BEAM * L_BEAM * tens);
        double thermal_load = TENSION_0 * M_th / D_th;

        performStep(wi_thermo, wim1_thermo, tau, eps_t, thermal_load);

        time += tau;

        if (counter % 1000 == 0) {
            printf("t=%.3e  T_top=%.1f  T_bot=%.1f  M=%.3e  D=%.3e  w_mid=%.3e\n",
                   time, cells[1].T, cells[NC].T, M_th, D_th, wi_thermo[M_BEAM]);
        }
    }

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
        fprintf(ft, "%.6e\t%.6e\n", (double)i * H_BEAM, wi_thermo[i]);
    fclose(ft);

    /* === Результаты === */
    double M_fin = computeThermalMoment(cells);
    double D_fin = computeBendingStiffness(cells);
    printf("\n=== Термоупругость ===\n");
    printf("M (тепловой момент)   = %.6e Н\n",    M_fin);
    printf("D (изгибная жёсткость)= %.6e Н·м²\n", D_fin);
    printf("w_thermo[M_BEAM]      = %.6e\n",       wi_thermo[M_BEAM]);
    printf("w_mech  [M_BEAM]      = %.6e\n",       w_mech[M_BEAM]);
    printf("T2 = %f K\n", cells[NC].T);

    /* п.8: конвективный теплообмен с газом */
    double q_laser_final = LASER_A * (exp(INTEGRATE_TIME) - 1.0);
    analyzeConvection(cells[NC].T, q_laser_final);

    return 0;
}
