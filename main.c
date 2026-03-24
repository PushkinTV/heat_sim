#include "heat.h"
#include "settings.h"
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

    makeMesh(cells, faces, L);
    init(cells, T_INIT);
    applyBoundaries(cells, T1, T2, 0.0);

    /* Таймеры (наносекундное разрешение) */
    struct timespec ts_start, ts_end;
    double ms_fluxes = 0, ms_integrate = 0, ms_sources = 0;
    double ms_temp = 0, ms_params = 0, ms_bc = 0;

    struct timespec ts_total_start;
    clock_gettime(CLOCK_MONOTONIC, &ts_total_start);

    double time = 0.0;
    for (int counter = 0; time <= INTEGRATE_TIME; ++counter) {
        clock_gettime(CLOCK_MONOTONIC, &ts_start);
        computeFluxes(cells, faces);
        clock_gettime(CLOCK_MONOTONIC, &ts_end);
        ms_fluxes += elapsed_ms(ts_start, ts_end);

        clock_gettime(CLOCK_MONOTONIC, &ts_start);
        integrate(cells, TAU);
        clock_gettime(CLOCK_MONOTONIC, &ts_end);
        ms_integrate += elapsed_ms(ts_start, ts_end);

        clock_gettime(CLOCK_MONOTONIC, &ts_start);
        computeSources(cells, TAU);
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

        time += TAU;

        if (counter % 10 == 0) {
            printf("%f ...\n", time);
        }
    }

    struct timespec ts_total_end;
    clock_gettime(CLOCK_MONOTONIC, &ts_total_end);
    double ms_total = elapsed_ms(ts_total_start, ts_total_end);

    /* Вывод таймеров */
    printf("\n=== Profiling ===\n");
    printf("Total:           %.3f ms\n", ms_total);
    printf("computeFluxes:   %.3f ms\n", ms_fluxes);
    printf("integrate:       %.3f ms\n", ms_integrate);
    printf("computeSources:  %.3f ms\n", ms_sources);
    printf("updateTemp:      %.3f ms\n", ms_temp);
    printf("updateParams:    %.3f ms\n", ms_params);
    printf("applyBoundaries: %.3f ms\n", ms_bc);

    /* Вывод: T(x, t_end) — численный профиль */
    FILE *fp = fopen("T_profile.dat", "w");
    fprintf(fp, "# x  T_numerical\n");
    for (int i = 1; i <= NC; ++i) {
        fprintf(fp, "%f %f\n", cells[i].x, cells[i].T);
    }
    fclose(fp);

    /* Аналитическое решение в последний момент времени */
    double x_arr[NC], Ta[NC];
    for (int i = 0; i < NC; ++i) {
        x_arr[i] = cells[i + 1].x;
    }
    analyticalSolution(x_arr, Ta, NC, time);

    FILE *fa = fopen("T_analytical.dat", "w");
    fprintf(fa, "# x  T_analytical\n");
    for (int i = 0; i < NC; ++i) {
        fprintf(fa, "%f %f\n", x_arr[i], Ta[i]);
    }
    fclose(fa);

    return 0;
}
