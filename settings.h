#ifndef SETTINGS_H
#define SETTINGS_H
#pragma once

/* Режим: 1 — сплайн из .dat, 0 — полиномы SiO₂ */
#define USE_SPLINE 0

/* Файл свойств (для USE_SPLINE = 1) */
#define DATA_FILE "copper_properties.dat"

/* Сеточные параметры (тепловая задача по z) */
#define NC 20           /* число ячеек по z */
#define NN (NC + 1)     /* число узлов */
#define NG 2            /* призрачные ячейки */
#define NT (NC + NG)    /* общее число ячеек */
#define NF NN           /* число граней */

/* Геометрия мембраны */
#define H_THICK 1.0e-3          /* толщина по z [м] */
#define H       (H_THICK / NC)  /* шаг сетки по z [м] */

#define PI 3.1415926535897932384626433

/* ГУ 2 рода: лазерный поток */
#define USE_LASER_FLUX

/* Начальная температура */
#define T_INIT  573.0   /* 300°C по лекции 04.01 */
/* T1, T2 — не используются при USE_LASER_FLUX, оставлены для аналитики */
#define T1      400.0
#define T2      300.0

/* Лазерный поток: q₀(t) = LASER_A * (exp(t) - 1) */
#define LASER_A   1.1e-3  /* подобрано под T_верх ≈ 973 K при t=20 */
#define Q_LASER   LASER_A  /* алиас для applyBoundaries (cq) */

/* Каттанео: скорость фононов в SiO₂; τ_c = 3λ/(ρ·Cv·C̄²) считается на лету */
#define C_PHONON 5500.0   /* [м/с] */

/* Адаптивный шаг: CFL */
#define CFL_K          0.5
#define TAU            1e-10   /* начальный шаг [с], будет пересчитан */
#define INTEGRATE_TIME 20.0    /* время интегрирования [с] */

/* ===== Параметры балки (из NoStatic, переименованы) ===== */
#define N_BEAM    101           /* узлов по x (балка) */
#define M_BEAM    50            /* наблюдаемая точка (N_BEAM/2) */
#define L_BEAM    1e-3          /* длина мембраны по x [м] */
#define H_BEAM    (1.0 / (N_BEAM - 1))     /* нормированный шаг балки на [0,1] — как в NoStatic */
#define TAU_BEAM_MID 1e-2                  /* псевдо-шаг для calculateStaticDeflection */
#define TENSION_0 10.0          /* начальное натяжение [Н/м] */
#define R_TOL     1e-6          /* допуск сходимости балки */
#define ALPHA_LOAD 0.7          /* параметр нагрузки (из NoStatic) */

/* Расчёт газа */
#define D_GAP     1e-4          /* характерный размер канала [м] */

#endif // SETTINGS_H
