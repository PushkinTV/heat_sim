#ifndef SETTINGS_H
#define SETTINGS_H
#pragma once

/* Режим: 1 — сплайн из .dat, 0 — готовая формула */
#define USE_SPLINE 1

/* Файл свойств (для USE_SPLINE = 1) */
#define DATA_FILE "copper_properties.dat"

/* Сеточные параметры */
#define NC 100          /* число ячеек */
#define NN (NC + 1)     /* число узлов */
#define NG 2            /* призрачные ячейки */
#define NT (NC + NG)    /* общее число ячеек */
#define NF NN           /* число граней */

/* Физические параметры */
#define L   1.0         /* длина области [м] */
#define H   (L / NC)    /* шаг сетки */

#define PI 3.1415926535897932384626433

/* Тип ГУ: раскомментировать для ГУ 2 рода (поток) */
#define USE_LASER_FLUX

/* Начальная температура и ГУ 1 рода */
#define T_INIT 300.0    /* начальная T [K] */
#define T1     400.0    /* T на левой границе [K] */
#define T2     300.0    /* T на правой границе [K] */

/* Параметры лазерного потока (для USE_LASER_FLUX) */
#define Q_LASER 1.0e5   /* поток [Вт/м²] */

/* Временные параметры */
#define TAU            0.1     /* шаг по времени [с] */
#define INTEGRATE_TIME 100.0   /* время интегрирования [с] */

#endif // SETTINGS_H
