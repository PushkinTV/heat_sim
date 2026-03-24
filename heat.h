#ifndef HEAT_H
#define HEAT_H
#pragma once

#include "settings.h"

#define MAX_POINTS 50

/* --- Структуры данных --- */

struct cell_t {
    double rho;
    double lambda;
    double qv;
    double E;
    double Cv;
    double T;
    double x;
    double q;
};

struct face_t {
    double T;
    double lambda;
    double q;
    int cl;  /* индекс ячейки слева */
    int cr;  /* индекс ячейки справа */
};

/* --- Свойства вещества --- */

/* Чтение .dat файла */
int ReadProperties(const char* filename);

/* Линейная интерполяция */
double Interpolate(double T, const double Tdat[], const double val[], int n);

/* Свойства вещества */
double clambda(double T);
double cCv(double T);
double crho(double T);
double cqv(double T);
double cq(double t);

/* --- Функции решателя --- */
void makeMesh(struct cell_t *cells, struct face_t *faces, double length);
void init(struct cell_t cells[NT], double Tinit);
void applyBoundaries(struct cell_t cells[NT], double Tleft, double Tright, double t);
void computeFluxes(struct cell_t cells[NT], struct face_t faces[NF]);
void integrate(struct cell_t cells[NT], double tau);
void computeSources(struct cell_t cells[NT], double tau);
void updateTemp(struct cell_t cells[NT]);
void updateParams(struct cell_t cells[NT]);

/* --- Аналитическое решение (ряд Фурье) --- */
void analyticalSolution(double *x, double *Ta, int n, double t);

#endif // HEAT_H
