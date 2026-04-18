#ifndef STATIC_H
#define STATIC_H
#pragma once
double lambda(int i);
double phi(int i, int n);
void fourier(const double* input, double* output);
void calculateStaticDeflection(double* u, double tau);
double W(double x);
double Px(double i);
void performStep(double* wi, double* wim1, double tau, double eps, double thermal_load);
#endif // STATIC_H
