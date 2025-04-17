#pragma once


double Pf(double* res, double x, double y, double a, double c, double hx, double hy, int nx, int ny);
double calc_r1(double* res, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double));
double calc_r2(double* res, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double));
double calc_r3(double* res, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double)
        , double add_error);
double calc_r4(double* res, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double)
        , double add_error);
double find_max(double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double));
