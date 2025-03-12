#pragma once

void Newton_approximation(double* points_x, double* f_x, double* newton_coeff, int n);

double Newton_calck(double* newton_coeff, double* points_x, int n, double x);

void cube_three_points(double x, double y, double* b, double* res);

void cube_approximation(double* points_x, double* f_x, double* d, double* c, int n);

double cube_calc(double* points_x, double* c, double a, double b, int n, double x);

