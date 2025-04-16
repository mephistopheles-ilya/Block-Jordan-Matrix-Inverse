#pragma once


void ij2l(int nx, int /*ny*/, int i, int j, int& l);
void l2ij(int nx, int /*ny*/, int& i, int& j, int l);
int get_len_msr(int nx, int ny);
int get_off_diag(int nx, int ny, int i, int j, int* I_ij = nullptr);
int get_len_msr_off_diag(int nx, int ny);
int alocate_msr_matrix(int nx, int ny, double** p_A, int** p_I);
void fill_I(int nx, int ny, int* I);
void fill_A_ij(int nx, int ny, double hx, double hy, int i,int j, double* A_diag, double* A_offdiag);
void fill_A(int nx, int ny, double hx, double hy, int* I, double* A, int p, int k);
double F_IJ(int nx, int ny, double hx, double hy, double x0, double y0, int i, int j, double (*f)(double, double), double add_error);
void fill_B(double a, double c, int nx, int ny, double hx, double hy, double* b, int p, int k, double (*f)(double, double), double add_error);
