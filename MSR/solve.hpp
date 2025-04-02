#pragma once


double scalar_product(int n, double* x, double* y, int p, int k);
void matrix_mult_vector_msr(int n, double* A, int* I, double* x, double* y, int p, int k);
void mult_sub_vector(int n, double* x, double* y, double t, int p, int k);
void apply_preconditioner_msr_matrix(int n, double* A, int* I, double* v, double* r, int p, int k);
int min_residual_msr_matrix(int n, double* A, int* I, double* b, double* x, double* r, double* u, double* v
        , double eps, int maxit, int p, int k);
int min_residual_msr_matrix_full(int n, double* A, int* I, double* b, double* x, double* r, double* u, double* v
        , double eps, int maxit, int maxsteps, int p, int k);


