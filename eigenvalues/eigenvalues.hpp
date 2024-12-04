#pragma once

int matrix_to_almost_triangle(double* matrix, double* x, int n, double eps, double norm);
int find_eigenvalues(double* matrix, double* x1, double* x2, double* eigenvalues
        , int n, double eps, double norm, int& its); 
