#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include "eigenvalues.hpp"

int matrix_to_almost_triangle(double *matrix, double *x, int n, double eps, double norm) {
    double s_k = 0;
    double a_norm = 0;
    double x_norm = 0;
    double x_inv = 0;
    double dot_product = 0;
    for(int col = 0; col < (n - 2); ++col) {
        s_k = 0;
        a_norm = 0;
        x_norm = 0;
        x_inv = 0;
        for(int lin = col + 2; lin < n; ++lin) {
            s_k += matrix[n * lin + col] * matrix[n * lin + col];
        }
        if (s_k <= eps * norm) {
            continue;
        }
        a_norm = std::sqrt(s_k + matrix[n * (col + 1) + col] * matrix[n * (col + 1) + col]);
        x[col + 1] = matrix[n * (col + 1) + col] - a_norm;
        for(int i = col + 2; i < n; ++i) {
            x[i] = matrix[n * i + col];
        }
        x_norm = std::sqrt(x[col + 1] * x[col + 1] + s_k);
        x_inv = 1./x_norm;
        for(int i = col + 1; i < n; ++i) {
            x[i] *= x_inv;
        }

        matrix[n * (col + 1) + col] = a_norm;
        for(int i = col + 1; i < n; ++i) {
            dot_product = 0;
            for(int j = col + 1; j < n; ++j) {
                dot_product += x[j] * matrix[n * j + i];
            }
            dot_product *= 2;
            for(int j = col + 1; j < n; ++j) {
                matrix[n * j + i] -= x[j] * dot_product;
            }

        }
        for(int i = 0; i < n; ++i) {
            dot_product = 0;
            for(int j = col + 1; j < n; ++j) {
                dot_product += x[j] * matrix[n * i + j];
            }
            dot_product *= 2;
            for(int j = col + 1; j < n; ++j) {
                matrix[n * i + j] -= x[j] * dot_product;
            }
        }

    }
    return 0;
}
        
