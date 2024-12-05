#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include "eigenvalues.hpp"

inline int LIMIT_OF_ITERATIONS = 50000;

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
        if (std::sqrt(s_k) <= eps * norm) {
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

int find_eigenvalues(double *matrix, double *x1, double *x2, double* eigenvalues
        , int n, double eps, double norm, int& its) {
    double shift = 0;
    double s_k = 0;
    double a_norm = 0;
    double x_norm = 0;
    double dot_product = 0;
    for(int i = (n - 1); i >= 2; --i) {
        while(std::fabs(matrix[i * n + (i - 1)]) > eps * norm && (its < LIMIT_OF_ITERATIONS)) {
            ++its;
            if (its % 1000 == 0) {
                printf("Iterations = %d on step = %d\n", its, i);
            }
            shift = matrix[i * n + i] - 0.5 * matrix[i * n + (i - 1)];
            for(int j = 0; j <= i; ++j) {
                matrix[j * n + j] -= shift;
            }
            for(int k = 0; k < i; ++k) {
                if (std::fabs(matrix[n * (k + 1) + k]) <= eps * norm) {
                    x1[k] = 0;
                    x2[k] = 0;
                } else {
                    s_k = matrix[n * (k + 1) + k] * matrix[n * (k + 1) + k];
                    a_norm = std::sqrt(matrix[n * k + k] * matrix[n * k + k] + s_k);
                    x1[k] = matrix[k * n + k] - a_norm;
                    x2[k] = matrix[(k + 1) * n + k];
                    x_norm = std::sqrt(x1[k] * x1[k] + s_k);
                    x1[k] /= x_norm;
                    x2[k] /= x_norm;
                    matrix[k * n + k] = a_norm;
                    matrix[(k + 1) * n + k] = 0;

                    for(int l = k + 1; l <= i; ++l) {
                        dot_product = matrix[k * n + l] * x1[k];
                        dot_product += matrix[(k + 1) * n + l] * x2[k];
                        dot_product *= 2;
                        matrix[k * n + l] -= x1[k] * dot_product;
                        matrix[(k + 1)* n + l] -= x2[k] * dot_product;
                    }

                }
            }
            for(int k = 0; k < i; ++k) {
                for(int l = 0; l < k + 2; ++l) {
                    dot_product = matrix[l * n + k] * x1[k];
                    dot_product += matrix[l * n + k + 1] * x2[k];
                    dot_product *= 2;
                    matrix[l * n + k] -= x1[k] * dot_product;
                    matrix[l * n + k + 1] -= x2[k] * dot_product;
                }
            }
           for(int j = 0; j <= i; ++j) {
                matrix[j * n + j] += shift;
           }
        }
        eigenvalues[i] = matrix[i * n + i];
    }
    if (its >= LIMIT_OF_ITERATIONS) {
        printf("The limit of iterations has been reached\n");
        return 1;
    }
    if (n > 1) {
        double b = -(matrix[0] + matrix[n + 1]);
        double c = (matrix[0] * matrix[n + 1] - matrix[n] * matrix[1]);
        double discr = b * b - 4 * c;
        double val1 = 0, val2 = 0;

        if (discr < 0) {
            printf("No real solutions\n");
            return 1;
        } else {
            discr = std::sqrt(discr);
        }
        if (std::fabs(b) <= eps * norm && std::fabs(c) < eps * norm) {
            val1 = 0;
            val2 = 0;
        } else if (b > 0) {
            val1 = (-b - discr) * 0.5;
            val2 = c / val1;
        } else {
            val1 = (-b + discr) * 0.5;
            val2 = c / val1;
        }
        eigenvalues[1] = std::max(val1, val2);
        eigenvalues[0] = std::min(val1, val2);
    } 
    if (n == 1) {
        eigenvalues[0] = matrix[0];
    }
    return 0;
}










