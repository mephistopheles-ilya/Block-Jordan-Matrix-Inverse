#include "mpi.h"

#include <cmath>

#define EPS 1e-16

int inverse_block(double* block, double* inverse, int* permutations, int m, double norm) {
    int i = 0, j = 0;
    double max_elem = 0, tmp_d, inverse_elem;
    int max_col, tmp_i;

    for(i = 0; i < m; ++i) {
        permutations[i] = i;
    }

    for(j = 0; j < m; ++j) {
        max_elem = 0.0;
        max_col = j;
        for(i = j; i < m; ++i) {
            if(fabs(block[j * m + i]) > max_elem) {
                max_col = i;
                max_elem = fabs(block[j * m + i]);
            }
        }
        if (max_elem <= EPS * norm){
            return 1;
        }
        if(max_col != j) {
            //swap elements to remember permutation
            tmp_i = permutations[j];
            permutations[j] = permutations[max_col];
            permutations[max_col] = tmp_i;

            //swap columns in block
            for(int v = 0; v < m; ++v) {
                tmp_d = block[v * m + j];
                block[v * m + j] = block[v * m + max_col];
                block[v * m + max_col] = tmp_d;
            }
        }

        inverse_elem = 1. / block[j * m + j];

        inverse[j * m + j] = inverse_elem;

        for(int v = j + 1; v < m; ++v) {
            block[j * m + v] *= inverse_elem;
        }
        for(int v = 0; v < j; ++v) {
            inverse[j * m + v] *= inverse_elem;
        }
        for(int v = 0; v < j; ++v) {
            for(int u = j + 1; u < m; ++u) {
                block[v * m + u] -=  block[j * m + u] * block[v * m + j];
            }
            for(int u = 0; u <= j; ++u) {
                inverse[v * m + u] -= inverse[j * m + u] * block[v * m + j];
            }
        }
        for(int v = j + 1; v < m; ++v) {
            for(int u = j + 1; u < m; ++u) {
                block[v * m + u] -= block[j * m + u] * block[v * m + j];
            }
            for(int u = 0; u <= j; ++u) {
                inverse[v * m + u] -= inverse[j * m + u] * block[v * m + j];
            }
        }
    }

    for(i = 0; i < m;) {
        if(permutations[i] != i) {
            for(int v = 0; v < m; ++v) {
                tmp_d = inverse[i * m + v];
                inverse[i * m + v] = inverse[permutations[i] * m + v];
                inverse[permutations[i] * m + v] = tmp_d;
            }

            tmp_i = permutations[i];
            permutations[i] = permutations[tmp_i];
            permutations[tmp_i] = tmp_i;
        } else {
            ++i;
        }
    }
    return 0;
}

int inverse_matrix(double* matrix, double* inverse, int n, int m, int* permutations
        , double* block_m, double* inv_block_m, int* permutations_m, double matrix_norm
        , double* matrix_buf, double* inverse_buf, int proc_num, int p, MPI_Comm com) {
    int k = n / m;
    int l = n - m * k;
    int i, j, u, v;
    double block_min_norm, current_block_norm;
    int min_col, tmp_i;
    int el_in_block_line = k * m * m + l * m;
    bool inv_exist;

    return 0;
}


