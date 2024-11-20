#include <cmath>
#include <math.h>
#include <limits>
#include <string.h>

#include "matrix.hpp"
#include "thread.hpp"
#include "reduce_sum.hpp"

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

int inverse_matrix(Arg* a, int* permutations, double* block_m, double* inv_block_m, int* permutations_m
        , double matrix_norm) {
    int n = a->n;
    int m = a->m;
    int thread_number = a->thread_number;
    int p = a->p;
    double* matrix = a->matrix;
    double* inverse = a->inverse;

    int k = n / m;
    int l = n - m * k;
    int i, j, u, v;
    double block_min_norm, current_block_norm;
    int min_col, tmp_i;
    int el_in_block_line = k * m * m + l * m;
    bool inv_exist;

    for(j = 0; j < k; ++j) {
        min_col = j;
        block_min_norm = std::numeric_limits<double>::max(); 
        inv_exist = false;

        i = thread_number;
        while (i < j ) { i += p; }

        for(; i < k; i += p) {
            memcpy(block_m, matrix + j * el_in_block_line + i * m * m, m * m * sizeof(double));
            fill_id_block(inv_block_m, m);
            if (inverse_block(block_m, inv_block_m, permutations_m, m, matrix_norm) == 0) {
                current_block_norm = norm_block(inv_block_m, m);
                if (current_block_norm < block_min_norm) {
                    block_min_norm = current_block_norm;
                    min_col = i;
                }
                inv_exist = true;
            }
        }
        reduce_main_element(p, block_min_norm, min_col, inv_exist);

        if (inv_exist == false) {
            a -> error = 1;
            return 1;
        }


        if (min_col != j) {
            swap_block_col_th(a, min_col, j, block_m);
            tmp_i = permutations[j];
            permutations[j] = permutations[min_col];
            permutations[min_col] = tmp_i;
        }

        pthread_barrier_wait(a->barrier);

        memcpy(block_m, matrix + j * el_in_block_line + j * m * m, m * m * sizeof(double));
        fill_id_block(inv_block_m, m);
        inverse_block(block_m, inv_block_m, permutations_m, m, matrix_norm);

        u = thread_number;
        while (u < (j + 1) ) { u += p; }

        for(; u < k; u += p) {
            block_mult(inv_block_m, m, m, matrix + j * el_in_block_line 
                    + u * m * m, m, block_m);
            memcpy(matrix + j * el_in_block_line + u * m * m, block_m, m * m * sizeof(double));
        }
        if (l != 0 && (k % p == thread_number)) {
            block_mult(inv_block_m, m, m, matrix + j * el_in_block_line + k * m * m, l, block_m); 
            memcpy(matrix + j * el_in_block_line + k * m * m, block_m, m * l * sizeof(double));
        }

        for(u = thread_number; u < j; u += p) {
            block_mult(inv_block_m, m, m, inverse + j * el_in_block_line 
                    + u * m * m, m, block_m);
            memcpy(inverse + j * el_in_block_line + u * m * m, block_m, m * m * sizeof(double));
        }
        if (j % p == thread_number) {
            memcpy(inverse + j * el_in_block_line + j * m * m, inv_block_m, m * m * sizeof(double));
        }

        pthread_barrier_wait(a->barrier);
        
        for(v = thread_number; v < j; v += p) {
            for(u = j + 1; u < k; ++u) {
                block_mult_sub(matrix + v * el_in_block_line + j * m * m
                        , m , m, matrix + j * el_in_block_line + u * m * m, m
                        , matrix + v * el_in_block_line + u * m * m);
            }
            block_mult_sub(matrix + v * el_in_block_line +  j * m * m
                    , m , m, matrix + j * el_in_block_line + u * m * m, l
                    , matrix + v * el_in_block_line + u * m * m);

            for(u = 0; u <= j; ++u) {
                block_mult_sub(matrix + v * el_in_block_line +  j * m * m
                        , m , m, inverse + j * el_in_block_line + u * m * m, m
                        , inverse + v * el_in_block_line + u * m * m);
            }
        }

        v = thread_number;
        while (v < (j + 1) ) { v += p; }

        for(; v < k; v += p) {
            for(u = j + 1; u < k; ++u) {
                block_mult_sub(matrix + v * el_in_block_line + j * m * m
                        , m , m, matrix + j * el_in_block_line + u * m * m, m
                        , matrix + v * el_in_block_line + u * m * m);
            }
            block_mult_sub(matrix + v * el_in_block_line +  j * m * m
                    , m , m, matrix + j * el_in_block_line + u * m * m, l
                    , matrix + v * el_in_block_line + u * m * m);
            for(u = 0; u <= j; ++u) {
                block_mult_sub(matrix + v * el_in_block_line +  j * m * m
                        , m , m, inverse + j * el_in_block_line + u * m * m, m
                        , inverse + v * el_in_block_line + u * m * m);
            }
        }
        if (l != 0 && (k % p == thread_number)) {
            for(u = j + 1; u < k; ++u) {
                block_mult_sub(matrix + k * el_in_block_line + j * m * l
                        , l , m, matrix + j * el_in_block_line + u * m * m, m
                        , matrix + k * el_in_block_line + u * m * l);
            }
            block_mult_sub(matrix + k * el_in_block_line +  j * l * m
                    , l , m, matrix + j * el_in_block_line + u * m * m, l
                    , matrix + k * el_in_block_line + u * m * l);
            for(u = 0; u <= j; ++u) {
                block_mult_sub(matrix + k * el_in_block_line +  j * l * m
                        , l , m, inverse + j * el_in_block_line + u * m * m, m
                        , inverse + k * el_in_block_line + u * m * l);
            }
        }

        pthread_barrier_wait(a->barrier);
    }
    if (l != 0) {
        memcpy(block_m, matrix + j * el_in_block_line + j * l * m, l * l * sizeof(double));
        fill_id_block(inv_block_m, l);
        if (inverse_block(block_m, inv_block_m, permutations_m, l, matrix_norm) != 0) {
            a->error = 1;
            return -1;
        }
        for(u = thread_number; u < j; u += p) {
            block_mult(inv_block_m, l, l, inverse + j * el_in_block_line 
                    + u * m * l, m, block_m);
            memcpy(inverse + j * el_in_block_line + u * m * l, block_m, l * m * sizeof(double));
        }
        if (k % p == thread_number) {
            memcpy(inverse + j * el_in_block_line + j * m * l, inv_block_m, l * l * sizeof(double));
        }

        pthread_barrier_wait(a->barrier);

        for(v = thread_number; v < j; v += p) {
            for(u = 0; u < j; ++u) {
                block_mult_sub(matrix + v * el_in_block_line +  j * m * m
                        , m , l, inverse + j * el_in_block_line + u * l * m, m
                        , inverse + v * el_in_block_line + u * m * m);
            }
            block_mult_sub(matrix + v * el_in_block_line +  j * m * m
                    , m , l, inverse + j * el_in_block_line + u * l * m, l
                    , inverse + v * el_in_block_line + u * m * m);
        }
    }

    pthread_barrier_wait(a->barrier);

    if (thread_number == 0) {
        for(i = 0; i < k;) {
            if(permutations[i] != i) {
                swap_block_lin(inverse, n, m, i, permutations[i], block_m);
                tmp_i = permutations[i];
                permutations[i] = permutations[tmp_i];
                permutations[tmp_i] = tmp_i;
            } else {
                ++i;
            }
        }
    }
    return 0;   
}

#if 0
int inverse_matrix(Arg* a, int* permutations, double* block_m, double* inv_block_m, int* permutations_m
        , double matrix_norm) {
    int n = a->n;
    int m = a->m;
    int thread_number = a->thread_number;
    int p = a->p;
    double* matrix = a->matrix;
    double* inverse = a->inverse;

    int k = n / m;
    int l = n - m * k;
    int i, j, u, v;
    double block_min_norm, current_block_norm;
    int min_col, tmp_i;
    int el_in_block_line = k * m * m + l * m;
    bool inv_exist;

    for(j = 0; j < k; ++j) {
        min_col = j;
        block_min_norm = std::numeric_limits<double>::max(); 
        inv_exist = false;

        i = thread_number;
        while (i < j ) { i += p; }

        for(; i < k; i += p) {
            memcpy(block_m, matrix + j * el_in_block_line + i * m * m, m * m * sizeof(double));
            fill_id_block(inv_block_m, m);
            if (inverse_block(block_m, inv_block_m, permutations_m, m, matrix_norm) == 0) {
                current_block_norm = norm_block(inv_block_m, m);
                if (current_block_norm < block_min_norm) {
                    block_min_norm = current_block_norm;
                    min_col = i;
                }
                inv_exist = true;
            }
        }
        reduce_main_element(p, block_min_norm, min_col, inv_exist);

        if (inv_exist == false) {
            a -> error = 1;
            return 1;
        }


        if (min_col != j) {
            swap_block_col_th(a, min_col, j, block_m);
            tmp_i = permutations[j];
            permutations[j] = permutations[min_col];
            permutations[min_col] = tmp_i;
        }

        pthread_barrier_wait(a->barrier);

        memcpy(block_m, matrix + j * el_in_block_line + j * m * m, m * m * sizeof(double));
        fill_id_block(inv_block_m, m);
        inverse_block(block_m, inv_block_m, permutations_m, m, matrix_norm);

        u = thread_number;
        while (u < (j + 1) ) { u += p; }

        for(; u < k; u += p) {
            memcpy(tmp_block, matrix + j * el_in_block_line + u * m * m, m * m * sizeof(double));
            block_mult(inv_block_m, m, m, tmp_block, m, block_m);
            memcpy(matrix + j * el_in_block_line + u * m * m, block_m, m * m * sizeof(double));
        }
        if (l != 0 && (k % p == thread_number)) {
            block_mult(inv_block_m, m, m, matrix + j * el_in_block_line + k * m * m, l, block_m); 
            memcpy(matrix + j * el_in_block_line + k * m * m, block_m, m * l * sizeof(double));
        }

        for(u = thread_number; u < j; u += p) {
            memcpy(tmp_block, inverse + j * el_in_block_line + u * m * m, m * m * sizeof(double));
            block_mult(inv_block_m, m, m, tmp_block, m, block_m);
            memcpy(inverse + j * el_in_block_line + u * m * m, block_m, m * m * sizeof(double));
        }
        if (j % p == thread_number) {
            memcpy(inverse + j * el_in_block_line + j * m * m, inv_block_m, m * m * sizeof(double));
        }

        pthread_barrier_wait(a->barrier);
        
        for(v = thread_number; v < j; v += p) {
            for(u = j + 1; u < k; ++u) {
                memcpy(tmp_block, matrix + j * el_in_block_line + u * m * m, m * m * sizeof(double));
                block_mult_sub(matrix + v * el_in_block_line + j * m * m
                        , m , m, tmp_block, m
                        , matrix + v * el_in_block_line + u * m * m);
            }
            block_mult_sub(matrix + v * el_in_block_line +  j * m * m
                    , m , m, matrix + j * el_in_block_line + u * m * m, l
                    , matrix + v * el_in_block_line + u * m * m);

            for(u = 0; u <= j; ++u) {
                memcpy(tmp_block, inverse + j * el_in_block_line + u * m * m, m * m * sizeof(double));
                block_mult_sub(matrix + v * el_in_block_line +  j * m * m
                        , m , m, tmp_block, m
                        , inverse + v * el_in_block_line + u * m * m);
            }
        }

        v = thread_number;
        while (v < (j + 1) ) { v += p; }

        for(; v < k; v += p) {
            for(u = j + 1; u < k; ++u) {
                memcpy(tmp_block, matrix + j * el_in_block_line + u * m * m, m * m * sizeof(double));
                block_mult_sub(matrix + v * el_in_block_line + j * m * m
                        , m , m, tmp_block, m
                        , matrix + v * el_in_block_line + u * m * m);
            }
            block_mult_sub(matrix + v * el_in_block_line +  j * m * m
                    , m , m, matrix + j * el_in_block_line + u * m * m, l
                    , matrix + v * el_in_block_line + u * m * m);
            for(u = 0; u <= j; ++u) {
                memcpy(tmp_block, inverse + j * el_in_block_line + u * m * m, m * m * sizeof(double));
                block_mult_sub(matrix + v * el_in_block_line +  j * m * m
                        , m , m, tmp_block, m
                        , inverse + v * el_in_block_line + u * m * m);
            }
        }
        if (l != 0 && (k % p == thread_number)) {
            for(u = j + 1; u < k; ++u) {
                block_mult_sub(matrix + k * el_in_block_line + j * m * l
                        , l , m, matrix + j * el_in_block_line + u * m * m, m
                        , matrix + k * el_in_block_line + u * m * l);
            }
            block_mult_sub(matrix + k * el_in_block_line +  j * l * m
                    , l , m, matrix + j * el_in_block_line + u * m * m, l
                    , matrix + k * el_in_block_line + u * m * l);
            for(u = 0; u <= j; ++u) {
                block_mult_sub(matrix + k * el_in_block_line +  j * l * m
                        , l , m, inverse + j * el_in_block_line + u * m * m, m
                        , inverse + k * el_in_block_line + u * m * l);
            }
        }

        pthread_barrier_wait(a->barrier);
    }
    if (l != 0) {
        memcpy(block_m, matrix + j * el_in_block_line + j * l * m, l * l * sizeof(double));
        fill_id_block(inv_block_m, l);
        if (inverse_block(block_m, inv_block_m, permutations_m, l, matrix_norm) != 0) {
            a->error = 1;
            return -1;
        }
        for(u = thread_number; u < j; u += p) {
            block_mult(inv_block_m, l, l, inverse + j * el_in_block_line 
                    + u * m * l, m, block_m);
            memcpy(inverse + j * el_in_block_line + u * m * l, block_m, l * m * sizeof(double));
        }
        if (k % p == thread_number) {
            memcpy(inverse + j * el_in_block_line + j * m * l, inv_block_m, l * l * sizeof(double));
        }

        pthread_barrier_wait(a->barrier);

        for(v = thread_number; v < j; v += p) {
            for(u = 0; u < j; ++u) {
                block_mult_sub(matrix + v * el_in_block_line +  j * m * m
                        , m , l, inverse + j * el_in_block_line + u * l * m, m
                        , inverse + v * el_in_block_line + u * m * m);
            }
            block_mult_sub(matrix + v * el_in_block_line +  j * m * m
                    , m , l, inverse + j * el_in_block_line + u * l * m, l
                    , inverse + v * el_in_block_line + u * m * m);
        }
        pthread_barrier_wait(a->barrier);
    }


    if (thread_number == 0) {
        for(i = 0; i < k;) {
            if(permutations[i] != i) {
                swap_block_lin(inverse, n, m, i, permutations[i], block_m);
                tmp_i = permutations[i];
                permutations[i] = permutations[tmp_i];
                permutations[tmp_i] = tmp_i;
            } else {
                ++i;
            }
        }
    }
    return 0;   
}
#endif
