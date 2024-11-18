#include <cmath>
#include <math.h>
#include <limits>
#include <string.h>

#include "matrix.hpp"
#include "thread.hpp"

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
            return -1;
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



int inverse_matrix(Arg* a, double* matrix, double* inverse, int n, int m, int* permutations
        , double* block_m, double* inv_block_m, int* permutations_m) {
    double matrix_norm = norm_matrix(matrix, n, m);
    int k = n / m;
    int l = n - m * k;
    int i, j, u, v;
    double block_min_norm, current_block_norm;
    int min_col, tmp_i;
    int el_in_block_line = k * m * m + l * m;
    bool inv_exist;
    for(int i = 0; i < k; ++i) {
        permutations[i] = i;
    }

    for(j = 0; j < k; ++j) {
        min_col = j;
        block_min_norm = std::numeric_limits<double>::max(); 
        inv_exist = false;
        for(i = j; i < k; ++i) {
            memcpy(block_m, matrix + j * el_in_block_line + i * m * m, m * m * sizeof(double));
            fill_id_block(inv_block_m, m);
            if (inverse_block(block_m, inv_block_m, permutations_m, m, matrix_norm)
                    == 0) {
                current_block_norm = norm_block(inv_block_m, m);
                if (current_block_norm < block_min_norm) {
                    block_min_norm = current_block_norm;
                    min_col = i;
                }
                inv_exist = true;
            }
        }
        if (!inv_exist) {
            return -1;
        }
        if (min_col != j) {
            swap_block_col(matrix, n, m, min_col, j, block_m);
            tmp_i = permutations[j];
            permutations[j] = permutations[min_col];
            permutations[min_col] = tmp_i;
        }
        memcpy(block_m, matrix + j * el_in_block_line + j * m * m, m * m * sizeof(double));
        fill_id_block(inv_block_m, m);
        inverse_block(block_m, inv_block_m, permutations_m, m, matrix_norm);
        for(u = j + 1; u < k; ++u) {
            block_mult(inv_block_m, m, m, matrix + j * el_in_block_line 
                    + u * m * m, m, block_m);
            memcpy(matrix + j * el_in_block_line + u * m * m, block_m, m * m * sizeof(double));
        }
        block_mult(inv_block_m, m, m, matrix + j * el_in_block_line + u * m * m, l, block_m); 
        memcpy(matrix + j * el_in_block_line + u * m * m, block_m, m * l * sizeof(double));
        for(u = 0; u < j; ++u) {
            block_mult(inv_block_m, m, m, inverse + j * el_in_block_line 
                    + u * m * m, m, block_m);
            memcpy(inverse + j * el_in_block_line + u * m * m, block_m, m * m * sizeof(double));
        }
        memcpy(inverse + j * el_in_block_line + u * m * m, inv_block_m, m * m * sizeof(double));
        
        for(v = 0; v < j; ++v) {
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
        for(v = j + 1; v < k; ++v) {
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
        if (l != 0) {
            for(u = j + 1; u < k; ++u) {
                block_mult_sub(matrix + v * el_in_block_line + j * m * l
                        , l , m, matrix + j * el_in_block_line + u * m * m, m
                        , matrix + v * el_in_block_line + u * m * l);
            }
            block_mult_sub(matrix + v * el_in_block_line +  j * l * m
                    , l , m, matrix + j * el_in_block_line + u * m * m, l
                    , matrix + v * el_in_block_line + u * m * l);
            for(u = 0; u <= j; ++u) {
                block_mult_sub(matrix + v * el_in_block_line +  j * l * m
                        , l , m, inverse + j * el_in_block_line + u * m * m, m
                        , inverse + v * el_in_block_line + u * m * l);
            }
        }
    }
    if (l != 0) {
        memcpy(block_m, matrix + j * el_in_block_line + j * l * m, l * l * sizeof(double));
        fill_id_block(inv_block_m, l);
        if (inverse_block(block_m, inv_block_m, permutations_m, l, matrix_norm) 
                != 0) {
                return -1;
        }
        for(u = 0; u < j; ++u) {
            block_mult(inv_block_m, l, l, inverse + j * el_in_block_line 
                    + u * m * l, m, block_m);
            memcpy(inverse + j * el_in_block_line + u * m * l, block_m, l * m * sizeof(double));
        }
        memcpy(inverse + j * el_in_block_line + u * m * l, inv_block_m, l * l * sizeof(double));
        
        for(v = 0; v < j; ++v) {
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
    return 0;   
}

