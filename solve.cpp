#include <cmath>
#include <math.h>
#include <limits>
#include <string.h>

#include "matrix.h"

#define EPS 1e-15

double norm_block(double* block, int m) {
    double max = 0;
    double current = 0;
    for(int i = 0; i < m; ++i) {
        for(int j = 0; j < m; ++j) {
            current += fabs(block[i * m + j]);
        }
        max = (current > max) ? current : max;
        current = 0;
    }
    return max;
}

double norm_matrix(double* matrix, int n, int m) {
    double max = 0, current = 0;
    int num_in_line1; //amount of lines in block
    int num_in_line2; // size of lines in blcok
    double* cur_pointer;
    int tmp;
    int k = n / m;
    int l = n - m * k;
    int elements_in_blcok_line = k * m * m + l * m;
    int line_of_blocks, block_in_line, line_in_block, i;

    for(line_of_blocks =  0; line_of_blocks <= k; ++line_of_blocks) {
        num_in_line1 = (line_of_blocks < k) ? m : l;
        tmp = line_of_blocks * elements_in_blcok_line;
        for(line_in_block = 0; line_in_block < num_in_line1; ++line_in_block) {
            for(block_in_line = 0; block_in_line <= k; ++block_in_line) {
                num_in_line2 = (block_in_line < k) ? m : l;
                cur_pointer = matrix + tmp + block_in_line * m * num_in_line1
                    + line_in_block * num_in_line2;
                for(i = 0; i < num_in_line2; ++i) {
                    current += fabs(*(cur_pointer + i));
                }
            }
            max = (current > max) ? current : max;
            current = 0;
        }
    }
    return max;
}

int inverse_block(double* block, double* inverse, int* permutations, int m, double norm) {
    int i = 0, j = 0;
    double max_elem = 0, tmp_d, inverse_elem;
    int max_col, tmp_i;

    for(i = 0; i < m; ++i) {
        permutations[i] = i;
    }

    for(j = 0; j < m; ++j) {
        max_elem = 0;
        max_col = j;
        for(i = j; i < m; ++i) {
            if(fabs(block[j * m + i]) > max_elem) {
                max_col = j;
                max_elem = fabs(block[j * m + i]);
            }
        }
        if (max_elem < EPS * norm) {
            return -1;
        }
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

        inverse_elem = 1 / block[j * m + j];

        for(int v = j + 1; v < m; ++v) {
            block[j * m + v] *= inverse_elem;
        }
        for(int v = 0; v < j; ++v) {
            inverse[j * m + v] *= inverse_elem;
        }
        inverse[j * m + j] = inverse_elem;
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


void swap_block_col(double* matrix, int n, int m, int i, int j, double* tmp_block_m) {
    int k = n / m;
    int l = n - k * m;
    int u;
    int el_in_block_line = k * m * m + m * l, imm = i * m * m, jmm = j * m * m;
    for(u = 0; u < k; ++u) {
        memcpy(tmp_block_m, matrix + u * el_in_block_line + imm, m * m * sizeof(double));
        memcpy(matrix + u * el_in_block_line + imm
                , matrix + u * el_in_block_line + jmm, m * m * sizeof(double));
        memcpy(matrix + u * el_in_block_line + jmm, tmp_block_m, m * m * sizeof(double));
    }
    if(l != 0) {
        memcpy(tmp_block_m, matrix + u * el_in_block_line + i * m * l, l * m * sizeof(double));
        memcpy(matrix + u * el_in_block_line + i * m * l
                , matrix + u * el_in_block_line + j * m * l, l * m * sizeof(double));
        memcpy(matrix + u * el_in_block_line + j * l * m, tmp_block_m, l * m * sizeof(double));
    }
} 

void swap_block_lin(double* matrix, int n, int m, int i, int j, double* tmp_block_m) {
    int k = n / m;
    int l = n - k * m;
    int u;
    int iline = i * (k * m * m + m * l), jline = j * (k * m * m + m * l);
    for(u = 0; u < k; ++u) {
        memcpy(tmp_block_m, matrix + iline + u * m * m, m * m * sizeof(double));
        memcpy(matrix + iline + u * m * m
                , matrix + jline + u * m * m, m * m * sizeof(double));
        memcpy(matrix + jline + u * m * m, tmp_block_m, m * m * sizeof(double));
    }
    if(l != 0) {
        memcpy(tmp_block_m, matrix + iline + u * m * m, l * m * sizeof(double));
        memcpy(matrix + iline + u * m * m
                , matrix + jline + u * m * m, l * m * sizeof(double));
        memcpy(matrix + jline + u * m * m, tmp_block_m, l * m * sizeof(double));
    }
} 


void fill_id_block(double* block, int m) {
    for(int i = 0; i < m; ++i) {
        for(int j = 0; j < m; ++j) {
            block[m * i + j] = (i == j) ? 1 : 0;
        }
    }
}


int inverse_matrix(double* matrix, double* inverse, int n, int m, int* permutations
        , double* block_m, double* inv_block_m, int* permutations_m) {
    double matrix_norm = norm_matrix(matrix, n, m);
    int k = n / m;
    int l = n - m * k;
    int i, j, u, v;
    double block_min_norm, current_block_norm;
    int min_col, tmp_i;
    int el_in_block_line = k * m * m + l * m;
    for(int i = 0; i < k; ++i) {
        permutations[i] = i;
    }

    for(j = 0; j < k; ++j) {
        min_col = j;
        block_min_norm = std::numeric_limits<double>::max(); 
        for(i = j; i < k; ++i) {
            memcpy(block_m, matrix + j * el_in_block_line + i * m * m, m * m * sizeof(double));
            fill_id_block(inv_block_m, m);
            if (inverse_block(block_m, inv_block_m, permutations_m, m, matrix_norm) != 0) {
                continue;
            }
            current_block_norm = norm_block(inv_block_m, m);
            if (current_block_norm < block_min_norm) {
                block_min_norm = current_block_norm;
                min_col = i;
            }
        }
        if (block_min_norm >= std::numeric_limits<double>::max()) {
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
        if (inverse_block(block_m, inv_block_m, permutations_m, l, matrix_norm) != 0) {
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


        

            






        
















