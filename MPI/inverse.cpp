#include "mpi.h"
#include "matrix.hpp"
#include "inverse.hpp"
#include "read_print_fill.hpp"

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
        , double* matrix_buf, double* inverse_buf, int proc_num, int p, MPI_Comm comm) {
    int k = n / m;
    int l = n - m * k;
    int i_glob = 0, i_loc = 0;
    int el_in_block_line = k * m * m + l * m;
    int inv_exist_loc = 0, inv_exist_glob = 0;
    int send_count_matrix = 0, send_count_inv = 0;
    int min_col_loc = 0, min_col_glob = 0;
    double current_block_norm = 0, block_min_norm_loc = 0, block_min_norm_glob = 0;
    int u = 0, v = 0;

    for(i_glob = 0; i_glob < k; ++i_glob) {
        i_loc = i_glob / p;
        send_count_matrix = (k - i_glob + p - 1) / p;
        send_count_inv = (i_glob + 1 + p - 1) / p;
        MPI_Scatter(matrix + i_loc * el_in_block_line + i_glob * m * m, send_count_matrix, MPI_DOUBLE
                , matrix_buf, send_count_matrix, MPI_DOUBLE, proc_num, comm);
        MPI_Scatter(inverse + i_loc * el_in_block_line, send_count_inv, MPI_DOUBLE
                , inverse_buf, send_count_inv, MPI_DOUBLE, proc_num, comm);
        inv_exist_loc = 0;
        inv_exist_glob = 0;
        min_col_loc = 0;
        block_min_norm_loc = std::numeric_limits<double>::max();
        for(int q = 0; q < send_count_matrix; ++q) {
            memcpy(block_m, matrix_buf, m * m * sizeof(double)); 
            fill_id_block(inv_block_m, m);
            if (inverse_block(block_m, inv_block_m, permutations_m, m, matrix_norm) == 0) {
                current_block_norm = norm_block(inv_block_m, m);
                if (current_block_norm < block_min_norm_loc) {
                    block_min_norm_loc = current_block_norm;
                    min_col_loc = q;
                }
                inv_exist_loc = 1;
            }
        }
        double_int di_loc, di_glob;
        if(inv_exist_loc == 1) {
            di_loc.x = 1/block_min_norm_loc;
        } else {
            di_loc.x = 0;
        }
        di_loc.y = min_col_loc + send_count_matrix * p + i_glob;

        MPI_Allreduce(&di_loc, &di_glob, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
        if (di_glob.x <= 0) {
            return 1;
        }

        if (proc_num == (di_glob.y - i_glob) / send_count_matrix) {
            memcpy(block_m, matrix_buf + m * m * (di_glob.y % send_count_matrix), m * m * sizeof(double));;
        }
        MPI_Bcast(block_m, m * m, MPI_DOUBLE, proc_num, comm);
        fill_id_block(block_m, m);
        inverse_block(block_m, inv_block_m, permutations_m, m, matrix_norm);

        for(int q = 0; q < send_count_matrix; ++q) {
            block_mult(inv_block_m, m, m, matrix_buf + m * m * q, m, block_m);
            memcpy(matrix_buf + m * m * q, block_m, m * m * sizeof(double));
        }

        for(int q = 0; q < send_count_inv; ++q) {
            block_mult(inv_block_m, m, m, inverse_buf + m * m * q, m, block_m);
            memcpy(inverse_buf + m * m * q, block_m, m * m * sizeof(double));
        }

        MPI_Allgather(matrix_buf, send_count_matrix, MPI_DOUBLE, matrix_buf, send_count_matrix, MPI_DOUBLE, comm);
        MPI_Allgather(inverse_buf,send_count_inv, MPI_DOUBLE, inverse_buf, send_count_inv, MPI_DOUBLE, comm);

        

        int rows = get_rows(n, m, p, proc_num);
        for(v = 0; v < rows; ++v) {
            for(u = i_glob + 1; u < k; ++u) {
                block_mult_sub(matrix + v * el_in_block_line + i_glob * m * m
                        , m , m, matrix + i_glob * el_in_block_line + u * m * m, m
                        , matrix + v * el_in_block_line + u * m * m);
            }
            for(u = 0; u <= i_glob; ++u) {
                block_mult_sub(matrix + v * el_in_block_line +  i_glob * m * m
                        , m , m, inverse + i_glob * el_in_block_line + u * m * m, m
                        , inverse + v * el_in_block_line + u * m * m);
            }
        }
    }

    return 0;
}


