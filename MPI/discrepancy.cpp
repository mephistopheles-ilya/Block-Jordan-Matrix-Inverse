#include "mpi.h"

#include "read_print_fill.hpp"
#include "matrix.hpp"

#include <cmath>

double calculate_discrepancy(double* matrix, double* inverse, int n, int m, double* tmp_block_m
        , double* tmp_line_n, int proc_num, int p, MPI_Comm comm, double* inverse_buf) {
    int k = n / m;
    int l = n - m * k;
    int el_in_blcok_line = k * m * m + l * m;
    int i_glob = 0, j = 0, q = 0;
    bool small_row = (l != 0 && k%p == proc_num) ? true : false;
    double sum = 0;
    int max_rows = get_max_rows(n, m, p);
    int rows = get_rows(n, m, p, proc_num);
    rows = (small_row == true) ? rows - 1 : rows;
    memset(tmp_line_n, 0, n * sizeof(double));
    for(i_glob = 0; i_glob < k; ++i_glob) {
        for(j = 0; j < rows; ++j) {
            memcpy(inverse_buf + j * m * m, inverse + j * el_in_blcok_line + i_glob * m * m, m * m * sizeof(double));
        }
        if (small_row == true) {
            memcpy(inverse_buf  + j * m * m, inverse + j * el_in_blcok_line + i_glob * m * l, m * l * sizeof(double));
        }
        MPI_Allgather(inverse_buf, max_rows * m * m, MPI_DOUBLE, matrix + max_rows * el_in_blcok_line
                , max_rows * m * m, MPI_DOUBLE, comm);
        for(q = 0; q < rows; ++q) {
            memset(tmp_block_m, 0, m * m * sizeof(double));
            for(j = 0; j < k; ++j) {
                int proc_j = j % p;
                int j_loc = j / p; 
                block_mult_add(matrix + q * el_in_blcok_line + j * m * m, m, m
                        , matrix  + max_rows * el_in_blcok_line + proc_j * max_rows * m * m + j_loc * m * m
                        , m, tmp_block_m);
            }
            int proc_j = j % p;
            int j_loc = j / p; 
            block_mult_add(matrix + q * el_in_blcok_line + j * m * m, m, l
                    , matrix + max_rows * el_in_blcok_line + proc_j * max_rows * m * m + j_loc * m * m
                    , m, tmp_block_m); 
            if ((q * p + proc_num) == i_glob) {
                for(int s = 0; s < m; ++s) {
                    tmp_block_m[s * m + s] -= 1;
                }
            }
            for(int w = 0; w < m; ++w) {
                sum = 0;
                for(int v = 0; v < m; ++v) {
                    sum += fabs(tmp_block_m[v * m + w]);
                }
                tmp_line_n[i_glob * m + w] += sum;
            }
        }
        if (small_row == true) {
            memset(tmp_block_m, 0, m * m * sizeof(double));
            for(j = 0; j < k; ++j) {
                int proc_j = j % p;
                int j_loc = j / p; 
                block_mult_add(matrix + q * el_in_blcok_line + j * m * l, l, m
                        , matrix  + max_rows * el_in_blcok_line + proc_j * max_rows * m * m + j_loc * m * m
                        , m, tmp_block_m);
        
            }
            int proc_j = j % p;
            int j_loc = j / p; 
            block_mult_add(matrix + q * el_in_blcok_line + j * m * l, l, l
                    , matrix + max_rows * el_in_blcok_line + proc_j * max_rows * m * m + j_loc * m * m
                    , m, tmp_block_m); 
            for(int w = 0; w < m; ++w) {
                sum = 0;
                for(int v = 0; v < l; ++v) {
                    sum += fabs(tmp_block_m[v * m + w]);
                }
                tmp_line_n[i_glob * m + w] += sum;
            }
        }
    }

    if (l != 0) {
        for(j = 0; j < rows; ++j) {
            memcpy(inverse_buf + j * l * m, inverse + j * el_in_blcok_line + i_glob * m * m, l * m * sizeof(double));
        }
        if (small_row == true) {
            memcpy(inverse_buf  + j * l * m, inverse + j * el_in_blcok_line + i_glob * m * l, l * l * sizeof(double));
        }

        MPI_Allgather(inverse_buf, max_rows * m * l, MPI_DOUBLE, matrix + max_rows * el_in_blcok_line
                , max_rows * m * l, MPI_DOUBLE, comm);
        for(q = 0; q < rows; ++q) {
            memset(tmp_block_m, 0, m * m * sizeof(double));
            for(j = 0; j < k; ++j) {
                int proc_j = j % p;
                int j_loc = j / p; 
                block_mult_add(matrix + q * el_in_blcok_line + j * m * m, m, m
                        , matrix  + max_rows * el_in_blcok_line + proc_j * max_rows * m * l + j_loc * l * m
                        , l, tmp_block_m);
            }
            int proc_j = j % p;
            int j_loc = j / p; 
            block_mult_add(matrix + q * el_in_blcok_line + j * m * m, m, l
                    , matrix + max_rows * el_in_blcok_line + proc_j * max_rows * l * m + j_loc * l * m
                    , l, tmp_block_m); 
            for(int w = 0; w < l; ++w) {
                sum = 0;
                for(int v = 0; v < m; ++v) {
                    sum += fabs(tmp_block_m[v * l + w]);
                }
                tmp_line_n[i_glob * m + w] += sum;
            }
        }
        if (small_row == true) {
            memset(tmp_block_m, 0, m * m * sizeof(double));
            for(j = 0; j < k; ++j) {
                int proc_j = j % p;
                int j_loc = j / p; 
                block_mult_add(matrix + q * el_in_blcok_line + j * m * l, l, m
                        , matrix  + max_rows * el_in_blcok_line + proc_j * max_rows * l * m + j_loc * l * m
                        , l, tmp_block_m);
        
            }
            int proc_j = j % p;
            int j_loc = j / p; 
            block_mult_add(matrix + q * el_in_blcok_line + j * m * l, l, l
                    , matrix + max_rows * el_in_blcok_line + proc_j * max_rows * l * m + j_loc * l * m
                    , l, tmp_block_m); 
            if ((q * p + proc_num) == i_glob) {
                for(int s = 0; s < l; ++s) {
                    tmp_block_m[s * l + s] -= 1;
                }
            }
            for(int w = 0; w < l; ++w) {
                sum = 0;
                for(int v = 0; v < l; ++v) {
                    sum += fabs(tmp_block_m[v * l + w]);
                }
                tmp_line_n[i_glob * m + w] += sum;
            }
        }
    }
    MPI_Allreduce(tmp_line_n, inverse_buf, n, MPI_DOUBLE, MPI_SUM, comm);

    double max = inverse_buf[0];
    for(int i = 0; i < n; ++i) {
        if(inverse_buf[i] > max) max = inverse_buf[i];
    }
    return max;
}

