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
    int inv_exist_loc = 0;
    int send_count_matrix = 0, send_count_inv = 0;
    int min_col_loc = 0, min_col_glob = 0;
    double current_block_norm = 0, block_min_norm_loc = 0;
    int u = 0, v = 0, tmp_i = 0, q = 0;
    int real_num_in_my_part = 0;
    int blocks_in_line = (l == 0) ? k : (k + 1);
    bool small_row = (l != 0 && k%p == proc_num) ? true : false;

    for(i_glob = 0; i_glob < k; ++i_glob) {
        i_loc = i_glob / p;
        send_count_matrix = (blocks_in_line - i_glob + p - 1) / p;
        send_count_inv = (i_glob + 1 + p - 1) / p;
        MPI_Scatter(matrix + i_loc * el_in_block_line + i_glob * m * m, send_count_matrix * m * m, MPI_DOUBLE
                , matrix_buf, send_count_matrix * m * m, MPI_DOUBLE, i_glob%p, comm);
        MPI_Scatter(inverse + i_loc * el_in_block_line, send_count_inv * m * m, MPI_DOUBLE
                , inverse_buf, send_count_inv * m * m, MPI_DOUBLE, i_glob%p, comm);
        inv_exist_loc = 0;
        min_col_loc = 0;
        block_min_norm_loc = std::numeric_limits<double>::max();
        real_num_in_my_part = (blocks_in_line - i_glob - send_count_matrix * (proc_num + 1));
        real_num_in_my_part = (real_num_in_my_part <= 0 && real_num_in_my_part > -send_count_matrix && l != 0) ? 
            (real_num_in_my_part - 1) : real_num_in_my_part;
        real_num_in_my_part = (real_num_in_my_part > 0) ? send_count_matrix
            : send_count_matrix + real_num_in_my_part;
        for(int q = 0; q < real_num_in_my_part; ++q) {
            memcpy(block_m, matrix_buf + q * m * m, m * m * sizeof(double)); 
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
        di_loc.y = min_col_loc + send_count_matrix * proc_num + i_glob;

        MPI_Allreduce(&di_loc, &di_glob, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
        if (di_glob.x <= 0) {
            return 1;
        }
        min_col_glob = di_glob.y;

        if (proc_num == (di_glob.y - i_glob) / send_count_matrix) {
            memcpy(block_m, matrix_buf + m * m * ((di_glob.y - i_glob) % send_count_matrix), m * m * sizeof(double));;
        }
        MPI_Bcast(block_m, m * m, MPI_DOUBLE, (di_glob.y - i_glob)/send_count_matrix, comm);
        fill_id_block(inv_block_m, m);
        inverse_block(block_m, inv_block_m, permutations_m, m, matrix_norm);

        for(q = 0; q < real_num_in_my_part; ++q) {
            block_mult(inv_block_m, m, m, matrix_buf + m * m * q, m, block_m);
            memcpy(matrix_buf + m * m * q, block_m, m * m * sizeof(double));
        }
        real_num_in_my_part = (blocks_in_line - i_glob - send_count_matrix * (proc_num + 1));
        if(real_num_in_my_part <= 0 && real_num_in_my_part > -send_count_matrix && l != 0) {
            block_mult(inv_block_m, m, m, matrix_buf + m * m * q, l, block_m);
            memcpy(matrix_buf + m * m * q, block_m, l * m * sizeof(double));
        }

        real_num_in_my_part = (i_glob + 1 - send_count_inv * (proc_num + 1));
        real_num_in_my_part = (real_num_in_my_part <= 0 && real_num_in_my_part > -send_count_inv && l != 0) ? 
            (real_num_in_my_part - 1) : real_num_in_my_part;
        real_num_in_my_part = (real_num_in_my_part > 0) ? send_count_inv
            : send_count_inv + real_num_in_my_part;

        for(q = 0; q < real_num_in_my_part; ++q) {
            block_mult(inv_block_m, m, m, inverse_buf + m * m * q, m, block_m);
            memcpy(inverse_buf + m * m * q, block_m, m * m * sizeof(double));
        }
        real_num_in_my_part = (i_glob + 1 - send_count_inv * (proc_num + 1));
        if (real_num_in_my_part <= 0 && real_num_in_my_part > -send_count_matrix) {
            memcpy(inverse_buf + m * m * q, inv_block_m, m * m * sizeof(double));
        }


        int max_rows = get_max_rows(n, m, p);
        MPI_Allgather(matrix_buf, send_count_matrix * m * m, MPI_DOUBLE
                , matrix + max_rows * el_in_block_line, send_count_matrix * m * m, MPI_DOUBLE, comm);
        MPI_Allgather(inverse_buf, send_count_inv * m * m, MPI_DOUBLE
                , inverse + max_rows * el_in_block_line, send_count_inv * m * m, MPI_DOUBLE, comm);
        if (i_glob % p == proc_num) {
            memcpy(matrix + i_loc * el_in_block_line + i_glob * m * m, matrix + max_rows * el_in_block_line
                    , ((k - i_glob) * m * m + l * m) * sizeof(double));
            memcpy(inverse + i_loc * el_in_block_line, inverse + max_rows * el_in_block_line
                    , (i_glob + 1) * m * m * sizeof(double));
        }

        if (i_glob != min_col_glob) {
            if (proc_num == 0) {
                printf("SWAP\n");
            }
            swap_block_col(matrix, n, m, i_glob, min_col_glob, proc_num, p, block_m);
            memcpy(matrix + max_rows * el_in_block_line + m * m * (min_col_glob - i_glob)
                    , matrix + max_rows * el_in_block_line, m * m * sizeof(double)); 
            memcpy(matrix + max_rows * el_in_block_line, block_m, m * m * sizeof(double));
            tmp_i = permutations[i_glob];
            permutations[i_glob] = permutations[min_col_glob];
            permutations[min_col_glob] = tmp_i;
        }

        int rows = get_rows(n, m, p, proc_num);
        rows = (small_row == true) ? (rows - 1) : rows;
        for(v = 0; v < rows; ++v) {
            if (i_glob % p != proc_num || i_glob / p != v) {
                for(u = i_glob + 1; u < k; ++u) {
                    block_mult_sub(matrix + v * el_in_block_line + i_glob * m * m
                            , m , m, matrix + max_rows * el_in_block_line + (u - i_glob) * m * m, m
                            , matrix + v * el_in_block_line + u * m * m);
                }
                block_mult_sub(matrix + v * el_in_block_line +  i_glob * m * m
                        , m , m, matrix + max_rows * el_in_block_line + (u - i_glob) * m * m, l
                        , matrix + v * el_in_block_line + u * m * m);
                for(u = 0; u <= i_glob; ++u) {
                    block_mult_sub(matrix + v * el_in_block_line +  i_glob * m * m
                            , m , m, inverse + max_rows * el_in_block_line + u * m * m, m
                            , inverse + v * el_in_block_line + u * m * m);
                }
            }
        }
        if (small_row == true) {
            for(u = i_glob + 1; u < k; ++u) {
                block_mult_sub(matrix + v * el_in_block_line + i_glob * m * l
                        , l , m, matrix + max_rows * el_in_block_line + (u - i_glob) * m * m, m
                        , matrix + v * el_in_block_line + u * m * l);
            }
            block_mult_sub(matrix + v * el_in_block_line +  i_glob * l * m
                    , l , m, matrix + max_rows * el_in_block_line + (u - i_glob) * m * m, l
                    , matrix + v * el_in_block_line + u * m * l);
            for(u = 0; u <= i_glob; ++u) {
                block_mult_sub(matrix + v * el_in_block_line +  i_glob * l * m
                        , l , m, inverse + max_rows * el_in_block_line + u * m * m, m
                        , inverse + v * el_in_block_line + u * m * l);
            }
        }
    }
    if (l != 0) {
        if (proc_num == (k % p)) {
            int rows = get_rows(n, m, p, proc_num);
            memcpy(block_m, matrix + (rows - 1) * el_in_block_line + k * m * l, l * l * sizeof(double));
        }
        MPI_Bcast(block_m, l * l, MPI_DOUBLE, k % p, comm);
        fill_id_block(inv_block_m, l);
        if (inverse_block(block_m, inv_block_m, permutations_m, l, matrix_norm) != 0) {
            return 1;
        }
        i_loc = i_glob / p;
        send_count_inv = (i_glob + 1 + p - 1) / p;

        MPI_Scatter(inverse + i_loc * el_in_block_line, send_count_inv * l * m, MPI_DOUBLE
                , inverse_buf, send_count_inv * l * m, MPI_DOUBLE, i_glob%p, comm);

        real_num_in_my_part = (i_glob + 1 - send_count_inv * (proc_num + 1));
        real_num_in_my_part = (real_num_in_my_part <= 0 && real_num_in_my_part > -send_count_inv && l != 0) ? 
            (real_num_in_my_part - 1) : real_num_in_my_part;
        real_num_in_my_part = (real_num_in_my_part > 0) ? send_count_inv
            : send_count_inv + real_num_in_my_part;


        for(q = 0; q < real_num_in_my_part; ++q) {
            block_mult(inv_block_m, l, l, inverse_buf + l * m * q, m, block_m);
            memcpy(inverse_buf + l * m * q, block_m, l * m * sizeof(double));
        }
        real_num_in_my_part = (i_glob + 1 - send_count_inv * (proc_num + 1));
        if(real_num_in_my_part <= 0 && real_num_in_my_part > -send_count_inv) {
            memcpy(inverse_buf + l * m * q, inv_block_m, l * l * sizeof(double));
        }

        int max_rows = get_max_rows(n, m, p);
        MPI_Allgather(inverse_buf, send_count_inv * m * l, MPI_DOUBLE
                , inverse + max_rows * el_in_block_line, send_count_inv * m * l, MPI_DOUBLE, comm);
        if (i_glob % p == proc_num) {
            memcpy(inverse + i_loc * el_in_block_line, inverse + max_rows * el_in_block_line
                    , (k * l * m + l * m) * sizeof(double));
        }

        int rows = get_rows(n, m, p, proc_num);
        rows = (small_row == true) ? (rows - 1) : rows;
        for(v = 0; v < rows; ++v) {
            for(u = 0; u < i_glob; ++u) {
                block_mult_sub(matrix + v * el_in_block_line +  i_glob * m * m
                        , m , l, inverse + max_rows * el_in_block_line + u * l * m, m
                        , inverse + v * el_in_block_line + u * m * m);
            }
            block_mult_sub(matrix + v * el_in_block_line +  i_glob * m * m
                    , m , l, inverse + max_rows * el_in_block_line + u * l * m, l
                    , inverse + v * el_in_block_line + u * m * m);
        }
    }
    for(int i = 0; i < k;) {
        if(permutations[i] != i) {
            swap_block_line(inverse, n, m, i, permutations[i], proc_num, p, inverse_buf, comm);  
            tmp_i = permutations[i];
            permutations[i] = permutations[tmp_i];
            permutations[tmp_i] = tmp_i;
        } else {
            ++i;
        }
    }

    return 0;
}


