#include "matrix.hpp"
#include "read_print_fill.hpp"

#include <cmath>
#include <cstring>

void block_mult(double* A, int av, int ag, double* B, int bg, double* C) {
    int av_reminder = av % 3;
    int bg_reminder = bg % 3;
    int _av = av - av_reminder;
    int _bg = bg - bg_reminder;
    int i = 0, j = 0, q = 0;
    double c00, c01, c02, c10, c11, c12, c20, c21, c22;
    for(i = 0; i < _av; i += 3) {
        for(j = 0; j < _bg; j+= 3) {
            c00 = 0; c01 = 0; c02 = 0;
            c10 = 0; c11 = 0; c12 = 0;
            c20 = 0; c21 = 0; c22 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
                c01 += A[ag * (i) + q] * B[bg * q + (j + 1)];
                c02 += A[ag * (i) + q] * B[bg * q + (j + 2)];
                c10 += A[ag * (i + 1) + q] * B[bg * q + (j)];
                c11 += A[ag * (i + 1) + q] * B[bg * q + (j + 1)];
                c12 += A[ag * (i + 1) + q] * B[bg * q + (j + 2)];
                c20 += A[ag * (i + 2) + q] * B[bg * q + (j)];
                c21 += A[ag * (i + 2) + q] * B[bg * q + (j + 1)];
                c22 += A[ag * (i + 2) + q] * B[bg * q + (j + 2)];
            }
            C[bg * (i) + (j)] = c00;
            C[bg * (i) + (j + 1)] = c01;
            C[bg * (i) + (j + 2)] = c02;
            C[bg * (i + 1) + (j)] = c10;
            C[bg * (i + 1) + (j + 1)] = c11;
            C[bg * (i + 1) + (j + 2)] = c12;
            C[bg * (i + 2) + (j)] = c20;
            C[bg * (i + 2) + (j + 1)] = c21;
            C[bg * (i + 2) + (j + 2)] = c22;
        }
        for(;j < bg; ++j) {
            c00 = 0;
            c10 = 0;
            c20 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
                c10 += A[ag * (i + 1) + q] * B[bg * q + (j)];
                c20 += A[ag * (i + 2) + q] * B[bg * q + (j)];
            }
            C[bg * (i) + (j)] = c00;
            C[bg * (i + 1) + (j)] = c10;
            C[bg * (i + 2) + (j)] = c20;
        }
    }
    for(; i < av; ++i) {
        for(j = 0; j < _bg; j+= 3) {
            c00 = 0; c01 = 0; c02 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
                c01 += A[ag * (i) + q] * B[bg * q + (j + 1)];
                c02 += A[ag * (i) + q] * B[bg * q + (j + 2)];
            }
            C[bg * (i) + (j)] = c00;
            C[bg * (i) + (j + 1)] = c01;
            C[bg * (i) + (j + 2)] = c02;
        }
        for(;j < bg; ++j) {
            c00 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
            }
            C[bg * (i) + (j)] = c00;
        }
    }
}

void block_mult_sub(double* A, int av, int ag, double* B, int bg, double* C) {
    int av_reminder = av % 3;
    int bg_reminder = bg % 3;
    int _av = av - av_reminder;
    int _bg = bg - bg_reminder;
    int i = 0, j = 0, q = 0;
    double c00, c01, c02, c10, c11, c12, c20, c21, c22;
    for(i = 0; i < _av; i += 3) {
        for(j = 0; j < _bg; j+= 3) {
            c00 = 0; c01 = 0; c02 = 0;
            c10 = 0; c11 = 0; c12 = 0;
            c20 = 0; c21 = 0; c22 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
                c01 += A[ag * (i) + q] * B[bg * q + (j + 1)];
                c02 += A[ag * (i) + q] * B[bg * q + (j + 2)];
                c10 += A[ag * (i + 1) + q] * B[bg * q + (j)];
                c11 += A[ag * (i + 1) + q] * B[bg * q + (j + 1)];
                c12 += A[ag * (i + 1) + q] * B[bg * q + (j + 2)];
                c20 += A[ag * (i + 2) + q] * B[bg * q + (j)];
                c21 += A[ag * (i + 2) + q] * B[bg * q + (j + 1)];
                c22 += A[ag * (i + 2) + q] * B[bg * q + (j + 2)];
            }
            C[bg * (i) + (j)] -= c00;
            C[bg * (i) + (j + 1)] -= c01;
            C[bg * (i) + (j + 2)] -=c02;
            C[bg * (i + 1) + (j)] -= c10;
            C[bg * (i + 1) + (j + 1)] -= c11;
            C[bg * (i + 1) + (j + 2)] -= c12;
            C[bg * (i + 2) + (j)] -= c20;
            C[bg * (i + 2) + (j + 1)] -= c21;
            C[bg * (i + 2) + (j + 2)] -= c22;
        }
        for(;j < bg; ++j) {
            c00 = 0;
            c10 = 0;
            c20 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
                c10 += A[ag * (i + 1) + q] * B[bg * q + (j)];
                c20 += A[ag * (i + 2) + q] * B[bg * q + (j)];
            }
            C[bg * (i) + (j)] -= c00;
            C[bg * (i + 1) + (j)] -= c10;
            C[bg * (i + 2) + (j)] -= c20;
        }
    }
    for(; i < av; ++i) {
        for(j = 0; j < _bg; j+= 3) {
            c00 = 0; c01 = 0; c02 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
                c01 += A[ag * (i) + q] * B[bg * q + (j + 1)];
                c02 += A[ag * (i) + q] * B[bg * q + (j + 2)];
            }
            C[bg * (i) + (j)] -= c00;
            C[bg * (i) + (j + 1)] -= c01;
            C[bg * (i) + (j + 2)] -= c02;
        }
        for(;j < bg; ++j) {
            c00 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
            }
            C[bg * (i) + (j)] -= c00;
        }
    }
}


void block_mult_add(double* A, int av, int ag, double* B, int bg, double* C) {
    int av_reminder = av % 3;
    int bg_reminder = bg % 3;
    int _av = av - av_reminder;
    int _bg = bg - bg_reminder;
    int i = 0, j = 0, q = 0;
    double c00, c01, c02, c10, c11, c12, c20, c21, c22;
    for(i = 0; i < _av; i += 3) {
        for(j = 0; j < _bg; j+= 3) {
            c00 = 0; c01 = 0; c02 = 0;
            c10 = 0; c11 = 0; c12 = 0;
            c20 = 0; c21 = 0; c22 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
                c01 += A[ag * (i) + q] * B[bg * q + (j + 1)];
                c02 += A[ag * (i) + q] * B[bg * q + (j + 2)];
                c10 += A[ag * (i + 1) + q] * B[bg * q + (j)];
                c11 += A[ag * (i + 1) + q] * B[bg * q + (j + 1)];
                c12 += A[ag * (i + 1) + q] * B[bg * q + (j + 2)];
                c20 += A[ag * (i + 2) + q] * B[bg * q + (j)];
                c21 += A[ag * (i + 2) + q] * B[bg * q + (j + 1)];
                c22 += A[ag * (i + 2) + q] * B[bg * q + (j + 2)];
            }
            C[bg * (i) + (j)] += c00;
            C[bg * (i) + (j + 1)] += c01;
            C[bg * (i) + (j + 2)] +=c02;
            C[bg * (i + 1) + (j)] += c10;
            C[bg * (i + 1) + (j + 1)] += c11;
            C[bg * (i + 1) + (j + 2)] += c12;
            C[bg * (i + 2) + (j)] += c20;
            C[bg * (i + 2) + (j + 1)] += c21;
            C[bg * (i + 2) + (j + 2)] += c22;
        }
        for(;j < bg; ++j) {
            c00 = 0;
            c10 = 0;
            c20 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
                c10 += A[ag * (i + 1) + q] * B[bg * q + (j)];
                c20 += A[ag * (i + 2) + q] * B[bg * q + (j)];
            }
            C[bg * (i) + (j)] += c00;
            C[bg * (i + 1) + (j)] += c10;
            C[bg * (i + 2) + (j)] += c20;
        }
    }
    for(; i < av; ++i) {
        for(j = 0; j < _bg; j+= 3) {
            c00 = 0; c01 = 0; c02 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
                c01 += A[ag * (i) + q] * B[bg * q + (j + 1)];
                c02 += A[ag * (i) + q] * B[bg * q + (j + 2)];
            }
            C[bg * (i) + (j)] += c00;
            C[bg * (i) + (j + 1)] += c01;
            C[bg * (i) + (j + 2)] += c02;
        }
        for(;j < bg; ++j) {
            c00 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
            }
            C[bg * (i) + (j)] += c00;
        }
    }
}

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

double norm_matrix(double* matrix, int n, int m, int p, int proc_num, MPI_Comm comm) {
    int rows = get_rows(n, m, p, proc_num);

    int k = n/m;
    int l = n - m * k;
    int b = (n + m - 1) / m;

    int line_of_blocks = 0, line_in_block = 0, block_in_line = 0;
    double norm_loc = 0, tmp = 0, norm_glob = 0;

    for(line_of_blocks = 0; line_of_blocks < rows; ++line_of_blocks) {
        int num_line_in_block = m;
        if(l != 0 && (b - 1)%p == proc_num && line_of_blocks == (rows - 1)) {
            num_line_in_block = l;
        }
        for(line_in_block = 0; line_in_block < num_line_in_block; ++line_in_block) {
            tmp = 0;
            for(block_in_line = 0; block_in_line <= k; ++block_in_line) {
                int num_elem_in_line = (block_in_line < k) ? m : l;
                int shift = line_of_blocks * (m * m * k + l * m) + line_in_block * num_elem_in_line + 
                    block_in_line  * m * num_line_in_block;
                for(int q = 0; q < num_elem_in_line; ++q) {
                    tmp += std::fabs(*(matrix + shift + q));
                }
            }
            if (tmp > norm_loc) {
                norm_loc = tmp;
            }
        }
    }
    MPI_Allreduce(&norm_loc, &norm_glob, 1, MPI_DOUBLE, MPI_MAX, comm);
    return norm_glob;
}

void fill_id_block(double* block, int m) {
    for(int i = 0; i < m; ++i) {
        for(int j = 0; j < m; ++j) {
            block[m * i + j] = (i == j) ? 1 : 0;
        }
    }
}

void swap_block_col(double* matrix, int n, int m, int i, int j, int proc_num, int p, double* tmp_block) {
    int rows = get_rows(n, m, p, proc_num);
    int line_of_blocks = 0;
    int k = n/m;
    int l = n - m * k;
    int b = (n + m - 1) / m;
    int el_in_block_line = k * m * m + m * l;

    for(line_of_blocks = 0; line_of_blocks < rows; ++line_of_blocks) {
        int num_line_in_block = m;
        if(l != 0 && (b - 1)%p == proc_num && line_of_blocks == (rows - 1)) {
            num_line_in_block = l;
        }
        memcpy(tmp_block, matrix + line_of_blocks * el_in_block_line + i * m * num_line_in_block
                , m * num_line_in_block * sizeof(double));
        memcpy(matrix + line_of_blocks * el_in_block_line + i * m * num_line_in_block
                , matrix  + line_of_blocks * el_in_block_line + j * m * num_line_in_block
                , m * num_line_in_block * sizeof(double));
        memcpy(matrix + line_of_blocks * el_in_block_line + j * m * num_line_in_block
                , tmp_block, m * num_line_in_block * sizeof(double));
    }
}

void swap_block_line(double* matrix, int n, int m, int i_glob, int j_glob, int proc_num, int p, double* buf
        , MPI_Comm comm) {
    int k = n/m;
    int l = n - m * k;
    int el_in_block_line = m * m * k + l * m;
    if (proc_num == i_glob%p) {
        MPI_Status st;
        int i_loc = i_glob/p;
        MPI_Sendrecv(matrix + i_loc * el_in_block_line, el_in_block_line, MPI_DOUBLE
                , j_glob%p, 0, buf, el_in_block_line, MPI_DOUBLE, i_glob%p, 0, comm, &st);
        memcpy(matrix + i_loc * el_in_block_line, buf, el_in_block_line * sizeof(double));
    }
    if (proc_num == j_glob%p) {
        MPI_Status st;
        int j_loc = j_glob/p;
        MPI_Sendrecv(matrix + j_loc * el_in_block_line, el_in_block_line, MPI_DOUBLE
                , i_glob%p, 0, buf, el_in_block_line, MPI_DOUBLE, j_glob%p, 0, comm, &st);
        memcpy(matrix + j_loc * el_in_block_line, buf, el_in_block_line * sizeof(double));
    }
}







