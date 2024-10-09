#include <math.h>
#include <string.h>

//#pragma GCC push_options
//#pragma GCC optimize ("-ffast-math")
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
//#pragma GCC pop_options

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




