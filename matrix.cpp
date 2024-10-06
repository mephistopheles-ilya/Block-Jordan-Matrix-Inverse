#include <string.h>
#include <math.h>

void block_mult(double* A, int av, int ag, double* B, int bg, double* C) {
#if 1
    int av_reminder = av % 3;
    int bg_reminder = bg % 3;
    int _av = av - av_reminder;
    int _bg = bg - bg_reminder;
    int i = 0, j = 0, q = 0, tmp_j = 0;
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
            C[bg * (i) + (j + 2)]=  c02;
            C[bg * (i + 1) + (j)] = c10;
            C[bg * (i + 1) + (j + 1)] = c11;
            C[bg * (i + 1) + (j + 2)] = c12;
            C[bg * (i + 2) + (j)] = c20;
            C[bg * (i + 2) + (j + 1)] = c21;
            C[bg * (i + 2) + (j + 2)] = c22;
        }
    }
    tmp_j = j;

    for(;i < av; ++i) {
        for(j = 0; j < bg; ++j) {
            c00 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * i + q] * B[bg * q + j];
            }
            C[bg * i + j] = c00;
        }
    }

    for(i = 0; i < _av; ++i) {
        for(j = tmp_j; j < bg; ++j) {
            c00 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * i + q] * B[bg * q + j];
            }
            C[bg * i + j] = c00;
        }
    }
#endif
#if 0
    for(int i = 0; i < av; ++i) {
        for(int j = 0; j < bg; ++j) {
            double c00 = 0;
            for(int q = 0; q < ag; ++q) {
                c00 += A[ag * i + q] * B[bg * q + j];
            }
            C[bg * i + j] = c00;
        }
    }
#endif

}

void block_mult_add(double* A, int av, int ag, double* B, int bg, double* C) {
#if 1
    int av_reminder = av % 3;
    int bg_reminder = bg % 3;
    int _av = av - av_reminder;
    int _bg = bg - bg_reminder;
    int i = 0, j = 0, q = 0, tmp_j = 0;
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
            C[bg * (i) + (j + 2)] += c02;
            C[bg * (i + 1) + (j)] += c10;
            C[bg * (i + 1) + (j + 1)] += c11;
            C[bg * (i + 1) + (j + 2)] += c12;
            C[bg * (i + 2) + (j)] += c20;
            C[bg * (i + 2) + (j + 1)] += c21;
            C[bg * (i + 2) + (j + 2)] += c22;
        }
    }
    tmp_j = j;

    for(;i < av; ++i) {
        for(j = 0; j < bg; ++j) {
            c00 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * i + q] * B[bg * q + j];
            }
            C[bg * i + j] += c00;
        }
    }

    for(i = 0; i < _av; ++i) {
        for(j = tmp_j; j < bg; ++j) {
            c00 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * i + q] * B[bg * q + j];
            }
            C[bg * i + j] += c00;
        }
    }
#endif
#if 0
    for(int i = 0; i < av; ++i) {
        for(int j = 0; j < bg; ++j) {
            double c00 = 0;
            for(int q = 0; q < ag; ++q) {
                c00 += A[ag * i + q] * B[bg * q + j];
            }
            C[bg * i + j] += c00;
        }
    }
#endif

}

void block_mult_sub(double* A, int av, int ag, double* B, int bg, double* C) {
#if 1
    int av_reminder = av % 3;
    int bg_reminder = bg % 3;
    int _av = av - av_reminder;
    int _bg = bg - bg_reminder;
    int i = 0, j = 0, q = 0, tmp_j = 0;
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
            C[bg * (i) + (j + 2)] -= c02;
            C[bg * (i + 1) + (j)] -= c10;
            C[bg * (i + 1) + (j + 1)] -= c11;
            C[bg * (i + 1) + (j + 2)] -= c12;
            C[bg * (i + 2) + (j)] -= c20;
            C[bg * (i + 2) + (j + 1)] -= c21;
            C[bg * (i + 2) + (j + 2)] -= c22;
        }
    }
    tmp_j = j;

    for(;i < av; ++i) {
        for(j = 0; j < bg; ++j) {
            c00 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * i + q] * B[bg * q + j];
            }
            C[bg * i + j] -= c00;
        }
    }

    for(i = 0; i < _av; ++i) {
        for(j = tmp_j; j < bg; ++j) {
            c00 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * i + q] * B[bg * q + j];
            }
            C[bg * i + j] -= c00;
        }
    }
#endif
#if 0
    for(int i = 0; i < av; ++i) {
        for(int j = 0; j < bg; ++j) {
            double c00 = 0;
            for(int q = 0; q < ag; ++q) {
                c00 += A[ag * i + q] * B[bg * q + j];
            }
            C[bg * i + j] -= c00;
        }
    }
#endif

}


double calculate_discrepancy(double* matrix, double* inverse, int n, int m, double* tmp_block_m
        , double* tmp_line_n) {
    int k = n / m;
    int l = n - m * k;
    int i, j, q;
    int real_sz = (l != 0) ? (k + 1) : k;
    int elements_in_blcok_line = k * m * m + l * m;
    double sum = 0;
    double max = 0;
    memset(tmp_line_n, 0, n * sizeof(double)); 
    for(i = 0; i < k; ++i) { // columns in inverse
        for(j = 0; j < k; ++j) { // rows in matrix
            memset(tmp_block_m, 0, sizeof(double) * m * m);
            for(q = 0; q < k; ++q) {
                block_mult_add(matrix + j * elements_in_blcok_line + q * m * m
                        , m, m, inverse + q * elements_in_blcok_line + i * m * m, m, tmp_block_m);
            }
            //last block
            block_mult_add(matrix + j * elements_in_blcok_line + q * m * m
                    , m, l, inverse + q * elements_in_blcok_line + i * l * m, m, tmp_block_m);
            if (i == j) {
                for(int s = 0; s < m; ++s) {
                    tmp_block_m[s * m + s] -= 1;
                }
            }
            for(int w = 0; w < m; ++w) {
                sum = 0;
                for(int v = 0; v < m; ++v) {
                    sum += fabs(tmp_block_m[v * m + w]);
                }
                tmp_line_n[i * m + w] += sum;
            }   
        }
        //last row in matrix
        for(; j < real_sz; ++j) {
            memset(tmp_block_m, 0, sizeof(double) * m * m);
            for(q = 0; q < k; ++q) {
                block_mult_add(matrix + j * elements_in_blcok_line + q * l * m
                        , l , m, inverse + q * elements_in_blcok_line + i * m * m, m, tmp_block_m);
            }
            //last block
            block_mult_add(matrix + j * elements_in_blcok_line + q * l * m
                    , l, l, inverse + q * elements_in_blcok_line + i * l * m, m, tmp_block_m);
            for(int w = 0; w < m; ++w) {
                sum = 0;
                for(int v = 0; v < l; ++v) {
                    sum += fabs(tmp_block_m[v * m + w]);
                }
                tmp_line_n[i * m + w] += sum;
            }   

        }
    }
    //last colomn in inverse
    for(; i < real_sz; ++i) {
        for(j = 0; j < k; ++j) {
            memset(tmp_block_m, 0, sizeof(double) * m * m);
            for(q = 0; q < k; ++q) {
                block_mult_add(matrix + j * elements_in_blcok_line + q * m * m
                        , m, m, inverse + q * elements_in_blcok_line + i * m * m, l, tmp_block_m);
            }
            //last block
            block_mult_add(matrix + j * elements_in_blcok_line + q * m * m
                    , m, l, inverse + q * elements_in_blcok_line + i * l * m, l, tmp_block_m);
            for(int w = 0; w < l; ++w) {
                sum = 0;
                for(int v = 0; v < m; ++v) {
                    sum += fabs(tmp_block_m[v * l + w]);
                }
                tmp_line_n[i * m + w] += sum;
            }   
        }
        //last row in matrix
        for(; j < real_sz; ++j) {
            memset(tmp_block_m, 0, sizeof(double) * m * m);
            for(q = 0; q < k; ++q) {
                block_mult_add(matrix + j * elements_in_blcok_line + q * l * m
                        , l , m, inverse + q * elements_in_blcok_line + i * m * m, l, tmp_block_m);
            }
            //last block
            block_mult_add(matrix + j * elements_in_blcok_line + q * l * m
                    , l, l, inverse + q * elements_in_blcok_line + i * l * m, l, tmp_block_m);
            if (i == j) {
                for(int s = 0; s < l; ++s) {
                    tmp_block_m[s * l + s] -= 1;
                }
            }
            for(int w = 0; w < l; ++w) {
                sum = 0;
                for(int v = 0; v < l; ++v) {
                    sum += fabs(tmp_block_m[v * l + w]);
                }
                tmp_line_n[i * m + w] += sum;
            }   
        }
    }
    max = tmp_line_n[0];
    for(i = 0; i < n; ++i) {
        if(tmp_line_n[i] > max) max = tmp_line_n[i];
    }
    return max;
}









