#include "matrix.hpp"
#include "thread.hpp"
#include <string.h>
#include <math.h>

#include "reduce_sum.hpp"

#if 1
double calculate_discrepancy(Arg* a, double* matrix, double* inverse, double* tmp_block_m, double* tmp_line_n) {
    int n = a->n;
    int m = a->m;
    int p = a->p;
    int thread_number = a->thread_number;
    int k = n / m;
    int l = n - m * k;
    int i, j, q;
    int real_sz = (l != 0) ? (k + 1) : k;
    int elements_in_blcok_line = k * m * m + l * m;
    double sum = 0;
    double max = 0;
    memset(tmp_line_n, 0, n * sizeof(double)); 
    for(i = 0; i < k; ++i) { // columns in inverse
        for(j = thread_number; j < k; j+=p) { // rows in matrix
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
        if (l!= 0 && (k % p == thread_number)) {
            j = k;
        //for(; j < real_sz; ++j) {
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
        for(j = thread_number; j < k; j+=p) {
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
        if (l!=0 && (k % p == thread_number)) {
            j = k;
        //for(; j < real_sz; ++j) {
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
    reduce_sum_fabs(p, tmp_line_n, n);
    max = tmp_line_n[0];
    for(i = 0; i < n; ++i) {
        if(tmp_line_n[i] > max) max = tmp_line_n[i];
    }
    return max;
}
#endif

#if 0
double calculate_discrepancy(Arg* a, double* matrix, double* inverse, double* tmp_block_m, double* tmp_line_n) {
    int n = a->n;
    int m = a->m;
    int p = a->p;
    int thread_number = a->thread_number;
    int k = n / m;
    int l = n - m * k;
    int i, j, q;
    int real_sz = (l != 0) ? (k + 1) : k;
    int elements_in_blcok_line = k * m * m + l * m;
    double sum = 0;
    double max = 0;
    memset(tmp_line_n, 0, n * sizeof(double)); 
    for(i = thread_number; i < k; i+=p) {
        for(j = 0; j < k; ++j) {
            memset(tmp_block_m, 0, m * m * sizeof(double));
            for(q = 0; q < k; ++q) {
                block_mult_add(matrix + i * elements_in_blcok_line + q * m * m
                        , m, m, inverse + q * elements_in_blcok_line + j * m * m, m, tmp_block_m);
            }
            block_mult_add(matrix + i * elements_in_blcok_line + q * m * m
                    , m, l, inverse + q * elements_in_blcok_line + j * l * m, m, tmp_block_m);
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
                tmp_line_n[j * m + w] += sum;
            }
        }
        for(;j < real_sz; ++j) {
            memset(tmp_block_m, 0, m * m * sizeof(double));
            for(q = 0; q < k; ++q) {
                block_mult_add(matrix + i * elements_in_blcok_line + q * m * m
                        , m, m, inverse + q * elements_in_blcok_line + j * m * m, l, tmp_block_m);
            }
            block_mult_add(matrix + i * elements_in_blcok_line + q * m * m
                    ,m, l, inverse + q * elements_in_blcok_line + j * l * m, l, tmp_block_m);
            for(int w = 0; w < l; ++w) {
                sum = 0;
                for(int v = 0; v < m; ++v) {
                    sum += fabs(tmp_block_m[v * m + w]);
                }
                tmp_line_n[j * m + w] += sum;
            }
        }
    }
    if (l!= 0 && (k % p == thread_number)) {
        i = k;
        for(j = 0; j < k; ++j) {
            memset(tmp_block_m, 0, m * m * sizeof(double));
            for(q = 0; q < k; ++q) {
                block_mult_add(matrix + i * elements_in_blcok_line + q * l * m
                        , l, m, inverse + q * elements_in_blcok_line + j * m * m, m, tmp_block_m);
            }
            block_mult_add(matrix + i * elements_in_blcok_line + q * l * m
                    , l, l, inverse + q * elements_in_blcok_line + j * l * m, m, tmp_block_m);
            for(int w = 0; w < m; ++w) {
                sum = 0;
                for(int v = 0; v < l; ++v) {
                    sum += fabs(tmp_block_m[v * m + w]);
                }
                tmp_line_n[j * m + w] += sum;
            }
        }
        memset(tmp_block_m, 0, m * m * sizeof(double));
        for(q = 0; q < k; ++q) {
            block_mult_add(matrix + i * elements_in_blcok_line + q * l * m
                    , l, m, inverse + q * elements_in_blcok_line + j * m * m, l, tmp_block_m);
        }
        block_mult_add(matrix + i * elements_in_blcok_line + q * l * m
                ,l, l, inverse + q * elements_in_blcok_line + j * l * m, l, tmp_block_m);
        if (j == i) {
            for(int s = 0; s < l; ++s) {
                tmp_block_m[s * l + s] -= 1;
            }
        }
        for(int w = 0; w < l; ++w) {
            sum = 0;
            for(int v = 0; v < l; ++v) {
                sum += fabs(tmp_block_m[v * m + w]);
            }
            tmp_line_n[j * m + w] += sum;
        }
    }
    reduce_sum_fabs(p, tmp_line_n, n);
    max = tmp_line_n[0];
    for(i = 0; i < n; ++i) {
        if(tmp_line_n[i] > max) max = tmp_line_n[i];
    }
    return max;
}
#endif

#if 0
double calculate_discrepancy(Arg* a, double* matrix, double* inverse, double* tmp_block_m, double* tmp_line_n) {
    int n = a->n;
    int m = a->m;
    int p = a->p;
    int thread_number = a->thread_number;
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
#endif
