#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <cmath>

#include "func.hpp"


double f(int k, int n, int i, int j) {
    ++i;
    ++j;
    if (k == 1) {
        return n - std::max(i, j) + 1;
    }
    if (k == 2) {
        if (i == j) {
            return 2;
        }
        if (std::fabs(i - j) == 1) {
            return -1;
        }
        return 0;
    }
    if (k == 3) {
        if (i == j && i < n) {
            return 1;
        }
        if (j == n) {
            return i;
        }
        if (i == n) {
            return j;
        }
        return 0;
    }
    if (k == 4) {
        return 1./(i + j - 1);
    }
    return 0;
}


io_status read_matrix_from_file(double* matrix, int n, char* filename) {
    FILE* f = fopen(filename, "r");
    if (f == nullptr) {
        return io_status::cannot_open_file;
    }

    int ret_val = 0;
    int i = 0;
    for(i = 0; i < n * n; ++i) {
        if((ret_val = fscanf(f, "%lf", matrix + i)) != 1) {
            break;
        }
    }
    fclose(f);
    if (i == n * n) {
        return io_status::succes;
    }
    if (ret_val == 0) {
        return io_status::error_read;
    }
    if (ret_val == EOF) {
        return io_status::not_enough_elements;
    }
    return io_status::succes;
}

io_status fill_matrix_with_formula(double* matrix, int n, int k) {
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            matrix[i * n + j] = f(k, n, i, j);
        }
    }
    return io_status::succes;
}

double matrix_norm(double* matrix, int n) {
    double norm = 0;
    double sum = 0;
    for(int i = 0; i < n; ++i) {
        sum = 0;
        for(int j = 0; j < n; ++j) {
            sum += std::fabs(matrix[j * n + i]);
        }
        if (sum > norm) {
            norm = sum;
        }
    }
    return norm;
}

double matrix_trace(double* matrix, int n) {
    double trace = 0;
    for(int i = 0; i < n; ++i) {
        trace += matrix[i * n + i];
    }
    return trace;
}

double matrix_length(double* matrix, int n) {
    double length = 0;
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            length += matrix[i * n + j] * matrix[j * n + i];
        }
    }
    length = (length < 0) ? 0 : length;
    return std::sqrt(length);
}

void print_matrix(double* matrix, int a, int b, int m, int n) {
    a = (a > m) ? m : a;
    b = (b > m) ? m : b;
    for(int i = 0; i < a; ++i) {
        for(int j = 0; j < b; ++j) {
            printf(" %10.3e", matrix[i * n + j]);
        }
        printf("\n");
    }
}




