#include "mpi.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
double calculate_discrepancy(double* matrix, double* inverse, int n, int m, double* tmp_block_m
        , double* tmp_line_n, int proc_num, int p, MPI_Comm comm) {
    int k = n / m;
    int l = n - m * k;
    int i, j, q;
    int real_sz = (l != 0) ? (k + 1) : k;
    int elements_in_blcok_line = k * m * m + l * m;
    double sum = 0;
    double max = 0;

    return 0;
}

