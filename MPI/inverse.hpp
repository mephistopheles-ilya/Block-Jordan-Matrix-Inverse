#pragma once

#include "mpi.h"

int inverse_block(double* block, double* inverse, int* permutations, int m, double norm);
int inverse_matrix(double* matrix, double* inverse, int n, int m, int* permutations
        , double* block_m, double* inv_block_m, int* permutations_m, double matrix_norm
        , double* matrix_buf, double* inverse_buf, int proc_num, int p, MPI_Comm com);

struct double_int {
    double x = 0;
    int y = 0;
};

