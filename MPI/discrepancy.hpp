#pragma once

#include "mpi.h"
double calculate_discrepancy(double* matrix, double* inverse, int n, int m, double* tmp_block_m
        , double* tmp_line_n, int proc_num, int p, MPI_Comm comm);

