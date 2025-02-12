#pragma once

#include "mpi.h"

int get_max_rows(int n, int m, int p);
int get_rows(int n, int m, int p, int proc_num);
int fill_matrix(double* matrix, int n, int m, int s, char* file_name, int proc_num, int p
        , double* buf, MPI_Comm comm );
int ltg(int m, int p, int k, int i_loc);
int fill_id_matrix(double* inverse, int n, int m, int proc_num, int p);
int print_matrix(double* matrix, int n, int m, int proc_num, int p, int r, double* buf, MPI_Comm comm);
