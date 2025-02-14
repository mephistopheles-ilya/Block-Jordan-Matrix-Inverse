#pragma once

#include "mpi.h"

void block_mult(double* A, int av, int ag, double* B, int bg, double* C); 
void block_mult_sub(double* A, int av, int ag, double* B, int bg, double* C); 
void block_mult_add(double* A, int av, int ag, double* B, int bg, double* C); 
double norm_block(double* block, int m); 
double norm_matrix(double* matrix, int n, int m, int p, int proc_num, MPI_Comm comm);
void fill_id_block(double* block, int m);
void swap_block_col(double* matrix, int n, int m, int i, int j, int proc_num, int p, double* tmp_block);
void swap_block_line(double* matrix, int n, int m, int i_glob, int j_glob, int proc_num, int p, double* buf
        , MPI_Comm comm);
