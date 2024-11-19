#pragma once

#include "thread.hpp"

void block_mult(double* A, int av, int ag, double* B, int bg, double* C); 
void block_mult_sub(double* A, int av, int ag, double* B, int bg, double* C); 
void block_mult_add(double* A, int av, int ag, double* B, int bg, double* C); 
double norm_block(double* block, int m); 
double norm_matrix(double* matrix, int n, int m);
void swap_block_col_th(Arg* a, int i, int j, double* tmp_block_m);
void swap_block_lin(double* matrix, int n, int m, int i, int j, double* tmp_block_m);
void fill_id_block(double* block, int m);
