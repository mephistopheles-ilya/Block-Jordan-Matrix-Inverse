#pragma once

void block_mult(double* A, int av, int ag, double* B, int bg, double* C); 
double calculate_discrepancy(double* matrix, double* inverse, int n, int m, double* tmp_block_m
        , double* tmp_line_n);
void block_mult_sub(double* A, int av, int ag, double* B, int bg, double* C); 
