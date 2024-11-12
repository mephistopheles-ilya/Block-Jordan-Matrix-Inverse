#pragma once

int inverse_block(double* block, double* inverse, int* permutations, int m, double norm);
int inverse_matrix(double* matrix, double* inverse, int n, int m, int* permutations
        , double* block_m, double* inv_block_m, int* permutations_m);

