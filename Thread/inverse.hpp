#pragma once

#include "thread.hpp"

int inverse_block(double* block, double* inverse, int* permutations, int m, double norm);
int inverse_matrix(Arg* a, double* matrix, double* inverse, int n, int m, int* permutations
        , double* block_m, double* inv_block_m, int* permutations_m);

