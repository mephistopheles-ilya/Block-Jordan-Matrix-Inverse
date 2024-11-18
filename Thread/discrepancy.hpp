#pragma once

#include "thread.hpp"

double calculate_discrepancy(Arg* a, double* matrix, double* inverse, int n, int m, double* tmp_block_m
        , double* tmp_line_n);

