#pragma once

enum class io_status {
    undef,
    error_read,
    not_enough_elements,
    cannot_open_file,
    succes
};

io_status read_matrix_from_file(double* matrix, int n, char* filename);
io_status fill_matrix_with_formula(double* matrix, int n, int k);
double f(int k, int n, int i, int j);
double matrix_norm(double* matrix, int n);
double matrix_trace(double* matrix, int n);
double matrix_length(double* matrix, int n);
void print_matrix(double* matrix, int a, int b, int m);

