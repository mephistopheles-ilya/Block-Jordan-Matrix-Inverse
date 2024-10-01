#pragma once

void fill_matrix_with_formula(double* matrix, int n, int m, int s);
int read_matrix_from_file(double* matrix, int n, int m, char* filename);
void print_matrix(double* matrix, int n, int m, int r);
void create_id_matrix(double* matrix, int n, int m);
double* get_block(double* matrix, int n, int m, int i, int j);
double calculate_discrepancy(double*, double*, int, int);
