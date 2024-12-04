#include <stdio.h>
#include <stdlib.h>
#include <new>
#include <ctime>
#include <cmath>

#include "func.hpp"
#include "eigenvalues.hpp"

int main(int argc, char* argv[]) {
    int n = 0, m = 0, k = 0;
    double eps = 0;
    double trace1 = 0, length1 = 0, norm = 0;
    double t1 = 0, t2 = 0, res1 = 0, res2 = 0;
    int its = 0;
    if (argc != 5 && argc != 6) {
        printf("Usage <n> <m> <eps> <k> <filename>\n");
        return 1;
    }
    if (sscanf(argv[1], "%d", &n) != 1 || sscanf(argv[2], "%d", &m) != 1
            || sscanf(argv[3], "%lf", &eps) != 1 || sscanf(argv[4], "%d", &k) != 1) {
        printf("Usage <n> <m> <eps> <k> <filename>\n");
        return 1;
    }
    if (n <= 0 || m < 0 || k < 0 || k > 4 || eps < 0) {
        printf("Wrong type of arguments\n");
        return 1;
    }

    double* matrix = new (std::nothrow) double[n * n];
    double* x = new (std::nothrow) double[n];
    if (matrix == nullptr || x == nullptr) {
        printf("Not enough memmory\n");
        delete[] matrix;
        delete[] x;
        return 2;
    }
    io_status status = io_status::undef;
    if (k == 0) {
        status = read_matrix_from_file(matrix, n, argv[5]);
    } else {
        status = fill_matrix_with_formula(matrix, n, k);
    }

    if (status == io_status::cannot_open_file) {
        printf("Cannot open file\n");
        delete[] matrix;
        delete[] x;
        return 3;
    }
    if (status == io_status::error_read) {
        printf("Cannot read elements from file\n");
        delete[] matrix;
        delete[] x;
        return 4;
    }
    if (status == io_status::not_enough_elements) {
        printf("Not enough elements in file\n");
        delete[] matrix;
        delete[] x;
        return 5;
    }
    norm = matrix_norm(matrix, n);
    trace1 = matrix_trace(matrix, n);
    length1 = matrix_length(matrix, n);

    printf("Initial Matrix:\n");
    print_matrix(matrix, n, n, m);

    t1 = clock();
    matrix_to_almost_triangle(matrix, x, n, eps, norm);
    t1 = (clock() - t1)/CLOCKS_PER_SEC;

    double trace2 = matrix_trace(matrix, n);
    printf("DEBUG TRACE_DIFF: %e %e \n", std::fabs(trace1 - trace2), std::fabs(trace1 - trace2)/norm);

    printf ("%s : Residual1 = %e Residual2 = %e Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n",
            argv[0], res1, res2, its, its / n, t1, t2);
    print_matrix(matrix, n, n, m);


    delete[] matrix;
    delete[] x;
    return 0;
}


