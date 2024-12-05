#include <stdio.h>
#include <stdlib.h>
#include <new>
#include <ctime>
#include <cmath>
#include <string.h>

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
    double* x1 = new (std::nothrow) double[n];
    double* x2 = new (std::nothrow) double[n];
    double* eigenvalues = new (std::nothrow) double[n];
    if (matrix == nullptr || x1 == nullptr || x2 == nullptr || eigenvalues == nullptr) {
        printf("Not enough memmory\n");
        delete[] matrix;
        delete[] x1;
        delete[] x2;
        delete[] eigenvalues;
        return 2;
    }

    memset(matrix, 0, n * n * sizeof(double));
    memset(x1, 0, n * sizeof(double));
    memset(x2, 0, n * sizeof(double));
    memset(eigenvalues, 0, n * sizeof(double));

    io_status status = io_status::undef;
    if (k == 0) {
        status = read_matrix_from_file(matrix, n, argv[5]);
    } else {
        status = fill_matrix_with_formula(matrix, n, k);
    }

    if (status == io_status::cannot_open_file) {
        printf("Cannot open file\n");
        delete[] matrix;
        delete[] x1;
        delete[] x2;
        delete[] eigenvalues;
        return 3;
    }
    if (status == io_status::error_read) {
        printf("Cannot read elements from file\n");
        delete[] matrix;
        delete[] x1;
        delete[] x2;
        delete[] eigenvalues;
        return 4;
    }
    if (status == io_status::not_enough_elements) {
        printf("Not enough elements in file\n");
        delete[] matrix;
        delete[] x1;
        delete[] x2;
        delete[] eigenvalues;
        return 5;
    }
    norm = matrix_norm(matrix, n);
    trace1 = matrix_trace(matrix, n);
    length1 = matrix_length(matrix, n);

    printf("Initial Matrix:\n");
    print_matrix(matrix, n, n, m, n);

    t1 = clock();
    matrix_to_almost_triangle(matrix, x1, n, eps, norm);
    t1 = (clock() - t1)/CLOCKS_PER_SEC;

    int error_code = 0;
    t2 = clock();
    error_code = find_eigenvalues(matrix, x1, x2, eigenvalues, n, eps, norm, its);
    t2 = (clock() - t2)/CLOCKS_PER_SEC; 

    if (error_code == 0) {
        double trace2 = 0;
        for(int i = 0; i < n; ++i) {
            trace2 += eigenvalues[i];
        }
        double length2 = 0;
        for(int i = 0; i < n; ++i) {
            length2 += eigenvalues[i] * eigenvalues[i];
        }
        length2 = std::sqrt(length2);

        if (norm > 0) {
            res1 = std::fabs(trace1 - trace2) / norm;
            res2 = std::fabs(length1 - length2) / norm;
        } else {
            res1 = 0;
            res2 = 0;
        }

        printf("Eigenvalues :\n");
        print_matrix(eigenvalues, 1, n, m, n);
    } else {
        res1 = -1;
        res2 = -1;
    }

    printf ("%s : Residual1 = %e Residual2 = %e Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n",
            argv[0], res1, res2, its, its / n, t1, t2);

    delete[] matrix;
    delete[] x1;
    delete[] x2;
    delete[] eigenvalues;
    return 0;
}



