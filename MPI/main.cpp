#include "mpi.h"
#include "read_print_fill.hpp"
#include "matrix.hpp"
#include "inverse.hpp"
#include "discrepancy.hpp"

#include <stdio.h>
#include <new>
#include <string.h>


#include <fenv.h>


int main(int argc, char* argv[]) {
    int n = 0, m = 0, r = 0, s = 0, p = 0, proc_num = 0;
    int task = 19;
    double r1 = -1, r2 = -1, t1 = 0, t2 = 0;
    int error_loc = 0, error_glob = 0;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);

    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &proc_num);

    if (argc != 5 && argc != 6) {
        error_loc = 1;
        if (proc_num == 0) {
            printf("Wrong amount of arguments\n");
        }
    } else if (!(sscanf(argv[1], "%d", &n) == 1 
        && sscanf(argv[2], "%d", &m) == 1
        && sscanf(argv[3], "%d", &r) == 1
        && sscanf(argv[4], "%d", &s) == 1)) {
        error_loc = 1;
        if (proc_num == 0) {
            printf("Wrong type of arguments\n");
        }
    } else if (n <= 0 || m <= 0 || r < 0 || s < 0 || s > 4 || (s==0 && argc==5) || n < m || p <= 0) {
        error_loc = 1;
        if (proc_num == 0) {
            printf("Wrong arguments\n");
        }
    }

    MPI_Allreduce(&error_loc, &error_glob, 1, MPI_INT, MPI_MAX, comm);
    if(error_glob != 0) {
        if (proc_num == 0) {
            printf("Errors with arguments\n");
        }
        MPI_Finalize();
        return 0;
    }

    if (proc_num == 0) {
        printf("Usage %s %d %d %d %d %s\n", argv[0], n, m, r, s, argv[5]);
        printf("Amount of processes = %d\n", p);
    }
    
    int max_rows = get_max_rows(n, m, p);

    double* matrix = new(std::nothrow) double[max_rows * m * n];
    double* inverse = new(std::nothrow) double[max_rows * m * n];
    double* tmp_row_matrix = new(std::nothrow) double[m * n];
    double* tmp_row_inverse = new(std::nothrow) double[m * n];
    int* permutations = new(std::nothrow) int[n/m];
    int* permutations_m = new(std::nothrow) int[m];
    double* block_m = new(std::nothrow) double[m * m];
    double* inv_block = new(std::nothrow) double[m * m];

    if (matrix == nullptr || inverse == nullptr || tmp_row_inverse == nullptr || tmp_row_matrix == nullptr
            || permutations == nullptr || permutations_m == nullptr || block_m == nullptr || inv_block == nullptr) {
        error_loc = 1;
    }

    MPI_Allreduce(&error_loc, &error_glob, 1, MPI_INT, MPI_MAX, comm);
    if (error_glob != 0) {
        if (proc_num == 0) {
            printf("Not enough memmory\n");
        }
        delete[] matrix;
        delete[] inverse;
        delete[] tmp_row_inverse;
        delete[] tmp_row_matrix;
        delete[] permutations;
        delete[] permutations_m;
        delete[] block_m;
        delete[] inv_block;
        MPI_Finalize();
        return 0;
    }
    memset(matrix, 0, max_rows * m * n * sizeof(double));
    memset(inverse, 0, max_rows * m * n * sizeof(double));
    memset(tmp_row_matrix, 0, n * m * sizeof(double));
    memset(tmp_row_inverse, 0, n * m * sizeof(double));
    memset(permutations, 0, n/m * sizeof(int));
    memset(permutations_m, 0, m * sizeof(int));
    memset(block_m, 0, m * m * sizeof(double));
    memset(inv_block, 0, m * m * sizeof(double));

    fill_id_matrix(inverse, n, m, proc_num, p);
    error_loc = fill_matrix(matrix, n, m, s, argv[5], proc_num, p, tmp_row_matrix, comm);
    if (error_loc != 0) {
        if (proc_num == 0) {
            if (error_loc == 1) {
                printf("Cannot open file\n");
            }
            if (error_loc == 2) {
                printf("Not enough elements in file or cannot read element\n");
            }
        }
        delete[] matrix;
        delete[] inverse;
        delete[] tmp_row_inverse;
        delete[] tmp_row_matrix;
        delete[] permutations;
        delete[] permutations_m;
        delete[] block_m;
        delete[] inv_block;
        MPI_Finalize();
        return 0;
    }
    if (proc_num == 0) {
        printf("Initial Matrix :\n");
    }
    print_matrix(matrix, n, m, proc_num, p, r, tmp_row_matrix, comm);

    double matrix_norm = norm_matrix(matrix, n, m, p, proc_num, comm);
    MPI_Barrier(comm);
    t1 = get_full_time();
    error_loc = inverse_matrix(matrix, inverse, n, m, permutations, block_m, inv_block
            , permutations_m, matrix_norm, tmp_row_matrix, tmp_row_inverse, proc_num, p, comm); 
    t1 = get_full_time() - t1;
    MPI_Allreduce(&error_loc, &error_glob, 1, MPI_INT, MPI_MAX, comm);
    if (error_glob != 0) {
        if (proc_num == 0) {
            printf("This algorithm can't be applied\n");
        }
    }

    if (error_glob == 0 && n <= 11000) {
        fill_matrix(matrix, n, m, s, argv[5], proc_num, p, tmp_row_matrix, comm);

        MPI_Barrier(comm);
        t2 = get_full_time();
        r1 = calculate_discrepancy(matrix, inverse, n, m, tmp_row_matrix, tmp_row_inverse, proc_num, p, comm);
        r2 = calculate_discrepancy(inverse, matrix, n, m, tmp_row_matrix, tmp_row_inverse, proc_num, p, comm);
        t2 = get_full_time() - t2;
    }

    if(error_glob == 0) {
        if (proc_num == 0) {
            printf("Inverse matrix :\n");
        }
        print_matrix(inverse, n, m, proc_num, p, r, tmp_row_matrix, comm);
    }

    if(proc_num == 0) {
        printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
                argv[0], task, r1, r2, t1, t2, s, n, m, p);
    }


    delete[] matrix;
    delete[] inverse;
    delete[] tmp_row_inverse;
    delete[] tmp_row_matrix;
    delete[] permutations;
    delete[] permutations_m;
    delete[] block_m;
    delete[] inv_block;
    MPI_Finalize();
    return 0;
}

