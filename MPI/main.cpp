#include "mpi.h"
#include <stdio.h>
#include <new>
#include <string.h>


#include <fenv.h>
#include "read_print_fill.hpp"

#define LEN 256

int main(int argc, char* argv[]) {
    int n = 0, m = 0, r = 0, s = 0, p = 0, proc_num = 0;
    int error_loc = 0, error_glob = 0;
    char buf[256];
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Status status;
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
    int rows = get_rows(n, m, p, proc_num);

    double* matrix = new(std::nothrow) double[max_rows * m * n];
    double* inverse = new(std::nothrow) double[max_rows * m * n];
    double* tmp_row_matrix = new(std::nothrow) double[m * n];
    double* tmp_row_inverse = new(std::nothrow) double[m * n];

    if (matrix == nullptr || inverse == nullptr || tmp_row_inverse == nullptr || tmp_row_matrix == nullptr) {
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
        MPI_Finalize();
        return 0;
    }
    memset(matrix, 0, max_rows * m * n);
    memset(inverse, 0, max_rows * m * n);
    memset(tmp_row_matrix, 0, n * m);
    memset(tmp_row_inverse, 0, n * m);

    snprintf(buf, LEN, "max_rows = %d, rows = %d, proc_num = %d\n", max_rows, rows, proc_num);
    if(proc_num != 0) {
        MPI_Send(buf, strlen(buf) + 1, MPI_CHAR, 0, 0, comm);
    } else {
        printf("%s", buf);
        for(int i = 1; i < p; ++i) {
            MPI_Recv(buf, LEN, MPI_CHAR, i, 0, comm, &status);
            printf("%s", buf);
        }
    }

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
        MPI_Finalize();
        return 0;
    }
    if (proc_num == 0) {
        printf("Initial Matrix :\n");
    }
    print_matrix(matrix, n, m, proc_num, p, r, tmp_row_matrix, comm);


    delete[] matrix;
    delete[] inverse;
    delete[] tmp_row_inverse;
    delete[] tmp_row_matrix;
    MPI_Finalize();
    return 0;
}

