#include "mpi.h"
#include <stdio.h>
#include <new>

#include <fenv.h>

int main(int argc, char* argv[]) {
    int n = 0, m = 0, r = 0, s = 0, p = 0, proc_num = 0;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);

    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &proc_num);

    if (argc != 5 && argc != 6) {
        if (proc_num == 0) {
            printf("Wrong amount of arguments\n");
        }
        MPI_Finalize();
        return 0;
    }
    if (!(sscanf(argv[1], "%d", &n) == 1 
        && sscanf(argv[2], "%d", &m) == 1
        && sscanf(argv[3], "%d", &r) == 1
        && sscanf(argv[4], "%d", &s) == 1)) {
        if (proc_num == 0) {
            printf("Wrong type of arguments\n");
        }
        MPI_Finalize();
        return 0;
    } 
    if (n <= 0 || m <= 0 || r < 0 || s < 0 || s > 4 || (s==0 && argc==5) || n < m || p <= 0) {
        if (proc_num == 0) {
            printf("Wrong arguments\n");
        }
        MPI_Finalize();
        return 0;
    }

    if (proc_num == 0) {
        printf("Usage %s %d %d %d %d %s\n", argv[0], n, m, r, s, argv[5]);
        printf("Amount of processes = %d\n", p);
    }
    
    //but actually it is columns
    int num_of_rows = (n % m == 0) ? n/m : (n/m + 1); 
    int local_num_rows = (num_of_rows % p == 0) ? num_of_rows/p : (num_of_rows/p + 1);

    double* matrix = new(std::nothrow) double[local_num_rows * m * n];
    double* inverse = new(std::nothrow) double[local_num_rows * m * n];

    if (matrix == nullptr || inverse == nullptr) {
        if (proc_num == 0) {
            printf("Not enough memmory\n");
        }
        delete[] matrix;
        delete[] inverse;
    }

    MPI_Finalize();
    return 0;
}

