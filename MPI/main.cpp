#include "mpi.h"
#include <stdio.h>
#include <fenv.h>

int main(int argc, char* argv[]) {
    int p, k;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);

    //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &k);
    if (k == 0) {
        printf("p = %d\n", p);
        double a = 1.0/0;
    }
    MPI_Finalize();
    return 0;
}

