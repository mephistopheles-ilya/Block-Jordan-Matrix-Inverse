#include <ctime>
#include <new>
#include <stdio.h>
#include <sched.h>
#include <sys/sysinfo.h>
#include <unistd.h>

#include "functions.h"

int main(int argc, char* argv[]) {
    int n, m, r, s;

//check command line arguments
    if (argc != 5 && argc != 6) {
        printf("Wrong amount of arguments\n");
        return 1;
    }

    if (!(sscanf(argv[1], "%d", &n) == 1 
        && sscanf(argv[2], "%d", &m) == 1
        && sscanf(argv[3], "%d", &r) == 1 
        && sscanf(argv[4], "%d", &s) == 1)) {
        printf("Wrong type of arguments\n");
        return 1;
    } 
    if (n <= 0 || m <= 0 || r < 0 || s < 0 || s > 4 || (s==0 && argc==5) || n < m) {
        printf("Wrong arguments\n");
        return 1;
    }

//stick this thread to core
    int nproc = get_nprocs();
    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    CPU_SET(nproc - 1, &cpu);
    sched_setaffinity(getpid(), sizeof(cpu), &cpu);

//allocate matrices
    double* matrix = new(std::nothrow) double[n * n];
    double* inverse = new(std::nothrow) double[n * n];
    if (matrix == nullptr || inverse == nullptr) {
        delete[] matrix;
        delete[] inverse;
        return 1;
    }

//read matrix 
    if (s != 0) {
        fill_matrix_with_formula(matrix, n, m, s);
    } else {
        int success = 0;
        success = read_matrix_from_file(matrix, n, m, argv[5]);
        if(success != 0) {
            delete[] matrix;
            delete[] inverse;
            printf("Error with file\n");
            return 1;
        }
    }


    double time1 = 0, time2 = 0;
//find inverse matrix
    time1 = clock();
    //TODO
    create_id_matrix(inverse, n, m);
    time2 = clock();
    double t1 = (time2 - time1)/CLOCKS_PER_SEC;

//compute discrepancy
    double r1 = 0, r2 = 0;
    double t2 = 0;
    if (n <= 11000) {
        time1 = clock();
        //TODO
        time2 = clock();
        t2 = (time2 - time1) / CLOCKS_PER_SEC;
    }

//print result
    print_matrix(inverse, n, m, r);
    printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n",
            argv[0], 19, r1, r2, t1, t2, s, n, m);
    delete[] matrix;
    delete[] inverse;
}

