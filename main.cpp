#include <ctime>
#include <new>
#include <stdio.h>
#include <sched.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#include <stdlib.h>

#include "read_print_fill.h"
#include "inverse.h"
#include "discrepancy.h"


int main(int argc, char* argv[]) {
    int n, m, r, s;
    int success1 = 0, success2 = 0;

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

    int nproc = get_nprocs();
    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    CPU_SET(nproc - 1, &cpu);
    sched_setaffinity(getpid(), sizeof(cpu), &cpu);

    double* matrix = new(std::nothrow) double[n * n];
    double* inverse = new(std::nothrow) double[n * n];
    int* permutations = new(std::nothrow) int[n / m];
    int* permutations_m = new(std::nothrow) int[m];  
    double* tmp_block_m = new(std::nothrow) double[m * m];
    double* tmp_inverse_m = new(std::nothrow) double[m * m];
    double* tmp_line_n = new(std::nothrow) double[n];

    if (matrix == nullptr || inverse == nullptr || permutations == nullptr
            || permutations_m == nullptr || tmp_block_m == nullptr
            || tmp_inverse_m == nullptr || tmp_line_n == nullptr) {
            delete[] matrix;
            delete[] inverse;
            delete[] permutations;
            delete[] permutations_m;
            delete[] tmp_block_m;
            delete[] tmp_inverse_m;
            delete[] tmp_line_n;
            return 1;
    }


    if (s != 0) {
        fill_matrix_with_formula(matrix, n, m, s);
    } else {
        success1 = read_matrix_from_file(matrix, n, m, argv[5]);
        if(success1 != 0) {
            delete[] matrix;
            delete[] inverse;
            delete[] permutations;
            delete[] permutations_m;
            delete[] tmp_block_m;
            delete[] tmp_inverse_m;
            delete[] tmp_line_n;
            printf("Error with file\n");
            return 1;
        }
    }

    create_id_matrix(inverse, n, m);

    double time1 = 0, time2 = 0;
    time1 = clock();
    success1 = inverse_matrix(matrix, inverse, n, m, permutations, tmp_block_m, tmp_inverse_m
            , permutations_m);
    time2 = clock();
    double t1 = (time2 - time1)/CLOCKS_PER_SEC;

    if (success1 != 0) {
        printf("This algorithm can't be applied\n");
    }
    if (s != 0) {
        fill_matrix_with_formula(matrix, n, m, s);
    } else {
        success2 = read_matrix_from_file(matrix, n, m, argv[5]);
    }

    double r1 = 0, r2 = 0;
    double t2 = 0;
    if (n <= 11000 && (success1 == 0) && (success2 == 0)) {
        time1 = clock();
        r1 = calculate_discrepancy(matrix, inverse, n, m, tmp_block_m, tmp_line_n);
        r2 = calculate_discrepancy(inverse, matrix, n, m, tmp_block_m, tmp_line_n);
        time2 = clock();
        t2 = (time2 - time1) / CLOCKS_PER_SEC;
    }

    print_matrix(inverse, n, m, r);
    printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n",
            argv[0], 19, r1, r2, t1, t2, s, n, m);
    delete[] matrix;
    delete[] inverse;
    delete[] permutations;
    delete[] permutations_m;
    delete[] tmp_block_m;
    delete[] tmp_inverse_m;
    delete[] tmp_line_n;
}
