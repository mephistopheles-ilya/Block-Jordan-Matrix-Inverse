#include <ctime>
#include <new>
#include <stdio.h>
#include <sched.h>
#include <sys/sysinfo.h>
#include <unistd.h>

#include "../functions.h"
#include "../matrix.h"

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

    block_mult(matrix, n, n, matrix, n, inverse);

    print_matrix(inverse, n, m, r);
    delete[] matrix;
    delete[] inverse;
}

