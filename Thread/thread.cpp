#include "thread.hpp"
#include "reduce_sum.hpp"
#include "read_print_fill.hpp"
#include "matrix.hpp"

#include <sys/sysinfo.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


void stick_matrices_to_core(Arg* a) {
    int n = a->n;
    int m = a->m;
    int p = a->p;
    int thread_number = a->thread_number;
    double* matrix = a->matrix;
    double* inverse = a->inverse;
    int k = n / m;
    int l = n - k * m;
    for(int i = thread_number; i < k; i += p) {
        memset(matrix + i * n * m, 0, n * m * sizeof(double));
        memset(inverse + i * n * m, 0, n * m * sizeof(double));
    }
    if (l != 0 && ((k + 1) % p == thread_number)) {
        memset(matrix + k * n * m, 0, n * l * sizeof(double));
        memset(inverse + k * n * m, 0, n * l * sizeof(double));
    }
    pthread_barrier_wait(a->barrier);
}

void fill_matrices(Arg* a) {
    int n = a->n;
    int m = a->m;
    int s = a->s;
    int p = a->p;
    int thread_number = a->thread_number;
    char* filename = a->filename;
    double* matrix = a->matrix;
    double* inverse = a->inverse;

    int error = 0;

    if (thread_number == 0) {
        if (filename == nullptr) {
            error = read_matrix_from_file(matrix, n, m, filename);
        } else {
            fill_matrix_with_formula(matrix, n, m, s);
        }
        create_id_matrix(inverse, n, m);
        a->error = error;
    }
    reduce_sum(p, &(a->error), 1);
}


void* thread_func(void* arg) {
    Arg* a = (Arg*)arg;
    int n = a->n;
    int m = a->m;
    int p = a->p;
    int r = a->r;
    int s = a->s;
    int thread_number = a->thread_number;
    char* filename = a->filename;
    double* matrix = a->matrix;
    double* inverse = a->inverse;

    double norm = 0;
    double local_time = 0;
    double global_time = 0;

    int nproc = get_nprocs();
    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    CPU_SET(nproc - 1 - (thread_number % (nproc)), &cpu);
    pthread_setaffinity_np(pthread_self(), sizeof(cpu), &cpu);

    stick_matrices_to_core(a);
    fill_matrices(a);
    if (a->error > 0) {
        if (thread_number == 0) {
            printf("Bad file\n");
            //TODO free extra arrays
        }
        return nullptr;
    }

    if (thread_number == 0) {
        norm = norm_matrix(matrix, n, m);
    }
    pthread_barrier_wait(a->barrier);


    return nullptr;
}

