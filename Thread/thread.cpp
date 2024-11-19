#include "thread.hpp"
#include "reduce_sum.hpp"
#include "read_print_fill.hpp"
#include "matrix.hpp"
#include "inverse.hpp"
#include "discrepancy.hpp"

#include <sys/sysinfo.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <new>

double get_cpu_time() {
    struct rusage buf;
    getrusage(RUSAGE_THREAD, &buf);
    return (double)(buf.ru_utime.tv_sec) + (double)(buf.ru_utime.tv_usec)/1000000.;
}
double get_full_time() {
    struct timeval buf;
    gettimeofday(&buf, 0);
    return (double)(buf.tv_sec) + (double)(buf.tv_usec)/1000000.;
}


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
    if (l != 0 && (k % p == thread_number)) {
        memset(matrix + k * n * m, 0, n * l * sizeof(double));
        memset(inverse + k * n * m, 0, n * l * sizeof(double));
    }
    pthread_barrier_wait(a->barrier);
}

void fill_matrices(Arg* a, bool id) {
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
        if (s == 0) {
            error = read_matrix_from_file(matrix, n, m, filename);
        } else {
            fill_matrix_with_formula(matrix, n, m, s);
        }
        if(id == true) {
            create_id_matrix(inverse, n, m);
        }
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

    double matrix_norm = 0;
    double local_time = 0;
    double global_time = 0;
    double r1 = 0, r2 = 0;

    int nproc = get_nprocs();
    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    CPU_SET(nproc - 1 - (thread_number % (nproc)), &cpu);
    pthread_setaffinity_np(pthread_self(), sizeof(cpu), &cpu);

    int* permutations = new(std::nothrow) int[n / m];
    int* permutations_m = new(std::nothrow) int[m];  
    double* block_m = new(std::nothrow) double[m * m];
    double* inv_block_m = new(std::nothrow) double[m * m];
    double* line_n = new(std::nothrow) double[n];

    if (permutations == nullptr || permutations_m == nullptr 
        || block_m == nullptr || inv_block_m == nullptr || line_n == nullptr) {
        delete[] permutations;
        delete[] permutations_m;
        delete[] block_m;
        delete[] inv_block_m;
        delete[] line_n;
        if (thread_number == 0) {
            printf("Not enough memmory in thread number %d\n", thread_number);
        }
        return nullptr;
    }

    memset(permutations, 0, (n/m) * sizeof(int));
    memset(permutations_m, 0, m * sizeof(int));
    memset(block_m, 0, m * m * sizeof(double));
    memset(inv_block_m, 0, m * m * sizeof(double));
    memset(line_n, 0, n * sizeof(double));

    for(int i = 0; i < (n/m); ++i) {
        permutations[i] = i;
    }

    stick_matrices_to_core(a);
    fill_matrices(a, true);
    if (a->error > 0) {
        if (thread_number == 0) {
            printf("Bad file\n");
        }
        delete[] permutations;
        delete[] permutations_m;
        delete[] block_m;
        delete[] inv_block_m;
        delete[] line_n;
        return nullptr;
    }

    if (thread_number == 0) {
        printf("Initial matrix :\n");
        print_matrix(matrix, n, m, r);
        matrix_norm = norm_matrix(matrix, n, m);
    } else {
        matrix_norm = 0;
    }
    reduce_sum(p, &matrix_norm, 1);
    global_time = get_full_time();
    local_time = get_cpu_time();
    inverse_matrix(a, permutations, block_m, inv_block_m, permutations_m, matrix_norm);
    local_time = get_cpu_time() - local_time;
    reduce_sum(p, &(a->error), 1);
    global_time = get_full_time() - global_time;
    a -> local_time_inv = local_time;
    a -> global_time_inv = global_time;
    if (a-> error > 0) {
        if (thread_number == 0) {
            printf("This algorithm can't be applied\n");
        }
        a -> r1 = -1;
        a -> r2 = -1;
        delete[] permutations;
        delete[] permutations_m;
        delete[] block_m;
        delete[] inv_block_m;
        delete[] line_n;
        return nullptr;
    }
    fill_matrices(a, false);
    if (a->error > 0) {
        if (thread_number == 0) {
            printf("Bad file\n");
        }
        delete[] permutations;
        delete[] permutations_m;
        delete[] block_m;
        delete[] inv_block_m;
        delete[] line_n;
        return nullptr;
    }

    if (n < 11000) {
        global_time = get_full_time();
        local_time = get_cpu_time();
        if (thread_number == 0) {
            r1 = calculate_discrepancy(a, matrix, inverse, n, m, block_m, line_n);  
        }
        a -> local_time_dis = get_cpu_time() - local_time;
        pthread_barrier_wait(a->barrier);
        local_time = get_cpu_time();
        if (thread_number == 0) {
            r2 = calculate_discrepancy(a, inverse, matrix, n, m, block_m, line_n);
        }
        a -> local_time_dis += get_cpu_time() - local_time;
        pthread_barrier_wait(a->barrier);
        global_time = get_full_time() - global_time;
        a -> global_time_dis = global_time;
        a -> r1 = r1;
        a -> r2 = r2;
    }
    if (thread_number == 0) {
        printf("Inverse matrix :\n");
        print_matrix(inverse, n, m, r);
    }

    delete[] permutations;
    delete[] permutations_m;
    delete[] block_m;
    delete[] inv_block_m;
    delete[] line_n;
    return nullptr;
}

