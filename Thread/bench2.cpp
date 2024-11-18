#include <iostream>
#include <pthread.h>
#include <sched.h>
#include <sys/sysinfo.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <mutex>
#include <ctime>
#include <stdio.h>
#include <sched.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#include <stdlib.h>


std::mutex mut;

struct Args {
    int n = 0;
    int m = 0;
    int p = 0;
    int k = 0;
    double* matrix1 = nullptr; 
    double* matrix2 = nullptr;
    pthread_t tid = 0;
    pthread_barrier_t* barier = nullptr;
    double time = 0;
};

void block_mult(double* A, int av, int ag, double* B, int bg, double* C) {
    int av_reminder = av % 3;
    int bg_reminder = bg % 3;
    int _av = av - av_reminder;
    int _bg = bg - bg_reminder;
    int i = 0, j = 0, q = 0;
    double c00, c01, c02, c10, c11, c12, c20, c21, c22;
    for(i = 0; i < _av; i += 3) {
        for(j = 0; j < _bg; j+= 3) {
            c00 = 0; c01 = 0; c02 = 0;
            c10 = 0; c11 = 0; c12 = 0;
            c20 = 0; c21 = 0; c22 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
                c01 += A[ag * (i) + q] * B[bg * q + (j + 1)];
                c02 += A[ag * (i) + q] * B[bg * q + (j + 2)];
                c10 += A[ag * (i + 1) + q] * B[bg * q + (j)];
                c11 += A[ag * (i + 1) + q] * B[bg * q + (j + 1)];
                c12 += A[ag * (i + 1) + q] * B[bg * q + (j + 2)];
                c20 += A[ag * (i + 2) + q] * B[bg * q + (j)];
                c21 += A[ag * (i + 2) + q] * B[bg * q + (j + 1)];
                c22 += A[ag * (i + 2) + q] * B[bg * q + (j + 2)];
            }
            C[bg * (i) + (j)] = c00;
            C[bg * (i) + (j + 1)] = c01;
            C[bg * (i) + (j + 2)] = c02;
            C[bg * (i + 1) + (j)] = c10;
            C[bg * (i + 1) + (j + 1)] = c11;
            C[bg * (i + 1) + (j + 2)] = c12;
            C[bg * (i + 2) + (j)] = c20;
            C[bg * (i + 2) + (j + 1)] = c21;
            C[bg * (i + 2) + (j + 2)] = c22;
        }
        for(;j < bg; ++j) {
            c00 = 0;
            c10 = 0;
            c20 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
                c10 += A[ag * (i + 1) + q] * B[bg * q + (j)];
                c20 += A[ag * (i + 2) + q] * B[bg * q + (j)];
            }
            C[bg * (i) + (j)] = c00;
            C[bg * (i + 1) + (j)] = c10;
            C[bg * (i + 2) + (j)] = c20;
        }
    }
    for(; i < av; ++i) {
        for(j = 0; j < _bg; j+= 3) {
            c00 = 0; c01 = 0; c02 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
                c01 += A[ag * (i) + q] * B[bg * q + (j + 1)];
                c02 += A[ag * (i) + q] * B[bg * q + (j + 2)];
            }
            C[bg * (i) + (j)] = c00;
            C[bg * (i) + (j + 1)] = c01;
            C[bg * (i) + (j + 2)] = c02;
        }
        for(;j < bg; ++j) {
            c00 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
            }
            C[bg * (i) + (j)] = c00;
        }
    }
}

void block_mult_add(double* A, int av, int ag, double* B, int bg, double* C) {
    int av_reminder = av % 3;
    int bg_reminder = bg % 3;
    int _av = av - av_reminder;
    int _bg = bg - bg_reminder;
    int i = 0, j = 0, q = 0;
    double c00, c01, c02, c10, c11, c12, c20, c21, c22;
    for(i = 0; i < _av; i += 3) {
        for(j = 0; j < _bg; j+= 3) {
            c00 = 0; c01 = 0; c02 = 0;
            c10 = 0; c11 = 0; c12 = 0;
            c20 = 0; c21 = 0; c22 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
                c01 += A[ag * (i) + q] * B[bg * q + (j + 1)];
                c02 += A[ag * (i) + q] * B[bg * q + (j + 2)];
                c10 += A[ag * (i + 1) + q] * B[bg * q + (j)];
                c11 += A[ag * (i + 1) + q] * B[bg * q + (j + 1)];
                c12 += A[ag * (i + 1) + q] * B[bg * q + (j + 2)];
                c20 += A[ag * (i + 2) + q] * B[bg * q + (j)];
                c21 += A[ag * (i + 2) + q] * B[bg * q + (j + 1)];
                c22 += A[ag * (i + 2) + q] * B[bg * q + (j + 2)];
            }
            C[bg * (i) + (j)] += c00;
            C[bg * (i) + (j + 1)] += c01;
            C[bg * (i) + (j + 2)] +=c02;
            C[bg * (i + 1) + (j)] += c10;
            C[bg * (i + 1) + (j + 1)] += c11;
            C[bg * (i + 1) + (j + 2)] += c12;
            C[bg * (i + 2) + (j)] += c20;
            C[bg * (i + 2) + (j + 1)] += c21;
            C[bg * (i + 2) + (j + 2)] += c22;
        }
        for(;j < bg; ++j) {
            c00 = 0;
            c10 = 0;
            c20 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
                c10 += A[ag * (i + 1) + q] * B[bg * q + (j)];
                c20 += A[ag * (i + 2) + q] * B[bg * q + (j)];
            }
            C[bg * (i) + (j)] += c00;
            C[bg * (i + 1) + (j)] += c10;
            C[bg * (i + 2) + (j)] += c20;
        }
    }
    for(; i < av; ++i) {
        for(j = 0; j < _bg; j+= 3) {
            c00 = 0; c01 = 0; c02 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
                c01 += A[ag * (i) + q] * B[bg * q + (j + 1)];
                c02 += A[ag * (i) + q] * B[bg * q + (j + 2)];
            }
            C[bg * (i) + (j)] += c00;
            C[bg * (i) + (j + 1)] += c01;
            C[bg * (i) + (j + 2)] += c02;
        }
        for(;j < bg; ++j) {
            c00 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
            }
            C[bg * (i) + (j)] += c00;
        }
    }
}



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


void stick_matrix_to_core(Args* a) {
    int lines = a->n / a->m;
    for(int i = a->k; i < lines; i+=a->p) {
        memset(a->matrix1 + i * a->n * a->m, 0, a->n * a->m * sizeof(double));
        memset(a->matrix2 + i * a->n * a->m, 0, a->n * a->m * sizeof(double));
    }
    pthread_barrier_wait(a->barier);
}

void bench_mult(Args* a, double* block1) {
    double* matrix1 = a->matrix1;
    double* matrix2 = a->matrix2;
    int n = a->n;
    int m = a->m;
    int p = a->p;
    int lines = n / m;
    for(int i = a->k; i < lines; i+=(a->p)) {
        for(int j = 0; j < lines; ++j) {
            block_mult_add(matrix1 + i * n * m + j * m * m, m, m, matrix2 + j * n * m
                    , m, block1);
        }
    }
    //pthread_barrier_wait(a->barier);
    return;
}


void* thread_func(void* arg) {
    Args* a = (Args*)(arg);
    int m = a->m;

    double all_time = 0;
    double time_for_thread = 0;


    int nproc = get_nprocs();
    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    CPU_SET(nproc - 1 - a->k, &cpu);
    pthread_setaffinity_np(pthread_self(), sizeof(cpu), &cpu);
    //sched_setaffinity(pthread_self(), sizeof(cpu), &cpu);

    double* block1 = new double[m * m];
    memset(block1, 0, sizeof(double) * m * m);


    stick_matrix_to_core(a);

    all_time = get_full_time();
    time_for_thread = get_cpu_time();
    bench_mult(a, block1);
    time_for_thread = get_cpu_time() - time_for_thread;
    all_time = get_full_time() - all_time;
    a->time = time_for_thread;

    {
        std::lock_guard<std::mutex> lk(mut);
    //if (a -> k == 0) {
        std::cout << "Full time = " << all_time << std::endl;
    //}
    }
    delete[] block1;

    return nullptr;
}


int main(int argc, char* argv[]) {
    int n, m, p, s;
    pthread_barrier_t barrier;

    if (argc != 5 && argc != 6) {
        printf("Wrong amount of arguments\n");
        return 1;
    }
    if (!(sscanf(argv[1], "%d", &n) == 1 
        && sscanf(argv[2], "%d", &m) == 1
        && sscanf(argv[3], "%d", &p) == 1 
        && sscanf(argv[4], "%d", &s) == 1)) {
        printf("Wrong type of arguments\n");
        return 1;
    } 
    if (n <= 0 || m <= 0 || p < 0 || s < 0 || s > 4 || (s==0 && argc==5) || n < m) {
        printf("Wrong arguments\n");
        return 1;
    }
    if (n % m != 0) {
        std::cout << "n % m != 0" << std::endl;
        return 1;
    }
    if((n/m) % p != 0) {
        std::cout << "(n/m) % p != 0 " << std::endl;
        return 1;
    }

    double* matrix1 = new double[n * n];
    double* matrix2 = new double[n * n];
    Args* args = new Args[p];

    pthread_barrier_init(&barrier, 0, p);
    for(int i = 1; i < p; ++i) {
        args[i].p = p;
        args[i].k = i;
        args[i].matrix1 = matrix1;
        args[i].matrix2 = matrix2;
        args[i].n = n;
        args[i].m = m;
        args[i].barier = &barrier;
    }
    args[0].p = p;
    args[0].k = 0;
    args[0].matrix1 = matrix1;
    args[0].matrix2 = matrix2;
    args[0].n = n;
    args[0].m = m;
    args[0].barier = &barrier;
    for(int i = 1; i < p; ++i) {
        pthread_create(&(args[i].tid), nullptr, thread_func, args + i);
    }
    args[0].tid = pthread_self();
    thread_func(args);

    for(int i = 1; i < p; ++i) {
        pthread_join(args[i].tid, nullptr);
    }

    for(int i = 0; i < p; ++i) {
        std::cout << "Threda number " << i << " time = " << args[i].time << std::endl;
    }

    delete[] matrix1;
    delete[] matrix2;
    delete[] args;
    pthread_barrier_destroy(&barrier);
    return 0;
}



