#include <ctime>
#include <stdio.h>
#include <sched.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#include <stdlib.h>
#include <new>

#include "thread.hpp"

int main(int argc, char* argv[]) {
    int n = 0, m = 0, p = 0, r = 0, s = 0;
    pthread_barrier_t barrier;


    if (argc != 6 && argc != 7) {
        printf("Wrong amount of arguments\n");
        return 1;
    }
    if (!(sscanf(argv[1], "%d", &n) == 1 
        && sscanf(argv[2], "%d", &m) == 1
        && sscanf(argv[3], "%d", &p) == 1 
        && sscanf(argv[4], "%d", &r) == 1
        && sscanf(argv[5], "%d", &s) == 1)) {
        printf("Wrong type of arguments\n");
        return 1;
    } 
    if (n <= 0 || m <= 0 || r < 0 || s < 0 || s > 4 || (s==0 && argc==6) || n < m || p <= 0) {
        printf("Wrong arguments\n");
        return 1;
    }

    printf("Usage %s %d %d %d %d %d %s\n", argv[0], n, m, p, r, s, argv[6]);


    double* matrix = new(std::nothrow) double[n * n];
    double* inverse = new(std::nothrow) double[n * n];
    Arg* args = new(std::nothrow) Arg[p];
    if (matrix == nullptr || inverse == nullptr || args == nullptr) {
        printf("Not enough memmory\n");
        delete[] matrix;
        delete[] inverse;
        delete[] args;
        return 2;
    }

    pthread_barrier_init(&barrier, nullptr, p);

    for(int thread_number = 0; thread_number < p; ++thread_number) {
        args[thread_number].n = n;
        args[thread_number].m = m;
        args[thread_number].p = p;
        args[thread_number].r = r;
        args[thread_number].s = s;
        args[thread_number].filename = argv[6];
        args[thread_number].thread_number = thread_number; 
        args[thread_number].matrix = matrix;
        args[thread_number].inverse = inverse;
        args[thread_number].barrier = &barrier;
    }

    for(int thread_number = 1; thread_number < p; ++thread_number) {
        if (pthread_create(&(args[thread_number].tid), nullptr, thread_func, args + thread_number) != 0) {
            printf("Cannot create new thread\n");
            for(int l = 1; l < thread_number; ++l) {
                pthread_join(args[thread_number].tid, nullptr);
            }
            pthread_barrier_destroy(&barrier);
            delete[] matrix;
            delete[] inverse;
            delete[] args;
            return 3;
        }
    }
    args[0].tid = pthread_self();
    thread_func(args + 0);
    for(int thread_number = 1; thread_number < p; ++thread_number) {
        pthread_join(args[thread_number].tid, nullptr);
    }
    int task = 19;
    double r1 = args[0].r1;
    double r2 = args[0].r2;
    double t1 = args[0].global_time_inv;
    double t2 = args[0].global_time_dis;
    printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
            argv[0], task, r1, r2, t1, t2, s, n, m, p);

#if 0
    for(int i = 0; i < p; ++i) {
        printf("Time for inv thread %2d : %lf\n", i, args[i].local_time_inv);
    }
    for(int i = 0; i < p; ++i) {
        printf("Time for dis thread %2d : %lf\n", i, args[i].local_time_dis);
    }

#endif

#if 0
    double sum1 = 0;
    for(int i = 0; i < p; ++i) {
        sum1 += args[i].local_time_inv * args[i].local_time_inv;
    }
    sum1 /= p;
    double sum2 = 0;
    for(int i = 0; i < p; ++i) {
        sum2 += args[i].local_time_inv;
    }
    sum2 /= p;
    printf("D : %.20lf\n", sum1 - sum2 * sum2);
#endif

    pthread_barrier_destroy(&barrier);
    delete[] matrix;
    delete[] inverse;
    delete[] args;

    return 0;
}
