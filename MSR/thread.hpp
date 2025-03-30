#pragma once

#include <pthread.h>

struct Arg {
    double a = -1., b = -1., c = -1., d = -1.;
    int nx = -1, ny = -1;
    int k = -1;
    double eps = -1;
    int mi = -1;
    int p = -1;
    int thr_num = -1;
    
    pthread_t tid = 0;

    double r1 = -1., r2 = -1., r3 = -1., r4 = -1.;
    double t1 = -1., t2 = -1.;
    int it = -1;

};

void* thread_func(void* arg);


