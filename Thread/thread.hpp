#pragma once

#include <pthread.h>

enum class io_status {
    undef,
    error_open,
    error_read,
    not_enough_elements,
    succes
};

struct Arg {
    int n = 0;
    int m = 0;
    int p = 0;
    int r = 0;
    int s = 0;
    char* filename = nullptr;

    int thread_number = 0;
    double* matrix = nullptr;
    double* inverse = nullptr;

    int error = 0;

    pthread_t tid = 0;
    pthread_barrier_t* barrier = nullptr;
};


void * thread_func(void* arg);


    
