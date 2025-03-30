#pragma once

#include <pthread.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <vector>

template <typename T = int>
void reduce_sum(int p, T* a = nullptr, int n = 0) {
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static T* r = nullptr;
    int i;
    if(p <= 1) return;
    pthread_mutex_lock(&m);
    if(r == nullptr) {
        r = a;
    } else {
        for(i = 0; i < n; ++i) r[i] += a[i];
    }
    ++t_in;
    if(t_in >= p) {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    } else {
        while(t_in < p) {
            pthread_cond_wait(&c_in, &m);
        }
    }
    if(r != a) {
        for(i = 0; i < n; ++i) {
            a[i] = r[i];
        }
    }
    ++t_out;
    if(t_out >= p) {
        t_in = 0;
        r = nullptr;
        pthread_cond_broadcast(&c_out);
    } else {
        while (t_out < p) {
            pthread_cond_wait(&c_out, &m);
        }
    }
    pthread_mutex_unlock(&m);
}

inline double reduce_sum_double_det(int p, int k, double s) {

    struct alignas(128) alignad_double {
        double num = -1.;

    };

    static std::vector<alignad_double> results(p);

    double sum = 0;
    results[k].num = s;
    reduce_sum(p);

    for(int l = 0; l < p; ++l) {
        sum += results[l].num;
    }

    return sum;
}

inline double get_cpu_time() {
    struct rusage buf;
    getrusage(RUSAGE_THREAD, &buf);
    return (double)(buf.ru_utime.tv_sec) + (double)(buf.ru_utime.tv_usec)/1000000.;
}

inline double get_full_time() {
    struct timeval buf;
    gettimeofday(&buf, 0);
    return (double)(buf.tv_sec) + (double)(buf.tv_usec)/1000000.;
}

inline void thread_rows(int n, int p, int k, int& i1, int& i2) {
    i1 = n * k;
    i1 /= p;
    i2 = n * (k + 1);
    i2 /= p;
}
