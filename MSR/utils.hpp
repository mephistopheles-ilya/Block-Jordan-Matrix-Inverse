#pragma once

#include <pthread.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <vector>
#include <tuple>
#include <utility>


template <typename... Args>
class memmory_manager {

private:

    std::tuple<Args*...> pointers;
    int nullptr_counter = 0;

    struct expand_t {
        expand_t(int...) {}
    };

    template<typename T>
    void is_nullptr(T* ptr) {
        if (ptr == nullptr) ++nullptr_counter;
    }

    template<size_t... N>
    void deleter(std::integer_sequence<size_t, N...>) {
        expand_t((delete[] std::get<N>(pointers), 0)...);
        expand_t((std::get<N>(pointers) = nullptr, 0)...);
    }

    template<size_t... N>
    void count_nullptr(std::integer_sequence<size_t, N...>) {
        expand_t((is_nullptr(std::get<N>(pointers)), 0)...);
    }


public:

    memmory_manager(Args*... args): pointers(args...) {}

    bool check_nullptr() {
        count_nullptr(std::make_index_sequence<sizeof...(Args)>{});
        if (nullptr_counter != 0) return true;
        return false;
    }

    ~memmory_manager() {
        deleter(std::make_index_sequence<sizeof...(Args)>{});
    }
};

        
        
            


template <typename T>
void reduce_sum(int p, T* a, int n) {
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

inline void barrier(int p) {

    struct closure_barrier {
        pthread_barrier_t barrier;
        closure_barrier(int p) {
            pthread_barrier_init(&barrier, nullptr, p);
        }
        ~closure_barrier() {
            pthread_barrier_destroy(&barrier);
        }
    };

    static closure_barrier bar(p);

    pthread_barrier_wait(&bar.barrier);

}

inline void reduce_sum_double_det(int p, int k, double& s) {

    struct alignas(128) alignad_double {
        double num = -1.;

    };

    static std::vector<alignad_double> results(p);

    double sum = 0;
    results[k].num = s;
    barrier(p);

    for(int l = 0; l < p; ++l) {
        sum += results[l].num;
    }

    s = sum;
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
