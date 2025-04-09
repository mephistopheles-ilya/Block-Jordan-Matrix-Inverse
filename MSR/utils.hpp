#pragma once

#include <pthread.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <vector>
#include <tuple>
#include <utility>
#include <algorithm>
#include <cstring>
#include <stdio.h>
#include <unistd.h>


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

    struct no_false_sharing_double {
        double nums[16] = {0};

    };

    static std::vector<no_false_sharing_double> results(p);

    double sum = 0;
    results[k].nums[0] = s;

    barrier(p);

    for(int l = 0; l < p; ++l) {
        sum += results[l].nums[0];
    }

    s = sum;

    barrier(p);
}



template <typename T>
void reduce_max(int p, T* a, int n) {
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
        for(i = 0; i < n; ++i) r[i] = std::max(r[i], a[i]);
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

template <typename... Args>
class fill_with_zeros {

private :

    std::tuple<Args*...> pointers;

    template <size_t... Is, typename... Nums>
    void apply_memset_impl(std::index_sequence<Is...>, Nums... nums) {
        expand_t((memset(std::get<Is>(pointers), 0, nums * 
                        sizeof(std::remove_pointer_t< std::remove_reference_t< decltype(std::get<Is>(pointers)) > >)), 0)...);
    }


public: 

    fill_with_zeros(Args*... args): pointers(args...) {}

    struct expand_t {
        expand_t(int...) {}
    };

    template <typename... Nums>
    void apply_memset(Nums... nums) {
        apply_memset_impl(std::make_index_sequence<sizeof...(Nums)>{}, nums...);
    }
};            


#if 0
inline void reduce_sum_two_double_det(int p, int k, double& s1, double& s2) {
    struct no_false_sharing_double {
        double nums[16] = {0};

    };

    static std::vector<no_false_sharing_double> results1(p);

    double sum1 = 0, sum2 = 0;
    results1[k].nums[0] = s1;
    results1[k].nums[1] = s2;

    barrier(p);

    for(int l = 0; l < p; ++l) {
        sum1 += results1[l].nums[0];
        sum2 += results1[l].nums[1];
    }

    s1 = sum1;
    s2 = sum2;

    barrier(p);
}

inline void make_fpe(double x = 0) {
    printf("%lf\n", 1./x);
}

inline void handler(int /*n*/) {
    char msg[] = "got signal\n";
    auto ret = write(1, msg, sizeof(msg) - 1);
    (void)ret;
}




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




template <typename T, size_t alignment> 
struct aligned_allocator {
    typedef T value_type;
    typedef T *pointer;

    pointer allocate(size_t n) {
        void *chunk = std::aligned_alloc(alignment,  n * sizeof(value_type));
        return static_cast<pointer>(chunk);
    }

    void deallocate(pointer p, size_t) {
        free(p);
    }

    template <typename U> 
    aligned_allocator(const aligned_allocator<U, alignment> &) {}

    aligned_allocator() {}
    aligned_allocator(const aligned_allocator &) {}

    template <typename U> 
    struct rebind { typedef aligned_allocator<U, alignment> other; };

};

#endif

