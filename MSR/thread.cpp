#include "thread.hpp"
#include "utils.hpp"

void* thread_func(void* arg ) {
    Arg* a = (Arg *)arg;
    int p = a->p;
    int thr_num = a->thr_num;
    barrier(p);
    reduce_sum_double_det(p, thr_num, 1);
    barrier(p);
    reduce_sum_double_det(p, thr_num, 1);
    barrier(p);
    reduce_sum_double_det(p, thr_num, 1);
    barrier(p);
    reduce_sum_double_det(p, thr_num, 1);

    return nullptr;
};
