#include "utils.hpp"
#include <immintrin.h>


double saclar_product(int n, double* x, double* y, int p, int k) {
    int i1 = -1, i2 = -1, i = -1;
    double sum = 0;
    double res12[] = {0, 0};
    thread_rows(n, p, k, i1, i2);

    __m128d xmm0 = _mm_setzero_pd();
    for(i = i1; i <= i2 - 2; i += 2) {
        __m128d v1 = _mm_loadu_pd(x + i);
        __m128d v2 = _mm_loadu_pd(y + i);
        __m128d prod = _mm_mul_pd(v1, v2);
        xmm0 = _mm_add_sd(xmm0, prod);
    }

    if(i < n) {
        __m128d v1 = _mm_load_sd(x + i);
        __m128d v2 = _mm_load_sd(y + i);
        __m128d prod = _mm_mul_sd(v1, v2);
        xmm0 = _mm_add_sd(xmm0, prod);
    }

    _mm_store_pd(res12, xmm0);

    sum = res12[0] + res12[1];

    sum = reduce_sum_double_det(p, k, sum);

    return sum;
}
