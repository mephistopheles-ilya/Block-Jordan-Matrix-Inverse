#include "utils.hpp"
#include "solve.hpp"
#include <immintrin.h>

int min_residual_msr_matrix(int n, double* A, int* I, double* b, double* x, double* r, double* u, double* v
        , double eps, int maxit, int p, int k) {
    double prec = 0., b_norm2 = 0., tau = 0., c1 = 0., c2 = 0.;
    int it = 0;
    b_norm2 = scalar_product(n, b, b, p, k);
    prec = b_norm2 * eps * eps;
    matrix_mult_vector_msr(n, A, I, x, r, p, k);
    mult_sub_vector(n, r, b, 1., p, k);
    for(it = 0; it < maxit; ++it) {
        apply_preconditioner_msr_matrix(n, A, I, v, r, p, k);
        matrix_mult_vector_msr(n, A, I, v, u, p, k);
        c1 = scalar_product(n, u, r, p, k);
        c2 = scalar_product(n, u, u, p, k);
        if (c1 < prec || c2 < prec) {
            break;
        }
        tau = c1 / c2;
        mult_sub_vector(n, x, v, tau, p, k);
        mult_sub_vector(n, r, u, tau, p, k);
    }
    if(it > maxit) {
        return -1;
    }
    return it;
}

int min_residual_msr_matrix_full(int n, double* A, int* I, double* b, double* x, double* r, double* u, double* v
        , double eps, int maxit, int maxsteps, int p, int k) {
    int step = 0, ret = 0, its = 0;
    for(step = 0; step < maxsteps; ++step) {
        ret = min_residual_msr_matrix(n, A, I, b, x, r, u, v, eps, maxit, p, k);
        if (ret >= 0) {
            its += ret;
            break;
        }
        its += maxit;
    }
    if (step >= maxsteps) {
        return -1;
    }
    return its;
}

void matrix_mult_vector_msr(int n, double* A, int* I, double* x, double* y, int p, int k) {
    int i = 0, i1 = 0, i2 = 0, l = 0, J = 0;
    double s = 0;
    thread_rows(n, p, k, i1, i2);
    for(i = i1; i < i2; ++i) {
        s = A[i] * x[i];
        l = I[i + 1] - I[i];
        J = I[i];
        for(int j = 0; j < l; ++j) {
            s += A[J + j] * x[I[J + j]];
        }
        y[i] = s;
    }
    barrier(p);
}

double scalar_product(int n, double* x, double* y, int p, int k) {
    int i = 0, i1 = 0, i2 = 0;
    double s = 0;
    thread_rows(n, p, k, i1, i2);
    for(i = i1; i < i2; ++i) {
        s += x[i] * y[i];
    }
    reduce_sum_double_det(p, k, s);
    return s;
}

void mult_sub_vector(int n, double* x, double* y, double t, int p, int k) {
    int i = 0, i1 = 0, i2 = 0;
    thread_rows(n, p, k, i1, i2);
    for(i = i1; i < i2; ++i) {
        x[i] -= t * y[i];
    }
    barrier(p);
}

void apply_preconditioner_msr_matrix(int n, double* A, int* I, double* v, double* r, int p, int k) {
    int i = 0, i1 = 0, i2 = 0, l = 0, J = 0;
    int num = 0;
    double s = 0;
    thread_rows(n, p, k, i1, i2);
    for(i = i1; i < i2; ++i) {
        l = I[i + 1] - I[i];
        J = I[i];
        s = r[i];
        for(int j = 0; j < l; ++j) {
            num = I[J + j];
            if(num >= i1 && num < i) {
                s -= A[J + j] * v[num];
            }
        }
        v[i] = s / A[i];
    }

    for(i = i1; i < i2; ++i) {
        v[i] *= A[i];
    }
    
    for(i = i2 - 1; i >= i1; ++i) {
        l = I[i + 1] - I[i];
        J = I[i];
        s = v[i];
        for(int j = 0; j < l; ++j) {
            num = I[J + j];
            if (num > i && num < i2) {
                s -= A[J + j] * v[num];
            }
        }
        v[i] = s / A[i];
    }
    barrier(p);
}







#if 0
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
#endif 
