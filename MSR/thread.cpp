#include <cmath>

#include "thread.hpp"
#include "utils.hpp"
#include "fill_msr.hpp"
#include "solve.hpp"
#include "residuals.hpp"

void* thread_func(void* argument ) {
    Arg* arg = (Arg *)argument;
    const int maxsteps = 100;
    double a = arg->a;
    double b = arg->b;
    double c = arg->c;
    double d = arg->d;
    int maxit = arg->mi;
    int nx = arg->nx;
    int ny = arg->ny;
    double eps = arg->eps;
    int p = arg->p;
    int thr_num = arg->thr_num;
    int func_num = arg->k;
    double* A = arg->A;
    double* B = arg->B;
    int* I = arg->I;
    double* x = arg->x;
    double* r = arg->r;
    double* u = arg->u;
    double* v = arg->v;
    double (*f)(double, double) = get_funk(func_num);
    double t1 = -1, t2 = -1, r1 = -1, r2 = -1, r3 = -1, r4 = -1;
    int it = 0;

    double hx = (b - a) / nx;
    double hy = (d - c) / ny;

    int len_diag = (nx + 1) * (ny + 1);
    if (thr_num == 0) {
        fill_I(nx, ny, I);
    }
    barrier(p);

    fill_A(nx, ny, hx, hy, I, A, p, thr_num);
    fill_B(a, c, nx, ny, hx, hy, B, p, thr_num, f);

    barrier(p);
    t1 = get_full_time();
    it = min_residual_msr_matrix_full(len_diag, A, I, B, x, r, u, v, eps, maxit, maxsteps, p, thr_num); 
    barrier(p);
    t1 = get_full_time() - t1;

    barrier(p);
    t2 = get_full_time();
    r1 = calc_r1(x, a, c, hx, hy, nx, ny, p, thr_num, f);
    r2 = calc_r2(x, a, c, hx, hy, nx, ny, p, thr_num, f);
    r3 = calc_r3(x, a, c, hx, hy, nx, ny, p, thr_num, f);
    r4 = calc_r4(x, a, c, hx, hy, nx, ny, p, thr_num, f);
    barrier(p);
    t2 = get_full_time() - t2;

    arg->t1 = t1;
    arg->t2 = t2;
    arg->r1 = r1;
    arg->r2 = r2;
    arg->r3 = r3;
    arg->r4 = r4;
    arg->it = it;

    return nullptr;
};


double (*get_funk(int num))(double, double) {
    static auto f0 = [](double /*x*/, double /*y*/) {
        return 1.;
    };
    static auto f1 = [](double x, double /*y*/) {
        return x;
    };
    static auto f2 = [](double /*x*/, double y) {
        return y;
    };
    static auto f3 = [](double x, double y) {
        return x + y;
    };
    static auto f4 = [](double x, double y) {
        return std::sqrt(x * x + y * y);
    };
    static auto f5 = [](double x, double y) {
        return x * x + y * y;
    };
    static auto f6 = [](double x, double y) {
        return std::exp(x * x - y * y);
    };
    static auto f7 = [](double x, double y) {
        return 1./(25. * (x * x + y * y) + 1);
    };

    switch (num) {
        case 0: 
            return f0;
        case 1:
            return f1;
        case 2:
            return f2;
        case 3:
            return f3;
        case 4:
            return f4;
        case 5:
            return f5;
        case 6:
            return f6;
        case 7:
            return f7;
        default :
            return nullptr;
    }
    return nullptr;
}
