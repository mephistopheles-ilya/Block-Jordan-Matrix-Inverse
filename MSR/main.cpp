#include <stdio.h>
#include <new>

#include "thread.hpp"
#include "fill_msr.hpp"
#include "utils.hpp"


#include <fenv.h>




int main(int argc, char* argv[]) {

    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

#if 0
    struct sigaction sa = {};
    sa.sa_handler = handler;
    sa.sa_flags = SA_RESTART | SA_NODEFER;
    sigaction(SIGABRT, &sa, NULL);
#endif

    double a = -1., b = -1., c = -1., d = -1.;
    int nx = -1, ny = -1;
    int k = -1;
    double eps = -1;
    int mi = -1;
    int p = -1;


    if (argc != 11 
        || sscanf(argv[1], "%lf", &a) != 1 || sscanf(argv[2], "%lf", &b) != 1 
        || sscanf(argv[3], "%lf", &c) != 1 || sscanf(argv[4], "%lf", &d) != 1
        || sscanf(argv[5], "%d", &nx) != 1 || sscanf(argv[6], "%d", &ny) != 1
        || sscanf(argv[7], "%d", &k) != 1 || sscanf(argv[8], "%lf", &eps) != 1
        || sscanf(argv[9], "%d", &mi) != 1 || sscanf(argv[10], "%d", &p) != 1) {

        printf("Usage %s a b c d nx ny k eps mi p\n", argv[0]);
        return 1;
    }

    if (b <= a || d <= c || nx < 2 || ny < 2 || k < 0 || k > 7 || eps < 0 || mi < 0 || p <= 0) {
        printf("Wrong parametrs of commad line\n");
        return 2;
    }

    const int len_msr = get_len_msr(nx, ny) + 1;
    const int len_diag = (nx + 1) * (ny + 1);
    double* A = new (std::nothrow) double[len_msr];
    double* B = new (std::nothrow) double[len_diag];
    int* I  = new (std::nothrow) int[len_msr];
    double* x = new (std::nothrow) double[len_diag];
    double* r = new (std::nothrow) double[len_diag];
    double* u = new (std::nothrow) double[len_diag];
    double* v = new (std::nothrow) double[len_diag];
    Arg* args = new (std::nothrow) Arg[p];

    memmory_manager<double, double,  int, double, double, double, double, Arg> mm(A, B, I, x, r, u, v, args);

    if (mm.check_nullptr() == true) {
        printf("Not enough memmory\n");
        return 3;
    }

#if 0
    fill_with_zeros<double, double, int, double, double, double, double> fwz(A, B, I, x, r, u, v);
    fwz.apply_memset(len_msr, len_diag, len_msr, len_diag, len_diag, len_diag, len_diag);
#endif
    memset(x, 0, len_diag * sizeof(double));

    for(int thr_num = 0; thr_num < p; ++thr_num) {
        args[thr_num].a = a;
        args[thr_num].b = b;
        args[thr_num].c = c;
        args[thr_num].d = d;
        args[thr_num].nx = nx;
        args[thr_num].ny = ny;
        args[thr_num].k = k;
        args[thr_num].eps = eps;
        args[thr_num].mi = mi;
        args[thr_num].p = p;
        args[thr_num].thr_num = thr_num;
        args[thr_num].A = A;
        args[thr_num].B = B;
        args[thr_num].I = I;
        args[thr_num].x = x;
        args[thr_num].r = r;
        args[thr_num].u = u;
        args[thr_num].v = v;
    }

    for(int thr_num = 1; thr_num < p; ++thr_num) {
        if (pthread_create(&(args[thr_num].tid), nullptr, thread_func, args + thr_num) != 0) {
            printf("Cannot create new thread\n");
            for(int i = 1; i < thr_num; ++i) {
                pthread_join(args[thr_num].tid, nullptr);
            }
            return 4;
        }
    }

    args[0].tid = pthread_self();
    thread_func(args + 0);

    for(int thr_num = 1; thr_num < p; ++thr_num) {
        pthread_join(args[thr_num].tid, nullptr);
    }

    const int task = 2;
    double r1 = args[0].r1;
    double r2 = args[0].r2;
    double r3 = args[0].r3;
    double r4 = args[0].r4;
    double t1 = args[0].t1;
    double t2 = args[0].t2;
    int it = args[0].it;
    
    printf ("%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f \n\t It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
            argv[0], task, r1, r2, r3, r4, t1, t2, it, eps, k, nx, ny, p);

    return 0;
}
