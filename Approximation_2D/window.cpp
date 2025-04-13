#include "window.hpp"

#include "fill_msr.hpp"

#include "utils.hpp"

#include "solve.hpp"

#include "residuals.hpp"

#include <cstring>
#include <sys/sysinfo.h>
#include <QTimer>
#include <QCloseEvent>

MainWindow::MainWindow(double a, double b, double c, double d, int nx, int ny
        ,int mx, int my, int k, double eps, int mi, int p): 
    a(a), b(b), c(c), d(d), nx(nx), ny(ny), mx(mx), my(my), k(k), eps(eps), mi(mi), p(p) {
    
    int len_msr = get_len_msr(nx, ny) + 1;
    int len_diag = (nx + 1) * (ny + 1);
    args = new (std::nothrow) Arg[p];
    allocate_memmory(len_msr, len_diag);
    if(check_nullptr() == true) {
        deallocate_memmory();
        QTimer::singleShot(0, this, [this]() {
                this->close();
                }
        );
    }

    timer = new QTimer(this);
    connect(timer, &QTimer::timeout, this, &MainWindow::time_out_checker);
    timer -> start(100);
}

void MainWindow::time_out_checker() {
    pthread_mutex_lock(&p_mutex);
    if (calculations_done == true) {
        memcpy(x_gui, x_msr, (nx + 1) * (ny + 1) * sizeof(double));
    }
    pthread_mutex_unlock(&p_mutex);
}



void MainWindow::create_threads() {
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
        args[thr_num].x = x_gui;
        args[thr_num].r = r;
        args[thr_num].u = u;
        args[thr_num].v = v;
        args[thr_num].window = this;
    }
    for(int thr_num = 0; thr_num < p; ++thr_num) {
        if (pthread_create(&(args[thr_num].tid), nullptr, create_msr_approximation, args + thr_num) != 0) {
            printf("Cannot create new thread\n");
            abort();
        }
    }
    for(int thr_num = 0; thr_num < p; ++thr_num) {
        pthread_join(args[thr_num].tid, nullptr);
    }
    this->close();
}



void MainWindow::allocate_memmory(int len_msr, int len_diag) {
    A = new (std::nothrow) double[len_msr];
    B = new (std::nothrow) double[len_diag];
    I  = new (std::nothrow) int[len_msr];
    x_gui = new (std::nothrow) double[len_diag];
    x_msr = new (std::nothrow) double[len_diag];
    r = new (std::nothrow) double[len_diag];
    u = new (std::nothrow) double[len_diag];
    v = new (std::nothrow) double[len_diag];
    memset(x_gui, 0, len_diag * sizeof(double));
}

void MainWindow::deallocate_memmory() {
    delete[] A;
    delete[] B;
    delete[] I;
    delete[] x_gui;
    delete[] x_msr;
    delete[] r;
    delete[] u;
    delete[] v;
}

bool MainWindow::check_nullptr() const {
    if (A == nullptr || B == nullptr || I == nullptr || x_gui == nullptr || x_msr == nullptr
            || r == nullptr || u == nullptr || v == nullptr || args == nullptr) {
        return true;
    } 
    return false;
}

void MainWindow::closeEvent(QCloseEvent* event) {
    int local_caluclations_done = true;
    pthread_mutex_lock(&p_mutex);
    msr_continue = false;
    if (calculations_done == false) {
        local_caluclations_done = false;
    }
    pthread_mutex_unlock(&p_mutex);    

    if (local_caluclations_done == false) {
        event->ignore();
        return;
    }
    pthread_cond_broadcast(&p_cond);

    deallocate_memmory();
    event->accept();
}




void* MainWindow::create_msr_approximation(void* argument) {
    Arg* arg = (Arg *)argument;
    const int maxit_no_restart = 50;
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

    int nproc = get_nprocs();
    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    CPU_SET(nproc - 1 - (thr_num % (nproc)), &cpu);
    pthread_setaffinity_np(pthread_self(), sizeof(cpu), &cpu);

    double hx = (b - a) / nx;
    double hy = (d - c) / ny;

    int len_diag = (nx + 1) * (ny + 1);
    int t_in = 0;
    static bool local_msr_constinue = true;
    while(true) {

        barrier(p);
        t1 = get_full_time();

        if (thr_num == 0) {
            fill_I(nx, ny, I);
        }
        barrier(p);

        fill_A(nx, ny, hx, hy, I, A, p, thr_num);
        fill_B(a, c, nx, ny, hx, hy, B, p, thr_num, f);

        it = min_residual_msr_matrix_full(len_diag, A, I, B, x, r, u, v, eps, maxit, maxit_no_restart, p, thr_num);
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

        pthread_mutex_lock(&arg->window->p_mutex);
        ++t_in;
        if (t_in >= p) {
            t_in = 0;
            if (arg->window->msr_continue == false) {
                local_msr_constinue = false;
            }
            printf ("%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f \n\t It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
                    "a.out", 3, r1, r2, r3, r4, t1, t2, it, eps, func_num, nx, ny, p);
            arg->window->calculations_done = true;

            pthread_cond_wait(&arg->window->p_cond, &arg->window->p_mutex);
        } else {
            pthread_cond_wait(&arg->window->p_cond, &arg->window->p_mutex);
        }

        if (local_msr_constinue == false) {
            break;
        }
    }


    return nullptr;
}


