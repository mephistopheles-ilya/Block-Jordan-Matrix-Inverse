#include "window.hpp"
#include "fill_msr.hpp"
#include "utils.hpp"
#include "solve.hpp"
#include "residuals.hpp"

#include <iostream>
#include <limits>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <pthread.h>
#include <sys/sysinfo.h>
#include <QTimer>
#include <QCloseEvent>
#include <QPainter>
#include <QMessageBox>
#include <QMenuBar>

MainWindow::MainWindow(double a, double b, double c, double d, int nx, int ny
        ,int mx, int my, int k, double eps, int mi, int p, const char* prog_name) {
    plot_task.a = a;
    plot_task.b = b;
    plot_task.c = c;
    plot_task.d = d;
    plot_task.nx = nx;
    plot_task.ny = ny;
    plot_task.mx = mx;
    plot_task.my = my;
    plot_task.func_num = k;
    plot_task.f = get_funk(k);
    plot_task.s_a = a;
    plot_task.s_b = b;
    plot_task.s_c = c;
    plot_task.s_d = d;
    plot_task.current_paint = what_to_paint::function;
    plot_task.x_approximation = new double[(nx + 1) * (ny + 1)];

    msr_task.nx = nx;
    msr_task.ny = ny;
    msr_task.func_num = k;
    msr_task.f = get_funk(k);
    msr_task.condition = msr_condition::has_task;

    args = new Arg[p];
    for(int thr_num = 0; thr_num < p; ++thr_num) {
        args[thr_num].thr_num = thr_num;
        args[thr_num].eps = eps;
        args[thr_num].maxit = mi;
        args[thr_num].p = p;
        args[thr_num].a = a;
        args[thr_num].b = b;
        args[thr_num].c = c;
        args[thr_num].d = d;
        args[thr_num].prog_name = prog_name;
        args[thr_num].task = &msr_task;
        args[thr_num].p_mutex = &p_mutex;
        args[thr_num].p_cond = &p_cond;
    }

    setWindowTitle("Approximation 2D");
    //QMenuBar* menu = new QMenuBar(this);
    //setMenuBar(menu);
    //info = menu->addMenu("Calculatiing...");
    status_label = new QLabel(this);
    status_label->setAlignment(Qt::AlignCenter);
    menuBar()->setCornerWidget(status_label, Qt::TopLeftCorner);
    status_label->setText("Calculatiing...");
    status_label->setMinimumWidth(500);
    status_label->setStyleSheet("QLabel { padding: 3px; background: #f0f0f0; }");
        

    timer = new QTimer(this);
    connect(timer, &QTimer::timeout, this, &MainWindow::time_out_checker);
    timer -> start(200);

    for(int thr_num = 0; thr_num < p; ++thr_num) {
        if (pthread_create(&(args[thr_num].tid), nullptr, create_msr_approximation, args + thr_num) != 0) {
            printf("Cannot create new thread\n");
            abort();
        }
    }
}

void MainWindow::updateMenuTitle() {
    get_max_min_value(plot_task);
    QString title;
    int func_num = plot_task.func_num;
    int nx = plot_task.nx, ny = plot_task.ny;
    int mx = plot_task.mx, my = plot_task.my;
    int s = plot_task.s;
    int p = plot_task.inaccuracy;
    double f_abs_max = plot_task.f_abs_max;

    QString mode;
    switch (plot_task.current_paint) {
        case what_to_paint::function:
            mode = "Function";
            break;
        case what_to_paint::approximation:
            mode = "Approximation";
            break;
        case what_to_paint::residual:
            mode = "Residual";
            break;
    }

    switch (func_num) {
        case 0: 
            title = QString("k = 0 f(x, y) = 1, nx = %1, ny = %2, mx = %3, my = %4, s = %5, p = %6, max{|F_min|, |F_max|} = %7, Mode : %8")
                .arg(nx).arg(ny).arg(mx).arg(my).arg(s).arg(p).arg(f_abs_max).arg(mode);
            break;
        case 1: 
            title = QString("k = 1 f(x, y) = x, nx = %1, ny = %2, mx = %3, my = %4, s = %5, p = %6, max{|F_min|, |F_max|} = %7, Mode : %8")
                .arg(nx).arg(ny).arg(mx).arg(my).arg(s).arg(p).arg(f_abs_max).arg(mode);
           break;
        case 2: 
            title = QString("k = 2 f(x, y) = y, nx = %1, ny = %2, mx = %3, my = %4, s = %5, p = %6, max{|F_min|, |F_max|} = %7, Mode : %8")
                .arg(nx).arg(ny).arg(mx).arg(my).arg(s).arg(p).arg(f_abs_max).arg(mode);
           break;
        case 3: 
            title = QString("k = 3 f(x, y) = x + y, nx = %1, ny = %2, mx = %3, my = %4, s = %5, p = %6, max{|F_min|, |F_max|} = %7, Mode : %8")
                .arg(nx).arg(ny).arg(mx).arg(my).arg(s).arg(p).arg(f_abs_max).arg(mode);
           break;
        case 4: 
            title 
                = QString("k = 4 f(x, y) = sqrt(x^2 + y^2), nx = %1, ny = %2, mx = %3, my = %4, s = %5, p = %6, max{|F_min|, |F_max|} = %7, Mode : %8")
                .arg(nx).arg(ny).arg(mx).arg(my).arg(s).arg(p).arg(f_abs_max).arg(mode);
           break;
        case 5: 
            title = QString("k = 5 f(x, y) = x^2 + y^2, nx = %1, ny = %2, mx = %3, my = %4, s = %5, p = %6, max{|F_min|, |F_max|} = %7, Mode : %8")
                .arg(nx).arg(ny).arg(mx).arg(my).arg(s).arg(p).arg(f_abs_max).arg(mode);
           break;
        case 6: 
            title = QString("k = 6 f(x, y) = exp(x^2 - y^2), nx = %1, ny = %2, mx = %3, my = %4, s = %5, p = %6, max{|F_min|, |F_max|} = %7, Mode : %8")
                .arg(nx).arg(ny).arg(mx).arg(my).arg(s).arg(p).arg(f_abs_max).arg(mode);
           break;
        case 7: 
            title 
            = QString("k = 7 f(x) = 1/(25*(x^2 + y^2) + 1), nx = %1, ny = %2, mx = %3, my = %4, s = %5, p = %6, max{|F_min|, |F_max|} = %7, Mode : %8")
                .arg(nx).arg(ny).arg(mx).arg(my).arg(s).arg(p).arg(f_abs_max).arg(mode);
           break;
    }
    //info->setTitle(title);
    status_label->setText(title);
    status_label->adjustSize();
    update();
}

MainWindow::~MainWindow() {
    int p = args[0].p;
    for(int thr_num = 0; thr_num < p; ++thr_num) {
        pthread_join(args[thr_num].tid, nullptr);
    }
    delete[] plot_task.x_approximation;
    delete[] args;
}

QSize MainWindow::minimumSizeHint() const {
    return QSize(100, 100);
}

QSize MainWindow::sizeHint() const {
    return QSize(1000, 1000);
}

void MainWindow::time_out_checker() {
    bool need_update = false;
    pthread_mutex_lock(&p_mutex);
    if (msr_task.condition == msr_condition::task_is_ready) {
        need_update = true;
        msr_task.condition = msr_condition::no_task;
        //copy answr and task parametrs
        int current_len = (plot_task.nx + 1) * (plot_task.ny + 1);
        int new_len = (msr_task.nx + 1) * (msr_task.ny + 1);
        if (new_len > current_len) {
            delete[] plot_task.x_approximation;
            plot_task.x_approximation = new double[new_len];
        }
        plot_task.nx = msr_task.nx;
        plot_task.ny = msr_task.ny;
        plot_task.func_num = msr_task.func_num;
        plot_task.f = msr_task.f;
        plot_task.inaccuracy = msr_task.inaccuracy;
        memcpy(plot_task.x_approximation, msr_task.x_approximation, new_len * sizeof(double)); 
        plot_task.is_ready = true;
    }
    pthread_mutex_unlock(&p_mutex);
    if (need_update == true) {
        updateMenuTitle();
    }
}



void MainWindow::closeEvent(QCloseEvent* event) {
    bool accept = false;
    pthread_mutex_lock(&p_mutex);
    if (msr_task.condition == msr_condition::no_task 
            || msr_task.condition == msr_condition::task_is_ready) {
        msr_task.condition = msr_condition::quit_app;
        pthread_cond_broadcast(&p_cond);
        accept = true;
    }
    pthread_mutex_unlock(&p_mutex);
    if (accept == true) {
        event->accept();
    } else {
        QMessageBox::information(this, "Calculations in progress...", "Wait until the calculations end");
        event->ignore();
    }
}




void* MainWindow::create_msr_approximation(void* argument) {
    const int maxit_no_restart = 50;
    static data_to_msr glob_task;
    data_to_msr task;

    static double* glob_A = nullptr, *glob_B = nullptr, *glob_u = nullptr, *glob_v = nullptr, *glob_r = nullptr, *glob_x = nullptr;
    static int* glob_I = nullptr;

    Arg* arg = (Arg *)argument;
    const int p = arg->p;
    const int thr_num = arg->thr_num;
    const int maxit = arg->maxit;
    const double eps = arg->eps;
    const double a = arg->a, b = arg->b, c = arg->c, d = arg->d;
    const char* prog_name = arg->prog_name;
    data_to_msr* const pointer_task = arg->task;

    double* A = nullptr, *B = nullptr, *u = nullptr, *v = nullptr, *r = nullptr, *x = nullptr;
    int* I = nullptr;



    pthread_mutex_t* p_mutex = arg->p_mutex;
    pthread_cond_t* p_cond = arg->p_cond;

    double t1 = 0, t2 = 0, r1 = 0, r2 = 0, r3 = 0, r4 = 0;
    double hx = 0, hy = 0;
    int it = 0, len_diag = 0;

    static bool quit_app = false;
    static int t_in = 0;
    while(true) {
        if (thr_num == 0) {
            //get_task and signal to stop
            pthread_mutex_lock(p_mutex);
            glob_task = *pointer_task;
            pthread_mutex_unlock(p_mutex);
            if (glob_task.condition == msr_condition::has_task) {
                //allocate memmory for task
                delete[] glob_A;
                delete[] glob_B;
                delete[] glob_u;
                delete[] glob_v;
                delete[] glob_r;
                delete[] glob_x;
                delete[] I;
                int len_msr = get_len_msr(glob_task.nx, glob_task.ny) + 1;
                int len_diag = (glob_task.nx + 1) * (glob_task.ny + 1);
                glob_A = new (std::nothrow) double[len_msr];
                glob_B = new (std::nothrow) double[len_diag];
                glob_I  = new (std::nothrow) int[len_msr];
                glob_x = new (std::nothrow) double[len_diag];
                glob_r = new (std::nothrow) double[len_diag];
                glob_u = new (std::nothrow) double[len_diag];
                glob_v = new (std::nothrow) double[len_diag];
                memset(glob_x, 0, len_diag * sizeof(double));
            }

        }
        barrier(p);

        //get_task for every thread
        task = glob_task;
        A = glob_A;
        B = glob_B;
        I = glob_I;
        x = glob_x;
        r = glob_r;
        u = glob_u;
        v = glob_v;

        hx = (b - a) / task.nx;
        hy = (d - c) / task.ny;
        len_diag = (task.nx + 1) * (task.ny + 1);


        if (task.condition == msr_condition::quit_app) {
            if (thr_num == 0) {
                //free memmory
                delete[] glob_A;
                delete[] glob_B;
                delete[] glob_u;
                delete[] glob_v;
                delete[] glob_r;
                delete[] glob_x;
                delete[] I;
            }
            break;
        }

        if (task.condition == msr_condition::has_task) {

            barrier(p);
            t1 = get_full_time();

            if (thr_num == 0) {
                fill_I(task.nx, task.ny, I);
            }
            barrier(p);

            double max_f_value = 0;
            fill_A(task.nx, task.ny, hx, hy, I, A, p, thr_num);
            max_f_value = find_max(a, c, hx, hy, task.nx, task.ny, p, thr_num, task.f);
            double add_error = task.inaccuracy * 0.1 * max_f_value;
            fill_B(a, c, task.nx, task.ny, hx, hy, B, p, thr_num, task.f, add_error);

            it = min_residual_msr_matrix_full(len_diag, A, I, B, x, r, u, v, eps, maxit, maxit_no_restart, p, thr_num);
            barrier(p);
            t1 = get_full_time() - t1;

            barrier(p);
            t2 = get_full_time();
            r1 = calc_r1(x, a, c, hx, hy, task.nx, task.ny, p, thr_num, task.f);
            r2 = calc_r2(x, a, c, hx, hy, task.nx, task.ny, p, thr_num, task.f);
            r3 = calc_r3(x, a, c, hx, hy, task.nx, task.ny, p, thr_num, task.f, 0);
            r4 = calc_r4(x, a, c, hx, hy, task.nx, task.ny, p, thr_num, task.f, 0);
            barrier(p);
            t2 = get_full_time() - t2;


            if (thr_num == 0) {
                pthread_mutex_lock(p_mutex);
                if (pointer_task->condition == msr_condition::quit_app) {
                    quit_app = true;
                }
                pointer_task->condition = msr_condition::task_is_ready;
                pointer_task->x_approximation = x;
                printf ("%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f \n\t It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
                        prog_name, 3, r1, r2, r3, r4, t1, t2, it, eps, task.func_num, task.nx, task.ny, p);
                pthread_mutex_unlock(p_mutex);
            }
            barrier(p);
            if (quit_app == true) {
                if (thr_num == 0) {
                    //free memmory
                    delete[] glob_A;
                    delete[] glob_B;
                    delete[] glob_u;
                    delete[] glob_v;
                    delete[] glob_r;
                    delete[] glob_x;
                    delete[] I;
                }
                break;
            }
        }

        pthread_mutex_lock(p_mutex);
        ++t_in;
        if (t_in >= p) {
            t_in = 0;
            if (pointer_task->condition == msr_condition::no_task || pointer_task->condition == msr_condition::task_is_ready ) {
                pthread_cond_wait(p_cond, p_mutex);
            }
        } else {
            pthread_cond_wait(p_cond, p_mutex);
        }
        pthread_mutex_unlock(p_mutex);
    }


    return nullptr;
}

void MainWindow::keyPressEvent(QKeyEvent* event) {
    if (event->key() == Qt::Key_0) {
        pthread_mutex_lock(&p_mutex);
        msr_condition cond = msr_task.condition;
        pthread_mutex_unlock(&p_mutex);
        if (cond != msr_condition::no_task) {
            QMessageBox::information(this, "Calculations in progress...", "Wait until the calculations end");
            return;
        }
        pthread_mutex_lock(&p_mutex);
        msr_task.func_num = (msr_task.func_num + 1) % 8;
        msr_task.f = get_funk(msr_task.func_num);
        msr_task.condition = msr_condition::has_task;
        pthread_mutex_unlock(&p_mutex);
        pthread_cond_broadcast(&p_cond);

    } else if (event->key() == Qt::Key_1) {
        pthread_mutex_lock(&p_mutex);
        switch (plot_task.current_paint) {
            case what_to_paint::function:
                plot_task.current_paint = what_to_paint::approximation;
                break;
            case what_to_paint::approximation:
                plot_task.current_paint = what_to_paint::residual;
                break;
            case what_to_paint::residual:
                plot_task.current_paint = what_to_paint::function;
                break;
        }
        pthread_mutex_unlock(&p_mutex);
        updateMenuTitle();

    } else if (event->key() == Qt::Key_2) {
        pthread_mutex_lock(&p_mutex);
        double s_a = plot_task.s_a, s_b = plot_task.s_b;
        double s_c = plot_task.s_c, s_d = plot_task.s_d;
        double tmp = (3 * s_a + s_b) / 4;
        s_b = (3 * s_b + s_a) / 4;
        s_a = tmp;
        tmp = (3 * s_c + s_d) / 4;
        s_d = (3 * s_d + s_c) / 4;
        s_c = tmp;

        if (s_b - s_a < 1e-6 || s_d - s_c < 1e-6) {
            pthread_mutex_unlock(&p_mutex);
            return;
        }
        plot_task.s += 1;
        plot_task.s_a = s_a;
        plot_task.s_b = s_b;
        plot_task.s_c = s_c;
        plot_task.s_d = s_d;
        pthread_mutex_unlock(&p_mutex);
        updateMenuTitle();
    } else if (event->key() == Qt::Key_3) {
        pthread_mutex_lock(&p_mutex);
        if (plot_task.s == 0) {
            pthread_mutex_unlock(&p_mutex);
            return;
        }
        plot_task.s -= 1;
        if (plot_task.s == 0) {
            plot_task.s_a = plot_task.a;
            plot_task.s_b = plot_task.b;
            plot_task.s_c = plot_task.c;
            plot_task.s_d = plot_task.d;
        } else {
            double s_a = plot_task.s_a, s_b = plot_task.s_b;
            double s_c = plot_task.s_c, s_d = plot_task.s_d;
            double tmp = (3 * s_a - s_b) / 2;
            s_b = (3 * s_b - s_a) / 2;
            s_a = tmp;
            tmp = (3 * s_c - s_d) / 2;
            s_d = (3 * s_d - s_c) / 2;
            s_c = tmp;
            plot_task.s_a = s_a;
            plot_task.s_b = s_b;
            plot_task.s_c = s_c;
            plot_task.s_d = s_d;
        }
        pthread_mutex_unlock(&p_mutex);
        updateMenuTitle();
    } else if (event->key() == Qt::Key_4) {
        pthread_mutex_lock(&p_mutex);
        msr_condition cond = msr_task.condition;
        pthread_mutex_unlock(&p_mutex);
        if (cond != msr_condition::no_task) {
            QMessageBox::information(this, "Calculations in progress...", "Wait until the calculations end");
            return;
        }
        pthread_mutex_lock(&p_mutex);
        msr_task.nx = msr_task.nx * 2;
        msr_task.ny = msr_task.ny * 2;
        msr_task.condition = msr_condition::has_task;
        pthread_mutex_unlock(&p_mutex);
        pthread_cond_broadcast(&p_cond);
    } else if (event->key() == Qt::Key_5) {
        pthread_mutex_lock(&p_mutex);
        msr_condition cond = msr_task.condition;
        pthread_mutex_unlock(&p_mutex);
        if (cond != msr_condition::no_task) {
            QMessageBox::information(this, "Calculations in progress...", "Wait until the calculations end");
            return;
        }
        pthread_mutex_lock(&p_mutex);
        msr_task.nx = msr_task.nx / 2;
        if (msr_task.nx < 5) {
            msr_task.nx = 5;
        }
        msr_task.ny = msr_task.ny / 2;
        if (msr_task.ny < 5) {
            msr_task.ny = 5;
        }
        msr_task.condition = msr_condition::has_task;
        pthread_mutex_unlock(&p_mutex);
        pthread_cond_broadcast(&p_cond);
    } else if (event->key() == Qt::Key_6) {
        pthread_mutex_lock(&p_mutex);
        msr_condition cond = msr_task.condition;
        pthread_mutex_unlock(&p_mutex);
        if (cond != msr_condition::no_task) {
            QMessageBox::information(this, "Calculations in progress...", "Wait until the calculations end");
            return;
        }
        pthread_mutex_lock(&p_mutex);
        msr_task.inaccuracy += 1;
        msr_task.condition = msr_condition::has_task;
        pthread_mutex_unlock(&p_mutex);
        pthread_cond_broadcast(&p_cond);
    } else if (event->key() == Qt::Key_7) {
        pthread_mutex_lock(&p_mutex);
        msr_condition cond = msr_task.condition;
        pthread_mutex_unlock(&p_mutex);
        if (cond != msr_condition::no_task) {
            QMessageBox::information(this, "Calculations in progress...", "Wait until the calculations end");
            return;
        }
        pthread_mutex_lock(&p_mutex);
        msr_task.inaccuracy -= 1;
        msr_task.condition = msr_condition::has_task;
        pthread_mutex_unlock(&p_mutex);
        pthread_cond_broadcast(&p_cond);
    } else if (event->key() == Qt::Key_8) {
        pthread_mutex_lock(&p_mutex);
        plot_task.mx = plot_task.mx * 2;
        plot_task.my = plot_task.my * 2;
        pthread_mutex_unlock(&p_mutex);
        updateMenuTitle();
    } else if (event->key() == Qt::Key_9) {
        pthread_mutex_lock(&p_mutex);
        plot_task.mx = plot_task.mx / 2;
        if (plot_task.mx < 5) {
            plot_task.mx = 5;
        }
        plot_task.my = plot_task.my / 2;
        if (plot_task.my < 5) {
            plot_task.my = 5;
        }
        pthread_mutex_unlock(&p_mutex);
        updateMenuTitle();
    }
}

void MainWindow::get_rgb_color(double value, double max_value, double min_value, double& R, double& G, double& B) {
    static double R1 = 0, R2 = 255;
    static double G1 = 255, G2 = 0;
    static double B1 = 0, B2 = 0;
    double normal_mean = 0;
    if (max_value - min_value > 1e-16) {
        normal_mean = (value - min_value)/(max_value - min_value);
    }
    R = R1 * (1 - normal_mean) + R2 * normal_mean;
    G = G1 * (1 - normal_mean) + G2 * normal_mean;
    B = B1 * (1 - normal_mean) + B2 * normal_mean;
}

void MainWindow::get_max_min_value(data_to_plot& data) {
    double local_max = -1, local_min = std::numeric_limits<double>::max();
    double a = data.a, b = data.b;
    double c = data.c, d = data.d;
    double s_a = data.s_a, s_b = data.s_b;
    double s_c = data.s_c, s_d = data.s_d;
    int nx = data.nx, ny =  data.ny;
    int mx = data.mx, my = data.my;
    double (*f)(double, double) = data.f;
    double* x_approximation = data.x_approximation;
    double hmx = (s_b - s_a) / mx;
    double hmy = (s_d - s_c) / my;
    double hnx = (b - a) / nx;
    double hny = (d - c) / ny;
    auto current_paint = data.current_paint;
    if (current_paint == what_to_paint::function) {
        for(int j = 0; j < my; ++j) {
            for(int i = 0; i < mx; ++i) {
                double value_low_triangle =  f(s_a + hmx * (i + 2./3), s_c + hmy * (j + 1./3));
                double value_up_triangle = f(s_a + hmx * (i + 1./3), s_c + hmy * (j + 2./3));
                local_max = std::max({local_max, value_low_triangle, value_up_triangle});
                local_min = std::min({local_min, value_low_triangle, value_up_triangle});
            }
        }
        data.max_value = local_max;
        data.min_value = local_min;
        data.f_abs_max = std::max(std::fabs(local_max), std::fabs(local_min));
        return;
    }
    if (current_paint == what_to_paint::approximation) {
        for(int j = 0; j < my; ++j) {
            for(int i = 0; i < mx; ++i) {
                double value_low_triangle =  Pf(x_approximation, s_a + hmx * (i + 2./3), s_c + hmy * (j + 1./3), a, c, hnx, hny, nx, ny); 
                double value_up_triangle = Pf(x_approximation, s_a + hmx * (i + 1./3), s_c + hmy * (j + 2./3), a, c, hnx, hny, nx, ny);
                local_max = std::max({local_max, value_low_triangle, value_up_triangle});
                local_min = std::min({local_min, value_low_triangle, value_up_triangle});
            }
        }
        data.max_value = local_max;
        data.min_value = local_min;
        data.f_abs_max = std::max(std::fabs(local_min), std::fabs(local_max));
        return;
    }
    if (current_paint == what_to_paint::residual) {
        for(int j = 0; j < my; ++j) {
            for(int i = 0; i < mx; ++i) {
                double value_low_triangle = std::fabs(f(s_a + hmx * (i + 2./3), s_c + hmy * (j + 1./3))
                             - Pf(x_approximation, s_a + hmx * (i + 2./3), s_c + hmy * (j + 1./3), a, c, hnx, hny, nx, ny));
                double value_up_triangle = std::fabs(f(s_a + hmx * (i + 1./3), s_c + hmy * (j + 2./3))
                             - Pf(x_approximation, s_a + hmx * (i + 1./3), s_c + hmy * (j + 2./3), a, c, hnx, hny, nx, ny));
                local_max = std::max({local_max, value_low_triangle, value_up_triangle});
                local_min = std::min({local_min, value_low_triangle, value_up_triangle});
            }
        }
        data.max_value = local_max;
        data.min_value = local_min;
        data.f_abs_max = std::max(std::fabs(local_max), std::fabs(local_min));
        return;
    }
}

QPointF MainWindow::l2g (double x_loc, double y_loc, double x_min, double x_max, double y_min, double y_max) {
    int h = menuBar()->height();
    double x_gl = (x_loc - x_min) / (x_max - x_min) * width ();
    double y_gl = (y_max - y_loc) / (y_max - y_min) * (height() - h) + h;
    return QPointF(x_gl, y_gl);
}

void MainWindow::paint_graph(data_to_plot& data) {
    double max_value = data.max_value, min_value = data.min_value;
    double R = 0, G = 0, B = 0;
    QPainter painter(this);
    QBrush brush(Qt::SolidPattern);
    QPen pen(Qt::black, 0, Qt::SolidLine);
    QPointF triangle[3];
    double a = data.a, b = data.b;
    double c = data.c, d = data.d;
    double s_a = data.s_a, s_b = data.s_b;
    double s_c = data.s_c, s_d = data.s_d;
    int nx = data.nx, ny =  data.ny;
    int mx = data.mx, my = data.my;
    double (*f)(double, double) = data.f;
    double* x_approximation = data.x_approximation;
    double hmx = (s_b - s_a) / mx;
    double hmy = (s_d - s_c) / my;
    double hnx = (b - a) / nx;
    double hny = (d - c) / ny;
    auto current_paint = data.current_paint;
    if (current_paint == what_to_paint::function) {
        for(int j = 0; j < my; ++j) {
            for(int i = 0; i < mx; ++i) {
                double value_low_triangle =  f(s_a + hmx * (i + 2./3), s_c + hmy * (j + 1./3));
                get_rgb_color(value_low_triangle, max_value, min_value, R, G, B);
                brush.setColor(QColor(R, G, B));
                pen.setColor(QColor(R, G, B));
                painter.setPen(pen);
                painter.setBrush(brush);
                triangle[0] = l2g(s_a + hmx * i, s_c + hmy * j, s_a, s_b, s_c, s_d);
                triangle[1] = l2g(s_a + hmx * (i + 1), s_c + hmy * j, s_a, s_b, s_c, s_d);
                triangle[2] = l2g(s_a + hmx * (i + 1), s_c + hmy * (j + 1), s_a, s_b, s_c, s_d);
                painter.drawPolygon(triangle, 3);

                double value_up_triangle = f(s_a + hmx * (i + 1./3), s_c + hmy * (j + 2./3));
                get_rgb_color(value_up_triangle, max_value, min_value, R, G, B);
                brush.setColor(QColor(R, G, B));
                pen.setColor(QColor(R, G, B));
                painter.setPen(pen);
                painter.setBrush(brush);
                triangle[1] = l2g(s_a + hmx * (i), s_c + hmy * (j + 1), s_a, s_b, s_c, s_d);
                painter.drawPolygon(triangle, 3);
            }
        }
        return;
    }
    if (current_paint == what_to_paint::approximation) {
        for(int j = 0; j < my; ++j) {
            for(int i = 0; i < mx; ++i) {
                double value_low_triangle =  Pf(x_approximation, s_a + hmx * (i + 2./3), s_c + hmy * (j + 1./3), a, c, hnx, hny, nx, ny); 
                get_rgb_color(value_low_triangle, max_value, min_value, R, G, B);
                brush.setColor(QColor(R, G, B));
                pen.setColor(QColor(R, G, B));
                painter.setPen(pen);
                painter.setBrush(brush);
                triangle[0] = l2g(s_a + hmx * i, s_c + hmy * j, s_a, s_b, s_c, s_d);
                triangle[1] = l2g(s_a + hmx * (i + 1), s_c + hmy * j, s_a, s_b, s_c, s_d);
                triangle[2] = l2g(s_a + hmx * (i + 1), s_c + hmy * (j + 1), s_a, s_b, s_c, s_d);
                painter.drawPolygon(triangle, 3);

                double value_up_triangle = Pf(x_approximation, s_a + hmx * (i + 1./3), s_c + hmy * (j + 2./3), a, c, hnx, hny, nx, ny);
                get_rgb_color(value_up_triangle, max_value, min_value, R, G, B);
                brush.setColor(QColor(R, G, B));
                pen.setColor(QColor(R, G, B));
                painter.setPen(pen);
                painter.setBrush(brush);
                triangle[1] = l2g(s_a + hmx * (i), s_c + hmy * (j + 1), s_a, s_b, s_c, s_d);
                painter.drawPolygon(triangle, 3);
            }
        }
        return;
    }
    if (current_paint == what_to_paint::residual) {
        for(int j = 0; j < my; ++j) {
            for(int i = 0; i < mx; ++i) {
                double value_low_triangle = std::fabs(f(s_a + hmx * (i + 2./3), s_c + hmy * (j + 1./3))
                             - Pf(x_approximation, s_a + hmx * (i + 2./3), s_c + hmy * (j + 1./3), a, c, hnx, hny, nx, ny));
                get_rgb_color(value_low_triangle, max_value, min_value, R, G, B);
                brush.setColor(QColor(R, G, B));
                pen.setColor(QColor(R, G, B));
                painter.setPen(pen);
                painter.setBrush(brush);
                triangle[0] = l2g(s_a + hmx * i, s_c + hmy * j, s_a, s_b, s_c, s_d);
                triangle[1] = l2g(s_a + hmx * (i + 1), s_c + hmy * j, s_a, s_b, s_c, s_d);
                triangle[2] = l2g(s_a + hmx * (i + 1), s_c + hmy * (j + 1), s_a, s_b, s_c, s_d);
                painter.drawPolygon(triangle, 3);

                double value_up_triangle = std::fabs(f(s_a + hmx * (i + 1./3), s_c + hmy * (j + 2./3))
                             - Pf(x_approximation, s_a + hmx * (i + 1./3), s_c + hmy * (j + 2./3), a, c, hnx, hny, nx, ny));
                get_rgb_color(value_up_triangle, max_value, min_value, R, G, B);
                brush.setColor(QColor(R, G, B));
                pen.setColor(QColor(R, G, B));
                painter.setPen(pen);
                painter.setBrush(brush);
                triangle[1] = l2g(s_a + hmx * (i), s_c + hmy * (j + 1), s_a, s_b, s_c, s_d);
                painter.drawPolygon(triangle, 3);
            }
        }
        return;
    }

}

void MainWindow::paintEvent(QPaintEvent*) {
    if (plot_task.is_ready == true) {
        paint_graph(plot_task);
        pthread_mutex_lock(&p_mutex);
        printf("max{|F_max|, |F_min|} = %e\n", plot_task.f_abs_max);
        pthread_mutex_unlock(&p_mutex);
        //updateMenuTitle();
    }
}

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
