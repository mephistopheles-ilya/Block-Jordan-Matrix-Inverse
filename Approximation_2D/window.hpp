#pragma once

#include <QMainWindow>
#include <QTimer>

#include "thread.hpp"

class MainWindow;

struct Arg {
    double a = -1., b = -1., c = -1., d = -1.;
    int nx = -1, ny = -1;
    int k = -1;
    double eps = -1;
    int mi = -1;
    int p = -1;
    int thr_num = -1;
    
    pthread_t tid = 0;

    double r1 = -1., r2 = -1., r3 = -1., r4 = -1.;
    double t1 = -1., t2 = -1.;
    int it = -1;

    double* A = nullptr;
    double* B = nullptr;
    int* I = nullptr;
    double* x = nullptr;
    double* r = nullptr;
    double* u =  nullptr;
    double* v = nullptr;

    MainWindow* window = nullptr;

    void get() const;

};


class MainWindow: public QMainWindow {
    Q_OBJECT

    public:

        MainWindow(double a, double b, double c, double d, int nx, int ny
                , int mx, int my, int k, double eps, int mi, int p);
        QSize minimumSizeHint() const override;
        QSize sizeHint() const override;
        ~MainWindow();

    protected:

        void keyPressEvent(QKeyEvent* event) override;
        void paintEvent(QPaintEvent* ) override;
        void closeEvent(QCloseEvent* event ) override;

    private:
 
        void allocate_memmory(int len_msr, int len_diag);
        void deallocate_memmory();
        bool check_nullptr() const;
        void create_threads();
        static void* create_msr_approximation(void* argument);
        void time_out_checker();

    private:

        Arg* args = nullptr;
        double* A = nullptr;
        double* B = nullptr;
        int* I = nullptr;
        double* x_gui = nullptr;
        double* x_msr = nullptr;
        double* r = nullptr;
        double* u = nullptr;
        double* v = nullptr;

        double a = -1, b = -1, c = -1, d = -1;
        int nx = 0, ny = 0, mx = 0, my = 0;
        int k = 0;
        double eps = -1;
        int mi = -1;
        int p = -1;

        volatile bool msr_continue = true;
        bool calculations_done = false;

        QTimer* timer = nullptr;
        pthread_mutex_t p_mutex = PTHREAD_MUTEX_INITIALIZER;
        pthread_cond_t p_cond = PTHREAD_COND_INITIALIZER;

};
