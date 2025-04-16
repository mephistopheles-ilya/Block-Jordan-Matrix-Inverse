#pragma once

#include <QMainWindow>
#include <QTimer>


class MainWindow;
double (*get_funk(int num))(double, double);



class MainWindow: public QMainWindow {
    Q_OBJECT

    public:

        MainWindow(double a, double b, double c, double d, int nx, int ny
                , int mx, int my, int k, double eps, int mi, int p);
        void create_threads();
        void wait_threads();
        QSize minimumSizeHint() const override;
        QSize sizeHint() const override;
        ~MainWindow();

    protected:

        void keyPressEvent(QKeyEvent* event) override;
        void paintEvent(QPaintEvent* ) override;
        void closeEvent(QCloseEvent* event ) override;

    private:
        enum class what_to_paint {
            function,
            approximation,
            residual,
            udefine
        };
        enum class msr_condition {
            no_task,
            has_task,
            task_is_ready,
            quit_app
        };

         struct data_to_plot {
            double a = 0, b = 0, c = 0, d = 0;
            int nx = 0, ny = 0, mx = 0, my = 0;
            int func_num = 0;
            double (*f)(double, double) = nullptr;
            double* x_approximation = nullptr;
            what_to_paint current_paint = what_to_paint::udefine;
            int s = 0;
            double s_a = 0, s_b = 0, s_c = 0, s_d = 0;
            int inaccuracy = 0;
            bool is_reday = false;
        };

        struct data_to_msr {
            int nx = 0, ny = 0;
            int func_num = 0;
            double (*f)(double, double) = nullptr;
            double inaccuracy = 0;
            double* x_approximation = nullptr;
            msr_condition condition = msr_condition::no_task; 
        };

        struct Arg {
            int thr_num = 0;
            double eps = 0;
            int maxit = 0;
            int p = 0;
            double a = 0, b = 0, c = 0, d = 0;
            pthread_t tid = 0;
            data_to_msr* task = nullptr;
            pthread_mutex_t* gui_mutex = nullptr;
            pthread_cond_t* gui_cond = nullptr;
        };


        QPointF l2g (double x_loc, double y_loc, double x_min, double x_max, double y_min, double y_max);
        static void* create_msr_approximation(void* argument);
        void time_out_checker();
        void paint_graph(const data_to_plot& data);
        void get_rgb_color(double value, double max_value, double min_value, double& R, double& G, double& B);
        void get_max_min_value(const data_to_plot& data, double& max_value, double& min_value);

    private:

        Arg* args = nullptr;

        data_to_plot plot_task;
        data_to_msr msr_task;


        QTimer* timer = nullptr;
        pthread_mutex_t p_mutex = PTHREAD_MUTEX_INITIALIZER;
        pthread_cond_t p_cond = PTHREAD_COND_INITIALIZER;

};
