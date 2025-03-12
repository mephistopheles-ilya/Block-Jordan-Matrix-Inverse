#pragma once

#include <QMainWindow>


class MainWindow: public QMainWindow {
    Q_OBJECT

    private:

        QMenu *functionsMenu = nullptr;
        
        double (*f)(double);

        double a = 0, b = 0;
        int n = 0, k = 0;
        int s = 0, p = 0;

        double F_abs_max = 0;

        double* points_x = nullptr;
        double* f_x = nullptr;
        double* newton_coeff = nullptr;
        double* d = nullptr;
        double* c = nullptr;

        int allocated_n = 0;
        int plot_number = 0;
        bool is_plot_func = false;
        bool is_plot_apr1 = false;
        bool is_plot_apr2 = false;
        bool is_plot_resid = false;

    public:

        MainWindow(double a, double b, int n, int k, QWidget* parent = nullptr);

        QSize minimumSizeHint() const override;

        QSize sizeHint() const override;

    protected:

        void fill_func_array();
        void fill_points_array();
        void delete_allocate();
        void add_error();

        void keyPressEvent(QKeyEvent *event) override;
        QPointF l2g (double x_loc, double y_loc, double y_min, double y_max);


        void paintEvent(QPaintEvent* ) override;

        void set_func_number(int new_k);

        void next_function();
        void updateMenuTitle();
        void update_function();

    public:

        ~MainWindow();
};

