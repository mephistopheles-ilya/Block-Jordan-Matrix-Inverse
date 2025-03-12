#include <QMenuBar>
#include <QKeyEvent>
#include <QPainter>
#include <QLabel>
#include <QHBoxLayout>
#include <cmath>
#include <algorithm>

#include "window.hpp"
#include "approximations.hpp"

#define L2G(X,Y) (l2g ((X), (Y), min_y, max_y))

inline double f0(double ) {
    return 1;
}
inline double f1(double x) {
    return x;
}
inline double f2(double x) {
    return x * x;
}
inline double f3(double x) {
    return x * x * x;
}
inline double f4(double x) {
    return x * x * x * x;
}
inline double f5(double x) {
    return std::exp(x);
}
inline double f6(double x) {
    return 1./(25 * x * x + 1);
}


MainWindow::MainWindow(double a, double b, int n, int k, QWidget* parent):
    QMainWindow(parent), a(a), b(b), n(n), k(k) {

        newton_coeff = new double[50];

        points_x = new double[n];
        f_x = new double[n];
        d = new double[n];
        c = new double[4 * (n - 1)];

        update_function();
        fill_points_array();
        fill_func_array();
        add_error();

        allocated_n = n;

        plot_number = 1;
        is_plot_func = true;
        is_plot_apr1 = true;

        if (n <= 50) {
            Newton_approximation(points_x, f_x, newton_coeff, n);
        }

        setWindowTitle("Function Approximation");

        QMenuBar *menuBar = new QMenuBar(this);
        setMenuBar(menuBar);

        functionsMenu = menuBar->addMenu("Functions");

        QAction *func0 = functionsMenu->addAction(QIcon::fromTheme("k0"), "k = 0 f(x) = 1");
        QAction *func1 = functionsMenu->addAction(QIcon::fromTheme("k1"), "k = 1 f(x) = x");
        QAction *func2 = functionsMenu->addAction(QIcon::fromTheme("k2"), "k = 2 f(x) = x^2");
        QAction *func3 = functionsMenu->addAction(QIcon::fromTheme("k3"), "k = 3 f(x) = x^3");
        QAction *func4 = functionsMenu->addAction(QIcon::fromTheme("k4"), "k = 4 f(x) = x^4");
        QAction *func5 = functionsMenu->addAction(QIcon::fromTheme("k5"), "k = 5 f(x) = exp(x)");
        QAction *func6 = functionsMenu->addAction(QIcon::fromTheme("k6"), "k = 6 f(x) = 1/(25*x^2 + 1)");

        functionsMenu->addSeparator();

        connect(func0, &QAction::triggered, [this]() {set_func_number(0);});
        connect(func1, &QAction::triggered, [this]() {set_func_number(1);});
        connect(func2, &QAction::triggered, [this]() {set_func_number(2);});
        connect(func3, &QAction::triggered, [this]() {set_func_number(3);});
        connect(func4, &QAction::triggered, [this]() {set_func_number(4);});
        connect(func5, &QAction::triggered, [this]() {set_func_number(5);});
        connect(func6, &QAction::triggered, [this]() {set_func_number(6);});

        updateMenuTitle();
    }

void MainWindow::updateMenuTitle() {
    QString title;
    switch (k) {
        case 0: 
            title = QString("k = 0 f(x) = 1, n = %1, s = %2, [a, b] = [%3, %4], max{|F_min|, |F_max|} = %5, p = %6")
                .arg(n).arg(s).arg(a).arg(b).arg(F_abs_max).arg(p);
            break;
        case 1: 
            title = QString("k = 1 f(x) = x, n = %1, s = %2, [a, b] = [%3, %4], max{|F_min|, |F_max|} = %5, p = %6")
                .arg(n).arg(s).arg(a).arg(b).arg(F_abs_max).arg(p);
            break;
        case 2: 
            title = QString("k = 2 f(x) = x^2, n = %1, s = %2, [a, b] = [%3, %4], max{|F_min|, |F_max|} = %5, p = %6")
                .arg(n).arg(s).arg(a).arg(b).arg(F_abs_max).arg(p);
            break;
        case 3: 
            title = QString("k = 3 f(x) = x^3, n = %1, s = %2, [a, b] = [%3, %4], max{|F_min|, |F_max|} = %5, p = %6")
                .arg(n).arg(s).arg(a).arg(b).arg(F_abs_max).arg(p);
            break;
        case 4: 
            title = QString("k = 4 f(x) = x^4, n = %1, s = %2, [a, b] = [%3, %4], max{|F_min|, |F_max|} = %5, p = %6")
                .arg(n).arg(s).arg(a).arg(b).arg(F_abs_max).arg(p);
            break;
        case 5: 
            title = QString("k = 5 f(x) = exp(x), n = %1, s = %2, [a, b] = [%3, %4], max{|F_min|, |F_max|} = %5, p = %6")
                .arg(n).arg(s).arg(a).arg(b).arg(F_abs_max).arg(p);
            break;
        case 6: 
            title = QString("k = 6 f(x) = 1/(25*x^2 + 1), n = %1, s = %2, [a, b] = [%3, %4], max{|F_min|, |F_max|} = %5, p = %6")
                .arg(n).arg(s).arg(a).arg(b).arg(F_abs_max).arg(p);
            break;
    }
    functionsMenu->setTitle(title);
}

void MainWindow::update_function() {
    switch (k) {
        case 0:
            f = f0;
            break;
        case 1:
            f = f1;
            break;
        case 2:
            f = f2;
            break;
        case 3:
            f = f3;
            break;
        case 4:
            f = f4;
            break;
        case 5:
            f = f5;
            break;
        case 6:
            f = f6;
    }
}

void MainWindow::next_function() {
    k = (k + 1) % 7;
}

QSize MainWindow::minimumSizeHint() const {
    return QSize(100, 100);
}

QSize MainWindow::sizeHint() const {
    return QSize(1000, 1000);
}

void MainWindow::fill_func_array() {
    for(int i = 0; i < n; ++i) {
        f_x[i] = f(points_x[i]);
    }
}

void MainWindow::add_error() {
    double x1 = 0, y1 = 0;
    double delta_x = (b - a)/width();
    double max_f = 0;
    for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
        y1 = std::fabs(f(x1));
        if (y1 > max_f)
            max_f = y1;
    }
    f_x[n/2] = f(points_x[n/2]) + p * 0.1 * max_f;

}

void MainWindow::delete_allocate() {
    delete[] points_x;
    delete[] f_x;
    delete[] d;
    delete[] c;

    points_x = new double[n];
    f_x = new double[n];
    d = new double[n];
    c = new double[4 * (n - 1)];
}

void MainWindow::fill_points_array() {
    double delta_x = (b - a) / (n - 1);
    for(int i = 0; i < n; ++i) {
        points_x[i] = a + i * delta_x;
    }
}

QPointF MainWindow::l2g (double x_loc, double y_loc, double y_min, double y_max) {
    int h = menuBar()->height();
    double x_gl = (x_loc - a) / (b - a) * width ();
    double y_gl = (y_max - y_loc) / (y_max - y_min) * (height() - h) + h;
    return QPointF (x_gl, y_gl);
}

MainWindow::~MainWindow() {
    delete[] points_x;
    delete[] f_x;
    delete[] newton_coeff;
    delete[] d;
    delete[] c;
}

void MainWindow::keyPressEvent(QKeyEvent *event) {
    if (event->key() == Qt::Key_0) {
        next_function();
        update_function();
        fill_func_array();
        add_error();

        if ((is_plot_apr1 == true || is_plot_resid == true) && n <= 50) {
            Newton_approximation(points_x, f_x, newton_coeff, n);
        }
        if (is_plot_apr2 == true || is_plot_resid == true) {
            cube_approximation(points_x, f_x, d, c, n);
        }
        update();
    } else if (event->key() == Qt::Key_4) {
        n *= 2;
        if (allocated_n < n) {
            delete_allocate();
        }
        fill_points_array();
        fill_func_array();
        add_error();
        if ((is_plot_apr1 == true || is_plot_resid == true) && n <= 50) {
            Newton_approximation(points_x, f_x, newton_coeff, n);
        }
        if (is_plot_apr2 == true || is_plot_resid == true) {
            cube_approximation(points_x, f_x, d, c, n);
        }

        update();
    } else if (event->key() == Qt::Key_5) {
        if (n >= 10) {
            n /= 2;
        } else {
            n = 5;
        }
        fill_points_array();
        fill_func_array();
        add_error();
        if ((is_plot_apr1 == true || is_plot_resid == true) && n <= 50) {
            Newton_approximation(points_x, f_x, newton_coeff, n);
        }
        if (is_plot_apr2 == true || is_plot_resid == true) {
            cube_approximation(points_x, f_x, d, c, n);
        }

        update();
    } else if (event->key() == Qt::Key_2) {
        ++s;
        if (std::fabs(a) >= std::numeric_limits<double>::max()/2 
                || std::fabs(b) >= std::numeric_limits<double>::max()/2) {
            --s;
            return;
        }
        double tmp = (3 * a - b) / 2;
        b = (3 * b - a) / 2;
        a = tmp;
        fill_points_array();
        fill_func_array();
        add_error();
        if ((is_plot_apr1 == true || is_plot_resid == true) && n <= 50) {
            Newton_approximation(points_x, f_x, newton_coeff, n);
        }
        if (is_plot_apr2 == true || is_plot_resid == true) {
            cube_approximation(points_x, f_x, d, c, n);
        }

        update();
    } else if (event->key() == Qt::Key_3) {
        double tmp_a = a, tmp_b = b;
        --s;
        double tmp = (3 * a + b) / 4;
        b = (3 * b + a) / 4;
        a = tmp;
        if (a >= b) {
            a = tmp_a;
            b = tmp_b;
            ++s;
            return;
        }
        fill_points_array();
        fill_func_array();
        add_error();
        if ((is_plot_apr1 == true || is_plot_resid == true) && n <= 50) {
            Newton_approximation(points_x, f_x, newton_coeff, n);
        }
        if (is_plot_apr2 == true || is_plot_resid == true) {
            cube_approximation(points_x, f_x, d, c, n);
        }

        update();
    } else if (event->key() == Qt::Key_6) {
        ++p;
        add_error();
        if ((is_plot_apr1 == true || is_plot_resid == true) && n <= 50) {
            Newton_approximation(points_x, f_x, newton_coeff, n);
        }
        if (is_plot_apr2 == true || is_plot_resid == true) {
            cube_approximation(points_x, f_x, d, c, n);
        }

        update();
    } else if (event->key() == Qt::Key_7) {
        --p;
        add_error();
        if ((is_plot_apr1 == true || is_plot_resid == true) && n <= 50) {
            Newton_approximation(points_x, f_x, newton_coeff, n);
        }
        if (is_plot_apr2 == true || is_plot_resid == true) {
            cube_approximation(points_x, f_x, d, c, n);
        }

        update();
    } else if (event->key() == Qt::Key_1) {
        ++plot_number;
        if (plot_number > 4) {
            plot_number = plot_number % 4;
        }
        printf("%d\n", plot_number);
        if (plot_number == 1 || plot_number == 2 || plot_number == 3) {
            is_plot_func = true;
        } else {
            is_plot_func = false;
        }
        if ((plot_number == 1 || plot_number == 3) && n <= 50) {
            is_plot_apr1 = true;
        } else {
            is_plot_apr1 = false;
        }
        if (plot_number == 2 || plot_number == 3) {
            is_plot_apr2 = true;
        } else {
            is_plot_apr2 = false;
        }
        if (plot_number == 4) {
            is_plot_resid = true;
        } else {
            is_plot_resid = false;
        }
        if ((is_plot_apr1 == true || is_plot_resid == true) && n <= 50) {
            Newton_approximation(points_x, f_x, newton_coeff, n);
        }
        if (is_plot_apr2 == true || is_plot_resid == true) {
            cube_approximation(points_x, f_x, d, c, n);
        }
        update();
    }
}

void MainWindow::paintEvent(QPaintEvent* ) {
    QPainter painter (this);
    QPen pen_blue(Qt::blue, 2, Qt::SolidLine); 
    QPen pen_green(Qt::green, 2, Qt::SolidLine); 
    QPen pen_red(Qt::red, 2, Qt::SolidLine); 
    QPen pen_black(Qt::black, 2, Qt::SolidLine);


    double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
    double delta_y = 0, delta_x = (b - a)/width();
    double min_y = 0, max_y = 0;

    if (is_plot_func == true) {
        for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
            y1 = f(x1);
            if (y1 < min_y)
                min_y = y1;
            if (y1 > max_y)
                max_y = y1;
        }
    }
    if (is_plot_apr1 == true && n <= 50) {
        for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
            y1 = Newton_calck(newton_coeff, points_x, n, x1);
            if (y1 < min_y)
                min_y = y1;
            if (y1 > max_y)
                max_y = y1;
        }
    }
    if (is_plot_apr2 == true) {
        for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
            y1 = cube_calc(points_x, c, a, b, n, x1);
            if (y1 < min_y)
                min_y = y1;
            if (y1 > max_y)
                max_y = y1;
        }
    }
    if (is_plot_resid == true) {
        for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
            y1 = std::fabs(f(x1) - cube_calc(points_x, c, a, b, n, x1));
            if (y1 < min_y)
                min_y = y1;
            if (y1 > max_y)
                max_y = y1;
        }
        if (n <= 50) {
            for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
                y1 = std::fabs(f(x1) - Newton_calck(newton_coeff, points_x, n, x1));
                if (y1 < min_y)
                    min_y = y1;
                if (y1 > max_y)
                    max_y = y1;
            }
        }
    }

    F_abs_max = std::max(std::fabs(min_y), std::fabs(max_y));

    delta_y = 0.01 * (max_y - min_y);
    min_y -= delta_y;
    max_y += delta_y;

    painter.setPen(pen_black);
    painter.drawLine(L2G(a, 0), L2G(b, 0));
    painter.drawLine(L2G(0, min_y), L2G(0, max_y));

    if (is_plot_func == true) {
        painter.setPen(pen_blue);
        x1 = a;
        y1 = f(x1);
        for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
            y2 = f(x2);
            painter.drawLine(L2G(x1, y1), L2G(x2, y2));
            x1 = x2;
            y1 = y2;
        }
        x2 = b;
        y2 = f(x2);
        painter.drawLine(L2G(x1, y1), L2G(x2, y2));
    }

    if (is_plot_apr1 == true && n <= 50) {
        painter.setPen(pen_green);
        x1 = a;
        y1 = Newton_calck(newton_coeff, points_x, n, x1);
        for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
            y2 = Newton_calck(newton_coeff, points_x, n, x1);
            painter.drawLine(L2G(x1, y1), L2G(x2, y2));
            x1 = x2;
            y1 = y2;
        }
        x2 = b;
        y2 = Newton_calck(newton_coeff, points_x, n, x1);
        painter.drawLine(L2G(x1, y1), L2G(x2, y2));
    }
    if (is_plot_apr2 == true) {
        painter.setPen(pen_red);
        x1 = a;
        y1 = cube_calc(points_x, c, a, b, n, x1);
        for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
            y2 = cube_calc(points_x, c, a, b, n, x1);
            painter.drawLine(L2G(x1, y1), L2G(x2, y2));
            x1 = x2;
            y1 = y2;
        }
        x2 = b;
        y2 = cube_calc(points_x, c, a, b, n, x1);
        painter.drawLine(L2G(x1, y1), L2G(x2, y2));
    }
    if (is_plot_resid == true) {
        painter.setPen(pen_red);
        x1 = a;
        y1 = std::fabs(f(x1) - cube_calc(points_x, c, a, b, n, x1));
        for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
            y2 = std::fabs(f(x1) - cube_calc(points_x, c, a, b, n, x1));
            painter.drawLine(L2G(x1, y1), L2G(x2, y2));
            x1 = x2;
            y1 = y2;
        }
        x2 = b;
        y2 = std::fabs(f(x1) - cube_calc(points_x, c, a, b, n, x1));
        painter.drawLine(L2G(x1, y1), L2G(x2, y2));
        if (n <= 50) {
            painter.setPen(pen_green);
            x1 = a;
            y1 = std::fabs(f(x1) - Newton_calck(newton_coeff, points_x, n, x1));
            for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
                y2 = std::fabs(f(x1) - Newton_calck(newton_coeff, points_x, n, x1));
                painter.drawLine(L2G(x1, y1), L2G(x2, y2));
                x1 = x2;
                y1 = y2;
            }
            x2 = b;
            y2 = std::fabs(f(x1) - Newton_calck(newton_coeff, points_x, n, x1));
            painter.drawLine(L2G(x1, y1), L2G(x2, y2));
        }
    }
    updateMenuTitle();
}


void MainWindow::set_func_number(int new_k) {
    k = new_k;
    update_function();
    fill_func_array();
    if (is_plot_apr1 == true || (is_plot_resid == true && n <= 50)) {
        Newton_approximation(points_x, f_x, newton_coeff, n);
    }
    if (is_plot_apr2 == true || is_plot_resid == true) {
        cube_approximation(points_x, f_x, d, c, n);
    }
    update();
}



