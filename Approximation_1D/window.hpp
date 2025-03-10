#include <QApplication>
#include <QMainWindow>
#include <QMenuBar>
#include <QKeyEvent>
#include <QPainter>
#include <cmath>
#include <algorithm>

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


class MainWindow: public QMainWindow {
    Q_OBJECT

    private:

        QMenu *functionsMenu = nullptr;
        double (*f)(double);
        double a = 0, b = 0;
        int n = 0, k = 0;
        int s = 0, p = 0;
        double F_max = 0, F_min = 0;
        double* points_x = nullptr;
        double* f_x = nullptr;
        double* newton_coeff = nullptr;
        double* d = nullptr;
        double* c = nullptr;
        double max_f = 0;

    public:

        MainWindow(double a, double b, int n, int k, QWidget* parent = nullptr):
            QMainWindow(parent), a(a), b(b), n(n), k(k) {

                setWindowTitle("Math Function Plotter");

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
                updateFunction();

            }

        QSize minimumSizeHint() const override {
            return QSize(100, 100);
        }

        QSize sizeHint() const override {
            return QSize(1000, 1000);
        }

    protected:

        void memmory_refilling() {
            delete[] points_x;
            delete[] f_x;
            delete[] newton_coeff;
            delete[] d;
            delete[] c;

            points_x = new double[n];
            f_x = new double[n];
            newton_coeff = new double[n];
            d = new double[n];
            c = new double[4 * (n - 1)];

            double delta_x = (b - a) / (n - 1);
            for(int i = 0; i < n; ++i) {
                points_x[i] = a + i * delta_x;
            }

            for(int i = 0; i < n; ++i) {
                f_x[i] = f(points_x[i]);
            }
            f_x[n/2] += p * 0.1 * max_f;
        }

        void Newton_approximation() {
            for(int i = 0; i < n; ++i) {
                newton_coeff[i] = f_x[i];
            }

            for(int j = 1; j < n; ++j) {
                for(int i = n - 1; i >= j; --i) {
                    double delta_f = (newton_coeff[i] - newton_coeff[i - 1]) / (points_x[i] - points_x[i - j]);
                    newton_coeff[i] = delta_f;
                }
            }
        }

        double Newton_calck(double x) {
            double res = newton_coeff[n - 1];
            for(int i = n - 2; i >= 0; --i) {
                res *= (x - points_x[i]);
                res += newton_coeff[i];
            }
            return res;
        }

        int sign(double x) {
            if (x > 0) return 1;
            if (x < 0) return -1;
            return 0;
        }

        void cube_three_points(double x, double y, double* b, double* res) {
            double xx = x * x, xxx = x * x * x, yy = y * y, yyy = y * y * y;
            double den = xx * yyy - xxx * yy;
            double f = b[0], q = b[1], p = b[2], k = b[3];
            res[0] = p;
            res[1] = k;
            double tmp = f * yyy + k * (xxx * y - x * yyy) + p * (xxx - yyy) - q * xxx;
            res[2] = tmp/den;
            tmp = -f * yy + k * (x * yy - xx * y) + p * (yy - xx) + q * xx;
            res[3] = tmp/den;
        }

        void cube_approximation() {
            for(int i = 1; i < n - 1; ++i) {
                double df1 = (f_x[i] - f_x[i - 1])/(points_x[i] - points_x[i - 1]);
                double df2 = (f_x[i + 1] - f_x[i])/(points_x[i + 1] - points_x[i]);
                if (df1 * df2 > 0) {
                    d[i] = sign(df1) * std::min(std::fabs(df1), std::fabs(df2));
                } else {
                    d[i] = 0;
                }
            }

            for(int i = 2; i < n - 3; ++i) {
                c[4 * i + 0] = f_x[i];
                c[4 * i + 1] = d[i];
                double delta_x = points_x[i + 1] - points_x[i];
                double df = (f_x[i + 1] - f_x[i])/delta_x;
                c[4 * i + 2] = (3 * df - 2 * d[i] - d[i + 1])/delta_x;
                c[4 * i + 3] = (d[i] + d[i + 1] - 2 * df)/(delta_x * delta_x);
            }
            double b1[] = {f_x[0], f_x[1], f_x[2], d[2]};
            cube_three_points(points_x[0] - points_x[2], points_x[1] - points_x[2], b1, c + 0);
            cube_three_points(points_x[0] - points_x[2], points_x[1] - points_x[2], b1, c + 4);

            double b2[] = {f_x[n - 1], f_x[n - 2], f_x[n - 3], d[n - 3]};
            cube_three_points(points_x[n - 1] - points_x[n - 3], points_x[n - 2] - points_x[n - 3], b2, c + 4 * (n - 1) - 4);
            cube_three_points(points_x[n - 1] - points_x[n - 3], points_x[n - 2] - points_x[n - 3], b2, c + 4 * (n - 1) - 8);
        }

        double cube_calc(double x) {
            int i = (x - a) * (n - 1)/ (b - a);
            double x_i = 0;
            x_i = points_x[i];
            if (x <= points_x[2]) {
                x_i = points_x[2];
            }
            if (x >= points_x[n - 3]) {
                x_i = points_x[n - 3];
                if (i == n - 1) {
                    i = n - 2;
                }
            }
            return c[4 * i + 0] + c[4 * i + 1] * (x - x_i) + c[4 * i + 2] * (x - x_i) * (x - x_i) + c[4 * i + 3] * (x - x_i) * (x - x_i) * (x - x_i);
        }




        void keyPressEvent(QKeyEvent *event) override {
            if (event->key() == Qt::Key_0) {
                next_function();
                updateMenuTitle();
                updateFunction();
                update();
            } else if (event->key() == Qt::Key_4) {
                n *= 2;
                updateMenuTitle();
                update();
            } else if (event->key() == Qt::Key_5) {
                if (n >= 8) {
                    n /= 2;
                } else {
                    n = 3;
                }
                updateMenuTitle();
                update();
            } else if (event->key() == Qt::Key_2) {
                ++s;
                double tmp = (3 * a - b) / 2;
                b = (3 * b - a) / 2;
                a = tmp;
                updateMenuTitle();
                update();
            } else if (event->key() == Qt::Key_3) {
                --s;
                double tmp = (3 * a + b) / 4;
                b = (3 * b + a) / 4;
                a = tmp;
                updateMenuTitle();
                update();
            } else if (event->key() == Qt::Key_6) {
                ++p;
                updateMenuTitle();
                update();
            } else if (event->key() == Qt::Key_7) {
                --p;
                updateMenuTitle();
                update();
            }
        }
        QPointF l2g (double x_loc, double y_loc, double y_min, double y_max) {
            int h = menuBar()->height();
            double x_gl = (x_loc - a) / (b - a) * width ();
            double y_gl = (y_max - y_loc) / (y_max - y_min) * (height() - h) + h;
            return QPointF (x_gl, y_gl);
        }

        void clalc_max_f() {
            double x1 = 0, y1 = 0;
            double delta_x = (b - a)/width();
            max_f = 0;
            for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
                y1 = std::fabs(f(x1));
                if (y1 > max_f)
                    max_f = y1;
            }
        }


        void paintEvent(QPaintEvent* ) override {
            QPainter painter (this);
            QPen pen_blue(Qt::blue, 4, Qt::SolidLine); 
            QPen pen_green(Qt::green, 4, Qt::SolidLine); 
            QPen pen_red(Qt::red, 4, Qt::SolidLine); 
            QPen pen_black(Qt::black, 4, Qt::SolidLine);


            clalc_max_f();
            memmory_refilling();
            //Newton_approximation();
            cube_approximation();

            double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
            double delta_y = 0, delta_x = (b - a)/width();
            double min_y = 0, max_y = 0;
            for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
                y1 = f(x1);
                if (y1 < min_y)
                    min_y = y1;
                if (y1 > max_y)
                    max_y = y1;

            }

#if 0
            for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
                y1 = Newton_calck(x1);
                if (y1 < min_y)
                    min_y = y1;
                if (y1 > max_y)
                    max_y = y1;
            }

#endif
            for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
                y1 = cube_calc(x1);
                if (y1 < min_y)
                    min_y = y1;
                if (y1 > max_y)
                    max_y = y1;
            }
            for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
                y1 = std::fabs(f(x1) - cube_calc(x1));
                if (y1 < min_y)
                    min_y = y1;
                if (y1 > max_y)
                    max_y = y1;
            }

            delta_y = 0.01 * (max_y - min_y);
            min_y -= delta_y;
            max_y += delta_y;

            painter.setPen(pen_black);
            painter.drawLine(L2G(a, 0), L2G(b, 0));
            painter.drawLine(L2G(0, min_y), L2G(0, max_y));

#if 0
            painter.setPen(pen_blue);
            x1 = a;
            y1 = std::fabs(f(x1) - cube_calc(x1));
            for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
                y2 = std::fabs(f(x2) - cube_calc(x2));
                painter.drawLine(L2G(x1, y1), L2G(x2, y2));
                x1 = x2;
                y1 = y2;
            }
            x2 = b;
            y2 = std::fabs(f(x2) - cube_calc(x2));
            painter.drawLine(L2G(x1, y1), L2G(x2, y2));
#endif
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
#if 0
            painter.drawLine(L2G(x1, y1), L2G(x2, y2));
            painter.setPen(pen_green);
            x1 = a;
            y1 = Newton_calck(x1);
            for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
                y2 = Newton_calck(x2);
                painter.drawLine(L2G(x1, y1), L2G(x2, y2));
                x1 = x2;
                y1 = y2;
            }
            x2 = b;
            y2 = Newton_calck(x2);
            painter.drawLine(L2G(x1, y1), L2G(x2, y2));
#endif
            painter.setPen(pen_red);
            x1 = a;
            y1 = cube_calc(x1);
            for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
                y2 = cube_calc(x2);
                painter.drawLine(L2G(x1, y1), L2G(x2, y2));
                x1 = x2;
                y1 = y2;
            }
            x2 = b;
            y2 = cube_calc(x2);
            painter.drawLine(L2G(x1, y1), L2G(x2, y2));

            F_min = min_y;
            F_max = max_y;
            updateMenuTitle();

        }

        void set_func_number(int new_k) {
            k = new_k;
            updateMenuTitle();
            updateFunction();
            update();
        }
        void next_function() {
            k = (k + 1) % 7;
        }

        void updateMenuTitle() {
            QString title;
            switch (k) {
                case 0: 
                    title = QString("k = 0 f(x) = 1, n = %1, s = %2, [a, b] = [%3, %4], [F_min, F_max] = [%5, %6], p = %7")
                        .arg(n).arg(s).arg(a).arg(b).arg(F_min).arg(F_max).arg(p);
                    break;
                case 1: 
                    title = QString("k = 1 f(x) = x, n = %1, s = %2, [a, b] = [%3, %4], [F_min, F_max] = [%5, %6], p = %7")
                        .arg(n).arg(s).arg(a).arg(b).arg(F_min).arg(F_max).arg(p);
                    break;
                case 2: 
                    title = QString("k = 2 f(x) = x^2, n = %1, s = %2, [a, b] = [%3, %4], [F_min, F_max] = [%5, %6], p = %7")
                        .arg(n).arg(s).arg(a).arg(b).arg(F_min).arg(F_max).arg(p);
                    break;
                case 3: 
                    title = QString("k = 3 f(x) = x^3, n = %1, s = %2, [a, b] = [%3, %4], [F_min, F_max] = [%5, %6], p = %7")
                        .arg(n).arg(s).arg(a).arg(b).arg(F_min).arg(F_max).arg(p);
                    break;
                case 4: 
                    title = QString("k = 4 f(x) = x^4, n = %1, s = %2, [a, b] = [%3, %4], [F_min, F_max] = [%5, %6], p = %7")
                        .arg(n).arg(s).arg(a).arg(b).arg(F_min).arg(F_max).arg(p);
                    break;
                case 5: 
                    title = QString("k = 5 f(x) = exp(x), n = %1, s = %2, [a, b] = [%3, %4], [F_min, F_max] = [%5, %6], p = %7")
                        .arg(n).arg(s).arg(a).arg(b).arg(F_min).arg(F_max).arg(p);
                    break;
                case 6: 
                    title = QString("k = 6 f(x) = 1/(25*x^2 + 1), n = %1, s = %2, [a, b] = [%3, %4], [F_min, F_max] = [%5, %6], p = %7")
                        .arg(n).arg(s).arg(a).arg(b).arg(F_min).arg(F_max).arg(p);
                    break;
            }
            functionsMenu->setTitle(title);
        }
        void updateFunction() {
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

    public:

        ~MainWindow() {
            delete[] points_x;
            delete[] f_x;
            delete[] newton_coeff;
            delete[] d;
            delete[] c;
        };
};

