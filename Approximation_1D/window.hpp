#include <QApplication>
#include <QMainWindow>
#include <QMenuBar>
#include <QKeyEvent>
#include <QPainter>
#include <cmath>

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
        int s = 0;

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
                if (n > 1) {
                    n /= 2;
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
            }
        }
        QPointF l2g (double x_loc, double y_loc, double y_min, double y_max) {
            int h = menuBar()->height();
            double x_gl = (x_loc - a) / (b - a) * width ();
            double y_gl = (y_max - y_loc) / (y_max - y_min) * (height() - h) + h;
            return QPointF (x_gl, y_gl);
        }

        void paintEvent(QPaintEvent* ) override {
            QPainter painter (this);
            QPen pen_blue(Qt::blue, 0, Qt::SolidLine); 
            QPen pen_black(Qt::black, 0, Qt::SolidLine); 

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

            delta_y = 0.01 * (max_y - min_y);
            min_y -= delta_y;
            max_y += delta_y;

            painter.setPen(pen_black);
            painter.drawLine(L2G(a, 0), L2G(b, 0));
            painter.drawLine(L2G(0, min_y), L2G(0, max_y));

            painter.setPen(pen_blue);
            x1 = a;
            y1 = f (x1);
            for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
                y2 = f (x2);
                painter.drawLine(L2G(x1, y1), L2G(x2, y2));
                x1 = x2;
                y1 = y2;
            }
            x2 = b;
            y2 = f(x2);
            painter.drawLine(L2G(x1, y1), L2G(x2, y2));

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
                    title = QString("k = 0 f(x) = 1, n = %1, s = %2").arg(n).arg(s); 
                    break;
                case 1: 
                    title = QString("k = 1 f(x) = x, n = %1, s = %2").arg(n).arg(s); 
                    break;
                case 2: 
                    title = QString("k = 2 f(x) = x^2, n = %1, s = %2").arg(n).arg(s); 
                    break;
                case 3: 
                    title = QString("k = 3 f(x) = x^3, n = %1, s = %2").arg(n).arg(s); 
                    break;
                case 4: 
                    title = QString("k = 4 f(x) = x^4, n = %1, s = %2").arg(n).arg(s); 
                    break;
                case 5: 
                    title = QString("k = 5 f(x) = exp(x), n = %1, s = %2").arg(n).arg(s); 
                    break;
                case 6: 
                    title = QString("k = 6 f(x) = 1/(25*x^2 + 1), n = %1, s = %2").arg(n).arg(s); 
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

        ~MainWindow() {};
};

