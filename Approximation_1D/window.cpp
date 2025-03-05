#include <QPainter>
#include <stdio.h>

#include "window.hpp"

#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_N 10
#define L2G(X,Y) (l2g ((X), (Y), min_y, max_y))

static double f_0 (double x) {
    return x;
}

static double f_1 (double x) {
    return x * x * x;
}

Window::Window (QWidget *parent): QWidget (parent) {
    a = DEFAULT_A;
    b = DEFAULT_B;
    n = DEFAULT_N;

    func_id = 0;

    change_func ();
}

QSize Window::minimumSizeHint () const {
    return QSize (100, 100);
}

QSize Window::sizeHint () const {
    return QSize (1000, 1000);
}

int Window::parse_command_line (int argc, char *argv[]) {
    if (argc == 1)
        return 0;

    if (argc == 2)
        return -1;

    if (   sscanf (argv[1], "%lf", &a) != 1
            || sscanf (argv[2], "%lf", &b) != 1
            || b - a < 1.e-6
            || (argc > 3 && sscanf (argv[3], "%d", &n) != 1)
            || n <= 0)
        return -2;

    return 0;
}

/// change current function for drawing
void Window::change_func () {
    func_id = (func_id + 1) % 2;

    switch (func_id) {
        case 0:
            f_name = "f (x) = x";
            f = f_0;
            break;
        case 1:
            f_name = "f (x) = x * x * x";
            f = f_1;
            break;
    }
    update ();
}

QPointF Window::l2g (double x_loc, double y_loc, double y_min, double y_max) {
    double x_gl = (x_loc - a) / (b - a) * width ();
    double y_gl = (y_max - y_loc) / (y_max - y_min) * height ();
    return QPointF (x_gl, y_gl);
}

/// render graph
void Window::paintEvent (QPaintEvent * /* event */) {  
    QPainter painter (this);
    double x1, x2, y1, y2;
    double max_y, min_y;
    double delta_y, delta_x = (b - a) / n;
    QPen pen_black(Qt::black, 0, Qt::SolidLine); 
    QPen pen_red(Qt::red, 0, Qt::SolidLine); 

    painter.setPen (pen_black);

    // calculate min and max for current function
    max_y = min_y = 0;
    for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
        y1 = f (x1);
        if (y1 < min_y)
            min_y = y1;
        if (y1 > max_y)
            max_y = y1;
    }

    delta_y = 0.01 * (max_y - min_y);
    min_y -= delta_y;
    max_y += delta_y;

    // draw approximated line for graph
    x1 = a;
    y1 = f (x1);
    for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
        y2 = f (x2);
        // local coords are converted to draw coords
        painter.drawLine (L2G(x1, y1), L2G(x2, y2));

        x1 = x2, y1 = y2;
    }
    x2 = b;
    y2 = f (x2);
    painter.drawLine (L2G(x1, y1), L2G(x2, y2));

    // draw axis
    painter.setPen (pen_red);
    painter.drawLine (L2G(a, 0), L2G(b, 0));
    painter.drawLine (L2G(0, min_y), L2G(0, max_y));

    // render function name
    painter.setPen ("blue");
    painter.drawText (0, 20, f_name);

}
