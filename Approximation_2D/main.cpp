#include <QApplication>
#include <QMessageBox>

#include <stdio.h>

#include "window.hpp"

#include <fenv.h>

int main(int argc, char* argv[]) {

    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

    QApplication app(argc, argv);

    double a = -1., b = -1., c = -1., d = -1.;
    int nx = -1, ny = -1;
    int mx = -1, my = -1;
    int k = -1;
    double eps = -1;
    int mi = -1;
    int p = -1;


    if (argc != 13 
        || sscanf(argv[1], "%lf", &a) != 1 || sscanf(argv[2], "%lf", &b) != 1 
        || sscanf(argv[3], "%lf", &c) != 1 || sscanf(argv[4], "%lf", &d) != 1
        || sscanf(argv[5], "%d", &nx) != 1 || sscanf(argv[6], "%d", &ny) != 1
        || sscanf(argv[7], "%d", &mx) != 1 || sscanf(argv[8], "%d", &my) != 1
        || sscanf(argv[9], "%d", &k) != 1 || sscanf(argv[10], "%lf", &eps) != 1
        || sscanf(argv[11], "%d", &mi) != 1 || sscanf(argv[12], "%d", &p) != 1) {

        std::string message = std::string("Usage ") + std::string(argv[0]) + std::string(" a b c d nx ny mx my k eps mi p\n");
        QMessageBox::critical(nullptr, "Error", message.data());
        return 1;
    }

    if (b <= a || d <= c || nx < 5 || ny < 5 || mx < 5 || my < 5 || k < 0 || k > 7 || eps < 0 || mi < 0 || p <= 0) {
        QMessageBox::critical(nullptr, "Error", "Wrong parametrs of commad line\n");
        return 2;
    }

    MainWindow window(a, b, c, d, nx, ny, mx, my, k, eps, mi, p, argv[0]);
    window.show();
    app.exec();

 
    return 0 ;
}
