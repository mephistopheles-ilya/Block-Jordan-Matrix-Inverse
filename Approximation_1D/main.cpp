#include <QApplication>
#include <QMessageBox>

#include <stdio.h>

//#include <fenv.h>

#include "window.hpp"

int main(int argc, char *argv[]) {

    //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

    QApplication app(argc, argv);

    double a = 0, b = 0;
    int k = 0, n = 0;
    if (argc != 5 || sscanf(argv[1], "%lf", &a) != 1 || sscanf(argv[2], "%lf", &b) != 1 
                || sscanf(argv[3], "%d", &n) != 1 || sscanf(argv[4], "%d", &k) != 1) {
        QMessageBox::critical(nullptr, "Error", "Wrong amoun of arguments");
        return 1;
    }
    if (a > b || k < 0 || k > 6 || n < 5) {
        QMessageBox::critical(nullptr, "Error", "Wrong arguments");
        return 2;
    }

    MainWindow window(a, b, n, k);
    window.show();

    int ret_code = 0;
    ret_code = app.exec();
    printf("QApplication return code %d\n", ret_code);
    return 0;
}
