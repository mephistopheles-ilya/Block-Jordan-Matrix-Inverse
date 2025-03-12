#include <algorithm>
#include <cmath>

void Newton_approximation(double* points_x, double* f_x, double* newton_coeff, int n) {
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

double Newton_calck(double* newton_coeff, double* points_x, int n, double x) {
    double res = newton_coeff[n - 1];
    for(int i = n - 2; i >= 0; --i) {
        res *= (x - points_x[i]);
        res += newton_coeff[i];
    }
    return res;
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

void cube_approximation(double* points_x, double* f_x, double* d, double* c, int n) {
    for(int i = 1; i < n - 1; ++i) {
        double df1 = (f_x[i] - f_x[i - 1])/(points_x[i] - points_x[i - 1]);
        double df2 = (f_x[i + 1] - f_x[i])/(points_x[i + 1] - points_x[i]);
        if (df1 * df2 > 0) {
            double sign = (df1 > 0) ? 1 : -1;
            d[i] = sign * std::min(std::fabs(df1), std::fabs(df2));
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

double cube_calc(double* points_x, double* c, double a, double b, int n, double x) {
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

