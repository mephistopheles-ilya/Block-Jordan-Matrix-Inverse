#include "fill_msr.hpp"
#include "utils.hpp"

#include <cmath>
 
/*__________
< govnokod >
 ----------
        \   ^__^
         \  (oo)\_______
            (__)\       )\/\
                ||----w |
                ||     ||
*/

double Pf(double* res, double x, double y, double a, double c, double hx, double hy, int nx, int ny) {
    int i = (x - a) / hx;
    int j = (y - c) / hx;
    double k[8] = {0};
    double x0 = 0., y0 = 0., z0 = 0.;
    double x1 = 0., y1 = 0., z1 = 0.;
    double x2 = 0., y2 = 0., z2 = 0.;
    int l0 = 0, l1 = 0, l2 = 0;

    ij2l(nx, ny, i + 0, j + 0, l0);
    ij2l(nx, ny, i + 1, j + 1, l2);
    x0 = a + (i + 0) * hx;
    y0 = c + (j + 0) * hy;
    z0 = res[l0];
    x2 = a + (i + 1) * hx;
    y2 = c + (j + 1) * hy;
    z2 = res[l2];

    if (hy/hx * (y - y0) >= (x - x0)) {
        ij2l(nx, ny, i + 1, j + 0, l1);
        x1 = a + (i + 1) * hx;
        y1 = c + (j + 0) * hy;
        z1 = res[l1];
    } else {
        ij2l(nx, ny, i + 0, j + 1, l1);
        x1 = a + (i + 0) * hx;
        y1 = c + (j + 1) * hy;
        z1 = res[l1];
    }


    k[0] = x - x0;
    k[1] = y - y0;
    k[2] = x1 - x0;
    k[3] = y1 - y0;
    k[4] = z1 - z0;
    k[5] = x2 - x0;
    k[6] = y2 - y0;
    k[7] = z2 - z0;

    double num = k[3] * k[5] * z0 - k[2] * k[6] * z0 + k[1] * k[4] * k[5] - k[0] * k[4] * k[6] 
        - k[1] * k[2] * k[7] + k[0] * k[3] * k[7];
    double den = k[3] * k[5] - k[2] * k[6];

    return num / den;
}

double calc_r1(double* res, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double)) {
    int j = 0, j1 = 0, j2 = 0;
    double max = -1, cur_val = 0;
    thread_rows(ny - 1, p, k, j1, j2);
    for(j = j1; j < j2; ++j) {
        for(int i = 0; i < nx; ++i) {
            cur_val = std::fabs(f(a + (i + 2./3) * hx, c + (j + 1./3) * hy) 
                    - Pf(res, a + (i + 2./3) * hx, c + (j + 1./3) * hy, a, c, hx, hy, nx, ny));
            if (cur_val > max) {
                max = cur_val;
            }
            cur_val = std::fabs(f(a + (i + 1./3) * hx, c + (j + 2./3) * hy) 
                    - Pf(res, a + (i + 1./3) * hx, c + (j + 2./3) * hy, a, c, hx, hy, nx, ny));
            if (cur_val > max) {
                max = cur_val;
            }

        }
    }
    reduce_max(p, &max, 1);
    return max;
}

double calc_r2(double* res, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double)) {
    int j = 0, j1 = 0, j2 = 0;
    double sum = 0;
    thread_rows(ny - 1, p, k, j1, j2);
    for(j = j1; j < j2; ++j) {
        for(int i = 0; i < nx; ++i) {
            sum += std::fabs(f(a + (i + 2./3) * hx, c + (j + 1./3) * hy) 
                    - Pf(res, a + (i + 2./3) * hx, c + (j + 1./3) * hy, a, c, hx, hy, nx, ny));
            sum += std::fabs(f(a + (i + 1./3) * hx, c + (j + 2./3) * hy) 
                    - Pf(res, a + (i + 1./3) * hx, c + (j + 2./3) * hy, a, c, hx, hy, nx, ny));

        }
    }
    sum *= hx * hx * 0.5;
    reduce_sum_double_det(p, k, sum);
    return sum;
}

double calc_r3(double* res, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double)) {
    int j = 0, j1 = 0, j2 = 0, l = 0;
    double max = -1, cur_val = 0;
    thread_rows(ny, p, k, j1, j2);
    for(j = j1; j < j2; ++j) {
        for(int i = 0; i <= nx; ++i) {
            ij2l(nx, ny, i, j, l);
            cur_val = std::fabs(f(a + i * hx, c + j * hy) - res[l]);
            if (cur_val > max) {
                max = cur_val;
            }
        }
    }
    reduce_max(p, &max, 1);
    return max;
}

double calc_r4(double* res, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double)) {
    int j = 0, j1 = 0, j2 = 0, l = 0;
    double sum = 0;
    thread_rows(ny, p, k, j1, j2);
    for(j = j1; j < j2; ++j) {
        for(int i = 0; i <= nx; ++i) {
            ij2l(nx, ny, i, j, l);
            sum += std::fabs(f(a + i * hx, c + j * hy) - res[l]);
        }
    }
    sum *= hx * hy;
    reduce_sum_double_det(p, k, sum);
    return sum;
}
           








