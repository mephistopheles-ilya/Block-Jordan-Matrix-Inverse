#include "fill_msr.hpp"
#include "utils.hpp"

#include <new>
#include <cstdlib>

void ij2l(int nx, int /*ny*/, int i, int j, int& l) {
    l = i + j * (nx + 1);
}

void l2ij(int nx, int /*ny*/, int& i, int& j, int l) {
    j = l / (nx + 1);
    i = l - j * (nx + 1);
}

int get_len_msr(int nx, int ny) {
    return (nx + 1) * (ny + 1) 
        + 6 * (nx - 1) * (ny - 1) + (2 * (nx - 1) + 2 * (ny - 1)) * 4 + 2 * 3 + 2 * 2;
}

#define Fa(I, J) do{ij2l(nx, ny, I, J, k); if(I_ij) {I_ij[m] = k;} ++m;} while(0)

int get_off_diag(int nx, int ny, int i, int j, int* I_ij) {
    int m = 0, k = 0;
    if (i < nx) { Fa(i + 1, j); }
    if (j > 0) { Fa(i, j - 1); }
    if (i > 0 && j > 0) { Fa(i -1 , j - 1); }
    if (i > 0) { Fa(i - 1, j); }
    if (j < ny) { Fa(i, j + 1); }
    if (i < nx && j < ny) { Fa(i + 1, j + 1); }
    return m;
}

int get_len_msr_off_diag(int nx, int ny) {
    int m = 0, i = 0, j = 0;
    for(i = 0; i <= nx; ++i) {
        for (j = 0; j <= ny; ++j) {
            m += get_off_diag(nx, ny, i, j);
        }
    }
    return m;
}

int alocate_msr_matrix(int nx, int ny, double** p_A, int** p_I) {
    int diag_len = (nx + 1) * (ny + 1);
    int off_diag = get_len_msr_off_diag(nx, ny);
    int len = diag_len + off_diag + 1;
    double* A = nullptr;
    int* I = nullptr;
    A = new (std::nothrow) double[len];
    I = new (std::nothrow) int[len];
    if (A == nullptr || I == nullptr) {
        delete[] A;
        delete[] I;
        return 1;
    }
    *p_A = A;
    *p_I = I;
    return 0;
}


void fill_I(int nx, int ny, int* I) {
    int N = (nx + 1) * (ny + 1);
    int l = 0, i = 0, j = 0, r = N + 1, m = 0;
    for(l = 0; l < N; ++l) {
        l2ij(nx, ny, i, j, l);
        I[l] = r;
        m = get_off_diag(nx, ny, i, j, I + r);
        r += m;
    }
    I[l] = r;
}

void fill_A_ij(int nx, int ny, double hx, double hy, int i,int j, double* A_diag, double* A_offdiag) {
    double s = hx * hy;
    int l = 0;
    if (i > 0 && i < nx && j > 0 && j < ny) {
        A_diag[0] = 6 * s / 12.;
        for(l = 0; l < 6; ++l) {
            A_offdiag[l] = s / 12.;
        }
        return;
    } else if (j == 0 && i > 0 && i < nx) {
        A_diag[0] = 3 * s / 12.;
        A_offdiag[0] = 1 * s / 24.;
        A_offdiag[1] = 1 * s / 24.;
        A_offdiag[2] = 2 * s / 24.;
        A_offdiag[3] = 2 * s / 24.;
        return;
    } else if (j == ny && i > 0 && i < nx) {
        A_diag[0] = 3 * s / 12.;
        A_offdiag[0] = 1 * s / 24.;
        A_offdiag[1] = 2 * s / 24.;
        A_offdiag[2] = 2 * s / 24.;
        A_offdiag[3] = 1 * s / 24.;
        return;
    } else if (i == 0 && j > 0 && j < ny) {
        A_diag[0] = 3 * s / 12.;
        A_offdiag[0] = 2 * s / 24.;
        A_offdiag[1] = 1 * s / 24.;
        A_offdiag[2] = 1 * s / 24.;
        A_offdiag[3] = 2 * s / 24.;
        return;
    } else if (i == nx && j > 0 && j < ny) {
        A_diag[0] = 3 * s / 12.;
        A_offdiag[0] = 1 * s / 24.;
        A_offdiag[1] = 2 * s / 24.;
        A_offdiag[2] = 2 * s / 24.;
        A_offdiag[3] = 1 * s / 24.;
        return;
    } else if (i == 0 && j == 0) {
        A_diag[0] = 2 * s / 12.;
        A_offdiag[0] = 1 * s / 24.;
        A_offdiag[1] = 1 * s / 24.;
        A_offdiag[2] = 2 * s / 24.;
        return;
    } else if (i == nx && j == ny) {
        A_diag[0] = 2 * s / 24.;
        A_offdiag[0] = 1 * s / 24.;
        A_offdiag[1] = 2 * s / 24.;
        A_offdiag[2] = 1 * s / 24.;
        return;
    } else if (i == 0 && j == ny) {
        A_diag[0] = 1 * s / 12.;
        A_offdiag[0] = 1 * s / 24.;
        A_offdiag[1] = 1 * s / 24.;
        return;
    } else if (i == nx && j == 0) {
        A_diag[0] = 1 * s / 12.;
        A_offdiag[0] = 1 * s / 24.;
        A_offdiag[1] = 1 * s / 24.;
        return;
    }

    make_fpe();
}

void fill_A(int nx, int ny, double hx, double hy, int* I, double* A, int p, int k) {
    int l1 = 0, l2 = 0, l = 0;
    int N = (nx + 1) * (ny + 1);
    int i = 0, j = 0;
    l1 = N * k;
    l1 /= p;
    l2 = N * (k + 1);
    l2 /= p;
    for(l = l1; l < l2; ++l) {
        l2ij(nx, ny, i, j, l);
        double* A_diag = A + l;
        double* A_offdiag = A + I[l];
        fill_A_ij(nx, ny, hx, hy, i, j, A_diag, A_offdiag);
    }
    barrier(p);
}

#define Fb(I, J) (f(x0 + (I) * hx, y0 + (J) * hy))

double F_IJ(int nx, int ny, double hx, double hy, double x0, double y0, int i, int j, double (*f)(double, double)) {
    double w = hx * hy / 192;
    if (i > 0 && i < nx && j > 0 && j < ny) {
        return w * (
                36 * Fb(i, j)
                + 20 * (
                    Fb(i + 0.5, j) + Fb(i, j - 0.5) + Fb(i - 0.5, j - 0.5)
                    + Fb(i - 0.5, j) + Fb(i, j + 0.5) + Fb(i + 0.5, j + 0.5)
                    )
                + 4 * (
                    Fb(i + 0.5, j - 0.5) + Fb(i - 0.5, j - 1) + Fb(i - 1, j - 0.5)
                    + Fb(i - 0.5, j + 0.5) + Fb(i + 0.5, j + 1) + Fb(i + 1, j + 0.5)
                    )
                + 2 * (
                    Fb(i + 1, j) + Fb(i, j - 1) + Fb(i - 1, j - 1) 
                    + Fb(i - 1, j) + Fb(i, j + 1) + Fb(i + 1, j + 1)
                    )
                );
    }
    if (i > 0 && i < nx && j == 0) {
        return w * (
                18 * Fb(i, j)
                + 10 * (Fb(i + 0.5, j) + Fb(i - 0.5, j))
                + 20 * (Fb(i, j + 0.5) + Fb(i + 0.5, j + 0.5))
                + 4 * (Fb(i - 0.5, j + 0.5) + Fb(i + 0.5, j + 1) + Fb(i + 1, j + 0.5))
                + 1 * (Fb(i - 1, j) + Fb(i + 1, j))
                + 2 * (Fb(i, j + 1) + Fb(i + 1, j + 1))
                );
    }
    if (i > 0 && i < nx && j == ny) {
        return w * (
                18 * Fb(i, j)
                + 10 * (Fb(i - 0.5, j) + Fb(i + 0.5, j))
                + 20 * (Fb(i, j - 0.5) + Fb(i - 0.5, j - 0.5))
                + 4 * (Fb(i + 0.5, j - 0.5) + Fb(i - 0.5, j - 1) + Fb(i - 1, j - 0.5))
                + 1 * (Fb(i - 1, j) + Fb(i + 1, j))
                + 2 * (Fb(i, j - 1) + Fb(i - 1, j - 1))
                );
    }
    if (i == 0 && j > 0 && j < ny) {
        return w * (
                18 * Fb(i, j)
                + 10 * (Fb(i, j - 0.5) + Fb(i, j + 0.5))
                + 20 * (Fb(i + 0.5, j) + Fb(i + 0.5, j + 0.5))
                + 4 * (Fb(i + 0.5, j - 0.5) + Fb(i + 0.5, j + 1) + Fb(i + 1, j + 0.5))
                + 1 * (Fb(i, j - 1) + Fb(i, j + 1))
                + 2 * (Fb(i + 1, j) + Fb(i + 1, j + 1))
                );
    }
    if (i == nx && j > 0 && j < ny) {
        return w * (
                18 * Fb(i, j)
                + 10 * (Fb(i, j - 0.5) + Fb(i, j + 0.5))
                + 20 * (Fb(i - 0.5, j) + Fb(i - 0.5, j + 0.5))
                + 4 * (Fb(i - 0.5, j) + Fb(i - 1, j - 0.5) + Fb(i - 0.5, j + 0.5))
                + 1 * (Fb(i, j - 1) + Fb(i, j + 1))
                + 2 * (Fb(i - 1, j) + Fb(i - 1, j - 1))
                );
    }
    if (i == 0 && j == 0) {
        return w * (
                12 * Fb(i, j)
                + 10 * (Fb(i + 0.5, j) + Fb(i, j + 0.5))
                + 20 * (Fb(i + 0.5, j + 0.5))
                + 4 * (Fb(i + 1, j + 0.5) + Fb(i + 0.5, j + 1))
                + 1 * (Fb(i + 1, j) + Fb(i, j + 1))
                + 2 * (Fb(i + 1, j + 1))
                );
    }
    if (i == nx && j == ny) {
        return w * (
                12 * Fb(i, j)
                + 10 * (Fb(i - 0.5, j) + Fb(i, j - 0.5))
                + 20 * (Fb(i - 0.5, j - 0.5))
                + 4 * (Fb(i - 0.5, j - 1) + Fb(i - 1, j - 0.5))
                + 1 * (Fb(i, j - 1) + Fb(i - 1, j))
                + 2 * (Fb(i - 1, j - 1))
                );
    }
    if (i == 0 && j == ny) {
        return w * (
                6 * Fb(i, j)
                + 10 * (Fb(i + 0.5, j) + Fb(i, j - 0.5))
                + 4 * (Fb(i + 0.5, j - 0.5))
                + 1 * (Fb(i + 1, j) + Fb(i, j - 1))
                );
    }
    if (i == nx && j == 0) {
        return w * (
                6 * Fb(i, j)
                + 10 * (Fb(i - 0.5, j) + Fb(i, j + 0.5))
                + 4 * (Fb(i - 0.5, j + 0.5))
                + 1 * (Fb(i - 1, j) + Fb(i, j + 1))
                );
    }

    make_fpe();
    return 1e308;
}

void fill_B(double a, double c, int nx, int ny, double hx, double hy, double* b, int p, int k, double (*f)(double, double)) {
    int l1 = 0, l2 = 0, l = 0;
    int N = (nx + 1) * (ny + 1);
    int i = 0, j = 0;
    l1 = N * k;
    l1 /= p;
    l2 = N * (k + 1);
    l2 /= p;
    for(l = l1; l < l2; ++l) {
        l2ij(nx, ny, i, j, l);
        b[l] = F_IJ(nx, ny, hx, hy, a, c, i, j, f); 
    }
    barrier(p);
}


    


        






