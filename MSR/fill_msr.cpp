#include "fill_msr.hpp"

#include <new>

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

#define F(I, J) do{ij2l(nx, ny, I, J, k); if(I_ij) {I_ij[m] = k;} ++m;} while(0)

int get_off_diag(int nx, int ny, int i, int j, int* I_ij) {
    int m = 0, k = 0;
    if (i < nx) { F(i + 1, j); }
    if (j > 0) { F(i, j - 1); }
    if (i > 0 && j > 0) { F(i -1 , j - 1); }
    if (i > 0) { F(i - 1, j); }
    if (j < ny) { F(i, j + 1); }
    if (i < nx && i < ny) { F(i + 1, j + 1); }
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
        






