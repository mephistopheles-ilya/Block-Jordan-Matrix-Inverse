#pragma once


void ij2l(int nx, int /*ny*/, int i, int j, int& l);
void l2ij(int nx, int /*ny*/, int& i, int& j, int l);
int get_len_msr(int nx, int ny);
int get_off_diag(int nx, int ny, int i, int j, int* I_ij = nullptr);
int get_len_msr_off_diag(int nx, int ny);
