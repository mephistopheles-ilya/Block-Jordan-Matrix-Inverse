#include <stdio.h>
#include <algorithm>
#include <cstring>

#include "read_print_fill.hpp"

int get_max_rows(int n, int m, int p) {
    int b = (n + m - 1)/m;
    return (b + p - 1)/p;
}

int get_rows(int n, int m, int p, int proc_num) {
    int b = (n + m - 1) / m;
    return (proc_num >= b%p ? b/p : (b + p - 1)/p);
}

int ltg(int m, int p, int k, int i_loc) {
    int i_loc_m = i_loc/m;
    int i_glob_m = i_loc_m * p + k;
    return i_glob_m * m + i_loc % m;
}

double f(int s, int n, int i_glob, int j_glob) {
    switch (s) {
        case 1: return n - std::max(i_glob, j_glob) + 1;
        case 2: return std::max(i_glob, j_glob);
        case 3: return std::abs(i_glob - j_glob);
        case 4: return 1./(i_glob + j_glob - 1);
        default: return 0;
    }
}

int fill_matrix(double* matrix, int n, int m, int s, char* file_name, int proc_num, int p, double* buf
        , MPI_Comm comm) {
    FILE* file = nullptr;
    int err = 0;
    if (s == 0) {
        if (proc_num == 0) {
            file = fopen(file_name, "r");
            if (!file) {
                err = 1;
            }
        }
        MPI_Bcast(&err, 1, MPI_INT, 0, comm);
        if (err != 0) {
            return err;
        }
        memset(buf, 0, n * m * sizeof(double));
        int counter = 0;
        int k = n/m;
        int l = n - m * k;
        int line_of_blocks = 0, line_in_block = 0, block_in_line = 0, num_line_in_block = 0, num_elem_in_line = 0;
        for(line_of_blocks = 0; line_of_blocks <= k; ++line_of_blocks) {
            num_line_in_block = (line_of_blocks < k) ? m : l;
            for(line_in_block = 0; line_in_block < num_line_in_block; ++line_in_block) {
                for(block_in_line = 0; block_in_line <= k; ++block_in_line) {
                    num_elem_in_line = (block_in_line < k) ? m : l;
                    int shift = line_in_block * num_line_in_block + block_in_line * m * m;
                    if (proc_num == 0) {
                        for(int q = 0; q < num_elem_in_line; ++q) {
                            if (fscanf(file, "%lf", buf + shift + q) == 1) {
                                ++counter;
                            } else {
                                err = 2;
                                goto lable;
                            }
                        }
                    }
                }
            }
lable:
            MPI_Bcast(&err, 1, MPI_INT, 0, comm);
            MPI_Bcast(&counter, 1, MPI_INT, 0, comm);
            if (err != 0 && counter != n * n) {
                if (proc_num == 0) {
                    fclose(file);
                }
                return err;
            }
            int owner = line_of_blocks % p;
            int line_of_blocks_loc = line_of_blocks / p;
            if (proc_num == 0) {
                if (owner == 0) {
                    memcpy(matrix + line_of_blocks_loc * (m * n * k + l * m), buf, k * num_line_in_block * m + l * 
                            num_line_in_block);
                } else {
                    MPI_Send(buf, k * num_line_in_block * m + l * num_line_in_block, MPI_DOUBLE, owner, 0, comm);
                }
            } else {
                if (proc_num == owner) {
                    MPI_Status st;
                    MPI_Recv(matrix + line_of_blocks_loc * (m * n * k + l * m), k * num_line_in_block * m + 
                            l * num_line_in_block, MPI_DOUBLE, 0, 0, comm, &st);
                }
            }        
        }
        if (proc_num == 0) {
            fclose(file);
            return 0;
        }
    } else {
        int rows = get_rows(n, m, p, proc_num);

        int k = n/m;
        int l = n - m * k;
        int b = (n + m - 1) / m;

        int line_of_blocks = 0, line_in_block = 0, block_in_line = 0;
        int i_loc = 0, j_loc = 0;
        int i_glob = 0, j_glob = 0;

        for(line_of_blocks = 0; line_of_blocks < rows; ++line_of_blocks) {
            int num_line_in_block = m;
            if(l != 0 && (b - 1)%p == proc_num && line_of_blocks == (rows - 1)) {
                num_line_in_block = l;
            }
            for(line_in_block = 0; line_in_block < num_line_in_block; ++line_in_block) {
                i_loc = line_of_blocks * m + line_in_block;
                for(block_in_line = 0; block_in_line <= k; ++block_in_line) {
                    int num_elem_in_line = (block_in_line < k) ? m : l;
                    j_loc = block_in_line * m + 1;
                    int shift = line_of_blocks * n * m + line_in_block * num_line_in_block + 
                        block_in_line  * m * m;
                    for(int q = 0; q < num_elem_in_line; ++q) {
                        j_glob = j_loc;
                        i_glob = ltg(m, p, proc_num, i_loc);
                        *(matrix + shift + q) = f(s, n, i_glob + 1, j_glob + 1);
                        ++j_loc;
                    }
                }
            }
        }
    }
    return 0;
}

int fill_id_matrix(double* inverse, int n, int m, int proc_num, int p) {
    int rows = get_rows(n, m, p, proc_num);
    int k = n/m;
    int l = n - m * k;
    int b = (n + m - 1)/m;
    int line_of_blocks = 0, line_in_block = 0, block_in_line = 0;
    int i_loc = 0, j_loc = 0;
    int i_glob = 0, j_glob = 0;

    for(line_of_blocks = 0; line_of_blocks < rows; ++line_of_blocks) {
        int num_line_in_block = m;
        if(l != 0 && (b - 1)%p == proc_num && line_of_blocks == (rows - 1)) {
            num_line_in_block = l;
        }
        for(line_in_block = 0; line_in_block < num_line_in_block; ++line_in_block) {
            i_loc = line_of_blocks * m + line_in_block;
            for(block_in_line = 0; block_in_line <= k; ++block_in_line) {
                int num_elem_in_line = (block_in_line < k) ? m : l;
                j_loc = block_in_line * m + 1;
                int shift = line_of_blocks * n * m + line_in_block * num_line_in_block + 
                    block_in_line  * m * m;
                for(int q = 0; q < num_elem_in_line; ++q) {
                    j_glob = j_loc;
                    i_glob = ltg(m, p, proc_num, i_loc);
                    *(inverse + shift + q) =  (i_glob + 1 == j_glob + 1) ? 1 : 0;
                    ++j_loc;
                }
            }
        }
    }
    return 0;
}

int print_matrix(double* matrix, int n, int m, int proc_num, int p, int r, double* buf, MPI_Comm comm) {
    if (r > n) r = n;
    int k = n / m;
    int l = n - m * k;
    int line_of_blocks = 0, block_in_line = 0, line_in_block = 0;
    int num_line_in_block = 0, num_elem_in_line = 0;
    int r_num_line_in_block = 0, r_num_elem_in_line = 0;

    int restrict_k = r / m;
    int restrict_l = r - m * restrict_k;

    for(line_of_blocks =  0; line_of_blocks <= restrict_k; ++line_of_blocks) {
        num_line_in_block = (line_of_blocks < k) ? m : l;
        r_num_line_in_block = (line_of_blocks < restrict_k) ? m : restrict_l;
        int owner = line_of_blocks % p;
        int line_of_blocks_loc = line_of_blocks / p;
        if (proc_num == 0) {
            if (owner == 0) {
                memcpy(buf, matrix + line_of_blocks_loc*(m * m * k + l * m), k * num_line_in_block * m 
                        + num_line_in_block * l);
            } else {
                MPI_Status st;
                MPI_Recv(buf, k * num_line_in_block * m + num_line_in_block * l, MPI_DOUBLE, owner, 0, comm, &st);
            }
        } else {
            if (owner == proc_num) {
                MPI_Send(matrix + line_of_blocks_loc*(m * m * k + l * m)
                        ,  k * num_line_in_block * m + num_line_in_block * l, MPI_DOUBLE, 0, 0, comm);
            }
        }

        if (proc_num == 0) {
            for(line_in_block = 0; line_in_block < r_num_line_in_block; ++line_in_block) {
                for(block_in_line = 0; block_in_line <= restrict_k; ++block_in_line) {
                    num_elem_in_line = (block_in_line < k) ? m : l;
                    r_num_elem_in_line = (block_in_line < restrict_k) ? m : restrict_l;
                    int shift = block_in_line * m * num_line_in_block + line_in_block * num_elem_in_line;
                    for(int q = 0; q < r_num_elem_in_line; ++q) {
                        printf(" %10.3e", *(buf + shift + q));
                    }
                }
            }
        }
    }
    return 0;
}



