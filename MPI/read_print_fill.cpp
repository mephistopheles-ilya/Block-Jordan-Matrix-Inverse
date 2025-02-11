#include <stdio.h>

int get_max_rows(int n, int m, int p) {
    int b = (n + m - 1)/m;
    return (b + p - 1)/p;
}

int get_rows(int n, int m, int p, int proc_num) {
    int b = (n + m - 1) / m;
    //return (b%p <= proc_num ? b/p : b/p + 1);
    return (proc_num > b%p ? b/p : (b + p - 1)/p);
}

int fill_matrix(double* matrix, int n, int m, int s, char* file_name, int proc_num, int p) {
    FILE* file = nullptr;
    if (s == 0) {
        file = fopen(file_name, "r");
        if (!file) {
            return 1;
        }
    } else {

        int max_rows = get_max_rows(n, m, p);
        int rows = get_rows(n, m, p, proc_num);

        int k = n/m;
        int l = n - m * k;

        int line_of_blocks = 0, line_in_block = 0, block_in_line = 0;
        int i_loc = 0, j_loc = 0;

        for(line_of_blocks = 0; line_of_blocks <= k; ++line_of_blocks) {
            int num_line_in_block = (line_of_blocks < k) ? m : l;
            for(line_in_block = 0; line_in_block < num_line_in_block; ++line_in_block) {
                i_loc = line_of_blocks * m + line_in_block + 1;
                for(block_in_line = 0; block_in_line <= max_rows; ++block_in_line) {
                    int num_elem_in_line = m;
                    if (block_in_line == max_rows) {
                        int b = (n + m - 1) / m;
                        if (rows < max_rows) {
                            num_elem_in_line = 0;
                        } else {
                            if (b%p == proc_num && n%m != 0) {
                                num_elem_in_line = l;
                            }
                        }
                    }
                    j_loc = block_in_line * m + 1;
                    int shift = 0; 
                    for(int q = 0; q < num_elem_in_line; ++q) {
                        *(matrix + shift) = 0;
                        ++j_loc;
                    }
                }
            }
        }


    }
    return 0;
}
