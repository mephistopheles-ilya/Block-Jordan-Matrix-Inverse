#include "functions.h"
#include <stdio.h>

int read_matrix_from_file(double *matrix, int n, int m, char *filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        return 1;
    }
    /* 
     * n = m * k + l
     */
    int counter = 0;
    int num_in_line1; //amount of lines in block
    int num_in_line2; // size of lines in blcok
    double* cur_pointer;
    int tmp;
    int k = n / m;
    int l = n - m * k;
    int elements_in_blcok_line = k * m * m + l * m;
    int line_of_blocks, block_in_line, line_in_block, i;

    for(line_of_blocks =  0; line_of_blocks <= k; ++line_of_blocks) {
        num_in_line1 = (line_of_blocks < k) ? m : l;
        tmp = line_of_blocks * elements_in_blcok_line;
        for(line_in_block = 0; line_in_block < num_in_line1; ++line_in_block) {
            for(block_in_line = 0; block_in_line <= k; ++block_in_line) {
                num_in_line2 = (block_in_line < k) ? m : l;
                cur_pointer = matrix + tmp + block_in_line * m * num_in_line1
                    + line_in_block * num_in_line2;
                for(i = 0; i < num_in_line2; ++i) {
                    if(fscanf(file, "%lf", (cur_pointer + i)) == 1) {
                        ++counter;
                    } else {
                        goto lable;
                    }
                }
            }
        }
    }

    lable:
    fclose(file);
    if (counter != n * n) {
        return 1;
    } else {
        return 0;
    }
}

void create_id_matrix(double *matrix, int n, int m) {
    int num_in_line1; //amount of lines in block
    int num_in_line2; // size of lines in blcok
    double* cur_pointer;
    int tmp;
    int k = n / m;
    int l = n - m * k;
    int elements_in_blcok_line = k * m * m + l * m;
    int line_of_blocks, block_in_line, line_in_block, q;
    int i, j;

    for(line_of_blocks =  0; line_of_blocks <= k; ++line_of_blocks) {
        num_in_line1 = (line_of_blocks < k) ? m : l;
        tmp = line_of_blocks * elements_in_blcok_line;
        for(line_in_block = 0; line_in_block < num_in_line1; ++line_in_block) {
            i = line_of_blocks * m + line_in_block + 1;
            for(block_in_line = 0; block_in_line <= k; ++block_in_line) {
                num_in_line2 = (block_in_line < k) ? m : l;
                cur_pointer = matrix + tmp + block_in_line * m * num_in_line1
                    + line_in_block * num_in_line2;
                j = block_in_line * m + 1;
                for(q = 0; q < num_in_line2; ++q) {
                    *(cur_pointer + q) = (i == j) ? 1 : 0;
                    ++j;
                }
            }
        }
    }
}


double f2(int i, int j) {
    return (i > j) ? i : j;
}
double f3(int i, int j ) {
    return ((i - j) < 0) ? (j - i) : (i - j);
}
double f4(int i, int j) {
    return 1./(i + j - 1);
}


void fill_matrix_with_formula(double *matrix, int n, int m, int s) {
    int num_in_line1; //amount of lines in block
    int num_in_line2; // size of lines in blcok
    double* cur_pointer;
    int k = n / m;
    int l = n - m * k;
    int elements_in_blcok_line = k * m * m + l * m;
    int line_of_blocks, block_in_line, line_in_block, q;
    int i, j, tmp;
    double (*p)(int, int) = f2;

    if (s == 1) {
        for(line_of_blocks =  0; line_of_blocks <= k; ++line_of_blocks) {
            num_in_line1 = (line_of_blocks < k) ? m : l;
            tmp = line_of_blocks * elements_in_blcok_line;
            for(line_in_block = 0; line_in_block < num_in_line1; ++line_in_block) {
                i = line_of_blocks * m + line_in_block + 1;
                for(block_in_line = 0; block_in_line <= k; ++block_in_line) {
                    num_in_line2 = (block_in_line < k) ? m : l;
                    cur_pointer = matrix + tmp + block_in_line * m * num_in_line1
                        + line_in_block * num_in_line2;
                    j = block_in_line * m + 1;
                    for(q = 0; q < num_in_line2; ++q) {
                        *(cur_pointer + q) = n - ((i > j) ? i : j) + 1;
                        ++j;
                    }
                }
            }
        }
        return;
    }
    switch(s) {
        case 2: 
            p = f2;
            break;
        case 3:
            p = f3;
            break;
        case 4:
            p = f4;
            break;
        default:
            break;
    }
    for(line_of_blocks =  0; line_of_blocks <= k; ++line_of_blocks) {
        num_in_line1 = (line_of_blocks < k) ? m : l;
        tmp = line_of_blocks * elements_in_blcok_line;
        for(line_in_block = 0; line_in_block < num_in_line1; ++line_in_block) {
            i = line_of_blocks * m + line_in_block + 1;
            for(block_in_line = 0; block_in_line <= k; ++block_in_line) {
                num_in_line2 = (block_in_line < k) ? m : l;
                cur_pointer = matrix + tmp + block_in_line * m * num_in_line1
                    + line_in_block * num_in_line2;
                j = block_in_line * m + 1;
                for(q = 0; q < num_in_line2; ++q) {
                    *(cur_pointer + q) = p(i, j);
                    ++j;
                }
            }
        }
    }
}


void print_matrix(double *matrix, int n, int m, int r) {
    /* 
     * n = m * k + l
     */
    if (r > n) r = n;
    int num_in_line1; //amount of lines in block
    int num_in_line2; // size of lines in blcok
    double* cur_element;
    int k = n / m;
    int l = n - m * k;
    int elements_in_blcok_line = k * m * m + l * m;
    int line_of_blocks, block_in_line, line_in_block, i;

    int restrict_k = r / m;
    int restrict_l = r - m * restrict_k;
    int restrict_num_in_line1;
    int restrict_num_in_line2;

    for(line_of_blocks =  0; line_of_blocks <= restrict_k; ++line_of_blocks) {
        num_in_line1 = (line_of_blocks < k) ? m : l;
        restrict_num_in_line1 = (line_of_blocks < restrict_k) ? m : restrict_l;
        for(line_in_block = 0; line_in_block < restrict_num_in_line1; ++line_in_block) {
            for(block_in_line = 0; block_in_line <= restrict_k; ++block_in_line) {
                num_in_line2 = (block_in_line < k) ? m : l;
                restrict_num_in_line2 = (block_in_line < restrict_k) ? m : restrict_l;
                cur_element = matrix + line_of_blocks * elements_in_blcok_line + 
                    block_in_line * m * num_in_line1 + line_in_block * num_in_line2;
                for(i = 0; i < restrict_num_in_line2; ++i) {
                    printf("%lf ", *(cur_element + i));
                }
            }
            printf("\n");
        }
    }
}



