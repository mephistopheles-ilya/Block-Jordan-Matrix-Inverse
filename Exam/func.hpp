#pragma once 

#include "mpi.h"

int read_array_from_file(double* array, int n, int p, int k, MPI_Comm comm, char* filame);

void do_my_task(double* array, int n, int p, int k, MPI_Comm comm);

void get_indexes(int p, int k, int n, int& l1, int& l2);
