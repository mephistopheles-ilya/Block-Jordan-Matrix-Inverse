#include "func.hpp"
#include <cstring>
#include <cstdio>

int read_array_from_file(double* array, int n, int p, int k, MPI_Comm comm, char* filename) {
	int err = 0;
	FILE* f = nullptr;
	if (k == 0) {
		f = fopen(filename, "r");
		if (f == nullptr) {
			err = 1;
		}
	}
	MPI_Bcast(&err, 1, MPI_INT, 0, comm);
	if (err != 0) {
		return err;
	}
	memset(array, 0, n * sizeof(double));
	int i = 0;
	if (k == 0) {
		for(i = 0; i < n; ++i) {
			if (fscanf(f, "%lf", array + i) != 1) {
				i = -1;
				break;
			}
		}
	}
	if (i < 0) {
		err = 2;
	}
	MPI_Bcast(&err, 1, MPI_INT, 0, comm);
	if (err != 0) {
		return err;
	}
	MPI_Bcast(array, n, MPI_DOUBLE, 0, comm);
	if (k == 0) {
		fclose(f);
	}
	return 0;
}

void get_indexes(int p, int k, int n, int& l1, int& l2) {
	int chank = n/p;
	if (chank == 0) {
		chank = 1;
	}
	l1 = k * chank;
	l2 = (k + 1) * chank;
	if (k == (p - 1) && l2 < n) {
		l2 = n;
	}
	if (l1 >= n) {
		l1 = -1;
		l2 = -1;
	}

}


void do_my_task(double* array, int n, int p, int k, MPI_Comm comm) {
	int l1 = 0, l2 = 0;
	get_indexes(p, k, n, l1, l2);
	if (l1 >= 0 && l2 >= 0) {
		for(int i = 0; i < n; ++i) {
			array[i] = k;
		}
	}
}
