#include "func.hpp"
#include <stdio.h>

int main(int argc, char* argv[]) {
	int p = 0, k = 0;
	int err_loc = 0, err_glob = 0;
	int n = 0;
	char* filename = nullptr;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Init(&argc, &argv);

	MPI_Comm_size(comm, &p);
	MPI_Comm_rank(comm, &k);

	if (argc != 3) {
		err_loc = 1;
	} else if (sscanf(argv[1], "%d", &n) != 1) {
		err_loc = 1;
	} else if (n <= 0) {
		err_loc = 1;
	}

	MPI_Allreduce(&err_loc, &err_glob, 1, MPI_INT, MPI_MAX, comm);
	if (err_glob != 0) {
		if (k == 0) {
			printf("Wring command line argumants\n");
		}
		MPI_Finalize();
		return 0;
	}
	filename = argv[2];

	double* array = new (std::nothrow) double[n];
	if (array == nullptr) {
		err_loc = 1;
	}
	MPI_Allreduce(&err_loc, &err_glob, 1, MPI_INT, MPI_MAX, comm);
	if (err_glob != 0) {
		if (k == 0) {
			printf("Not enough memmory\n");
		}
		MPI_Finalize();
		return 0;
	}

	err_loc = read_array_from_file(array, n, p, k, comm, filename);
	MPI_Allreduce(&err_loc, &err_glob, 1, MPI_INT, MPI_MAX, comm);
	if (err_glob != 0) {
		if (k == 0) {
			printf("Errors with file\n");
		}
		MPI_Finalize();
		return 0;
	}
	if (k == 0) {
		printf(" %d", n);
		for(int i = 0; i < n; ++i) {
			printf(" %lf", array[i]);
		}
		printf("\n");
	}
	
	double time = MPI_Wtime();
	do_my_task(array, n, p, k, comm);
	time = MPI_Wtime() - time;

	int src, dst, l1, l2;
	MPI_Status st;
	for(int i = 1; i < p; ++i) {
		src = i;
		dst = 0;
		get_indexes(p, src, n, l1, l2);
		if (l1 >= 0 && l2 >= 0) {
			if (k == 0) {
				MPI_Recv(array + l1, l2 - l1, MPI_DOUBLE, src, 0, comm, &st);
			} else if (k == src) {
				MPI_Send(array + l1, l2 - l1, MPI_DOUBLE, dst, 0, comm);
			}
		}
	}
		
	if (k == 0) {
		for(int i = 0; i < n; ++i) {
			printf("%lf ", array[i]);
		}
		printf("\n");
	}

	double all_time = 0;
	MPI_Reduce(&time, &all_time, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
	if (k == 0) {
		printf("%s : %lf\n", "All time ", all_time);
	}

	MPI_Gather(&time, sizeof(double), MPI_BYTE, array, sizeof(double), MPI_BYTE, 0, comm);

	if (k == 0) {
		for(int i = 0; i < p; ++i) {
			printf("Time for process %d : %.10lf\n", i, array[i]);
		}
	}
	
	delete[] array;	
	MPI_Finalize();
	return 0;
}
