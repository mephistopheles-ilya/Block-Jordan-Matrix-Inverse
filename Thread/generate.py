import numpy as np


def generate_matrix(size, a, b):
    # Generate a random matrix
    matrix = np.random.randint(a, b, (size, size))
    return matrix

def generate_matrix_non_zero_determinant(size, a, b):
    while True:
        matrix = np.random.normal(size/2, size/2, (size, size))
        #if np.linalg.det(matrix) != 0:
        if np.linalg.matrix_rank(matrix) == size:
            return matrix

def write_matrix_to_file(matrix, filename):
    np.savetxt(filename, matrix, fmt='%.6f')

def main():
    size = 2000  # Size of the matrix (10x10)
    a = -size     # Lower bound
    b = size   # Upper bound

    matrix = generate_matrix_non_zero_determinant(size, a, b)
    write_matrix_to_file(matrix, "matrix.txt");

    inverse = np.linalg.inv(matrix);
    write_matrix_to_file(inverse, "inverse.txt");
    
    prod = np.dot(matrix, inverse)  - np.eye(size);
    col_sums = np.abs(prod).sum(axis=0)
    max_sum = np.max(col_sums)
    print(max_sum);

    prod = np.dot(inverse, matrix)  - np.eye(size);
    col_sums = np.abs(prod).sum(axis=0)
    max_sum = np.max(col_sums)
    print(max_sum);

if __name__ == "__main__":
    main()
