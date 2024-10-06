import numpy as np


def generate_matrix(size, a, b):
    # Generate a random matrix
    matrix = np.random.randint(a, b, (size, size)).astype(np.float32)
    return matrix


def write_matrix_to_file(matrix, filename):
    np.savetxt(filename, matrix, fmt='%.6f')

def main():
    size = 10  # Size of the matrix (10x10)
    a = -10     # Lower bound
    b = 10    # Upper bound

    matrix = generate_matrix(size, a, b);
    write_matrix_to_file(matrix, "matrix.txt")
    matrix2 = matrix
    write_matrix_to_file(np.dot(matrix, matrix2), "resmatrix.txt");


if __name__ == "__main__":
    main()
