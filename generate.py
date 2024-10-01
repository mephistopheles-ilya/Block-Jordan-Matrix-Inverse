import numpy as np

def generate_non_degenerate_matrix(size, a, b):
    while True:
        # Generate a random matrix
        matrix = np.random.uniform(a, b, (size, size))
        # matrix = np.random.uniform(a, b, (size, size)).astype(np.float128)
        # Check if the determinant is non-zero
        #if np.linalg.matrix_rank(matrix) == size:
        det  = np.linalg.det(matrix);
        if det != 0:
            print("det = ", det)
            return matrix

def generate_matrix(size, a, b):
    # Generate a random matrix
    matrix = np.random.randint(a, b, (size, size))
    return matrix


def write_matrix_to_file(matrix, filename, index):
    np.savetxt(filename, matrix[:index, :index], fmt='%.6f')

def main():
    size = 10  # Size of the matrix (10x10)
    a = 0     # Lower bound
    b = 10    # Upper bound

    matrix = generate_matrix(size, a, b)

    for i in range(1,size + 1):
        filename = f"matrices/matrix_{i}.txt"
        write_matrix_to_file(matrix, filename, i)


if __name__ == "__main__":
    main()
