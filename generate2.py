import numpy as np

def f(i, j):
    return (j - i) if (i - j) < 0 else (i - j)

def generate_matrix(size, a, b):
    matrix = np.random.uniform(a, b, (size, size))
    '''
    matrix = np.zeros((size, size))
    for i in range(size):
        for j in range(size):
            matrix[i, j] = f(i, j)
    '''
    return matrix

def generate_matrix_non_zero_determinant(size, a, b):
    while True:
        matrix = np.random.normal(size/2, size/2, (size, size))
        #if np.linalg.det(matrix) != 0:
        if np.linalg.matrix_rank(matrix) == size:
            return matrix

def write_matrix_to_file(matrix, filename):
    #np.savetxt(filename, matrix, fmt='%.6f')
    with open(filename, 'w') as f:
        for row in matrix:
            for element in row:
                if abs(element) < 1e-5:
                    element = 0
                f.write(f"{element:.6f}\n")  # Write each element on a new line


def main():
    size = 2000  # Size of the matrix (10x10)
    a = -size     # Lower bound
    b = size   # Upper bound

    matrix = generate_matrix(size, a, b)
    write_matrix_to_file(matrix, "matrix.txt");
    inverse = np.linalg.inv(matrix);
    write_matrix_to_file(inverse, "inverse.txt");
    
'''
    prod = np.dot(matrix, inverse)  - np.eye(size);
    col_sums = np.abs(prod).sum(axis=0)
    max_sum = np.max(col_sums)
    print(max_sum);

    prod = np.dot(inverse, matrix)  - np.eye(size);
    col_sums = np.abs(prod).sum(axis=0)
    max_sum = np.max(col_sums)
    print(max_sum);
'''

if __name__ == "__main__":
    main()
