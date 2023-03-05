import numpy as np

def matrix_determinant(A):
    if A.shape[0] == 2:
        return A[0, 0] * A[1,1] - A[0, 1] * A[1, 0]
    else:
        det = 0
        for i in range(A.shape[0]):
            det += A[0, i] * matrix_determinant(np.delete(np.delete(A, 0, 0), i, 1)) * (-1) ** i
        return det


def matrix_transpose(A):
    n = A.shape[0]
    m = A.shape[1]
    B = np.zeros((m, n))
    for i in range(n):
        for j in range(m):
            B[j, i] = A[i, j]
    return B

def matrix_multiplication(A, B):
    return A.dot(B)


A = np.matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
B = np.matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
print(matrix_determinant(A))
print(np.delete(A, 2, 0))

print(A)
print(matrix_transpose(A))
print(matrix_multiplication(A, B))

print(np.prod([1,2,3,4,5,6]))

print(A[:,0])
