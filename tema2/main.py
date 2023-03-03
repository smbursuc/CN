import numpy as np


def generate_random_matrix(n):
    while True:
        A = np.random.rand(n, n)
        if np.linalg.matrix_rank(A) == n:
            break
    return A @ A.T


def compute_cholesky(A):
    A_init = np.copy(A)
    n = A.shape[0]
    d = np.zeros(n)

    for p in range(n):
        # calculeaza elementul dp
        dp = A_init[p, p]
        for k in range(p):
            dp = dp - d[k] * A[p, k] ** 2
        d[p] = dp

        # calculeaza elementele coloanei p
        for i in range(p + 1, n):
            s = 0
            for k in range(p):
                s = s + d[k] * A_init[i, k] * A[p, k]
            A[i, p] = (A_init[i, p] - s) / dp
    print("A:\n", A)
    print("d:\n", d)


action = input("Introduceti 1 pentru a genera o matrice random, introduceti 2 pentru a folosi matricea predefinita:")
if action == "1":
    n = 150
    A = generate_random_matrix(n)
    compute_cholesky(A)
elif action == "2":
    A = np.array([[1, 2.5, 3], [2.5, 8.25, 15.5], [3, 15.5, 43]])
    compute_cholesky(A)
else:
    print("Actiune invalida")
