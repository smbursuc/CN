import numpy as np
import scipy.linalg as la


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
    return A, d


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

def get_L_matrix(A):
    n = A.shape[0]
    L = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                L[i, j] = 1
            elif i > j:
                L[i, j] = A[i, j]
    return L


def solve_equation_system(A, b):
    n = A.shape[0]
    A, D = compute_cholesky(A)
    L = get_L_matrix(A)
    LT = matrix_transpose(L)
    y = np.zeros(n)
    x = np.zeros(n)
    z = np.zeros(n)

    # Ly = b
    for i in range(n):
        s = 0
        for j in range(i):
            s = s + L[i, j] * y[j]
        y[i] = (b[i] - s) / L[i, i]

    # Dz = y
    for i in range(n):
        z[i] = y[i] / D[i]

    # LTx = z for x
    for i in range(n - 1, -1, -1):
        s = 0
        for j in range(i + 1, n):
            s = s + LT[i, j] * x[j]
        x[i] = (z[i] - s) / LT[i, i]
    

    return x


def solve_equation_system_numpy(A, b):
    LU, piv = la.lu_factor(A)
    x = la.lu_solve((LU, piv), b)
    return x
    

action = input("Introduceti 1 pentru a genera o matrice random, introduceti 2 pentru a folosi matricea predefinita:")
if action == "1":
    n = 150
    A = generate_random_matrix(n)
    A, D = compute_cholesky(A)
    print(A)
    print(D)
elif action == "2":
    A = np.array([[1, 2.5, 3], [2.5, 8.25, 15.5], [3, 15.5, 43]])
    print("Matricea initiala:\n",A)
    print("Descompunerea Cholesky a matricei predefinite este:\n")
    A_Cholesky, D = compute_cholesky(A)
    print("A:\n",A_Cholesky)
    print("D:\n",D)

    L = get_L_matrix(A_Cholesky)
    print("L\n",L)
    LT = matrix_transpose(L)
    print("LT:\n",LT)


    # proba pentru a vedea daca merge descompunerea Cholesky
    # D_matrix_form = np.zeros((3, 3))
    # for i in range(3):
    #     D_matrix_form[i, i] = D[i]
    # A_composed = matrix_multiplication(matrix_multiplication(L, D_matrix_form), LT)
    # print(A_composed)

    diagonal_matrix_determiant = np.prod(D)
    determinant = matrix_determinant(L) * diagonal_matrix_determiant * matrix_determinant(LT)

    A = np.array([[1, 2.5, 3], [2.5, 8.25, 15.5], [3, 15.5, 43]])
    print("Determinantul matricei predefinite fara descompunere este:\n",matrix_determinant(A))
    print("Determinantul matricei predefinite folosind descompunerea este:\n", determinant)

    print("Solutie pentru A si b = [12, 38, 68]:\n")
    print(solve_equation_system(A, np.array([12, 38, 68])))

    print("Solutie pentru A si b = [12, 38, 68] folosind numpy:\n")

    P, L, U = la.lu(A)
    print("Descompunerea LU a matricei predefinite este:\n")
    print("L:\n", L)
    print("U:\n", U)

    A = np.array([[1, 2.5, 3], [2.5, 8.25, 15.5], [3, 15.5, 43]])
    print("Solutia folosind numpy:\n", solve_equation_system_numpy(A, np.array([12, 38, 68])))

    A = np.array([[1, 2.5, 3], [2.5, 8.25, 15.5], [3, 15.5, 43]])
    print("Verificare calcul norma rezultat:\n")
    X_Cholesky = solve_equation_system(A, np.array([12, 38, 68]))
    A = np.array([[1, 2.5, 3], [2.5, 8.25, 15.5], [3, 15.5, 43]])
    print(la.norm(A.dot(X_Cholesky) - np.array([12, 38, 68])))

    # verify if the solution is correct
    print("Verificare calcul norma rezultat folosind numpy:\n")
    A = np.array([[1, 2.5, 3], [2.5, 8.25, 15.5], [3, 15.5, 43]])
    print(la.norm(A.dot(solve_equation_system_numpy(A, np.array([12, 38, 68]))) - np.array([12, 38, 68])))



else:
    print("Actiune invalida")
