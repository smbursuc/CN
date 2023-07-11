import numpy as np
import math
def generate_I_matrix(n):
    I = np.zeros((n, n))
    for i in range(n):
        I[i][i] = 1
    return I


def householder(A, s):
    n = len(A)
    Q_ = generate_I_matrix(n)
    epsilon = 1e-6
    b = calculate_b(A, s)
    for r in range(n-1):
        sigma = 0
        for i in range(r, n):
            sigma += A[i][r]**2
        if sigma <= epsilon:
            break
        k = math.sqrt(sigma)
        if(A[r][r]>0):
            k = -k
        beta = sigma - A[r][r]*k
        u = np.zeros(n)
        u[r] = A[r][r] - k
        for i in range(r+1, n):
            u[i] = A[i][r]
        for j in range(r+1,n):
            gamma = 0
            for i in range(r, n):
                gamma += u[i]*A[i][j]
            gamma = gamma/beta
            for i in range(r, n):
                A[i][j] -= gamma*u[i]
        A[r][r] = k
        for i in range(r+1, n):
            A[i][r] = 0
        gamma = 0
        for i in range(r, n):
            gamma += u[i]*b[i]
        gamma = gamma/beta
        for i in range(r, n):
            b[i] -= gamma*u[i]
        for i in range(r,n):
            b[i] = b[i] - gamma*u[i]
        for i in range(n):
            gamma = 0
            for j in range(r, n):
                gamma += u[j]*Q_[i][j]
            gamma = gamma/beta
            for j in range(r, n):
                Q_[i][j] -= gamma*u[j]
    return A, Q_


def solve_RxQTb(R, Q, b):
    QTb = np.dot(Q.T, b)
    n = len(R)
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        s = 0
        for j in range(i + 1, n):
            s = s + R[i, j] * x[j]
        x[i] = (QTb[i] - s) / R[i, i]
    return x

def solve_upper_triangular(A, b):
    n = len(A)
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        s = 0
        for j in range(i + 1, n):
            s = s + A[i, j] * x[j]
        x[i] = (b[i] - s) / A[i, i]
    return x


def inverse_of_matrix_householder(Q,R):
    A_inverse_columns = np.zeros((len(Q), len(Q)))
    for i in range(len(Q)):
        b = Q[i]
        x = solve_upper_triangular(R, b)
        A_inverse_columns[:, i] = x
    return A_inverse_columns

def calculate_b(A,s):
    b = np.zeros(len(A))
    for i in range(len(A)):
        for j in range(len(A)):
            b[i] += A[i][j]*s[j]
    return b

def test_with_n(n):
    A = np.random.rand(n, n)
    for i in range(n):
        if(A[i][i] == 0):
            break
    A_copy = A.copy()
    #print(A,"\n")
    s = np.random.rand(n)
    b = calculate_b(A,s)
    R, Q = householder(A, b)
    x_householder = solve_RxQTb(R, Q, b)
    Q, R = np.linalg.qr(A_copy)
    x_qr = solve_RxQTb(R, Q, b)
    print("Norma verificare solutie propie vs numpy\n", np.linalg.norm(x_qr - x_householder))
    print("\nNormele cerute la punctul 4 din tema:")
    print("\n|| A_init dot x_householder - b_init ||\n",np.linalg.norm(A.dot(x_householder)-b))

    print("\n|| A_init dot x_QR - b_init ||\n",np.linalg.norm(A.dot(x_qr)-b))

    print("\n|| x_householder - s || / || s ||\n",np.linalg.norm(x_householder-s)/np.linalg.norm(s))

    print("\n|| x_QR - s || / || s ||\n",np.linalg.norm(x_qr-s)/np.linalg.norm(s))


action = input("Pentru a testa random tasta 1, pentru cazul exemplu din pdf tasta 2")
if(action == "1"):
    for n in range(10,100,10):
        print("Test pentru matrice cu n=",n,"\n")
        test_with_n(n)
else:
    print("\nA initial:\n",np.array([[0, 0, 4], [1, 2, 3], [0, 1, 2]]))

    R, Q = householder(np.array([[0, 0, 4], [1, 2, 3], [0, 1, 2]]),np.array([3, 2, 1]))
    print("\nQ folosind householder:\n",Q)
    print("\nR folosind householder:\n",R)
    A = np.array([[0, 0, 4], [1, 2, 3], [0, 1, 2]])
    #b = np.dot(A, np.array([3, 2, 1]))
    b = calculate_b(A,np.array([3, 2, 1]))
    x_householder = solve_RxQTb(R, Q, b)
    print("\nSolutie folosind Householder:\n",x_householder)

    # QR decomposition using numpy
    Q, R = np.linalg.qr(A)
    print("\nQ folosind numpy:\n",Q)
    print("\nR folosind numpy:\n",R)
    x_qr = solve_RxQTb(R, Q, b)
    print("\nSolutie folosind numpy:\n",x_qr)
    print("\nNorma verificare solutie propie vs numpy\n",np.linalg.norm(x_qr-x_householder))


    R, Q = householder(np.array([[0, 0, 4], [1, 2, 3], [0, 1, 2]]),np.array([3, 2, 1]))
    print("\nInversa lui A folosing descompunerea Householder\n",inverse_of_matrix_householder(Q,R))

    A_inverse = np.linalg.inv(A)
    print("\nInversa lui A folosind numpy\n",A_inverse)

    print("\nNorma verificare inversa proprie vs numpy\n",np.linalg.norm(A_inverse-inverse_of_matrix_householder(Q,R)))

    print("\nNormele cerute la punctul 4 din tema:")
    A = np.array([[0, 0, 4], [1, 2, 3], [0, 1, 2]])
    print("\n|| A_init dot x_householder - b_init ||\n",np.linalg.norm(A.dot(x_householder)-b))

    print("\n|| A_init dot x_QR - b_init ||\n",np.linalg.norm(A.dot(x_qr)-b))

    print("\n|| x_householder - s || / || s ||\n",np.linalg.norm(x_householder-np.array([3, 2, 1]))/np.linalg.norm(np.array([3, 2, 1])))

    print("\n|| x_QR - s || / || s ||\n",np.linalg.norm(x_householder-np.array([3, 2, 1]))/np.linalg.norm(np.array([3, 2, 1])))







        

        
        
