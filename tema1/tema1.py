import math
import pprint
import numpy as np

def ex1():
    m = 0
    while math.pow(10, -m) + 1 != 1:
        m += 1
    m-=1
    return math.pow(10, -m)

print("u = ",ex1())
print()

def ex2():
    x = 1.0
    u = ex1()
    print("(x + u) + u ------> ",(x + u) + u)
    print("x + (u + u) ------> ",x + (u + u))
    print("(x + u) + u == x + (u + u) ------> ",(x + u) + u == x + (u + u))
    print()
    x = 0.1
    y = 0.2
    z = 0.3
    print("(x * u) * u ------> ",(x * y) * z)
    print("x * (u * u) ------> ",x * (y * z))
    print("(x * u) * u == x * (u * u) ------> ",(x * y) * z == x * (y * z))

ex2()

print()

def matrix_addition(A, B):
    n = len(A)
    C = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            C[i][j] = A[i][j] + B[i][j]
    return C

def strassen(A, B, n):
    if n==2:
        P1 = (A[0][0] + A[1][1]) * (B[0][0] + B[1][1])
        P2 = (A[1][0] + A[1][1]) * B[0][0]
        P3 = A[0][0] * (B[0][1] - B[1][1])
        P4 = A[1][1] * (B[1][0] - B[0][0])
        P5 = (A[0][0] + A[0][1]) * B[1][1]
        P6 = (A[1][0] - A[0][0]) * (B[0][0] + B[0][1])
        P7 = (A[0][1] - A[1][1]) * (B[1][0] + B[1][1])
        C = [[0, 0], [0, 0]]
        C[0][0] = P1 + P4 - P5 + P7
        C[0][1] = P3 + P5
        C[1][0] = P2 + P4
        C[1][1] = P1 + P3 - P2 + P6
        return C
    else:
        middle = int(n/2)
        A11 = [[A[i][j] for j in range(middle)] for i in range(middle)]
        A12 = [[A[i][j] for j in range(middle, n)] for i in range(middle)] 
        A21 = [[A[i][j] for j in range(middle)] for i in range(middle, n)] 
        A22 = [[A[i][j] for j in range(middle, n)] for i in range(middle, n)] 
        B11 = [[B[i][j] for j in range(middle)] for i in range(middle)]
        B12 = [[B[i][j] for j in range(middle, n)] for i in range(middle)]
        B21 = [[B[i][j] for j in range(middle)] for i in range(middle, n)] 
        B22 = [[B[i][j] for j in range(middle, n)] for i in range(middle, n)]

        # pprint.pprint(A11)
        # pprint.pprint(A12)
        # pprint.pprint(A21)
        # pprint.pprint(A22)
        # pprint.pprint(B11)
        # pprint.pprint(B12)
        # pprint.pprint(B21)
        # pprint.pprint(B22)
        # return 0
        
        C11 = matrix_addition(strassen(A11, B11, middle),strassen(A12, B21, middle))
        C12 = matrix_addition(strassen(A11, B12, middle),strassen(A12, B22, middle))
        C21 = matrix_addition(strassen(A21, B11, middle),strassen(A22, B21, middle))
        C22 = matrix_addition(strassen(A21, B12, middle),strassen(A22, B22, middle))

        C = [[C11,C12],[C21,C22]]
        return C
    
# A = [[1,2],[3,4]]
# B = [[5,6],[7,8]]

# print(strassen(A,B,2))

A = [[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]]
B = [[4,3,2,1],[8,7,6,5],[12,11,10,9],[16,15,14,13]]

result = strassen(A,B,4)
pprint.pprint(result)
