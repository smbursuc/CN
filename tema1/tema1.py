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
    return np.add(A, B)

def matrix_subtraction(A, B):
    return np.subtract(A, B)

def strassen(A, B, n):
    if n==2:
        P1 = (A[0,0] + A[1,1]) * (B[0,0] + B[1,1])
        P2 = (A[1,0] + A[1,1]) * B[0,0]
        P3 = A[0,0] * (B[0,1] - B[1,1])
        P4 = A[1,1] * (B[1,0] - B[0,0])
        P5 = (A[0,0] + A[0,1]) * B[1,1]
        P6 = (A[1,0] - A[0,0]) * (B[0,0] + B[0,1])
        P7 = (A[0,1] - A[1,1]) * (B[1,0] + B[1,1])
        C = np.matrix(np.zeros((2,2)))
        C[0,0] = P1 + P4 - P5 + P7
        C[0,1] = P3 + P5
        C[1,0] = P2 + P4
        C[1,1] = P1 + P3 - P2 + P6
        return C
    else:
        middle = int(n/2)
        # print(A)
        # print(B)
        A11 = np.matrix([[A[i,j] for j in range(middle)] for i in range(middle)])
        A12 = np.matrix([[A[i,j] for j in range(middle, n)] for i in range(middle)]) 
        A21 = np.matrix([[A[i,j] for j in range(middle)] for i in range(middle, n)]) 
        A22 = np.matrix([[A[i,j] for j in range(middle, n)] for i in range(middle, n)]) 
        B11 = np.matrix([[B[i,j] for j in range(middle)] for i in range(middle)])
        B12 = np.matrix([[B[i,j] for j in range(middle, n)] for i in range(middle)])
        B21 = np.matrix([[B[i,j] for j in range(middle)] for i in range(middle, n)]) 
        B22 = np.matrix([[B[i,j] for j in range(middle, n)] for i in range(middle, n)])

        # pprint.pprint(A11)
        # pprint.pprint(A12)
        # pprint.pprint(A21)
        # pprint.pprint(A22)
        # pprint.pprint(B11)
        # pprint.pprint(B12)
        # pprint.pprint(B21)
        # pprint.pprint(B22)
        # return 0

        # print(A11)
        # print(A22)
        # print(B11)
        # print(B22)
        P1 = strassen(matrix_addition(A11,A22), matrix_addition(B11,B22), middle)
        P2 = strassen(matrix_addition(A21,A22), B11, middle)
        P3 = strassen(A11, matrix_subtraction(B12,B22), middle)
        P4 = strassen(A22, matrix_subtraction(B21,B11), middle)
        P5 = strassen(matrix_addition(A11,A12), B22, middle)
        P6 = strassen(matrix_subtraction(A21,A11), matrix_addition(B11,B12), middle)
        P7 = strassen(matrix_subtraction(A12,A22), matrix_addition(B21,B22), middle)
        
        # C11 = matrix_addition(strassen(A11, B11, middle),strassen(A12, B21, middle))
        # C12 = matrix_addition(strassen(A11, B12, middle),strassen(A12, B22, middle))
        # C21 = matrix_addition(strassen(A21, B11, middle),strassen(A22, B21, middle))
        # C22 = matrix_addition(strassen(A21, B12, middle),strassen(A22, B22, middle))

        C11 = matrix_addition(matrix_subtraction(matrix_addition(P1,P4),P5),P7)
        C12 = matrix_addition(P3,P5)
        C21 = matrix_addition(P2,P4)
        C22 = matrix_addition(matrix_subtraction(matrix_addition(P1,P3),P2),P6)

        C = np.concatenate((np.concatenate((C11,C12),axis=1),np.concatenate((C21,C22),axis=1)),axis=0)
        return C
    
# A = [[1,2],[3,4]]
# B = [[5,6],[7,8]]

# print(strassen(A,B,2))

A = [[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]]
B = [[4,3,2,1],[8,7,6,5],[12,11,10,9],[16,15,14,13]]

A = np.matrix(A)
B = np.matrix(B)

# result = strassen(A,B,4)
# pprint.pprint(result)

A = [[1,2,3,4,5,6,7,8],[9,10,11,12,13,14,15,16],[17,18,19,20,21,22,23,24],[25,26,27,28,29,30,31,32],[33,34,35,36,37,38,39,40],[41,42,43,44,45,46,47,48],[49,50,51,52,53,54,55,56],[57,58,59,60,61,62,63,64]]
B = [[1,2,3,4,5,6,7,8],[9,10,11,12,13,14,15,16],[17,18,19,20,21,22,23,24],[25,26,27,28,29,30,31,32],[33,34,35,36,37,38,39,40],[41,42,43,44,45,46,47,48],[49,50,51,52,53,54,55,56],[57,58,59,60,61,62,63,64]]

A = np.matrix(A)
B = np.matrix(B)


result = strassen(A,B,8)
pprint.pprint(result)


