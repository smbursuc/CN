import numpy as np
import math
import pprint

def parse_sparse_matrix_file(file):
    f = open(file, 'r')
    lines = f.readlines()
    n = int(lines[0].rsplit('\n')[0])
    A = []
    lines.remove(lines[0])
    last_i = int(lines[0].replace('\n', '').split(",")[1])
    x = 0
    line_vals = []
    while x <= len(lines):
        if(x == len(lines)):
            A.append(line_vals.copy())
            break
        line = lines[x].replace('\n', '')
        line = line.split(",")
        val = float(line[0])
        i = int(line[1])
        j = int(line[2])
        if(i==j):
            if(val==0):
                print("0 pe diagonala")
                break
        if i!=last_i:
            A.append(line_vals.copy())
            line_vals.clear()
            tupla = (val, j)
            line_vals.append(tupla)
            last_i = i
        else:
            tupla = (val, j)
            line_vals.append(tupla)
        x+=1
    f.close()
    return A


def sparse_matrix_multiplication(A, x):
    y = np.zeros(len(x))
    for i in range(len(A)):
        sum = 0.0
        for j in range(len(A[i])):
            sum += A[i][j][0]*x[A[i][j][1]]
        y[i] = sum
    return y


def gauss_seidel_2_vectors(A, b, kmax):
    #x_c = np.zeros(len(b))
    x_gs = np.zeros(len(b))
    # x_c = x_k+1 si x_p = x_k 
    k = 0
    epsilon = 1e-10
    delta_x = 0
    while True:
        # compute new x_c using x_p
        for i in range(len(A)):
            sum = 0.0
            elem_diag = 0
            x_gs_copy = x_gs.copy()
            for j in range(len(A[i])):
                if(A[i][j][1]==i):
                    elem_diag = A[i][j][0]
            for j in range(len(A[i])):
                if(A[i][j][0]!=elem_diag):
                    sum += A[i][j][0]*x_gs_copy[A[i][j][1]]
                    #print(A[i][j][0], x_p[A[i][j][1]], elem_diag)
            x_gs[i] = (b[i]-sum)/elem_diag
        # compute delta_x
        delta_x = np.linalg.norm(x_gs-x_gs_copy)
        k+=1
        print("Pas: ",k)
        if(delta_x < epsilon or delta_x>math.pow(10,8) or k > kmax):
            break
        # if(k > kmax):
        #     break
    if(delta_x<epsilon):
        print("x:",x_gs)
        print("k: ",k)
    else:
        print("Divergence")
    print("Norma ceruta:",np.linalg.norm(sparse_matrix_multiplication(A, x_gs)-b))

def gauss_seidel(A, b, kmax):
    #x_c = np.zeros(len(b))
    x_gs = np.zeros(len(b))
    # x_c = x_k+1 si x_p = x_k 
    k = 0
    epsilon = 1e-10
    delta_x = 0
    while True:
        # compute new x_c using x_p
        norm_partial = 0
        for i in range(len(A)):
            aux_prev_value = x_gs[i]
            sum = 0.0
            elem_diag = 0
            for j in range(len(A[i])):
                if(A[i][j][1]==i):
                    elem_diag = A[i][j][0]
            for j in range(len(A[i])):
                if(A[i][j][0]!=elem_diag):
                    sum += A[i][j][0]*x_gs[A[i][j][1]]
                    #print(A[i][j][0], x_p[A[i][j][1]], elem_diag)
            x_gs[i] = (b[i]-sum)/elem_diag

            norm_partial += (x_gs[i] - aux_prev_value) ** 2
        norm_partial = math.sqrt(norm_partial)
        # compute delta_x
        #delta_x = np.linalg.norm(x_gs-x_gs_copy)
        k+=1
        print("Pas: ",k)
        if(norm_partial < epsilon or norm_partial>math.pow(10,8) or k > kmax):
            break
        # if(k > kmax):
        #     break
    if(norm_partial<epsilon):
        print("x:",x_gs)
        print("k: ",k)
    else:
        print("Divergence")
    print("Norma ceruta:",np.linalg.norm(sparse_matrix_multiplication(A, x_gs)-b))

def parse_b_file(file):
    f = open(file, 'r')
    lines = f.readlines()
    b = np.zeros(int(lines[0].replace('\n', '')))
    lines.remove(lines[0])
    for i in range(len(lines)):
        b[i] = float(lines[i].replace('\n', ''))
    f.close()
    return b


action = input("Tasta 1 pentru exemplu pdf tasta 2 pentru fisierele date\n")

if action=='1':
    A = parse_sparse_matrix_file('fisiere/a_example.txt')
    print(pprint.pprint(A))

    b = [6.0, 7.0, 8.0, 9.0, 1.0]
    gauss_seidel(A, b, 10)
else:
    root = "fisiere/"
    iteratii = 100
    for i in range(1, 6):
        print("Exemplul a_" + str(i) + " b_" + str(i))
        A = parse_sparse_matrix_file(root+'a_'+str(i)+'.txt')
        b = parse_b_file(root+'b_'+str(i)+'.txt')
        gauss_seidel(A, b, iteratii)
        print("Solutie PDF:")
        if(i==1):
            sol = np.zeros(len(b))
            for j in range(len(b)):
                sol[j] = 1.0
            print("x_i = 1 oricare ar fi i")
            print("solutia: ",sol)
        elif i==2:
            sol = np.zeros(len(b))
            for j in range(len(b)):
                sol[j] = 4.0/3.0
            print("x_i = 4.0/3.0 oricare ar fi i")
            print("solutia: ",sol)
        elif i==3:
            sol = np.zeros(len(b))
            for j in range(len(b)):
                sol[j] = 0.4*(j+1)
            print("x_i = 0.4*(i+1) oricare ar fi i")
            print("solutia: ",sol)
        elif i==4:
            sol = np.zeros(len(b))
            for j in range(len(b)):
                sol[j] = j/7.0
            print("x_i = i/7.0 oricare ar fi i")
            print("solutia: ",sol)
        elif i==5:
            sol = np.zeros(len(b))
            for j in range(len(b)):
                sol[j] = 2.0
            print("x_i = 2.0 oricare ar fi i")
            print("solutia: ",sol)
        print("")