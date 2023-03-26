import numpy as np
import pprint

def generate_random_sparse_matrix(n):
    A = []
    for i in range(n):
        A.append([])
    visited = []
    for i in range(n):
        for j in range(n):
            if (i,j) not in visited:
                if i==j:
                    val = round(np.random.rand(),2) * 100
                    A[i].append((val, j))
                else:
                    chance = np.random.rand()
                    if(chance<0.30):
                        val = round(np.random.rand(),2) * 100
                        A[i].append((val, j))
                        A[j].append((val, i))
                visited.append((j,i))
    return A



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
            line_vals.insert(j, tupla)
            last_i = i
        else:
            tupla = (val, j)
            line_vals.insert(j, tupla)
        x+=1
    f.close()
    return A


def is_symmetric(A):
    for i in range(len(A)):
        for j in range(len(A[i])):
            index = A[i][j][1]
            if i!=index:
                #fac rounding pt ca in fisier nu a dat exact egale valorile
                val1 = round(A[i][j][0], 4)
                #print(val1)
                for k in range(len(A[index])):
                    if A[index][k][1] == i:
                        val2 = round(A[index][k][0], 4)
                        if val1 != val2:
                            return False
                        break
    return True

def parse_files():
    # 2023 nu merge?
    ends = ["512","1024"]
    base_name = "m_rar_sim_2023_"
    for end in ends:
        file_name = base_name + end + ".txt"
        A = parse_sparse_matrix_file(file_name)
        print("E simetrica: " + str(is_symmetric(A)))
        print("Power method:")
        eingenvalue, eigenvector = power_method(A, len(A), 10000)
        print(eingenvalue, eigenvector[:5])
        print("Singular value decomposition:")
        singular_value_decomposition(A)

def sparse_matrix_multiplication(A, x):
    y = np.zeros(len(x))
    for i in range(len(A)):
        sum = 0.0
        for j in range(len(A[i])):
            sum += A[i][j][0]*x[A[i][j][1]]
        y[i] = sum
    return y


def power_method(A, n, kmax):
    v = np.random.random(n)
    # normalizare vector asa incat norma lui sa fie 1
    v = v/np.linalg.norm(v)
    w = sparse_matrix_multiplication(A, v)
    lambda_ = np.dot(w, v)
    k = 0
    while True:
        k += 1
        v = w/np.linalg.norm(w)
        w = sparse_matrix_multiplication(A, v)
        if(np.linalg.norm(w-np.dot(v,lambda_) <= n*10**(-9)) or k>kmax):
            break
        lambda_ = np.dot(w, v)
    return lambda_, v

def sparse_to_dense(A):
    n = len(A)
    dense = np.zeros((n,n))
    for i in range(n):
        for j in range(len(A[i])):
            dense[i][A[i][j][1]] = A[i][j][0]
    return dense

def singular_value_decomposition(A):
    A_dense = sparse_to_dense(A)
    #print("A_dense:\n",A_dense)
    u, s, vt = np.linalg.svd(A_dense, full_matrices=True)
    print("Singular values:",s)
    rank = np.sum(s > 10**(-9))
    print("Rank:",rank)
    condition_number = np.max(s)/np.min(s)
    print("Conditioning number:",condition_number)
    v = vt.T
    s_matrix_form = np.zeros((len(A_dense), len(A_dense)))
    for i in range(len(s)):
        s_matrix_form[i][i] = s[i]
    mp_inverse = v.dot(s_matrix_form).dot(u.T)
    print("Moore-Penroe inverse",mp_inverse)
    b = np.random.random(len(mp_inverse))
    x_i = np.linalg.solve(mp_inverse, b)
    print("x_i:",x_i)
    print("||b-Ax||:",np.linalg.norm(b-np.dot(A_dense, x_i)))
    least_squares_pseudo_inverse = np.dot(np.linalg.inv(np.dot(A_dense.T, A_dense)), A_dense.T)
    print("Least squares pseudo inverse:",least_squares_pseudo_inverse)
    print("||A_I - A_J||:",np.linalg.norm(mp_inverse-least_squares_pseudo_inverse))


    


action = input("Tasta 1 pentru a genera random, tasta 2 pentru a parsa fisierele\n")
if action == "1":
    n = 100
    random_matrix = generate_random_sparse_matrix(n)
    pp = pprint.PrettyPrinter(width=100, compact=True)
    pp.pprint(random_matrix)
    print("E simetrica: " + str(is_symmetric(random_matrix)))
    print("Power method:")
    eingenvalue, eigenvector = power_method(random_matrix, n, 10000)
    print(eingenvalue, eigenvector[:5])
    print("Singular value decomposition:")
    singular_value_decomposition(random_matrix)
else:
    parse_files()