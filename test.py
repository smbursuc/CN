import numpy as np

A = np.matrix([[1, 2], [3, 4]])
B = np.matrix([[5, 6], [7, 8]])

print(np.add(A, B))

print(np.concatenate((A, B), axis=0))