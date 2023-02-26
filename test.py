import numpy as np

def recursion(a):
    if a == 1:
        return 1
    else:
        b = recursion(a-1)
        print(b)
        b = 3

recursion(5)