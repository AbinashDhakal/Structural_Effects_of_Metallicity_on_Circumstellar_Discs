import sympy as smp
from sympy import erf
from scipy import special
from sympy import diff
import numpy as np

def f(t, x, **params):

    a = params['a']
    c = params['c']

    f1 = a * (x[0] - x[0] * x[1] +x[2])
    f2 = -c * (x[1] - x[0] * x[1] + x[2])  
    #f3 = a*(x[0] - x[0] * x[1] + x[2] )

    return np.array([f1, f2], dtype = np.float)

def df(t, x, **params):

    eps = 1e-10
    J = np.zeros([len(x), len(x)], dtype = np.float)

    for i in range(len(x)):
        x1 = x.copy()
        x2 = x.copy()

        x1[i] += eps
        x2[i] -= eps

        f1 = f(t, x1, **params)
        print(f1)
        f2 = f(t, x2, **params)
        print(f2)
        #f3 = f(t, x2, **params)


        J[ : , i] = (f1 - f2) / (2 * eps)

    return J

t = 0
x = np.array([1, 2,1], dtype = np.float)
print( df(t, x, a = 1, c = 1))