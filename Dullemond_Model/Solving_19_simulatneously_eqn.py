import numpy as np
from scipy.optimize import fsolve, minimize
np.set_printoptions(suppress=True)

def f(x) :
    y = np.zeros(np.size(x))
    y[0] = -29.79 + 1.795*x[1] + 6.33*x[0] + -0.305*x[3] - 0.024*x[1]*x[3] + 0.222*x[0]*x[3] -2.063*x[4] + 0.122*x[1]*x[4] + 0.864*x[0]*x[4] + x[5] - x[6]
    y[1] = -23.269 + 1.795*x[0] + 4.547*x[2] - 0.342*x[1] - 0.024*x[0]*x[3] + 0.072*x[2]*x[3] - 0.052*x[1]*x[3] - 1.545*x[4] + 0.122*x[0]*x[4] + 0.367*x[2]*x[4] + 0.114*x[1]*x[4] + x[7] - x[8]
    y[2] = 0.587 + 4.547*x[1] + 0.046*x[3] + 0.072*x[1]*x[3] - 0.211*x[4] + 0.367*x[1]*x[4] + x[9] - x[10]
    y[3] = -1.657 - 0.305*x[0] - 0.230*x[1] + 0.046*x[2] - 0.024*x[0]*x[1] + 0.072*x[1]*x[2] + 0.111*x[0]**2.0 - 0.026*x[1]**2.0 + 0.78981 + x[11]**2.0
    y[4] = -0.106 - 2.063*x[0] - 1.545*x[1] - 0.211*x[2] + 0.122*x[0]*x[1] + 0.367*x[1]*x[2] + 0.432*x[0]**2.0 + 0.057*x[1]**2.0 + 1.34744 + x[12]**2.0
    y[5] = x[0] - 0.3 + x[13]**2.0
    y[6] = -x[0] + 0.1 + x[14]**2.0
    y[7] = x[1] - 80 + x[15]**2.0
    y[8] = -x[1] + 30 + x[16]**2.0
    y[9] = x[2] - 230 + x[17]**2.0
    y[10] = -x[2] + 200 + x[18]**2.0
    y[11] = 2*x[3]*x[11]
    y[12] = 2*x[4]*x[12]
    y[13] = 2*x[5]*x[13]
    y[14] = 2*x[6]*x[14]
    y[15] = 2*x[7]*x[15]
    y[16] = 2*x[8]*x[16]
    y[17] = 2*x[9]*x[17]
    y[18] = 2*x[10]*x[18]
    return y

x0 = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
x =fsolve(f, x0)

print('fsolve')
print ("x:")
print(x)
print ("|x|_2:")
print(np.linalg.norm(x))

""" Try more starting-points fsolve """
N = 1000
best_obj = np.inf
best_x = 0
for i in range(N):
    x0 = np.random.uniform(size=x0.shape)
    x = fsolve(f, x0)
    fun = np.linalg.norm(x)
    if fun < best_obj:
        best_obj = fun
        best_x = x

print('-----')
print('multi-start fsolve')
print('x: ')
print(best_x)
print('|x|_2:')
print(np.linalg.norm(best_x))