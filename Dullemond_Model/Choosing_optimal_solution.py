from scipy.optimize import fsolve
import numpy as np
import random

def eq(x):
     return x[0] + x[1]**2 - 4 , np.sin(x[0]) + x[0]*x[1] - 3
 
x0 = np.array([1,1])
z= fsolve(eq, x0)
print(fsolve(eq,x0))
print(np.linalg.norm(z))

N = 1000
best_obj = np.inf
bext_z =0
for i in range(N):
    x0 = np.random.uniform(size=x0.shape)
    z= fsolve(eq, x0)
    print(z)
    fun = np.linalg.norm(z)
    if fun < best_obj:
        best_obj = fun
        best_z = z

print('-----')
print('multi-start fsolve')
print('x: ')
print(best_z)
print('|x|_2:')
print(np.linalg.norm(best_z))