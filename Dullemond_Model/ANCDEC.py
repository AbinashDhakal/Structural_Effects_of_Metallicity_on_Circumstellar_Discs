import numpy as np
#import sympy as smp
#from sympy import diff
#from sympy import erf
import math

#Natural constant
GG =  6.67430e-8  # U# gravitational constant  [dyn cm^2 g^-2]
kb =  1.380649e-16  # Boltzmann constant defined [erg K^-1]
hp =  6.62607015e-27  # Planck constant [erg s] defined
clight = 2.99792458e10  # light speed [cm s^-1] defined 
NA = 6.02214076e23  # Avogadro constant defined 
mp = 1.672621898e-24  # proton mass [g]; not 1.0/NA anymore. 
ss = 5.6703e-5  # Stefan-Boltzmann const  [erg/cm^2/K^4/s]

#Astronomical constant
AU = 1.495978707e13  # Distance between Earth and Sun [cm] defined 
pc = 3.085677581e18  # Parsec [ cm] 
ms = 1.9884e33  # Sun mass [g] 
ts = 5.777e3  # effective temperature of the sun [K] 
ls = 3.828e33  # solar intensity [erg/s] defined by IAU 
rs = 6.96e10  # solar radius [cm] 

#Gas constant
mu = 2.353  # Average molecular weight  (H2 + He + metals)

# grid parameters 
nr = 256  # Number of radial grids. 
ntheta = 257  # Number of grids in elevation.
nz = 256 # Number of elevation grids. 

#Stellar Parameter
mstar = 2.4 * ms
rstar = 6.4 * rs
tstar = 9500
lstar = 47* ls

#Disk parameter
T_rim= t_sublimation = 1500.0  # Temperature at which dust start to condense
sigma0 = sigmao = 2*10 **3  # Surface density at 1 AU  Σ0 [g/cm^2] 
kp_stellar = 400.0  # Monochromatic opacity = 400cm^2 g^-1

nu = 5.879 * (10 ** 10) * tstar  #Wien law
t_virial = GG * mstar * mu * mp / (kb * rstar)  # Virial Temp
gamma = 2.0
beta = -1.5


def rim_function(variable ):
    x1, x2, x3, x4, x5 = variable
    #R_rim, h_rim, H_rim, X_rim, sigma_rim
    f =[((lstar/(np.pi* 4*T_rim**4 *ss))**0.5 *(1+ x3/x1) -x1),
           (((kb*T_rim*x1**3)/(mu*mp*GG*mstar))**0.5 -x2),
           (4*x5*kp_stellar*(1-math.erf(x4)) - 1),
           (x2*x4 - x3),
           (sigma0*(x1/AU)**beta - x5)]
    return f
pi = np.pi
def jacobian_rim(variable):
    x1, x2, x3, x4, x5 = variable
    J = [[-0.25*x3*(lstar/(T_rim**4*ss))**0.5/(pi**0.5*x1**2*(x3/x1 + 1)**0.5) - 1,       0,          0.25*(lstar/(T_rim**4*ss))**0.5/(pi**0.5*x1*(x3/x1 + 1)**0.5),       0,     0],
		[1.5*(x1**3*T_rim*kb/(GG*mp*mstar*mu))**0.5/x1 ,                         -1    , 0,0,0],
		[0 ,      x4,    -1,    x2,0],
		[0,0,0,-8*kp_stellar*x5*np.exp(-x4**2)/np.sqrt(pi),4*kp_stellar*(1 - math.erf(x4))],
        [beta*sigma0*(x1/AU)**beta/x1, 0, 0, 0, -1]]
    return J


"""
def iter_newton(X,function,jacobian,imax = 1e6,tol = 1e-5):
    for i in range(int(imax)):
        J = jacobian(X) # calculate jacobian J = df(X)/dY(X) 
        Y = function(X) # calculate function Y = f(X)
        dX = np.linalg.solve(J,Y) # solve for increment from JdX = Y 
        X -= dX # step X by dX 
        print(dX)
        if np.linalg.norm(dX)<tol: # break if converged
            print('converged.')
            break
    return X
"""

"""

def x_delta_by_gauss(J,b):

    return np.linalg.solve(J,b)

bezao = (rim_function(1,2,3,4,5))
jotinha  = (jacobian_rim(1,2,3,2,5))


print (x_delta_by_gauss(jotinha, bezao))
x_delta_test = x_delta_by_gauss(jotinha,bezao)

def x_plus_1(x_delta,x_previous):

    x_next = x_previous + x_delta

    return x_next

print (x_plus_1(x_delta_test,[1,2,3,4,5]))



def newton_method(x_init):

    first = x_init[0]

    second = x_init[1]

    third = x_init[2]
    fourth = x_init[3]
    fifth = x_init[4]

    jacobian = jacobian_rim(first, second, third, fourth, fifth)

    vector_b_f_output = rim_function(first, second, third, fourth, fifth)

    x_delta = x_delta_by_gauss(jacobian, vector_b_f_output)

    x_plus_1 = x_delta + x_init

    return x_plus_1

def iterative_newton(x_init):

    counter = 0

    x_old = list(map(abs,x_init))
   #print ("x_old", x_old)

    x_new =list(map(abs, newton_method(x_old)))
   #print ("x_new", x_new)

    diff = np.linalg.norm(x_old-x_new)
   #print (diff)

    while diff>0.0000000000001:

        counter += 1

       #print ("x_old", x_old)
        x_new = list(map(abs, newton_method(x_old)))
       #print ("x_new", x_new)

        diff = np.linalg.norm(x_old-x_new)
       #print (diff)

        x_old = x_new

    convergent_val = x_new
   #print (counter)

    return convergent_val

#print (iterative_newton([1,2]))
print("Hello")
print (list(map(float,(iterative_newton([100,200,3,4,5])))))

"""


def iterative_newton(fun, x_init, jacobian):
    max_iter = 50
    epsilon = 100000000

    x_last = x_init

    for k in range(max_iter):
        # Solve J(xn)*( xn+1 - xn ) = -F(xn):
        J = np.array(jacobian(x_last))
        print(J)
        F = np.array(fun(x_last))
        print(F)
        diff = np.linalg.solve( J, -F )
        x_last = x_last + diff

        # Stop condition:
        if np.linalg.norm(diff) < epsilon:
            print('convergence!, nre iter:', k )
            break

    else: # only if the for loop end 'naturally'
        print('not converged')

    return x_last


X_guess_rim =abs(np.array( [0.4904092456436323*AU, 0.39336608131611195*AU,0.13113940101171412*AU,4.512601124660966*AU,1650.5889740859373], dtype =float))
x_sol = iterative_newton(rim_function,  [0.4904092456436323*AU, 0.39336608131611195*AU,0.13113940101171412*AU,4.512601124660966*AU,1650.5889740859373], jacobian_rim)
print('solution exercice:', x_sol )
print('F(sol)', rim_function(x_sol) )



#x_sol = iterative_newton(rim_function, X_guess_rim, jacobian_rim)
#print('solution exercice:', x_sol )
#print('F(sol)', rim_function(x_sol) )