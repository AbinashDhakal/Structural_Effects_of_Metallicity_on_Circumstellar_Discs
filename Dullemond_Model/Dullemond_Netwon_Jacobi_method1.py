import numpy as np
import sympy as smp
from sympy import diff
from sympy import erf



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
sigma0 = 2*10 **3  # Surface density at 1 AU  Σ0 [g/cm^2] 
kp_stellar = 400.0  # Monochromatic opacity = 400cm^2 g^-1

nu = 5.879 * (10 ** 10) * tstar  #Wien law
t_virial = GG * mstar * mu * mp / (kb * rstar)  # Virial Temp
gamma = 2.0
beta = -1.5


pi = np.pi
"""
def eqn_Rim(x):
    eqn1 = (lstar/(4*np.pi * T_rim**4 *ss))**0.5 * (1 + x[2]/x[0])**0.5 - x[0]
    #Define h_rim
    eqn2 =(kb*T_rim*x[0]**3/(mu*mp*GG*mstar))**0.5 - x[1]
    #Define H_rim
    eqn3 = x[3] *x[1] -x[2]
    #Define X_rim
    #integral[H/h = X  --- inf] of e^-x^2 dx =(Pi^0.5 /2)(1- erf(1/X_rim))
    # Tau(zo) = 8*sigma* kp_stellar* integral/(pi^0.5) 
    eqn4 = (1-smp.erf(x[0]))*(4*sigma_rim*kp_stellar) -1
    #Define surface density
    eqn5 = sigma0*(x[0]/AU)**beta -x[4]
    return eqn1 ,eqn2, eqn3, eqn4, eqn5
"""

#import symengine
#lstar, rstar, T_rim,kb, mu,mp,GG,mstar, kp_stellar,beta,AU, sigma0, ss = smp.symbols('lstar, rstar, T_rim,kb, mu,mp,GG,mstar, kp_stellar,beta,AU, sigma0 ss')

"""
eqn1 = (lstar/(4*np.pi * T_rim**4 *ss))**0.5 * (1 + H_rim/R_rim)**0.5 - R_rim
#Define h_rim
eqn2 =(kb*T_rim*R_rim**3/(mu*mp*GG*mstar))**0.5 - h_rim
#Define H_rim
eqn3 = X_rim * h_rim -H_rim
#Define X_rim
#integral[H/h = X  --- inf] of e^-x^2 dx =(Pi^0.5 /2)(1- erf(1/X_rim))
# Tau(zo) = 8*sigma* kp_stellar* integral/(pi^0.5) 
eqn4 = (1-smp.erf(X_rim))*(4*sigma_rim*kp_stellar) -1
#Define surface density
eqn5 = sigma0*(R_rim/AU)**beta -sigma_rim
"""
"""
vars = smp.symbols('R_rim h_rim H_rim X_rim sigma_rim') # Define x and y variables


f = smp.sympify(['(lstar/(4*pi * T_rim**4 *ss))**0.5 * (1 + H_rim/R_rim)**0.5 - R_rim',
                 '(kb*T_rim*R_rim**3/(mu*mp*GG*mstar))**0.5 - h_rim',
                 'X_rim * h_rim -H_rim',
                 '((1-erf(X_rim))*(4*sigma_rim*kp_stellar) -1)',
                 'sigma0*(R_rim/AU)**beta -sigma_rim']) # Define function
J = smp.zeros(len(f),len(vars)) # Initialise Jacobian matrix
"""


"""
vars = smp.symbols('R_rim h_rim H_rim X_rim sigma_rim', positive=True, real=True) # Define x and y variables


f = (['pow((lstar/(4*pi*(pow(T_rim,4))*ss)), 0.5)*pow((1+ H_rim/R_rim),0.5) -R_rim',
                 'pow((kb*T_rim*pow(R_rim,3)/(mu*mp*GG*mstar)),0.5) - h_rim',
                 'X_rim * h_rim -H_rim',
                 '((1-erf(X_rim))*(4*sigma_rim*kp_stellar) -1)',
                 'sigma0*pow((R_rim/AU),beta) -sigma_rim']) # Define function
print(f)
J = smp.zeros(len(f),len(vars)) # Initialise Jacobian matrix

# Fill Jacobian matrix with entries
for i, fi in enumerate(f):
    for j, s in enumerate(vars):
        J[i,j] = smp.diff(fi, s)
        print(J[i,j])

print (J)

"""
vars = smp.symbols('R_rim h_rim H_rim X_rim sigma_rim') # Define x and y variables

def function(vars):
    f = (['pow((lstar/(4*pi*(pow(T_rim,4))*ss)), 0.5)*pow((1+ H_rim/R_rim),0.5) -R_rim',
                     'pow((kb*T_rim*pow(R_rim,3)/(mu*mp*GG*mstar)),0.5) - h_rim',
                     'X_rim * h_rim -H_rim',
                     '((1-erf(X_rim))*(4*sigma_rim*kp_stellar) -1)',
                     'sigma0*pow((R_rim/AU),beta) -sigma_rim']) # Define function
    return f


def Jacobian(vars):
    f = function(vars)
    J = smp.zeros(len(f),len(vars)) # Initialise Jacobian matrix

    # Fill Jacobian matrix with entries
    for i, fi in enumerate(f):
        for j, s in enumerate(vars):
            J[i,j] = smp.diff(fi, s)
            print(J[i,j])
    return J
    
def iterative_newton(fun, x_init, jacobian):
    max_iter = 50
    epsilon = 100

    x_last = x_init

    for k in range(max_iter):
        # Solve J(xn)*( xn+1 - xn ) = -F(xn):
        J = np.array(jacobian(x_last))
        F = np.array(fun(x_last))

        diff = np.linalg.solve( J, -F )
        x_last = x_last + diff

        # Stop condition:
        if np.linalg.norm(diff) < epsilon:
            print('convergence!, nre iter:', k )
            break

    else: # only if the for loop end 'naturally'
        print('not converged')

    return x_last
"""
def iterative_newton(fun, x_init, jacobian):
    max_iter = 50
    epsilon = 100

    x_last = x_init

    for k in range(max_iter):
        # Solve J(xn)*( xn+1 - xn ) = -F(xn):
        J = np.array(jacobian(x_last))
        F = np.array(fun(x_last))

        diff = np.linalg.solve( J, -F )
        x_last = x_last + diff

        # Stop condition:
        if np.linalg.norm(diff) < epsilon:
            print('convergence!, nre iter:', k )
            break

    else: # only if the for loop end 'naturally'
        print('not converged')

    return x_last
"""
x_sol = iterative_newton(function, [0.4904092456436323*AU, 0.39336608131611195*AU,0.13113940101171412*AU,4.512601124660966*AU,1650.5889740859373], Jacobian)
print('solution exercice:', x_sol )
