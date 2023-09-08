# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 21:14:20 2022

@author: abudh
"""

import numpy as np
#import sympy as smp
#from sympy import diff
#from sympy import erf
import math
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt



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
psi_s, psi_i = 1,1

def rim_function(variable):
    x1, x2, x3, x4, x5 = variable
    #R_rim, h_rim, H_rim, X_rim, sigma_rim
    f =[((lstar/(np.pi* 4*T_rim**4 *ss))**0.5 *(1+ x3/x1)**0.5 -x1),
           (((kb*T_rim*x1**3)/(mu*mp*GG*mstar))**0.5 -x2),
           (x2*x4 - x3),
           (4*x5*kp_stellar*(1-math.erf(x4)) - 1),
           (sigma0*(x1/AU)**beta - x5)]
    return f


pi = np.pi
def jacobian_rim(variable):
    x1, x2, x3, x4, x5 = variable
    #R_rim, h_rim, H_rim, X_rim, sigma_rim

    J = [[-0.25*x3*(lstar/(T_rim**4*ss))**0.5/(pi**0.5*x1**2*(x3/x1 + 1)**0.5) - 1,       0,          0.25*(lstar/(T_rim**4*ss))**0.5/(pi**0.5*x1*(x3/x1 + 1)**0.5),       0,     0],
		[1.5*(x1**3*T_rim*kb/(GG*mp*mstar*mu))**0.5/x1 ,                         -1    , 0,0,0],
		[0 ,      x4,    -1,    x2,0],
		[0,0,0,-8*kp_stellar*x5*np.exp(-x4**2)/np.sqrt(pi),4*kp_stellar*(1 - math.erf(x4))],
        [beta*sigma0*(x1/AU)**beta/x1, 0, 0, 0, -1]]
    return J



#Definig grid system in radial direction
rin = 0.52 *AU # Sublimation radius  # inner rim of disk. Density is assumed to
rout = 100*AU
ri = np.linspace(rin,rout,nr+1) 
R   = 0.5 * ( ri[0:nr] + ri[1:nr+1] )             # Take the in-"between values
R = np.array(R)

#sigma[g/cm^2] 
def Sigma(R):
    return sigma0 * (R/AU)**beta

sigma = Sigma(R)
sigma = np.array(sigma)


def disk_function(variable,R, psi_s, psi_i):
    x1, x2, x3, x4, x5 = variable
    #X_cg, H_cg, hcg, Ti, alpha
    f =[(1-math.erf(x1/np.sqrt(2)) -2* x5/(sigma* kp_stellar)),
       (x1*x3 -x2),
       ((x4/t_virial)**0.5 * (R/rstar)**0.5 -x3/R),
       (x5*psi_s/psi_i/2)**0.25 * (rstar/R)**0.5 *tstar -x4,
       (0.4*rstar/R) +(gamma-1)*x2/R -x5]
    return f
    

def jacobian_disk(variable, R, psi_s, psi_i):
    x1, x2, x3, x4, x5 = variable
    #X_cg, H_cg, hcg, Ti, alpha

    J = [[-np.sqrt(2)*np.exp(-x1**2/2)/np.sqrt(pi) , 0 ,0, 0 ,-2/(kp_stellar*sigma)],
		[x3, -1 ,x1 ,0, 0],
		[0, 0, -1/R ,0.5*(R/rstar)**0.5*(x4/t_virial)**0.5/x4, 0],
		[0,   0,    0, -1, tstar*(rstar/R)**0.5*(psi_s*x5/psi_i/2)**0.25/x5],
        [0, (gamma - 1)/R ,0, 0, -1]]
    return J



def iterative_newton(fun, x_init, jacobian):
    max_iter = 50
    epsilon = 1e-8

    x_last = x_init
    listt = []
    for k in range(max_iter):
        # Solve J(xn)*( xn+1 - xn ) = -F(xn):
        J = np.array(jacobian(x_last))
        F = np.array(fun(x_last))

        diff = np.linalg.solve( J, -F )
        x_last = x_last + diff
        #print("I am trying to converge:  ",np.linalg.norm(diff) ,k, "iteration" )
        listt.append(round(np.linalg.norm(diff), 1))

        # Stop condition:
        if np.linalg.norm(diff) < epsilon:
            print('convergence!, nre iter:', k )
            break

    else: # only if the for loop end 'naturally'
        print('not converged')
    
    print(listt)
    return x_last



sample_data = pd.read_csv('Ti_dataset.csv')
sample_data2 = pd.read_csv('hR_datasets.csv')
sample_data3= pd.read_csv('HHHR_datasets.csv')
fig, ax = plt.subplots(figsize = (9, 6))


R_iterpolate_range = np.linspace(rin/AU,rout/AU,nr)
Temp_interpolate_f =interp1d(sample_data.R, sample_data.Ti, 'quadratic')
hR_interpolate_f =interp1d(sample_data2.R, sample_data2.h, 'quadratic')
HR_interpolate_f =interp1d(sample_data3.R, sample_data3.H, 'quadratic')

Ti_interpolate = Temp_interpolate_f(R_iterpolate_range)
hR_interpolate = hR_interpolate_f(R_iterpolate_range)
h_interpolate = hR_interpolate_f(R_iterpolate_range)*AU
HR_interpolate = np.array(HR_interpolate_f(R_iterpolate_range))
H_interpolate = np.array(HR_interpolate_f(R_iterpolate_range))*AU

X_interpolate = HR_interpolate/hR_interpolate
alpha_interpolate = (0.4*rstar/R_iterpolate_range/AU) + (gamma -1)*HR_interpolate
plt.semilogx(R_iterpolate_range, Ti_interpolate)
plt.show


x_init_rim =[0.539*AU, 0.019*AU,0.09344*AU,1.53,4500]
x_sol_rim = iterative_newton(rim_function, x_init_rim, jacobian_rim)
print('solution exercice:', x_sol_rim )
print('F(sol)', rim_function(x_sol_rim) )


X_cg = np.zeros(nr)
impinge_angle = np.zeros(nr)
H_cg = np.zeros(nr)
h_cg = np.zeros(nr)
Ti = np.zeros(nr)

"""
for i, (r_i, s_i) in enumerate(zip(R, sigma)):
    x_init_disk = [X_interpolate[i] , H_interpolate[i], h_interpolate[i],Ti_interpolate[i],alpha_interpolate[i]]
    #X_cg, H_cg, hcg, Ti, alpha
    x_sol_disk = iterative_newton(disk_function, x_init_disk, jacobian_disk, args = (r_i, s_i))
    X_cg[i], H_cg[i], h_cg[i],Ti[i], impinge_angle[i] = x_sol_disk[0], x_sol_disk[1], x_sol_disk[2],x_sol_disk[3],x_sol_disk[4]
    #print(x_init_disk)
    #print(x_sol_disk)
"""


