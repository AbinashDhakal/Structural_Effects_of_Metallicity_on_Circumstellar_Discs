import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import fsolve
import pandas as pd
import math
from scipy.special import erfinv

import sympy as smp

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
nr = 1000  # Number of radial grids. 
ntheta = 257  # Number of grids in elevation.
nz = 1000 # Number of elevation grids. 

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

#R_rim, h_rim, H_rim, X_rim,sigma_rim = smp.symbols('R_rim h_rim H_rim X_rim sigma_rim')
R_rim, h_rim, H_rim, X_rim,sigma_rim = smp.symbols('R_rim, h_rim, H_rim, X_rim,sigma_rim')

#Define R_rim
eqn1 = smp.Eq( (lstar/(4*smp.pi * T_rim**4 *ss))**0.5 * (1 + H_rim/R_rim)**0.5,R_rim)
#Define h_rim
eqn2 =smp.Eq((kb*T_rim*R_rim**3/(mu*mp*GG*mstar))**0.5, h_rim)
#Define H_rim
eqn3 = smp.Eq(X_rim * h_rim,H_rim)
#Define X_rim
#integral[H/h = X  --- inf] of e^-x^2 dx =(Pi^0.5 /2)(1- erf(X_rim))
# Tau(zo) = 8*sigma* kp_stellar* integral/(pi^0.5) 
eqn4 = smp.Eq((1-smp.erf(X_rim))*(4*sigma_rim*kp_stellar), 1)
#Define surface density
eqn5 = smp.Eq(sigma0*(R_rim/AU)**beta, sigma_rim)

print(smp.solve([eqn1, eqn2,eqn3,eqn4, eqn5], (R_rim, h_rim, H_rim, X_rim,sigma_rim)))
"""
def Rim(vars):
    R_rim, h_rim, H_rim, X_rim,sigma_rim = vars
    #Define R_rim
    eqn1 = (lstar/(4*np.pi * T_rim**4 *ss))**0.5 * (1 + H_rim/R_rim)**0.5 - R_rim
    #Define h_rim
    eqn2 =(kb*T_rim*R_rim**3/(mu*mp*GG*mstar))**0.5 - h_rim
    #Define H_rim
    eqn3 = X_rim * h_rim -H_rim
    #Define X_rim
    #integral[H/h = X  --- inf] of e^-x^2 dx =(Pi^0.5 /2)(1- erf(X_rim))
    # Tau(zo) = 8*sigma* kp_stellar* integral/(pi^0.5) 
    eqn4 = (1-math.erf(X_rim))*(4*sigma_rim*kp_stellar) -1
    #Define surface density
    eqn5 = sigma0*(R_rim/AU)**beta -sigma_rim
    return [eqn1, eqn2,eqn3,eqn4, eqn5]

R_rim, h_rim, H_rim,X_rim, sigma_rim = fsolve(Rim, (0.52*AU, 0.1*AU,0.057*AU,5.3,6207))
#print(R_rim/AU,H_rim/R_rim, X_rim)
"""