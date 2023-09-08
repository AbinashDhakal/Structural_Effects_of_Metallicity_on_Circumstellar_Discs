import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import fsolve
import pandas as pd
import math
from scipy.special import erfinv
import scipy.optimize
import itertools
import random 

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


def Rim(vars):
    R_rim, h_rim, H_rim, X_rim,sigma_rim = vars
    #Define R_rim
    eqn1 = ((lstar/(4*np.pi * T_rim**4 *ss))**0.5 * (1 + H_rim/R_rim)**0.5 - R_rim)
    #Define h_rim
    eqn2 =((kb*T_rim*R_rim**3/(mu*mp*GG*mstar))**0.5 - h_rim)
    #Define H_rim
    eqn3 = (X_rim * h_rim -H_rim)
    #Define X_rim
    #integral[H/h = X  --- inf] of e^-x^2 dx =(Pi^0.5 /2)(1- erf(X_rim))
    # Tau(zo) = 8*sigma* kp_stellar* integral/(pi^0.5) 
    eqn4 = erfinv(1-1/(4*sigma_rim*kp_stellar)) -X_rim
    #eqn4 = ((1-math.erf(1/X_rim))*(4*sigma_rim*kp_stellar) -1)
    #Define surface density
    eqn5 = (sigma0*(R_rim/AU)**beta -sigma_rim)
    return [eqn1, eqn2,eqn3,eqn4, eqn5]

X_guess_rim = np.array([0.4904092456436323*AU, 0.39336608131611195*AU,0.13113940101171412*AU,4.512601124660966*AU,1650.5889740859373])
#X_guess_rim = np.array([random.uniform(0.01*AU, AU), random.uniform(0.01*AU, AU),random.uniform(0.001*AU, 0.2*AU),random.uniform(0.1,5),random.uniform(1000, 2000)])

R_rim, h_rim, H_rim,X_rim, sigma_rim = fsolve(Rim,X_guess_rim)

"""
variable_rim = R_rim, h_rim, H_rim,X_rim, sigma_rim = fsolve(Rim, X_guess_rim)

def Optimise(variable, X_guess, function):
    N = 10
    best_obj = np.inf
    best_variable = 0
    for i in range(N):
        X_guess = np.random.uniform(size =X_guess.shape) *variable
        variable =fsolve(Rim, X_guess)
        #print(variable)
        fun = np.linalg.norm(variable)
        #positive_Check = (all(number > 0 for number in variable))
        #print(positive_Check) # Check whether all the variable is positive or not
        if (fun < best_obj):
            best_obj =fun
            best_variable =variable
    #print(variable)
    return best_variable

            
variable_rim = Optimise(variable_rim, X_guess_rim, Rim)
R_rim, h_rim, H_rim,X_rim, sigma_rim = variable_rim[0], variable_rim[1], variable_rim[2], variable_rim[3], variable_rim[4]
#print(variable_rim)

"""
#Definig grid system in radial direction
rin = R_rim  # Sublimation radius  # inner rim of disk. Density is assumed to be 0 inside [cm] 
rout = 100*AU
#rout = 2.7 * (10 ** 2) * AU  # Outer edge of disk. Density is set to 0 outside this [cm] 
ri = np.linspace(rin,rout,nr+1) 
#R = np.array([0.5 * (ri[i:] + ri[:i]) for i in range(len(ri))])  # Take the in-"between values
#R   = np. array(0.5 * ( ri[::-1] + ri) ) 
R   = 0.5 * ( ri[0:nr] + ri[1:nr+1] )             # Take the in-"between values
R = np.array(R)

#sigma[g/cm^2] 
def Sigma(R):
    return sigma0 * (R/AU)**beta

sigma = Sigma(R)
sigma = np.array(sigma)

#def Disk(X_cg, impinge_angle, H_cg, h_cg, Ti, R,  sigma, psi_i , psi_s):
def Disk(vars,  R,  sigma, psi_i , psi_s):
    X_cg, impinge_angle, H_cg, h_cg, Ti =vars
   #Define X_cg
    #eqn1 = 1-math.erf(X_cg/np.sqrt(2)) -2* impinge_angle/(sigma* kp_stellar)
    eqn1 = (2**0.5)*erfinv(1- 2*impinge_angle/kp_stellar/sigma) -X_cg
    #Define H_cg
    eqn2 = X_cg*h_cg -H_cg
    #Define h_cg
    eqn3 = (Ti/t_virial)**0.5 * (R/rstar)**0.5 -h_cg/R
    #Define Ti
    eqn4 = (impinge_angle*psi_s/psi_i)**0.25 * (rstar/R)**0.5 *tstar -Ti
    #Define Impinge angle
    eqn5 = (0.4*rstar/R) +(gamma-1)*H_cg/R -impinge_angle
    return [eqn1, eqn2,eqn3,eqn4, eqn5]

X_cg = np.zeros(nr)
impinge_angle = np.zeros(nr)
H_cg = np.zeros(nr)
h_cg = np.zeros(nr)
Ti= np.zeros(nr)
psi_i =np.array(list(itertools.repeat(1, nr)))
psi_s = np.array(list(itertools.repeat(1, nr)))
X_guess_disk = np.array([X_rim, 0.15, H_rim, h_rim, T_rim])

for i,( r_i,s_i, p_i,p_s), in enumerate(zip(R,sigma, psi_i, psi_s)):
    X_cg[i], H_cg[i], h_cg[i], Ti[i],impinge_angle[i] =fsolve(Disk,X_guess_disk, args = (r_i, s_i, p_i, p_s))

fig, ax = plt.subplots()
#ax.plot(R, Ti)
plt.plot(R/AU, Ti)

ax.set_xlabel("Radius(AU)")
ax.set_ylabel("Temperature")
ax.set_title("Protoplanetary System Temp vs. Radius")

fig, ax = plt.subplots()
ax.plot(R/AU, h_cg)
ax.set_xlabel("Radius")
ax.set_ylabel("Surface_height")
ax.set_title("Surface_height vs. Radius")
plt.show()


#Solving for Rim and Disk for 1000 times 
Iter =10000
for i in range(Iter):
    ri = np.linspace(R_rim,rout,nr+1) 
    R   = 0.5 * ( ri[0:nr] + ri[1:nr+1] )             # Take the in-"between values
    R = np.array(R)
    #print(R_rim)
    R_rim, h_rim, H_rim,X_rim, sigma_rim = fsolve(Rim,X_guess_rim)
    X_guess_rim =  np.array([R_rim, h_rim, H_rim,X_rim, sigma_rim])*random.uniform(0.95,1)
    #updat new value with each element on ith column and  
    #Standardised the list and output the std as best value
    sigma = Sigma(R)
    for i,( r_i,s_i, p_i,p_s), in enumerate(zip(R,sigma, psi_i, psi_s)):
        X_cg[i], H_cg[i], h_cg[i], Ti[i],impinge_angle[i] =fsolve(Disk,X_guess_disk, args = (r_i, s_i, p_i, p_s))
        X_guess_disk =np.array([X_cg[i], H_cg[i], h_cg[i], Ti[i],impinge_angle[i]])
        print(X_guess_disk)
ax.set_xlabel("Radius(AU)")
ax.set_ylabel("Temperature")
ax.set_title("Protoplanetary System Temp vs. Radius")

fig, ax = plt.subplots()
ax.plot(R/AU, h_cg)
ax.set_xlabel("Radius")
ax.set_ylabel("Surface_height")
ax.set_title("Surface_height vs. Radius")
plt.show()