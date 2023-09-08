import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import fsolve
import pandas as pd
import math
from scipy.special import erfinv

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


#Definig grid system in radial direction
rin = R_rim =0.4904092456436323* AU # Sublimation radius  # inner rim of disk. Density is assumed to be 0 inside [cm] 
rout = 2*AU
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


df = pd.read_csv("Planck_data.csv")
filtering_data =( df[
    (df["T_stellar"] == 3000) 
    & (df["Metallicity"] == -0.3 )
   ])
T_data = np.array(list(filtering_data["T_gas"]))
kp_data = np.array(list(filtering_data["Planck_opacity"]))
kp_t_array = np.zeros((2, len(kp_data)))
kp_t_array[0] =T_data
kp_t_array[1] =kp_data
#change_array = T_kp_data.to_numpy()

impinge_angle = pd.read_csv("alpha.csv")
alpha =np.array(list(impinge_angle["aalpha"]))

def kp_interp(T_guess):
    """Interpolate Temperature to find accurate Kp(T) value."""
    kp_t_guess = np.interp(T_guess, T_data, kp_data)
    #kp_t_guess = np.interp(T_guess, T_data, Kp_data)

    return kp_t_guess

def t_s_eq(T_iter, const):
    """Solve iteration for Ts."""
    #print(kp_interp(T_iter))
    eps_iter = (kp_interp(tstar) / kp_interp(T_iter)) ** 0.25
    #print(kp_interp(T_iter))
    return eps_iter * const - T_iter


def t_i_eq(T_iter, const):
    """Solve iteration for Ti."""
    psi_i_iter = kp_interp(T_iter)
    return const * (psi_i_iter ** -0.25) - T_iter
"""
def Temp_disk(vars,R, sigma):
    Psi_s, Psi_i,Ti,Ts,epsilon =vars
    eqn1 = (alpha*Psi_s/Psi_i)**0.25 *(rstar/R)**0.5 *tstar -Ti
    eqn2 = (epsilon)**0.25 * (rstar/2/ R)**0.5 *tstar -Ts
    eqn3 =  (kp_interp(tstar) / kp_interp(Ts)) -epsilon
    eqn4 = (sigma*kp_interp(Ti))-Psi_i
    eqn5 = (sigma*kp_interp(Ts)) - Psi_s
    return [eqn1,eqn2, eqn3, eqn4, eqn5]

Psi_s= np.zeros(nr)
Psi_i= np.zeros(nr)
Ti= np.zeros(nr)
Ts= np.zeros(nr)
epsilon = np.zeros(nr)

X_guess = np.array([3,3,t_sublimation, t_sublimation,3])
for i, (r_i,s_i) in enumerate(zip(R,sigma)):
    Psi_s[i], Psi_i[i] ,Ti[i] ,Ts[i],epsilon[i] = fsolve(Temp_disk,X_guess, args = (r_i,s_i))
    print(Psi_s[i])





#R_rim, h_rim, H_rim,X_rim, sigma_rim = fsolve(Rim, (0.52*AU, 0.1*AU,0.057*AU,5.3,6207))
X_guess_rim = np.array([3,3,t_sublimation, t_sublimation,3])
#X_guess_rim = np.array([random.uniform(0.01*AU, AU), random.uniform(0.01*AU, AU),random.uniform(0.001*AU, 0.2*AU),random.uniform(0.1,5),random.uniform(1000, 2000)])

#variable_Temp = Psi_s, Psi_i, Ti, Ts, epsilion= fsolve(Rim, X_guess_rim)
"""
"""
def Optimise(variable, X_guess, function):
    N = 1000
    best_obj = np.inf
    best_variable = 0
    for i in range(N):
        X_guess = np.array([random.uniform(0.01*AU, AU), random.uniform(0.01*AU, AU),random.uniform(0.001*AU, 0.2*AU),random.uniform(0.1,5),random.uniform(1000, 2000)])
        variable =  fsolve(Rim, X_guess)
        #print(variable)
        #X_guess = np.random.uniform(size =X_guess.shape)
        #print(X_guess)
        fun = np.linalg.norm(variable)
        if fun < best_obj:
            best_obj =fun
            best_variable =variable
    return variable

            
variable_rim = Optimise(variable_rim, X_guess_rim, Rim)
#R_rim, h_rim, H_rim,X_rim, sigma_rim = variable_rim[0], variable_rim[1], variable_rim[2], variable_rim[3], variable_rim[4]
#print(variable_rim)

fig, ax = plt.subplots()
ax.plot(R/AU, Ti)
ax.set_xlabel("Radius")
ax.set_ylabel("Surface_Height")
ax.set_title("Surface_Height vs. Radius")
plt.show()
"""

"""
"""
#alpha = ([])

kp_tstar = kp_interp(tstar)
T_s = np.zeros(len(R))
T_s_0 = 3000.0
Psi_s = np.zeros(len(R))
#alpha = impinge_angle
T_i = np.zeros(len(R))
T_i_0 = 3000.0
Psi_i = np.zeros(len(R))
T_prev = t_sublimation
for i, (r_i,s_i) in enumerate(zip(R,sigma)):
   # print(i)
    # Numerically solve for T_s
    const_1 = np.sqrt(0.5 * rstar / r_i) * tstar
   # print(const_1)
    T_s[i] = fsolve(t_s_eq, t_sublimation, args=const_1)
    #print(T_s[i])
    # Numerical solution for T_s generates Psi_s
    Psi_s[i] = kp_interp(T_s[i]) 
    #print(Psi_s[i])
   # print(Psi_s[i])
    # Simplify constant elements of equation for T_i
    #const_2 = ( np.sqrt(rstar / r_i)) * tstar
    #const_2 = ((alpha[i] * Psi_s[i]/Psi_i[i]) ** 0.25) * np.sqrt(rstar / r_i) * tstar
    # Numerically solve for T_i via Psi_i
    #T_i[i] = fsolve(t_i_eq, T_i_0, args=const_2)
    # Numerical solution for T_i generates Psi_i
   # Psi_i[i] = kp_interp(T_i[i])*s_i
   # Psi_frac = Psi_s[i]/Psi_i[i]
    #print(Psi_frac)

for i, (r_i,s_i,Ps_s,a_i) in enumerate(zip(R,sigma, Psi_s,alpha)):

    #print(Psi_s[i])
   # print(Psi_s[i])
    # Simplify constant elements of equation for T_i
    #const_2 = ( np.sqrt(rstar / r_i)) * tstar
    #print(a_i)
    const_2 = ((alpha[i] * Psi_s[i]/Psi_i[i]) ** 0.25) * np.sqrt(rstar / r_i) * tstar
    #print(const_2)
    # Numerically solve for T_i via Psi_i
    T_i[i] = fsolve(t_i_eq, T_i_0, args=const_2)
    print(const_2)
    # Numerical solution for T_i generates Psi_i
   # Psi_i[i] = kp_interp(T_i[i])*s_i
   # Psi_frac = Psi_s[i]/Psi_i[i]
    #print(Psi_frac)

fig, ax = plt.subplots()
ax.plot(R/AU, T_s)
ax.set_xlabel("Radius")
ax.set_ylabel("Surface_Height")
ax.set_title("Surface_Height vs. Radius")
plt.show()
