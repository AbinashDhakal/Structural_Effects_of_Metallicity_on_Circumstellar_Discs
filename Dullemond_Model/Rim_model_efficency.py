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


#R_rim, h_rim, H_rim,X_rim, sigma_rim = fsolve(Rim, (0.52*AU, 0.1*AU,0.057*AU,5.3,6207))
X_guess_rim = np.array([0.4904092456436323*AU, 0.39336608131611195*AU,0.13113940101171412*AU,4.512601124660966,1650.5889740859373])
#X_guess_rim = np.array([random.uniform(0.01*AU, AU), random.uniform(0.01*AU, AU),random.uniform(0.001*AU, 0.2*AU),random.uniform(0.1,5),random.uniform(1000, 2000)])

variable_rim = R_rim, h_rim, H_rim,X_rim, sigma_rim = fsolve(Rim, X_guess_rim)
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
        fun = np.linalg.norm(variabl2e)
        if fun < best_obj:
            best_obj =fun
            best_variable =variable
    return variable

            
variable_rim = Optimise(variable_rim, X_guess_rim, Rim)
#R_rim, h_rim, H_rim,X_rim, sigma_rim = variable_rim[0], variable_rim[1], variable_rim[2], variable_rim[3], variable_rim[4]
#print(variable_rim)
"""

#Definig grid system in radial direction
rin = R_rim  # Sublimation radius  # inner rim of disk. Density is assumed to be 0 inside [cm] 
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


#def Disk(X_cg, impinge_angle, H_cg, h_cg, Ti, R,  sigma, psi_i , psi_s):
def Disk(vars,R,  sigma, psi_i , psi_s):
    X_cg, H_cg, h_cg, Ti,impinge_angle   = vars
    #Define X_cg (Appendix A2)
    #Integral of exp(-z^2 /2H_cvg^2)dz = (pi/2)^0.5 x H_ccg (1- erf(1/2^2))
    #1 -erf(1/2^0.5) = 1- erf(X_cg/ 2^0.5) = 2 impingement_angle(X_cg)/ (sigma x Kp_stellar)
    eqn1 = 1-math.erf(X_cg/np.sqrt(2)) -2* impinge_angle/(sigma* kp_stellar)
    #Define H_cg
    eqn2 = X_cg*h_cg -H_cg
    #Define h_cg
    eqn3 = (Ti/t_virial)**0.5 * (R/rstar)**0.5 -h_cg/R
    #Define Ti
    eqn4 = (impinge_angle*psi_s/psi_i)**0.25 * (rstar/R)**0.5 *tstar -Ti
    #Define Impinge angle
    eqn5 = (0.4*rstar/R) +(gamma-1)*H_cg/R -impinge_angle
    return [eqn1, eqn2,eqn3,eqn4, eqn5]

psi_i, psi_s = [1,1]
X_cg = np.zeros(nr)
impinge_angle = np.zeros(nr)
H_cg = np.zeros(nr)
h_cg = np.zeros(nr)
Ti= np.zeros(nr)
for i,( r_i,s_i), in enumerate(zip(R,sigma)):
    X_cg[i], H_cg[i], h_cg[i], Ti[i],impinge_angle[i] =fsolve(Disk,(X_rim, 0.15, H_rim, h_rim, T_rim), args = (r_i, s_i, psi_i, psi_s))

print("Hello")

"""
psi_i, psi_s = [1,1]
X_guess2 = np.array([X_rim, H_rim, h_rim, T_rim,0.15])
variable_disk = X_cg, H_cg, h_cg, Ti,impinge_angle =fsolve(Disk,X_guess2, args = (R[0],sigma[0], psi_i, psi_s))

X_cg = np.zeros(nr)
impinge_angle = np.zeros(nr)
H_cg = np.zeros(nr)
h_cg = np.zeros(nr)
Ti= np.zeros(nr)
psi_i, psi_s = [1,1]
X_guess2 = np.array([X_rim, H_rim, h_rim, T_rim,0.15])
variable_disk = X_cg, H_cg, h_cg, Ti,impinge_angle =fsolve(Disk,X_guess2, args = (R[0],sigma[0], psi_i, psi_s))

def Optimise2(variable, X_guess2, function, R, sigma):
    X_cg = np.zeros(nr)
    impinge_angle = np.zeros(nr)
    H_cg = np.zeros(nr)
    h_cg = np.zeros(nr)
    Ti= np.zeros(nr)
    #print(X_cg)
    N = 10
    best_obj = np.inf
    best_variable = 0
    for ii in range(N):
        X_guess2 = np.array([random.uniform(0.1,7),random.uniform(0.01*AU, 5*AU),random.uniform(0.01*AU, AU),random.uniform(100,2000),random.uniform(0.001, 100)])
        #print(variable)
        #X_guess = np.random.uniform(size =X_guess.shape)
        #print(X_guess)
        for i,( r_i,s_i), in enumerate(zip(R,sigma)):
            X_cg[i], H_cg[i], h_cg[i], Ti[i],impinge_angle[i] =fsolve(Disk,(X_rim, 0.15, H_rim, h_rim, T_rim), args = (r_i, s_i, psi_i, psi_s))
            print(X_cg, i)
            fun = np.linalg.norm(variable)
            if fun < best_obj:
                best_obj =fun
                best_variable =variable
           # print(variable[0],i, "Abinash")
        #print(variable[0],ii, "Dhakal")
    return variable

variable_disk = Optimise2(variable_disk, X_guess2, Disk, R, sigma)

X_cg, H_cg, h_cg, Ti,impinge_angle = variable_disk[0], variable_disk[1], variable_disk[2], variable_disk[3], variable_disk[4]

#for i,( r_i,s_i), in enumerate(zip(R,sigma)):
#    X_cg[i], H_cg[i], h_cg[i], Ti[i],impinge_angle[i] =fsolve(Disk,(X_rim, 0.15, H_rim, h_rim, T_rim), args = (r_i, s_i, psi_i, psi_s))
"""

"""
def Z_grid(H_cg):
    zi = np.linspace(0.0, H_cg, nz +1)
    z_elevation   = 0.5 * (zi[0:nz] + zi[1:nz+1])          # Take the in-between values
    z_elevation = np.array(z_elevation)
    #z_elevation   = np.array([0.5 * ( zi[i] + zi[i+1])   for i in range(len(zi) - 1)])  # Take the in-"between values       # Take the in-between values
    return z_elevation
z_elevation = Z_grid(H_cg)
"""
zi = np.linspace(0.0, H_cg, nz +1)
z_elevation   = 0.5 * (zi[0:nz] + zi[1:nz+1])          # Take the in-between values
z_elevation = np.array(z_elevation)
    #z_elevation   = np.array([0.5 * ( zi[i] + zi[i+1])   for i in range(len(zi) - 1)])  # Take the in-"between values       # Take the in-between values
#rhod is spatial density [g/cm^3], sigmad is [g/cm^2]. Finally, when entering rhod, a vertical structure was assumed. 
# This time, we assumed a Gaussian distribution with scale height hh in the height direction.




def Rho(z_elevation, h_cg, sigma):
    return sigma*np.exp(- z_elevation**2 / 2/h_cg**2) / (np.sqrt(2*np.pi) *h_cg)
rho = Rho(z_elevation, h_cg, sigma)




Simulation =100
for i in range (Simulation):
    ri = np.linspace(rin,rout,nr+1) 
    R   = 0.5 * ( ri[0:nr] + ri[1:nr+1] )             # Take the in-"between values
    R = np.array(R)
    X_guess_rim = np.array([R_rim, h_rim, H_rim,X_rim, sigma_rim])
    R_rim, h_rim, H_rim,X_rim, sigma_rim = fsolve(Rim, X_guess_rim)
    sigma = Sigma(R)
    sigma = np.array(sigma)
    print("Hello", i)
    for i,( r_i,s_i), in enumerate(zip(R,sigma)):
        X_guess_disk = np.array([X_cg[i], H_cg[i], h_cg[i], Ti[i],impinge_angle[i]])
        print(X_guess_disk[0])
        X_cg[i], H_cg[i], h_cg[i], Ti[i],impinge_angle[i] =fsolve(Disk,X_guess_disk, args = (r_i, s_i, psi_i, psi_s))
        x_guess_last_disk = np.array([X_cg[i], H_cg[i], h_cg[i], Ti[i],impinge_angle[i]])
        
        print(x_guess_last_disk[0])
        #print(np.linalg.norm(x_guess_last_disk-X_guess_disk))
        epsilonnnn = 1
        
        #######################HELP####################################
        if (np.lina.norm(x_guess_last_disk - X_guess_disk)) < epsilonnnn:
            x_guess_last_disk = X_guess_disk
    zi = np.linspace(0.0, H_cg, nz +1)
    z_elevation   = 0.5 * (zi[0:nz] + zi[1:nz+1])          # Take the in-between values
    z_elevation = np.array(z_elevation)
    #z_elevation   




fig, ax = plt.subplots()
ax.plot(R, H_cg)
ax.set_xlabel("Radius")
ax.set_ylabel("Surface_Height")
ax.set_title("Surface_Height vs. Radius")
plt.show()


fig, ax = plt.subplots()
ax.plot(R, Ti)
ax.set_xlabel("Radius")
ax.set_ylabel("Temp")
ax.set_title("Temp vs. Radius")
plt.show()




#Pressure [dyn/cm^2] = [g cm/s^2/cm^2] = [g/s^2/cm]
def Pressure(rho, T_disk):
    return (kb*T_disk*rho/mu/mp)
pres =Pressure(rho, Ti)

def T_surface(R, epsilon):
    return epsilon**-0.25 * (rstar/2/R)**0.5 * tstar
Ts = T_surface(R,1)







"""
def closest(lst, K):
     lst = np.asarray(lst)
     idx = (np.abs(lst - K)).argmin()
     return lst[idx]
  
#Planck opacity [cm^2/g]
def planckOpacity(rho,T_disk,H_cg ,pres):
    df = pd.read_csv("Planck_data.csv")
    planck_opacity = []
    planck_opacity = np.array(planck_opacity)
    for i,(T_i,H_i) in enumerate(zip(T_disk, H_cg)):
        print(T_i)
        filtering_data = ( df[
            (df["T_stellar"] == (closest(df["T_stellar"], tstar))) 
            & (df["T_gas"] == (closest(df["T_gas"],T_disk[T_i])))
           # & (df["Rho_gas"] == (closest(df["Rho_gas"],rho[T_i][H_cg] )))
           # & (df["P_gas"] == (closest(df["P_gas"],  pres)[T_i][H_i]))
           ])
        filtering_data_array = filtering_data.to_numpy()
        planck_opacity_data = filtering_data_array[0][5]
        planck_opacity = np.append(planck_opacity, planck_opacity_data)
    return planck_opacity

planck_opacity = planckOpacity(rho,Ts,H_cg ,pres)
   
 """


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



def kp_interp(T_guess):
    """Interpolate Temperature to find accurate Kp(T) value."""
    kp_t_guess = np.interp(T_guess, T_data, kp_data)
    #kp_t_guess = np.interp(T_guess, T_data, Kp_data)

    return kp_t_guess


def t_s_eq(T_iter, const):
    """Solve iteration for Ts."""
    eps_iter = (kp_interp(tstar) / kp_interp(T_iter)) ** 0.25
    return eps_iter * const - T_iter


def t_i_eq(T_iter, const):
    """Solve iteration for Ti."""
    psi_i_iter = kp_interp(T_iter)
    return const * (psi_i_iter ** -0.25) - T_iter


kp_tstar = kp_interp(tstar)
T_s = np.zeros(len(R))
T_s_0 = 3000.0
Psi_s = np.zeros(len(R))
alpha = impinge_angle
T_i = np.zeros(len(R))
T_i_0 = 3000.0
Psi_i = np.zeros(len(R))
T_prev = t_sublimation
for i, r_i in enumerate(R):
    # Numerically solve for T_s
    const_1 = np.sqrt(0.5 * rstar / r_i) * tstar
    print(const_1)
    T_s[i] = fsolve(t_s_eq, T_prev, args=const_1)
    #print(T_s[i])
    T_prev= T_s[i]
    # Numerical solution for T_s generates Psi_s
    #Psi_s[i] = kp_interp(T_s[i])
    Psi_s[i] = (T_s[i]/const_1)**4 *kp_tstar
    
    # Simplify constant elements of equation for T_i
    const_2 = ((alpha * Psi_s[i]) ** 0.25) * np.sqrt(rstar / r_i) * tstar
    # Numerically solve for T_i via Psi_i
    T_i[i] = fsolve(t_i_eq, T_i_0, args=const_2)
    # Numerical solution for T_i generates Psi_i
    Psi_i[i] = kp_interp(T_i[i])


fig, ax = plt.subplots()
ax.plot(R, T_s)
ax.set_xlabel("Radius")
ax.set_ylabel("Surface_Height")
ax.set_title("Surface_Height vs. Radius")
plt.show()






