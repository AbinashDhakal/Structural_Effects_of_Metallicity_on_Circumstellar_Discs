import numpy as np
import sympy as smp
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.optimize import fsolve
import sympy as smp
import math
import pandas as pd


#Natural constant
GG      =  6.67430e-8      # Universal gravitational constant 6.67×10−8cm3g−1s−2
kb     =  1.380649e-16     # Boltzmann constant defined 
hp      =  6.62607015e-27  # Planck constant [erg s] defined
clight  =  2.99792458e10   # light speed [cm s^-1] defined 
NA     =  6.02214076e23    # Avogadro constant defined 
mp  =  1.672621898e-24     # proton mass [g]; not 1.0/NA anymore. 
ss     = 5.6703e-5         # Stefan-Boltzmann const  [erg/cm^2/K^4/s]

#Astronomical constant
AU     =  1.495978707e13   # Distance between Earth and Sun [cm] defined 
pc     =  3.085677581e18   # Parsec [ cm] 
ms     =  1.9884e33        # Sun mass [g] 
ts     =  5.777e3          # effective temperature of the sun [K] 
ls     =  3.828e33         # solar intensity [erg/s] defined by IAU 
rs     =  6.96e10          # solar radius [cm] 

#Gas constant
mu     =  2.353               # Average molecular weight  (H2 + He + metals)

# grid parameters 
nr         =  2560          # Number of radial grids. 
ntheta     =  257        # Number of grids in elevation.
nz         =  2560          # Number of elevation grids. 

#Stellar Parameter
mstar =0.5*ms
rstar =2.5*rs
tstar = 4_000
lstar = 0.7*ls

#Disk parameter
T_rim= t_sublimation = 1500.0    #Temperature at which dust start to condense
sigma0 = 10**3            #Surface density at 1 AU  Σ0 =10^3 gcm^1
opacity_const =400.0      #Monochromatic opacity = 400cm^2 g^-1
psi_i = 1.0
psi_s = 1.0
X_cg = 3.0                # Constant in range 2 - 6
sigma0 = 1e3             #Sigma dust at 1 AU [g/cm2^]

nu = 5.879*10**10 * tstar       #Wien law
t_virial = GG* mstar* mu* mp/ (kb* rstar) #Virial Temp
gamma =2

def Rim():
    def R_rimm(r_i):
        """Set equation for R_rim."""
        return (( np.sqrt(0.25 * lstar / (np.pi  * (t_sublimation ** 4) * ss))) * np.sqrt(1 + X_cg * np.sqrt((kb * t_sublimation / (mu * mp * GG * mstar)) * r_i))) - r_i
    
    # Create initial guess x0 = r_0
    r_0 = AU
    # Numerically solve for R_rim
    R_rim = fsolve(func=R_rimm, x0=r_0, factor=10)  
    
    # Calculate h_rim using eff_r_rim, and subsequent H_rim
    h_rim = np.sqrt(np.sqrt((kb * t_sublimation / (mu * mp * GG * mstar)) * (R_rim ** 3)))
    H_rim = X_cg * h_rim
    return h_rim,H_rim,R_rim

h_rim,H_rim,R_rim =Rim()

#Definig grid system in radial direction
rin        = R_rim     #Sublimation radius  # inner rim of disk. Density is assumed to be 0 inside [cm] 
rout       = 2.7*(10**2)*AU# Outer edge of disk. Density is set to 0 outside this [cm] 
ri        =  np.linspace ( rin , rout , nr + 1 )
#R   = np. array(0.5 * ( ri[:-1] + ri[1:] ) )             # Take the in-"between values
R = np.array([0.5 * (ri[i] + ri[i + 1]) for i in range(len(ri) - 1)])  # Take the in-"between values


#Calculating the ratio of planck mean opacity Eqn 8
def planck_mean_opacity(T):
    integrand1 = lambda x : ((2* nu**2 *opacity_const/ clight**2) * (hp * nu/(np.exp(hp*nu/(kb*T))-1)))
    lambda_rel = hp * clight / (kb * T)     #defining the reference lambda
    lambda_min = lambda_rel / 10  #Defining lower int bound
    lambda_max = lambda_rel * 10  #Defining upper int bound 
    integral,err1 = integrate.quad(integrand1, lambda_min, lambda_max, epsabs=0, limit=50)
    const = np.pi / (ss*T**4)
    planck_opacity = integral*const
    return planck_opacity

epsilon = planck_mean_opacity(tstar)/planck_mean_opacity(tstar)

#epsilion =(8*np.pi*kb*T*r_dust/hc)
#grain sphere with radius R =0.1* 10^-6m
#mass density rho_dust = 2 g cm^-3
#Optical depth: tau_v = 4 X 10^5 R_AU^-3/2 

#Calculating Dust temperature
Ts= (1/epsilon)**0.25 *(rstar/(2*R))**0.5 *tstar

#Calculating interior disk temperature
def disk():
    def Temp(Tii, R):
        """Calculate numerical result for temperature."""
        return (((0.4*rstar/R) + (gamma - 1)*(X_cg*(Tii/t_virial)**0.5  * (R/rstar)**0.5) /R)**0.25 * (rstar/R)**0.5 *(psi_s/psi_i)**0.25 * tstar) - Tii
    Ti = np.zeros(len(R))
    
    for i, r_i in enumerate(R):
        Ti[i] = fsolve(Temp, t_sublimation, args=r_i)
    Ti = np.array(Ti)
    
    #---------------------------------------------------------------------------------
    #**********************FIX IT**********************************
    h_cg = []
    h_cg = np.array(h_cg)
    
    h_cg = (Ti/t_virial)**0.25 * (R/rstar)**0.25
    columns = list(zip(*h_cg))
    h_cg = columns[0]
    h_cg = np.array(h_cg)
    
    H_cg = X_cg * h_cg
    impinge_angle = (0.4*rstar/R) + (gamma - 1) * H_cg/R
    return Ti, h_cg, H_cg, impinge_angle

Ti, h_cg, H_cg, impinge_angle = disk()

fig, ax = plt.subplots()
ax.plot(R, Ti)
ax.set_xlabel("Radius")
ax.set_ylabel("Temperature")
ax.set_title("Protoplanetary System Temp vs. Radius")

fig, ax = plt.subplots()
ax.plot(R, H_cg)
ax.set_xlabel("Radius")
ax.set_ylabel("Surface_height")
ax.set_title("Surface_height vs. Radius")
plt.show()

#surface_density =np.array(R)
#surface_density = (1/(2*planck_mean_opacity(Ts)))*(psi_i/psi_s)* (Ti/Ts)**4


zi        =  np.linspace ( 0 , H_cg , nz + 1 )
#z_elevation   = np.array(0.5 * ( zi[:-1] + zi[1:] )    )          # Take the in-between values
z_elevation   = np.array([0.5 * ( zi[i] + zi[i+1])   for i in range(len(zi) - 1)])  # Take the in-"between values       # Take the in-between values



#Inclination from the mid point
theta = ((R/R_rim -1)/(z_elevation/H_rim -1))
#Reference to the nearest star
delta = (R/R_rim -1)/ (z_elevation/H_rim -1)

#F_irr = impinge_angle * (lstar/(4*np.pi* R**2)) + 2*impinge_angle*ss*T_rim**4 *(R_rim/R)**2 /np.pi *np.cos(theta) * (delta * np.sqrt(1- delta) + np.arcsin(delta))


#Surface density
sigma = np.array(sigma0 *(R/AU)**-1.5)

rho = sigma * np.exp(z_elevation**2/(2*H_cg)) /(2*H_cg)


def closest(lst, K):
      
     lst = np.asarray(lst)
     idx = (np.abs(lst - K)).argmin()
     return lst[idx]
      

df = pd.read_csv("Planck_data.csv")

planck_opacity = []
planck_opacity = np.array(planck_opacity)
for T_i in range (len(Ti)):
    filtering_data = ( df[
        (df["T_stellar"] == (closest(df["T_stellar"], T_i))) 
        & (df["T_gas"] == (closest(df["T_gas"],600)))
        #& (df["Rho_gas"] == (closest(df["Rho_gas"],  3.7-20)))
        #& (df["P_gas"] == (closest(df["P_gas"],  1)))
       ])
    planck_opacity = filtering_data.to_numpy()[T_i][5]
