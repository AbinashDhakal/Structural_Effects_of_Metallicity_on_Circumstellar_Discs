import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfinv
import pandas as pd

#Natural constant
GG =  6.67430e-8  # U# gravitational constant  [dyn cm^2 g^-2]
kb =  1.380649e-16  # Boltzmann constant defined [erg K^-1]
hp =  6.62607015e-27  # Planck constant [erg s] defined
clight = 2.99792458e10  # light speed [cm s^-1] defined 
NA = 6.02214076e23  # Avogadro constant defined 
mp = 1.672621898e-24  # proton mass [g]; not 1.0/NA anymore. 
ss = 5.6703e-5  # Stefan-Boltzmann const  [erg/cm^2/K^4/s]
gasR = 	8.3143435* 10**7 #erg K-1 mole-1
a =  ss/0.25/clight
#Astronomical constant
AU = 1.495978707e13  # Distance between Earth and Sun [cm] defined 
pc = 3.085677581e18  # Parsec [ cm] 
ms = 1.9884e33  # Sun mass [g] 
ts = 5.777e3  # effective temperature of the sun [K] 
ls = 3.828e33  # solar intensity [erg/s] defined by IAU 
rs = 6.96e10  # solar radius [cm] 

#Gas constant
mu = 2.353  # Average molecular weight  (H2 + He + metals)
adia_gamma = 1.42 # adiabatic index
Pr = 1 #Prandtl number


# grid parameters 
nr = 1000  # Number of radial grids. 
ntheta = 257  # Number of grids in elevation.
nz = 256 # Number of elevation grids. 

#Stellar Parameter
mstar = 2.5 * ms
rstar = 2.5*rs
tstar = 9500
lstar = 47* ls

#Disk parameter
T_rim= t_sublimation = 1500.0  # Temperature at which dust start to condense
#sigma0 = 1000 **2  # Surface density at 1 AU  Σ0 [g/cm^2] 
#sigma0 = 17000
#sigma0 = 5*10**4
sigma0 = 2*10**3

kp_stellar = 400.0  # Monochromatic opacity = 400cm^2 g^-1

#kP(T*) = 2100 cm2 g-1;   https://iopscience.iop.org/article/10.3847/1538-4357/835/2/230/pdf 
#Kp(T_rim) = 700 cm2 g-1
#kp_stellar = 0.4
nu = 5.879 * (10 ** 10) * tstar  #Wien law
t_virial = GG * mstar * mu * mp / (kb * rstar)  # Virial Temp
gamma = 2.0
beta = -1.5
gray_const =1
visc_const = 10**-2

A_sil = 10**(-25.10)   #A_sil = 10^-25.10 cm^2.5
a_min = 5*10**-7
a_max = 2.5 * 10**-5
sil_density = 3.4 # g cm^-3
#################################]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]#############
#dust to gas mass ratios ζgrap = 0.0025 and ζsil = 0.0034
#(Draine & Lee 1984). 

num_sil_density = sil_density/mp *A_sil*(a_max**-2.5 - a_min**-2.5)/(-3.5)
kv = (3/4)*(1.5*10**4)*(10**(-25.10)) * ((a_max**-2.5 - a_min**-2.5)/(-3.5))
ni_const = (3/4)*(A_sil)*((a_max**-2.5 - a_min**-2.5)/(-3.5)) /(10**-6)

fd2g = 10**-2    # The gas: dust is assumed be 100 : 1.

import pandas as pd
def reading_dust_data():
    df = pd.read_csv("Planck_dust_data.csv")
    T_data = np.array(list(df["T_rad"]))
    Qabs_data = np.array(list(df["<Qabs>/a"]))
    kp_t_array = np.zeros((2, len(Qabs_data)))
    kp_t_array[0] =T_data
    kp_t_array[1] =Qabs_data
    #change_array = T_kp_data.to_numpy()
    return T_data, Qabs_data

T_data, Qabs_data = reading_dust_data()
def Qabs_interp(T_guess):
    T_data, Qabs_data = reading_dust_data()
    logT = np.log10(T_data)
    logQ = np.log10(Qabs_data)
    logT_interp = np.log10(T_guess)
    return np.power(10.0, np.interp(logT_interp, logT, logQ))
kp_stellar = Qabs_interp(tstar)*(3/4)*(1/sil_density)*10**4


def surface_density(R):
    return sigma0*(R/AU)**beta

def pressure_height(T, R):
    return ((T/t_virial)**0.5 * (R/rstar)**0.5 * R)

def scale_height(X_cg, h_cg):
    return (X_cg *h_cg)

def chi_cg(sigma,kp_stellar):
    return (erfinv(1-1/(4*sigma*kp_stellar)))
    #return 1-  (4*sigmaRim*kp_stellar)*(1-math.erf(X_rim_iter))

def Rim_radial(H_rim, Rrim_guess, gray_const):
    return gray_const**0.5*(lstar/(4*np.pi*T_rim**4 *ss))**0.5 *(1+ H_rim/Rrim_guess)**0.5


def density_distrib(z_elevation, h_cg, sigma):
    return sigma*np.exp(- z_elevation**2 / 2/h_cg**2) / (np.sqrt(2*np.pi) *h_cg)



R_rim =  0.5399999947476138 *AU
rin = R_rim # Sublimation radius  # inner rim of disk. Density is assumed to
rout = 100*AU
R =ri = np.linspace(rin,rout,nr) 
#R   = 0.5 * ( ri[0:nr] + ri[1:nr+1] )             # Take the in-"between values
R = np.array(R)



kp_stellar =  1586.4162432323064
kp_Trim =973.8276480930514


# Set maximum iteration
maxIter = 500

T = np.zeros(nr)
N = int((rout - rin)/nr)
#rout = 100
#rin = 0.5
#N measure number of grid?

#Running 500 iteration hoping for convergence
Q_visc = np.zeros(nr)
"""
dr = (R[1]-R[0]) *AU  
dr =30
for iteration in range(0, maxIter):
    T = (0.7)*(rstar/R)**0.5 *tstar
    for i in range(1,len(T)-1):
        T[0] = 1500
        print(i)

        #print(R[i]/AU)
        #Calcuating given paramater for ith element with T[i] and R[i]
        #T and R has same array size len(R) = len(T)
        hp = (T[i]/t_virial)**0.5 *(R[i]/rstar)**0.5 *R[i] #pressure scale
        cs = np.sqrt(kb*adia_gamma*T[i]/(mu*mp)) #sound speed
        omega = np.sqrt(GG*mstar/R[i]**3) #Keplerian velocity
        cV = kb/(mu*mp*(adia_gamma-1)) #specific const volume 
        cP = adia_gamma /cV #specific const pressure
        vt = visc_const*cs**2/omega #viscosity
        sigma = surface_density(R[i]) #surface density 
        rho = sigma/(np.sqrt(2*np.pi* hp**2)) #Density 
        cV = kb/(mu*mp*(adia_gamma-1)) #specific const volume 
        kt = rho*vt*adia_gamma*cV #Turblent thermal conductivity
        #Discretising laplace equation

        dT_dR_2 = (T[i+1] -2*T[i] +T[i-1])/(dr)**2
        dT_dR = (T[i+1] -T[i])/(dr)/2
        #Q_cool = 2 *ss* T
        #Q_condct = kt * ΔT
        #Q_cool = Q_conduct
        #T[i] =(kt*(dT_dR_2)/(2*ss))**0.25

        T[i] =(kt*(dT_dR_2 +dT_dR/R[i])/(2*ss))**0.25
        Q_visc[i] = rho*vt*(9*GG*mstar/4/R[i]**3)
    
"""   
#T = (0.7)*(rstar/R)**0.5 *tstar
dr = (R[1]-R[0])
  
def Q_conduction(R,T,dr):
    T = (0.7)*(rstar/R)**0.5 *tstar

    for i in range(1,len(T)-1):
        T[0] = 1500
        #print(R[i]/AU)
        #Calcuating given paramater for ith element with T[i] and R[i]
        #T and R has same array size len(R) = len(T)
        hp = (T[i]/t_virial)**0.5 *(R[i]/rstar)**0.5 *R[i] #pressure scale
        cs = np.sqrt(kb*adia_gamma*T[i]/(mu*mp)) #sound speed
        omega = np.sqrt(GG*mstar/R[i]**3) #Keplerian velocity
        #hp =cs/omega
        cV = kb/(mu*mp*(adia_gamma-1)) #specific const volume 
        cP = adia_gamma /cV #specific const pressure
        vt = visc_const*cs**2/omega #viscosity
        sigma = surface_density(R[i]) #surface density 
        rho = sigma/(np.sqrt(2*np.pi* hp**2)) #Density 
        cV = kb/(mu*mp*(adia_gamma-1)) #specific const volume 
        cP = adia_gamma /cV #specific const pressure
        kt = rho*vt*adia_gamma*cV #Turblent thermal conductivity
        #Discretising laplace equation

        dT_dR_2 = (T[i+1] -2*T[i] +T[i-1])/(dr)**2
        dT_dR = (T[i+1] -T[i])/dr/2
        Q =kt*(dT_dR_2 +dT_dR/R[i])
        #Q_cool = 2 *ss* T
        #Q_condct = kt * ΔT
        #Q_cool = Q_conduct
    
            
    return T
T_abc = np.zeros(nr)        
for iteration in range(0, maxIter):
    T_abc =Q_conduction(R,T,dr)



def Q_visc(R):
  
    yrs = 31557600 #s
    M_acc = 10**-8*ms/yrs
    T_visc = (3*GG*mstar*M_acc/(8*np.pi* ss* R**3))**0.25
    #print(T_visc)
    for iteration in range(0, maxIter):
        for i in range(1,len(T)-1):    
            sigma = surface_density(R[i])
            cs = np.sqrt(kb*T_visc[i]/mu/mp)
            Q_visc =(sigma/2/np.pi)*visc_const*cs*(GG*mstar/4/R[i]**3)
            T_visc[i] = (Q_visc/2/ss)**0.25
        print(iteration)
    return T_visc

T_cond =Q_conduction(R,T,dr)
T_visc = Q_visc(R)
fig, ax = plt.subplots()
ax.semilogx(R/AU, T_cond)


ax.set_xlabel("Radius")
ax.set_ylabel("Temperature")
ax.set_title("Temperature vs. Radius")
plt.show()   

fig, ax = plt.subplots()
ax.semilogx(R/AU, T_visc)
ax.set_xlabel("Radius")
ax.set_ylabel("Q_visc")
ax.set_title("Q_visc vs. Radius")
plt.show()   
