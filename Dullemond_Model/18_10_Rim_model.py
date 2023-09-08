import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import fsolve
import pandas as pd
import math

#Natural constant
GG =  6.67430e-8  # Universal gravitational constant 6.67×10−8cm3g−1s−2
kb =  1.380649e-16  # Boltzmann constant defined 
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
sigma0 = 2*10 **3  # Surface density at 1 AU  Σ0 =10^3 gcm^1
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

R_rim, h_rim, H_rim,X_rim, sigma_rim = fsolve(Rim, (0.52*AU, 0.1*AU,0.057*AU,5.3,6207))
print(R_rim/AU,H_rim/R_rim, X_rim)


#Definig grid system in radial direction
rin = R_rim  # Sublimation radius  # inner rim of disk. Density is assumed to be 0 inside [cm] 
rout = 2*AU
#rout = 2.7 * (10 ** 2) * AU  # Outer edge of disk. Density is set to 0 outside this [cm] 
ri = np.linspace(rin,rout,nr+1) 
#R = np.array([0.5 * (ri[i:] + ri[:i]) for i in range(len(ri))])  # Take the in-"between values
#R   = np. array(0.5 * ( ri[::-1] + ri) ) 
R   = 0.5 * ( ri[0:nr] + ri[1:nr+1] )             # Take the in-"between values
R = np.array(R)


def planckOpacity(rho,T_disk,H_cg pres):
    def closest(lst, K):
         lst = np.asarray(lst)
         idx = (np.abs(lst - K)).argmin()
         return lst[idx]
          
    df = pd.read_csv("Planck_data.csv")
    
    planck_opacity = []
    planck_opacity = np.array(planck_opacity)
    for i,(T_i,H_i) in enumerate(zip(T_disk, H_cg)):
        filtering_data = ( df[
            (df["T_stellar"] == (closest(df["T_stellar"], tstar))) 
            & (df["T_gas"] == (closest(df["T_gas"],T_disk[T_i])))
            & (df["Rho_gas"] == (closest(df["Rho_gas"],rho[T_i][H_cg] )))
            & (df["P_gas"] == (closest(df["P_gas"],  pres)[T_i][H_i]))
           ])
        filtering_data_array = filtering_data.to_numpy()
        planck_opacity_data = filtering_data_array[0][5]
        planck_opacity = np.append(planck_opacity, planck_opacity_data)
    return planck_opacity


def Disk(R, vars):
    X_cg, impinge_angle, sigma, H_cg, h_cg, Ti, Ts, phi_i, phi_s, epsilon =vars
    eqn1 = 1-math.erf(X_cg/np.sqrt(2)) -2* impinge_angle/(sigma* kp_stellar)
    eqn2 = X_cg*h_cg -H_cg
    eqn3 = sigma0*(R/AU)**beta - sigma
    eqn4 = (Ti/t_virial)**0.5 * (R/rstar)**0.5 -h_cg/R
    eqn5 = (impinge_angle*psi_s/psi_i)**0.25 * (rstar/R)**0.5 *tstar -Ti
    eqn6 = (1/epsilon)**0.25 *(rstar/(2R))**0.5 *tstar -Ts
    eqn7 = (0.4*rstar/R) +(gamma-1)*H_cg/R -impinge_angle
    eqn8 = kp_stellar/
    
    
    
def Disk(R,psi_i, psi_s):

    def Temp(Tii, R):
       #Calculate numerical result for temperature.
        return ((((0.4 * rstar / R) + ((gamma - 1) * ( X_cg* ((kb*Tii*R**3/(mu*mp*GG*mstar))**0.5)) / R))*psi_s/2/psi_i)**0.25 *(rstar/R)**0.5)*tstar -Tii
       #return (((0.4*rstar/R) + (gamma - 1)*(X_cg*(Tii/t_virial)**0.5  * (R/rstar)**0.5) /R)**0.25 * (rstar/R)**0.5 *(psi_s/psi_i)**0.25 * tstar) - Tii
    
    Ti = np.zeros(len(R))
    for i, r_i in enumerate(R):
        Ti[i] = fsolve(Temp, t_sublimation, args=r_i)
    Ti = np.array(Ti)
    
    h_cg = (Ti * R / (rstar * t_virial)) ** 0.5 *R
    H_cg = X_cg * h_cg
    impinge_angle = (0.4 * rstar / R) + ((gamma - 1) * H_cg / R)
    
    
    return Ti, h_cg, H_cg, impinge_angle

Ti, h_cg, H_cg, impinge_angle = Disk(R,psi_i,psi_s)

#surface_density =np.array(R)
surface_density = (1/2*(psi_i/psi_s)* (Ti/Ts)**4)

zi = np.linspace(0.0, H_cg, nz +1)
z_elevation   = 0.5 * (zi[0:nz] + zi[1:nz+1])          # Take the in-between values
z_elevation = np.array(z_elevation)
#z_elevation   = np.array([0.5 * ( zi[i] + zi[i+1])   for i in range(len(zi) - 1)])  # Take the in-"between values       # Take the in-between values



"""
fig, ax = plt.subplots()
#ax.plot(R, Ti)
plt.plot(R/AU, Ti)

ax.set_xlabel("Radius")
ax.set_ylabel("Temperature")
ax.set_title("Protoplanetary System Temp vs. Radius")

fig, ax = plt.subplots()
ax.plot(R/AU, H_cg)
ax.set_xlabel("Radius")
ax.set_ylabel("Surface_height")
ax.set_title("Surface_height vs. Radius")
plt.show()
"""
#R_fl =np.where(H_cg/R < H_rim/R_rim)     
#print(R[R_fl])

#R_flare = H_cg*(R_rim/ H_rim)

#Inclination from the mid point
theta = ((R/R_rim -1)/(z_elevation/H_rim -1))
theta = np.array(theta)
#Reference to the nearest star
delta = (R/R_rim -1)/ (z_elevation/H_rim -1)
delta = np.array(delta)



#Surface density
sigma = np.array(sigma0 *(R/AU)**-1.5)
"""
fig, ax = plt.subplots()
ax.plot(R/AU, sigma)
ax.set_xlabel("Radius")
ax.set_ylabel("Surface_density")
ax.set_title("Surface_density vs. Radius")
plt.show()
"""
rho = sigma * np.exp(-z_elevation**2/(2*H_cg**2)) /((np.sqrt(2*np.pi) * H_cg))

#Pressure
P = R*rho*Ti/mu
Ra, Za = np.meshgrid(R,z_elevation)
rhoo = (sigma0 *(Ra/AU)**-1.5) * np.exp(-Za**2/(2*(X_cg * (((((1 / epsilon) ** 0.25) * np.sqrt(rstar / (2 * Ra)) * tstar) * Ra / (rstar * t_virial)) ** 0.5 *Ra)
)**2)) /((np.sqrt(2*np.pi) * (  X_cg * (((((1 / epsilon) ** 0.25) * np.sqrt(rstar / (2 * Ra)) * tstar) * Ra / (rstar * t_virial)) ** 0.5 *Ra)
)))
"""
plt.contour(Za,Ra,rhoo)
plt.colorbar()
plt.ylim(0,0.5* 1e13)

plt.show()
"""
#F_irr = impinge_angle *(lstar/(4*np.pi* R**2)) + 2*impinge_angle*ss*T_rim**4 *(R_rim/R)**2 /np.pi *np.cos(theta) * (delta * np.sqrt(1- delta) + np.arcsin(delta))

edge_angle_i = 67    # Viewing angle for observer
dist =144 * AU* 206265   #Distamce from the stellar to observer

def Blackbody(nu, T):
    return (2*hp*nu**3/clight**2)*(np.exp(hp*nu/(kb*T))-1)**-1

F_rim = 4*(R_rim* H_rim/ dist**2)*np.sin(math.radians(edge_angle_i))*Blackbody(nu, T_rim)

def luminosity_func(nu, Ti, R):
    for i, T_i in enumerate(Ti):
        print(Blackbody(nu, T_i))

#luminosity_func(nu, Ti)  
    

"""
planck_opacity = np.array(nr+1)
for i, (T_i, T_s, rho_i) in enumerate(zip(Ti, Ts,rho)):
    
    planck_opacity[i]= planckOpacity(rho_i,T_i)
    epsilion = planckOpacity(rho,Ts)
psi_ii = []
psi_ii = np.array(psi_ii)
psi_ss = []
psi_ss = np.array(psi_ss)
psi_ii = planck_opacity*sigma
psi_ss = planck_opacity
"""

#Appendix A2
#Integral of exp(-z^2 /2H_cvg^2)dz = (pi/2)^0.5 x H_ccg (1- erf(1/2^2))
#1 -erf(1/2^0.5) = 1- erf(X_cg/ 2^0.5) = 2 impingement_angle(X_cg)/ (sigma x Kp_stellar)
H_ccg = np.zeros(nr)
def H_cg_refine(H_ccg, R, sigma):
    #Calculate numerical iteration for solving H_ccg
    analytical_sol =(kp_stellar * R/(0.4 * rstar + (gamma -1) * H_ccg) * (sigma/2)) * (1- math.erf(1/2**0.5)) - 1
    return analytical_sol - H_ccg

for i, r_i, in enumerate(R):
    H_ccg[i]= fsolve(H_cg_refine, H_rim, args=(r_i, sigma[i]))


fig, ax = plt.subplots()
ax.plot(R/AU, H_cg)
ax.set_xlabel("Radius")
ax.set_ylabel("Surface_Heightttttttttt")
ax.set_title("Surface_Height vs. Radius")
plt.show()

def intergrand(x):
    return np.exp(-x**2)

def Zo(z,sigma, h_rim):
    integral, err =quad(intergrand,z/h_rim, np.inf)
    const = 8* sigma* kp_stellar/(np.pi)**0.5
    return 1-integral*const

zo = np.zeros(nz)
for i,s_i in enumerate(sigma):
    zo[i] = fsolve(Zo,H_rim, args = (s_i, h_rim))
    

