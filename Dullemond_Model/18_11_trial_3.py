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
rstar = 0.5 * rs
tstar = 9500
lstar = 47* ls

#Disk parameter
T_rim= t_sublimation = 1500.0  # Temperature at which dust start to condense
sigma0 = 2*10 **3  # Surface density at 1 AU  Σ0 [g/cm^2] 
kp_stellar = 2100.0  # Monochromatic opacity = 400cm^2 g^-1

#kP(T*) = 2100 cm2 g-1;   https://iopscience.iop.org/article/10.3847/1538-4357/835/2/230/pdf 
#Kp(T_rim) = 700 cm2 g-1

nu = 5.879 * (10 ** 10) * tstar  #Wien law
t_virial = GG * mstar * mu * mp / (kb * rstar)  # Virial Temp
gamma = 2.0
beta = -1.5
gray_const =1
yrs = 31557600 #s
M_acc = 10**-8*ms/yrs

A_sil = 10**(-25.10)   #A_sil = 10^-25.10 cm^2.5
a_min = 5*10**-7
a_max = 2.5 * 10**-5
sil_density = 3.4 # g cm^-3
#################################]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]#############
#dust to gas mass ratios ζgrap = 0.0025 and ζsil = 0.0034
#(Draine & Lee 1984). 

nH = 0.22 #0.22 cm**-3  number density of H nuclei
dust_radius = 1000    #cm

num_sil_density = sil_density/mp *A_sil*(a_max**-2.5 - a_min**-2.5)/(-3.5)
kv = (3/4)*(1.5*10**4)*(10**(-25.10)) * ((a_max**-2.5 - a_min**-2.5)/(-3.5))
ni_const = (3/4)*(A_sil)*((a_max**-2.5 - a_min**-2.5)/(-3.5)) 



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

def Rho(z_elevation, h_cg, sigma):
    return sigma*np.exp(- z_elevation**2 / 2/h_cg**2) / (np.sqrt(2*np.pi) *h_cg)

def converging_rim(gray_const):
    R_guess = 0.53*AU
    R_rim = 0.54*AU
    H_rim =0.1162*AU
    H_guess =0.11*R_rim
    iteration =0
    kp_stellar = 1200
    while np.round(H_guess/R_rim,2) != np.round(H_rim/R_rim,2):
        H_guess =        H_rim
        w = 0.9999999
        error = 1*AU
        while error >10e-5 *AU:
            temp_r = R_rim
            R_rim = (w*R_rim +(1-w)*Rim_radial(H_rim, R_rim,gray_const))
            error = np.max(np.abs(temp_r-R_rim))
            #print(R_rim/AU, temp_r/AU)
        iteration +=1
        #print("I am running iteration", iteration)
        #print("This is dust_const", gray_const)
        sigma_rim = surface_density(R_rim)
        h_rim = ((T_rim/t_virial)**0.5 * (R_rim/rstar)**0.5 * R_rim)
        X_rim = chi_cg(sigma_rim,kp_stellar)
        H_rim = scale_height(X_rim, h_rim)
        #print(H_rim/R_rim, H_guess/R_rim)
        Rho_rim = Rho(0, h_rim, sigma_rim)
        kp_Trim = Qabs_interp(T_rim)*(1/Rho_rim)*ni_const
        kp_stellar = 3*(kp_Trim)
    return R_rim, H_rim, h_rim, sigma_rim, X_rim, kp_stellar



R_rim, H_rim, h_rim, sigma_rim, X_rim, kp_stellar = converging_rim(gray_const)
#print(R_rim/AU,H_rim/R_rim)   

#Definig grid system in radial direction
rin = R_rim # Sublimation radius  # inner rim of disk. Density is assumed to
rout = 100*AU
R =ri = np.linspace(rin,rout,nr) 
#R   = 0.5 * ( ri[0:nr] + ri[1:nr+1] )             # Take the in-"between values
R = np.array(R)


def viscous_temp(R):
    return (3*GG*mstar*M_acc/(8*np.pi* ss* R**3))**0.25

dr = R[1] - R[0]
def Rim_heating(R,kp_rim, kp_T, dr):
    return T_rim*kp_rim*np.exp(-(1-0.5**(0.25))*(R_rim - R)/(dr))*kp_T

def alpha_angle(H_cg,R):
    return (0.4*rstar/R   +(gamma - 1)*H_cg/R)

def Rho(z_elevation, h_cg, sigma):
    return sigma*np.exp(- z_elevation**2 / 2/h_cg**2) / (np.sqrt(2*np.pi) *h_cg)

#Pressure [dyn/cm^2] = [g cm/s^2/cm^2] = [g/s^2/cm]
def GasPressure(rho, T_disk):
    return (kb*T_disk*rho/mu/mp)

#fraction of stellar flux that will be absorbed by the interior
def flux_absorbed_fraction(sigma, kp_T): #psi_s
    if (sigma*kp_T) < 1:
        return kp_T
    else:
        return 1
         
#The possiblity that the disk interior is not fully optically thick to its own emission
def self_absorbed_fraction(sigma,kp_T): #psi_i
    if (sigma*kp_T) < 1:
        return sigma*kp_T
    else:
        return 1

def irr_disk_temp(impinge_angle, R,psi_i, psi_s): #Ti
    return (impinge_angle*psi_s/psi_i)**0.25 * (rstar/R)**0.5 *tstar


def irr_surface_Temp(R, kp_Ts, kp_stellar): #Ts
    epsilon =kp_Ts/kp_stellar
    T_s = (1/epsilon)**0.25 * (rstar/2/R)**0.5 *tstar
    return T_s



rho_rim = Rho(0, h_rim, sigma_rim)
kp_rim = Qabs_interp(T_rim)*(2/rho_rim)*(ni_const)



X_cg = np.zeros(nr)
impinge_angle = np.zeros(nr)
H_cg = np.zeros(nr)
h_cg = np.zeros(nr)
sigma = surface_density(R)
Ti = np.zeros(nr)


#convering iteration such that the it iterate until the error is smaller than 0.00001
#THis only check for a radial point
alpha_guess = 0.2

def converging_disk(alpha_guess, r_ii):
    a_i = 0.22
    psi_i = 1
    psi_s = 1 
    kp_Ts = (1/3)*kp_stellar
      
    while np.round(a_i, 5) != np.round(alpha_guess ,5):
    #while alpha_guess !=a_i:
        alpha_guess = a_i
        s_i = surface_density(r_ii)
        T_i =irr_disk_temp(alpha_guess, r_ii,psi_i, psi_s)
        T_s = irr_surface_Temp(R, kp_Ts, kp_stellar)
        X_i = chi_cg(s_i,kp_stellar)
        h_i = pressure_height(T_i,r_ii)
        H_i = X_i*h_i
        rho_dust = Rho(0, h_i, s_i)
        kp_T = Qabs_interp(T_i)*(1/rho_dust)*ni_const
        kp_Ts =Qabs_interp(T_s)*(1/rho_dust)*ni_const
        psi_i = self_absorbed_fraction(s_i,kp_T)
        psi_s = flux_absorbed_fraction(s_i, kp_T)
        a_i =alpha_angle(H_i, r_ii)
    return a_i, T_i


def converging_disk2(R):
    for i, (r_i) in enumerate (R):
        impinge_angle[i], Ti[i] = converging_disk(alpha_guess, r_i)
        print(Ti[i])
        sigma[i] = surface_density(r_i)
        X_cg[i] = chi_cg(sigma[i],kp_stellar)
        h_cg[i] = pressure_height(Ti[i],r_i)
        H_cg[i] = X_cg[i]*h_cg[i]
    return impinge_angle,Ti,sigma,X_cg,h_cg,H_cg

impinge_angle,Ti,sigma,X_cg,h_cg,H_cg =converging_disk2(R)

z_elevation =zi = np.linspace(0.0, H_cg, nz)
#z_elevation = np.array(z_elevation)
#z_elevation   = np.array([0.5 * ( zi[i] + zi[i+1])   for i in range(len(zi) - 1)])  # Take the in-"between values       # Take the in-between values    
Z = np.linspace(0.0, 50*AU, nz)

h_cg[0] = h_rim 
H_cg[0] = H_rim
X_cg[0] = X_rim



Rflare = np.where (R >H_cg*R_rim/H_rim  )#and H_cg*R_rim/R/H_rim >1)
Rshadow = np.where (R <H_cg*R_rim/H_rim)

Ti[Rflare] = 0

h_cg[Rflare] = 0
H_cg[Rflare] = 0
"""

T_effective = (Ti**4 )**0.25
h_cg = pressure_height(T_effective, R)
X_cg = chi_cg(sigma,kp_stellar)
H_cg = scale_height(X_cg, h_cg)
"""

h_cg[0] = h_rim 
H_cg[0] = H_rim
X_cg[0] = X_rim




















fig, ax = plt.subplots()
ax.semilogx(R/AU, Ti)
ax.set_xlabel("Radius")
ax.set_ylabel("Interior_Temperature")
ax.set_title("Interior_Temperature vs. Radius")
plt.show()   




fig, ax = plt.subplots()
ax.semilogx(R/AU, impinge_angle)


ax.set_xlabel("Radius")
ax.set_ylabel("Impinge_angle")
ax.set_title("Impinge_angle vs. Radius")
plt.show()   





fig, ax = plt.subplots()
ax.semilogx(R/AU, H_cg/R)
ax.axhline(y=H_rim/R_rim, xmin=0, xmax=100, color = 'b', linestyle = '--')

ax.set_xlabel("Radius")
ax.set_ylabel("Surface_Height")
ax.set_title("Surface_Height vs. Radius")
#ax.set_ylim([0, 0.5])

plt.show()

fig, ax = plt.subplots()
ax.semilogx(R/AU, h_cg/R)

ax.set_xlabel("Radius")
ax.set_ylabel("Pressure_Height")
ax.set_title("pressure_Height vs. Radius")
ax.set_ylim([0, 0.5])

plt.show()



fig, ax = plt.subplots()
ax.semilogx(R/AU, X_cg)


ax.set_xlabel("Radius")
ax.set_ylabel("X_cg")
ax.set_title("x_cg vs. Radius")
plt.show()   

 


fig, ax = plt.subplots()
ax.plot(R/AU, H_cg*R_rim/R/H_rim)


ax.set_xlabel("R_cg")
ax.set_ylabel("Flaring_point")
ax.set_title("Flaring_point vs. Radius")
plt.show()   
 
    




