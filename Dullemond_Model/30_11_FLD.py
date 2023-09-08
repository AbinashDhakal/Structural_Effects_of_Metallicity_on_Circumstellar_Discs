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

# grid parameters 
nr = 150  # Number of radial grids. 
ntheta = 257  # Number of grids in elevation.
nz = 256 # Number of elevation grids. 

#Stellar Parameter
mstar = 2.5 * ms
rstar = 2.5*rs
tstar = 4000
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

def converging_rim(gray_const):
    R_guess = 0.53*AU
    R_rim = 0.54*AU
    H_rim =0.1162*AU
    H_guess =0.11*R_rim
    iteration =0
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
    return R_rim, H_rim, h_rim, sigma_rim, X_rim, 


R_rim, H_rim, h_rim, sigma_rim, X_rim = converging_rim(gray_const)
#print(R_rim/AU,H_rim/R_rim)   



#Definig grid system in radial direction
rin = R_rim # Sublimation radius  # inner rim of disk. Density is assumed to
rout = 1000*AU
R =ri = np.linspace(rin,rout,nr) 
#R   = 0.5 * ( ri[0:nr] + ri[1:nr+1] )             # Take the in-"between values
R = np.array(R)

kp_Trim = Qabs_interp(T_rim)*(3/4)*(1/sil_density)*10**4
kp_const=kp_stellar/kp_Trim
#gray_const = kp_stellar/kp_Trim*(2/np.exp(1) +kp_stellar/ kp_Trim*(1+2*kp_stellar/kp_Trim))
gray_const = (2*kp_const +1+ 2/kp_const/np.exp(-1))

def viscous_temp(R):
    yrs = 31557600 #s
    M_acc = 10**-8*ms/yrs
    return (3*GG*mstar*M_acc/(8*np.pi* ss* R**3))**0.25
    #return 107*(R/AU)**-0.75 * (mstar/2.5/ms)**0.25 *(M_acc/(10**-8 *ms/yrs))**0.25
T_visc = viscous_temp(R)

betaaa = 1
def Rim_heating(R):
    #return T_rim*np.exp(-1*(1-0.5**(1/(4+betaaa)))  *(R-R_rim)) 
    return T_rim*np.exp(-1*(1-0.5**(1/(4+betaaa)))  *(R-R_rim)/(5834903615946.781))
RimHeating = Rim_heating(R)


"""
dr = R[1] - R[0]
def Rim_heating(R,dr, kp_Trim, kp_T):
    return (kp_Trim/kp_T)**0.25* np.exp(-(1-(0.5)**0.25)*(R-R_rim)/dr)*T_rim


T_guess = 1500
def converging_rimHeating(Ri,T_guess):
    T_i = 1560
    kp_Trim = Qabs_interp(T_rim)*(3/4)*(1/sil_density)*10**4
    print(kp_Trim)
    while np.round(T_i, 5) != np.round(T_guess,5):
        print(T_i, T_guess)
    #while alpha_guess !=a_i:
        T_guess = T_i
        kp_T = Qabs_interp(T_guess)*(3/4)*(1/sil_density)*10**4
        T_i = Rim_heating(Ri,dr, kp_Trim,kp_T)

    print(T_i, T_guess)
    return T_i

"""


def alpha_angle(H_cg,R):
    return (0.4*rstar/R   +(gamma - 1)*H_cg/R)


def Rho(z_elevation, h_cg, sigma):
    return sigma*np.exp(- z_elevation**2 / 2/h_cg**2) / (np.sqrt(2*np.pi) *h_cg)

#Pressure [dyn/cm^2] = [g cm/s^2/cm^2] = [g/s^2/cm]
def GasPressure(rho, T_disk):
    return (kb*T_disk*rho/mu/mp)


def closest(lst, K):
     lst = np.asarray(lst)
     idx = (np.abs(lst - K)).argmin()
     return lst[idx]


#Planck opacity [cm^2/g]
def planckOpacity(Metallicity,T_stellar,T_gas,P_gas,Rho_gas):
    df = pd.read_csv("Planck_data.csv")
    planck_opacity = []
    planck_opacity = np.array(planck_opacity)
    filtering_data = ( df[
        (df["Metallicity"] == (closest(df["Metallicity"], Metallicity))) 
        &(df["T_stellar"] == (closest(df["T_stellar"], T_stellar))) 
        & (df["T_gas"] == (closest(df["T_gas"],T_gas)))
        #& (df["Rho_gas"] == (closest(df["Rho_gas"],Rho_gas)))
        & (df["P_gas"] == (closest(df["P_gas"],  P_gas)))
       ])
    #print( closest(df["T_stellar"], T_stellar), closest(df["T_gas"],T_gas),(closest(df["Rho_gas"],Rho_gas)),(closest(df["P_gas"],  P_gas)))
    #print(filtering_data)
    filtering_data_array = filtering_data.to_numpy()
    #print(filtering_data_array)
    planck_opacity_data = filtering_data_array
    #print(planck_opacity_data[0][5])
    
    return planck_opacity_data[0][5]



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

psi_i, psi_s = 1, 1
def irr_disk_temp(impinge_angle, R,psi_i, psi_s):
    return ((impinge_angle*psi_s/psi_i)**0.25 * (rstar/R)**0.5 *tstar)


def irr_surface_Temp(R, kp_Ts):
    epsilon =kp_Ts/kp_stellar
    T_s = (1/epsilon)**0.25 * (rstar/2/R)**0.5 *tstar
    return T_s
"""

alpha_guess = 0.2

Metallicity = 0.3
T_stellar = 5000

#convering iteration such that the it iterate until the error is smaller than 0.00001
#THis only check for a radial point
def converging_disk(alpha_guess, r_ii):
    a_i = 0.22
    psi_i = 1
    psi_s = 1 
    while np.round(a_i, 5) != np.round(alpha_guess,5):
    #while alpha_guess !=a_i:
        alpha_guess = a_i
        T_i =irr_disk_temp(alpha_guess, r_ii,psi_i, psi_s)
        h_i= pressure_height(T_i,r_ii)
        s_i = surface_density(r_ii)
        X_i = chi_cg(s_i,kp_stellar)
        H_i = X_i*h_i
        a_i =alpha_angle(H_i, r_ii)
        T_s= irr_surface_Temp(r_ii, kp_Ts =3.2)

        w = 0.9999999
        error = 10
        #This check Ti multiple times such that error is less than 0.01
        while error >0.01:
            Temp = T_i
            Rho_gas =Rho(0, h_i, s_i)
            P_gas =GasPressure(Rho_gas, T_i)
            kp_Ts = planckOpacity(Metallicity,T_stellar,T_s,P_gas,Rho_gas)
            kp_Ti = planckOpacity(Metallicity,T_stellar,T_i,P_gas,Rho_gas)
            psi_i = self_absorbed_fraction(s_i, kp_Ti)
            psi_s = flux_absorbed_fraction(s_i, kp_Ts)
            #print("Printing , psi value, ",psi_i, psi_s)
            T_i =  w*Temp +(1-w)*irr_disk_temp(alpha_guess, r_ii,psi_i, psi_s)
            T_s= irr_surface_Temp(r_ii, kp_Ts)
            error = np.max(np.abs(Temp - T_i))
            print(T_i, T_s)
        #print(alpha_guess, a_i)
        
    return a_i



X_cg = np.zeros(nr)
impinge_angle = np.zeros(nr)
H_cg = np.zeros(nr)
h_cg = np.zeros(nr)
sigma = surface_density(R)
Ti = np.zeros(nr)
Ts = np.zeros(nr)
   

#This iterate through every radial points
def converging_disk2(R, psi_i, psi_s,alpha_guess):
    for i, r_i in enumerate(R):
        impinge_angle_iter =converging_disk(alpha_guess, r_i)
        alpha_guess = impinge_angle_iter
        sigma_iter = surface_density(r_i)
        Ti_iter = irr_disk_temp(alpha_guess, r_i,psi_i, psi_s)

        h_cg_iter = pressure_height(Ti_iter,r_i)
        X_cg_iter= chi_cg(sigma_iter,kp_stellar)
        H_cg_iter = X_cg_iter*h_cg_iter
        #impinge_angle_iter =converging(alpha_guess, r_i)    
        #np.append(Ti, Ti_iter)
        Ti[i] = Ti_iter
        h_cg[i] =h_cg_iter
        sigma[i] =sigma_iter
        X_cg[i] =X_cg_iter
        H_cg[i] =H_cg_iter
        impinge_angle[i] =impinge_angle_iter
    return Ti,h_cg,sigma,X_cg,H_cg,impinge_angle







"""















from scipy.integrate import quad
def intergrand(R, Tgas,kp_T):
    return kp_T*(sigma0*(R/AU)**beta)/np.sqrt(2*np.pi)/((Tgas/t_virial)**0.5 * (R/rstar)**0.5 * R)
"""
def Tau(R, Tgas,kp_T):
    return quad(intergrand,np.inf,R,args=(Tgas,kp_T))[0]
"""
"""
def Tau(R, Tgas, kp_stellar):
    return -0.25*np.sqrt(2)*kp_stellar*sigma0/(np.sqrt(np.pi)*(R/AU)**1.5*(R/rstar)**0.5*(Tgas/t_virial)**0.5) + 0.25*np.sqrt(2)*kp_stellar*sigma0/(np.sqrt(np.pi)*(rstar/AU)**1.5*(Tgas/t_virial)**0.5)
"""
def Tau(R, Tgas, kp_stellar):
    ABC =  -0.25*np.sqrt(2)*kp_stellar*sigma0/(np.sqrt(np.pi)*(R/AU)**1.5*(R/rstar)**0.5*(Tgas/t_virial)**0.5) + 0.25*np.sqrt(2)*kp_stellar*sigma0/(np.sqrt(np.pi)*(rstar/AU)**1.5*(Tgas/t_virial)**0.5)
    return np.exp(-(R/AU +7)/10)
def stellar_irradation(R,kp_stellar,Tgas):
    Tau_i = Tau(R, Tgas,kp_stellar)
    F_star = ss*tstar**4 **(rstar/R)*np.exp(-Tau_i)
    return F_star



kp_stellar =  1586.4162432323064
kp_Trim =973.8276480930514
R_rim =  0.5399999947476138 *AU
print( stellar_irradation(R_rim,kp_stellar,T_rim))
Tn =1520
#delta_t = 2.5*10**-15
delta_t = 10**-2
adabatic_index = 5/3
cV = 3/2 *gasR
cP = 5/2 *gasR
Tgas = np.zeros(nr) 
F_irr = np.zeros(nr) 
rad_E = np.zeros(nr) 
fc = np.zeros(nr) 
Qheat = np.zeros(nr)

def Flux_limited_diffustion(R,kp_stellar,kp_Trim):
    Tgas[0] = tstar

    for i,Ri in enumerate(R):
        print(i)
        F_irr[i] = stellar_irradation(Ri,kp_stellar,Tgas) #Eqn 17
        rad_E[i] = a*Tgas[i]**4 + (kp_stellar/kp_Trim)*(abs(F_irr[i]))/clight  #eqn16
        h_cg = pressure_height(Tgas[i], Ri)
        sigma = surface_density(Ri)
        #print(sigma, h_cg)
        Rho_i =Rho(0, h_cg, sigma) #Deriving midplane density
        
        #fc = deriving eqn 9 
        fc[i] = (cV*Rho_i/(4*a*Tgas[i]**3) +1)**(-1)
        
        cs = np.sqrt(kb*Tgas[i]/(mu*mp))
        vt =cs*h_cg*10**-2
        Qheat[i] = Rho_i *vt*9*GG*mstar/(4*Ri**3)
        
        
        kt = Rho_i*vt*cP
        
        
        

        #Eqn 15
        #∇.F = dE/dt
        #∇.F* = -kp(T*)* rho * F* (F* =eqn 17)
        rad_E[i+1]  = rad_E[i] +(delta_t)*(Rho_i*kp_stellar*F_irr[i] -Qheat[i])/(1+fc[i])
        
        #Eqn 16
        Tgas[i+1] = (((kp_stellar/kp_Trim)*(abs(F_irr[i]))/clight +rad_E[i+1])/clight)**0.25
        print(Tgas[i], Tgas[i+1])
    return Tgas

print(Flux_limited_diffustion(R,kp_stellar,kp_Trim))

Tgass= Flux_limited_diffustion(R,kp_stellar,kp_Trim)



fig, ax = plt.subplots()
ax.semilogx(R/AU, Tgas)

ax.set_xlabel("Radius")
ax.set_ylabel("Temp_profile")
ax.set_title("Temp_profile vs. Radius")

plt.show()


fig, ax = plt.subplots()
ax.semilogx(R/AU, Tau(R, T_rim, kp_stellar))

ax.set_xlabel("Radius")
ax.set_ylabel("Surface Density")
ax.set_title("Surface Density vs. Radius")

plt.show()












alpha_guess = 0.2


#convering iteration such that the it iterate until the error is smaller than 0.00001
#THis only check for a radial point
def converging_disk(alpha_guess, r_ii, psi_i, psi_s):
    a_i = 0.22
    while np.round(a_i, 5) != np.round(alpha_guess,5):
    #while alpha_guess !=a_i:
        alpha_guess = a_i
        T_i =irr_disk_temp(alpha_guess, r_ii,psi_i, psi_s)
        h_i= pressure_height(T_i,r_ii)
        s_i = surface_density(r_ii)
        X_i = chi_cg(s_i,kp_stellar)
        H_i = X_i*h_i
        a_i =alpha_angle(H_i, r_ii)
        
    return a_i



X_cg = np.zeros(nr)
impinge_angle = np.zeros(nr)
H_cg = np.zeros(nr)
h_cg = np.zeros(nr)
sigma = np.zeros(nr)
Ti = np.zeros(nr)
   
Metallicity =0.3
T_stellar = tstar
psi_i, psi_s =1, 1
#This iterate through every radial points

def converging_disk2(R, psi_i, psi_s,alpha_guess):
        
    for i, r_i in enumerate(R):
        impinge_angle_iter =converging_disk(alpha_guess, r_i,psi_i, psi_s)
        alpha_guess = impinge_angle_iter
        sigma_iter = surface_density(r_i)
        w = 0.9999999
        T_s = 1200
        h_cg_iter = 0.3
        error = 10
        Ti_iter = irr_disk_temp(alpha_guess, r_i,psi_i, psi_s)
        #This check Ti multiple times such that error is less than 0.01
        while error >0.01:
            Temp = Ti_iter
            Rho_gas =Rho(0, h_cg_iter, sigma_iter)
            P_gas =GasPressure(Rho_gas, Temp)
            kp_Ts = planckOpacity(Metallicity,T_stellar,T_s,P_gas,Rho_gas)
            T_s= irr_surface_Temp(r_i, kp_Ts)

            kp_Ti = planckOpacity(Metallicity,T_stellar,Temp,P_gas,Rho_gas)
            psi_i = self_absorbed_fraction(sigma_iter, kp_Ti)
            psi_s = flux_absorbed_fraction(sigma_iter, kp_Ts)
            #print("Printing , psi value, ",psi_i, psi_s)
            Ti_iter =  w*Temp +(1-w)*irr_disk_temp(alpha_guess, r_i,psi_i, psi_s)
            error = np.max(np.abs(Temp - Ti_iter))
      
            

        
        h_cg_iter = pressure_height(Ti_iter,r_i)
        X_cg_iter = chi_cg(sigma_iter,kp_stellar)
        H_cg_iter = X_cg_iter*h_cg_iter
        #impinge_angle_iter =converging(alpha_guess, r_i)    
        #np.append(Ti, Ti_iter)
        Ti[i] = Ti_iter
        h_cg[i] =h_cg_iter
        sigma[i] =sigma_iter
        X_cg[i] =X_cg_iter
        H_cg[i] =H_cg_iter
        impinge_angle[i] =impinge_angle_iter
    return Ti,h_cg,sigma,X_cg,H_cg,impinge_angle

Ti,h_cg,sigma,X_cg,H_cg,impinge_angle =converging_disk2(R, psi_i, psi_s,alpha_guess)
T_interior = Ti
RimHeating =0
T_effective = (T_visc**4  +Ti**4 +RimHeating**4)**0.25


z_elevation =zi = np.linspace(0.0, H_cg, nz)
#z_elevation = np.array(z_elevation)
#z_elevation   = np.array([0.5 * ( zi[i] + zi[i+1])   for i in range(len(zi) - 1)])  # Take the in-"between values       # Take the in-between values    
Z = np.linspace(0.0, 50*AU, nz)




Rho_rim = Rho(0, h_rim, sigma_rim)
P_rim =GasPressure(Rho_rim, T_rim)


Rho_gas = Rho(z_elevation, h_cg, sigma)
Gas_pressure =GasPressure(Rho_gas, T_effective)
P_gas = GasPressure(Rho_gas, T_effective)




Rflare = np.where (R <H_cg*R_rim/H_rim  )#and H_cg*R_rim/R/H_rim >1)
Rshadow = np.where (R >=H_cg*R_rim/H_rim)



Ti[Rshadow] = 0
T_effective = (T_visc**4  +Ti**4 +RimHeating**4)**0.25
h_cg = pressure_height(T_effective, R)
X_cg = chi_cg(sigma,kp_stellar)
H_cg = scale_height(X_cg, h_cg)

#h_cg[0] = h_rim 
#H_cg[0] = H_rim
#X_cg[0] = X_rim


z_elevation =zi = np.linspace(0.0, H_cg, nz)
#z_elevation = np.array(z_elevation)
#z_elevation   = np.array([0.5 * ( zi[i] + zi[i+1])   for i in range(len(zi) - 1)])  # Take the in-"between values       # Take the in-between values    
Z = np.linspace(0.0, 50*AU, nz)




Rho_rim = Rho(0, h_rim, sigma_rim)
P_rim =GasPressure(Rho_rim, T_rim)


Rho_gas = Rho(z_elevation, h_cg, sigma)
Gas_pressure =GasPressure(Rho_gas, T_effective)
P_gas = GasPressure(Rho_gas, T_effective)






def Stability(h_cg,sigma,R,Ti):
    cs = np.sqrt(kb*Ti/(mu*mp))
    kp = np.sqrt(GG*mstar/R**3)
    return cs*kp/(np.pi * GG*sigma)

stability_cretrion = Stability(h_cg,sigma,R,T_effective)


SD = 9.03*10**4*(R/AU)**(-1.75)
To = (300882 * (R/AU)**(-2*0.5) + 10**2)**0.5
cs = np.sqrt(kb*To/(mu*mp))
kp = np.sqrt(GG*0.5*ms/R**3)
Q_stab=cs*kp/(3.14*GG*SD)


Q_stability = (np.sqrt(GG*mstar/R**3)*np.sqrt(kb*T_effective/mu/mp))/(np.pi*GG*sigma)
Q_stability_visc = (np.sqrt(GG*mstar/R**3)*np.sqrt(kb*T_visc/mu/mp))/(np.pi*GG*sigma)

    

fig, ax = plt.subplots()
ax.semilogx(R/AU, Ti, marker = '*', markevery = (Q_stability <1))
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
ax.semilogx(R/AU, H_cg/R,marker = '*', markevery = (Q_stability <1))
ax.axhline(y=H_rim/R_rim, xmin=0, xmax=100, color = 'b', linestyle = '--')

ax.set_xlabel("Radius")
ax.set_ylabel("Surface_Height")
ax.set_title("Surface_Height vs. Radius")
#ax.set_ylim([0, 0.5])

plt.show()

fig, ax = plt.subplots()
ax.semilogx(R/AU, h_cg/R,marker = '*', markevery = (Q_stability <1))

ax.set_xlabel("Radius")
ax.set_ylabel("Pressure_Height")
ax.set_title("pressure_Height vs. Radius")
ax.set_ylim([0, 0.5])

plt.show()



fig, ax = plt.subplots()
ax.semilogx(R/AU, surface_density(R),marker = '*', markevery = (Q_stability <1))

ax.set_xlabel("Radius")
ax.set_ylabel("Surface Density")
ax.set_title("Surface Density vs. Radius")

plt.show()





fig, ax = plt.subplots()
ax.semilogx(R/AU, X_cg)


ax.set_xlabel("Radius")
ax.set_ylabel("X_cg")
ax.set_title("x_cg vs. Radius")
plt.show()   

 


fig, ax = plt.subplots()
ax.plot(R/AU, H_cg*R_rim/R/H_rim,marker = '*', markevery = (Q_stability <1))
ax.set_xlabel("R_cg")
ax.set_ylabel("Flaring_point")
ax.set_title("Flaring_point vs. Radius")
plt.show()   
 
R_inst = (3*10**(-5)  *(mstar/ms)**(13/3) *(lstar/ls)**(-4) * (H_cg/h_cg)**(-4)* (sigma**(14/3))) *AU  



fig, ax = plt.subplots()
ax.semilogx(R/AU, Q_stability)

ax.set_xlabel("Radius")
ax.set_ylabel("Toomre criterion")
ax.set_title("Ioomre criterion vs. Radius")
plt.show()   


fig, ax = plt.subplots()
ax.loglog(R/AU, Ti, marker = '*', markevery = (Q_stability <1))
ax.set_xlabel("Radius")
ax.set_ylabel("Interior_Temperature")
ax.set_title("Interior_Temperature vs. Radius")
plt.show()   


Rim_tempeature = (h_rim/R_rim)**2*(GG*mstar/R_rim)*(mu/gasR)
print("Rim temperature", Rim_tempeature)


inner_temp = (2/3/np.pi)**0.25 *(rstar/R_rim)**.75 *tstar