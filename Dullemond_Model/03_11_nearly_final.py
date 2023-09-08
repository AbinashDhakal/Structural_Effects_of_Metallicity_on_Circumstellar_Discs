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
cV = kb/(mu*mp*(adia_gamma-1)) #specific const volume 
cP = adia_gamma *cV #specific const pressure

# grid parameters 
nr = 256  # Number of radial grids. 
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


psi_i, psi_s = 1, 1



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

def irr_disk_temp(impinge_angle, R,psi_i, psi_s):
    return ((impinge_angle*psi_s/psi_i)**0.25 * (rstar/R)**0.5 *tstar)


def irr_surface_Temp(R, kp_Ts):
    epsilon =kp_Ts/kp_stellar
    T_s = (1/epsilon)**0.25 * (rstar/2/R)**0.5 *tstar
    return T_s

def alpha_angle(H_cg,R):
    return (0.4*rstar/R   +(gamma - 1)*H_cg/R)



#Q_visc = Rho* vt*(r*dr_omega)**2
#Q_visc = (sigma/(2*pi*h))*(visc_const)*(np.sqrt(kb*T/mu/mh))*(9*G*Mstar/4/R**3)
def Q_visc(R,T):
    sigma = surface_density(R)
    cs = np.sqrt(kb*T/mu/mp)
    return (sigma/2/np.pi)*visc_const*cs*(GG*mstar/4/R**3)

#Q_cond = kt∇^2T
#kt : thermal conductivity
def Q_cond(R,T,dr):
    #Defining a potential form of T
    
    T = (0.7)*(rstar/R)**0.5 *tstar
    Q = np.zeros(nr)
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
        vt = visc_const*cs**2/omega #viscosity
        sigma = surface_density(R[i]) #surface density 
        rho = sigma/(np.sqrt(2*np.pi* hp**2)) #Density 
        kt = rho*vt*adia_gamma*cV #Turblent thermal conductivity
        #Discretising laplace equation

        dT_dR_2 = (T[i+1] -2*T[i] +T[i-1])/(dr)**2
        dT_dR = (T[i+1] -T[i])/dr/2
        Q[i] =kt*(dT_dR_2 +dT_dR/R[i])
        #Q_cool = 2 *ss* T
        #Q_condct = kt * ΔT
        #Q_cool = Q_conduct
    return Q
    

alpha_guess = 0.2
#convering iteration such that the it iterate until the error is smaller than 0.00001
#THis only check for a radial point
def converging_alpha(alpha_guess, r_ii, psi_i, psi_s):
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
    
    
 #alpha*lstar/(4*np.pi*R**2)
def Q_irr(R,alpha_guess):
    T = (0.5)*(rstar/R)**0.5 *tstar
    for i in range(1,nr-1):
        alpha = converging_alpha(alpha_guess, R[i], psi_i, psi_s)
        T[i] = ((alpha*lstar/(4*np.pi*R[i]**2))/(2*ss))**0.25#alpha*lst00000ar/(4*np.pi*R**2)
    return T
    
R_rim =  0.5399999947476138 *AU
H_rim = 0.0593989*AU
rin = R_rim # Sublimation radius  # inner rim of disk. Density is assumed to
rout = 100*AU
R =ri = np.linspace(rin,rout,nr) 
#R   = 0.5 * ( ri[0:nr] + ri[1:nr+1] )             # Take the in-"between values
R = np.array(R)

print(Q_irr(R,alpha_guess))

kp_stellar =  1586.4162432323064
kp_Trim =973.8276480930514


# Set maximum iteration
maxIter = 500

#THis is collecting everything and converging and outputing T[i]
#Q_heat = Q_cool 
#Q_heat = Q_visc + Q_irr+ Q_cond
def converging(R,alpha_guess):
    T = (0.7)*(rstar/R)**0.5 *tstar

    Q_heat= np.zeros(nr)
    dr = (R[1]-R[0])

    for iteration in range(0, maxIter):
        for i in range(1,len(T)):
            alpha = converging_alpha(alpha_guess, R[i], psi_i, psi_s)
            alpha_guess =alpha
            Q_heat[i] = Q_irr(R[i],T[i]) + Q_visc(R[i],T[i]) + Q_cond(R[i],T[i],dr)
            #Q_cool = 2 ss* T[i]**4
            T[i] = (Q_heat/2/ss)**0.25 
    return T
            

T_effective = converging(R,alpha_guess)
h_cg = pressure_height(T_effective, R)
sigma = surface_density(R)
X_cg = chi_cg(sigma,kp_stellar)
H_cg = scale_height(X_cg, h_cg)
alpha = alpha_angle(H_cg, R)





#########################Modelling shadowed region
#THis should be defububg what happen in the region its in shaded region



h_cg[0] = h_rim  = 0.01932884992839139*AU
H_cg[0] = H_rim =0.074*AU
X_cg[0] = X_rim = 4.294689941474017

Rshadow = np.where (R >=H_cg*R_rim/H_rim)
Rflare = np.where(R>10*AU)


Q_irr(alpha,R)[Rflare]= 0   #This should set Q_irr =0 in shadowed region
T_effective = converging(R)
h_cg = pressure_height(T_effective, R)
sigma = surface_density(R)
X_cg = chi_cg(sigma,kp_stellar)
H_cg = scale_height(X_cg, h_cg)
impinge_angle = alpha_angle(H_cg, R)



h_cg[0] = h_rim  = 0.01932884992839139*AU
H_cg[0] = H_rim =0.074*AU
X_cg[0] = X_rim = 4.294689941474017









#######################Graphing#######################################################




def Stability(h_cg,sigma,R,Ti):
    cs = np.sqrt(kb*Ti/(mu*mp))
    kp = np.sqrt(GG*mstar/R**3)
    return cs*kp/(np.pi * GG*sigma)



    
Q_stability = (np.sqrt(GG*mstar/R**3)*np.sqrt(kb*T_effective/mu/mp))/(np.pi*GG*sigma)

fig, ax = plt.subplots()
ax.semilogx(R/AU, T_effective, marker = '*', markevery = (Q_stability <1))
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



