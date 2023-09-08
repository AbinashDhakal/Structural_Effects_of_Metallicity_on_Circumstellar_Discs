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
kp_stellar = 400.0  # Monochromatic opacity = 400cm^2 g^-1

#kP(T*) = 2100 cm2 g-1;   https://iopscience.iop.org/article/10.3847/1538-4357/835/2/230/pdf 
#Kp(T_rim) = 700 cm2 g-1

nu = 5.879 * (10 ** 10) * tstar  #Wien law
t_virial = GG * mstar * mu * mp / (kb * rstar)  # Virial Temp
gamma = 2.0
beta = -1.5
gray_const =1


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
"""
def converging_rim(gray_const):
    R_guess = 0.53*AU
    R_rim = 0.54*AU
    H_rim =0.1162*AU
    iteration =0
    while np.round(R_guess/AU, 8) != np.round(R_rim/AU,8):
        R_guess =        R_rim
        w = 0.9999999
        error = 1*AU
        while error >10e-9 *AU:
            temp_r = R_rim
            R_rim = (w*R_rim +(1-w)*Rim_radial(H_rim, R_rim,gray_const))
            error = np.max(np.abs(temp_r-R_rim))
            print(R_rim/AU, temp_r/AU)
        iteration +=1
        print("I am running iteration", iteration)
        sigma_rim = surface_density(R_rim)
        h_rim = ((T_rim/t_virial)**0.5 * (R_rim/rstar)**0.5 * R_rim)
        X_rim = chi_cg(sigma_rim,kp_stellar)
        H_rim = scale_height(X_rim, h_rim)

    return R_rim, H_rim, h_rim, sigma_rim, X_rim

"""

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
        #print(H_rim/R_rim, H_guess/R_rim)

    return R_rim, H_rim, h_rim, sigma_rim, X_rim, 




R_rim, H_rim, h_rim, sigma_rim, X_rim = converging_rim(gray_const)
#print(R_rim/AU,H_rim/R_rim)   


#Definig grid system in radial direction
rin = R_rim # Sublimation radius  # inner rim of disk. Density is assumed to
rout = 100*AU
R =ri = np.linspace(rin,rout,nr) 
#R   = 0.5 * ( ri[0:nr] + ri[1:nr+1] )             # Take the in-"between values
R = np.array(R)



def viscous_temp(R):
    yrs = 31557600 #s
    M_acc = 10**-8*ms/yrs
    return (3*GG*mstar*M_acc/(8*np.pi* ss* R**3))**0.25
    #return 107*(R/AU)**-0.75 * (mstar/2.5/ms)**0.25 *(M_acc/(10**-8 *ms/yrs))**0.25

T_visc = viscous_temp(R)

betaaa = 1
def Rim_heating(R):
    #return 1.2*T_rim*np.exp(-0.12944943670387588  *(R-R_rim))
    
    
    #Rimheating = np.zeros(nr)
    #for i,(ri) in enumerate(R):
    #    Rimheating[i] =1.2*T_rim*np.exp((-1)*(1-0.5**(1/(4+betaaa)))  *(ri/AU-R_rim/AU))
    return T_rim*np.exp(-1*(1-0.5**(1/(4+betaaa)))  *(R-R_rim)/(5834903615946.781))
    #return Rimheating
RimHeating = Rim_heating(R)

def surface_temp(R):
    return (0.25)**(1/(4+beta)) * (rstar/R)**(2/(4+beta)) * tstar
Ts = surface_temp(R)

def alpha_angle(H_cg,R):
    return (0.4*rstar/R   +(gamma - 1)*H_cg/R)

psi_i, psi_s = 1, 1
def irr_disk_temp(impinge_angle, R):
    return ((impinge_angle*psi_s/psi_i)**0.25 * (rstar/R)**0.5 *tstar)


alpha_guess = 0.2
#convering iteration such that the it iterate until the error is smaller than 0.00001
#THis only check for a radial point
def converging_disk(alpha_guess, R):
    a_i = 0.22
    for i,r_ii in enumerate(R):
       while np.round(a_i, 5) != np.round(alpha_guess,5):
        #while alpha_guess !=a_i:
            alpha_guess = a_i
            T_i =irr_disk_temp(alpha_guess, r_ii)
            T_e = (Rim_heating(r_ii)**4 + viscous_temp(r_ii)**4 +T_i**4 )**0.25
            h_i= pressure_height(T_e,r_ii)
            s_i = surface_density(r_ii)
            X_i = chi_cg( s_i,kp_stellar)
            H_i = X_i*h_i
            a_i =alpha_angle(H_i, r_ii)
        

            Ti[i] = T_i
            T_effective[i] = T_e
            h_cg[i] =h_i
            sigma[i] =s_i
            X_cg[i] =X_i
            H_cg[i] =H_i
            impinge_angle[i] =a_i
    return X_cg,impinge_angle,H_cg,h_cg,sigma,Ti,T_effective

X_cg = np.zeros(nr)
impinge_angle = np.zeros(nr)
H_cg = np.zeros(nr)
h_cg = np.zeros(nr)
sigma = np.zeros(nr)
Ti = np.zeros(nr)
T_effective = np.zeros(nr)
X_cg,impinge_angle,H_cg,h_cg,sigma,Ti,T_effective =converging_disk(alpha_guess, R)
   




#This iterate through every radial points
"""
for i, r_i in enumerate(R):
    impinge_angle_iter =converging_disk(alpha_guess, r_i)
    alpha_guess = impinge_angle_iter
    sigma_iter = surface_density(r_i)
    Ti_iter = irr_disk_temp(alpha_guess, r_i)
    Teff_iter = (Rim_heating(r_i)**4 + viscous_temp(r_i)**4+Ti_iter**4 )**0.25
    h_cg_iter = pressure_height(Teff_iter,r_i)
    X_cg_iter= chi_cg(sigma_iter, kp_stellar)
    H_cg_iter = X_cg_iter*h_cg_iter
    #impinge_angle_iter =converging(alpha_guess, r_i)    
    #np.append(Ti, Ti_iter)
    Ti[i] = Ti_iter
    T_effective[i] = Teff_iter
    h_cg[i] =h_cg_iter
    sigma[i] =sigma_iter
    X_cg[i] =X_cg_iter
    H_cg[i] =H_cg_iter
    impinge_angle[i] =impinge_angle_iter

#T_effective = (RimHeating**4 + T_visc**4 +Ts**4  +Ti**4)**0.25
"""

T_effective = (T_visc**4  +Ti**4 +RimHeating**4)**0.25


"""
R_rim = 0.47*AU
H_rim - 0.11*R_rim
"""

Rflare = np.where (R <H_cg*R_rim/H_rim  )#and H_cg*R_rim/R/H_rim >1)
Rshadow = np.where (R >=H_cg*R_rim/H_rim)



Ti[Rshadow] = 0
T_effective = (T_visc**4  +Ti**4 +RimHeating**4)**0.25
h_cg = pressure_height(T_effective, R)
X_cg = chi_cg(sigma,kp_stellar)
H_cg = scale_height(X_cg, h_cg)

h_cg[0] = h_rim 
H_cg[0] = H_rim
X_cg[0] = X_rim


z_elevation =zi = np.linspace(0.0, H_cg, nz)
#z_elevation = np.array(z_elevation)
#z_elevation   = np.array([0.5 * ( zi[i] + zi[i+1])   for i in range(len(zi) - 1)])  # Take the in-"between values       # Take the in-between values    
Z = np.linspace(0.0, 50*AU, nz)


def Rho(z_elevation, h_cg, sigma):
    return sigma*np.exp(- z_elevation**2 / 2/h_cg**2) / (np.sqrt(2*np.pi) *h_cg)

Rho_gas = Rho(z_elevation, h_cg, sigma)
#Pressure [dyn/cm^2] = [g cm/s^2/cm^2] = [g/s^2/cm]

def GasPressure(rho, T_disk):
    return (kb*T_disk*rho/mu/mp)

#Force due to gravity
#g(z) = GMz/R^3

#Gas pressure distribution
#dP/dz = -rho F_gravity
# dP [Pc,Pz] = (mu*G*M/(RR*T R**3)) ,RR is ideal gas consntat RR = 8.3145\,  J \cdot mol^{-1}\cdot K^{-1}
RR = 8.31446261815324*10**7	#erg⋅K−1⋅mol−1
#|Pz = Pc* exp(-mu GM z^2/(2RR T R^3))

#Surface density: Σ = Integral of ( rho dz) between [-inf, inf]
#surface density: Σ = Σ0(R/AU)^-1.5

#Σ = integral[-inf, inf] of (mu*Pc/(RR T)) exp(-mu GG M/(2RR T R^3)) dz
#Σ/ Pc = (2 pi mu R^3/(RR TGM)) , Pc = central pressure
def central_pressure(sigma, R, Tdisk):
    return (2* np.pi* mu* R**3/(RR*Tdisk*GG*mstar))**-0.5 * sigma
#Pc = central_pressure(sigma, R, T_effective)

    
def pressure_distribution(z_elevation, Ti, R, Pc):
    return (Pc* np.exp(-mu*GG*mstar*z_elevation**2/(2*RR*Ti/R**3)))
#Pz = pressure_distribution(z_elevation, Ti, R, Pc)





Rho_rim = Rho(0, h_rim, sigma_rim)
Pc_rim = central_pressure(sigma_rim, Rho_rim, T_rim)
P_rim =GasPressure(Rho_rim, T_rim)

Metallicity = -0.3
T_stellar = tstar
print(Metallicity, T_stellar)
Gas_pressure =GasPressure(Rho_gas, T_effective)
P_gas = GasPressure(Rho_gas, T_effective)


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
    
    mean_stellar_opacity =  df[
        (df["T_stellar"] == (closest(df["T_stellar"], T_stellar))) 
        & (df["Metallicity"] == (closest(df["Metallicity"], Metallicity))) ]
    kp_stellar = mean_stellar_opacity["Planck_opacity"].mean()

    print(kp_stellar, planck_opacity_data[0][5])
    return kp_stellar, planck_opacity_data[0][5]
kp_stellar, Planck_rim =  planckOpacity(Metallicity,T_stellar,T_rim,P_rim,Rho_rim)

kp_stellar, Planck_rim = kp_stellar/10, Planck_rim *100
print(R_rim/AU, H_rim/R_rim)   
gray_const =( 2* Planck_rim/kp_stellar + 1 + 2*kp_stellar/Planck_rim * np.exp(-1))



#gray_const = 4
R_rimm, H_rimm, h_rimm, sigma_rimm, X_rimm = converging_rim(gray_const)
print(R_rimm/AU, H_rimm/R_rimm)   






kg =10**-4
planck_Td = 700
kp_stellar = 2100

def dust_gas_ratio(T_dust,planck_Td,kp_stellar,R):
    const_top = kg*(4*R**2 *T_dust**4 - R**2*tstar**4)
    const_bottom = kp_stellar*R**2*tstar**4 - 4*planck_Td*R**2*T_dust**4
    return const_top/const_bottom

fd2g = dust_gas_ratio(T_rim,planck_Td,kp_stellar,R_rim)

def epsilon(planck_Td, kp_stellar, fd2g):
    const_top = kg +fd2g*planck_Td
    const_bottom = kg + fd2g*kp_stellar
    return const_top/const_bottom

#e = epsilon(planck_Td, kp_stellar, fd2g)


















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
ax.semilogx(R/AU, Ts)

ax.set_xlabel("Radius")
ax.set_ylabel("Surface_Temperature")
ax.set_title("surface_Temperature vs. Radius")
plt.show()   
    



fig, ax = plt.subplots()
ax.semilogx(R/AU, T_visc)

ax.set_xlabel("Radius")
ax.set_ylabel("Viscous_Temperature")
ax.set_title("Viscous_Temperature vs. Radius")
plt.show()   
    

fig, ax = plt.subplots()
ax.semilogx(R/AU, RimHeating)

ax.set_xlabel("Radius")
ax.set_ylabel("RimHeating")
ax.set_title("RimHeating vs. Radius")

plt.show()




fig, ax = plt.subplots()
ax.semilogx(R/AU, T_effective)

ax.set_xlabel("Radius")
ax.set_ylabel("Effective_temperature")
ax.set_title("Effective_temperature vs. Radius")

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
 
    
