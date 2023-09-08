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
tstar = 9000
lstar = 47* ls

#Disk parameter
T_rim= t_sublimation = 1500.0  # Temperature at which dust start to condense
sigma0 = 2*10**3
sigmagas0 =2
psi_i, psi_s = 1, 1 # Absorption and emmission factor 

#Inner gas hole parameter
Metallicity =0.3
T_stellar =tstar


nu = 5.879 * (10 ** 10) * tstar  #Wien law
t_virial = GG * mstar * mu * mp / (kb * rstar)  # Virial Temp
gamma = 2.0  #Accounting circumstellar rim 
beta = -1.5  #density decay factor 
gray_const =1
visc_const = 10**-2 #Shakura & Sunyaev turblent viscosity constant

#Astro-silicate paramater
A_sil = 10**(-25.10)   #A_sil = 10^-25.10 cm^2.5
a_min = 5*10**-7
a_max = 2.5 * 10**-5
sil_density = 3.4 # g cm^-3
#dust to gas mass ratios ζgrap = 0.0025 and ζsil = 0.0034
#(Draine & Lee 1984). 
num_sil_density = sil_density/mp *A_sil*(a_max**-2.5 - a_min**-2.5)/(-3.5)
kv = (3/4)*(1.5*10**4)*(10**(-25.10)) * ((a_max**-2.5 - a_min**-2.5)/(-3.5))
ni_const = (3/4)*(A_sil)*((a_max**-2.5 - a_min**-2.5)/(-3.5)) /(10**-6)
fd2g = 10**-2    # The gas: dust is assumed be 100 : 1.



# Set maximum iteration
maxIter = 500





################################*Defning different paramater *#################

#Surface density profile
def surface_density(R):
    return sigma0*(R/AU)**beta

def pressure_height(T, R):
    return ((T/t_virial)**0.5 * (R/rstar)**0.5 * R)

def scale_height(X_cg, h_cg):
    return (X_cg *h_cg)

def chi_cg(sigma,kp_stellar):
    return (erfinv(1-1/(4*sigma*kp_stellar)))
    #return 1-  (4*sigmaRim*kp_stellar)*(1-math.erf(X_rim_iter))


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

#Defining interior heating temperature
def irr_disk_temp(impinge_angle, R,psi_i, psi_s):
    return ((impinge_angle*psi_s/psi_i)**0.25 * (rstar/R)**0.5 *tstar)

#Definging surface tempeature
def irr_surface_Temp(R, kp_Ts):
    epsilon =kp_Ts/kp_stellar
    T_s = (1/epsilon)**0.25 * (rstar/2/R)**0.5 *tstar
    return T_s


#Direct stellar radiation impinges onto the disk at an angle
def alpha_angle(H_cg,R):
    return (0.4*rstar/R   +(gamma - 1)*H_cg/R)



def density_distribution(z_elevation, h_cg, sigma):
    return sigma*np.exp(- z_elevation**2 / 2/h_cg**2) / (np.sqrt(2*np.pi) *h_cg)

#Pressure [dyn/cm^2] = [g cm/s^2/cm^2] = [g/s^2/cm]
def GasPressure(rho, Tgas):
    return (kb*Tgas*rho/mu/mp)


from scipy.integrate import quad
def intergrand(R, Tgas,kp_T):
    return kp_T*(sigma0*(R/AU)**beta)/np.sqrt(2*np.pi)/((Tgas/t_virial)**0.5 * (R/rstar)**0.5 * R)

def Tau(R, Tgas,kp_T):
    return quad(intergrand,np.inf,R,args=(Tgas,kp_T))[0]

def closest(lst, K):
     lst = np.asarray(lst)
     idx = (np.abs(lst - K)).argmin()
     return lst[idx]

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



#############################*Modelling Kp_stellar*############################

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



#################################*Modelling Inner hole*########################

from scipy.integrate import quad
def intergrand(R, Tgas,kp_T,sigma0):
    return kp_T*(sigma0*(R/AU)**beta)/np.sqrt(2*np.pi)/((Tgas/t_virial)**0.5 * (R/rstar)**0.5 * R)

def Tau(R, Tgas,kp_T,sigma0):
    return quad(intergrand,R,(rstar,R),args=(Tgas,kp_T))[0]

#Tau[R*, R] = integral of rho as a function of (sigma/np.sqrt(2pih**2))
#where h = ((Tgas/t_virial)**0.5 * (R/rstar)**0.5 * R)
#where sigma= sigma0 *(R/AU)**-1.5
def optical_depth(Tgas,R,kp_T,sigma0):
    return -0.25*np.sqrt(2)*kp_T*sigma0/(np.sqrt(np.pi)*(R/AU)**1.5*(R/rstar)**0.5*(Tgas/t_virial)**0.5) + 0.25*np.sqrt(2)*kp_T*sigma0/(np.sqrt(np.pi)*(rstar/AU)**1.5*(Tgas/t_virial)**0.5)




#where sigma0 = 1000gcm^-2
def gas_density(R):
    return sigmagas0*(R/AU)**(-1.5)


def innerHole(Ri):
    Tguess = T_rim
    sigma0 =2

    Tgas = 9200
    #Load original data
    data = np.loadtxt('planck_data.csv', delimiter=',', skiprows=1)
    #Load without kp(T) column
    modified_data = np.loadtxt('planck_data.csv', delimiter=',', skiprows=1,  usecols=(0,1, 2,3,4))
    #Load without Kp(T) and T*
    #THis is to find kp(Tgas)
    gas_data = np.loadtxt('Planck_gas.csv', delimiter=',', skiprows=1)
    modified_gas_data = np.loadtxt('Planck_gas.csv', delimiter=',', skiprows=1, usecols=(0,1, 2,3))

    #print(Tgas, Tguess)
    iteration = 0
    while np.round(Tguess, 2) != np.round(Tgas,2):
        Tguess = Tgas
        s_gas = gas_density(Ri)
        #s_gas =surface_density(Ri) #Suurface density profile
        h_gas =pressure_height(Tgas, Ri)# gas pressure scale height
        
        Rho_gas = density_distribution(0, h_gas, s_gas) # Density profile in midplane z = 0
        Rho_gas = 3.16e-14
        P_gas = GasPressure(Rho_gas, Tgas)
        #print("printing pressute",P_gas,Rho_gas)
        # Vector we try to match with data
                #Defining Kp(Tgas, Trad)
        
        search_vect =  np.array([Metallicity, T_stellar, Tgas ,P_gas,Rho_gas])
        #search_vect =  np.array([-0.3,	3000,	700,	1.00E-09,	3.99E-20])

        # normalize data in each column with the column mean
        norm_data = modified_data /np.mean(modified_data, axis=0)
        #print("printing normalised vecotr",norm_data)
        # normalize the search vector with the same values
        s_vec_n = search_vect/np.mean(modified_data , axis=0)
        idx = np.argmin(np.linalg.norm((s_vec_n - norm_data), axis=1))
        kp_TgasTrad = data[idx][5]
        #print(idx)
        
        #Defining Kp(Tgas, Tgas)
        #I have removed Tstellar column as stellar temperature only consist
        # of T>3000K so I thought there would be any convergence since
        #it cannot find closest value 
        
        search_vect2 =  np.array([Metallicity, Tgas ,P_gas,Rho_gas])
        # normalize data in each column with the column mean
        norm_data2 = modified_gas_data /np.mean(modified_gas_data, axis=0)
        #print("printing normalised vecotr",norm_data)
        # normalize the search vector with the same values
        s_vec_n2 = search_vect2/np.mean(modified_gas_data , axis=0)
        idx2 = np.argmin(np.linalg.norm((s_vec_n2 - norm_data2), axis=1))
        kp_Tgas = gas_data[idx2][5]
        
        
        #This is method 2: where it find closest of each row
        #and take 5th row as kp(T) values
        kp_TgasTrad2 =planckOpacity(Metallicity,T_stellar,Tgas,P_gas,Rho_gas)
        kp_Tgas2 =planckOpacity(Metallicity,Tgas,Tgas,P_gas,Rho_gas)
        eps2 = kp_Tgas2/kp_TgasTrad2
        #print(idx, idx2)
        #print(kp_TgasTrad2,kp_Tgas2)
        #Tgas2 =(0.5*(1+3*Tauuuuu/2))**0.5*tstar
        
        #Correction factor
        #Close to star we cannt use R>>R*
        C = 2*(1-np.sqrt(1-(rstar/Ri)**2))*(Ri/rstar)**2 
        eps = kp_Tgas/kp_TgasTrad
        #Tgas = (C/eps)**(0.25)* (rstar/Ri/2)**0.5 *tstar
        Tgas = (C/eps)**(0.25)* (rstar/Ri/2)**0.5 *tstar
        #print(kp_TgasTrad, kp_Tgas)
        #print("\n",Tgas)
        
        
        #The solution doesnt converge
        iteration +=1
        if iteration == 50:
            break
    return kp_Tgas,Rho_gas, Tgas

def innerHole_Temp(R_hole):
    for i, ri in enumerate(R_hole):
        Temp_innerhole[i] = innerHole(ri)[2]
    return Temp_innerhole



################################*Modelling Puffed up rim*#####################


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
            kp_TgasTrad,Rho_gas, Tgas = innerHole(R_rim)
            opticalll =  kp_TgasTrad*Rho_gas*(R_rim - rstar)
            gray_const =(1-np.exp(-opticalll))
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
print(R_rim/AU,H_rim/R_rim)   


#calculating rim temperature
R_hole = np.linspace(rstar,R_rim,nr)
#Rim_temp = innerHole(R_rim)
#print(innerHole(R_rim))
#Temp_innerhole = np.zeros(nr)
sigma_gas = gas_density(R_hole)



#calculating rim temperature
Metallicity = -0.3
R_rim, H_rim, h_rim, sigma_rim, X_rim = converging_rim(gray_const)
print(R_rim/AU,H_rim/R_rim,Metallicity) 
#0.5399999942393016 0.14096094205659093 -0.3
#H_rim =0.07611890789852563
Metallicity = 0
R_rim, H_rim, h_rim, sigma_rim, X_rim = converging_rim(gray_const)
print(R_rim/AU,H_rim/R_rim,Metallicity) 
#0.5399999943692885 0.1409609420735568 0
#H_rim =0.07611890792601027

Metallicity = 0.3
R_rim, H_rim, h_rim, sigma_rim, X_rim = converging_rim(gray_const)
print(R_rim/AU,H_rim/R_rim,Metallicity) 
#0.5399999944663407 0.140960942086224 0.3
#H_rim =0.07611890794653114

