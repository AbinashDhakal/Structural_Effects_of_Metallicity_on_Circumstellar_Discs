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
rstar = 2.5 * rs
tstar = 9000
lstar = 47* ls

#Disk parameter
T_rim= t_sublimation = 1500.0  # Temperature at which dust start to condense
sigma0 = 2*10 **3  # Surface density at 1 AU  Σ0 [g/cm^2] 
kp_stellar = 400.0  # Monochromatic opacity = 400cm^2 g^-1

sigma0gas = 1000 #https://www.ita.uni-heidelberg.de/~dullemond/lectures/radtrans_2012/Chapter_8.pdf
sigma0=2

#kP(T*) = 2100 cm2 g-1;   https://iopscience.iop.org/article/10.3847/1538-4357/835/2/230/pdf 
#Kp(T_rim) = 700 cm2 g-1
#kp_stellar = 0.4
nu = 5.879 * (10 ** 10) * tstar  #Wien law
t_virial = GG * mstar * mu * mp / (kb * rstar)  # Virial Temp
gamma = 2.0
beta = -1.5
gray_const =1


def surface_density(R):
    return sigma0*(R/AU)**beta


def pressure_height(T, R):
    return ((T/t_virial)**0.5 * (R/rstar)**0.5 * R)

def Rho(z_elevation, h_cg, sigma):
    return sigma*np.exp(- z_elevation**2 / 2/h_cg**2) / (np.sqrt(2*np.pi) *h_cg)

#Pressure [dyn/cm^2] = [g cm/s^2/cm^2] = [g/s^2/cm]
def GasPressure(rho, Tgas):
    return (kb*Tgas*rho/mu/mp)

R_rim = 0.5399999944538275*AU
H_rim = 0.07581315724362893*AU
Metallicity =0.3
T_stellar =tstar
T_rim = 1500
R =np.linspace(rstar,R_rim,nr) 
#R   = 0.5 * ( ri[0:nr] + ri[1:nr+1] )             # Take the in-"between values
R = np.array(R)


from scipy.integrate import quad
def intergrand(R, Tgas,kp_T):
    return kp_T*(sigma0*(R/AU)**beta)/np.sqrt(2*np.pi)/((Tgas/t_virial)**0.5 * (R/rstar)**0.5 * R)

def Tau(R, Tgas,kp_T):
    return quad(intergrand,R,(rstar,R),args=(Tgas,kp_T))[0]

#Tau[R*, R] = integral of rho as a function of (sigma/np.sqrt(2pih**2))
#where h = ((Tgas/t_virial)**0.5 * (R/rstar)**0.5 * R)
#where sigma= sigma0 *(R/AU)**-1.5
def optical_depth(Tgas,R,kp_T):
    return -0.25*np.sqrt(2)*kp_T*sigma0/(np.sqrt(np.pi)*(R/AU)**1.5*(R/rstar)**0.5*(Tgas/t_virial)**0.5) + 0.25*np.sqrt(2)*kp_T*sigma0/(np.sqrt(np.pi)*(rstar/AU)**1.5*(Tgas/t_virial)**0.5)


#This function find the closest value of my input from the given table
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

#where sigma0 = 1000gcm^-2
def gas_density(R):
    return sigma0*(R/AU)**(-1.5)


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
        
        Rho_gas = Rho(0, h_gas, s_gas) # Density profile in midplane z = 0
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



#calculating rim temperature
Metallicity = 0
R_hole = np.linspace(rstar,R_rim,nr)
kp_TgasTrad,Rho_gas, Tgas = innerHole(R_rim)
print(innerHole(R_rim))
opticalll =  kp_TgasTrad*Rho_gas*(R_rim - rstar)
print(np.exp(-opticalll))
sublimation_temp= ((1-np.exp(-opticalll))*ss*tstar**4*(rstar/R_rim)**2/ss)**0.25
print(sublimation_temp,Metallicity )


Metallicity = 0.3

kp_TgasTrad,Rho_gas, Tgas = innerHole(R_rim)
print(innerHole(R_rim))
opticalll =  kp_TgasTrad*Rho_gas*(R_rim - rstar)
print(np.exp(-opticalll))
sublimation_temp= ((1-np.exp(-opticalll))*ss*tstar**4*(rstar/R_rim)**2/ss)**0.25
print(sublimation_temp,Metallicity )

Metallicity = -0.3
kp_TgasTrad,Rho_gas, Tgas = innerHole(R_rim)
print(innerHole(R_rim))
opticalll =  kp_TgasTrad*Rho_gas*(R_rim - rstar)
print(np.exp(-opticalll))
sublimation_temp= ((1-np.exp(-opticalll))*ss*tstar**4*(rstar/R_rim)**2/ss)**0.25
print(sublimation_temp,Metallicity )





Temp_innerhole = np.zeros(nr)
sigma_gas = gas_density(R_hole)

#Graphing inner hole temp
fig, ax = plt.subplots()
ax.semilogx(R_hole/AU,sigma_gas, label= Metallicity)
ax.set_xlabel("Radius")
ax.set_ylabel("sigma_gas")
ax.set_title("sigma_gas vs. Radius")
plt.show()   


Temp_innerhole = innerHole_Temp(R_hole)

#Graphing inner hole temp
fig, ax = plt.subplots()
ax.semilogx(R_hole/AU,Temp_innerhole, label= Metallicity)
ax.set_xlabel("Radius")
ax.set_ylabel("Interior_Temperature")
ax.set_title("Interior_Temperature vs. Radius")
plt.show()   



