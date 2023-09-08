import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfinv
import pandas as pd
from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator

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

sigmagas0 = 1.7*2*10**3 #https://www.ita.uni-heidelberg.de/~dullemond/lectures/radtrans_2012/Chapter_8.pdf
#sigma0=1700
sigma0 = 2*10**3
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

#Midplane density profile :Rho(0, h_cg, sigma)
def Rho(z_elevation, h_cg, sigma):
    return sigma*np.exp(- z_elevation**2 / 2/h_cg**2) / (np.sqrt(2*np.pi) *h_cg)

#Pressure [dyn/cm^2] = [g cm/s^2/cm^2] = [g/s^2/cm]
def GasPressure(rho, Tgas):
    return (kb*Tgas*rho/mu/mp)

#where sigma0 = 1000gcm^-2
#Surface density
def gas_density(R):
    return sigmagas0*(R/AU)**(-1.5)


R_rim = 0.5399999944538275*AU
H_rim = 0.07581315724362893*AU
Metallicity =-0.3
T_stellar =tstar 
T_Stellar =9000
T_rim = 1500
R =np.linspace(rstar,R_rim,nr) 
#R   = 0.5 * ( ri[0:nr] + ri[1:nr+1] )             # Take the in-"between values
R = np.array(R)



#THis should interpolate my Kp(Tgas,Trad)
#Its treating input data with 3D array and 4th dimension as a result (kp(Trad, Tgas))

def kp_TgasTrad_interpolation(Input_array,Output_array ,request):
    Interp= LinearNDInterpolator(Input_array, Output_array, rescale = True)
    output = Interp(request)
    #It linearly interpolate from pre recorded data and try to find closest Match
    iteration = 0
    request_inital = request
    
    #IF its returning noting, they try to re-evaluate until it output something
    while output != output:
        request_again =  np.round(request*np.random.uniform(0.9, 1.1),1)
        print( "kp_TgasTrad_interpolation",request_again)
        request = request_again
        output = Interp(request_again)
        #print(output)
        iteration += 1
        
        #IF it cant find anythig, then try to find someething very close to that value
        #Find the index of close value and then usee that index to calculate kp

        if iteration ==10:
        # Vector we try to match with data
            search_vect = request_inital
            # normalize data in each column with the column mean
            norm_data = Input_array/np.mean(Input_array, axis=0)
            # normalize the search vector with the same values
            s_vec_n = search_vect/np.mean(Input_array, axis=0)
            # Define closest fit as the smallest norm of the difference vector
            idx = np.argmin(np.linalg.norm((s_vec_n - norm_data), axis=1))
            output = Output_array[idx]

    return output


#THis uses same procedure as Kp_Tgas_Trad but the data frame is different
def kp_Tgas_interpolation(Input_array,Output_array ,request):
    #Rescale = True : IT rescale all the element interm of [0,1]
    linInter= LinearNDInterpolator(Input_array, Output_array,rescale = True)
    output = linInter(request)
    request_inital = request
 
    iteration = 0
    while output != output:
        request_again =  np.round(request*np.random.uniform(0.9, 1.1),1)
        print( "kp_Tgas_interpolation",request_again)
        request = request_again
        output = linInter(request_again)
        #print(output)
        iteration += 1
        if iteration ==10:
        # Vector we try to match with data
            search_vect = request_inital
            # normalize data in each column with the column mean
            norm_data = Input_array/np.mean(Input_array, axis=0)
            # normalize the search vector with the same values
            s_vec_n = search_vect/np.mean(Input_array, axis=0)
            # Define closest fit as the smallest norm of the difference vector
            idx = np.argmin(np.linalg.norm((s_vec_n - norm_data), axis=1))
            output = Output_array[idx]

    return output


df = pd.read_csv("Planck_data.csv")
df = df.loc[df['T_stellar']==T_Stellar]
df = df.loc[df['Metallicity']==Metallicity]
#Rho and P gas changes exponentially 
#When interpolating linearly, its harder to interpolate
df['Log_Rho_gas'] = np.log10(np.array(df['Rho_gas'])) 
df['Log_P_gas'] = np.log10(np.array(df['P_gas']))    
Input_array = np.array(df[['T_gas','Log_P_gas', 'Log_Rho_gas']])
Output_array = np.array(df["kp_stellar"])

#Md = Meta data
#Input array = [Metallicity(Implicitly applied),
                                #Tgas, log10(Pgas), log10(Rho_gas)]
#Output array correspond to kp(Input)
md = pd.read_csv("Planck_gas.csv")
md= md.loc[md['Metallicity']==Metallicity]
md['Log_Rho_gas'] = np.log10(np.array(md['Rho_gas'])) 
md['Log_P_gas'] = np.log10(np.array(md['P_gas'])) 
Input_array_gas = np.array(md[['T_gas','Log_P_gas', 'Log_Rho_gas']])
Output_array_gas = np.array(md["kp_stellar"])


request1= np.round(np.array([700	,-8.71444,	-19.1129]),2)
request2 =  np.round(np.array([700	,2.28556,	-8.11295]),2)
request3 = np.round(np.array([11007,0.285557,	-11.8633]),2)
request4 =  np.round(np.array([82540	,8.57171,	-4.46344]),2)

print(kp_Tgas_interpolation(Input_array_gas,Output_array_gas ,request1))
print(kp_Tgas_interpolation(Input_array_gas,Output_array_gas ,request2))
print(kp_Tgas_interpolation(Input_array_gas,Output_array_gas ,request3))
print(kp_Tgas_interpolation(Input_array_gas,Output_array_gas ,request4))

"""

request1= np.round(np.array([700	,-9	,-19.399]),2)
request2 =  np.round(np.array([7100,	-0.428291,	-12.3143]),2)
request3 = np.round(np.array([74989,	8.42813	,-4.56864]),2)
request4 =  np.round(np.array([82540	,8.57171,	-4.46344]),2)
print(kp_TgasTrad_interpolation(Input_array,Output_array ,request1))
print(kp_TgasTrad_interpolation(Input_array,Output_array ,request2))
print(kp_TgasTrad_interpolation(Input_array,Output_array ,request3))
print(kp_TgasTrad_interpolation(Input_array,Output_array ,request4))
"""