import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import fsolve
import itertools


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
rstar = 6.4 * rs
tstar = 9500
lstar = 47* ls

#Disk parameter
T_rim= t_sublimation = 1500.0  # Temperature at which dust start to condense
sigma0 = 2*10 **3  # Surface density at 1 AU  Σ0 [g/cm^2] 
kp_stellar = 400.0  # Monochromatic opacity = 400cm^2 g^-1

nu = 5.879 * (10 ** 10) * tstar  #Wien law
t_virial = GG * mstar * mu * mp / (kb * rstar)  # Virial Temp
gamma = 2.0
beta = -1.5

psi_i =np.array(list(itertools.repeat(1, nr)))
psi_s = np.array(list(itertools.repeat(1, nr)))


"""
#Reading Planck_data file and converting as a list
df = pd.read_csv("Planck_data.csv")
T_data = np.array(list(df["T_gas"]))
kp_data = np.array(list(df["Planck_opacity"]))
kp_t_array = np.zeros((2, len(kp_data)))
kp_t_array[0] =T_data
kp_t_array[1] =kp_data
#change_array = T_kp_data.to_numpy()

# Assuming kp_t_array is a 2 x n array, with two columns.  Column 1 contains
# Temperature values, and column 2 contains corresponding Kp(T) values.
kp_t_array = np.zeros((2, 1000))


def kp_interp(T_guess):
    #Interpolate Temperature to find accurate Kp(T) value.
    kp_t_guess = np.interp(T_guess, kp_t_array[0], kp_t_array[1])
    return kp_t_guess


def t_i_eq(T_iter, const):
    #Solve iteration for Ti.
    psi_i_iter = kp_interp(T_iter)
    return const * (psi_i_iter ** -0.25) - T_iter


def t_s_eq(T_iter, const):
    #Solve iteration for Ts.
    eps_iter = (kp_interp(tstar) / kp_interp(T_iter)) ** 0.25
    return eps_iter * const - T_iter
"""


#Solve for angle of Direct stellar radiation impinges onto the disk 
def alpha_angle(alpha, H_cg, const):
    # const = (0.4*rstar + (gamma- 1))(1/r_i)
    return (const* H_cg) - alpha

const_alpha = (0.4*rstar + (gamma- 1))(1/r_i)
impinge_angle= fsolve(alpha_angle, 0.2, args=const_alpha)



#Solving for interior temperature with constant phi_i and phi_s
def t_i_eq(T_iter, const):
    #const_temp = (a_i*p_s/2/p_i)**0.25 *(rstar/r_i)**0.5 *tstar
    return const * (psi_i_iter ** -0.25) - T_iter

const_temp = (a_i*p_s/2/p_i)**0.25 *(rstar/r_i)**0.5 *tstar
Ti= fsolve(t_i_eq,T_rim, args=const_temp)
