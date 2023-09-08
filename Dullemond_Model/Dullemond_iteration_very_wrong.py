import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import fsolve
import itertools
import math
from scipy.optimize import fsolve



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

#psi_i =np.array(list(itertools.repeat(1, nr)))
#psi_s = np.array(list(itertools.repeat(1, nr)))

psi_i = 1
psi_s = 1
def Rim(vars):
    R_rim, h_rim, H_rim, X_rim,sigma_rim = vars

    #Define R_rim
    eqn1 = np.abs((lstar/(4*np.pi * T_rim**4 *ss))**0.5 * (1 + H_rim/R_rim)**0.5 - R_rim)
    #Define h_rim
    eqn2 =np.abs((kb*T_rim*R_rim**3/(mu*mp*GG*mstar))**0.5 - h_rim)
    #Define H_rim
    eqn3 =np.abs( X_rim * h_rim -H_rim)

    #Define X_rim
    #integral[H/h = X  --- inf] of e^-x^2 dx =(Pi^0.5 /2)(1- erf(X_rim))
    # Tau(zo) = 8*sigma* kp_stellar* integral/(pi^0.5) 
    eqn4 = np.abs((1-math.erf(X_rim))*(4*sigma_rim*kp_stellar) -1)
    #Define surface density
    eqn5 = np.abs(sigma0*(R_rim/AU)**beta -sigma_rim)
    return [eqn1, eqn2,eqn3,eqn4, eqn5]


R_rim, h_rim, H_rim,X_rim, sigma_rim = fsolve(Rim, (0.52*AU, 0.1*AU,0.057*AU,5.3,6207))


#Definig grid system in radial direction
rin = 0.52 *AU # Sublimation radius  # inner rim of disk. Density is assumed to
rout = 100*AU
#R = ri = np.linspace(rin,rout,nr+1) 
R = ri = np.linspace(rin,rout,nr) 

#R   = 0.5 * ( ri[0:nr] + ri[1:nr+1] )             # Take the in-"between values
R = np.array(R)

#sigma[g/cm^2] 
def Sigma(R):
    return sigma0 * (R/AU)**beta

sigma = Sigma(R)
sigma = np.array(sigma)

import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


sample_data = pd.read_csv('Ti_dataset.csv')
sample_data2 = pd.read_csv('hR_datasets.csv')
sample_data3= pd.read_csv('HHHR_datasets.csv')


R_iterpolate_range = R/AU
Temp_interpolate_f =interp1d(sample_data.R, sample_data.Ti, 'quadratic')
hR_interpolate_f =interp1d(sample_data2.R, sample_data2.h, 'quadratic')
HR_interpolate_f =interp1d(sample_data3.R, sample_data3.H, 'quadratic')

Ti_interpolate = Temp_interpolate_f(R_iterpolate_range)
hR_interpolate = hR_interpolate_f(R_iterpolate_range)
h_interpolate = hR_interpolate_f(R_iterpolate_range)*AU
HR_interpolate = np.array(HR_interpolate_f(R_iterpolate_range))
H_interpolate = np.array(HR_interpolate_f(R_iterpolate_range))*AU

X_interpolate = HR_interpolate/hR_interpolate
alpha_interpolate = (0.4*rstar/R_iterpolate_range/AU) + (gamma -1)*HR_interpolate
plt.semilogx(R_iterpolate_range, Ti_interpolate)



def rim_point(R_iter, const):
    #const_Rrim = (lstar/(4*np.pi*T_rim**4 *ss))**0.5  *(1+ H_i/r_i)
    return np.abs(const - R_iter)

#const_Rrim = (lstar/(4*np.pi*T_rim**4 *ss))**0.5  *(1+ H_i/r_i)
#R_rim = fsolve(rim_point, 0.1*AU,const_Rrim)

def surface_density(const_sigma):
    return const_sigma

def chi_Rim(X_iter, sigma):
    return np.abs((4*sigma*kp_stellar)*(1-math.erf(X_iter))) -1

#Solve for angle of Direct stellar radiation impinges onto the disk 
def alpha_angle(alpha_iter, const):
    # const = (0.4*rstar + (gamma- 1)*H_i)(1/r_i)
    return np.abs(const - alpha_iter)

def alpha_angle_matrix(alpha_iter, const):
    # const = (0.4*rstar + (gamma- 1)*H_i)(1/r_i)
    return [(const - alpha_iter)]

#Solving for interior temperature with constant phi_i and phi_s
def temp_interior(T_iter, const):
    #const_temp = (a_i*p_s/2/p_i)**0.25 *(rstar/r_i)**0.5 *tstar
    return np.abs(const  - T_iter)

def temp_interior_matrix(T_iter, const):
    #const_temp = (a_i*p_s/2/p_i)**0.25 *(rstar/r_i)**0.5 *tstar
    return [(const  - T_iter)]

def pressure_height(hcg_iter,const):
    #const_hcg = (T_i/t_virial)**0.5 * (r_i/rstar)**0.5*r_i
    return np.abs(const - hcg_iter)


def pressure_height_matrix(hcg_iter,const):
    #const_hcg = (T_i/t_virial)**0.5 * (r_i/rstar)**0.5*r_i
    return [(const - hcg_iter)]


def surface_height(Hcg_iter, const):
    #const_Hcg = X_i * h_i
    return np.abs(const - Hcg_iter)

def surface_height_matrix(Hcg_iter, const):
    #const_Hcg = X_i * h_i
    return [(const - Hcg_iter)]

def jacobian_(Iter):
    #R_rim, h_rim, H_rim, X_rim, sigma_rim
    J = [-1]
   
    return J


def chi_disk(X_iter, const):
    #const_X = 2*a_i/kp_stellar/s_i
    #const_X = 2*a_i/kp_stellar/s_i
    return np.abs(1- math.erf(X_iter/(np.sqrt(2))) - const)

def chi_disk_matrix(X_iter, const):
    #const_X = 2*a_i/kp_stellar/s_i
    #const_X = 2*a_i/kp_stellar/s_i
    return [(1- math.erf(X_iter/(np.sqrt(2))) - const)]

def jacobian_chi(X_iter):
    return[(-np.sqrt(2)*np.exp(-X_iter**2/2)/np.sqrt(np.pi))]



def iterative_newton(func, x_init, jacobian, const):
    max_iter = 50
    epsilon = 1e-8

    x_last = x_init

    for k in range(max_iter):
        # Solve J(xn)*( xn+1 - xn ) = -F(xn):
        J = np.array(jacobian(x_last))
        F = np.array(func(x_last, const))

        diff = np.abs( J -F )
        x_last = x_last + diff
        print("I am trying to converge:  ",np.linalg.norm(diff) ,k )

        # Stop condition:
        if np.abs(diff) < epsilon:
            print('convergence!, nre iter:', k )
            print(x_last)
            break

    else: # only if the for loop end 'naturally'
        print('not converged')

    return x_last








X_cg = np.zeros(nr)
impinge_angle = np.zeros(nr)
H_cg = np.zeros(nr)
h_cg = np.zeros(nr)
sigma = np.zeros(nr)
Ti = np.zeros(nr)
abcd = np.zeros(nr)

psi_i = 1
psi_s = 1

Simulation = 5
#print("I am this value outside the boxk :",X_cg_trial[0])





#This doesnt update after 1st simulation
for z in range(Simulation):
    for i, (r_i, H_i, h_i, T_i,X_i, a_i, s_i) in enumerate(zip(R, H_interpolate, h_interpolate, Ti_interpolate, X_interpolate,alpha_interpolate, sigma)):
       #print(i)
        #Simulation -=1
        const_alpha = (0.4*rstar + (gamma- 1)*H_i)*(1/r_i)
        x_sol= iterative_newton(alpha_angle_matrix,[a_i], jacobian_,const_alpha)
        #print(const_alpha)
        impinge_angle[i]= fsolve(alpha_angle,x_sol , args=const_alpha)
        #I have multiplied argument with 0.95 as it would be different from measured value
        #and potentially produce new answer that converges
        

        const_temp = (a_i*psi_s/2/psi_i)**0.25 *(rstar/r_i)**0.5 *tstar
        #print("Const temperaureate is ",const_temp)
        x_sol= iterative_newton(temp_interior_matrix,[T_i], jacobian_,const_temp)

        Ti[i]= fsolve(temp_interior,x_sol, args=const_temp)


        const_hcg = (T_i/t_virial)**0.5 * (r_i/rstar)**0.5 *r_i
        x_sol= iterative_newton(pressure_height_matrix,[h_i], jacobian_,const_hcg)
        h_cg[i] = fsolve(pressure_height,x_sol, args = const_hcg)


        const_Hcg = X_i * h_i
        x_sol= iterative_newton(surface_height_matrix,[H_i], jacobian_,const_Hcg)

        H_cg[i] = fsolve(surface_height,x_sol, args =const_Hcg )
        #print(4)
       #print("Help4\n\n\n\n")
        const_X = 2*a_i/kp_stellar/s_i
        x_sol= iterative_newton(chi_disk_matrix,[X_i], jacobian_chi,const_X)

        X_cg[i] = fsolve(chi_disk, x_sol, args = const_X)
        
        
       #print("Hello, I want to check this :  ", X_cg_trial[0], X_cg[0])
        
    #I dont know why but this doesnt over lay is list
    X_interpolate = X_cg.copy()
        #print(X_cg_trial[i])
    alpha_interpolate= impinge_angle.copy()
    H_interpolate = H_cg.copy()

    h_interpolate = h_cg.copy()
    Ti_interpolate =  Ti.copy()


"""

#This doesnt update after 1st simulation
for z in range(Simulation):
    for i, (r_i, H_i, h_i, T_i,X_i, a_i, s_i) in enumerate(zip(R, H_interpolate, h_interpolate, Ti_interpolate, X_interpolate,alpha_interpolate, sigma)):
       #print(i)
        #Simulation -=1
        #print(Simulation)
       #print("Hello\n\n\n\n")
        #print(r_i, H_i, h_i, T_i,X_i, a_i, s_i )
       #print("Help\n\n\n\n")
        const_alpha = (0.4*rstar + (gamma- 1)*H_i)*(1/r_i)
        x_sol= iterative_newton(alpha_angle_matrix,[a_i], jacobian_,const_alpha)
        print(x_sol)
        #print(const_alpha)
        impinge_angle[i]= fsolve(alpha_angle,a_i , args=const_alpha)
        #I have multiplied argument with 0.95 as it would be different from measured value
        #and potentially produce new answer that converges
        
        
        #print(impinge_angle)
       #print("Help1\n\n")
        #print(1)
        const_temp = (a_i*psi_s/2/psi_i)**0.25 *(rstar/r_i)**0.5 *tstar
        #print("Const temperaureate is ",const_temp)
        print(T_i)
        Ti[i]= fsolve(temp_interior,T_i, args=const_temp)
        print(T_i, "Afterrrrr")
        #print(2)
       #print("Help2\n\n")

        const_hcg = (T_i/t_virial)**0.5 * (r_i/rstar)**0.5 *r_i
        h_cg[i] = fsolve(pressure_height,h_i, args = const_hcg)
        #print(3)
       #print("Help3\n\n")

        const_Hcg = X_i * h_i
        #const_X = 2*a_i/kp_stellar/ 

        H_cg[i] = fsolve(surface_height,H_i, args =const_Hcg )
        #print(4)
       #print("Help4\n\n\n\n")
        const_X = 2*a_i/kp_stellar/s_i
        X_cg[i] = fsolve(chi_disk, X_i, args = const_X)
       #print("Hiiiiiiiiiiiiiiiiiiiiiiiiii")
        
       #print("Hello, I want to check this :  ", X_cg_trial[0], X_cg[0])
        
    #I dont know why but this doesnt over lay is list
    X_interpolate = X_cg.copy()
        #print(X_cg_trial[i])
    alpha_interpolate= impinge_angle.copy()
    H_interpolate = H_cg.copy()

    h_interpolate = h_cg.copy()
    Ti_interpolate =  Ti.copy()

"""



fig, ax = plt.subplots()
ax.semilogx(R/AU, H_cg/R)

ax.set_xlabel("Radius")
ax.set_ylabel("Surface_Height")
ax.set_title("Surface_Height vs. Radius")
ax.set_ylim([0, 0.5])

plt.show()


fig, ax = plt.subplots()
ax.semilogx(R/AU, h_cg/R)

ax.set_xlabel("Radius")
ax.set_ylabel("Pressure_Height")
ax.set_title("pressure_Height vs. Radius")
ax.set_ylim([0, 0.5])

plt.show()


fig, ax = plt.subplots()
ax.semilogx(R/AU, Ti)

ax.set_xlabel("Radius")
ax.set_ylabel("Temperature")
ax.set_title("Temperature vs. Radius")
plt.show()   
    
    