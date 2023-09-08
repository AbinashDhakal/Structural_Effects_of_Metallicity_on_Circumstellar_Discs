
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

#kP(T*) = 2100 cm2 g-1;   https://iopscience.iop.org/article/10.3847/1538-4357/835/2/230/pdf 
#Kp(T_rim) = 700 cm2 g-1




nu = 5.879 * (10 ** 10) * tstar  #Wien law
t_virial = GG * mstar * mu * mp / (kb * rstar)  # Virial Temp
gamma = 2.0
beta = -1.5

#psi_i =np.array(list(itertools.repeat(1, nr)))
#psi_s = np.array(list(itertools.repeat(1, nr)))

psi_i = 1
psi_s = 1


def alpha_angle(H_cg,R):
    return (0.4*rstar/R   +(gamma - 1)*H_cg/R)

def disk_temp(impinge_angle, R, psi_i, psi_s):
    return ((impinge_angle*psi_s/2/psi_i)**0.25 * (rstar/R)**0.5 *tstar)

def pressure_height(Ti, R):
    return ((Ti/t_virial)**0.5 * (R/rstar)**0.5 * R)

def scale_height(X_cg, h_cg):
    return (X_cg *h_cg)

def sublimation_point(H_cg, R_rim_iter):
    return (lstar/(4*np.pi*T_rim**4 * ss))**0.5 *(1+ H_cg/R_rim_iter)**0.5 * -R_rim_iter
    
#using fslove and nonlinear problem to sovle the solution for H_cg

def surface_density(R):
    return sigma0*(R/AU)**beta





"""
from scipy.special import erfinv
def chi_disk(X_cg_iter, *data):
    impinge_angle, sigma = data
    return    ( 1- math.erf(X_cg_iter) - 2*impinge_angle/(sigma*kp_stellar))
"""

from scipy.special import erfinv

def chi_disk(impinge_angle, sigma):
    return np.sqrt(2)*erfinv(1- 2*impinge_angle/(sigma*kp_stellar))


def chi_rim(sigmaRim):
    return (erfinv(1-1/(4*sigmaRim*kp_stellar)))
    #return 1-  (4*sigmaRim*kp_stellar)*(1-math.erf(X_rim_iter))


def rim_function(variable):
    x1, x2, x3, x4, x5 = variable
    #R_rim, h_rim, H_rim, X_rim, sigma_rim
    f =[((lstar/(np.pi* 4*T_rim**4 *ss))**0.5 *(1+ x3/x1)**0.5 -x1),
           (((kb*T_rim*x1**3)/(mu*mp*GG*mstar))**0.5 -x2),
           (x2*x4 - x3),
           (4*x5*kp_stellar*(1-math.erf(x4)) - 1),
           (sigma0*(x1/AU)**beta - x5)]
    return f


pi = np.pi
def jacobian_rim(variable):
    x1, x2, x3, x4, x5 = variable
    #R_rim, h_rim, H_rim, X_rim, sigma_rim

    J = [[-0.25*x3*(lstar/(T_rim**4*ss))**0.5/(pi**0.5*x1**2*(x3/x1 + 1)**0.5) - 1,       0,          0.25*(lstar/(T_rim**4*ss))**0.5/(pi**0.5*x1*(x3/x1 + 1)**0.5),       0,     0],
		[1.5*(x1**3*T_rim*kb/(GG*mp*mstar*mu))**0.5/x1 ,                         -1    , 0,0,0],
		[0 ,      x4,    -1,    x2,0],
		[0,0,0,-8*kp_stellar*x5*np.exp(-x4**2)/np.sqrt(pi),4*kp_stellar*(1 - math.erf(x4))],
        [beta*sigma0*(x1/AU)**beta/x1, 0, 0, 0, -1]]
    return J


def iterative_newton(fun, x_init, jacobian):
    max_iter = 5000
    epsilon = 1e-8

    x_last = x_init
    listt = []
    for k in range(max_iter):
        # Solve J(xn)*( xn+1 - xn ) = -F(xn):
        J = np.array(jacobian(x_last))
        F = np.array(fun(x_last))

        diff = np.linalg.solve( J, -F )
        #listt.append(round(diff,2))
   
        x_last = x_last + diff
        #print("I am trying to converge:  ",np.linalg.norm(diff) ,k )
        listt.append(round(np.linalg.norm(diff), 2))
        # Stop condition:
        if np.linalg.norm(diff) < epsilon:
            #print('convergence!, nre iter:', k )
            break

    else: # only if the for loop end 'naturally'
        print('not converged')
    
    #print(listt)
    return x_last

x_init_rim =[0.539*AU, 0.019*AU,0.09344*AU,1.53,4500]
x_sol_rim = iterative_newton(rim_function, x_init_rim, jacobian_rim)

R_rim, h_rim, H_rim, X_rim, sigmaRim = x_sol_rim[0], x_sol_rim[1], x_sol_rim[2], x_sol_rim[3], x_sol_rim[4]

#Definig grid system in radial direction
rin = R_rim # Sublimation radius  # inner rim of disk. Density is assumed to
rout = 100*AU
R =ri = np.linspace(rin,rout,nr) 
#R   = 0.5 * ( ri[0:nr] + ri[1:nr+1] )             # Take the in-"between values
R = np.array(R)




alpha_guess = 0.2

def converging_disk(alpha_guess, r_ii):
    a_i = 0.22
    while np.round(a_i, 5) != np.round(alpha_guess,5):
    #while alpha_guess !=a_i:
        alpha_guess = a_i
        T_i =disk_temp(alpha_guess, r_ii,psi_i, psi_s)
        h_i= pressure_height(T_i,r_ii)
        s_i = surface_density(r_ii)
        X_i = chi_disk(alpha_guess, s_i)
        H_i = X_i*h_i
        a_i =alpha_angle(H_i, r_ii)
        #print(alpha_guess, a_i)
        
    return a_i



X_cg = np.zeros(nr)
impinge_angle = np.zeros(nr)
H_cg = np.zeros(nr)
h_cg = np.zeros(nr)
sigma = np.zeros(nr)
Ti = np.zeros(nr)
   

for i, r_i in enumerate(R):
    impinge_angle_iter =converging_disk(alpha_guess, r_i)
    alpha_guess = impinge_angle_iter
    Ti_iter = disk_temp(alpha_guess, r_i,psi_i, psi_s)
    h_cg_iter = pressure_height(Ti_iter,r_i)
    sigma_iter = surface_density(r_i)
    X_cg_iter= chi_disk(alpha_guess, sigma_iter)
    H_cg_iter = X_cg_iter*h_cg_iter
    #impinge_angle_iter =converging(alpha_guess, r_i)    
    #np.append(Ti, Ti_iter)
    Ti[i] = Ti_iter
    h_cg[i] =h_cg_iter
    sigma[i] =sigma_iter
    X_cg[i] =X_cg_iter
    H_cg[i] =H_cg_iter
    impinge_angle[i] =impinge_angle_iter




"""

for i, r_i in enumerate(R):
    Ti_iter = disk_temp(alpha_guess, r_i)
    h_cg_iter = pressure_height(Ti_iter,r_i)
    sigma_iter = surface_density(r_i)
    X_cg_iter= chi_disk(alpha_guess, sigma_iter)
    #print(X_cg_iter)
    #data = (alpha_guess, sigma_iter)
    #X_cg = fsolve(chi_disk,0.7,args = data)
    #print(h_cg_iter)
    H_cg_iter = X_cg_iter*h_cg_iter
    impinge_angle_iter =converging_disk(alpha_guess, r_i)
    print(alpha_guess, impinge_angle_iter)

    #impinge_angle_iter = alpha_angle(H_cg, r_i)
    alpha_guess = impinge_angle_iter
    
    Ti_iter = disk_temp(alpha_guess, r_i)
    h_cg_iter = pressure_height(Ti_iter,r_i)
    sigma_iter = surface_density(r_i)
    X_cg_iter= chi_disk(alpha_guess, sigma_iter)
    H_cg_iter = X_cg_iter*h_cg_iter
    #impinge_angle_iter =converging(alpha_guess, r_i)    
    #np.append(Ti, Ti_iter)
    Ti[i] = Ti_iter
    h_cg[i] =h_cg_iter
    sigma[i] =sigma_iter
    X_cg[i] =X_cg_iter
    H_cg[i] =H_cg_iter
    impinge_angle[i] =impinge_angle_iter
"""    
 
z_elevation =zi = np.linspace(0.0, H_cg, nz)
#z_elevation = np.array(z_elevation)
#z_elevation   = np.array([0.5 * ( zi[i] + zi[i+1])   for i in range(len(zi) - 1)])  # Take the in-"between values       # Take the in-between values    
Z = np.linspace(0.0, 50*AU, nz)
def Rho(R, Z):

    impinge_angle = (0.4*rstar/R   +(gamma - 1)*H_cg/R)
    Ti = ((impinge_angle*psi_s/2/psi_i)**0.25 * (rstar/R)**0.5 *tstar)
    h_cg = ((Ti/t_virial)**0.5 * (R/rstar)**0.5 * R)

    return (sigma0*(R/AU)**beta)*np.exp(- Z**2 / 2/h_cg**2) / (np.sqrt(2*np.pi) *h_cg)
rho = Rho(R,z_elevation)




"""
def Rho(z_elevation, h_cg, sigma):
    return sigma*np.exp(- z_elevation**2 / 2/h_cg**2) / (np.sqrt(2*np.pi) *h_cg)

"""










"""

def Rho(z_elevation, h_cg, R):
    return sigma0*(R/AU)**beta*np.exp(-z_elevation**2/2/h_cg**2)/(np.sqrt(2*np.pi)*h_cg)
    

H_grid, R_grid = np.meshgrid(H_cg, R)
h_grid, R_grid = np.meshgrid(h_cg, R)
sigma_grid, R_grid  = np.meshgrid(sigma, R)
#rho = Rho(z_elevation, h_cg, sigma_grid)
rho = Rho(z_elevation, h_grid,R_grid)
plt.pcolormesh( R_grid,z_elevation,  rho, shading='auto')
plt.colorbar()

plt.show()


"""











R_grid = np.meshgrid(R)
H_grid = np.meshgrid(H_cg)
h_grid = np.meshgrid(h_cg)
sigma_grid = np.meshgrid(sigma)
#rho = Rho(z_elevation, h_cg, sigma_grid)
rho = Rho(R_grid, H_grid)

R_grid = np.meshgrid(R)
plt.pcolormesh( R_grid,z_elevation,  rho, shading='auto')
plt.colorbar()

plt.show()

#Pressure [dyn/cm^2] = [g cm/s^2/cm^2] = [g/s^2/cm]
def GasPressure(rho, T_disk):
    return (kb*T_disk*rho/mu/mp)
Gas_pressure =GasPressure(rho, Ti)



#Force due to gravity
#g(z) = GMz/R^3

def F_gravity(z,R):
    return GG*mstar*z/R**3



#Gas pressure distribution
#dP/dz = -rho F_gravity
# dP [Pc,Pz] = (mu*G*M/(RR*T R**3)) ,RR is ideal gas consntat RR = 8.3145\,  J \cdot mol^{-1}\cdot K^{-1}
RR = 8.31446261815324*10**7	#erg⋅K−1⋅mol−1
#|Pz = Pc* exp(-mu GM z^2/(2RR T R^3))

#Surface density: Σ = Integral of ( rho dz) between [-inf, inf]
#surface density: Σ = Σ0(R/AU)^-1.5

#Σ = integral[-inf, inf] of (mu*Pc/(RR T)) exp(-mu GG M/(2RR T R^3)) dz
#Σ/ Pc = (2 pi mu R^3/(RR TGM)) , Pc = central pressure
def central_pressure(sigma, R, Ti):
    return (2* np.pi* mu* R**3/(RR*Ti*GG*mstar))**-0.5 * sigma
Pc = central_pressure(sigma, R, Ti)

    
def pressure_distribution(z_elevation, Ti, R, Pc):
    return (Pc* np.exp(-mu*GG*mstar*z_elevation**2/(2*RR*Ti/R**3)))
Pz = pressure_distribution(z_elevation, Ti, R, Pc)



import pandas as pd
df = pd.read_csv("Planck_data.csv")
filtering_data =( df[
    (df["T_stellar"] == 3000) 
    & (df["Metallicity"] == -0.3 )
   ])
T_data = np.array(list(filtering_data["T_gas"]))
kp_data = np.array(list(filtering_data["Planck_opacity"]))
kp_t_array = np.zeros((2, len(kp_data)))
kp_t_array[0] =T_data
kp_t_array[1] =kp_data
#change_array = T_kp_data.to_numpy()


def kp_interp(T_guess):
    """Interpolate Temperature to find accurate Kp(T) value."""
    kp_t_guess = np.interp(T_guess, T_data, kp_data)
    #kp_t_guess = np.interp(T_guess, T_data, Kp_data)

    return kp_t_guess


kp_tstar = kp_interp(tstar)

def surface_Temp(epsilon,R):
    return (1/epsilon)**0.25 * (rstar/2/R)**0.5 *tstar


w = 0.999999

    
def converging_temperature(r_ii, psi_s, psi_i, angle_i):
    T_guess = 1200
    r_ii = 0.52*AU
    T_s = 1200
    Tii = 1200
    error = 10

    while error >0.0001:
        T_guess = T_s
        epsilon = kp_interp(T_guess)
        T_s = surface_Temp(epsilon, r_ii)
        psi_s= kp_interp(T_s)
        s_i = surface_density(r_ii)
        psi_i = s_i *kp_interp(Tii)
        Tii = disk_temp(0.19459029132621866, r_ii, psi_i, psi_s)
        error = np.abs(T_s -T_guess)
        print(T_guess, T_s,Tii)
    return  Tii,T_s




Ts = np.zeros(nr)

for i, r_i in enumerate(R):
    impinge_angle_iter =converging_disk(alpha_guess, r_i)
    alpha_guess = impinge_angle_iter
    Ti_iter, Ts_iter = converging_temperature(r_i, psi_s, psi_i, impinge_angle_iter)
    h_cg_iter = pressure_height(Ti_iter,r_i)
    sigma_iter = surface_density(r_i)
    X_cg_iter= chi_disk(alpha_guess, sigma_iter)
    H_cg_iter = X_cg_iter*h_cg_iter
    #impinge_angle_iter =converging(alpha_guess, r_i)    
    #np.append(Ti, Ti_iter)
    Ti[i] = Ti_iter
    Ts[i] = Ts_iter
    h_cg[i] =h_cg_iter
    sigma[i] =sigma_iter
    X_cg[i] =X_cg_iter
    H_cg[i] =H_cg_iter
    impinge_angle[i] =impinge_angle_iter
#converging_temperature(1400, R[0], 1, 1, impinge_angle[0])
fig, ax = plt.subplots()
ax.semilogx(R/AU, impinge_angle)

ax.set_xlabel("Radius")
ax.set_ylabel("angle")
ax.set_title("angle vs. Radius")
plt.show()   
    



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
    

fig, ax = plt.subplots()
ax.semilogx(R/AU, Ts)


ax.set_xlabel("Radius")
ax.set_ylabel("Disk Temperature")
ax.set_title("Disk Temperature vs. Radius")
plt.show()   
    

fig, ax = plt.subplots()
ax.semilogx(R/AU, X_cg)


ax.set_xlabel("Radius")
ax.set_ylabel("X_cg")
ax.set_title("x_cg vs. Radius")
plt.show()   


