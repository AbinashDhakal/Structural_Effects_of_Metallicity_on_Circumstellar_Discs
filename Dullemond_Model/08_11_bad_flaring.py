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



def alpha_angle(H_cg,R):
    return (0.4*rstar/R   +(gamma - 1)*H_cg/R)


def irr_surface_Temp(R, Ts_guess):
    es =  40*(Ts_guess/tstar)**1.5
    if es >1:
        epsilon = 1
    else:
        epsilon = es
     
    return (1/epsilon)**0.25 *(rstar/2/R)* tstar 


Ts_guess = T_rim
#convering iteration such that the it iterate until the error is smaller than 0.00001
#THis only check for a radial point
"""
def converging_Temp(Ts_guess, r_ii):
    T_s = 30
    while np.round(T_s, 5) != np.round(Ts_guess,5):
    #while alpha_guess !=a_i:
        Ts_guess = T_s
        T_s =irr_surface_Temp(r_ii, Ts_guess)
        
    return T_s



Ts = converging_Temp(R, T_rim)
"""


#fraction of stellar flux that will be absorbed by the interior
def flux_absorbed_fraction(sigma, Ts_guess): #psi_s
    psi_s = 1

    if (sigma*kp_interp(Ts_guess)) < 1:
        return kp_interp(Ts_guess)
    else:
        return 1
         

#The possiblity that the disk interior is not fully optically thick to its own emission
def self_absorbed_fraction(sigma, Ti_guess): #psi_i
    if (sigma*kp_interp(Ti_guess)) < 1:
        return kp_interp(Ti_guess)*sigma
    else:
        return 1

def irr_disk_temp(impinge_angle, R, Ti_guess, Ts_guess):
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


import pandas as pd
def reading_data():

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
    return T_data, kp_data

T_data, kp_data = reading_data()
kp_t_array = np.zeros((2, len(kp_data)))
kp_t_array[0], kp_t_array[1] =T_data, kp_data

def kp_interp(T_guess):
    """Interpolate Temperature to find accurate Kp(T) value."""
    kp_t_guess = np.interp(T_guess, T_data, kp_data)
    #kp_t_guess = np.interp(T_guess, T_data, Kp_data)
    return kp_t_guess


kp_tstar = kp_interp(tstar)

"""
def irr_surface_Temp(R, T_guess):
    epsilon = kp_interp(T_guess)/kp_interp(tstar)
    T_s = (1/epsilon)**0.25 * (rstar/2/R)**0.5 *tstar
    w = 0.9999999
    error = 10
    while error >0.01:
        Temp = T_s
        epsilon = kp_interp(Temp)/kp_interp(tstar)
        #x[k+1] = w x[k] + (1-w) f(x[k+1])
        T_s = w*Temp +(1-w)*(1/epsilon)**0.25 * (rstar/2/R)**0.5 *tstar
        error = np.max(np.abs(Temp - T_s))

        
    return T_s
"""

def irr_surface_Temp(R, T_guess):
    epsilon = kp_interp(T_guess)/kp_interp(tstar)
    T_s = (1/epsilon)**0.25 * (rstar/2/R)**0.5 *tstar
    w = 0.9999999
    error = 10
    while error >0.01:
        Temp = T_s
        epsilon = kp_interp(Temp)/kp_interp(tstar)
        #x[k+1] = w x[k] + (1-w) f(x[k+1])
        T_s = w*Temp +(1-w)*(1/epsilon)**0.25 * (rstar/2/R)**0.5 *tstar
        error = np.max(np.abs(Temp - T_s))

        
    return T_s

T_s = irr_surface_Temp(R, T_rim)


from scipy.special import erfinv
def chi_disk(impinge_angle, sigma):
    return np.sqrt(2)*erfinv(1- 2*impinge_angle/(sigma*kp_stellar))


def chi_rim(sigmaRim):
    return (erfinv(1-1/(4*sigmaRim*kp_stellar)))
    #return 1-  (4*sigmaRim*kp_stellar)*(1-math.erf(X_rim_iter))


fig, ax = plt.subplots()
ax.semilogx(R/AU, T_s)

ax.set_xlabel("Radius")
ax.set_ylabel("Surface_temperature")
ax.set_title("Surface_tempeatureeeeeee vs. Radius")

plt.show()


alpha_guess = 0.2


#convering iteration such that the it iterate until the error is smaller than 0.00001
#THis only check for a radial point
def converging_disk(alpha_guess, r_ii):
    a_i = 0.22
    while np.round(a_i, 5) != np.round(alpha_guess,5):
    #while alpha_guess !=a_i:
        alpha_guess = a_i
        T_i =irr_disk_temp(alpha_guess, r_ii,psi_i, psi_s)
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
   


#This iterate through every radial points

for i, r_i in enumerate(R):
    impinge_angle_iter =converging_disk(alpha_guess, r_i)
    alpha_guess = impinge_angle_iter
    sigma_iter = surface_density(r_i)
    w = 0.9999999
    error = 10
    Ti_iter = irr_disk_temp(alpha_guess, r_i,psi_i, psi_s)
    #This check Ti multiple times such that error is less than 0.01
    while error >0.01:
        Temp = Ti_iter
        psi_i = self_absorbed_fraction(sigma_iter, Temp)
        psi_s = flux_absorbed_fraction(sigma_iter, T_s[i])
        print("Printing , psi value, ",psi_i, psi_s)
        Ti_iter =  w*Temp +(1-w)*irr_disk_temp(alpha_guess, r_i,psi_i, psi_s)
        error = np.max(np.abs(Temp - Ti_iter))
    #Use that value of T to compute h_cg and others
    
    h_cg_iter = pressure_height(Ti_iter,r_i)
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
    impinge_angle_iter =converging_disk(alpha_guess, r_i)
    alpha_guess = impinge_angle_iter
    Ti_iter = irr_disk_temp(alpha_guess, r_i,psi_i, psi_s)
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


def Rho(z_elevation, h_cg, sigma):
    return sigma*np.exp(- z_elevation**2 / 2/h_cg**2) / (np.sqrt(2*np.pi) *h_cg)

rho = Rho(z_elevation, h_cg, sigma)
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

R_flare = H_cg*R_rim/H_rim

def viscous_temp(R):
    
    yrs = 31557600 #s
    M_acc = 10**-8*ms/yrs
    return (3*GG*mstar*M_acc/(8*np.pi* ss* R**3))**0.25
    #return 107*(R/AU)**-0.75 * (mstar/2.5/ms)**0.25 *(M_acc/(10**-8 *ms/yrs))**0.25

T_visc = viscous_temp(R)

#This create a log-grid scale so if I were to integrate between [z, inf]
#This should allow me to focus on bottom point more when adding up

log_R      = np.logspace(np.log10(rout/AU),np.log10(rin/AU), nr,endpoint=True)

def Optical_depthhh(H_cg,h_cg, sigma, impinge_angle):
    I= 0 # Inegral
    Tau = np.zeros((nr,nz))
    for i,(H_i) in enumerate(H_cg):
        up_lim = H_i
        low_lim = H_i - i*nz
        Hz = ((up_lim - low_lim)/nz) #creating grid along vertical axisz
        for j in range (nz):
            zi = low_lim+ Hz/2 +j*Hz
            rho = kp_stellar/impinge_angle[i]*sigma[i]/(np.sqrt(2*np.pi)*h_cg[i]) * np.exp(-zi**2/h_cg[i]**2)
            I +=rho
        Tau[i][j] = I

    
    return Tau

    
    

TTau = Optical_depthhh(H_cg,h_cg, sigma, impinge_angle)


#Method 2
#This integrand should supposedly create a cuntion as const*ext(-z**2)
#then use bottom function to integrate
#I wnat 2D table containing tau value 
def intergrand(z):
    for i, (h_i, s_i) in enumerate(h_cg, sigma):
        rho = s_i/(np.sqrt(2*np.pi)*h_i) * np.exp(-z**2/h_i**2)*kp_stellar/impinge_angle[i]
    return rho


def Tau(H_cg):
    for j,(H_i) in enumerate(H_cg):
        integral, err =quad(intergrand,H_i -j*nz, H_i)
    return integral






Rflare = np.where (R <H_cg*R_rim/H_rim)
Rshadow = np.where (R >H_cg*R_rim/H_rim)

Ti[Rshadow] = 0

T_effective = ( T_visc**4  +Ti**4)**0.25

h_cg = pressure_height(T_effective, R)
H_cg = scale_height(X_cg, h_cg)
h_cg[0] = h_rim
H_cg[0] = H_rim
X_cg[0] = X_rim
T_effective[0] = T_rim





"""
def optical_depth(h_cg, z_elevation, impinge_angle,sigma):
    return sigma*kp_stellar/(2*impinge_angle) *(1-z_elevation/h_cg / np.sqrt(2))
    
Tau = optical_depth(h_cg, z_elevation,impinge_angle,sigma)
"""

T_effective = np.sqrt(T_visc**2 + T_s**2+ Ti**2)
#T_net    = T_viscousity + T_surface + T_inerior

listtt =np.zeros(len(H_cg))
for i, (H_i) in enumerate (H_cg):
    listtt[i] = H_i*R_rim/H_rim
    


h_redefine = pressure_height(T_effective, R)

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
ax.semilogx(R/AU, T_s)


ax.set_xlabel("Radius")
ax.set_ylabel("Disk Temperature")
ax.set_title("Disk Temperature vs. Radius")
plt.show()   

    
fig, ax = plt.subplots()
ax.semilogx(R/AU, viscous_temp(R))


ax.set_xlabel("Radius")
ax.set_ylabel("viscous_temp(R)")
ax.set_title("viscous_temp vs. Radius")
plt.show()   
    
fig, ax = plt.subplots()
ax.semilogx(R/AU, T_effective)


ax.set_xlabel("Radius")
ax.set_ylabel("Effective temp")
ax.set_title("Effective vs. Radius")
plt.show()   
    


fig, ax = plt.subplots()
ax.semilogx(R/AU, X_cg)


ax.set_xlabel("Radius")
ax.set_ylabel("X_cg")
ax.set_title("x_cg vs. Radius")
plt.show()   

fig, ax = plt.subplots()
ax.plot(R/AU, 150*(R/AU)**-(3/7))
ax.plot(R/AU,  550*(R/AU)**-(2/5))
   
"""

fig, ax = plt.subplots()
ax.semilogx(R/AU, impinge_angle)
ax.semilogx(R/AU, H_cg/R)
ax.semilogx(R/AU, 0.4* rstar/R)

"""

ax.set_xlabel("Radius")
ax.set_ylabel("Impinge_angle")
ax.set_title("Impinge_angle vs. Radius")
plt.show()   

fig, ax = plt.subplots()
ax.plot(H_cg/R, Ti)


ax.set_xlabel("Hcg/R")
ax.set_ylabel("Temp")
ax.set_title("Temp vs. scale_heigh/ Radius")
plt.show()   


R_rim = 0.54*AU
H_rim= 0.3*R_rim
fig, ax = plt.subplots()
ax.plot(R/AU, H_cg*R_rim/R/H_rim)


ax.set_xlabel("R_cg")
ax.set_ylabel("Flaring_point")
ax.set_title("Flaring_point vs. Radius")
plt.show()   
