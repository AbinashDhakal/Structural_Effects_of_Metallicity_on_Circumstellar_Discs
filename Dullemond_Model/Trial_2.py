
import numpy as np
import sympy as smp
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.optimize import fsolve


AU     =  1.495978707e13             # Distance between Earth and Sun [cm] defined 
GG      =  6.67430e-8                # Universal gravitational constant changing.... 
mu     =  2.353                      # Average molecular weight 
kb     =  1.380649e-16               # Boltzmann constant defined 
hp      =  6.62607015e-27            # Planck constant [erg s] defined
clight  =  2.99792458e10             # light speed [cm s^-1] defined 
NA     =  6.02214076e23              # Avogadro constant defined 
mP  =  1.672621898e-24               # proton mass [g]; not 1.0/NA anymore. 
pc     =  3.085677581e18             # Parsec [ cm] 
ms     =  1.9884e33                  # Sun mass [g] 
ts     =  5.777e3                    # effective temperature of the sun [K] 
ls     =  3.828e33                   # solar intensity [erg/s] defined by IAU 
rs     =  6.96e10                    # solar radius [cm] 
ss     = 5.6703e-5                   # Stefan-Boltzmann const  [erg/cm^2/K^4/s]



# grid parameters 
nr         =  2560          # Number of radial grids. 
ntheta     =  257        # Number of grids in elevation.
nz         =  2560          # Number of elevation grids. 

rin        =  0.1* AU     #Sublimation radius  # inner rim of disk. Density is assumed to be 0 inside [cm] 
rout       = 2.7*(10**2)*AU# Outer edge of disk. Density is set to 0 outside this [cm] 
theta_min = np.pi/2-0.18;           #Polar boundaries in radian
theta_max = np.pi/2+0.18;

#Stellar Parameter
mstar =0.5*ms
rstar =2.5*rs
tstar = 4000
lstar = 0.7*ls

#Disk parameter
t_sublimation = 1500            # Temperature at which dust start to condense
sigma0 = 10**3                  #Surface density at 1 AU  Σ0 =10^3 gcm^1
opacity_const =400             #Monochromatic opacity = 400cm^2 g^-1
psi_i = 1
psi_s = 1
X_cg = 2
nu = 5.879*10**10 * tstar       #Wien law
t_virial = GG* mstar* mu* mP/ (kb* rstar)
alpha =2

ri        =  np.linspace ( rin , rout , nr + 1 )
R   = 0.5 * ( ri[:-1] + ri[1:] )              # Take the in-"between values



def Temp(Ti, R):
    """Set equation for R_rim."""
    return ((0.4*rstar/R) + (alpha - 1)*(X_cg*(Ti/t_virial)**0.5  * (R/rstar)**0.5) /R)**0.25 * (rstar/R)**0.5 * tstar

t_0 = t_sublimation

eff_t_effective =  np.empty(len(R))
for i in range(len(eff_t_effective)):
    eff_t_effective[i] = fsolve(Temp,t_0,  args= R[i])
    h_rim = (eff_t_effective/t_virial)**0.25 * (R/rstar)**0.25
    H_cg = X_cg * h_rim
    xab = (0.4*rstar/R) + (alpha - 1)*H_cg/R

Tdust = np.empty(len(R))* 0.1
surface_density = (1/(2*opacity_const)*eff_t_effective/Tdust)
"""
# Numerically solve for R_rim
eff_r_rim = fsolve(func=r_rim, x0=r_0, factor=100)

# Calculate h_rim using eff_r_rim, and subsequent H_rim
h_rim = np.sqrt(np.sqrt((kb * t_sublimation / (mu * mP * GG * mstar)) * (eff_r_rim ** 3)))
H_rim = X_cg * h_rim

"""

"""
def impingement_angle (H_cg, R):
    gamma = 2
    impinging_angle = (2*rstar/(5*R)) + (gamma -1)(H_cg/R)
    return impinging_angle

impinging_angle = impingement_angle(H_cg, R)

T_inerior = (impinging_angle*psi_s/ (2*psi_i))**0.25 * (rstar/R)**0.5 *tstar

T_virial = GG*mstar*mu*mP/(kb*rstar)      #Given =8 x 10^6

def pressure_scale_height(T_inerrior,R):
    h_cg = (T_interior/ T_virial)**0.5 * (R/rstar)**0.5 * R
    return h_cg
h_cg= pressure_scale_height(T_inerrior, R)

H_cg = X_cg * h_cg


"""

T = range(500, 6000)
def planck_mean_opacity(T):
    planck_opacity = []

    for t in T:
        lambda_rel = hp * clight / (kb * t)     #defining the reference lambda
        lambda_min = lambda_rel / 10  #Defining lower int bound
        lambda_max = lambda_rel * 10  #Defining upper int bound        
        integrand1 = lambda x : ((2* nu**2 *opacity_const/ clight**2) * (hp * nu/(np.exp(hp*nu/(kb*t))-1)))
        integrand2 = lambda x : ((2* nu**2 / clight**2) * (hp * nu/(np.exp(hp*nu/(kb*t))-1)))
        result1,err1 = integrate.quad(integrand1, lambda_min, lambda_max, epsabs=0, limit=50)
        result2,err2= integrate.quad(integrand2, lambda_min, lambda_max, epsabs=0, limit=50)
        planck_opac = result1/result2
        planck_opacity.append(planck_opac)
    return planck_opacity

plt.plot(T,planck_mean_opacity(T))

"""
# Generate the symbolic variables
H_rim, X_rim, h_rim, Kb, T_rim, R_rim, mu, mp,GG, M, L ,ss = smp.symbols('H_rim X_ri h_rim Kb T_rim R_rim mu mp GG M L ss')

h_rim = smp.sqrt(Kb*T_rim* R_rim**3 / (mu*mp*GG, M))

H_rim = X_rim * h_rim
ABC= (L/(4* smp.pi * T_rim* ss))**0.5 *(1 + H_rim/ R_rim)**0.5 -R_rim
X_rim, Kb, T_rim, mu, mp,GG, M, L ,ss = [2,kb,1500,mu,mP,GG,0.5*ms,0.7*ls,ss]

smp.solve(ABC,R_rim)

h_rim_f = smp.lambdify([kb,T_rim, mu,mp,GG, M],h_rim)
H_rim_f = smp.lambdify([X_rim],H_rim)
R_rim_f = smp.lambdify([L,T_rim, ss])



h_rim_ans = h_rim_f(kb,T_rim, mu,mp,GG, M)
H_rim_ans = H_rim_f(X_rim)
R_rim_ans = R_rim_f(L,T_rim, ss)
"""