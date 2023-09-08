import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import fsolve
import itertools
import math
from scipy.optimize import fsolve
from scipy.special import erfinv


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


gray_const= 46.77977670461847




def rim_function(variable):
    x1, x2, x3, x4, x5 = variable
    #R_rim, h_rim, H_rim, X_rim, sigma_rim
    f =[(gray_const*(lstar/(np.pi* 4*T_rim**4 *ss))**0.5 *(1+ x3/x1)**0.5 -x1),
           (((kb*T_rim*x1**3)/(mu*mp*GG*mstar))**0.5 -x2),
           (x2*x4 - x3),
           (4*x5*kp_stellar*(1-math.erf(x4)) - 1),
           (sigma0*(x1/AU)**beta - x5)]
    return f


pi = np.pi
def jacobian_rim(variable):
    x1, x2, x3, x4, x5 = variable
    #R_rim, h_rim, H_rim, X_rim, sigma_rim

    J = [[-0.25*x3*gray_const*(lstar/(T_rim**4*ss))**0.5/(pi**0.5*x1**2*(x3/x1 + 1)**0.5) - 1,       0,          gray_const*0.25*(lstar/(T_rim**4*ss))**0.5/(pi**0.5*x1*(x3/x1 + 1)**0.5),       0,     0],
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















def surface_density(R):
    return sigma0*(R/AU)**beta

sigma = surface_density(R)



def viscous_temp(R):
    yrs = 31557600 #s
    M_acc = 10**-8*ms/yrs

    return (3*GG*mstar*M_acc/(8*np.pi* ss* R**3))**0.25
    #return 107*(R/AU)**-0.75 * (mstar/2.5/ms)**0.25 *(M_acc/(10**-8 *ms/yrs))**0.25

T_visc = viscous_temp(R)

R_rim = 0.54*AU
H_rim= 0.3*R_rim
T_rim = 1500
beta = 1
betaaa = 1
def Rim_heating(R):
    return 1.2*T_rim*np.exp(-1*(1-0.5**(1/(4+betaaa)))  *(R-R_rim))
RimHeating = Rim_heating(R)


def pressure_height(Ti, R):
    return ((Ti/t_virial)**0.5 * (R/rstar)**0.5 * R)

def scale_height(X_cg, h_cg):
    return (X_cg *h_cg)


def chi_cg(sigma):
    return (erfinv(1-1/(4*sigma*kp_stellar)))
    #return 1-  (4*sigmaRim*kp_stellar)*(1-math.erf(X_rim_iter))
 
def surface_temp(R):
    return (0.25)**(1/(4+beta)) * (rstar/R)**(2/(4+beta)) * tstar


Ts = surface_temp(R)
T_effective = (RimHeating**4 + T_visc**4 +Ts**4 )**0.25
h_cg = pressure_height(T_effective, R)
X_cg = chi_cg(sigma)
H_cg = scale_height(X_cg, h_cg)

def alpha_angle(H_cg,R):
    return (0.4*rstar/R   +(gamma - 1)*H_cg/R)

impinge_angle= alpha_angle(H_cg, R)

psi_i, psi_s = 1, 1
def irr_disk_temp(impinge_angle, R):
    return ((impinge_angle*psi_s/psi_i)**0.25 * (rstar/R)**0.5 *tstar)

Ti = irr_disk_temp(impinge_angle, R)

#T_effective = (RimHeating**4 + T_visc**4 +Ts**4  +Ti**4)**0.25
T_effective = (RimHeating**4 + T_visc**4  +Ti**4)**0.25

h_cg = pressure_height(T_effective, R)
X_cg = chi_cg(sigma)
H_cg = scale_height(X_cg, h_cg)

Rflare = np.where (R <H_cg*R_rim/H_rim )#and H_cg*R_rim/R/H_rim >1)
Rshadow = np.where (R >H_cg*R_rim/H_rim)

Ti[Rshadow] = 0
T_effective = (RimHeating**4 + T_visc**4  +Ti**4)**0.25
h_cg = pressure_height(T_effective, R)
X_cg = chi_cg(sigma)
H_cg = scale_height(X_cg, h_cg)




h_cg[0] = h_rim
H_cg[0] = H_rim
X_cg[0] = X_rim
T_effective[0] = T_rim

fig, ax = plt.subplots()
ax.semilogx(R/AU, Ti)
Rflare = np.where (R <H_cg*R_rim/H_rim)

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
ax.semilogx(R, RimHeating)

ax.set_xlabel("Radius")
ax.set_ylabel("RimHeating")
ax.set_title("RimHeating vs. Radius")
ax.set_ylim([0, 0.5])

plt.show()




fig, ax = plt.subplots()
ax.semilogx(R, T_effective)

ax.set_xlabel("Radius")
ax.set_ylabel("Effective_temperature")
ax.set_title("Effective_temperature vs. Radius")

plt.show()



fig, ax = plt.subplots()
ax.semilogx(R/AU, H_cg/R)

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


