import numpy as np
import sympy as smp
import matplotlib.pyplot as plt


# ------------------------------------------------- ------------------- 
# Make the grid 
# Create a 3D grid and write it to a file called amr_grid.inp. 
# ------------------------------------------------- -------------------
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

# 
# grid parameters 

nr         =  2560          # Number of radial grids. 
ntheta     =  257        # Number of grids in elevation.
nz         =  2560          # Number of elevation grids. 

rin        =  0.1 * AU     #Sublimation radius  # inner rim of disk. Density is assumed to be 0 inside [cm] 
rout       = 4* AU      # Outer edge of disk. Density is set to 0 outside this [cm] 
theta_min = np.pi/2-0.18;           #Polar boundaries in radian
theta_max = np.pi/2+0.18;

#Stellar Parameter
mstar =0.5*ms
rstar =2.5*rs
tstar = 4000
lstar = 0.7*ls

#Disk parameter
t_sublimation = 1500













def func(H_rim, x_dummy):
    R_rim = (lstar/(4*np.pi*t_sublimation* ss))**0.5 *(1+ H_rim/x_dummy)**0.5
    return R_rim

Z=np.linspace(-30,30,200)
R_rim = np.zeros(x.shape)
w = 1
d = 10 #error cretria
track_d = []

while d > 0.001:
    track_d.append(d)
    temp = R_rim
    R_rim = w * R_rim + (1-w) * func(const, H_rim,R_rim, x)
    d = np.max(np.abs(temp-R_rim))

y=func(x, R_rim,A,B,C,D)
plt.plot(x,y)
plt.show()
# look at the convergence
plt.plot(track_d)
plt.show()