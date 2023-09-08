import numpy as np
import sympy as smp
from sympy import diff
from sympy import erf


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


vars = smp.symbols('R_rim h_rim H_rim X_rim sigma_rim') #Defining variable to calculate partial differentation
f = (['pow((lstar/(4*pi*(pow(T_rim,4))*ss)), 0.5)*pow((1+ H_rim/R_rim),0.5) -R_rim', #R_rim 
                     'pow((kb*T_rim*pow(R_rim,3)/(mu*mp*GG*mstar)),0.5) - h_rim', #h_rim
                     'X_rim * h_rim -H_rim', #H_rim
                     '(1-erf(X_rim))*(4*sigma_rim*kp_stellar) -1', #X_rim
                     'sigma0*pow((R_rim/AU),beta) -sigma_rim'])

J = smp.zeros(len(f),len(vars)) # Initialise Jacobian matrix
for i, fi in enumerate(f):
    for j, s in enumerate(vars):
        print(s)
        J[i,j] = smp.diff(fi, s)
        print(J[i,j])
        
def jacobian_example(xy):
    x, y = xy
    return [[1, 2],
            [2*x, 8*y]]

vars= ('x1', 'x2', 'x3', 'x4' , 'x5')
x1 ,x2, x3,x4,x5 = vars
f = ([((lstar/(4*np.pi*((T_rim)**4)*ss))**0.5)*((1+ x2/x1)**0.5) -x1, #R_rim 
                     ((kb*T_rim*(x1**3)/(mu*mp*GG*mstar))**0.5) - x3, #h_rim
                     x4 * x3 -x2, #H_rim
                     (1-erf(x4))*(4*x5*kp_stellar) -1, #X_rim
                     sigma0*((x1/AU)**beta) -x5])

J = smp.zeros(len(f),len(vars)) # Initialise Jacobian matrix
for i, fi in enumerate(f):
    for j, s in enumerate(vars):
        print(s)
        J[i,j] = smp.diff(fi, s)
        print(J[i,j])