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
sigmao = 2*10 **3  # Surface density at 1 AU  Σ0 [g/cm^2] 
kp_stellar = 400.0  # Monochromatic opacity = 400cm^2 g^-1

nu = 5.879 * (10 ** 10) * tstar  #Wien law
t_virial = GG * mstar * mu * mp / (kb * rstar)  # Virial Temp
gamma = 2.0
beta = -1.5


pi = np.pi
#vars = smp.symbols('R_rim h_rim H_rim X_rim sigma_rim') # Define x and y variables
#R_rim, h_rim ,H_rim,X_rim, sigma_rim =vars # Define x and y variables


#This should return function that describe Rim property
#def function_Rim(R_rim, h_rim ,H_rim,X_rim, sigma_rim):
def function_Rim(vars):
    #f(R_rim, h_rim, H_rim, X_rim, sigma_rim) =0
    vars = smp.symbols('R_rim h_rim H_rim X_rim sigma_rim') #Defining variable to calculate partial differentation
    f = (['pow((lstar/(4*pi*(pow(T_rim,4))*ss)), 0.5)*pow((1+ H_rim/R_rim),0.5) -R_rim', #R_rim 
                         'pow((kb*T_rim*pow(R_rim,3)/(mu*mp*GG*mstar)),0.5) - h_rim', #h_rim
                         'X_rim * h_rim -H_rim', #H_rim
                         '(1-erf(X_rim))*(4*sigma_rim*kp_stellar) -1', #X_rim
                         'sigmao*pow((R_rim/AU),beta) -sigma_rim']) #sigma_rim
    return f

    #vars = smp.symbols('R_rim h_rim H_rim X_rim sigma_rim') # Define x and y variables
"""
    f = ([(lstar/(4*pi * T_rim**4 *ss))**0.5 * (1 + H_rim/R_rim)**0.5 - R_rim,
                     (kb*T_rim*R_rim**3/(mu*mp*GG*mstar))**0.5 - h_rim,
                     X_rim * h_rim -H_rim,
                     ((1-erf(X_rim))*(4*sigma_rim*kp_stellar) -1),
                     sigmao*(R_rim/AU)**beta -sigma_rim]) # Define function
    return f
"""
#This should return function that describe disk_cg property
def function_Disk(vars, R,sigma, psi_i, psi_s):
    #X_cg, impinge_angle, H_cg, h_cg, Ti =vars
    vars = smp.symbols('Ti impinge_angle h_cg H_cg X_cg') # Define x and y variables
    f = (['pow((impinge_angle*psi_s/psi_i),0.25) * pow((rstar/R), 0.5) *tstar -Ti',
          '(0.4*rstar/R) +(gamma-1)*H_cg/R -impinge_angle',
          'pow((Ti/t_virial),0.5) * pow((R/rstar),0.5) -h_cg/R',
          'X_cg*h_cg -H_cg',
          '1 - smp.erf(X_cg/(pow(2,0.5)))- (2*impinge_angle/sigma/kp_stellar)'])
    return f

"""


#This define function should call a function_rim or function_disk 
#and calculate the Jacobian matrix of the function
def Jacobian(vars,function):
    R_rim, h_rim ,H_rim,X_rim, sigma_rim =vars # Define x and y variables
    #R_rim, h_rim ,H_rim,X_rim, sigma_rim =vars # Define x and y variables
    vars = smp.symbols('R_rim h_rim H_rim X_rim sigma_rim') # Define x and y variables

    f = function(vars)
    J = smp.zeros(len(f),len(f)) # Initialise Jacobian matrix

    # Fill Jacobian matrix with entries
    for i, fi in enumerate(f):
        for j, s in enumerate(vars):
            J[i,j] = smp.diff(fi, s)
            print(J[i,j])
    return J
"""

def Jacobian(vars):
    #|df1/dx1  df1/dx2   ....... df1/dxn|
    #|df2/dx1  df2/dx2 ..........df2/dxn|
    #|df3/dx1  df3/dx2   ....... df3/dxn|
    #|dfn/dx1  dfn/dx2 ...........dfn/dxn/
    
    #where f represent function defining R_, H_, h_, T_, sigma_ X-
    #x represent partial variable for derivative
    
    
    #R_rim, h_rim ,H_rim,X_rim, sigma_rim =vars # Define x and y variables
    #vars = smp.symbols('R_rim h_rim H_rim X_rim sigma_rim') # Define x and y variables
    
    #vars = smp.symbols('Ti impinge_angle h_cg H_cg X_cg') # Define x and y variables
    #vars = smp.symbols('Ti impinge_angle h_cg H_cg X_cg') # Define x and y variables

    f = function_Rim(vars)
    J = smp.zeros(len(f),len(vars)) # Initialise Jacobian matrix

    # Fill Jacobian matrix with entries
    for i, fi in enumerate(f):
        for j, s in enumerate(vars):
            print(s)
            J[i,j] = smp.diff(fi, s)
            print(J[i,j])
    return J



def iterative_newton(function, x_inital, jacobian):
    max_iter = 500
    epsilon = 10**-6 # The modulus difference between |x[k+1] -x[k]|

    x_last = x_inital

    for k in range(max_iter):
        # Solve J(xn)*( xn+1 - xn ) = -F(xn):
        #Ja = lambdify(vars,jacobian, 'numpy')   #need to change it into numpy
        #J = Ja(x_last)
        
        #Fa = lambdify(vars, function, 'numpy')
        #F = Fa(x_last)
        J = np.array(jacobian(x_last))
        F = np.array(function(x_last))
        
        #Solving the value for delta x 
        #J(x[k])*delta x[k] = -f(x[k])
        #delta x[k] = x[k+1] - x[k]
        #x[k+1] = delta_X[k] + x[k]
        delta_x = np.linalg.solve( J, -F )
        #defining new value for x[k+1]
        x_last = x_last + delta_x

        # Stop condition:
            
        #np.linalg.norm(delta_x) calcualte the difference between x_old (x[k]) and x_new (x[k+1])
        if np.linalg.norm(delta_x) < epsilon:
            print('convergence!, nre iter:', k )
            break

    else: # only if the for loop end 'naturally'
        print('not converged')

    return x_last


#x_sol = iterative_newton(function_Rim, [0.4904092456436323*AU, 0.39336608131611195*AU,0.13113940101171412*AU,4.512601124660966*AU,1650.5889740859373], Jacobian)
#print('solution exercice:', x_sol )
X_guess_rim = [0.4904092456436323*AU, 0.39336608131611195*AU,0.13113940101171412*AU,4.512601124660966*AU,1650.5889740859373]

#Partial derivative with respective to variable for Rim
vars = smp.symbols('R_rim h_rim H_rim X_rim sigma_rim')
   
#vars = smp.symbols('Ti impinge_angle h_cg H_cg X_cg') # Define x and y variables
x_Rim_sol = iterative_newton(function_Rim,X_guess_rim, Jacobian)
R_rim,h_rim, H_rim, X_rim, sigma_rim = x_Rim_sol[0],x_Rim_sol[1],x_Rim_sol[2],x_Rim_sol[3],x_Rim_sol[4]



#Definig grid system in radial direction
rin = R_rim  # Sublimation radius  # inner rim of disk. Density is assumed to be 0 inside [cm] 
rout = 100*AU
#rout = 2.7 * (10 ** 2) * AU  # Outer edge of disk. Density is set to 0 outside this [cm] 
ri = np.linspace(rin,rout,nr+1) 
#R = np.array([0.5 * (ri[i:] + ri[:i]) for i in range(len(ri))])  # Take the in-"between values
#R   = np. array(0.5 * ( ri[::-1] + ri) ) 
R   = 0.5 * ( ri[0:nr] + ri[1:nr+1] )             # Take the in-"between values
R = np.array(R)



X_cg = np.zeros(nr)
impinge_angle = np.zeros(nr)
H_cg = np.zeros(nr)
h_cg = np.zeros(nr)
Ti= np.zeros(nr)
#psi_i =np.array(list(itertools.repeat(1, nr)))
#psi_s = np.array(list(itertools.repeat(1, nr)))
psi_i = 1
psi_s = 1


#Partial derivative with respective to variable for Disk function
vars = smp.symbols('Ti impinge_angle h_cg H_cg X_cg') 
x_Disk_sol = iterative_newton(function_Disk,X_guess_rim, Jacobian)
Ti, impinge_angle ,h_cg ,H_cg, X_cg = x_Rim_sol[0],x_Rim_sol[1],x_Rim_sol[2],x_Rim_sol[3],x_Rim_sol[4]

def Sigma(R):
    return sigmao * (R/AU)**beta

sigma = Sigma(R)
sigma = np.array(sigma)
X_guess_disk = [T_rim, 0.2,0.001*AU, 0.1*AU, 1,R_rim, sigma_rim]

for i,( r_i,s_i), in enumerate(zip(R,sigma)):
    x_Disk_sol = iterative_newton(function_Disk,X_guess_disk, Jacobian)
    Ti[i], impinge_angle[i] ,h_cg[i] ,H_cg[i], X_cg[i]  =x_Disk_sol[0],x_Disk_sol[1],x_Disk_sol[2],x_Disk_sol[3],x_Disk_sol[4]
    X_guess_disk = [Ti[i], impinge_angle[i] ,h_cg[i] ,H_cg[i], X_cg[i],r_i, s_i]


zi = np.linspace(0.0, H_cg, nz +1)
z_elevation   = 0.5 * (zi[0:nz] + zi[1:nz+1])          # Take the in-between values
z_elevation = np.array(z_elevation)
#z_elevation   = np.array([0.5 * ( zi[i] + zi[i+1])   for i in range(len(zi) - 1)])  # Take the in-"between values       # Take the in-between values



#Creating simulation to run same program such that solution converges with each simulation


Simulation =100
for i in range (Simulation):
    ri = np.linspace(rin,rout,nr+1) 
    R   = 0.5 * ( ri[0:nr] + ri[1:nr+1] )             # Take the in-"between values
    R = np.array(R)
    
    #Partial derivative with respective to variable for Rim
    vars = smp.symbols('R_rim h_rim H_rim X_rim sigma_rim')    
    x_Disk_sol = iterative_newton(function_Disk,X_guess_rim, Jacobian)
    Ti, impinge_angle ,h_cg ,H_cg, X_cg = x_Rim_sol[0],x_Rim_sol[1],x_Rim_sol[2],x_Rim_sol[3],x_Rim_sol[4]
    
    
    #Partial derivative with respective to variable for Disk function
    vars = smp.symbols('Ti impinge_angle h_cg H_cg X_cg') 
    for i,( r_i,s_i), in enumerate(zip(R,sigma)):
        X_guess_disk = [Ti[i], impinge_angle[i] ,h_cg[i] ,H_cg[i], X_cg[i],r_i, s_i]
        x_Disk_sol = iterative_newton(function_Disk,X_guess_disk, Jacobian)
        Ti[i], impinge_angle[i] ,h_cg[i] ,H_cg[i], X_cg[i]  =x_Disk_sol[0],x_Disk_sol[1],x_Disk_sol[2],x_Disk_sol[3],x_Disk_sol[4]


    zi = np.linspace(0.0, H_cg, nz +1)
    z_elevation   = 0.5 * (zi[0:nz] + zi[1:nz+1])          # Take the in-between values
    z_elevation = np.array(z_elevation)
    #z_elevation   