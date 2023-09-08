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
rout = 0.2*AU
ri = np.linspace(rin,rout,nr+1) 
R   = 0.5 * ( ri[0:nr] + ri[1:nr+1] )             # Take the in-"between values
R = np.array(R)

#sigma[g/cm^2] 
def Sigma(R):
    return sigma0 * (R/AU)**beta

sigma_trial = Sigma(R)
sigma_trial = np.array(sigma_trial)

#def Disk(X_cg, impinge_angle, H_cg, h_cg, Ti, R,  sigma, psi_i , psi_s):
def Disk(vars,R,  sigma, psi_i , psi_s):
    X_cg, H_cg, h_cg, Ti,impinge_angle   = vars
    #Define X_cg (Appendix A2)
    #Integral of exp(-z^2 /2H_cvg^2)dz = (pi/2)^0.5 x H_ccg (1- erf(1/2^2))
    #1 -erf(1/2^0.5) = 1- erf(X_cg/ 2^0.5) = 2 impingement_angle(X_cg)/ (sigma x Kp_stellar)
    eqn1 = np.abs(1-math.erf(X_cg/np.sqrt(2)) -2* impinge_angle/(sigma* kp_stellar))
    #Define H_cg
    eqn2 = np.abs(X_cg*h_cg -H_cg)
    #Define h_cg
    eqn3 = np.abs((Ti/t_virial)**0.5 * (R/rstar)**0.5 -h_cg/R)
    #Define Ti
    eqn4 = np.abs((impinge_angle*psi_s/psi_i)**0.25 * (rstar/R)**0.5 *tstar -Ti)
    #Define Impinge angle
    eqn5 = np.abs((0.4*rstar/R) +(gamma-1)*H_cg/R -impinge_angle)
    return [eqn1, eqn2,eqn3,eqn4, eqn5]

psi_i, psi_s = [1,1]
X_cg_trial = np.zeros(nr)
impinge_angle_trial = np.zeros(nr)
H_cg_trial = np.zeros(nr)
h_cg_trial = np.zeros(nr)
Ti_trial= np.zeros(nr)
for i,( r_i,s_i), in enumerate(zip(R,sigma_trial)):
    X_cg_trial[i], H_cg_trial[i], h_cg_trial[i], Ti_trial[i],impinge_angle_trial[i] =fsolve(Disk,(X_rim, 0.15, H_rim, h_rim, T_rim), args = (r_i, s_i, psi_i, psi_s))



fig, ax = plt.subplots()
ax.plot(R/AU, H_cg_trial/h_cg_trial)
ax.set_xlabel("Radius")
ax.set_ylabel("Surface_Heightttttttttt")
ax.set_title("Surface_Height vs. Radius")
plt.show()

fig, ax = plt.subplots()
ax.plot(R, Ti_trial)
ax.set_xlabel("Radius")
ax.set_ylabel("Temperature")
ax.set_title("Temperature vs. Radius")
plt.show()   
    
    

"""     
    diff = np.abs(X_cg_trial[i] *h_cg_trial  -H_cg_trial)
    print ("Hi      ",  diff)

    while diff>0.0000000000001:
        X_guess = (X_rim, 0.15, H_rim, h_rim, T_rim)*np.random.uniform(0.0001,1)
        print(X_guess)
        X_cg_trial[i], H_cg_trial[i], h_cg_trial[i], Ti_trial[i],impinge_angle_trial[i] =fsolve(Disk,X_guess, args = (r_i, s_i, psi_i, psi_s))
        diff = np.abs(X_cg_trial[i] *h_cg_trial  -H_cg_trial)
        #print(diff)

"""
    
    
    
    

#X_cg_trial, H_cg_trial, h_cg_trial, Ti_trial,impinge_angle_trial =fsolve(Disk,(X_rim, 0.15, H_rim, h_rim, T_rim), args = (R_rim, sigma_rim, psi_i, psi_s))
#X_cg_trial2, H_cg_trial2, h_cg_trial2, Ti_trial2,impinge_angle_trial2 = fsolve(Disk,(X_cg_trial, H_cg_trial, h_cg_trial, Ti_trial,impinge_angle_trial ),args = (R_rim, sigma_rim, psi_i, psi_s))


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

#Solving for interior temperature with constant phi_i and phi_s
def temp_interior(T_iter, const):
    #const_temp = (a_i*p_s/2/p_i)**0.25 *(rstar/r_i)**0.5 *tstar
    return np.abs(const  - T_iter)

def pressure_height(hcg_iter,const):
    #const_hcg = (T_i/t_virial)**0.5 * (r_i/rstar)**0.5*r_i
    return np.abs(const - hcg_iter)

def surface_height(Hcg_iter, const):
    #const_Hcg = X_i * h_i
    return np.abs(const - Hcg_iter)


def chi_disk(X_iter, const):
    #const_X = 2*a_i/kp_stellar/s_i
    #const_X = 2*a_i/kp_stellar/s_i

    return np.abs(1- math.erf(X_iter/(np.sqrt(2))) - const)



R_rim =0.47*AU
#Definig grid system in radial direction
rin = R_rim  # Sublimation radius  # inner rim of disk. Density is assumed to be 0 inside [cm] 
rout = 100*AU
#rout = 2.7 * (10 ** 2) * AU  # Outer edge of disk. Density is set to 0 outside this [cm] 
ri = np.linspace(rin,rout,nr+1) 
R   = 0.5 * ( ri[0:nr] + ri[1:nr+1] )             # Take the in-"between values
R = np.array(R)



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
    for i, (r_i, H_i, h_i, T_i,X_i, a_i, s_i) in enumerate(zip(R, H_cg_trial, h_cg_trial, Ti_trial, X_cg_trial,impinge_angle_trial, sigma_trial)):
       #print(i)
        #Simulation -=1
        #print(Simulation)
       #print("Hello\n\n\n\n")
        #print(r_i, H_i, h_i, T_i,X_i, a_i, s_i )
       #print("Help\n\n\n\n")
        const_alpha = (0.4*rstar + (gamma- 1)*H_i)*(1/r_i)
        #print(const_alpha)
        impinge_angle[i]= fsolve(alpha_angle,a_i *0.95 , args=const_alpha)
        #I have multiplied argument with 0.95 as it would be different from measured value
        #and potentially produce new answer that converges
        
        
        #print(impinge_angle)
       #print("Help1\n\n")
        #print(1)
        const_temp = (a_i*psi_s/2/psi_i)**0.25 *(rstar/r_i)**0.5 *tstar
        #print("Const temperaureate is ",const_temp)
        print(T_i)
        Ti[i]= fsolve(temp_interior,T_i*0.95, args=const_temp)
        print(T_i, "Afterrrrr")
        #print(2)
       #print("Help2\n\n")

        const_hcg = (T_i/t_virial)**0.5 * (r_i/rstar)**0.5 *r_i
        h_cg[i] = fsolve(pressure_height,h_i*0.95, args = const_hcg)
        #print(3)
       #print("Help3\n\n")

        const_Hcg = X_i * h_i
        #const_X = 2*a_i/kp_stellar/ 

        H_cg[i] = fsolve(surface_height,H_i*0.95, args =const_Hcg )
        #print(4)
       #print("Help4\n\n\n\n")
        const_X = 2*a_i/kp_stellar/s_i
        X_cg[i] = fsolve(chi_disk, X_i*0.95, args = const_X)
       #print("Hiiiiiiiiiiiiiiiiiiiiiiiiii")
        
       #print("Hello, I want to check this :  ", X_cg_trial[0], X_cg[0])
        
    #I dont know why but this doesnt over lay is list
    X_cg_trial = X_cg.copy()
        #print(X_cg_trial[i])
    impinge_angle_trial= impinge_angle.copy()
    print(H_cg[0])
    H_cg_trial = H_cg.copy()
    print(H_cg_trial[0])

    h_cg_trial = h_cg.copy()
    Ti_trial =  Ti.copy()
   #print("I am this hot: ",Ti_trial[0] )
        #print(H_cg_trial)
   #print("Hi  I am over here", z)
#   fig, ax = plt.subplots()
#   ax.plot(R, H_cg)
#    ax.set_xlabel("Radius")
#    ax.set_ylabel("Surface_Height")
#    ax.set_title("Surface_Height vs. Radius")
#    plt.show()
    #print(H_cg_trial)    
#print(Ti_trial)+
#print("NOooooooooooooooooooooooooooooo")
#print(Ti)
    

        

"""
for i,( r_i,s_i), in enumerate(zip(R,sigma_trial)):
    X_cg[i], H_cg[i], h_cg[i], Ti[i],impinge_angle[i] =fsolve(Disk,(X_rim, 0.15, H_rim, h_rim, T_rim), args = (r_i, s_i, psi_i, psi_s))


#This doesnt update after 1st simulation
for z in range(Simulation):
    for i, (r_i, H_i, h_i, T_i,X_i, a_i, s_i) in enumerate(zip(R, H_cg, h_cg, Ti, X_cg,impinge_angle, sigma)):
       #print(i)
        #Simulation -=1
        #print(Simulation)
       #print("Hello\n\n\n\n")
        #print(r_i, H_i, h_i, T_i,X_i, a_i, s_i )
       #print("Help\n\n\n\n")
        const_X = 2*a_i/kp_stellar/s_i
        X_cg[i] = fsolve(chi_disk, X_i*0.95, args = const_X)
        diff = np.abs(X_cg[i] - H_i/h_i)
        
        print("The difference is :", diff)
#        
#        while diff>0.0000000000001:
#            X_cg[i] = fsolve(chi_disk, X_i*0.95, args = const_X)
#            diff = np.abs(X_cg[i] - H_i/h_i)
#            #print(diff)
            
            
        const_alpha = (0.4*rstar + (gamma- 1)*H_i)*(1/r_i)
        #print(const_alpha)
        impinge_angle[i]= fsolve(alpha_angle,a_i *0.95 , args=const_alpha)
        #I have multiplied argument with 0.95 as it would be different from measured value
        #and potentially produce new answer that converges
        
        
        #print(impinge_angle)
       #print("Help1\n\n")
        #print(1)
        const_temp = (a_i*psi_s/2/psi_i)**0.25 *(rstar/r_i)**0.5 *tstar
       #print("Const temperaureate is ",const_temp)
        Ti[i]= fsolve(temp_interior,T_i*0.95, args=const_temp)
        #print(2)
       #print("Help2\n\n")

        const_hcg = (T_i/t_virial)**0.5 * (r_i/rstar)**0.5 *r_i
        h_cg[i] = fsolve(pressure_height,h_i*0.95, args = const_hcg)
        #print(3)
       #print("Help3\n\n")

        const_Hcg = X_i * h_i
        #const_X = 2*a_i/kp_stellar/ 

        H_cg[i] = fsolve(surface_height,H_i*0.95, args =const_Hcg )
        #print(4)
       #print("Help4\n\n\n\n")

       #print("Hiiiiiiiiiiiiiiiiiiiiiiiiii")
        
       #print("Hello, I want to check this :  ", X_cg_trial[0], X_cg[0])

"""
fig, ax = plt.subplots()
ax.plot(R/AU, H_cg/h_cg)
ax.set_xlabel("Radius")
ax.set_ylabel("Surface_Heightttttttttt")
ax.set_title("Surface_Height vs. Radius")
plt.show()

fig, ax = plt.subplots()
ax.plot(R, Ti)
ax.set_xlabel("Radius")
ax.set_ylabel("Temperature")
ax.set_title("Temperature vs. Radius")
plt.show()   
    
    
    