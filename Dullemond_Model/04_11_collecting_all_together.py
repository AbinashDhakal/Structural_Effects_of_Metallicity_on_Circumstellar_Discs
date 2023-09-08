import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfinv
import pandas as pd

#Natural constant
GG =  6.67430e-8  # U# gravitational constant  [dyn cm^2 g^-2]
kb =  1.380649e-16  # Boltzmann constant defined [erg K^-1]
hp =  6.62607015e-27  # Planck constant [erg s] defined
clight = 2.99792458e10  # light speed [cm s^-1] defined 
NA = 6.02214076e23  # Avogadro constant defined 
mp = 1.672621898e-24  # proton mass [g]; not 1.0/NA anymore. 
ss = 5.6703e-5  # Stefan-Boltzmann const  [erg/cm^2/K^4/s]
gasR = 	8.3143435* 10**7 #erg K-1 mole-1
a =  ss/0.25/clight
#Astronomical constant
AU = 1.495978707e13  # Distance between Earth and Sun [cm] defined 
pc = 3.085677581e18  # Parsec [ cm] 
ms = 1.9884e33  # Sun mass [g] 
ts = 5.777e3  # effective temperature of the sun [K] 
ls = 3.828e33  # solar intensity [erg/s] defined by IAU 
rs = 6.96e10  # solar radius [cm] 

#Gas constant
mu = 2.353  # Average molecular weight  (H2 + He + metals)
adia_gamma = 1.42 # adiabatic index
Pr = 1 #Prandtl number
cV = kb/(mu*mp*(adia_gamma-1)) #specific const volume 
cP = adia_gamma *cV #specific const pressure

# grid parameters 
nr = 256  # Number of radial grids. 
ntheta = 257  # Number of grids in elevation.
nz = 256 # Number of elevation grids. 

#Stellar Parameter
mstar = 2.5 * ms
rstar = 2.5*rs
tstar = 5000
lstar = 47* ls

#Disk parameter
T_rim= t_sublimation = 1500.0  # Temperature at which dust start to condense
sigma0 = 2*10**3
psi_i, psi_s = 1, 1 # Absorption and emmission factor 

#Inner gas hole parameter
Metallicity =-0.3
T_stellar =tstar


nu = 5.879 * (10 ** 10) * tstar  #Wien law
t_virial = GG * mstar * mu * mp / (kb * rstar)  # Virial Temp
gamma = 2.0  #Accounting circumstellar rim 
beta = -1.5  #density decay factor 
gray_const =1
visc_const = 10**-2 #Shakura & Sunyaev turblent viscosity constant

#Astro-silicate paramater
A_sil = 10**(-25.10)   #A_sil = 10^-25.10 cm^2.5
a_min = 5*10**-7
a_max = 2.5 * 10**-5
sil_density = 3.4 # g cm^-3
#dust to gas mass ratios ζgrap = 0.0025 and ζsil = 0.0034
#(Draine & Lee 1984). 
num_sil_density = sil_density/mp *A_sil*(a_max**-2.5 - a_min**-2.5)/(-3.5)
kv = (3/4)*(1.5*10**4)*(10**(-25.10)) * ((a_max**-2.5 - a_min**-2.5)/(-3.5))
ni_const = (3/4)*(A_sil)*((a_max**-2.5 - a_min**-2.5)/(-3.5)) /(10**-6)
fd2g = 10**-2    # The gas: dust is assumed be 100 : 1.



# Set maximum iteration
maxIter = 500





################################*Defning different paramater *#################

#Surface density profile
def surface_density(R):
    return sigma0*(R/AU)**beta

def pressure_height(T, R):
    return ((T/t_virial)**0.5 * (R/rstar)**0.5 * R)

def scale_height(X_cg, h_cg):
    return (X_cg *h_cg)

def chi_cg(sigma,kp_stellar):
    return (erfinv(1-1/(4*sigma*kp_stellar)))
    #return 1-  (4*sigmaRim*kp_stellar)*(1-math.erf(X_rim_iter))


#fraction of stellar flux that will be absorbed by the interior
def flux_absorbed_fraction(sigma, kp_T): #psi_s
    if (sigma*kp_T) < 1:
        return kp_T
    else:
        return 1
         
#The possiblity that the disk interior is not fully optically thick to its own emission
def self_absorbed_fraction(sigma,kp_T): #psi_i
    if (sigma*kp_T) < 1:
        return sigma*kp_T
    else:
        return 1

#Defining interior heating temperature
def irr_disk_temp(impinge_angle, R,psi_i, psi_s):
    return ((impinge_angle*psi_s/psi_i)**0.25 * (rstar/R)**0.5 *tstar)

#Definging surface tempeature
def irr_surface_Temp(R, kp_Ts):
    epsilon =kp_Ts/kp_stellar
    T_s = (1/epsilon)**0.25 * (rstar/2/R)**0.5 *tstar
    return T_s


#Direct stellar radiation impinges onto the disk at an angle
def alpha_angle(H_cg,R):
    return (0.4*rstar/R   +(gamma - 1)*H_cg/R)


def Rho(z_elevation, h_cg, sigma):
    return sigma*np.exp(- z_elevation**2 / 2/h_cg**2) / (np.sqrt(2*np.pi) *h_cg)

#Pressure [dyn/cm^2] = [g cm/s^2/cm^2] = [g/s^2/cm]
def GasPressure(rho, Tgas):
    return (kb*Tgas*rho/mu/mp)


from scipy.integrate import quad
def intergrand(R, Tgas,kp_T):
    return kp_T*(sigma0*(R/AU)**beta)/np.sqrt(2*np.pi)/((Tgas/t_virial)**0.5 * (R/rstar)**0.5 * R)

def Tau(R, Tgas,kp_T):
    return quad(intergrand,np.inf,R,args=(Tgas,kp_T))[0]

def closest(lst, K):
     lst = np.asarray(lst)
     idx = (np.abs(lst - K)).argmin()
     return lst[idx]

def planckOpacity(Metallicity,T_stellar,T_gas,P_gas,Rho_gas):
    df = pd.read_csv("Planck_data.csv")
    planck_opacity = []
    planck_opacity = np.array(planck_opacity)
    filtering_data = ( df[
        (df["Metallicity"] == (closest(df["Metallicity"], Metallicity))) 
        &(df["T_stellar"] == (closest(df["T_stellar"], T_stellar))) 
        & (df["T_gas"] == (closest(df["T_gas"],T_gas)))
        #& (df["Rho_gas"] == (closest(df["Rho_gas"],Rho_gas)))
        & (df["P_gas"] == (closest(df["P_gas"],  P_gas)))
       ])
    #print( closest(df["T_stellar"], T_stellar), closest(df["T_gas"],T_gas),(closest(df["Rho_gas"],Rho_gas)),(closest(df["P_gas"],  P_gas)))
    #print(filtering_data)
    filtering_data_array = filtering_data.to_numpy()
    #print(filtering_data_array)
    planck_opacity_data = filtering_data_array
    #print(planck_opacity_data[0][5])
    
    return planck_opacity_data[0][5]



#############################*Modelling Kp_stellar*############################

def reading_dust_data():
    df = pd.read_csv("Planck_dust_data.csv")
    T_data = np.array(list(df["T_rad"]))
    Qabs_data = np.array(list(df["<Qabs>/a"]))
    kp_t_array = np.zeros((2, len(Qabs_data)))
    kp_t_array[0] =T_data
    kp_t_array[1] =Qabs_data
    #change_array = T_kp_data.to_numpy()
    return T_data, Qabs_data

T_data, Qabs_data = reading_dust_data()

def Qabs_interp(T_guess):
    T_data, Qabs_data = reading_dust_data()
    logT = np.log10(T_data)
    logQ = np.log10(Qabs_data)
    logT_interp = np.log10(T_guess)
    return np.power(10.0, np.interp(logT_interp, logT, logQ))
kp_stellar = Qabs_interp(tstar)*(3/4)*(1/sil_density)*10**4



#################################*Modelling Inner hole*########################
def innerHole(Ri):
    Tguess = T_rim
    Tgas = 9200
    data = np.loadtxt('planck_data.csv', delimiter=',', skiprows=1)
    modified_data = np.loadtxt('planck_data.csv', delimiter=',', skiprows=1,  usecols=(1, 2,3,4,5))
    #print(Tgas, Tguess)
    while np.round(Tguess, 2) != np.round(Tgas,2):
        Tguess = Tgas
        s_gas = surface_density(Ri)
        h_gas =pressure_height(Tgas, Ri)
        Rho_gas = 10**-13
        #Rho_gas = Rho(0, h_gas, s_gas)
        P_gas = GasPressure(Rho_gas, Tgas)
        #print("printing pressute",P_gas,Rho_gas)
        # Vector we try to match with data
        
        #Defining Kp(Tgas, Trad)
        search_vect =  np.array([Metallicity, T_stellar, Tgas ,P_gas, Rho_gas])
        # normalize data in each column with the column mean
        norm_data = modified_data /np.mean(modified_data, axis=0)
        #print("printing normalised vecotr",norm_data)
        # normalize the search vector with the same values
        s_vec_n = search_vect/np.mean(modified_data , axis=0)
        idx = np.argmin(np.linalg.norm((s_vec_n - norm_data), axis=1))
        kp_TgasTrad = data[idx][5]
        #print(kp_TgasTrad)
        
        #Defining Kp(Tgas, Tgas)
        search_vect2 =  np.array([Metallicity, T_stellar, T_stellar ,P_gas, Rho_gas])
        # normalize data in each column with the column mean
        norm_data2 = modified_data /np.mean(modified_data, axis=0)
        #print("printing normalised vecotr",norm_data)
        # normalize the search vector with the same values
        s_vec_n2 = search_vect2/np.mean(modified_data , axis=0)
        idx2 = np.argmin(np.linalg.norm((s_vec_n2 - norm_data2), axis=1))
        kp_Tgas = data[idx2][5]
        #print(kp_Tgas)
        Tgas =(kp_Tgas/kp_TgasTrad)**(-0.25)*(rstar/Ri)**(0.5) *tstar
        #Tgas = ((rstar/Ri)**2 *np.exp(-Tauuu )*Rho_gas * kp_T* tstar**4)**0.25
        
    return Tgas

R_rim = 0.47*AU
H_rim = 0.074*AU
R_hole = np.linspace(rstar,R_rim,nr)
Rim_temp = innerHole(R_rim)
print(innerHole(R_rim))
Temp_innerhole = np.zeros(nr)
def innerHole_Temp(R_hole):
    for i, ri in enumerate(R_hole):
        Temp_innerhole[i] = innerHole(ri)
    return Temp_innerhole
"""
Temp_innerhole = innerHole_Temp(R_hole)

#Graphing inner hole temp
fig, ax = plt.subplots()
ax.semilogx(R_hole/AU, innerHole_Temp(R_hole))
ax.set_xlabel("Radius")
ax.set_ylabel("Interior_Temperature")
ax.set_title("Interior_Temperature vs. Radius")
plt.show()   """


################################*Modelling Puffed up rim*#####################


def Rim_radial(H_rim, Rrim_guess, gray_const):
    return gray_const**0.5*(lstar/(4*np.pi*T_rim**4 *ss))**0.5 *(1+ H_rim/Rrim_guess)**0.5

def converging_rim(gray_const):
    R_guess = 0.53*AU
    R_rim = 0.54*AU
    H_rim =0.1162*AU
    H_guess =0.11*R_rim
    iteration =0
    while np.round(H_guess/R_rim,2) != np.round(H_rim/R_rim,2):
        H_guess =        H_rim
        w = 0.9999999
        error = 1*AU
        while error >10e-5 *AU:
            temp_r = R_rim
            R_rim = (w*R_rim +(1-w)*Rim_radial(H_rim, R_rim,gray_const))
            error = np.max(np.abs(temp_r-R_rim))
            #print(R_rim/AU, temp_r/AU)
        iteration +=1
        #print("I am running iteration", iteration)
        #print("This is dust_const", gray_const)

        sigma_rim = surface_density(R_rim)
        h_rim = ((T_rim/t_virial)**0.5 * (R_rim/rstar)**0.5 * R_rim)
        X_rim = chi_cg(sigma_rim,kp_stellar)
        H_rim = scale_height(X_rim, h_rim)
    return R_rim, H_rim, h_rim, sigma_rim, X_rim, 


R_rim, H_rim, h_rim, sigma_rim, X_rim = converging_rim(gray_const)
print(R_rim/AU,H_rim/R_rim)   


#############################*Defining Radial limits*##########################
rin = R_rim # Sublimation radius  # inner rim of disk. Density is assumed to
rout = 100*AU
R =ri = np.linspace(rin,rout,nr) 
R = np.array(R)



################################*Defining heating source*######################
#Q_visc = Rho* vt*(r*dr_omega)**2
#Q_visc = (sigma/(2*pi*h))*(visc_const)*(np.sqrt(kb*T/mu/mh))*(9*G*Mstar/4/R**3)
def Q_visc(R):
    yrs = 31557600 #s
    M_acc = 10**-8*ms/yrs
    T_visc = (3*GG*mstar*M_acc/(8*np.pi* ss* R**3))**0.25
    #print(T_visc)
    for iteration in range(0, maxIter):
        for i in range(1,len(T_visc)-1):    
            sigma = surface_density(R[i])
            cs = np.sqrt(kb*T_visc[i]/mu/mp)
            Q_visc =(sigma/2/np.pi)*visc_const*cs*(GG*mstar/4/R[i]**3)
            T_visc[i] = (Q_visc/2/ss)**0.25
    T_visc = (3*GG*mstar*M_acc/(8*np.pi* ss* R**3))**0.25

    return T_visc

T_visc = Q_visc(R)

fig, ax = plt.subplots()
ax.semilogx(R/AU, T_visc)
ax.set_xlabel("Radius")
ax.set_ylabel("T_visc")
ax.set_title("T_visc vs. Radius")
plt.show() 




T_inital_approx = (0.5)*(rstar/R)**0.5 *tstar

#Q_cond = kt∇^2T
#kt : thermal conductivity
dr = (R[1]-R[0]) *AU  
def Q_cond(R,T,dr):
    #Defining a potential form of T
    
    T = (0.5)*(rstar/R)**0.5 *tstar
    Q = np.zeros(nr)
    for i in range(1,len(T)-1):
        #T[0] = 1500
        #print(R[i]/AU)
        #Calcuating given paramater for ith element with T[i] and R[i]
        #T and R has same array size len(R) = len(T)
        hp = (T[i]/t_virial)**0.5 *(R[i]/rstar)**0.5 *R[i] #pressure scale
        cs = np.sqrt(kb*adia_gamma*T[i]/(mu*mp)) #sound speed
        omega = np.sqrt(GG*mstar/R[i]**3) #Keplerian velocity
        #hp =cs/omega
        cV = kb/(mu*mp*(adia_gamma-1)) #specific const volume 
        vt = visc_const*cs**2/omega #viscosity
        sigma = surface_density(R[i]) #surface density 
        rho = sigma/(np.sqrt(2*np.pi* hp**2)) #Density 
        kt = rho*vt*adia_gamma*cV #Turblent thermal conductivity
        #Discretising laplace equation

        dT_dR_2 = (T[i+1] -2*T[i] +T[i-1])/(dr)**2
        dT_dR = (T[i+1] -T[i])/dr/2
        Q[i] =kt*(dT_dR_2 +dT_dR/R[i])
        #Q_cool = 2 *ss* T
        #Q_condct = kt * ΔT
        #Q_cool = Q_conduct
    return T

T_cond =Q_cond(R,T_inital_approx,dr)

fig, ax = plt.subplots()
ax.semilogx(R/AU, T_cond)
ax.set_xlabel("Radius")
ax.set_ylabel("T_cond")
ax.set_title("T_cond vs. Radius")
plt.show()   



alpha_guess = 0.2
#convering iteration such that the it iterate until the error is smaller than 0.00001
#THis only check for a radial point
def converging_alpha(alpha_guess, r_ii, psi_i, psi_s):
    a_i = 0.22
    while np.round(a_i, 5) != np.round(alpha_guess,5):
    #while alpha_guess !=a_i:
        alpha_guess = a_i
        T_i =irr_disk_temp(alpha_guess, r_ii,psi_i, psi_s)
        h_i= pressure_height(T_i,r_ii)
        s_i = surface_density(r_ii)
        X_i = chi_cg(s_i,kp_stellar)
        H_i = X_i*h_i
        a_i =alpha_angle(H_i, r_ii)
        
    return a_i
    
X_cg = np.zeros(nr)
impinge_angle = np.zeros(nr)
H_cg = np.zeros(nr)
h_cg = np.zeros(nr)
sigma = np.zeros(nr)
T_irr2 = np.zeros(nr)
   
def Q_irr(R, psi_i, psi_s,alpha_guess):
    Ti = (0.5)*(rstar/R)**0.5 *tstar
        
    for i, r_i in enumerate(R):
        impinge_angle_iter =converging_alpha(alpha_guess, r_i,psi_i, psi_s)
        alpha_guess = impinge_angle_iter
        sigma_iter = surface_density(r_i)
        w = 0.9999999
        T_s = 1200
        h_cg_iter = 0.3
        error = 10
        Ti_iter = irr_disk_temp(alpha_guess, r_i,psi_i, psi_s)
        #This check Ti multiple times such that error is less than 0.01
        while error >0.01:
            Temp = Ti_iter
            Rho_gas =Rho(0, h_cg_iter, sigma_iter)
            P_gas =GasPressure(Rho_gas, Temp)
            kp_Ts = Qabs_interp(T_s)*(3/4)*(1/sil_density)*10**4
            T_s= irr_surface_Temp(r_i, kp_Ts)
            kp_stellar = Qabs_interp(tstar)*(3/4)*(1/sil_density)*10**4

            kp_Ti =Qabs_interp(Temp)*(3/4)*(1/sil_density)*10**4
            kp_Ti = planckOpacity(Metallicity,T_stellar,Temp,P_gas,Rho_gas)
            psi_i = self_absorbed_fraction(sigma_iter, kp_Ti)
            psi_s = flux_absorbed_fraction(sigma_iter, kp_Ts)
            #print("Printing , psi value, ",psi_i, psi_s)
            Ti_iter =  w*Temp +(1-w)*irr_disk_temp(alpha_guess, r_i,psi_i, psi_s)
            error = np.max(np.abs(Temp - Ti_iter))
      
          
        
        h_cg_iter = pressure_height(Ti_iter,r_i)
        X_cg_iter = chi_cg(sigma_iter,kp_stellar)
        H_cg_iter = X_cg_iter*h_cg_iter
        #impinge_angle_iter =converging(alpha_guess, r_i)    
        #np.append(Ti, Ti_iter)
        Ti[i] = Ti_iter
        h_cg[i] =h_cg_iter
        sigma[i] =sigma_iter
        X_cg[i] =X_cg_iter
        H_cg[i] =H_cg_iter
        impinge_angle[i] =impinge_angle_iter
    return Ti,h_cg,sigma,X_cg,H_cg,impinge_angle

T_irr,h_cg,sigma,X_cg,H_cg,impinge_angle =Q_irr(R, psi_i, psi_s,alpha_guess)

 #alpha*lstar/(4*np.pi*R**2)
def Q_irr2(R,alpha_guess):
    T = (0.5)*(rstar/R)**0.5 *tstar
    for i in range(1,nr-1):
        alpha = converging_alpha(alpha_guess, R[i], psi_i, psi_s)
        T[i] = ((alpha*lstar/(4*np.pi*R[i]**2))/(2*ss))**0.25#alpha*lst00000ar/(4*np.pi*R**2)
    return T

T_irr2 = Q_irr2(R,alpha_guess)

"""
fig, ax = plt.subplots()
ax.semilogx(R/AU, T_irr2)
ax.set_xlabel("Radius")
ax.set_ylabel("Ti")
ax.set_title("Ti vs. Radius")
plt.show()   
"""

fig, ax = plt.subplots()
ax.semilogx(R/AU, T_irr)
ax.set_xlabel("Radius")
ax.set_ylabel("T_irr")
ax.set_title("T_irr vs. Radius")
plt.show()


#####################*Without shadowed region*#################################
T_effective = (T_irr**4  +T_visc**4 +T_cond**4)**0.25
h_cg= pressure_height(T_effective, R)
sigma =surface_density(R)
X_cg =chi_cg(sigma, kp_stellar)
H_cg =scale_height(X_cg, h_cg)
impinge_angle =alpha_angle(H_cg, R)

fig, ax = plt.subplots()
ax.semilogx(R/AU, T_effective)
ax.set_xlabel("Radius")
ax.set_ylabel("T_effective")
ax.set_title("Assuming no shadowed region")
plt.show()

Rflare = np.where (R <H_cg*R_rim/H_rim  )#and H_cg*R_rim/R/H_rim >1)
Rshadow = np.where (R >=H_cg*R_rim/H_rim)

T_irr[Rshadow] = 0
T_effective = (T_irr**4  +T_visc**4 +T_cond**4)**0.25
h_cg = pressure_height(T_effective, R)
X_cg = chi_cg(sigma,kp_stellar)
H_cg = scale_height(X_cg, h_cg)

h_cg[0] = h_rim 
H_cg[0] = H_rim
X_cg[0] = X_rim
T_effective[0] =T_rim

z_elevation =zi = np.linspace(0.0, H_cg, nz)
#z_elevation = np.array(z_elevation)
#z_elevation   = np.array([0.5 * ( zi[i] + zi[i+1])   for i in range(len(zi) - 1)])  # Take the in-"between values       # Take the in-between values    
Z = np.linspace(0.0, 50*AU, nz)



#############################*Stability*cretrion*##############################

def Stability(h_cg,sigma,R,T_effective):
    Q_stability = (np.sqrt(GG*mstar/R**3)*np.sqrt(kb*T_effective/mu/mp))/(np.pi*GG*sigma)
    return Q_stability

Q_stability = Stability(h_cg,sigma,R,T_effective)



fig, ax = plt.subplots()
ax.semilogx(R/AU, T_effective, marker = '*', markevery = (Q_stability <1))
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
ax.semilogx(R/AU, H_cg/R,marker = '*', markevery = (Q_stability <1))
ax.axhline(y=H_rim/R_rim, xmin=0, xmax=100, color = 'b', linestyle = '--')

ax.set_xlabel("Radius")
ax.set_ylabel("Surface_Height")
ax.set_title("Surface_Height vs. Radius")
#ax.set_ylim([0, 0.5])

plt.show()


fig, ax = plt.subplots()
ax.semilogx(R/AU, h_cg/R,marker = '*', markevery = (Q_stability <1))

ax.set_xlabel("Radius")
ax.set_ylabel("Pressure_Height")
ax.set_title("pressure_Height vs. Radius")
ax.set_ylim([0, 0.5])

plt.show()



fig, ax = plt.subplots()
ax.semilogx(R/AU, surface_density(R),marker = '*', markevery = (Q_stability <1))

ax.set_xlabel("Radius")
ax.set_ylabel("Surface Density")
ax.set_title("Surface Density vs. Radius")

plt.show()



fig, ax = plt.subplots()
ax.semilogx(R/AU, X_cg)


ax.set_xlabel("Radius")
ax.set_ylabel("X_cg")
ax.set_title("x_cg vs. Radius")
plt.show()   

 

fig, ax = plt.subplots()
ax.plot(R/AU, H_cg*R_rim/R/H_rim,marker = '*', markevery = (Q_stability <1))
ax.set_xlabel("R_cg")
ax.set_ylabel("Flaring_point")
ax.set_title("Flaring_point vs. Radius")
plt.show()   
 


fig, ax = plt.subplots()
ax.semilogx(R/AU, Q_stability)
ax.set_xlabel("Radius")
ax.set_ylabel("Toomre criterion")
ax.set_title("Ioomre criterion vs. Radius")
plt.show()   


