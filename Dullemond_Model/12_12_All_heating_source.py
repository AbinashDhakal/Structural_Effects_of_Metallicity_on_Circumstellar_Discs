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
Metallicity =0.3
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



def density_distribution(z_elevation, h_cg, sigma):
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


"""
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
kp_stellar =1706.2401533687657
kp_Trim =973.8276480930514


#################################*Modelling Inner hole*########################
#https://arxiv.org/pdf/1905.00909.pdf 
def gas_density(R):
    return 2.5*(R/56/AU)**0.3 *np.exp(-(R/56/AU)**(2-0.3))





def innerHole(Ri):
    Tguess = T_rim
    Tgas = 9200
    data = np.loadtxt('planck_data.csv', delimiter=',', skiprows=1)
    modified_data = np.loadtxt('planck_data.csv', delimiter=',', skiprows=1,  usecols=(0,1, 2,3,4))
    #print(Tgas, Tguess)
    while np.round(Tguess, 2) != np.round(Tgas,2):
        Tguess = Tgas
        s_gas = gas_density(Ri)
        h_gas =pressure_height(Tgas, Ri)
        #Rho_gas = 10**-13
        Rho_gas = density_distribution(0, h_gas, s_gas)
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
        search_vect2 =  np.array([Metallicity, Tgas, Tgas ,P_gas, Rho_gas])
        # normalize data in each column with the column mean
        norm_data2 = modified_data /np.mean(modified_data, axis=0)
        #print("printing normalised vecotr",norm_data)
        # normalize the search vector with the same values
        s_vec_n2 = search_vect2/np.mean(modified_data , axis=0)
        idx2 = np.argmin(np.linalg.norm((s_vec_n2 - norm_data2), axis=1))
        kp_Tgas = data[idx2][5]
        
        kp_Tgas2 =planckOpacity(Metallicity,Tgas,Tgas,P_gas,Rho_gas)
        print(kp_Tgas,kp_TgasTrad, kp_Tgas2)
        print(idx, idx2)
        
        Tgas =(kp_Tgas/kp_TgasTrad)**(-0.25)*(rstar/Ri/2)**(0.5) *tstar
        #Tgas = ((rstar/Ri)**2 *np.exp(-Tauuu )*Rho_gas * kp_T* tstar**4)**0.25
        
    return Tgas


def innerHole_Temp(R_hole):
    for i, ri in enumerate(R_hole):
        Temp_innerhole[i] = innerHole(ri)
    return Temp_innerhole


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



#calculating rim temperature
R_hole = np.linspace(rstar,R_rim,nr)
Rim_temp = innerHole(R_rim)
print(innerHole(R_rim))
Temp_innerhole = np.zeros(nr)

sigma_gas = gas_density(R_hole)



"""



#############################*Defining Radial limits*##########################

R_rim, H_rim, X_rim , h_rim =8078284934830.44, 1134148689469.117, 3.9222797831411085, 289155479000.7479
kp_stellar =1706.2401533687657
kp_Trim =973.8276480930514

rin = R_rim # Sublimation radius  # inner rim of disk. Density is assumed to
rout = 100*AU
R =ri = np.linspace(rin,rout,nr) 
R = np.array(R)


#################################*Inital temperature profile*##################
#Inital approximation
T_inital_approx = (0.5*rstar/R)**0.5 * tstar
Qvisc= np.zeros(nr)
Qcond= np.zeros(nr)
Qirr= np.zeros(nr)
Qheat = np.zeros(nr)

################################*Defining heating source*######################
#Q_visc = Rho* vt*(r*dr_omega)**2
#Q_visc = (sigma/(2*pi*h))*(visc_const)*(np.sqrt(kb*T/mu/mh))*(9*G*Mstar/4/R**3)

#T_inital_approx = (0.5*rstar/R)**0.5 * tstar

#T_visc = (3*GG*mstar*M_acc/(8*np.pi* ss* R**3))**0.25
#(R*smp.diff(smp.sqrt(GG*mstar/R**3),R))**2 =9Gm*/4R^3


#Cant I make Qvisc as grobal variable but rather than calling Qvisc
def Q_visc(R,T,Qvisc):
    yrs = 31557600 #s
    M_acc = 10**-8*ms/yrs
    #print(T_visc)
    for i in range(0,len(T)):    
        sigma = surface_density(R[i])
        cs = np.sqrt(kb*T[i]/mu/mp)
        Qvisc[i] =(sigma/np.sqrt(2*np.pi)*visc_const*cs*(9*GG*mstar/4/R[i]**3))
    #T_visc = (3*GG*mstar*M_acc/(8*np.pi* ss* R**3))**0.25
    return Qvisc

#For some reason it doesnt include 
#Q_cond = kt∇^2T
#kt : thermal conductivity
def Q_cond(R,T,Qcond):
    #Defining a potential form of T
    dr = T[1]-T[0]
    for iteration in range(0, 10): 
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
            Qcond[i] =kt*(dT_dR_2 +dT_dR/R[i])
            #Q_cool = 2 *ss* T
            #Q_condct = kt * ΔT
            #Q_cool = Q_conduct
    return Qcond



X_cg = np.zeros(nr)
impinge_angle = np.zeros(nr)
H_cg = np.zeros(nr)
h_cg = np.zeros(nr)
sigma = np.zeros(nr)
T_irr2 = np.zeros(nr)


Tirr_guess = 800
#convering iteration such that the it iterate until the error is smaller than 0.00001
#THis only check for a radial point
#converging_alpha(T_inital_approx[0], R[0], psi_i, psi_s,Tirr_guess)
def converging_alpha(T_i, r_ii, psi_i, psi_s,Tirr_guess):
    T_i = 1550
    while np.round(T_i, 5) != np.round(Tirr_guess,5):
    #while alpha_guess !=a_i:
        Tirr_guess = T_i
        h_i= pressure_height(Tirr_guess,r_ii)
        s_i = surface_density(r_ii)
        X_i = chi_cg(s_i,kp_stellar)
        H_i = X_i*h_i
        a_i =alpha_angle(H_i, r_ii)
        T_i =irr_disk_temp(a_i, r_ii,psi_i, psi_s)
        #print(T_i,Tirr_guess)
    return a_i

#T_inital_approx = (0.5*rstar/R)**0.5 * tstar
#Q_irr(R,T_inital_approx,Tirr_guess,Qirr)
 #alpha*lstar/(4*np.pi*R**2)
def Q_irr(R,T,Tirr_guess,Qirr):
    for i in range(0,len(T)):
        alpha = converging_alpha(T[i], R[i], psi_i, psi_s,Tirr_guess)
        Qirr[i]= (alpha*lstar/(4*np.pi*R[i]**2))#alpha*lst00000ar/(4*np.pi*R**2)
    #print(Qirr)
    return Qirr



#T_inital_approx = (0.5*rstar/R)**0.5 * tstar

def Q_emit(R, T):
    for iteration in range(0, 1):
        #T_old = T
        Qheat = Q_irr(R,T,Tirr_guess,Qirr) +Q_cond(R,T,Qcond) +Q_visc(R,T,Qvisc)
        T = (Qheat/2/ss)**0.25
    #print(T)
    print(Qheat)
    return T
    
#With all heating source
T_effective = Q_emit(R, T_inital_approx)
h_cg= pressure_height(T_effective, R)
sigma =surface_density(R)
X_cg =chi_cg(sigma, kp_stellar)
H_cg =scale_height(X_cg, h_cg)

Rflare = np.where (R <H_cg*R_rim/H_rim  )#and H_cg*R_rim/R/H_rim >1)
Rshadow = np.where (R >=H_cg*R_rim/H_rim)
shadowed_length = R[len(np.where (R >=H_cg*R_rim/H_rim)[0]) +1] -R[0]

#Initally make shadowed region temperature 0
T_effective[Rshadow] = 0

def Q_secondaryheating(R, T):
    for iteration in range(0, 1):
        T_old = T
        T[0] =T_rim
        Qheat = Q_cond(R,T,Qcond) +Q_visc(R,T,Qvisc)
        T = (Qheat/2/ss)**0.25
    
    return T


T_secondary = Q_secondaryheating(R, T_inital_approx)

#Making flare part 0
T_secondary[Rflare] =0

#T final should include 
#   T_effective (all the heating source: This is used for flaring part)
#T_secondary (It include heating from viscous and conduction)

# T    = T[flare] +     T[shadowed]
T_final = T_effective + T_secondary
#I am just renaming as Teffective
T_effective = T_final
h_cg= pressure_height(T_effective, R)
sigma =surface_density(R)
X_cg =chi_cg(sigma, kp_stellar)
H_cg =scale_height(X_cg, h_cg)









"""
fig, ax = plt.subplots()
ax.semilogx(R/AU, T_irr2)
ax.set_xlabel("Radius")
ax.set_ylabel("Ti")
ax.set_title("Ti vs. Radius")
plt.show()   
"""




#####################*Without shadowed region*#################################
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
shadowed_length = R[len(np.where (R >=H_cg*R_rim/H_rim)[0]) +1] -R[0]
print("Shadowed_length and rim temp ",shadowed_length/AU)



