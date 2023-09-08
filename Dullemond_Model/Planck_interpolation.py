import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator
import pandas as pd
df = pd.read_csv("Planck_data.csv")
import os
 
print("File location using os.getcwd():", os.getcwd())
#Planck_gas
# put the available Metallicity, T_Stellar, Tgas,P_gas,Rho_gas data 
#                                                       as a numpy array
points = np.loadtxt('Planck_data.csv', delimiter=',', skiprows=1, 
                           usecols=(0,1, 2,6))


# and corresponding data values in a separate array:
kp_values = np.loadtxt('Planck_data.csv', delimiter=',', skiprows=1, 
                           usecols=(5))

# points to interpolate 
request1= ([-0.3,	3000,	7.00e+02,	-19.3990271])

request2 = np.array([-3.00e-01,  3.00e+03,  7.00e+02,   -19.3990271])
# First, define an interpolator function
linInter= LinearNDInterpolator(points, kp_values, rescale=True)
# Then, apply the function to one or more points
print( linInter(request1))

#  Defining each points
Metallicity = np.loadtxt('planck_data.csv', delimiter=',', skiprows=1, 
                           usecols=(0))
T_stellar = np.loadtxt('planck_data.csv', delimiter=',', skiprows=1, 
                           usecols=(1))

T_gas = np.loadtxt('planck_data.csv', delimiter=',', skiprows=1, 
                           usecols=(2))

P_gas = np.loadtxt('planck_data.csv', delimiter=',', skiprows=1, 
                           usecols=(3))

Rho_gas = np.loadtxt('planck_data.csv', delimiter=',', skiprows=1, 
                           usecols=(6))
# Create coordinate pairs
cartcoord = list(zip(Metallicity,T_stellar,T_gas,P_gas,Rho_gas))

interp= LinearNDInterpolator(cartcoord, kp_values, fill_value=0)
abc = interp(request1)
print( interp(request1))






