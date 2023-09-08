import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator
import pandas as pd
import os
import time

T_Stellar = 9000
Metallicity = -0.3
"""
df = pd.read_csv("Planck_data.csv")
df = df.loc[df['T_stellar']==T_Stellar]
df = df.loc[df['Metallicity']==Metallicity]

df['Log_Rho_gas'] = np.log10(np.array(df['Rho_gas'])) 
Input_array = np.array(df[['Metallicity','T_gas', 'Log_Rho_gas']])

Output_array = np.array(df["kp_stellar"])



interp = LinearNDInterpolator(Input_array, Output_array, rescale=True)

interp(-3.00000000e-01,  7.00000000e+02, -1.93990271e+01)


"""

print("File location using os.getcwd():", os.getcwd())
#Planck_gas
#os.chdir('\Users\abudh\OneDrive\Desktop\Python\Dullemond_Model')
#print("Current working directory: {0}".format(os.getcwd()))

# put the available Metallicity, T_Stellar, Tgas,P_gas,Rho_gas data 
#                                                       as a numpy array
points = np.loadtxt('Planck_0.csv', delimiter=',', skiprows=1, 
                           usecols=(1,3,4))


df = pd.read_csv("Planck_data.csv")
df = df.loc[df['T_stellar']==T_Stellar]
df = df.loc[df['Metallicity']==Metallicity]

df['Log_Rho_gas'] = np.log10(np.array(df['Rho_gas'])) 
Input_array = np.array(df[['Metallicity','T_gas', 'Log_Rho_gas']])
Output_array = np.array(df["kp_stellar"])


#max_rows= 500
#log_points = points.copy()
#log_points[:,(1,2,3)] = np.log(points[:,(1,2,3)])

log_points = np.log(points)


#print(log_points)
#for i in range(0, 3):
#    print(np.max(log_points[:,i]),np.min(log_points[:,i]) )

# and corresponding data values in a separate array:
kp_values = np.loadtxt('Planck_0.csv', delimiter=',', skiprows=1, 
                           usecols=(5))
lop_kp_values = np.log(kp_values)
#print(lop_kp_values)
print(np.max(lop_kp_values),np.min(lop_kp_values) )
#np.max(lop_kp_values)
#Out[3]: 15.432916640047551
#np.min(lop_kp_values)
#Out[4]: -1.987774353154012
# points to interpolate 


#request1= ([-0.3, 8.00637,6.55108, -20.7233,	-44.6679])#7.77E+00
request1= ([8.00637, -20.7233,	-44.6679])#7.77E+00

print(request1)


request2 = ([	30000,-5.101823517,	7.857332496]) #7.857332496


request3 = ([6000,398107,-7.596879479,6.1430148])
#request2 = np.array([-3.00e-01,  3.00e+03,  7.00e+02,   -19.3990271])
# First, define an interpolator function
start_time = time.time()

print("Hello")
linInter= LinearNDInterpolator(log_points, lop_kp_values)
#Then, apply the function to one or more points
print("Hi")
print( linInter(request1))
print("--- %s seconds ---" % (time.time() - start_time))
#[2.05027016] --- 286.8468542098999 seconds --- Take 5 min

"""
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




"""

