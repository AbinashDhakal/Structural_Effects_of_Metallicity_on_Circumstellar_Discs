import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator
import pandas as pd
import os
import time

T_Stellar = 9000
Metallicity = -0.3

df = pd.read_csv("Planck_data.csv")
df = df.loc[df['T_stellar']==T_Stellar]
df = df.loc[df['Metallicity']==Metallicity]
df['Log_Rho_gas'] = np.log10(np.array(df['Rho_gas'])) 
df['Log_P_gas'] = np.log10(np.array(df['P_gas']))    
Input_array = np.array(df[['T_gas','Log_P_gas', 'Log_Rho_gas']])
Output_array = np.array(df["kp_stellar"])

def kp_TgasTrad_interpolation(Input_array,Output_array ,request):
    Interp= LinearNDInterpolator(Input_array, Output_array)
    return Interp(request)


request = np.array([700,	-8.14267,-18.5421])
print(kp_TgasTrad_interpolation(Input_array,Output_array ,request))

def kp_Tgas_interpolation(Input_array,Output_array ,request):
    linInter= LinearNDInterpolator(Input_array, Output_array)
    
    output = linInter(request)
        #This calculate linear interpolation of a function
#Its not necessary, it will find a value so if it doesnt
#it return nan as a result  
    return output
    while np.isnan(output) != np.nan:
        #if interpolation return nothing, I want to re-evaluate
        #it choose different point (something close to re-evaluate)
        request = request*np.random.uniform(0.7, 1)
        print(output)
        print(request)
        output = linInter(request)
    return output
"""
    if np.isnan(output):
        return 1
   output = linInter(request)
    while output.isnan():
        #if interpolation return nothing, I want to re-evaluate
        #it choose different point (something close to re-evaluate)
        request = request*np.random.uniform(0.7, 1)
        output = linInter(request)
        
        #I want to loop until I get some number rather than nan
        print(request)
        iteration += 1
        print(iteration)"""

md = pd.read_csv("Planck_gas.csv")
md= md.loc[md['Metallicity']==Metallicity]
md['Log_T_gas'] = np.log10(np.array(md['T_gas'])) 

md['Log_Rho_gas'] = np.log10(np.array(md['Rho_gas'])) 
md['Log_P_gas'] = np.log10(np.array(md['P_gas'])) 
Input_array_gas = np.array(md[['Log_T_gas','Log_P_gas', 'Log_Rho_gas']])
Output_array_gas = np.array(md["kp_stellar"])
request1= np.round(np.array([2.8451,	-8.85699,	-19.2557]),3)
request2 = np.array([  6.        ,   8.71432976,  -5.42365865])
request3 =np.array([np.log10(5732)	,-0.28567,-11.8827])
request4 = np.array([np.log10(1e+06),	8.28556,	-5.85387])

print(kp_Tgas_interpolation(Input_array_gas,Output_array_gas ,request1))
print(kp_Tgas_interpolation(Input_array_gas,Output_array_gas ,request2))
print(kp_Tgas_interpolation(Input_array_gas,Output_array_gas ,request3))
print(kp_Tgas_interpolation(Input_array_gas,Output_array_gas ,request4))