import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator
import pandas as pd
import os
import time

T_Stellar = 9000
Metallicity = -0.3


def kp_Tgas_interpolation(Input_array,Output_array ,request):
    linInter= LinearNDInterpolator(Input_array, Output_array)
    output = linInter(request)
    iteration = 1
    while output != output:
        # Vector we try to match with data
        search_vect = request
        # normalize data in each column with the column mean
        norm_data = Input_array/np.mean(Input_array, axis=0)
        # normalize the search vector with the same values
        s_vec_n = search_vect/np.mean(Input_array, axis=0)
        # Define closest fit as the smallest norm of the difference vector
        idx = np.argmin(np.linalg.norm((s_vec_n - norm_data), axis=1))
        output = Output_array[idx]
    """
    iteration = 1
while output != output:
    request_again =  np.round(request*np.random.uniform(0.9, 1.1),1)
    print( "kp_TgasTrad_interpolation",request_again)
    request = request_again
    output = Interp(request_again)
    
    #print(output)
    iteration += 1
    if iteration ==10:
        #request =np.array([np.log10(5732)    ,-0.28567,-11.8827])
        #output = linInter(request)
        return 1"""
        
        
    return output


"""
def kp_Tgas_interpolation(Input_array,Output_array ,request):
    linInter= LinearNDInterpolator(Input_array, Output_array)
    try:
        output = linInter(request)
    except:
        request = request_again
        output = linInter(request_again)
    return output


"""



def kp_Tgas_interpolation(Input_array,Output_array ,request):
    linInter= LinearNDInterpolator(Input_array, Output_array)
    
    output = linInter(request)
        #This calculate linear interpolation of a function
#Its not necessary, it will find a value so if it doesnt
#it return nan as a result  
    while output != output:
        # Vector we try to match with data
        search_vect = request
        # normalize data in each column with the column mean
        norm_data = Input_array/np.mean(Input_array, axis=0)
        # normalize the search vector with the same values
        s_vec_n = search_vect/np.mean(Input_array, axis=0)
        # Define closest fit as the smallest norm of the difference vector
        idx = np.argmin(np.linalg.norm((s_vec_n - norm_data), axis=1))
        output = Output_array[idx]
        #print(idx)
        #print(output)
    return output





md = pd.read_csv("Planck_gas.csv")
md= md.loc[md['Metallicity']==Metallicity]
md['Log_T_gas'] = np.log10(np.array(md['T_gas'])) 

md['Log_Rho_gas'] = np.log10(np.array(md['Rho_gas'])) 
md['Log_P_gas'] = np.log10(np.array(md['P_gas'])) 
Input_array_gas = np.array(md[['Log_T_gas','Log_P_gas', 'Log_Rho_gas']])
Output_array_gas = np.array(md["kp_stellar"])
request1= np.round(np.array([2.8451,	-8.85699,	-19.54]),)
request2 = np.array([  6.        ,   8.71432976,  -5.42365865])
request3 =np.array([np.log10(5732)	,-0.28567,-11.8827])
request4 = np.array([np.log10(1e+06),	8.28556,	-5.85387])

print(kp_Tgas_interpolation(Input_array_gas,Output_array_gas ,request1))
print(kp_Tgas_interpolation(Input_array_gas,Output_array_gas ,request2))
print(kp_Tgas_interpolation(Input_array_gas,Output_array_gas ,request3))
print(kp_Tgas_interpolation(Input_array_gas,Output_array_gas ,request4))
