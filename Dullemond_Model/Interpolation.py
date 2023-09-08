import numpy as np
from scipy.interpolate import interp2d
import pandas as pd
from scipy import interpolate

"""
x = np.linspace(0., 1., num=5)
y = np.linspace(0., 1., num=5)
z = np.random.randint(1, 5, size=(5, 5))

f = interp2d(x, y, z, kind='linear', fill_value='-1')
print(f(0.111, 0.2533))
"""

df = pd.read_csv("Planck_data.csv")
filtering_data =( df[
    (df["T_stellar"] == 3000) 
    & (df["Metallicity"] == -0.3 )
   ])
Metallicity_data = np.array(list(filtering_data["Metallicity"]))

T_stellar_data = np.array(list(filtering_data["T_stellar"]))

T_data = np.array(list(filtering_data["T_gas"]))
P_data = np.array(list(filtering_data["P_gas"]))
Rho_data = np.array(list(filtering_data["Rho_gas"]))

kp_data = np.array(list(filtering_data["Planck_opacity"]))

x = np.linspace(1,6, 6)
y = np.linspace(1,len(filtering_data),len(filtering_data))
z = np.random.randint(1, 5, size=(len(x), len(y)))
z[:,0] =Metallicity_data
z[:,1] =T_stellar_data
z[:,2] =T_data
z[:,3] =P_data
z[:,4] = Rho_data
z[:5] = kp_data

f = interp2d(x, y, z, kind='linear', fill_value='-1')
print(f(-0.3	, 3000,700, 1.00E-09, 3.99E-20))








#f = interpolate.interp2d(x,y,test,kind='cubic')


f2 = interp2d(T_data, P_data,Rho_data, kp_data, kind='linear', fill_value='-1')
print(f2(700, 1.00E-09, 3.99E-20))

"""
kp_t_array = np.zeros((2, len(kp_data)))
kp_t_array[0] =T_data
kp_t_array[1] =kp_data
"""

