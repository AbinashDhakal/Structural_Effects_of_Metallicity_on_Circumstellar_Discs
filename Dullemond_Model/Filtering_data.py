import numpy as np
import pandas as pd

def closest(lst, K):
      
     lst = np.asarray(lst)
     idx = (np.abs(lst - K)).argmin()
     return lst[idx]
      
# Driver code

df = pd.read_csv("Planck_data.csv")
#df.to_numpy()[0].T  #create list of first row

#print(df.head(3))
#print(df[df["T_stellar"] >=4000])




closet_temp =(closest(df["T_stellar"], 5001))

#print(df.head(3))

#(df["Rho_gas"] ==(closest(df["Rho_gas"],  3.992E-20)))
Tstar = np.linspace(3000,5000,10)
Tstar = np.array(Tstar)
planck_opacity = []
planck_opacity = np.array(planck_opacity)
for T_i in range (len(Tstar)):
    filtering_data = ( df[
        (df["T_stellar"] == (closest(df["T_stellar"], T_i))) 
        & (df["T_gas"] == (closest(df["T_gas"],600)))
        #& (df["Rho_gas"] == (closest(df["Rho_gas"],  3.7-20)))
        #& (df["P_gas"] == (closest(df["P_gas"],  1)))
       ])
    planck_opacity = filtering_data.to_numpy()[T_i][5]

#planck_opacity = filtering_data.Planck_opacity
"""
for T_i in range(len(Tstar)):
    print(df[
        (df["T_stellar"] == (closest(df["T_stellar"], T_i))) 
        & (df["T_gas"] == (closest(df["T_gas"],700)))
        & (df["Rho_gas"] == (closest(df["Rho_gas"],  5.5900000e-20)))
        
        ])
"""

