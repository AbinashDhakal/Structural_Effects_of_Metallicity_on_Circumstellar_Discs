import pandas as pd
import numpy as np


import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfinv
import pandas as pd
from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator


df = pd.read_csv('Nearly_finarl_data.csv')

R_0 = np.array(list(df["R0"]))
R_1 = np.array(list(df["R1"]))
R_2 = np.array(list(df["R2"]))

SH_0 = np.array(list(df["SH0"]))
SH_1 = np.array(list(df["SH1"]))
SH_2 = np.array(list(df["SH2"]))

T0 = np.array(list(df["T0"]))
T1 = np.array(list(df["T1"]))
T2 = np.array(list(df["T2"]))

AU = 1.495978707e13  # Distance between Earth and Sun [cm] defined 

R =ri = np.linspace(0.4*AU,10*AU,500) 

fig, ax = plt.subplots()
ax.semilogx(R_0/AU, T0)
ax.semilogx(R_1/AU, T1)
ax.semilogx(R_2/AU, T2)

ax.set_xlabel("R (AU)")
ax.set_ylabel("T (K)")
plt.show()

fig, ax = plt.subplots()
ax.plot(R_0/AU, T0,label="-0.3")
ax.plot(R_1/AU, T1,label="0")
ax.plot(R_2/AU, T2,label="0.3")

plt.legend(['[Me/H] =-0.3', '[Me/H] =0', '[Me/H] =0.3'], loc='upper right')

#ax.set_xlim([0, 2])
ax.set_xlim([0, 1.5])

ax.set_xlabel("R(AU)")
ax.set_ylabel("T (K)")
plt.show()

fig, ax = plt.subplots()
ax.plot(R_0/AU, SH_0/R_0)
ax.plot(R_1/AU, SH_1/R_1)
ax.plot(R_2/AU, SH_2/R_2)
ax.set_xlim([0, 1.5])

ax.set_xlabel("Radius")
ax.set_ylabel("T_effective")
ax.set_title("Assuming no shadowed region")
plt.show()



fig, ax = plt.subplots()
ax.plot(R_0/AU, SH_0/1e13)
ax.plot(R_1/AU, SH_1/1e13)
ax.plot(R_2/AU, SH_2/1e13)
ax.set_xlim([0, 1.5])
#ax.set_ylim([0*1e13, 1.5*1e13])
ax.set_ylim([0, 0.2])

ax.set_xlabel("Radius")
ax.set_ylabel("T_effective")
ax.set_title("Assuming no shadowed region")
plt.show()