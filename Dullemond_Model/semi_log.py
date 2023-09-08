import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import curve_fit

x=np.linspace(0, 110,1000)
y=np.linspace(0,1500,1000)
a,b,c = 1,5,0.4
#plt.semilogx(5*(0.95*1/x+10 )**2 +0.5+ x**5,            y)
plt.semilogx((1+np.exp(5*(-50*x**2)/1000))/14, y)


plt.xlabel("Omega")
plt.ylabel("phase")
plt.show()

"""
from sympy import Eq, solve
from sympy.abc import x ,a,b,c,d
sol = solve([ Eq((a**2 + b*0.52 +c)**(-1) +d, 1500),
              Eq((a*10**2 + b*10+c)**(-1) +d, 200),
              Eq((a*100**2+ b*100+c)**(-1) +d, 190),
              Eq((a*0.6**2+b*0.6+c)**(-1) +d, 950)])

print({ s:sol[s].evalf() for s in sol })
"""

import pandas as pd
sample_data = pd.read_csv('Ti_dataset.csv')
sample_data2 = pd.read_csv('hR_datasets.csv')
sample_data3= pd.read_csv('HHHR_datasets.csv')
fig, ax = plt.subplots(figsize = (9, 6))

"""
ax.scatter(sample_data.R, sample_data.Ti, s =20, alpha =0.1, edgecolors = "k")


fig, ax = plt.subplots(figsize = (9, 6))
ax.scatter(sample_data.R, sample_data.Ti, s =20, alpha =0.1, edgecolors = "k")
ax.set_xscale("log")
plt.show()
"""
"""
def model_f(x,a,b,c):
    #return (1 + a*np.exp(b*x**2 +c))
    return a*b**x +c   #Out[25]: array([-6.86661811e-10,  3.00000000e+00, -2.83849861e+31])
    #return a*b**(x+c)

popt, pcov = curve_fit(model_f,sample_data.R,sample_data.Ti , p0 = [-6.86661811e-10,3,-2.83849861e+31], sigma =  None)

a_opt, b_opt, c_opt = popt
x_model = np.linspace(0,110,100000000)
y_model = model_f(x_model, a_opt,  b_opt, c_opt)

fig, ax = plt.subplots(figsize = (9, 6))
ax.scatter(sample_data.R, sample_data.Ti, s =20, alpha =0.1, edgecolors = "k")
ax.plot(x_model, y_model)
ax.set_xscale("log")
plt.show()
xData = np.array(sample_data.R, dtype=float)

import sys
import os
import numpy
import matplotlib.pyplot as plt
from pylab import *
from scipy.optimize import curve_fit
import scipy.optimize as optimization
from scipy.interpolate import interp1d
from scipy import interpolate


#plt.show()
#plt.semilogx(sample_data.R,sample_data.Ti)





"""
"""

fig, ax = plt.subplots(figsize = (9, 6))
ax.scatter(sample_data2.R, sample_data2.h, s =60, alpha =0.7, edgecolors = "k")
#ax.set_xscale("log")
plt.show()


fig, ax = plt.subplots(figsize = (9, 6))
ax.scatter(sample_data3.R, sample_data3.H, s =60, alpha =0.7, edgecolors = "k")
#ax.set_xscale("log")
plt.show()

"""



fig, ax = plt.subplots(figsize = (9, 6))
ax.scatter(sample_data2.R, sample_data2.h, s =60, alpha =0.7, edgecolors = "k")
ax.set_xscale("log")
plt.show()


#plt.semilogx(sample_data.R, sample_data.Ti)
#plt.scatter(sample_data2.R, sample_data2.h)
#plt.scatter(sample_data3.R, sample_data3.H)
from scipy.interpolate import interp1d
y_f =interp1d(sample_data2.R, sample_data2.h, 'quadratic')
x = np.linspace(min(sample_data2.R),max(sample_data2.R),10000)
y = y_f(x)
plt.semilogx(x, y)
plt.xlabel("R")
plt.ylabel("h")

#plt.scatter(x,y)
print("Hello")
plt.show()




fig, ax = plt.subplots(figsize = (9, 6))
ax.scatter(sample_data.R, sample_data.Ti, s =60, alpha =0.7, edgecolors = "k")
ax.set_xscale("log")
plt.show()


#plt.semilogx(sample_data.R, sample_data.Ti)
#plt.scatter(sample_data2.R, sample_data2.h)
#plt.scatter(sample_data3.R, sample_data3.H)
from scipy.interpolate import interp1d
y_f =interp1d(sample_data.R, sample_data.Ti, 'quadratic')
x = np.linspace(min(sample_data.R),max(sample_data.R),1000)
y = y_f(x)
plt.semilogx(x, y)
plt.xlabel("R")
plt.ylabel("h")

#plt.scatter(x,y)
print("Hello")
plt.show()
