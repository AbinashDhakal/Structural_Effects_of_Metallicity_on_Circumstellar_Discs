import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
from scipy.interpolate import interp1d


sample_data = pd.read_csv('Ti_dataset.csv')
sample_data2 = pd.read_csv('hR_datasets.csv')
sample_data3= pd.read_csv('HHHR_datasets.csv')
fig, ax = plt.subplots(figsize = (9, 6))


R_iterpolate_range = np.linspace(0.5,100,256)

Temp_interpolate_f =interp1d(sample_data.R, sample_data.Ti, 'quadratic')
hR_interpolate_f =interp1d(sample_data2.R, sample_data2.h, 'quadratic')
HR_interpolate_f =interp1d(sample_data3.R, sample_data3.H, 'quadratic')

Temp_interpolate = Temp_interpolate_f(R_iterpolate_range)
hR_interpolate = hR_interpolate_f(R_iterpolate_range)
HR_interpolate = HR_interpolate_f(R_iterpolate_range)

plt.semilogx(R_iterpolate_range, Temp_interpolate)
plt.show

