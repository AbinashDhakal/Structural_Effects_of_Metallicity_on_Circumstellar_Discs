import numpy, scipy, matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution
import warnings

#x and y are from experiment
#x=[0,1.778,2.921,3.302,6.317,9.524,10.54]
#y=[1,0.831763771,0.598411595,0.656145266,0.207014135,0.016218101,0.004102041]


import pandas as pd
sample_data = pd.read_csv('Ti_dataset.csv')
sample_data2 = pd.read_csv('hR_datasets.csv')
sample_data3= pd.read_csv('HHHR_datasets.csv')

x = numpy.array(sample_data.R)
y = numpy.array(sample_data.Ti)

# alias data to match previous example code
xData = numpy.array(x, dtype=float)
yData = numpy.array(y, dtype=float)


def func(x, a, b, c): # Sigmoidal Gompertz C from zunzun.com
    return (a*(x-b)**(-0.5) +c)
#a * numpy.exp(-b * numpy.exp(-c*x))


# function for genetic algorithm to minimize (sum of squared error)
def sumOfSquaredError(parameterTuple):
    warnings.filterwarnings("ignore") # do not print warnings by genetic algorithm
    val = func(xData, *parameterTuple)
    return numpy.sum((yData - val) ** 2.0)


def generate_Initial_Parameters():
    parameterBounds = []
    parameterBounds.append([2e5,2e7]) # search bounds for a
    parameterBounds.append([-3,3]) # search bounds for b
    parameterBounds.append([500, 1500]) # search bounds for c

    # "seed" the numpy random number generator for repeatable results
    result = differential_evolution(sumOfSquaredError, parameterBounds, seed=3)
    return result.x

# by default, differential_evolution completes by calling curve_fit() using parameter bounds
geneticParameters = generate_Initial_Parameters()

# now call curve_fit without passing bounds from the genetic algorithm,
# just in case the best fit parameters are aoutside those bounds
fittedParameters, pcov = curve_fit(func, xData, yData, geneticParameters)
print('Fitted parameters:', fittedParameters)
print()

modelPredictions = func(xData, *fittedParameters) 

absError = modelPredictions - yData

SE = numpy.square(absError) # squared errors
MSE = numpy.mean(SE) # mean squared errors
RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
Rsquared = 1.0 - (numpy.var(absError) / numpy.var(yData))

print()
print('RMSE:', RMSE)
print('R-squared:', Rsquared)

print()


##########################################################
# graphics output section
def ModelAndScatterPlot(graphWidth, graphHeight):
    f = plt.figure(figsize=(graphWidth/100.0, graphHeight/100.0), dpi=100)
    axes = f.add_subplot(111)

    # plot wuth log Y axis scaling
    plt.xscale('log')

    # first the raw data as a scatter plot
    axes.plot(xData, yData,  'D')

    # create data for the fitted equation plot
    xModel = numpy.linspace(min(xData), max(xData))
    yModel = func(xModel, *fittedParameters)

    # now the model as a line plot
    axes.plot(xModel, yModel)

    axes.set_xlabel('X Data') # X axis data label
    axes.set_ylabel('Y Data') # Y axis data label

    plt.show()
    plt.close('all') # clean up after using pyplot

graphWidth = 800
graphHeight = 600
ModelAndScatterPlot(graphWidth, graphHeight)