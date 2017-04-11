'''
Numerical Recipes Report 2 Code

The code reads in values from an input file, and performs a maximum likelihood fit to determine the decay times.

Author: Calum Henning
April 2017
'''

# Imports
import pylab
import numpy as np
import sys
import scipy.optimize as optimize
import matplotlib.pyplot as plt

# Define exponential function of form (1/tau)e^-(t/tau)
def f_exponential(t, tau):
    f = (1/tau)*np.exp(-t/tau)
    return f;

# Define Composite PDF for 2 different decay times
def composite_pdf(t, fraction1, tau1, tau2):
    f = fraction1*f_exponential(t, tau1) + (1-fraction1)*f_exponential(t, tau2)
    return f;
    
# Define negative log likelihood function
def Nll(c, decayTime):
    NLL = 0
    for i in range(len(decayTime)):
        if(composite_pdf(decayTime[i], c[0], c[1], c[2])>0):
            NLL = NLL + np.log(composite_pdf(decayTime[i], c[0], c[1], c[2]))
        else:
            NLL = NLL + np.log(0.0000001)
    NLL = NLL*-1.0
    return NLL;
    
#make plots of the NLL, when one variable is varied by +/- 3 sigma and the others are kept at the best fit values
def plotNLL(variables, option, sigma, data):
     
    #set lower and upper values for the variable being altered
    minX = np.clip(variables[option]-(3*sigma), 0, 1000)
    maxX = np.clip(variables[option]+(3*sigma), 0, 1000)
    #create an array of equally spaced values in this range to be used in the NLL function
    xVals = np.linspace(minX, maxX, num=50, endpoint=True)
    yVals = np.zeros(len(xVals))
        
    #For each new param value (xVals)
    for i in range(len(xVals)):
        temp = list(variables)
        temp[option] = xVals[i]
        yVals[i] = Nll(temp, data)
    
           
    if(option == 0):
        #create plot and set title
        plt.figure(0)
        plt.title('NLL Variation for Decay Fraction')
        plt.xlabel('decay fraction of peak 1')
    elif(option == 1):
        #create plot and set title
        plt.figure(1)
        plt.title('NLL Variation for tau1')
        plt.xlabel('tau1 [s]')
    elif(option == 2):
        #create plot and set title
        plt.figure(2)
        plt.title('NLL Variation for tau2')
        plt.xlabel('tau2 [s]')
    #set up the rest of the plot details
    plt.ylabel('NLL')
    plt.plot(xVals, yVals)
    pylab.savefig("experiment_plot_" + str(option) + ".png")
    #plt.show()
    
#make plots of the NLL, when one variable is varied by +/- 3 sigma and the others are kept at the best fit values
def errorsNLL(variables, option, sigma, data):
     
    #set lower and upper values for the variable being altered
    minX = np.clip(variables[option]-(sigma), 0, 1000)
    maxX = np.clip(variables[option]+(sigma), 0, 1000)
    #create an array of equally spaced values in this range to be used in the NLL function
    xVals = np.linspace(minX, maxX, num=100, endpoint=True)
    yVals = np.zeros(len(xVals))
    
    #bestNLL = Nll(variables, data)
    
    #For each new param value (xVals)
    for i in range(len(xVals)):
        temp = list(variables)
        temp[option] = xVals[i]
        yVals[i] = Nll(temp, data)

    plus_error = 0
    minus_error = 0
    yMin = np.amin(yVals)
    xMin = xVals[np.argmin(yVals)]
    
    #yCopy = [];
    xCopy = [];
    
    #sift out the values that are in the range between the minimum value of the NLL and NLL+0.5, around the minimum value
    for j in range(len(yVals)):
        if(yVals[j] <= yMin+0.5):
            #yCopy.append(yVals[j])
            xCopy.append(xVals[j])
    #locate minimum and maximum values of x in this range, leading to the values corresponding to y +0.5 at either side of the minimum        
    minus_error = xCopy[0]
    minus_error = abs(variables[option] - minus_error)
    plus_error = xCopy[(len(xCopy)-1)]
    plus_error = abs(variables[option] - plus_error)
    #set overall +/- error to be average of these two values
    calc_val = (minus_error + plus_error)/2
    
    #calc_vals = []
    
    #calc_vals.append(minus_error)
    #calc_vals.append(plus_error)
    #print("-ve")
    #print(abs(variables[option] - minus_error))
    #print("+ve")
    #print(abs(variables[option] - plus_error))
    
    return calc_val;
    
def main():
    # Open the file and read in the data
    #decayTime = np.loadtxt("DecayTimesData.txt")
    #decayTime = np.loadtxt("simulated.dat")
    decayTime  = np.loadtxt(sys.argv[1])

    #initial guesses:  (based on previous runs, but could be set to different values and would still function correctly)
    fraction1 = 0.2
    tau1 = 1.3
    tau2 = 0.2
    #calculate maximum likelihood values for the values by minimising the negative log likelihood
    result = optimize.minimize(Nll, [fraction1, tau1, tau2], decayTime, method='Nelder-Mead')

    #print result
    fraction1 = result.x[0]
    tau1 = result.x[1]
    tau2 = result.x[2]

    #Calculate errors in values
    #guess an initial value for sigma (that is at least larger than the largest error)
    sigma = 0.035
    variables = [fraction1, tau1, tau2]
    frac1_errors = errorsNLL(variables, 0, sigma, decayTime);
    tau1_errors = errorsNLL(variables, 1, sigma, decayTime);
    tau2_errors = errorsNLL(variables, 2, sigma, decayTime);

    #print out results and associated errors to the terminal
    print ("The fraction of peak one is: " + str("%.3f" % (fraction1)) + " +/- " + str("%.3f" % (frac1_errors)))
    print (" ")
    print ("The decay time of peak one is: " + str("%.2f" % (tau1))+"s" + " +/- " + str("%.2f" % (tau1_errors)))
    print (" ")
    print ("The decay time of peak two is: " + str("%.3f" % (tau2))+"s" + " +/- " + str("%.3f" % (tau2_errors)))
    print (" ")

    #Plot the variation of the NLL when each variable is varied by 3 sigma about the minimum
    variables = [fraction1, tau1, tau2]
    plotNLL(variables, 0, frac1_errors, decayTime);
    plotNLL(variables, 1, tau1_errors, decayTime);
    plotNLL(variables, 2, tau2_errors, decayTime);
    
main();