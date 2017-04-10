'''
Numerical Recipes Report 2 Code

This code is just a quick way to find out the maximum and minimum values in the 'DecayTimeData.txt' 
file so that reasonable values for 'a' and 'b' can be chosen to be used in the box method of simulation.

Author: Calum Henning
April 2017
'''
import sys
import numpy as np 

#load decay length data from text file
y = np.loadtxt(sys.argv[1])

#print amx and min values to terminal
print(np.amax(y))
print(np.amin(y))



