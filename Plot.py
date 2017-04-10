'''
Quick code to plot histogram of data for use in NUmerical Recipes Report 2
'''

import sys
import numpy as np 
import pylab as pl

#load decay length data from text file
y = np.loadtxt(sys.argv[1])

#histogram of data
pl.hist( y, 100)

# label the axis and add a title to the plot
pl.xlabel('x')
pl.ylabel('y')
pl.savefig('generated_data.png')
pl.show()


