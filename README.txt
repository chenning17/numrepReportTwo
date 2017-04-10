Numerical Recipes Code for Report 2 - Calum Henning
____________________________________________________________________________________________________________________

The python script 'MLL.py' carries out a Maximum Likelihood fit to the supplied data held in 'DecayTimesData.txt'
and calculates the associated errors in the values.

    It can be run using the command:                'python MLL.py DecayTimesData.txt'
(alternatively if the argument text file is changed, the program can fit to other data, e.g. data from a simulation)
    
Plots of NLL versus varying parameter values are also made.

--------------------------------------------------------------------------------------------------------------------

The python script 'MLLbox.py' generates simulated data for the parameters calculated in the previous code and
calculates the associated errors in the values.

    It can be run using the command:                'python MLLbox.py'

Plots of NLL versus varying parameter values are also made.

____________________________________________________________________________________________________________________

A few other scripts and files are present purely for self convenience or are output files from the main scripts:
    Plot.py     <dataFile.txt>              -  generates a histogram plot of the data file 
    
    maxVal.py   <dataFile.txt>              -  returns maximum and minimum values in data file
    
    DecayTimesData.txt                      -  contains the supplied data set of decay times
    
    run.cmd                                 -  contains the commands required to run the main scripts - if using 
                                               windows double click to run them for convenience
    
O   generated_data.txt                      -  contains values generated from running MLLbox.py
    
O   experiment_plot_<0,1 or 2>.png          -  output plot of NLL versus varying parameter. 0,1,2 values represent 
                                               fraction1, tau1,tau2 were parameter that was varied
    
O   sim_plot_vars_<number>_<0,1 or 2>.png   -  output plot of NLL versus varying parameter. 0,1,2 values again 
                                               represent fraction1, tau1,tau2 were parameter that was varied, number 
                                               is the value of this parameter
                                               
O = file is a generated output file

____________________________________________________________________________________________________________________
All code written and used in this project can be found here: https://github.com/chenning17/numrepReportTwo