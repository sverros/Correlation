import datetime
from neicio.readstation import readStation
from neicio.shake import ShakeGrid
import numpy as np
import time
from matplotlib import cm
import matplotlib.pyplot as plt
from neicio.gmt import GMTGrid
import sys
sys.path.append('/Users/sverros/Documents/Modules')
from Correlation.setup import initialize
from Correlation.loop import main
from Correlation.realizations import realizations
from Correlation.plotting import plot


# Variable of interest                                                                                                                                                                       
voi = 'PGA'
# Specify the radius of interest
r = 15

intensity_factor = 0.9
num_realizations = 1

vscorr = True
# Get shakemap for desired variable, PGA, uncertainty grid and stationdata                                                                                                                   
# Selected Stations: Units in pctg                                                                                                                                                           
shakemap = ShakeGrid('Inputs/grid.xml', variable = '%s' % voi)

# Uncertainty Data: Units in ln(pctg)                                                                                                                                                        
uncertainty = ShakeGrid('Inputs/uncertainty.xml', variable= 'STD%s' % voi)

# Station Data: Units in pctg                                                                                                                                                                
stationlist = 'Inputs/stationlist.xml'
stationdata = readStation(stationlist)

print 'Calling initialize'
variables = initialize(shakemap, uncertainty, stationdata, vscorr)

print variables['K'], ' stations', variables['N']*variables['M'], ' points'

print 'Calling main'
#rand = np.random.randn(variables['N']*variables['M'])
rand = []

with open('rand.txt') as f:
    for line in f:
        rand.append(float(line))
f.close()

out = main(variables, r, voi, rand, vscorr, intensity_factor)

realizations(num_realizations, variables['N'], variables['M'], out['grid_arr'], 
             out['mu_arr'], out['sigma_arr'], variables['uncertaintydata'], out['data'])

plot(out, variables, voi, shakemap, stationdata)


