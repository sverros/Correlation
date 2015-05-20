import datetime
from neicio.readstation import readStation
from neicio.shake import ShakeGrid
import numpy as np
import time
from matplotlib import cm
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
r = 6

intensity_factor = 0.9
num_realizations = 1
# Get shakemap for desired variable, PGA, uncertainty grid and stationdata                                                                                                                   
# Selected Stations: Units in pctg                                                                                                                                                           
shakemap = ShakeGrid('/Users/sverros/Documents/Northridge_Outputs/output.all.sta/grid.xml', variable = '%s' % voi)

# Uncertainty Data: Units in ln(pctg)                                                                                                                                                        
uncertainty = ShakeGrid('/Users/sverros/Documents/Northridge_Outputs/output.all.sta/uncertainty.xml', variable= 'STD%s' % voi)

# Station Data: Units in pctg                                                                                                                                                                
stationlist = '/Users/sverros/Documents/Northridge_Outputs/output.all.sta/stationlist.xml'
stationdata = readStation(stationlist)

print 'Calling initialize'
variables = initialize(shakemap, uncertainty, stationdata, True)
print 'Calling main'
rand = np.random.randn(variables['N']*variables['M'])
outputs = main(variables, r, voi, rand, intensity_factor)

ACCUM_ARRAY = realizations(num_realizations, variables['N'], variables['M'], outputs['grid_arr'], 
                           outputs['mu_arr'], outputs['sigma_arr'], variables['uncertaintydata'], outputs['data'])

plot(outputs, variables, voi, shakemap, stationdata, ACCUM_ARRAY)
