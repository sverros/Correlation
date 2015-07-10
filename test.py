#######################################################
# Code for computing the spatial correlation for a ShakeMap,
# adding to a ShakeMap grid, and computing multiple realizations
# VARIABLES:
#     voi - variable of interest, i.e. PGA
#     r - radius of influence
#     intensity_factor- factor for non-native data
#     num_realization- integer for desired number of realizations
#     vs_corr- Vs30 correlated bool, see JB2009
#     input data- grid.xml, uncertainty.xml, and stationlist.xml
#         stored in Inputs directory
#######################################################
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

voi = 'PGA'
r = np.array([15, 25, 35, 45, 55, 65, 75])
intensity_factor = 0.9
num_realizations = 100
vscorr = True

for R in range(0,np.size(r)):
    radius = r[R]

    # Get shakemap for desired variable, PGA, uncertainty grid and stationdata
    shakemap = ShakeGrid('Inputs/grid.xml', variable = '%s' % voi)

    # Uncertainty Data: Units in ln(pctg)
    uncertainty = ShakeGrid('Inputs/uncertainty.xml', variable= 'STD%s' % voi)

    # Station Data: Units in pctg
    stationlist = 'Inputs/stationlist.xml'
    stationdata = readStation(stationlist)

    print 'Calling initialize'
    variables = initialize(shakemap, uncertainty, stationdata, vscorr)

    print variables['K'], ' stations', variables['N']*variables['M'], ' points'

    # Calculate random vector
    rand = np.random.randn(variables['N']*variables['M'])
    # rand = []
    # 
    # with open('rand.txt') as f:
    #     for line in f:
    #         rand.append(float(line))
    # f.close()

    print 'Calling main'
    out = main(variables, radius, voi, rand, vscorr, intensity_factor)

    print 'Computing realizations'
    realizations(num_realizations, radius, variables['N'], variables['M'], out['grid_arr'], 
                 out['mu_arr'], out['sigma_arr'], variables['uncertaintydata'], out['data'])

    plot(out, variables, voi, shakemap, stationdata)


