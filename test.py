#######################################################
# Code for computing the spatial correlation for a ShakeMap,
# adding to a ShakeMap grid, and computing multiple realizations
# VARIABLES:
#     voi - variable of interest, i.e. PGA
#     r - radius of influence
#     num_realization- integer for desired number of realizations
#     corr_model- JB2009 or GA2010
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
r = [45]
num_realizations = 100
corr_model = 'JB2009'
vscorr = True
plot_on = False

for R in range(0,np.size(r)):
    radius = r[R]

    # Get shakemap for desired variable, PGA, uncertainty grid and stationdata
    shakemap = ShakeGrid('Inputs/grid.xml', variable = '%s' % voi)

    # Uncertainty Data: Units in ln(pctg)
    unc_INTRA = ShakeGrid('Inputs/uncertainty.xml', variable= 'GMPE_INTRA_STD%s' % voi)
    unc_INTER = ShakeGrid('Inputs/uncertainty.xml', variable= 'GMPE_INTER_STD%s' % voi)

    # Station Data: Units in pctg
    stationlist = 'Inputs/stationlist.xml'
    stationdata = readStation(stationlist)

    print 'Calling initialize'
    variables = initialize(shakemap, unc_INTRA, unc_INTRA, stationdata)

    print 'Radius: ', radius
    print variables['K'], ' stations', variables['N']*variables['M'], ' points'

    # Calculate random vector
    rand = np.random.randn(variables['N']*variables['M'])
    # Truncate to +/- 3sigma
#    for i in range(0,variables['N']*variables['M']):
#        if abs(rand[i]) > 3:
#            rand[i] = np.sign(rand[i])*(6-rand[i])

    print 'Calling main'
    out = main(variables, radius, voi, rand, corr_model, vscorr)

    if num_realizations > 0:
        print 'Computing realizations'
        realizations(num_realizations, radius, variables['N'], variables['M'], out['grid_arr'], 
                 out['mu_arr'], out['sigma_arr'], variables['uncertaintydata'], out['data'])

    if plot_on == True:
        print 'Plotting results'
        plot(out, variables, voi, shakemap, stationdata)


