#import matplotlib.pyplot as plt
import datetime
from neicio.readstation import readStation
from neicio.shake import ShakeGrid
import numpy as np
#from openquake.hazardlib.correlation import JB2009CorrelationModel
#from openquake.hazardlib.correlation import BaseCorrelationModel
#from openquake.hazardlib.geo.geodetic import geodetic_distance
#from openquake.hazardlib.imt import from_string
import time
#from mpl_toolkits.basemap import Basemap
from matplotlib import cm
from neicio.gmt import GMTGrid
#from openquake.hazardlib.imt import PGA as IMT
import sys
sys.path.append('/Users/sverros/Documents/Modules')

from Correlation.setup import initialize
from Correlation.loop import main



# Variable of interest                                                                                                                                                                       
voi = 'PGA'
# Specify the radius of interest
r = 6

# Get shakemap for desired variable, PGA, uncertainty grid and stationdata                                                                                                                   
# Selected Stations: Units in pctg                                                                                                                                                           
shakemap = ShakeGrid('/Users/sverros/Documents/output.all.sta/grid.xml', variable = '%s' % voi)

# Uncertainty Data: Units in ln(pctg)                                                                                                                                                        
uncertainty = ShakeGrid('/Users/sverros/Documents/output.all.sta/uncertainty.xml', variable= 'STD%s' % voi)

# Station Data: Units in pctg                                                                                                                                                                
stationlist = '/Users/sverros/Documents/output.all.sta/stationlist.xml'
stationdata = readStation(stationlist)

# Used for plotting 
topofile = '/Users/sverros/Documents/etopo1_bed_g_f4.grd'
print 'Calling initialize'
variables = initialize(shakemap, uncertainty, stationdata)
print 'Calling main'
rand = np.random.randn(variables['N']*variables['M'])
outputs = main(variables, r, voi, rand)


