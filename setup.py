import numpy as np
import math
from openquake.hazardlib.site import Site, SiteCollection
from openquake.hazardlib.geo import Point
import time


def initialize(shakemap, uncertainty, stationdata, vscorr, dm=1, dn=1):
#######################################################################
# Set up grid spacing, site collections, and other data values
#
# INPUTS: 
#     shakemap- shake grid of grid.xml
#     uncertainty- shake grid of uncertainty.xml
#     stationdata- data from stationlist.xml
#     vscorr - boolean, determines if Vs30 are correlation. See
#     JB2009
#     dm, dn- vertical and horozontal spacing, respectively. Default value is 1
# OUT:
#     Returns a dictionary with the following keys:    
#     M,N - number of points vertically and horozontally
#     K- number of stations
#     uncertaintydata- array of uncertainty data at each point
#     site_collection_SM- site collection for ShakeMap data
#     site_collection_station- site collection for station data
#     location_lat/lon_g- lat and lons for grid data, for plotting
#     data- ShakeMap data at each point
#     itensity- factor for non-native data
#     site_collection_both- site collection for both SM and stations
####################################################################### 

    start = time.time()

    attributes = shakemap.getAttributes()
    
    # Determine the size of the grid                                                                
    m = attributes['grid_specification']['nlat']
    n = attributes['grid_specification']['nlon']

    # M,N number of vertical, horizontal points considered
    M = int(math.floor(m/dm))
    N = int(math.floor(n/dn))

    # Allocate matrices for storing data and locations
    DATA = np.empty([M,N])
    location_lat_g = np.empty([M,N])
    location_lon_g = np.empty([M,N])
    site_SM = []
    uncertaintydata = np.empty([M,N])

    # Puts griddata into Site class, then turns the sites into a SiteCollection
    for i in range(0,M):
        for j in range(0,N):

            DATA[i,j] = shakemap.griddata[i*dm,j*dn] # UNITS pctg  

            lat,lon = shakemap.getLatLon(i*dm,j*dn)
            location_lat_g[i,j] = lat
            location_lon_g[i,j] = lon
            uncertaintydata[i,j] = uncertainty.griddata[i*dm,j*dn] # UNITS ln(pctg)            

            site_SM.append(Site(location = Point(location_lon_g[i,j], location_lat_g[i,j]), 
                                vs30 = 760, vs30measured = vscorr, z1pt0 = 100, z2pt5 = 1))

    site_collection_SM = SiteCollection(site_SM)

    # Determine number of stations
    K = np.size(stationdata['lat'])

    # Store lat lon points
    location_lat_s = np.empty([K])
    location_lon_s = np.empty([K])
    intensity = np.empty([K])
    site_station = []

    # Puts stationdata into Site class, then turns the sites into a SiteCollection
    for i in range(0,K):
        location_lat_s[i] = stationdata['lat'][i]
        location_lon_s[i] = stationdata['lon'][i]

        if stationdata['name'][i] == 'DERIVED':
            intensity[i] = 1
        else:
            intensity[i] = 0

        site = Site(location=Point(location_lon_s[i], location_lat_s[i]),
                    vs30 = 760, vs30measured = vscorr, z1pt0 = 100, z2pt5 = 1)
        site_station.append(site)
    
    site_collection_station = SiteCollection(site_station)

    site_both = site_station + site_SM
    site_collection_both = SiteCollection(site_both)
    
    end = time.time() - start
    print '\tInitialization Time:', end
    
    return {'M': M, 'N': N, 'K': K, 'uncertaintydata':uncertaintydata, \
                'site_collection_SM':site_collection_SM, \
                'site_collection_station':site_collection_station, \
                'location_lat_g': location_lat_g, 'location_lon_g':location_lon_g, \
                'data':DATA, 'intensity':intensity, \
                'site_collection_both':site_collection_both}
    
