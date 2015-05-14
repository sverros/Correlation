import matplotlib.pyplot as plt
import numpy as np
import scipy
from neicio.gmt import GMTGrid
from mpl_toolkits.basemap import Basemap
from matplotlib import cm


def plot(out, variables, voi, shakemap, stationdata, topofile, ACCUM_ARRAY):
    maxdata = np.amax(out['data_new'])
    attributes = shakemap.getAttributes()
    #station_lons = stationdata['lon']
    #station_lats = stationdata['lat']
    intensity = stationdata['name']
    SM = [i for i, value in enumerate(intensity) if value == 'UNK']
    IN = [i for i, value in enumerate(intensity) if value == 'DERIVED']
    sm_station_lons = [stationdata['lon'][j] for j in SM]
    sm_station_lats = [stationdata['lat'][j] for j in SM]
    in_station_lats = [stationdata['lat'][j] for j in IN]
    in_station_lons = [stationdata['lon'][j] for j in IN]
            

    xmin,xmax,ymin,ymax = shakemap.getRange()
    geodict = shakemap.getGeoDict()
    xmin -= (geodict['xdim']*2)
    xmax += (geodict['xdim']*2)
    ymin -= (geodict['ydim']*2)
    ymax += (geodict['ydim']*2)
    topo = GMTGrid(topofile,bounds=(xmin,xmax,ymin,ymax))
    topo.interpolateToGrid(geodict)
    clear_color = [0,0,0,0.0]

    CORdata = np.flipud(out['cor'].copy())
    ACCUMdata = np.flipud(ACCUM_ARRAY.copy())
    DATAdata = np.flipud(variables['data'].copy())
    DATA_NEWdata = np.flipud(out['data_new'].copy())
    topodata = np.flipud(topo.griddata.copy())
    iwater = np.where(topodata < 0)

    CORdata[iwater] = 0.0101
    ACCUMdata[iwater] = 0.0101
    DATAdata[iwater] = 0.0101
    DATA_NEWdata[iwater] = 0.0101
    
    palette = cm.jet
    COR_masked = np.ma.masked_equal(CORdata, 0.0101)
    ACCUM_masked = np.ma.masked_equal(ACCUMdata, 0.0101)
    DATA_masked = np.ma.masked_equal(DATAdata, 0.0101)
    DATA_NEW_masked = np.ma.masked_equal(DATA_NEWdata, 0.0101)
    
    clat = ymin + (ymax-ymin)/2.0
    clon = xmin + (xmax-xmin)/2.0
    m = Basemap(llcrnrlon=xmin,llcrnrlat=ymin,urcrnrlon=xmax,urcrnrlat=ymax,\
                    rsphere=(6378137.00,6356752.3142),\
                    resolution='i',area_thresh=1000.,projection='lcc',\
                    lat_1=clat,lon_0=clon)

    fig = plt.figure(figsize = (6,6))
    ax = fig.add_axes([0,0,1.0,1.0])
    shakemappable = m.imshow(COR_masked,cmap=palette)
    plt.autoscale(False)
    sta_x,sta_y = m(sm_station_lons, sm_station_lats)
    IN_x, IN_y = m(in_station_lons, in_station_lats)
    m.plot(IN_x, IN_y, 'g>', markersize = 6)
    m.plot(sta_x, sta_y, 'r^', markersize=6)

    m.drawparallels(np.arange(np.min(variables['location_lat_g']),np.max(variables['location_lat_g']),1),labels=[1,0,0,0], linewidth=0.0)
    m.drawmeridians(np.arange(np.min(variables['location_lon_g']),np.max(variables['location_lon_g']),1),labels=[0,0,0,1], linewidth=0.0)
    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    th = plt.title('Correlation Matrix for %s - %s M%.1f, (epsilon)' % (locstr,datestr,mag))
    ch=m.colorbar(mappable=shakemappable)
    plt.show()

    fig = plt.figure(figsize = (6,6))
    ax = fig.add_axes([0,0,1.0,1.0])
    shakemappable = m.imshow(ACCUM_masked,cmap=palette)
    plt.autoscale(False)
    sta_x,sta_y = m(sm_station_lons, sm_station_lats)
    IN_x, IN_y = m(in_station_lons, in_station_lats)
    m.plot(IN_x, IN_y, 'g>', markersize= 6)
    m.plot(sta_x, sta_y, 'r^', markersize=6)

    m.drawparallels(np.arange(np.min(variables['location_lat_g']),np.max(variables['location_lat_g']),1),labels=[1,0,0,0], linewidth=0.0)
    m.drawmeridians(np.arange(np.min(variables['location_lon_g']),np.max(variables['location_lon_g']),1),labels=[0,0,0,1], linewidth=0.0)
    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    th = plt.title('Average Adjusted Matrix for %s - %s M%.1f, (epsilon)' % (locstr,datestr,mag))
    ch=m.colorbar(mappable=shakemappable)
    
    
    fig = plt.figure(figsize = (6,6))
    ax = fig.add_axes([0,0,1.0,1.0])
    shakemappable = m.imshow(DATA_masked,cmap=palette)
    plt.autoscale(False)
    sta_x,sta_y = m(sm_station_lons, sm_station_lats)
    IN_x, IN_y = m(in_station_lons, in_station_lats)
    m.plot(IN_x, IN_y, 'g>', markersize= 6)
    m.plot(sta_x, sta_y, 'r^', markersize=6)

#plt.clim(0,maxdata)
    m.drawparallels(np.arange(np.min(variables['location_lat_g']),np.max(variables['location_lat_g']),1),labels=[1,0,0,0], linewidth=0.0)
    m.drawmeridians(np.arange(np.min(variables['location_lon_g']),np.max(variables['location_lon_g']),1),labels=[0,0,0,1], linewidth=0.0)
    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    th = plt.title('Shakemap %s for %s - %s M%.1f, (pctg)' % (voi,locstr,datestr,mag))
    ch=m.colorbar(mappable=shakemappable)
    
    
    fig = plt.figure(figsize = (6,6))
    ax = fig.add_axes([0,0,1.0,1.0])
    shakemappable = m.imshow(DATA_NEW_masked,cmap=palette)
    plt.autoscale(False)
    sta_x,sta_y = m(sm_station_lons, sm_station_lats)
    IN_x, IN_y = m(in_station_lons, in_station_lats)
    m.plot(IN_x, IN_y, 'g>', markersize= 6)
    m.plot(sta_x, sta_y, 'r^', markersize=6)

#plt.clim(0,maxdata)
    m.drawparallels(np.arange(np.min(variables['location_lat_g']),np.max(variables['location_lat_g']),1),labels=[1,0,0,0], linewidth=0.0)
    m.drawmeridians(np.arange(np.min(variables['location_lon_g']),np.max(variables['location_lon_g']),1),labels=[0,0,0,1], linewidth=0.0)
    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    th = plt.title('Adjusted %s for %s - %s M%.1f, (pctg)' % (voi,locstr,datestr,mag))
    ch=m.colorbar(mappable=shakemappable)
    return
