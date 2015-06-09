import cartopy
import matplotlib.pyplot as plt
import numpy as np
import scipy
from neicio.gmt import GMTGrid
from matplotlib import cm

WATER_COLOR = [.47,.60,.81]


def plot(out, variables, voi, shakemap, stationdata, ACCUM_ARRAY):
    maxdata = np.amax(out['data_new'])
    attributes = shakemap.getAttributes()
    intensity = stationdata['name']
    SM = [i for i, value in enumerate(intensity) if value == 'UNK']
    IN = [i for i, value in enumerate(intensity) if value == 'DERIVED']
    sm_station_lons = [stationdata['lon'][j] for j in SM]
    sm_station_lats = [stationdata['lat'][j] for j in SM]
    in_station_lats = [stationdata['lat'][j] for j in IN]
    in_station_lons = [stationdata['lon'][j] for j in IN]            
    palette = cm.jet

    fig = plt.figure(figsize=(10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    cartopy.feature.COASTLINE.scale = '50m'
    cartopy.feature.LAND.scale = '50m'
    cartopy.feature.OCEAN.scale = '50m'
    ax.add_feature(cartopy.feature.OCEAN,facecolor=WATER_COLOR)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':',zorder=10)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(shakemap.getRange())
    map = ax.imshow(out['cor'],extent=shakemap.getRange(), origin='upper',cmap=palette)
    plt.plot(sm_station_lons, sm_station_lats, 'g>', markersize = 6)
    plt.plot(in_station_lons, in_station_lats, 'r^', markersize = 6)
    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    th = plt.title('Correlation Matrix for %s - %s M%.1f, (epsilon)' % (locstr,datestr,mag), y = 1.08)
    ch=plt.colorbar(map, shrink=0.7)
    plt.show(map)

    fig = plt.figure(figsize = (10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    cartopy.feature.COASTLINE.scale = '50m'
    cartopy.feature.LAND.scale = '50m'
    cartopy.feature.OCEAN.scale = '50m'
    ax.add_feature(cartopy.feature.OCEAN,facecolor=WATER_COLOR)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':',zorder=10)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(shakemap.getRange())
    map = ax.imshow(out['data'],extent=shakemap.getRange(), origin='upper',cmap=palette)
    plt.plot(sm_station_lons, sm_station_lats, 'g>', markersize = 6)
    plt.plot(in_station_lons, in_station_lats, 'r^', markersize = 6)
    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    th = plt.title('ShakeMap for %s - %s M%.1f, (epsilon)' % (locstr,datestr,mag), y = 1.08)
    ch=plt.colorbar(map, shrink=0.7)
    plt.show(map)

    fig = plt.figure(figsize = (10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    cartopy.feature.COASTLINE.scale = '50m'
    cartopy.feature.LAND.scale = '50m'
    cartopy.feature.OCEAN.scale = '50m'
    ax.add_feature(cartopy.feature.OCEAN,facecolor=WATER_COLOR)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':',zorder=10)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(shakemap.getRange())
    map = ax.imshow(out['data_new'],extent=shakemap.getRange(), origin='upper',cmap=palette)
    plt.plot(sm_station_lons, sm_station_lats, 'g>', markersize = 6)
    plt.plot(in_station_lons, in_station_lats, 'r^', markersize = 6)
    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    th = plt.title('Avg Adj Matrix for %s - %s M%.1f, (epsilon)' % (locstr,datestr,mag), y = 1.08)
    ch=plt.colorbar(map, shrink=0.7)
    plt.show(map)


    fig = plt.figure(figsize = (10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    cartopy.feature.COASTLINE.scale = '50m'
    cartopy.feature.LAND.scale = '50m'
    cartopy.feature.OCEAN.scale = '50m'
    ax.add_feature(cartopy.feature.OCEAN,facecolor=WATER_COLOR)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':',zorder=10)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(shakemap.getRange())
    map = ax.imshow(ACCUM_ARRAY,extent=shakemap.getRange(), origin='upper',cmap=palette)
    plt.plot(sm_station_lons, sm_station_lats, 'g>', markersize = 6)
    plt.plot(in_station_lons, in_station_lats, 'r^', markersize = 6)
    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    th = plt.title('Accumulated Matrix for %s - %s M%.1f, (epsilon)' % (locstr,datestr,mag), y = 1.08)
    ch=plt.colorbar(map, shrink=0.7)
    plt.show(map)

    return
