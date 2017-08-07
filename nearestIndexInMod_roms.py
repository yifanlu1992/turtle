'''
Extract new data file named "ctd_good.csv "with new column "modNearestIndex" and "modDepthLayer"
'''
import pandas as pd
import numpy as np
import netCDF4
from datetime import datetime, timedelta
from turtleModule import str2ndlist, dist
import watertempModule as wtm   # a module that has all classes using ROMS and FVCOm model.
def nearest_point_index2(lon, lat, lons, lats):
    d = dist(lon, lat, lons ,lats)
    min_dist = np.min(d)
    index = np.where(d==min_dist)
    return index
obsData = pd.read_csv('ctd_extract_good.csv', index_col=0)
obsLon = obsData['LON']
obsLat = obsData['LAT']

starttime = datetime(2013,07,10) # starttime and endtime can be any time that included by model, we just want a url to get "lon_rho", "lat_rho", "h", "s_rho" in model.
endtime = starttime + timedelta(hours=1)
tempObj = wtm.water_roms()
url = tempObj.get_url(starttime, endtime)
modData = netCDF4.Dataset(url)
modLons = modData.variables['lon_rho'][:]
modLats = modData.variables['lat_rho'][:]
s_rho = modData.variables['s_rho'][:]
h = modData.variables['h'][:]
indexNotNull = obsLon[obsLon.isnull()==False].index # some obslat and obslon of point are empty, get rid of them.
                                                    # or this line can be the indices of TF which is less.
                                                    # indexTF = np.where(obsData['TF'].notnull())[0]

loc = []
for i in indexNotNull:
    ind = []
    lon = obsData['LON'][i]
    lat = obsData['LAT'][i]
    index = nearest_point_index2(lon, lat, modLons, modLats)
    ind.append(index[0][0])
    ind.append(index[1][0])
    loc.append(ind)
    print i
loc = pd.Series(loc, index=indexNotNull)
obsData['modNearestIndex'] = loc #add loc to obsData in case want to save it.

obsDepth = pd.Series(str2ndlist(obsData['TEMP_DBAR']), index=obsData.index)
layersAll = []
for i in indexNotNull:
    nearest_index = loc[i]
    layers = []
    depthLayers = h[nearest_index[0], nearest_index[1]] * s_rho
    for j in range(len(obsDepth[i])):
        # depthLayers = h[nearest_index[0], nearest_index[1]] * s_rho
        l = np.argmin(abs(depthLayers+obsDepth[i][j])) # obsDepth is positive and depthLayers is negitive. So the index of min sum is the layer
        layers.append(l)
        print i, j, l
    layersAll.append(layers)
layersAll = pd.Series(layersAll, index=indexNotNull)
obsData['modDepthLayer'] = layersAll
obsData.to_csv('ctd_good.csv')
