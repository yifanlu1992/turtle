'''
Extract new data file named "ctd_good_new.csv "with new column 'LayerDepth'and "modDepthLayer",meaning add the depth and layer of the mod
@author:yifan
'''
import pandas as pd
import numpy as np
import netCDF4
from datetime import datetime, timedelta
from turtleModule import str2ndlist, dist,closest_num

def nearest_point_index(lon, lat, lons, lats):
    d = dist(lon, lat, lons ,lats)
    min_dist = np.min(d)
    index = np.where(d==min_dist)
    return index
def get_url(starttime, endtime):
    url_oceantime = 'http://tds.marine.rutgers.edu/thredds/dodsC/roms/espresso/2013_da/his/ESPRESSO_Real-Time_v2_History_Best?time[0:1:31931]'
    oceantime = netCDF4.Dataset(url_oceantime).variables['time'][:]    #if url2006, ocean_time.
    t1 = (starttime - datetime(2013,5,18)).total_seconds()/3600 # for url2006 it's 2006,01,01; for url2013, it's 2013,05,18, and needed to be devide with 3600
    t2 = (endtime - datetime(2013,5,18)).total_seconds()/3600
    index1 = closest_num(t1, oceantime)
    index2 = closest_num(t2, oceantime)
    url = 'http://tds.marine.rutgers.edu/thredds/dodsC/roms/espresso/2013_da/his/ESPRESSO_Real-Time_v2_History_Best?h[0:1:81][0:1:129],s_rho[0:1:35],lon_rho[0:1:81][0:1:129],lat_rho[0:1:81][0:1:129],temp[{0}:1:{1}][0:1:35][0:1:81][0:1:129],time'
    url = url.format(index1, index2)
    return url
    
obsData = pd.read_csv('ctd_extract_good.csv', index_col=0)
obsLon = obsData['LON']
obsLat = obsData['LAT']
starttime = datetime(2013,07,10) # starttime and endtime can be any time that included by model, we just want a url to get "lon_rho", "lat_rho", "h", "s_rho" in model.
endtime = starttime + timedelta(hours=1)
url=get_url(starttime, endtime)
modData = netCDF4.Dataset(url)
modLons = modData.variables['lon_rho'][:]
modLats = modData.variables['lat_rho'][:]
s_rho = modData.variables['s_rho'][:]
h = modData.variables['h'][:]
indexNotNull = obsLon[obsLon.isnull()==False].index # some obslat and obslon of point are empty, get rid of them.
                                                    # or this line can be the indices of TF which is less.
                                                    # indexTF = np.where(obsData['TF'].notnull())[0]

loc = []  # get index of the nearest point from mod 
for i in indexNotNull:
    ind = []
    lon = obsData['LON'][i]
    lat = obsData['LAT'][i]
    index = nearest_point_index(lon, lat, modLons, modLats)
    #print index
    ind.append(index[0][0])
    ind.append(index[1][0])
    loc.append(ind)
    
loc = pd.Series(loc, index=indexNotNull)
obsData['modNearestIndex'] = loc #add loc to obsData in case want to save it.
obsDepth = pd.Series(str2ndlist(obsData['TEMP_DBAR']), index=obsData.index)

layersAll = []
LayerDepth=[]# depth of each layer
for i in indexNotNull:
    nearest_index = loc[i]
    layers = []
    layerdepth = h[nearest_index[0], nearest_index[1]] * s_rho
    for j in range(len(obsDepth[i])):
        l = np.argmin(abs(layerdepth+obsDepth[i][j])) # obsDepth is positive and layerdepth is negitive. So the index of min sum is the layer
        layers.append(l)
        #print i, j, l
    depth=[]
    for k in layers:
        la=-h[nearest_index[0], nearest_index[1]]*s_rho[k]
        depth.append(la)
    layersAll.append(layers)
    LayerDepth.append(depth)
layersAll = pd.Series(layersAll, index=indexNotNull)
LayerDepth = pd.Series(LayerDepth, index=indexNotNull)
obsData['modDepthLayer'] = layersAll
obsData['LayerDepth'] = LayerDepth
obsData.to_csv('ctd_good_new.csv')
