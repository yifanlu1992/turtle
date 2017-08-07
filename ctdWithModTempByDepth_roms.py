'''
get the roms model temp of different depth.

'''
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import watertempModule as wtm
from turtleModule import str2ndlist, closest_num, np_datetime, bottom_value, dist

def getModTemp(modTempAll, obsTime, modLayer, modNearestIndex, s_rho, waterDepth, starttime, oceantime):
    '''
    Return model temp based on observation layers or depth
    '''
    print "This program has assign the temp of island and land to 10000"
    indx = closest_num((starttime -datetime(2006,1,1)).total_seconds(), oceantime)
    modTemp = []
    l = len(modLayer.index)
    for i in modLayer.index:    # Loop the model location
        '''
        # For layers
        print i, l, 'getModTemp'
        timeIndex = closest_num((obsTime[i]-datetime(2006,01,01)).total_seconds(), oceantime)-ind
        modTempTime = modTempAll[timeIndex]
        modTempTime[modTempTime.mask] = 10000
        t = np.array([modTempTime[modLayer[i][j],modNearestIndex[i][0], modNearestIndex[i][1]] \
                          for j in range(len(modLayer[i]))])
        modTemp.append(t)
        '''
        # For depth
        print i, l, 'getModTemp'
        timeIndex1 = closest_num((obsTime[i]-datetime(2006,01,01)).total_seconds(), oceantime)
        timeIndex = timeIndex1 - indx
        temp = modTempAll[timeIndex]
        temp[temp.mask] = 10000 # Assign the temp of island and land to 10000
        a, b = int(modNearestIndex[i][0]), int(modNearestIndex[i][1]) # index of nearest model node
        t = []
        for depth in obsDepth[i]:
            depth = -depth
            locDepth = waterDepth[a, b]# Get the bottom depth of this location. waterDepth is 'h'
            lyrDepth = s_rho * locDepth# Depth of each layer
            if depth > lyrDepth[-1]: # Obs is shallower than last layer which is the surface.
                Temp = (temp[-2,a,b]-temp[-1,a,b])/(lyrDepth[-2]-lyrDepth[-1]) * \
                    (depth-lyrDepth[-1]) + temp[-1,a,b]
            elif depth < lyrDepth[0]: # Obs is deeper than first layer which is the bottom.
                Temp = (temp[1,a,b]-temp[0,a,b])/(lyrDepth[1]-lyrDepth[0]) * \
                    (depth-lyrDepth[0]) + temp[0,a,b]
            else:
                ind = closest_num(depth, lyrDepth)
                Temp = (temp[ind,a,b]-temp[ind-1,a,b])/(lyrDepth[ind]-lyrDepth[ind-1]) * \
                    (depth-lyrDepth[ind-1]) + temp[ind-1,a,b]
            t.append(Temp)
        modTemp.append(t)
    modTemp = np.array(modTemp)
    return modTemp
FONTSIZE = 25
obsData = pd.read_csv('ctd_good.csv')
tf_index = np.where(obsData['TF'].notnull())[0]
obsLon, obsLat = obsData['LON'][tf_index], obsData['LAT'][tf_index]
obsTime = pd.Series(np_datetime(obsData['END_DATE'][tf_index]), index=tf_index)
obsTemp = pd.Series(str2ndlist(obsData['TEMP_VALS'][tf_index]), index=tf_index)
# obsTemp = pd.Series(bottom_value(obs['TEMP_VALS'][tf_index]), index=tf_index)
obsDepth = pd.Series(str2ndlist(obsData['TEMP_DBAR'][tf_index]), index=tf_index)
modLayer = pd.Series(str2ndlist(obsData['modDepthLayer'][tf_index],bracket=True), index=tf_index)
modNearestIndex = pd.Series(str2ndlist(obsData['modNearestIndex'][tf_index], bracket=True), index=tf_index)

starttime = datetime(2009, 8, 24)
endtime = datetime(2013, 12, 13)
tempObj = wtm.waterCTD()
url = tempObj.get_url(starttime, endtime)
# modTemp1 = tempObj.watertemp(obsLon.values, obsLat.values, obsDepth.values, obsTime.values, url)
modDataAll = tempObj.get_data(url)
oceantime = modDataAll['ocean_time']
modTempAll = modDataAll['temp']
s_rho = modDataAll['s_rho']
waterDepth = modDataAll['h']
modTemp = getModTemp(modTempAll, obsTime, modLayer, modNearestIndex, s_rho, waterDepth, starttime, oceantime)
obsData['modTempByDepth'] = pd.Series(modTemp, index = tf_index)
obsData.to_csv('ctdWithModTempByDepth.csv')
