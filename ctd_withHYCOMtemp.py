# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 13:57:00 2015

@author: zhoabin
"""
'get temperature of HYCOM'
import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime,timedelta
from turtleModule import str2ndlist,np_datetime
def getHYcom(latpt,lonpt,time,depth):
    '''old version
    lonpt=360+lonpt+74.160034         #this is HYCOM`S LON range
    if time.month==1 and time.day==1:
        url='http://tds.hycom.org/thredds/dodsC/datasets/global/GLBa0.08_rect/data/temp/rarchv.'+str(time.year)+'_001_00_3zt.nc' 
    else:
        url='http://tds.hycom.org/thredds/dodsC/datasets/global/GLBa0.08_rect/data/temp/rarchv.'+str(time.year)+'_'+format(int(str(time-datetime(time.year,1,1,)).split(' ')[0])+1,'03')+'_00_3zt.nc'
        #use different days setting url
    '''   
    url='http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1/3hrly'
    nc = netCDF4.Dataset(url)
    lat = nc.variables['lat'][1400:1550] 
    lon = nc.variables['lon'][1250:1450]  # ROMS range 
    Depth = nc.variables['depth'][:]
    t=nc.variables['time'][:]
    layer=[]
    for i in depth[0:len(depth)+1]:
        layer.append(np.argmin(abs(i-Depth)))     #get layer of depth
    dist_sq=[]
    indlat,indlon=[],[]
    for k in range(len(lat)):
        for j in range(len(lon)):
            dist_sq.append((lat[k]-latpt)**2 + (lon[j]-lonpt)**2)  # find squared distance of every point on grid
            indlat.append(k)
            indlon.append(j)
    id = np.array(dist_sq).argmin() # 1D index of minimum dist_sq element
    t_diff=(time-datetime(2000,1,1)).total_seconds()/3600     #2000,01,01 is HYCOM`s start time.Unit is hour.
    TIME=np.argmin(abs(t_diff-np.array(t)))
    var=[]
    for i in layer:
        t = nc.variables['water_temp'][TIME,i,1400+indlat[id],1250+indlon[id]] #creates a "netCDF4 object"
        if t<100:   #change bad data to float.
            var.append(t)
        else:
            var.append(-100)
    nc.close()
    return var
##############################################################

obsData = pd.read_csv('ctd_good.csv')
tf_index = np.where(obsData['TF'].notnull())[0]
obsLon, obsLat = obsData['LON'][tf_index], obsData['LAT'][tf_index]
obsTime = pd.Series(np_datetime(obsData['END_DATE'][tf_index]), index=tf_index)
obsDepth = pd.Series(str2ndlist(obsData['TEMP_DBAR'][tf_index]), index=tf_index)
MODtemp=[]
for i in tf_index:
    t=getHYcom(obsLat[i],obsLon[i],obsTime[i],obsDepth[i])
    MODtemp.append(t)
    print i
obsData['modtemp_HYCOM']=pd.Series(MODtemp,index=tf_index)
obsData.to_csv('ctd_withHYCOMtemp.csv')
'''
shipData=pd.read_csv('ship06-08_ch.csv')       
shipLon=pd.Series(str2ndlist(shipData['lon'],bracket=True))
shipLat=pd.Series(str2ndlist(shipData['lat'],bracket=True))
shipTime=pd.Series(shipData['time'])
shipDepth=pd.Series(str2ndlist(shipData['depth'],bracket=True))

Modtemp=[]
for i in shipData.index:
    shipTime[i]=datetime.strptime(shipTime[i],"%Y-%m-%d %H:%M:%S")  # change str to datatime
    t=getHYcom(shipLat[i][0],shipLon[i][0],shipTime[i],shipDepth[i])
    Modtemp.append(t)
    print i
shipData['modtemp_HYCOM']=pd.Series(Modtemp,index=shipData.index)
shipData.to_csv('ship_withHYCOMtemp.csv')
'''
'''
Modtemp=[]
i=[1153, 1203, 1202, 336, 336, 336, 336, 331, 331, 1284, 1284, 695, 695, 695, 695, 695, 695, 695, 690, 690, 690, 925, 691, 692, 691, 692, 1067, 1067, 1067, 923, 924, 924, 924, 995, 919, 697, 698, 697, 698, 697, 698, 697, 698, 924, 924, 924, 924, 924, 925, 924, 925, 692, 919, 692, 692, 692, 692, 695, 695, 1265, 1265, 1265, 1265, 322, 322, 924, 925, 924, 925, 850, 850, 850, 850, 850, 850, 850, 690, 690, 690, 690, 1217, 1217, 40, 40, 41, 217, 41, 217, 41, 217, 21, 21, 21, 37, 37, 37, 324, 1160]
I=pd.Series(i).unique()
I.sort()
for i in I:
    shipTime[i]=datetime.strptime(shipTime[i],"%Y-%m-%d %H:%M:%S")  # change str to datatime
    t=getHYcom(shipLat[i][0],shipLon[i][0],shipTime[i],shipDepth[i])
    Modtemp.append(t)
    print i
shipData['modtemp_HYCOM']=pd.Series(Modtemp,index=I)
shipData.to_csv('ship_withHYCOMtemp1.csv')
'''