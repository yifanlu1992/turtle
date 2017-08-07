# -*- coding: utf-8 -*-
"""
Analysis data whether in FVCOM range or not ,calculate each dive`s nearest node,layer and time of model
"""

import pandas as pd
from datetime import datetime,timedelta
import netCDF4
from matplotlib.path import Path as mpPath
from turtleModule import str2ndlist,dist,np_datetime
def Closest_num(num, numlist, i=0):
    '''
    Return index of the closest number in the list
    '''
    index1, index2 = 0, len(numlist)
    indx = int(index2/2)
    if not numlist[0] < num < numlist[-1]:
        if numlist[0]>num:
            i=0
        else:
            i=len(numlist)-1
    else:        
        if index2 == 2:
            l1, l2 = num-numlist[0], numlist[-1]-num
            if l1 < l2:
                i = i
            else:
                i = i+1
        elif num == numlist[indx]:
            i = i + indx
        elif num > numlist[indx]:
            i = Closest_num(num, numlist[indx:],
                          i=i+indx)
        elif num < numlist[indx]:
            i = Closest_num(num, numlist[0:indx+1], i=i)
    return i
def point_in_poly(vertices, x, y):
    path =mpPath(vertices)
    point = (x,y)
    return path.contains_point(point)
#########################################################################
url='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?h[0:1:48450],lat[0:1:48450],lon[0:1:48450],siglay[0:1:44][0:1:48450],time[0:1:316008]'
data=netCDF4.Dataset(url)
h=data.variables['h'][:]
lat=data.variables['lat'][:]
lon=data.variables['lon'][:]
siglay=data.variables['siglay'][:]
time=data.variables['time'][:]
modData=pd.DataFrame(range(len(h)))
modData['lat']=pd.Series(lat)
modData['lon']=pd.Series(lon)
modData['h']=pd.Series(h)                    #get FVcom`s lat,lon,h
modtime=[]
for i in range(len(time)):
    t=timedelta(days=float(time[i])).total_seconds()
    modtime.append(t)                        #change time days to seconds
modTime=pd.Series(modtime)
depth=(-h*siglay).transpose()
print depth                #each layer`s depth
Depth=[]
for i in range(len(depth)):
    d=[]
    for j in depth[i]:
        d.append(round(j,3))
    Depth.append(d)
modData['depthBYlayer']=pd.Series(Depth)


obsData=pd.read_csv('ctd_extract_good.csv')
tf_index=np.where(obsData['TF'].notnull())[0]    #this is turtle data
obslat=pd.Series(obsData['LAT'][tf_index],index=tf_index)
obslon=pd.Series(obsData['LON'][tf_index],index=tf_index)
obsDepth = pd.Series(str2ndlist(obsData['TEMP_DBAR'][tf_index]), index=tf_index)
obsTime=pd.Series(np_datetime(obsData['END_DATE'][tf_index]),index=tf_index)      #this is CTD
'''
obsData=pd.read_csv('ship06-08_ch.csv')       #this is ship data
obslat=pd.Series(str2ndlist(obsData['lat'],bracket=True))
obslon=pd.Series(str2ndlist(obsData['lon'],bracket=True))
obsDepth = pd.Series(str2ndlist(obsData['depth'],bracket=True))   
obsTime=pd.Series(obsData['time'])
for i in range(len(obsTime)):
    obsTime[i]=datetime.strptime(obsTime[i], "%Y-%m-%d %H:%M:%S")  # change str to datatime   #this is ship
'''
pp=[(-75.6127,38.2253) ,(-75.6394,39.8235),(-74.3474,40.7558),(-69.9301,44.6404),(-61.717,46.9046),(-57.1221,41.8435)]            
for q in range(67,99):
    pp.append((modData['lon'][q],modData['lat'][q]))            
              #range of fvcom

FV_index=[]
for i in obslat.index:
    n=point_in_poly(pp,obslon[i],obslat[i])       #turtle
    'n=point_in_poly(pp,obslon[i][0],obslat[i][0])'   #ship
    if n==1:
        #print i
        FV_index.append(i)               #ensure turtle data in range of FVcom 
Depth=[]
node=[]
TIME=[]
for i in FV_index:
    distance=dist(obslon[i],obslat[i],modData['lon'],modData['lat'])
    k=np.argmin(distance)                 #get nearest node
    node.append(k)
    depth=[]
    #print i
    for j in range(len(obsDepth[i])):
        m=Closest_num(obsDepth[i][j],np.array(modData['depthBYlayer'][k]))   #in the nearest node get which layer it belongs to
        depth.append(int(m))
    Depth.append(depth)
    
    t_diff=(obsTime[i]-datetime(1858,11,17)).total_seconds()     #1858,11,17 is FVCOM`s start time
    T=Closest_num(t_diff,np.array(modTime))                      #get the nearest time
    TIME.append(T)
ctd_FV = pd.Series([True]*len(FV_index), index=FV_index)
obsData['in FVcom range'] = ctd_FV
obsData['modnode']=pd.Series(node,index=FV_index)
obsData['modtime']=pd.Series(TIME,index=FV_index)
obsData['modlayer']=pd.Series(Depth,index=FV_index)
obsData.to_csv('ctd_FVcom_node_time_layer.csv')        #CTD

