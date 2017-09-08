# -*- coding: utf-8 -*-
"""
Created on Tue Jan  6 14:59:03 2015
I use the Excel to convert the "ctd_data.dat" to the "ctd_data.csv"
change shipboard data:let data which is in one location put in one line
@author: yifan
"""
import pandas as pd 
import numpy as np
from datetime import datetime,timedelta
from turtleModule import mon_alpha2num
###############################################

def np_datetimes(m):
    'change shipdata time to datetime'
    dt = []
    for i in m:
        year = int('20'+i[7:])
        month = mon_alpha2num(i[3:6])
        day =  int(i[:2])
        temp = datetime(year,month,day)
        dt.append(temp)
        
    dt = np.array(dt)
    print len(dt)
    
    return dt
    

shipdata=pd.read_csv('ctd_data.csv')
day=pd.Series(np_datetimes(shipdata['gmt_date']))
shiptime=pd.Series(shipdata['gmt_time'])
time=[]
for i in shipdata.index:
    minute=round((shiptime[i]-int(shiptime[i]))*60,2)
    second=int((minute-int(minute))*60)
    hour=shiptime[i]/1
    time.append(day[i]+timedelta(hours=hour)+timedelta(minutes=int(minute))+timedelta(seconds=second))    
Time=pd.Series(time)
shipdata['Time']=pd.Series(Time,index=shipdata.index)   #change gmt to %Y,%M,%D %H:%M:%S
deepest=pd.Series(shipdata['bottom_depth'])
ids=pd.Series(shipdata['cruise_id'])
lat=pd.Series(shipdata['lat'])
lon=pd.Series(shipdata['lon'])
depth=pd.Series(shipdata['press'])
temp=pd.Series(shipdata['temp'])
num=[0]
for i in shipdata.index:
    if i+1>shipdata.index[-1]:
        break
    if Time[i]<>Time[i+1]:
        num.append(i+1)           # get the beginning index of next location
num.append(len(shipdata.index))    #use for the last line
Depths,Temps,IDs,bottom_depths,Lats,Lons=[],[],[],[],[],[]   #use for putting one location data in one line                                   
for i in range(len(num)):
    if i+1<len(num):
        Depth,Temp,ID,bottom_depth,Lat,Lon=[],[],[],[],[],[]   # use for each line       
        for j in range(num[i],num[i+1]):
            Depth.append(depth[j])
            Temp.append(round(temp[j],2))
            ID.append(ids[j])
            bottom_depth.append(deepest[j])
            Lat.append(lat[j])
            Lon.append(-lon[j])
        Depths.append(Depth)
        Temps.append(Temp)
        IDs.append(pd.Series(ID).unique())
        bottom_depths.append(pd.Series(bottom_depth).unique()) 
        Lats.append(pd.Series(Lat).unique())
        Lons.append(pd.Series(Lon).unique())
time=Time.unique()
Data=pd.DataFrame(range(len(Temps)))
Data['id']=pd.Series(IDs)
Data['bottom_depth']=pd.Series(bottom_depths)
Data['lat']=pd.Series(Lats)
Data['lon']=pd.Series(Lons)
Data['time']=pd.Series(time)
Data['depth']=pd.Series(Depths)
Data['temperature']=pd.Series(Temps)
Data.to_csv('shipdata.csv')  
