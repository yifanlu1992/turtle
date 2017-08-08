# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 13:48:41 2017
get a file including  detail information of turtle and ship after matched with 10 km and 3days 
@author: yifan
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from turtleModule import mon_alpha2num, np_datetime, dist,str2ndlist,colors,str2list
###################################################################################
r1,r2 = 0,10                # the obs position that has shipboard position within (r) kilometers might be considered as good data.
day = 3                # the obs time that has shipboard time within (day) days might be considered as good data.
obsData=pd.read_csv('ctdWithModTempByDepth.csv',index_col=0)
tf_index=np.where(obsData['TF'].notnull())[0]
obslat = pd.Series(obsData['LAT'][tf_index],index=tf_index)
obslon = pd.Series(obsData['LON'][tf_index],index=tf_index)
obstime = pd.Series(np_datetime(obsData['END_DATE'][tf_index]),index=tf_index)
obsDepth=pd.Series(str2ndlist(obsData['TEMP_DBAR'][tf_index]),index=tf_index)
obstemp=pd.Series(str2ndlist(obsData['TEMP_VALS'][tf_index]),index=tf_index)
turtle_id=obsData['PTT'][tf_index]

shipData=pd.read_csv('ship06-08_MODELtemp.csv',index_col=0)
shipid=shipData['id']
shiplat=shipData['LAT']
shiplon=shipData['LON']
shiptime=pd.Series((datetime.strptime(x, "%Y-%m-%d %H:%M:%S") for x in shipData['time']))
shipdepth=pd.Series(str2ndlist(shipData['depth'],bracket=True))
shiptemp=pd.Series(str2ndlist(shipData['temperature'],bracket=True))

index = []     #index of turtle
indx=[]      #index of shipboard 
for i in tf_index:
    for j in shipData.index:
        l = dist(obslon[i], obslat[i],shiplon[j],shiplat[j])
        if l<r2 and l>=r1:
            #print l        #distance
            maxtime = obstime[i]+timedelta(days=day)
            mintime = obstime[i]-timedelta(days=day)
            mx = shiptime[j]<maxtime
            mn = shiptime[j]>mintime
            TF = mx*mn  
            if TF==1:      #time
                index.append(i)     #turtle index
                indx.append(j)      #ship index

INDX=pd.Series(indx).unique() 
data=pd.DataFrame(range(len(indx)))
s_id,t_id,s_time,t_time,s_lat,s_lon,t_lat,t_lon=[],[],[],[],[],[],[],[]
s_depth,s_temp,t_depth,t_temp=[],[],[],[]
for i in INDX:
    for j in range(len(indx)):
        if indx[j]==i:
            s=indx[j]
            t=index[j]
            s_id.append(shipid[s])
            s_time.append(shipData['time'][s])
            s_lat.append(shiplat[s])
            s_lon.append(shiplon[s])
            s_depth.append(shipdepth[s])
            s_temp.append(shiptemp[s])
            t_id.append(turtle_id[t])
            t_time.append(obsData['END_DATE'][t])
            t_lat.append(obslat[t])
            t_lon.append(obslon[t])
            t_depth.append(obsDepth[t])
            t_temp.append(obstemp[t])
data['ship_id']=pd.Series(s_id)
data['ship_time']=pd.Series(s_time)
data['ship_lat']=pd.Series(s_lat)
data['ship_lon']=pd.Series(s_lon)
data['ship_depth']=pd.Series(s_depth)
data['ship_temp']=pd.Series(s_temp)
data['turtle_id']=pd.Series(t_id)
data['turtle_time']=pd.Series(t_time)
data['turtle_lat']=pd.Series(t_lat)
data['turtle_lon']=pd.Series(t_lon)
data['turtle_depth']=pd.Series(t_depth)
data['turtle_temp']=pd.Series(t_temp)
data.to_csv('matched_turtleVSship.csv') 
