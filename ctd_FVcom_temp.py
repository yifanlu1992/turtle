# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 11:23:07 2015

@author: zhaobin
"""
'get temperature from fvcom`s website and create ctd_FVcom_temp.csv or ship_FVcom_temp.csv'
import numpy as np
import pandas as pd
import netCDF4
from turtleModule import str2ndlist
##############################################################
url='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?temp[0:1:316008][0:1:44][0:1:48450]'
modData=netCDF4.Dataset(url)
modtemp=modData.variables['temp'] 

obsData=pd.read_csv('ctd_FVcom_node_time_layer.csv')     #ctd
#obsData=pd.read_csv('ship_FVcom_node_time_layer.csv')     #ship
tf_index=np.where(obsData['in FVcom range'].notnull())[0]
obsnode=pd.Series(obsData['modnode'][tf_index],index=tf_index)
obstime=pd.Series(obsData['modtime'][tf_index],index=tf_index)
obslayer=pd.Series(str2ndlist(obsData['modlayer'][tf_index],bracket=True),index=tf_index)
TEMP=[]

for i in tf_index:
    temp=[]
    for j in obslayer[i]:
        t=modtemp[int(obstime[i])][int(j)][int(obsnode[i])]
        temp.append(t)
    TEMP.append(temp)
    print i
obsData['modtempBYdepth']=pd.Series(TEMP,index=tf_index)
obsData.to_csv('ctd_FVcom_temp.csv')         #ctd
#obsData.to_csv('ship_FVcom_temp.csv')       #ship
    