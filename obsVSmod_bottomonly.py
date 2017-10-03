# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 10:24:20 2017

@author: zdong
"""

'''
#draw the correlation of the deepest observation(we assume it's the bottom of sea) and appropriate model data.
'''
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import netCDF4
import matplotlib.pyplot as plt
from matplotlib import path
from scipy import stats
from watertempModule import *
from turtleModule import mon_alpha2num, np_datetime, bottom_value
FONTSIZE = 25
obs = pd.read_csv('ctd_good.csv')
tf_index = np.where(obs['TF'].notnull())[0] #indices of True data.
obstime = pd.Series(np_datetime(obs['END_DATE'][tf_index]), index=tf_index)
tf_indexs=[]
for i in tf_index:
    if obstime[i]>=datetime(2009,10,11,2,0):
        tf_indexs.append(i)
        
obsLon, obsLat = obs['LON'][tf_indexs], obs['LAT'][tf_indexs]
obsTime = pd.Series(np_datetime(obs['END_DATE'][tf_indexs]), index=tf_indexs)
obsTemp = pd.Series(bottom_value(obs['TEMP_VALS'][tf_indexs]), index=tf_indexs)
obsDepth = obs['MAX_DBAR'][tf_indexs]
#obsdata = pd.DataFrame({'depth':obsdepth, 'temp':obstemp, 'lon':obslon,
                       # 'lat':obslat, 'time':obstime}).sort_index(by='depth')
time=[]
for i in obsTime.index:
    time.append(obsTime[i])
starttime = datetime(2009,10,11,2,0)#datetime(2009,8,24)
endtime = datetime(2013,12,13)


tempobj = water_roms()
url = tempobj.get_url(starttime, endtime)
temp = tempobj.watertemp(obsLon.values, obsLat.values, obsDepth.values, time, url)
temp = pd.Series(temp, index=tf_indexs)
i = temp[temp.isnull()==False].index

indexDeep = obs['MAX_DBAR'][i].index#[obs['MAX_DBAR'][i]>50]
indexShallow = obs['MAX_DBAR'][i][obs['MAX_DBAR'][i]<=50].index

tempModelDeep = temp[indexDeep]
tempModelShallow = temp[indexShallow]
tempObsDeep = obsTemp[indexDeep]
tempObsShallow = obsTemp[indexShallow]

x = np.arange(0.0, 30.0, 0.01)
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.scatter(tempObsDeep, tempModelDeep, s=15, c='b')
ax1.plot(x, x, 'r-')
fit1 = np.polyfit(tempObsDeep, tempModelDeep, 1)
fit_fn1 = np.poly1d(fit1)
ax1.plot(tempObsDeep, fit_fn1(tempObsDeep), 'y-')
plt.ylim([-5,35])
gradient1, intercept1, r_value1, p_value1, std_err1    = stats.linregress(tempObsDeep, tempModelDeep)
ax1.set_title('Deepest bottom, R-squard: %.4f' % r_value1**2, fontsize=18)
ax1.set_ylabel('Model temp', fontsize=16)
ax1.set_xlabel('OBS temp', fontsize=16)
plt.savefig('obsVSmodelBottomTemp1.png', dpi=200)
'''
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.scatter(tempObsShallow, tempModelShallow, s=50, c='b')
ax2.plot(x, x, 'r-')
fit2 = np.polyfit(tempObsShallow, tempModelShallow, 1)
fit_fn2 = np.poly1d(fit2)
ax2.plot(tempObsShallow, fit_fn2(tempObsShallow), 'y-')
gradient2, intercept2, r_value2, p_value2, std_err2\
    = stats.linregress(tempObsShallow, tempModelShallow)
ax2.set_title('Deepest bottom & <50m, R-squard: %.4f' % r_value2**2, fontsize=FONTSIZE)
ax2.set_ylabel('Model temp', fontsize=FONTSIZE)
ax2.set_xlabel('OBS temp', fontsize=FONTSIZE)
plt.savefig('obsVSmodelBottomTemp2.png', dpi=200)
plt.show()'''