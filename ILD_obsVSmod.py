# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 15:35:00 2017
draw temp change of one specific turtle data and model data.
@author: yifan
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime, timedelta
from turtleModule import str2list, str2ndlist, np_datetime, bottom_value, closest_num, dist

#########################################################
obsData = pd.read_csv('ctdWithModTempByDepth.csv', index_col=0)
tf_index = np.where(obsData['TF'].notnull())[0]
obsturtle_id=pd.Series(obsData['PTT'][tf_index],index=tf_index)

indx=[]  
for i in tf_index:
    if obsturtle_id[i]==118905:   #this turtle is same turtle with 4-second turtle
        indx.append(i)
Data = obsData.ix[indx]
                    
obsTime = pd.Series(np_datetime(Data['END_DATE'].values), index=indx)
obsTemp = pd.Series(str2ndlist(Data['TEMP_VALS'].values), index=indx)
obsDepth = pd.Series(str2ndlist(Data['TEMP_DBAR'].values), index=indx)
layers = pd.Series(str2ndlist(Data['modDepthLayer'], bracket=True), index=indx) # If str has '[' and ']', bracket should be True.
modTemp = pd.Series(str2ndlist(Data['modTempByDepth'].values,bracket=True), index=indx)

obsILD=[]
for i in obsTime.index:
    #if i!= 3 and i!= 44 and i!=50: # these profile have something wrong 
    try:   
        if obsTemp[i][0]==obsTemp[i][1]:
            min_slope1=1000  # 1000 is a large of random that represents infinity
        else:
            min_slope1=(obsDepth[i][1]-obsDepth[i][0])/(obsTemp[i][0]-obsTemp[i][1])
        for k in range(1,9):   # each profile have the 10 points record
            if obsTemp[i][k]==obsTemp[i][k+1] or obsDepth[i][k+1]==obsDepth[i][k]:
                break
            sl1=(obsDepth[i][k+1]-obsDepth[i][k])/(obsTemp[i][k]-obsTemp[i][k+1])
            #print sl
            if sl1<=min_slope1:
                min_slope1=sl1
            else:
                break
            ild1=(obsDepth[i][k+1]+obsDepth[i][k])/2 # ILD means isothermal layer depth
       
    except IndexError :
        continue
    obsILD.append(ild1)
modILD=[]
for i in obsTime.index:
    try:   
        if modTemp[i][0]==modTemp[i][1]:
            min_slope2=1000  # 1000 is a large of random that represents infinity
        else:
            min_slope2=(layers[i][-2]-layers[i][-1])/(modTemp[i][0]-modTemp[i][1])
        for k in range(1,9):   # each profile have the 10 points record
            if modTemp[i][k]==modTemp[i][k+1] or layers[i][k+1]==layers[i][k]:
                break
            sl2=(layers[i][k+1]-layers[i][k])/(modTemp[i][k]-modTemp[i][k+1])
            #print sl
            if sl2<=min_slope2:
                min_slope2=sl2
            else:
                break
            ild2=(layers[i][k+1]+layers[i][k])/2 # ILD means isothermal layer depth
       
    except IndexError :
        continue
    modILD.append(ild2)
   
data = pd.DataFrame({'obsTime':obsTime.values, 'obsILD':obsILD, 
                    'modILD': modILD}, index=range(len(obsTime)))
data = data.sort_index(by='obsTime')
# data['time'] = smooth(data['time'].values, timedelta(days=20))
Date=[]
for i in data.index:
    Date.append(data['obsTime'][i])
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(Date, data['obsILD'], color='b', linewidth=2,label='observed')

ax.plot(Date, data['modILD'], color='r', linewidth=2, label='modeled')

plt.legend()
ax.set_xlabel('Time', fontsize=20)
ax.set_ylabel('isothermal layer depth', fontsize=20)
dates = mpl.dates.drange(np.amin(obsTime), np.max(obsTime), timedelta(days=30))
dateFmt = mpl.dates.DateFormatter('%b,%Y')
ax.set_xticks(dates)
ax.xaxis.set_major_formatter(dateFmt)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.title('Time series of temp for turtle:{0}'.format(118905), fontsize=25)
#plt.savefig('timeSeries.png', pdi=200)
plt.show()


