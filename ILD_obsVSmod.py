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
from datetime import  timedelta
from turtleModule import  str2ndlist, np_datetime
obsData = pd.read_csv('ctdWithModTempByDepth.csv')
obs_data=pd.read_csv('ctd_good_new.csv')
tf_index1 = np.where(obsData['TF'].notnull())[0]
obsturtle_id=pd.Series(obsData['PTT'][tf_index1],index=tf_index1)

indx=[]  
for i in tf_index1:
    if obsturtle_id[i]==117170:   # we can get each turtle's ILD
        indx.append(i)
Data = obsData.ix[indx]
obs_data=obs_data.ix[indx]                    
obsTime = pd.Series(np_datetime(Data['END_DATE'].values), index=indx)
obsTemp = pd.Series(str2ndlist(Data['TEMP_VALS'].values), index=indx)
obsDepth = pd.Series(str2ndlist(Data['TEMP_DBAR'].values), index=indx)
modDepth = pd.Series(str2ndlist(obs_data['LayerDepth'], bracket=True), index=indx) # If str has '[' and ']', bracket should be True.
modTemp = pd.Series(str2ndlist(Data['modTempByDepth'].values,bracket=True), index=indx)
        
obsILD=[] # get the ILD of the observation
for i in obsTime.index: 
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
modILD=[]   #get the ILD of the mod
for i in obsTime.index:
    try:   
        if modTemp[i][0]==modTemp[i][1]:
            min_slope2=1000  # 1000 is a large of random that represents infinity
        else:
            min_slope2=(modDepth[i][1]-modDepth[i][0])/(modTemp[i][0]-modTemp[i][1])
        for k in range(1,9):   # each profile have the 10 points record
            if modTemp[i][k]==modTemp[i][k+1] or modDepth[i][k+1]==modDepth[i][k]:
                break
            sl2=(modDepth[i][k+1]-modDepth[i][k])/(modTemp[i][k]-modTemp[i][k+1])
            if sl2<=min_slope2:
                min_slope2=sl2
            else:
                break
            ild2=(modDepth[i][k+1]+modDepth[i][k])/2 
       
    except IndexError :
        continue
    modILD.append(ild2)
   
data = pd.DataFrame({'obsTime':obsTime.values, 'obsILD':obsILD, 
                    'modILD': modILD}, index=range(len(obsTime)))
data = data.sort_index(by='obsTime')
Date=[]
for i in data.index:
    Date.append(data['obsTime'][i])
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(Date, data['obsILD'], color='b', linewidth=2,label='observed')
ax.plot(Date, data['modILD'], color='r', linewidth=2, label='modeled')
plt.legend()
ax.set_xlabel('Time', fontsize=14)
ax.set_ylabel('isothermal layer depth', fontsize=14)
dates = mpl.dates.drange(np.amin(obsTime), np.max(obsTime), timedelta(days=30))
dateFmt = mpl.dates.DateFormatter('%b,%Y')
ax.set_xticks(dates)
ax.xaxis.set_major_formatter(dateFmt)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title('ILD for turtle:{0}'.format(117170), fontsize=15)
plt.savefig('ILD_obsVSmod%s.png'%(117170), pdi=200)
plt.show()