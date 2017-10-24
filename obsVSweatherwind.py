# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 11:59:50 2017

@author: yifan# -*- coding: utf-8 -*-
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import  timedelta
from turtleModule import  str2ndlist, np_datetime
from utilities import smooth , get_wind_fvcom
from get_ncep_wind_test import get_weather_wind
import math
example=[117170,129775,118905,129779] # those turtle's id is random  
T=[]  # get all information of six example turtle
for e in range(len(example)):
    obsData = pd.read_csv('ctdWithModTempByDepth.csv')
    obs_data=pd.read_csv('ctd_good_new.csv')
    tf_index = np.where(obsData['TF'].notnull())[0]
    obsturtle_id=pd.Series(obsData['PTT'][tf_index],index=tf_index)
    indx=[]  
    for i in tf_index:
        if obsturtle_id[i]==example[e]:   # we can get each turtle's ILD
            indx.append(i)
    Temp = pd.Series(str2ndlist(obsData['TEMP_VALS'][indx]), index=indx)
    Indx=[]
    for i in indx:
        if len(Temp[i])==10:
            Indx.append(i)
    Data = obsData.ix[Indx]
    obs_data=obs_data.ix[Indx]                    
    obsTime = pd.Series(np_datetime(Data['END_DATE'].values), index=Indx)
    obsTemp = pd.Series(str2ndlist(Data['TEMP_VALS'].values), index=Indx)
    obsDepth = pd.Series(str2ndlist(Data['TEMP_DBAR'].values), index=Indx)
    obsLon = pd.Series(Data['LON'],index=Indx)
    obsLat = pd.Series(Data['LAT'],index=Indx)
    modDepth = pd.Series(str2ndlist(obs_data['LayerDepth'], bracket=True), index=Indx) # If str has '[' and ']', bracket should be True.
    modTemp = pd.Series(str2ndlist(Data['modTempByDepth'].values,bracket=True), index=Indx)
    
            
    obsILD=[] # get the ILD of the observation
    for i in Indx: 
        try:   
            if obsTemp[i][0]==obsTemp[i][1]:
                min_slope=1000  # 1000 is a large of random that represents infinity
            else:
                min_slope=abs((obsDepth[i][1]-obsDepth[i][0])/(obsTemp[i][0]-obsTemp[i][1]))
            m=0
            for k in range(1,9):
                if obsTemp[i][k]==obsTemp[i][k+1]:  
                   obsTemp[i][k+1]-=0.000000000001
                sl=abs((obsDepth[i][k+1]-obsDepth[i][k])/(obsTemp[i][k]-obsTemp[i][k+1]))
                if sl==0:
                    sl=1000
                if sl<=min_slope:
                    min_slope=sl
                    m=k
            ild=(obsDepth[i][m+1]+obsDepth[i][m])/2
        except IndexError :
            continue
        obsILD.append(ild)
    
    wea_wind_speed=[]
    I=[]
    for i in obsTime.index:
        try:
            speed2=get_weather_wind(obsTime[i].year,obsLat[i],obsLon[i],obsTime[i])
            ws2=math.sqrt(speed2[0]**2+speed2[1]**2)
            wea_wind_speed.append(ws2)
        except IndexError :
            I.append(i)
            continue
    print I

    data = pd.DataFrame({'obsTime':obsTime.values, 'obsILD':obsILD,'wea_wind_speed':wea_wind_speed, 
                        }, index=range(len(obsTime)))
    data = data.sort_index(by='obsTime')
    data.index=range(len(obsTime))
    Date=[]
    for i in data.index:
        Date.append(data['obsTime'][i])
    ave_obs=round(np.mean(obsILD),1)
    ave_wea_wind_speed=round(np.mean(wea_wind_speed),1)
    t=[data,Date,ave_obs,ave_wea_wind_speed]
    T.append(t)
    print 'e'

O,W2=[],[]# smooth the model and observation ILD 
for i in range(4):
    num=6 # smooth by 6 is the best 
    ild2_smooth=smooth(T[i][0]['obsILD'],num,'hanning')
    difflen2=len(ild2_smooth)-len(T[i][0]['obsILD'])
    ilds2=ild2_smooth[difflen2/2:-difflen2/2]

    wea_wind_smooth=smooth(T[i][0]['wea_wind_speed'],num,'hanning')
    difflen4=len(wea_wind_smooth)-len(T[i][0]['wea_wind_speed'])
    wind2=wea_wind_smooth[difflen4/2:-difflen4/2]
    
    O.append(ilds2)
    W2.append(wind2)

fig=plt.figure()
for i in range(4):
    ax1 = fig.add_subplot(2,2,i+1,)
    ax1.plot(T[i][1], O[i], color='b', linewidth=1,label='obsmean: '+str(T[i][2]))#when we want to plot the smoothed 'obsILD' ,"T[i][0]['obsILD']" changed with "O[i]' 
    ax1.set_title('id: %s'%(example[i]), fontsize=8)
    ax1.set_ylim([35,0])
    plt.legend(loc='upper right',fontsize = 'xx-small')
    if i==1 or i==3: 
        plt.setp(ax1.get_yticklabels() ,visible=False)
    if i==0 or i==1:
        plt.setp(ax1.get_xticklabels() ,visible=False)
    ax2=ax1.twinx()
    ax2.plot(T[i][1], W2[i], color='m', linewidth=1,label='windmean: '+str(T[i][3]))
    ax2.set_ylim([13,-12]) 
    if i==0 or i==2:
        plt.setp(ax2.get_yticklabels() ,visible=False)
    if i==2:
       dates = mpl.dates.drange(np.amin(T[i][0]['obsTime']), np.max(T[i][0]['obsTime']+timedelta(days=30)), timedelta(days=30))
    else:
       dates = mpl.dates.drange(np.amin(T[i][0]['obsTime']), np.max(T[i][0]['obsTime']), timedelta(days=30))
    dateFmt = mpl.dates.DateFormatter('%b')
    ax1.set_xticks(dates)
    ax1.xaxis.set_major_formatter(dateFmt)
    ax2.set_xticks(dates)
    ax2.xaxis.set_major_formatter(dateFmt)

fig.text(0.5, 0.04, '2013', ha='center', va='center', fontsize=14)#  0.5 ,0.04 represent the  plotting scale of x_axis and y_axis
fig.text(0.06, 0.5, 'Isothermal Layer Depth(m)', ha='center', va='center',color='b', rotation='vertical',fontsize=12)
fig.text(0.98, 0.5, 'Wind Speed(m/s)', ha='center', va='center', rotation='vertical',color='m',fontsize=12)
plt.savefig('obsVSweather.png',dpi=200)
plt.show()
