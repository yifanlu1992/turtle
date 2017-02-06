# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 11:59:50 2017
the wind of weather and Fvcom
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
example=[129775,129779] # those turtle's id is random  
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
    Data = obsData.ix[indx]
    obs_data=obs_data.ix[indx]                    
    obsTime = pd.Series(np_datetime(Data['END_DATE'].values), index=indx)
    obsTemp = pd.Series(str2ndlist(Data['TEMP_VALS'].values), index=indx)
    obsDepth = pd.Series(str2ndlist(Data['TEMP_DBAR'].values), index=indx)
    obsLon = pd.Series(Data['LON'],index=indx)
    obsLat = pd.Series(Data['LAT'],index=indx)
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
            ild1=(obsDepth[i][k-1]+obsDepth[i][k])/2 # ILD means isothermal layer depth
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
            ild2=(modDepth[i][k-1]+modDepth[i][k])/2 
        except IndexError :
            continue
        modILD.append(ild2)
    
    wind_speed=[]
    wea_wind_speed=[]
    I=[]
    for i in obsTime.index:
        try:
            speed1=get_wind_fvcom(obsTime[i],obsTime[i],obsLat[i],obsLon[i])
            ws1=math.sqrt(speed1[1]**2+speed1[2]**2)
            wind_speed.append(ws1)
            
            speed2=get_weather_wind(obsTime[i].year,obsLat[i],obsLon[i],obsTime[i])
            ws2=math.sqrt(speed2[0]**2+speed2[1]**2)
            wea_wind_speed.append(ws2)
        except IndexError :
            I.append(i)
            continue
    print I
    
    
    data = pd.DataFrame({'obsTime':obsTime.values, 'obsILD':obsILD,'wea_wind_speed':wea_wind_speed, 
                        'modILD': modILD,'wind_speed':wind_speed}, index=range(len(obsTime)))
    data = data.sort_index(by='obsTime')
    data.index=range(len(obsTime))
    Date=[]
    for i in data.index:
        Date.append(data['obsTime'][i])
    ave_obs=round(np.mean(obsILD),1)
    ave_mod=round(np.mean(modILD),1)
    ave_wind_speed=round(np.mean(wind_speed),1)
    ave_wea_wind_speed=round(np.mean(wea_wind_speed),1)
    t=[data,Date,ave_obs,ave_mod,ave_wind_speed,ave_wea_wind_speed]
    T.append(t)
    print 'e'

M,O,W1,W2=[],[],[],[]# smooth the model and observation ILD 
for i in range(2):
    num=4 # smooth by 4 is the best 
    ild1_smooth=smooth(T[i][0]['modILD'],num,'hanning')
    difflen1=len(ild1_smooth)-len(T[i][0]['modILD'])
    ilds1=ild1_smooth[difflen1/2:-difflen1/2]
    
    ild2_smooth=smooth(T[i][0]['obsILD'],num,'hanning')
    difflen2=len(ild2_smooth)-len(T[i][0]['obsILD'])
    ilds2=ild2_smooth[difflen2/2:-difflen2/2]
    
    wind_smooth=smooth(T[i][0]['wind_speed'],num,'hanning')
    difflen3=len(wind_smooth)-len(T[i][0]['wind_speed'])
    wind1=wind_smooth[difflen3/2:-difflen3/2]
    
    wea_wind_smooth=smooth(T[i][0]['wea_wind_speed'],num,'hanning')
    difflen4=len(wea_wind_smooth)-len(T[i][0]['wea_wind_speed'])
    wind2=wea_wind_smooth[difflen4/2:-difflen4/2]
    
    M.append(ilds1)
    O.append(ilds2)
    W1.append(wind1)
    W2.append(wind2)

fig=plt.figure()
for i in range(2):
    ax1 = fig.add_subplot(2,1,i+1,)
    ax1.plot(T[i][1], O[i], color='b', linewidth=1)#when we want to plot the smoothed 'obsILD' ,"T[i][0]['obsILD']" changed with "O[i]' 
    ax1.plot(T[i][1], M[i], color='r', linewidth=1)#when we want to plot the smoothed 'modILD' ,"T[i][0]['modILD']" changed with "M[i]' 
    ax1.set_title('%s'%(example[i]), fontsize=8)
    if i==2:
       dates = mpl.dates.drange(np.amin(T[i][0]['obsTime']), np.max(T[i][0]['obsTime']+timedelta(days=30)), timedelta(days=30))
    else:
       dates = mpl.dates.drange(np.amin(T[i][0]['obsTime']), np.max(T[i][0]['obsTime']), timedelta(days=30))
    dateFmt = mpl.dates.DateFormatter('%b')
    ax1.set_xticks(dates)
    ax1.xaxis.set_major_formatter(dateFmt)
    ax1.set_ylim([35,2]) 
    plt.text(T[i][0]['obsTime'][1],31,'obsmean: '+str(T[i][2]))
    plt.text(T[i][0]['obsTime'][1],34,'modmean: '+str(T[i][3]))
    #plt.text(T[i][0]['obsTime'][1],34,'windmean: '+str(T[i][4]))
    #if i==1:
        #plt.setp(ax1.get_yticklabels() ,visible=False)
        #plt.setp(ax1.get_ylabel() ,visible=False)

    if i==0 :
        plt.setp(ax1.get_xticklabels() ,visible=False)
    ax2=ax1.twinx()
    ax2.plot(T[i][1], W1[i], color='g', linewidth=1)
    ax2.plot(T[i][1], W2[i], color='m', linewidth=1)
    ax2.set_ylim([15,-10]) 
    #if i==0:
        #plt.setp(ax2.get_yticklabels() ,visible=False)
fig.text(0.5, 0.04, '2013', ha='center', va='center', fontsize=14)#  0.5 ,0.04 represent the  plotting scale of x_axis and y_axis
fig.text(0.06, 0.5, 'Isothermal Layer Depth(m)', ha='center', va='center', rotation='vertical',fontsize=12)
fig.text(0.98, 0.5, 'Wind Speed(m/s)', ha='center', va='center', rotation='vertical',fontsize=12)
print 'Fvcom_windmean: '+str(T[i][4])
print 'Weather_windmean: '+str(T[i][5])
plt.savefig('obsVSmodVSweatherVSfvcom.png',dpi=200)
plt.show()


