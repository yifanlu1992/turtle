# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 15:35:00 2017
compare  isothermal layer depth(ILD) of 4 turtle data and model data.(turtles are random 4 turtle)
@author: yifan
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import  timedelta
from turtleModule import  str2ndlist, np_datetime

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
       
    data = pd.DataFrame({'obsTime':obsTime.values, 'obsILD':obsILD, 
                        'modILD': modILD}, index=range(len(obsTime)))
    data = data.sort_index(by='obsTime')
    data.index=range(len(obsTime))
    Date=[]
    for i in data.index:
        Date.append(data['obsTime'][i])
    ave_obs=round(np.mean(obsILD),1)
    ave_mod=round(np.mean(modILD),1)
    t=[data,Date,ave_obs,ave_mod]
    T.append(t)
    print e
    
fig=plt.figure()
for i in range(4):
    ax = fig.add_subplot(2,2,i+1,)
    ax.plot(T[i][1], T[i][0]['obsILD'], color='b', linewidth=1)
    ax.plot(T[i][1], T[i][0]['modILD'], color='r', linewidth=1)
    ax.set_title('%s'%(example[i]), fontsize=8)
    if i==2:
       dates = mpl.dates.drange(np.amin(T[i][0]['obsTime']), np.max(T[i][0]['obsTime']+timedelta(days=30)), timedelta(days=30))
    else:
       dates = mpl.dates.drange(np.amin(T[i][0]['obsTime']), np.max(T[i][0]['obsTime']), timedelta(days=30))
    dateFmt = mpl.dates.DateFormatter('%b')
    ax.set_xticks(dates)
    ax.xaxis.set_major_formatter(dateFmt)
    ax.set_ylim([35,2]) 
    plt.text(T[i][0]['obsTime'][1],31,'obsmean'+str(T[i][2]))
    plt.text(T[i][0]['obsTime'][1],34,'modmean'+str(T[i][3]))
    if i==1 or i==3: 
        plt.setp(ax.get_yticklabels() ,visible=False)
    if i==0 or i==1:
        plt.setp(ax.get_xticklabels() ,visible=False)
fig.text(0.5, 0.04, '2013', ha='center', va='center', fontsize=14)#  0.5 ,0.04 represent the  plotting scale of x_axis and y_axis
fig.text(0.06, 0.5, 'isothermal layer depth', ha='center', va='center', rotation='vertical',fontsize=14)
plt.savefig('ILD_obsVSmod.png',dpi=200)
plt.show()
