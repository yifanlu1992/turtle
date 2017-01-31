# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 12:27:49 2017
get the one turtle ILD by time
@author: yifan
"""
import matplotlib as mpl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from turtleModule import str2ndlist,np_datetime
from datetime import  timedelta

obsData = pd.read_csv('ctdWithModTempByDepth.csv') 
tf_index = np.where(obsData['TF'].notnull())[0]    # get the index of good data
obsturtle_id=pd.Series(obsData['PTT'][tf_index],index=tf_index)

indx=[]  # this indx is to get the specifical turtle all index in obsData ,if we use the "where" function ,we just get the length  of tf_index.
for i in tf_index:
    if obsturtle_id[i]==129777:   #we can change the turtle id we interest 
        indx.append(i)
Temp = pd.Series(str2ndlist(obsData['TEMP_VALS'][indx]), index=indx)
Time = pd.Series(np_datetime(obsData['END_DATE'][indx]), index=indx)
Time.sort()        
Indx=[]
for i in Time.index:
    if len(Temp[i])==10:
        Indx.append(i)
        
Data = obsData.ix[Indx]                  
obsTime = pd.Series(np_datetime(Data['END_DATE'].values), index=Indx)
obsTemp = pd.Series(str2ndlist(Data['TEMP_VALS'].values), index=Indx)
obsDepth = pd.Series(str2ndlist(Data['TEMP_DBAR'].values), index=Indx)
Index=range(len(Indx))
obsTime.index=Index
obsTemp.index=Index
obsDepth.index=Index

ILD=[] 
for i in obsTime.index:
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
    ILD.append(ild)
fig=plt.figure() 
ax = fig.add_subplot(111)   
plt.plot(obsTime,ILD,'b',linewidth=2)
plt.ylim([40, 2])
dates = mpl.dates.drange(np.amin(obsTime), np.max(obsTime), timedelta(days=30))
dateFmt = mpl.dates.DateFormatter('%b')
plt.xticks(dates,fontsize=10)    
ax.xaxis.set_major_formatter(dateFmt)
plt.xlabel('time', fontsize=10)
plt.ylabel('isothermal layer depth', fontsize=10)
plt.yticks(fontsize=10)
#plt.legend(loc='lower right')
plt.title('ILDvstime(129777)',fontsize=14)
#plt.text(1,0,'time:'+str(obsTime[i])+'')
#plt.text(1,1,'location:'+str(round(obsLon[i],2))+', '+str(round(obsLat[i],2))+'')
plt.savefig('129777_ILDvstime.png')
plt.show()
