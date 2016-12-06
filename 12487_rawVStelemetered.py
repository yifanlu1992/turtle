# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 21:20:06 2016
compare the telemetered data with raw upcast data for each profile
@author: yifan
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from turtleModule import str2ndlist,np_datetime
from gettopo import gettopo
def find_start(a,List):
    start=len(List)-3
    if List[a]<2:
        for i in range(a,len(List)-2):
            if List[i]<List[i+1] and List[i]<2:
                if List[i+1]<List[i+2] and List[i+1]>=2 :
                    start=i
                    break                            #find farest on surface before diving
        for i in range(start,len(List)-2):
            if List[i]>=2 and List[i]>List[i+1]:
               if List[i+1]<2 and List[i+1]>List[i+2]:
                   break                             #find nearest on surface after diving
    return [start,i+1]
def closest_time(time, timelist, i=0):
    '''
    Return index of the closest time in the list
    '''
    index = len(timelist)
    indx = int(index/2)
    if not timelist[0] < time < timelist[len(timelist)-1]:
        raise Exception('{0} is not in {1}'.format(str(time), str(timelist)))
    if index == 2:
        l1, l2 = time-timelist[0], timelist[len(timelist)-1]-time
        if l1 < l2:
            i = i
        else:
            i = i+1
    elif time == timelist[indx]:
        i = i + indx
    elif time > timelist[indx]:
        i = closest_time(time, timelist[indx:],
                          i=i+indx)
    elif time < timelist[indx]:
        i = closest_time(time, timelist[0:indx+1], i=i)
    return i
###########################################################################
obsData = pd.read_csv('ctdWithModTempByDepth.csv') 
tf_index = np.where(obsData['TF'].notnull())[0]    # get the index of good data
obsturtle_id=pd.Series(obsData['PTT'][tf_index],index=tf_index)
secondData=pd.read_csv('12487_location.csv')
tf_index1 = np.where(secondData['index'].notnull())[0]
time=pd.Series(secondData['time'],index=tf_index1)
depth=pd.Series(secondData['depth'],index=tf_index1)
temp=pd.Series(secondData['temp'],index=tf_index1)
lon=pd.Series(secondData['lon'],index=tf_index1)
lat=pd.Series(secondData['lat'],index=tf_index1)
inde=pd.Series(secondData['index'],index=tf_index1)
indx=[]
for i in tf_index:
    if obsturtle_id[i]==118905:   #this turtle is same turtle with 4-second turtle
        indx.append(i)
obsLon, obsLat = obsData['LON'][indx], obsData['LAT'][indx]
obsTime = pd.Series(np_datetime(obsData['END_DATE'][indx]), index=indx)
obsTemp = pd.Series(str2ndlist(obsData['TEMP_VALS'][indx]), index=indx)
obsDepth = pd.Series(str2ndlist(obsData['TEMP_DBAR'][indx]), index=indx)
for i in indx:
    waterdepth=-gettopo(obsLat[i],obsLon[i])
    print 'waterdepth: '+ str(waterdepth)
    Index_all=[]   #find indices which are in same area
    for j in tf_index1:
        if i==inde[j]:
           Index_all.append(j)
    dep=pd.Series(depth,index=Index_all)
    newdepth=dep[:]
    newtime=pd.Series(time,index=Index_all)[:]
    newtime=pd.to_datetime(newtime)
    newdepth.index=range(len(newdepth))
    newtime.index=range(len(newtime))
    Index=[]   # all dives for each profile
    for k in range(len(newdepth)):
        if newdepth[k]<2:
            I=find_start(k,newdepth)
            Index.append(I)
    Index = [list(x) for x in set(tuple(x) for x in Index)]
    Index.sort()
    INdex=[]  # all upcast index for each profile
    top_time=[]  #the time of end of one upcast
    for k in range(len(Index)):
        max_depth=max(newdepth[Index[k][0]:Index[k][1]+1])
        bottom=[]
        for j in range(Index[k][0],Index[k][1]+1):
            if newdepth[j]==max_depth:
                if newdepth[j]>=waterdepth*0.9:
                    bottom.append(j)
        if bottom==[]:
            pass
        else:
            INdex.append([bottom[-1],Index[k][1]])
            top_time.append(newtime[Index[k][1]])
    
    N=closest_time(obsTime[i],top_time)  #find the nearst index of upcast with profile.
    for k in range(len(dep)):
        if k==INdex[N][0]:
            down=Index_all[k]   #the clostest bottom index 
        if k==INdex[N][1]:
            up=Index_all[k]    #the clostest top index
    print 'i'
    plt.figure()     
    plt.plot(temp[down:up],depth[down:up],'ro-',linewidth=1)
    plt.plot(obsTemp[i],obsDepth[i],'bo-', label='telemetered',linewidth=3)
    plt.xlim([0, 30])
    plt.ylim([max(obsDepth[i])+3, -1])
    plt.xlabel('Temp', fontsize=10)
    plt.ylabel('Depth', fontsize=10)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    #ax.legend(loc='lower right')
    plt.title('profile',fontsize=25)
    plt.text(1,0,'time:'+str(obsTime[i])+'')
    plt.text(1,1,'location:'+str(round(obsLon[i],2))+', '+str(round(obsLat[i],2))+'')
    plt.show()
    