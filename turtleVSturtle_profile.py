# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 11:56:41 2017
compare the turtle's profile to each other within 1 km and 1 day  
@author: yifan
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import  timedelta
from turtleModule import np_datetime, dist,str2ndlist
from gettopo import gettopo

def small_dist_index(lon,lat,all_lon,all_lat):  # obtain the index according to the exacted positon within 1 km with all other turtle position
    for i in range(len(all_lon)):
        if all_lon[i]==lon:
            start=i+1
    index=[] # find the index within 1 kilometers 
    for i in range(start,len(all_lon)):
        l=dist(lon,lat,all_lon[i],all_lat[i])
        if l<=1:
            index.append(i)
    return index

Data=pd.read_csv('ctdWithModTempByDepth.csv',index_col=0)
tf_index=np.where(Data['TF'].notnull())[0]
turtle_id=Data['PTT'][tf_index]
lat = Data['LAT'][tf_index]
lon = Data['LON'][tf_index]
time = pd.Series(np_datetime(Data['END_DATE'][tf_index]),index=tf_index)
depth=pd.Series(str2ndlist(Data['TEMP_DBAR'][tf_index]),index=tf_index)
temp=pd.Series(str2ndlist(Data['TEMP_VALS'][tf_index]),index=tf_index)
indx=range(len(tf_index))
turtle_id.index=indx
lat.index=indx
lon.index=indx
time.index=indx
depth.index=indx
temp.index=indx

g_index=[] # each profile gives indexs including the matched profiles ,but most of them are empty.
for i in indx:
    if i%1000==0:
        print i
    l_index=small_dist_index(lon[i],lat[i],lon,lat)
    ind=[]
    for j in l_index:
        if turtle_id[j]!=turtle_id[i]:# we want to get the profiles from different turtle
            t=abs(time[i]-time[j])
            if t<=timedelta(days=1): # the time interval is one day
                ind.append(j)
    g_index.append(ind)

g_index=pd.Series(g_index,index=indx)
INDX1,INDX2=[],[] # get the raw index that profile have the matched profiles
for i in indx:
    if g_index[i]!=[]:
        INDX1.append(i)  #get rid of the empty ,there are 92 profiles in the last
        if len(g_index[i])!=1:
            INDX2.append(i)  # obtain the number of the turtle profile more than one
AV=[]
mean_STD=[]

n=0                              
for i in INDX1:#[1727,4011,4801,8038]:#[3274,4557,9106]:##[4801]:#[656,5327,1728,8699]:###
    ix=[] #g_index[i]
    
    cd=max (depth[i])
    
    for j in np.arange(len(g_index[i])): 
        d=max (depth[g_index[i][j]])
        if i==4801:
            ix.append(g_index[i][j])
        if abs(cd-d)<5:
            wd1=-gettopo(lat[i],lon[i])
            wd2=-gettopo(lat[g_index[i][j]],lon[g_index[i][j]])
            if cd>=wd1*0.9 and d>=wd2*0.9:
                ix.append(g_index[i][j])
    
    if ix!=[]:
        print i
        n+=1
        A=[]
        S=[]

        '''b=np.mean(temp[i])
        S.append(b)
        for k in range(len(ix)):
            c=np.mean(temp[ix[k]])
            a=b-c
            A.append(a)
            S.append(c)
        av=np.mean(np.array(A))
        std=np.std(np.array(S))

        print av ,std  
        AV.append(abs(av))
        mean_STD.append(std)'''
        b=temp[i][-1]
        S.append(b)
        for k in range(len(ix)):
            c=temp[ix[k]][-1]
            a=b-c
            A.append(a)
            S.append(c)
        av=np.mean(np.array(A))
        std=np.std(np.array(S))

        print av ,std  
        AV.append(abs(av))
        mean_STD.append(std)
        
        
        '''fig=plt.figure()
        plt.plot(temp[i],depth[i],'r',label='id: '+ str(turtle_id[i]))
        plt.plot(temp[ix[0]],depth[ix[0]],'b',label='id: '+ str(turtle_id[ix[0]]))
        for k in range(1,len(ix)):
           if turtle_id[ix[0]]==turtle_id[ix[k]]: 
               plt.plot(temp[ix[k]],depth[ix[k]],'b')
           elif k<=2:
               plt.plot(temp[ix[k]],depth[ix[k]],'g',label='id: '+str(turtle_id[ix[k]]))
           else:
               plt.plot(temp[ix[k]],depth[ix[k]],'g')

        plt.text(20,41,'mean: '+str(round(av,2))+u'°C',fontsize=10)
        plt.text(20,44,'std: '+str(round(std,2))+u'°C',fontsize=10)
        plt.ylim([45,0])
        plt.xlim([5,25])
        plt.ylabel('Depth(m)',fontsize=14)
        plt.xlabel('Temperature('+u'°C'+')',fontsize=12)
        plt.title('compare with 1 day and 1 km',fontsize=14)
        plt.legend(loc='upper left')
        #plt.savefig(str(turtle_id[i])+'new',dpi=200)
        plt.show()'''
print 'mean difference temp: ',np.mean(np.array(AV))
print 'mean std difference temp: ',np.mean(np.array(mean_STD))
print 'std mean difference temp: ',np.mean(np.std(AV))
print n
 