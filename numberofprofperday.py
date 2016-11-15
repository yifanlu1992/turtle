# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 13:45:38 2016
calculate number of profile per day for turtles
@author: yifan
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from turtleModule import np_datetime
obsData = pd.read_csv('ctdWithModTempByDepth.csv')
tf_index = np.where(obsData['TF'].notnull())[0]
obsturtle_id=pd.Series(obsData['REF'][tf_index],index=tf_index)
obsturtle_ids=obsturtle_id.unique()
length=len(obsturtle_ids)    #length is number of turtles 
ids=[]    #collect all indexes of each turtle
for i in range(length):
    ids.append([])
    for j in tf_index:
        if obsturtle_id[j]==obsturtle_ids[i]:
            ids[i].append(j)    
time_aves=[]    #collect each turtle`s average profile
for i in range(len(ids)):
    obs=obsData.ix[ids[i]]   
    d=pd.Series(np_datetime(obs['END_DATE']), index=ids[i]) 
    d=d.order()    
    d.index=range(len(d))
    diff_time=[]   #list the number of days between every two profile of each turtle
    for j in d.index:
        if j+1 == len(d):
            break
        t=(d[j+1]-d[j]).total_seconds()/3600/24   
        diff_time.append(t)
    time_ave=round((len(d)/sum(diff_time)),1)   #profile per day for each turtle (average)
    time_aves.append(time_ave)
ave_time=round(np.mean(np.array(time_aves)),2)
std_time=round(np.std(np.array(time_aves)),2)
time_aves.sort()
time_aves=pd.Series(time_aves)
y1=time_aves.value_counts()   #calculate amount of turtle of each average
y1=y1.sort_index()
x1=time_aves.unique()    
fig=plt.figure()
width=0.15
plt.bar(x1,y1,align="center",width=width,color='green')
plt.xlim([0,5])
plt.ylim([0,22]) 
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.xlabel('profile per day',fontsize=10)
plt.ylabel('Number of turtle',fontsize=10)
plt.title(str(ave_time)+' day average #profiles/day with standard deviation of '+str(std_time),fontsize=12)
plt.savefig('profile per day.png')
plt.show()
