# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 13:45:38 2016
calculate number of profile per day for turtles which data is extracted from 'ctdWithModTempByDepth.csv'
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
  
ids=[]    #collect all indexes of each turtle
for i in range(len(obsturtle_ids)):
    ids.append([])
    for j in tf_index:
        if obsturtle_id[j]==obsturtle_ids[i]:
            ids[i].append(j)    

days=[]    #collect each turtle`s number of days
for i in range(len(ids)):
    obs=obsData.ix[ids[i]]   
    d=pd.Series(np_datetime(obs['END_DATE']), index=ids[i]) 
    d=d.order()    
    d.index=range(len(d))
    n_day=0   #calculate the number of days for  each turtle
    for j in d.index:
        if j+1 == len(d):
            break
        t=(d[j+1]-d[j]).total_seconds()/3600/24   
        n_day+=t
    n_day=int(round(n_day))
    if i==10:  #  this turtle have the data from 2009-08-24 23:17:20 to 2010-07-09 19:17:20
        print 'max number of days is : '+ str(n_day)
    days.append(n_day)
    
ave_days=round(np.mean(np.array(days)),2)
std_days=round(np.std(np.array(days)),2)
days.sort()
days=pd.Series(days)
y1=days.value_counts()   #calculate amount of turtle of each average
y1=y1.sort_index()
x1=days.unique()    
fig=plt.figure()
width=0.45
plt.bar(x1,y1,align="center",width=width,color='green')
plt.xlim([0,320])
plt.ylim([0,15]) 
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.xlabel('number of days',fontsize=10)
plt.ylabel('Number of turtle',fontsize=10)
plt.title(str(ave_days)+'  average #day of profile with standard deviation of '+str(std_days),fontsize=12)
plt.savefig('days of profile(extracted data).png')
plt.show()
