# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 12:34:36 2016

@author: zdong
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 13:45:38 2016
calculate number of days for each turtle which data is all raw data
@author: yifan
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from turtleModule import np_datetime
obsData = pd.read_csv('ctdWithModTempByDepth.csv')
obsturtle_id=pd.Series(obsData['REF'])
obsturtle_ids=obsturtle_id.unique()   

ids=[]    #collect all indexes of each turtle
for i in range(len(obsturtle_ids)):
    ids.append([])
    for j in range(len(obsturtle_id)):
        if obsturtle_id[j]==obsturtle_ids[i]:
            ids[i].append(j)    

days=[]    #collect  number of days for each turtle
for i in range(len(ids)):
    obs=obsData.ix[ids[i]]   
    d=pd.Series(np_datetime(obs['END_DATE']), index=ids[i]) 
    d=d.order()    
    d.index=range(len(d))
    n_day=(d[len(d)-1]-d[0]).total_seconds()/3600/24   #calculate the number of days for  each turtle
    n_day=int(round(n_day))
    if i==82:  #  this turtle[tu68-H383-11] have the data from 2011-06-04 20:00:00 to 2013-04-17 12:00:00 
        print 'max number of days is : '+ str(n_day)
        print obsturtle_ids[i]
    days.append(n_day)
    
ave_day=round(np.mean(np.array(days)),2)
std_day=round(np.std(np.array(days)),2)
days.sort()
days=pd.Series(days)
y1=days.value_counts()   #calculate number of turtle 
y1=y1.sort_index()
x1=days.unique()
ynew=[]
xnew=[]
for a in np.arange(35): # make the x axis 20 equal parts
    sum0=0
    for i in np.arange(len(x1)):
        if x1[i]>=a*20 and x1[i]<a*20+20:
            sum0=sum0+y1[x1[i]]
    ynew.append(sum0)
    xnew.append(a*20)
    
fig=plt.figure()
width=10
plt.bar(xnew,ynew,align="center",width=width,color='green')
plt.xlim([0,700])
plt.ylim([0,12]) 
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.xlabel('number of days',fontsize=10)
plt.ylabel('Number of turtle',fontsize=10)
plt.title(str(ave_day)+'  average #day of profile with standard deviation of '+str(std_day),fontsize=12)
plt.savefig('days of profile_try(all data).png')
plt.show()
