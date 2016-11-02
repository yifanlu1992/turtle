# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 13:45:38 2016

@author: yifan
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from turtleModule import np_datetime
##################################main code####################################
obsData = pd.read_csv('ctdWithModTempByDepth.csv')
tf_index = np.where(obsData['TF'].notnull())[0]
obsturtle_id=pd.Series(obsData['REF'][tf_index],index=tf_index)
obsturtle_ids=obsturtle_id.unique()
length=len(obsturtle_ids)                #length is number of turtles 
ids=[]
for i in range(length):
    ids.append([])   
for i in range(length):
    for j in tf_index:
        if obsturtle_id[j]==obsturtle_ids[i]:
            ids[i].append(j)   #collect index of each turtle 
time_aves=[] #collect all average

for i in range(len(ids)):
    obs=obsData.ix[ids[i]]   
    goodTime=pd.Series(np_datetime(obs['END_DATE']), index=ids[i])
    d=pd.DataFrame({'Time':goodTime},index=ids[i])   #each turtle`s modtemp and obstemp
    d=d.sort(['Time'])
    d.index=range(len(d))
    diff_time=[]
    for j in d.index:
        if j+1 == len(d):
            break
        t=(d['Time'][j+1]-d['Time'][j]).total_seconds()/3600/24
        diff_time.append(t)
    time_ave=round(len(d)/sum(diff_time))#each turtle`s average
    time_aves.append(time_ave)
print time_aves
ave_time=round(np.mean(np.array(time_aves)),2)
std_time=round(np.std(np.array(time_aves)),2)
time_aves.sort()
time_aves=pd.Series(time_aves)
y1=time_aves.value_counts()
y1=y1.sort_index()
x1=time_aves.unique()
  #calculate quantity of each average number and use for plotting

fig=plt.figure()
plt.bar(x1,y1)
plt.title('Average profiles:'+str(ave_time)+' standard deviation:'+str(std_time))
plt.xlabel('profile per day')
plt.ylabel('Turtle numbers')
plt.ylim([0,100])
plt.xlim([0,5])
plt.savefig('profile per day.png')
plt.show()

