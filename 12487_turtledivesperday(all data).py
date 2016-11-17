# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 16:14:29 2016
calculate all dives one day for 4second turtle ,the turtle data we use all raw data.
@author: yifan
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from gettopo import gettopo

secondData=pd.read_csv('12487_location.csv')# /net/data5/jmanning/turtle/
index=pd.Series(secondData['index']) 
date=pd.Series(secondData['time'])
depth=pd.Series(secondData['depth'])
lon=pd.Series(secondData['lon']) 
lat=pd.Series(secondData['lat'])

index.fillna(0,inplace=True)  #replace missing value with 0
water_depth=[]  # get water depth according to lat and lon
for i in range(len(secondData)-1):
    if i%100000==0:   
        print(i)
    if i==0:
        wd1=-gettopo(lat[i],lon[i]) # this is a function to get water depth from a USGS database
    if index[i]!=index[i+1]: #  
       if index[i+1]==0:
          wd1=-gettopo(lat[i],lon[i])
       else:
           wd1=-gettopo(lat[i+1],lon[i+1])      
    water_depth.append(wd1)
wd=pd.Series(water_depth)

# the following code finds where the turtle passes the 10% depth (ie "start_dive")
start_dive=[]  #find all index of starting dive which define in 10% of water depth
for i in range(len(depth)-2):
    if depth[i]<wd[i]*0.1:  
        if depth[i]<depth[i+1]<depth[i+2] and depth[i+1]>=wd[i]*0.1:
                start_dive.append(i)
    if i%100000==0:   
        print(i)

# the following code finds where the turtle passes the 90% depth (ie "dives_index")
dives_index=[] # get the all index of dives which is selected
for i in range(len(start_dive)-1):
    max_depth=max(depth[start_dive[i]:start_dive[i+1]+1])
    for j in range(start_dive[i],start_dive[i+1]+1):   #find the max depth index in every dive
        if depth[j]==max_depth:
            if depth[j]>=wd[j]*0.9:
                dives_index.append(j)
                break

Date=[]  #remove hour, minute and second of date
for i in date:
    Date.append(i[0:10])
    
dive_date=[]  #get the date of every dive
for i in range(len(dives_index)):
    dive_date.append(Date[dives_index[i]])  # using dive's node index to find the dive date
dive_date=pd.Series(dive_date)
dive_dates=dive_date.unique()   #'list' object has no attribute 'unique',so need use pandas before

day_dives=[0]*len(dive_dates)  # in order to get all dives of each day
for i in range(len(dive_dates)):
    for j in range(len(dive_date)):
        if dive_date[j]==dive_dates[i]:
            day_dives[i]+=1

ave_dive=round(np.mean(np.array(day_dives)),2)
std_dive=round(np.std(np.array(day_dives)),2)    
day_dives.sort()
day_dives=pd.Series(day_dives)
y=day_dives.value_counts()
y=y.sort_index()
x=day_dives.unique()
fig=plt.figure()
width=0.4
plt.bar(x,y,align="center",width=width,color='green' ,label='turtle 12478 ')
plt.legend(loc='best')
plt.xlim([0,50])
plt.ylim([0,12]) 
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.xlabel('dives one day',fontsize=10)
plt.ylabel('amount of days',fontsize=10)
plt.title(str(ave_dive)+' average dives/day with standard deviation of '+str(std_dive),fontsize=12)
plt.savefig('12487_turtledives(all data).png')  # ,' average dives:'+str(ave_dive)
plt.show()

