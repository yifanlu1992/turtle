# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 12:08:52 2016

@author: yifan
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from turtleModule import str2ndlist
##################################main code####################################
secondData=pd.read_csv('turtle_13137_tdr.csv')
data=pd.read_csv('13137up_down.csv')
date=pd.Series(secondData['Date'])
Date=[]
for i in date:
    Date.append(i[0:10])#remove minute and second of date
down=pd.Series(str2ndlist(data['down'],bracket=True))
dive_date=[]
for i in range(len(down)):
    idx=down[i][0]
    dive_date.append(Date[int(idx)])#get the date of every dive
dive_date=pd.Series(dive_date)
dive_dates=dive_date.unique()#'list' object has no attribute 'unique',so need use pandas before
day_dives=[0]*len(dive_dates)
for i in range(len(dive_dates)):
    for j in range(len(dive_date)):
        if dive_date[j]==dive_dates[i]:
            day_dives[i]+=1#get all dives of each day
ave_dive=round(np.mean(np.array(day_dives)),2)
std_dive=round(np.std(np.array(day_dives)),2)    
day_dives.sort()
day_dives=pd.Series(day_dives)
y=day_dives.value_counts()
y=y.sort_index()
x=day_dives.unique()
fig=plt.figure()
width=0.7
plt.bar(x,y,align="center",width=width,color='green' ,label=( 'standard deviation:'+str(std_dive) +'\n'+'average dives:'+str(ave_dive)))
plt.legend(loc='best')
plt.xlim([0,35])
plt.ylim([0,3]) 
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlabel('dives one day',fontsize=18)
plt.ylabel('amount of dives',fontsize=18)
plt.title('Turtle 13137 dives',fontsize=20)
plt.savefig('profile per day.png')  # ,' average dives:'+str(ave_dive)
plt.show()


