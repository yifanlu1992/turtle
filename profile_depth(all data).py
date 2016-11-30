# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 10:27:13 2016
calculate the profile depth 
@author: yifan
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd        
obsData = pd.read_csv('ctd_good.csv') 
depth = obsData['MAX_DBAR'][:] # get deepest data file depth
obsID = obsData['PTT'][:]           # Get ID of turtle.

index1 = depth[depth>200].index
id = obsID[index1]
print 'ID of observations that depth>200:', id.drop_duplicates().values#find turle depth more than 200m 

depth.sort()
y=depth.unique()
x=depth.value_counts()
x=x.sort_index()

ynew,xnew=[],[]
for i in range(29):
    sum=0
    for j in range(len(y)):
        if y[j]>i*5 and y[j]<i*5+5:
            sum+=x[y[j]]
    xnew.append(sum)
    ynew.append(i*5)

fig = plt.figure()            
plt.barh(ynew, xnew,height=3,color='green')
plt.ylim(145,0)
plt.ylabel('profile Depth', fontsize=12)
plt.xlabel('Quantity', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(np.arange(0,145,5),fontsize=8)
plt.title('Profile  Depth ', fontsize=16)
plt.savefig('profile_depth(all data).png', dpi=100)
plt.show()