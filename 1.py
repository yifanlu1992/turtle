# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 23:55:09 2016

@author: JIM
"""


"""
Created on Tue Nov 29 23:46:38 2016
return ratio of the deepest depth
@author: yifan
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
        # A module of classes that using ROMS, FVCOM.
obsData = pd.read_csv('ctd_extract_good.csv') # From nearestIndexInMod.py
tf_index = np.where(obsData['TF'].notnull())[0] # Get  index of good data.
obsLat, obsLon = obsData['LAT'][tf_index], obsData['LON'][tf_index]
obsDeepest = obsData['MAX_DBAR'][tf_index] # get deepest data file depth
obsID = obsData['PTT'][tf_index]           # Get ID of turtle.




fig = plt.figure()
#ax = fig.add_subplot(111)
p = obsDeepest
'''index1 = p[p>1.5].index
id = obsID[index1]
print 'ID of observations that obsDeepest/newH>1.5:', id.drop_duplicates().values
'''

p.sort()
y=p.unique()
x=p.value_counts()
x=x.sort_index()

ynew,xnew=[],[]
for i in range(22):
    sum=0
    for j in range(len(y)):
        if y[j]>i*4 and y[j]<i*4+4:
            sum+=x[y[j]]
    xnew.append(sum)
    ynew.append(i*4)

            
plt.barh(ynew, xnew,height=3)
plt.ylim(88,0)
plt.ylabel('profile Depth', fontsize=12)
plt.xlabel('Quantity', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(np.arange(0,88,4),fontsize=8)
plt.title('Profile  Depth ', fontsize=20)
plt.savefig('deepestDepth.png', dpi=200)
plt.show()
