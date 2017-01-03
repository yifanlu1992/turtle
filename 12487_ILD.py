# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 11:57:08 2017

@author: yifan
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from turtleModule import str2ndlist,np_datetime
from gettopo import gettopo

###########################################################################
obsData = pd.read_csv('ctdWithModTempByDepth.csv') 
tf_index = np.where(obsData['TF'].notnull())[0]    # get the index of good data
obsturtle_id=pd.Series(obsData['PTT'][tf_index],index=tf_index)

indx=[]  
for i in tf_index:
    if obsturtle_id[i]==118905:   #this turtle is same turtle with 4-second turtle
        indx.append(i)

obsTime = pd.Series(np_datetime(obsData['END_DATE'][indx]), index=indx)
obsTemp = pd.Series(str2ndlist(obsData['TEMP_VALS'][indx]), index=indx)
obsDepth = pd.Series(str2ndlist(obsData['TEMP_DBAR'][indx]), index=indx)
Index=range(len(indx))
obsTime.index=Index
obsTemp .index=Index
obsDepth.index=Index

plt.figure() 
for i in range(3):
    for j in range(3):
        A=10
        if i==j: 
           plt.plot(obsTemp[i],obsDepth[i],'b',linewidth=2)
        plt.xlim([0, 30])
        plt.ylim([10, -1])
        plt.xlabel('Temp', fontsize=10)
        plt.ylabel('Depth', fontsize=10)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
    #plt.legend(loc='lower right')
    #plt.title(i,fontsize=14)
    #plt.text(1,0,'time:'+str(obsTime[i])+'')
    #plt.text(1,1,'location:'+str(round(obsLon[i],2))+', '+str(round(obsLat[i],2))+'')
    #plt.savefig('%s.png'%(i))
plt.show()
