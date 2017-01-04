# -*- coding: utf-8 -*-
"""
Created on Tue Jan 03 20:27:20 2017

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
for i in range(10):
    for k in range(1,9):
        min_slope=(obsDepth[i][1]-obsDepth[i][0])/(obsTemp[i][0]-obsTemp[i][1])
        if obsTemp[i][k]==obsTemp[i][k+1]:
            #ILD=[(obsTemp[i][k]+obsTemp[i][k+1])/2,(obsDepth[i][k+1]+obsDepth[i][k])/2]
            break
        sl=(obsDepth[i][k+1]-obsDepth[i][k])/(obsTemp[i][k]-obsTemp[i][k+1])
        #print sl
        if sl<=min_slope:
            min_slope=sl
        else:
            break
        ILD=[(obsTemp[i][k]+obsTemp[i][k+1])/2,(obsDepth[i][k+1]+obsDepth[i][k])/2]
    X=[]
    for j in range(10):  # seperated the profile by 2 degree difference
        x=obsTemp[i][j]+2*i
        X.append(x)
    plt.plot(X,obsDepth[i],'b',linewidth=2)
    plt.plot((ILD[0]-0.5+2*i,ILD[0]+0.5+2*i),(ILD[1],ILD[1]),linewidth=2,color='red')
    plt.xlim([18, 70])
    plt.ylim([45, -1])
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
