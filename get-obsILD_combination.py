# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 11:57:08 2017
calculate the "Isothermal Layer Depth" derived from turtle profiles
@author: yifan
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from turtleModule import str2ndlist,np_datetime

obsData = pd.read_csv('ctdWithModTempByDepth.csv') # has both observed and modeled profiles
tf_index = np.where(obsData['TF'].notnull())[0]    # get the index of good data
obsturtle_id=pd.Series(obsData['PTT'][tf_index],index=tf_index)
#indx=np.where(obsData['PTT'][tf_index]==ptt_of_interst)[0]
indx=[]  # this indx is to get the specifical turtle all index in obsData ,if we use the "where" function ,we just get the length  of tf_index.
for i in tf_index:
    if obsturtle_id[i]==118905:   #this turtle is same turtle with 4-second turtle
        indx.append(i)
       
obsTime = pd.Series(np_datetime(obsData['END_DATE'][indx]), index=indx)
obsTemp = pd.Series(str2ndlist(obsData['TEMP_VALS'][indx]), index=indx)
obsDepth = pd.Series(str2ndlist(obsData['TEMP_DBAR'][indx]), index=indx)
Index=range(len(indx))
obsTime.index=Index
obsTemp.index=Index
obsDepth.index=Index

ILD=[]#get the "isothermal layer depth" of all profile
for i in range(60):# just want plot six picture
    try:   
        if obsTemp[i][0]==obsTemp[i][1]:
            min_slope=1000  # 1000 is a large of random that represents infinity
        else:
            min_slope=(obsDepth[i][1]-obsDepth[i][0])/(obsTemp[i][0]-obsTemp[i][1])
        for k in range(1,9):   # each profile have the 10 points record
            if obsTemp[i][k]==obsTemp[i][k+1] or obsDepth[i][k+1]==obsDepth[i][k]:
                break
            sl=(obsDepth[i][k+1]-obsDepth[i][k])/(obsTemp[i][k]-obsTemp[i][k+1])
            #print sl
            if sl<=min_slope:
                min_slope=sl
            else:
                break
        ild=[(obsTemp[i][k]+obsTemp[i][k-1])/2,(obsDepth[i][k-1]+obsDepth[i][k])/2] # ILD means isothermal layer depth
    except IndexError :
        continue   
    ILD.append(ild)

X=[] # seperated the profile by 5 degree difference
for j in range(10):  
    x=obsTemp[i][j]+5*i
    X.append(x)

fig=plt.figure()
for i in range(6):
    ax=fig.add_subplot(2,3,i+1)
    for j in range(60):
        ax.plot(obsTemp[j]+5*j,obsDepth[j],'b',linewidth=1)
        ax.plot((ILD[j][0]-1.5+5*j,ILD[j][0]+1.5+5*j),(ILD[j][1],ILD[j][1]),linewidth=2,color='red')# remark the ILD by a short transverse line
        ax.set_ylim([35,2]) 

'''plt.plot(X,obsDepth[i],'b',linewidth=2)
plt.plot((ILD[0]-1.5+5*i,ILD[0]+1.5+5*i),(ILD[1],ILD[1]),linewidth=2,color='red')# remark the ILD by a short transverse line
#plt.xlim([15+n*55, 70+n*55])
plt.ylim([40, -1])
plt.xlabel('profile', fontsize=10)
plt.ylabel('Depth', fontsize=10)
#plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
#plt.legend(loc='lower right')
plt.title('isothermal layer depth',fontsize=14)
#plt.text(1,0,'time:'+str(obsTime[i])+'')
#plt.text(1,1,'location:'+str(round(obsLon[i],2))+', '+str(round(obsLat[i],2))+'')'''
plt.savefig('ILD.png')
    
plt.show()
