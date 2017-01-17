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

fig=plt.figure()
for i in range(60):
    ax1=fig.add_subplot(2,1,1)
    for j in range(0,30):
        ax1.plot(np.array(obsTemp[j])+5*j,obsDepth[j],'b',linewidth=1)
        ax1.plot((ILD[j][0]-1.5+5*j,ILD[j][0]+1.5+5*j),(ILD[j][1],ILD[j][1]),linewidth=2,color='red')# remark the ILD by a short transverse line
        ax1.set_ylim([40,-1])
        plt.setp(ax1.get_xticklabels() ,visible=False)
    ax2=fig.add_subplot(2,1,2)
    for j in range(30,60):
        ax2.plot(np.array(obsTemp[j])+5*j,obsDepth[j],'b',linewidth=1)
        ax2.plot((ILD[j][0]-1.5+5*j,ILD[j][0]+1.5+5*j),(ILD[j][1],ILD[j][1]),linewidth=2,color='red')# remark the ILD by a short transverse line
        ax2.set_ylim([40,-1])
        plt.setp(ax2.get_xticklabels() ,visible=False)

fig.text(0.5, 0.04, 'Temperature (20 degree of interval)', ha='center', va='center', fontsize=14)#  0.5 ,0.04 represent the  plotting scale of x_axis and y_axis
fig.text(0.06, 0.5, 'Depth(m)', ha='center', va='center', rotation='vertical',fontsize=14)
plt.savefig('obs_ILD.png',dpi=200)
plt.show()
