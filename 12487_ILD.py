# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 11:57:08 2017

@author: yifan
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from turtleModule import str2ndlist,np_datetime

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

for m in Index:
    if m%10==0: # every 10 profiles plotted in one picture
        plt.figure() 
        for i in range(m,m+10):
            if i==79: # total profile is 78
                break
            if i!= 3 and i!= 44 and i!=50: # these profile have something wrong 
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
                        ILD=[(obsTemp[i][k]+obsTemp[i][k+1])/2,(obsDepth[i][k+1]+obsDepth[i][k])/2] # ILD means isothermal layer depth
                    X=[]
                    for j in range(10):  # seperated the profile by 2 degree difference
                        x=obsTemp[i][j]+5*i
                        X.append(x)
                    plt.plot(X,obsDepth[i],'b',linewidth=2)
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
                    #plt.text(1,1,'location:'+str(round(obsLon[i],2))+', '+str(round(obsLat[i],2))+'')
                    plt.savefig('ILD_%s.png'%(m))
                except IndexError :
                    continue
        plt.show()
        