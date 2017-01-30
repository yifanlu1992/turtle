# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 12:27:49 2017
Calculates the "Isothermal Layer Depth" derived from turtle profiles
@author: yifan
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from turtleModule import str2ndlist,np_datetime

#############################################################################
#  HARDCODES
obsData = pd.read_csv('ctdWithModTempByDepth.csv') # has both observed and modeled profiles
ptt_of_interst=118905  # how do we know this is associated with 12487?
###############################################################################

tf_index = np.where(obsData['TF'].notnull())[0]    # get the index of good data
obsturtle_id=pd.Series(obsData['PTT'][tf_index],index=tf_index) # gets all the turtle ids
indx=np.where(obsData['PTT'][tf_index]==ptt_of_interst) # finds the index of the ppt_of_interest
        
obsTime = pd.Series(np_datetime(obsData['END_DATE'][indx]), index=indx)
obsTemp = pd.Series(str2ndlist(obsData['TEMP_VALS'][indx]), index=indx)
obsDepth = pd.Series(str2ndlist(obsData['TEMP_DBAR'][indx]), index=indx)

obsTime.sort()                    # sorts by time
Temp=obsTemp.ix[obsTime.index]    # makes a sorted Temp
Depth=obsDepth.ix[obsTime.index]  # makes a sorted Depth


plt.figure()
ILD=[] # list of isothermal layer depths for this turtle
for i in obsTime.index:
    #if i!= 3 and i!= 44 and i!=50: # these profile have something wrong 
    try:   
        if Temp[i][0]==Temp[i][1]:
            min_slope=1000  # 1000 is a large of random that represents infinity
        else:
            min_slope=(Depth[i][1]-Depth[i][0])/(Temp[i][0]-Temp[i][1])
        for k in range(1,9):   # each profile have the 10 points record
            if Temp[i][k]==Temp[i][k+1] or Depth[i][k+1]==Depth[i][k]:
                break # this breaks out of the for loop
            sl=(Depth[i][k+1]-Depth[i][k])/(Temp[i][k]-Temp[i][k+1])
            #print sl
            if sl<=min_slope:
                min_slope=sl
            else:
                break # this breaks out of the for loop
            ild=(Depth[i][k+1]+Depth[i][k])/2 # ILD means isothermal layer depth
       
    except IndexError :
        continue
    ILD.append(ild)
plt.plot(obsTime,ILD,'b',linewidth=2)
#plt.ylim([20, -1])
plt.xlabel('time', fontsize=10)
plt.ylabel('isothermal layer depth', fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
#plt.legend(loc='lower right')
plt.title('ILD vs time for '+str(ptt_of_interest),fontsize=14)
#plt.text(1,0,'time:'+str(obsTime[i])+'')
#plt.text(1,1,'location:'+str(round(obsLon[i],2))+', '+str(round(obsLat[i],2))+'')
plt.savefig('ILDvstime'+str(ptt_of_interest)+'.png')
plt.show()
