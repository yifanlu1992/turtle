# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 15:18:57 2017
draw the good positions of turtle and 70m contour line 
@author: yifan
"""
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import pandas as pd
from turtleModule import draw_basemap

ctd = pd.read_csv('ctd_extract_good.csv', index_col=0)
TF = ctd['TF']              # If True, data is good, if False, data is bad.
latGoodCTD, lonGoodCTD = ctd['LAT'][TF==True], ctd['LON'][TF==True]
obsData = pd.read_csv('ctdWithdepthofbottom_roms.csv')
obsLon, obsLat = obsData['LON'], obsData['LAT']     #use for plotting depth line
depthBottom = pd.Series(obsData['depth_bottom'],index=obsData.index)
for i in obsData.index:
    if depthBottom[i]>200:
        depthBottom[i]=200

lonsize = [-79.5, -71.5]
latsize = [33, 41]
#lonsize = [np.amin(lonGoodCTD), np.amax(lonGoodCTD)]
#latsize = [np.amin(latGoodCTD), np.amax(latGoodCTD)]
fig =plt.figure()
ax = fig.add_subplot(111)
plt.scatter(lonGoodCTD, latGoodCTD, color='y',s=1, label='Good Profiles')
draw_basemap(fig, ax, lonsize, latsize, interval_lon=2, interval_lat=2)

lon_is = np.linspace(lonsize[0],lonsize[1],150)
lat_is = np.linspace(latsize[0],latsize[1],150)  #use for depth line
depth_i=griddata(np.array(obsLon),np.array(obsLat),np.array(depthBottom),lon_is,lat_is,interp='linear')
cs=plt.contour(lon_is, lat_is,depth_i,levels=[70],colors = 'r',linewidths=2,linestyles='--')  #plot 100m depth
ax.annotate('70m water depth',xy=(-75.7089,34.5195),xytext=(-75.0034,34.0842),arrowprops=dict(facecolor='black'))
#plt.clabel(cs,'%.0f'%70.000,fmt='%s %m',inline=True,colors='k',fontsize=10)#fmt='%2.1d'
plt.title('Positions of turtle profiles after quality control checks', fontsize=10)
plt.legend(loc='lower right',fontsize = 'x-small')
plt.savefig('turtle_map_new',dpi=200)
plt.show()
