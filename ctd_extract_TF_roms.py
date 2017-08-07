'''
Extract data file ctd_extract_good.csv, add new column "TF".
If TF==True, data is good.
If TF==False, data is bad.
'''
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from turtleModule import mon_alpha2num, np_datetime, dist
r = 3                           # the ctd position that has gps position within (r) kilometers might be considered as good data.
hour = 3                        # the ctd time that has gps time within (hour) hours might be considered as good data.
ctd = pd.read_csv('2014_04_16_rawctd.csv') # original data file
ctdlat = ctd['LAT']
ctdlon = ctd['LON']
ctdtime = np_datetime(ctd['END_DATE'])
gps = pd.read_csv('2014_04_16_rawgps.csv') # orginal data file
gpslat = gps['LAT']
gpslon = gps['LON']
gpstime = np_datetime(gps['D_DATE'])
lonsize = [np.min(ctdlon), np.max(ctdlon)]
latsize = [np.min(ctdlat), np.max(ctdlat)]

index = []
i = 0
for lat, lon, ctdtm in zip(ctdlat, ctdlon, ctdtime):
    l = dist(lon, lat, gpslon, gpslat)
    p = np.where(l<r)
    maxtime = ctdtm+timedelta(hours=hour)
    mintime = ctdtm-timedelta(hours=hour)
    mx = gpstime[p[0]]<maxtime
    mn = gpstime[p[0]]>mintime
    TF = mx*mn
    if TF.any():
        index.append(i)
    i += 1
    #print(i)
ctd_TF = pd.Series([True]*len(index), index=index)
ctd['TF'] = ctd_TF
'''print(ctd)
print('{0} is OK(including "null" lon and lat values.).'.format(len(ctd_TF)/28975.0))
print('{0} is OK.'.format(len(ctd_TF)/15657.0))
print("save as 'ctd_extract_good.csv'")
ctd.to_csv('ctd_extract_good.csv')
'''