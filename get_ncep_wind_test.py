# program to test get_ncep_wind
import numpy as np
from datetime import datetime as dt
from datetime import timedelta as td
import netCDF4
import pandas as pd
def get_weather_wind(year=2016,latp=42.5,lonp=-70.5,datet=dt(2016,5,2)):
    option='one time' #'all year' or 'one time'
    #url_v='http://esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface/vwnd.sig995.'+str(year)+'.nc?vwnd,lat,lon,time'
    #url_u='http://esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface/uwnd.sig995.'+str(year)+'.nc?uwnd'
    #url_v='/net/data5/jmanning/wind/ncep/2016/vwnd.sig995.'+str(year)+'.nc'#?vwnd,lat,lon,time' # where I downloaded these from ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface/
    #url_u='/net/data5/jmanning/wind/ncep/2016/uwnd.sig995.'+str(year)+'.nc'#?uwnd'
    url_v='/home/zdong/yifan/wind/vwnd.sig995.'+str(year)+'.nc' # where I downloaded these from ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface/
    url_u='/home/zdong/yifan/wind/uwnd.sig995.'+str(year)+'.nc'
    ncv=netCDF4.Dataset(url_v)
    ncu=netCDF4.Dataset(url_u)
    lat=ncv['lat'][:]
    lon=ncv['lon'][:]
    # find the node of interest
    lat_i = list(np.arange(len(lat[0:])))
    lat_i.reverse()
    lat = list(lat[0:])
    lat.reverse()
    lon = list(lon[0:] - 360)
    lon_i = list(np.arange(len(lon[0:])))
    idlat = int(round(np.interp(latp,lat,lat_i)))
    idlon = int(round(np.interp(lonp,lon,lon_i)))
    # find the time index where there are estimates every 6 hours
    t=ncv['time'][:]
    moddate=[]
    for k in range(len(t)):
      moddate.append(dt(year,1,1,0,0,0)+td(hours=t[k]-t[0]))
      # datearray = np.array(pd.date_range(dt(year,1,1,0,0,0),freq='6H',periods=len(t)).tolist())
      # vitalii's way as follows:
      # datearray = np.arange(dt(year,1,1,0,0,0),dt(year,1,1,0,0,0)+td(hours=len(t)*4),td(hours=6)).astype(dt)
    if option=='one time':
      iddate=np.argmin(abs(np.array(moddate)-datet))
      uw=ncu['uwnd']
      u=uw[iddate,idlat,idlon]
      vw=ncv['vwnd']
      v=vw[iddate,idlat,idlon]
    elif option=='all year':
      uw=ncu['uwnd']
      u=uw[:,idlat,idlon]
      vw=ncv['vwnd']
      v=vw[:,idlat,idlon]
    #print u,v,t[iddate]
    return [u,v,t[iddate]]
