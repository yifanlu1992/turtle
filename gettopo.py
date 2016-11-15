# routine to get USGS global 30 sec arc topo data given lat & lon
import numpy as np
import matplotlib.pyplot as plt
import urllib
import netCDF4
#from mpl_toolkits.basemap import Basemap
#import sys
#sys.path.append('../modules')
#from conversions import m2fth

#lat=40.4
#lon=-70.03
isub=1
       # this function get nearby points and finds the closest value

def gettopo(lat,lon):
    base_url='http://coastwatch.pfeg.noaa.gov/erddap/griddap/usgsCeSrtm30v6.nc?'
    query='topo[(%f):%d:(%f)][(%f):%d:(%f)]' % (lat+.1,isub,lat-.1,lon-.1,1,lon+.1)
    url = base_url+query
    file='usgsCeSrtm30v6.nc'
    urllib.urlretrieve (url, file)
    nc = netCDF4.Dataset(file)
    ncv = nc.variables
    lons = ncv['longitude'][:]
    lats = ncv['latitude'][:]
    indlon,val=min(enumerate(lons),key=lambda x: abs(x[1]-lon))
    indlat,val=min(enumerate(lats),key=lambda x: abs(x[1]-lat))
    return ncv['topo'][indlat,indlon]

#
#depth=gettopo(lat,lon)
#print str(depth)+' meters or '+str(m2fth(depth))+' fathoms'
