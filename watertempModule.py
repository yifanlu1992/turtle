import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import netCDF4
import matplotlib.pyplot as plt
from matplotlib import path
import jata
class water(object):
    def __init__(self, startpoint):
        '''
        get startpoint of water, and the location of datafile.
        startpoint = [25,45]
        '''
        self.startpoint = startpoint
    def get_data(self, url):
        pass
    def bbox2ij(self, lons, lats, bbox):
        """
        Return tuple of indices of points that are completely covered by the 
        specific boundary box.
        i = bbox2ij(lon,lat,bbox)
        lons,lats = 2D arrays (list) that are the target of the subset, type: np.ndarray
        bbox = list containing the bounding box: [lon_min, lon_max, lat_min, lat_max]
    
        Example
        -------  
        >>> i0,i1,j0,j1 = bbox2ij(lat_rho,lon_rho,[-71, -63., 39., 46])
        >>> h_subset = nc.variables['h'][j0:j1,i0:i1]       
        """
        bbox = np.array(bbox)
        mypath = np.array([bbox[[0,1,1,0]],bbox[[2,2,3,3]]]).T
        p = path.Path(mypath)
        points = np.vstack((lons.flatten(),lats.flatten())).T
        tshape = np.shape(lons)
        # inside = p.contains_points(points).reshape((n,m))
        inside = []
        for i in range(len(points)):
            inside.append(p.contains_point(points[i]))
        inside = np.array(inside, dtype=bool).reshape(tshape)
        # ii,jj = np.meshgrid(xrange(m),xrange(n))
        index = np.where(inside==True)
        if not index[0].tolist():          # bbox covers no area
            raise Exception('no points in this area')
        else:
            # points_covered = [point[index[i]] for i in range(len(index))]
            # for i in range(len(index)):
                # p.append(point[index[i])
            # i0,i1,j0,j1 = min(index[1]),max(index[1]),min(index[0]),max(index[0])
            return index
    def nearest_point_index(self, lon, lat, lons, lats, length=(1, 1),num=4):
        '''
        Return the index of the nearest rho point.
        lon, lat: the coordinate of start point, float
        lats, lons: the coordinate of points to be calculated.
        length: the boundary box.
        '''
        bbox = [lon-length[0], lon+length[0], lat-length[1], lat+length[1]]
        # i0, i1, j0, j1 = self.bbox2ij(lons, lats, bbox)
        # lon_covered = lons[j0:j1+1, i0:i1+1]
        # lat_covered = lats[j0:j1+1, i0:i1+1]
        # temp = np.arange((j1+1-j0)*(i1+1-i0)).reshape((j1+1-j0, i1+1-i0))
        # cp = np.cos(lat_covered*np.pi/180.)
        # dx=(lon-lon_covered)*cp
        # dy=lat-lat_covered
        # dist=dx*dx+dy*dy
        # i=np.argmin(dist)
        # # index = np.argwhere(temp=np.argmin(dist))
        # index = np.where(temp==i)
        # min_dist=np.sqrt(dist[index])
        # return index[0]+j0, index[1]+i0
        index = self.bbox2ij(lons, lats, bbox)
        lon_covered = lons[index]
        lat_covered = lats[index]
        # if len(lat_covered) < num:
            # raise ValueError('not enough points in the bbox')
        # lon_covered = np.array([lons[i] for i in index])
        # lat_covered = np.array([lats[i] for i in index])
        cp = np.cos(lat_covered*np.pi/180.)
        dx = (lon-lon_covered)*cp
        dy = lat-lat_covered
        dist = dx*dx+dy*dy
        
        # get several nearest points
        dist_sort = np.sort(dist)[0:9]
        findex = np.where(dist==dist_sort[0])
        lists = [[]] * len(findex)
        for i in range(len(findex)):
            lists[i] = findex[i]
        if num > 1:
            for j in range(1,num):
                t = np.where(dist==dist_sort[j])
                for i in range(len(findex)):
                     lists[i] = np.append(lists[i], t[i])
        indx = [i[lists] for i in index]
        return indx, dist_sort[0:num]
        '''
        # for only one point returned
        mindist = np.argmin(dist)
        indx = [i[mindist] for i in index]
        return indx, dist[mindist]
        '''
    def nearest_point_index2(self, lon, lat, lons, lats):
        d = dist(lon, lat, lons ,lats)
        min_dist = np.min(d)
        index = np.where(d==min_dist)
        return index
class water_roms(water):
    '''
    ####(2009.10.11, 2013.05.19):version1(old) 2009-2013
    ####(2013.05.19, present): version2(new) 2013-present
    (2006.01.01 01:00, 2014.1.1 00:00)
    '''
    def __init__(self):
        pass
        # self.startpoint = lon, lat
        # self.dataloc = self.get_url(starttime)
    def get_url(self, starttime, endtime):
        '''
        get url according to starttime and endtime.
        '''
        '''
        self.starttime = starttime
        self.days = int((endtime-starttime).total_seconds()/60/60/24)+1 # get total days
        time1 = datetime(year=2009,month=10,day=11) # time of url1 that starts from
        time2 = datetime(year=2013,month=5,day=19)  # time of url2 that starts from
        url1 = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2009_da/avg?lon_rho[0:1:81][0:1:129],lat_rho[0:1:81][0:1:129],mask_rho[0:1:81][0:1:129],u[{0}:1:{1}][0:1:35][0:1:81][0:1:128],v[{0}:1:{1}][0:1:35][0:1:80][0:1:129]'
        url2 = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2013_da/avg_Best/ESPRESSO_Real-Time_v2_Averages_Best_Available_best.ncd?mask_rho[0:1:81][0:1:129],u[{0}:1:{1}][0:1:35][0:1:81][0:1:128],v[{0}:1:{1}][0:1:35][0:1:80][0:1:129],lon_rho[0:1:81][0:1:129],lat_rho[0:1:81][0:1:129]'
        if endtime >= time2:
            if starttime >=time2:
                index1 = (starttime - time2).days
                index2 = index1 + self.days
                url = url2.format(index1, index2)
            elif time1 <= starttime < time2:
                url = []
                index1 = (starttime - time1).days
                url.append(url1.format(index1, 1316))
                url.append(url2.format(0, self.days))
        elif time1 <= endtime < time2:
            index1 = (starttime-time1).days
            index2 = index1 + self.days
            url = url1.format(index1, index2)
        return url
        '''
        self.starttime = starttime
        # self.hours = int((endtime-starttime).total_seconds()/60/60) # get total hours
        # time_r = datetime(year=2006,month=1,day=9,hour=1,minute=0)
        if (starttime- datetime(2013,5,18)).total_seconds()/3600>25:
            #url_oceantime = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/hidden/2006_da/his?ocean_time'
            #url_oceantime = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2013_da/his_Best/ESPRESSO_Real-Time_v2_History_Best_Available_best.ncd?time'
            url_oceantime = 'http://tds.marine.rutgers.edu/thredds/dodsC/roms/espresso/2013_da/his/ESPRESSO_Real-Time_v2_History_Best?time[0:1:31931]'
            self.oceantime = netCDF4.Dataset(url_oceantime).variables['time'][:]    #if url2006, ocean_time.
            t1 = (starttime - datetime(2013,5,18)).total_seconds()/3600 # for url2006 it's 2006,01,01; for url2013, it's 2013,05,18, and needed to be devide with 3600
            t2 = (endtime - datetime(2013,5,18)).total_seconds()/3600
            self.index1 = self.closest_num(t1, self.oceantime)
            self.index2 = self.closest_num(t2, self.oceantime)
            # index1 = (starttime - time_r).total_seconds()/60/60
            # index2 = index1 + self.hours
            # url = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2006_da/his?h[0:1:81][0:1:129],s_rho[0:1:35],lon_rho[0:1:81][0:1:129],lat_rho[0:1:81][0:1:129],mask_rho[0:1:81][0:1:129],u[{0}:1:{1}][0:1:35][0:1:81][0:1:128],v[{0}:1:{1}][0:1:35][0:1:80][0:1:129]'
            #url = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/hidden/2006_da/his?s_rho[0:1:35],h[0:1:81][0:1:129],lon_rho[0:1:81][0:1:129],lat_rho[0:1:81][0:1:129],temp[{0}:1:{1}][0:1:35][0:1:81][0:1:129],ocean_time'
            #url = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2013_da/his_Best/ESPRESSO_Real-Time_v2_History_Best_Available_best.ncd?h[0:1:81][0:1:129],s_rho[0:1:35],lon_rho[0:1:81][0:1:129],lat_rho[0:1:81][0:1:129],temp[{0}:1:{1}][0:1:35][0:1:81][0:1:129],time' 
            url = 'http://tds.marine.rutgers.edu/thredds/dodsC/roms/espresso/2013_da/his/ESPRESSO_Real-Time_v2_History_Best?s_rho[0:1:35],lon_rho[0:1:81][0:1:129],lat_rho[0:1:81][0:1:129],time[0:1:32387],h[0:1:81][0:1:129],temp[0:1:32387][0:1:35][0:1:81][0:1:129]'
            url = url.format(self.index1, self.index2)
        else:
            #url_oceantime = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/hidden/2006_da/his?ocean_time'
            url_oceantime='http://tds.marine.rutgers.edu/thredds/dodsC/roms/espresso/2009_da/his?ocean_time'#[0:1:19145]
            self.oceantime = netCDF4.Dataset(url_oceantime).variables['ocean_time'][:]    #if url2006, ocean_time.
            t1 = (starttime - datetime(2006,1,1)).total_seconds() # for url2006 it's 2006,01,01; for url2013, it's 2013,05,18, and needed to be devide with 3600
            t2 = (endtime - datetime(2006,1,1)).total_seconds()
            self.index1 = self.closest_num(t1, self.oceantime)
            print 'self.index1' ,self.index1
            self.index2 = self.closest_num(t2, self.oceantime)
            #url = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/hidden/2006_da/his?s_rho[0:1:35],h[0:1:81][0:1:129],lon_rho[0:1:81][0:1:129],lat_rho[0:1:81][0:1:129],temp[{0}:1:{1}][0:1:35][0:1:81][0:1:129],ocean_time'
            url='http://tds.marine.rutgers.edu/thredds/dodsC/roms/espresso/2009_da/his?s_rho[0:1:35],h[0:1:81][0:1:129],lon_rho[0:1:81][0:1:129],lat_rho[0:1:81][0:1:129],ocean_time[0:1:19145],temp[0:1:19145][0:1:35][0:1:81][0:1:129]'
            url = url.format(self.index1, self.index2)
        return url
    def closest_num(self, num, numlist, i=0):
        '''
        Return index of the closest number in the list
        '''
        index1, index2 = 0, len(numlist)
        indx = int(index2/2)
        if not numlist[0] <= num < numlist[-1]:
            raise Exception('{0} is not in {1}'.format(str(num), str(numlist)))
        if index2 == 2:
            l1, l2 = num-numlist[0], numlist[-1]-num
            if l1 < l2:
                i = i
            else:
                i = i+1
        elif num == numlist[indx]:
            i = i + indx
        elif num > numlist[indx]:
            i = self.closest_num(num, numlist[indx:],
                              i=i+indx)
        elif num < numlist[indx]:
            i = self.closest_num(num, numlist[0:indx+1], i=i)
        return i
    def get_data(self, url):
        '''
        return the data needed.
        url is from water_roms.get_url(starttime, endtime)
        '''
        data = jata.get_nc_data(url, 'lon_rho', 'lat_rho', 'temp','h','s_rho', 'ocean_time')
        return data
    def watertemp(self, lon, lat, depth, time, url):
        #data = self.get_data(url)
        #url='http://tds.marine.rutgers.edu/thredds/dodsC/roms/espresso/2009_da/his'
        data =netCDF4.Dataset(url)
        lons = data['lon_rho'][:]
        lats = data['lat_rho'][:]
        #temp = data['temp']
        if type(lon) is list or type(lon) is np.ndarray:
            t = []
            for i in range(len(time)):
                #print(i)
                if i%100==0:
                    print i
                #print('depth: ', depth[i])
                watertemp = self.__watertemp(lon[i], lat[i], lons, lats, depth[i], time[i], data)
                t.append(watertemp)
                '''
                try:
                    print i, lon[i], lat[i], depth[i], time[i]
                    watertemp = self.__watertemp(lon[i], lat[i], lons, lats, depth[i], time[i], data)
                    t.append(watertemp)
                except Exception:
                    t.append(0)
                    continue
                '''
            t = np.array(t)
        else:
            #print('depth: ', depth)
            watertemp = self.__watertemp(lon, lat, lons, lats, depth, time, data)
            t = watertemp
        return t
    def __watertemp(self, lon, lat, lons, lats, depth, time, data):
        '''
        return temp
        '''
        index = self.nearest_point_index2(lon,lat,lons,lats)
        #print('index:',index)
        #print(time)
        depth_layers = data['h'][index[0][0]][index[1][0]]*data['s_rho']
        layer = np.argmin(abs(depth_layers+depth)) # Be careful, all depth_layers are negative numbers
        time_index = self.closest_num((time-datetime(2006,1,1,0,0,0)).total_seconds(),self.oceantime) - self.index1
        #print(time_index, layer, index[0][0], index[1][0])
        temp = data['temp'][time_index, layer, index[0][0], index[1][0]]
        return temp
    def layerTemp(self, layer, url):
        '''
        Get the temp of one whole specific layer.
        Only return temperature of the first 'time' index. 
        '''
        data = self.get_data(url)
        # lons, lats = data['lon_rho'][:], data['lat_rho'][:]
        # index = self.nearest_point_index2(lon, lat, lons, lats)
        # depth_layers = data['h'][index[0][0]][index[1][0]] * data['s_rho']
        # print depth_layers
        # layer = np.argmin(abs(depth_layers + depth))
        # layerTemp = data['temp'][0, layer]
        layerTemp = data['temp'][0, -layer]
        # depthRange = [depth_layers[-layer-1], depth_layers[-layer]]
        # depthRange = [depth_layers[layer-1], depth_layers[layer]]
        return layerTemp
    def depthTemp(self, depth, url):
        '''
        Return temp data of whole area in specific depth to draw contour
        '''
        data = self.get_data(url)
        temp = data['temp'][0]
        layerDepth = data['h']
        s_rho = data['s_rho']
        depthTemp = []
        for i in range(82):
            t = []
            for j in range(130):
                print(i, j, 'depthTemp')
                locDepth = layerDepth[i,j]  # The depth of this point
                lyrDepth = s_rho * locDepth
                if depth > lyrDepth[-1]: # Obs is shallower than last layer.
                    d = (temp[-2,i,j]-temp[-1,i,j])/(lyrDepth[-2]-lyrDepth[-1]) * \
                        (depth-lyrDepth[-1]) + temp[-1,i,j]
                elif depth < lyrDepth[0]: # Obs is deeper than first layer.
                    d = (temp[1,i,j]-temp[0,i,j])/(lyrDepth[1]-lyrDepth[0]) * \
                        (depth-lyrDepth[0]) + temp[0,i,j]
                else:
                    ind = self.closest_num(depth, lyrDepth)
                    d = (temp[ind,i,j]-temp[ind-1,i,j])/(lyrDepth[ind]-lyrDepth[ind-1]) * \
                        (depth-lyrDepth[ind-1]) + temp[ind-1,i,j]
                t.append(d)
            depthTemp.append(t)
        return np.array(depthTemp)

class waterCTD(water_roms):
    def watertemp(self, lon, lat, depth, time, url):
        data = self.get_data(url)
        lons = data['lon_rho'][:]
        lats = data['lat_rho'][:]
        t = []
        for i in range(len(time)):
            print(i, time[i])
            watertemp = self.__watertemp(lon[i], lat[i], lons, lats, depth[i], time[i], data)
            t.append(watertemp)
        return t
    def __watertemp(self, lon, lat, lons, lats, depth, time, data):
        index = self.nearest_point_index2(lon, lat, lons, lats)
        depth_layers = data['h'][index[0][0]][index[1][0]]*data['s_rho']
        t = []
        # depth = depth.split(',')
        time_index = self.closest_num((time-datetime(2006,1,1)).total_seconds(), self.oceantime) -\
                self.index1
        tem = data['temp'][time_index]
        tem[tem.mask] = 10000
        for dep in depth:
            layer = np.argmin(abs(depth_layers + dep))
            temp = tem[layer, index[0][0], index[1][0]]
            t.append(temp)
            # print time, dep, temp
        return t
class water_fvcom(water):
    def __init__(self, modelname='massbay'):
        self.modelname = modelname
    def get_url(self, starttime, endtime):
        '''
        get different url according to starttime and endtime.
        urls are monthly.
        '''
        self.hours = int((endtime-starttime).total_seconds()/60/60)
        if self.modelname is "30yr":
            url = []
            time1 = datetime(year=2011,month=1,day=1)      #all these datetime are made based on the model.
            time2 = datetime(year=2011,month=11,day=11)      #The model use different version data of different period.
            time3 = datetime(year=2013,month=5,day=9)
            time4 = datetime(year=2013,month=12,day=1)
            if endtime < time1:
                yearnum = starttime.year-1981
                standardtime = datetime.strptime(str(starttime.year)+'-01-01 00:00:00',
                                                 '%Y-%m-%d %H:%M:%S')
                index1 = int(26340+35112*(yearnum/4)+8772*(yearnum%4)+1+self.hours)
                index2 = index1 + self.hours
                furl = 'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?h[0:1:48450],lat[0:1:48450],latc[0:1:90414],lon[0:1:48450],lonc[0:1:90414],u[{0}:1:{1}][0:1:44][0:1:90414],v[{0}:1:{1}][0:1:44][0:1:90414],siglay'
                url.append(furl.format(index1, index2)) 
            elif time1 <= endtime < time2: # endtime is in GOM3_v11
                url.extend(self.__chain(starttime,endtime,time1,time2))
            elif time2 <= endtime < time3:  # endtime is in GOM3_v12
                url.extend(self.__chain(starttime,endtime,time2,time3))
            elif time3 <= endtime < time4:
                url.extend(self.__chain(starttime,endtime,time3,time4))
        elif self.modelname is "GOM3":
            url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc?lon[0:1:51215],lat[0:1:51215],lonc[0:1:95721],latc[0:1:95721],siglay[0:1:39][0:1:51215],h[0:1:51215],u[{0}:1:{1}][0:1:39][0:1:95721],v[{0}:1:{1}][0:1:39][0:1:95721]'
            period = starttime-\
                     (datetime.now().replace(hour=0,minute=0)-timedelta(days=3))
            index1 = int(period.total_seconds()/60/60)
            print('index1', index1)
            index2 = index1 + self.hours
            url = url.format(index1, index2)
        elif self.modelname is "massbay":
            url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc?lon[0:1:98431],lat[0:1:98431],lonc[0:1:165094],latc[0:1:165094],siglay[0:1:9][0:1:98431],h[0:1:98431],u[{0}:1:{1}][0:1:9][0:1:165094],v[{0}:1:{1}][0:1:9][0:1:165094]'
            period = starttime-\
                     (datetime.now().replace(hour=0,minute=0)-timedelta(days=3))
            index1 = int(period.total_seconds()/60/60)
            index2 = index1 + self.hours
            url = url.format(index1, index2)
        return url
    def __chain(self, starttime, endtime, time1, time2):
        if time1 <= endtime < time2:
            pass
        else:
            sys.exit('{0} not in the right period'.format(endtime))
        url = []
        if starttime >= time1:    #start time is from 2011.11.10 as v12
            if starttime.month == endtime.month:
                url.append(self.__url(starttime.year,starttime.month,
                                            [starttime.day,starttime.hour],
                                            [endtime.day,endtime.hour]))
            else:
                if starttime.year == endtime.year:
                    y = starttime.year
                    for i in range(starttime.month, endtime.month+1):
                        if i == starttime.month:
                            url.append(self.__url(y,i,
                                                  [starttime.day, starttime.hour],
                                                  [calendar.monthrange(y,i)[1],0]))
                        elif starttime.month < i < endtime.month:
                            url.append(self.__url(y,i,[1,0],
                                                  [calendar.monthrange(y,i)[1],0]))
                        elif i == endtime.month:
                            url.append(self.__url(y,i,[1,0],
                                                  [endtime.day,endtime.hour]))
                else:
                    for i in range(starttime.year, endtime.year+1):
                        if i == starttime.year:
                            url.extend(self.get_url(starttime,
                                               datetime(year=i,
                                                        month=12,day=31)))
                        elif i == endtime.year:
                            url.extend(self.get_url(datetime(year=i,month=1,day=1),
                                               endtime))
                        else:
                            url.extend(self.get_url(datetime(year=i,month=1,day=1),
                                               datetime(year=i,month=12,day=31)))
             
        else:
            url.extend(self.get_url(starttime,(time1-timedelta(minutes=1))))
            url.extend(self.get_url(time1,endtime))
        return url
    def __url(self, year, month, start_daytime, end_daytime):
        '''
        start_daytime,end_daytime: [day,hour]
        '''
        url_v11 = 'http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/NECOFS_GOM3_{0}/gom3v11_{0}{1}.nc?lon[0:1:48727],lat[0:1:48727],lonc[0:1:90997],latc[0:1:90997],h[0:1:48727],u[{2}:1:{3}][0:1:39][0:1:90997],v[{2}:1:{3}][0:1:39][0:1:90997],siglay[0:1:39][0:1:48727]'
        url_v12 = 'http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/NECOFS_GOM3_{0}/gom3v12_{0}{1}.nc?lon[0:1:48859],lat[0:1:48859],lonc[0:1:91257],latc[0:1:91257],h[0:1:48859],u[{2}:1:{3}][0:1:39][0:1:91257],v[{2}:1:{3}][0:1:39][0:1:91257],siglay[0:1:39][0:1:48859]'
        url_v13 = 'http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/NECOFS_GOM3_{0}/gom3v13_{0}{1}.nc?lon[0:1:51215],lat[0:1:51215],lonc[0:1:95721],latc[0:1:95721],h[0:1:51215],u[{2}:1:{3}][0:1:39][0:1:95721],v[{2}:1:{3}][0:1:39][0:1:95721],siglay[0:1:39][0:1:51215]'
        time1 = datetime(year=2011,month=1,day=1)      #all these datetime are made based on the model.
        time2 = datetime(year=2011,month=11,day=11)      #The model use different version data of different period.
        time3 = datetime(year=2013,month=5,day=9)
        time4 = datetime(year=2013,month=12,day=1)
        currenttime = datetime(year=year,month=month,day=start_daytime[0])
                                       
        if time1 <= currenttime < time2:
            version = '11'
        elif time2 <= currenttime < time3:
            version = '12'
        elif time3 <= currenttime < time4:
            version = '13'

        if year == 2011 and month == 11  and start_daytime[0] >10:
            start = str(24*(start_daytime[0]-1)+start_daytime[1]-240)
            end = str(24*(end_daytime[0]-1)+end_daytime[1]-240)
        elif year == 2013 and month == 5 and start_daytime[0] >8:
            start = str(24*(start_daytime[0]-1)+start_daytime[1]-192)
            end = str(24*(end_daytime[0]-1)+end_daytime[1]-192)
        else:
            start = str(24*(start_daytime[0]-1)+start_daytime[1])
            end = str(24*(end_daytime[0]-1)+end_daytime[1])
        year = str(year)
        month = '{0:02d}'.format(month)
        
        if version == '11':
            url = url_v11.format(year, month, start, end)
        elif version == '12':
            url = url_v12.format(year, month, start, end)
        elif version == '13':
            url = url_v13.format(year, month, start, end)
        return url
    def get_data(self,url):
        self.data = jata.get_nc_data(url,'lon','lat','latc','lonc',
                                     'u','v','siglay','h')
        return self.data
    def waternode(self, lon, lat, depth, url):
        if type(url) is str:
            nodes = dict(lon=[lon],lat=[lat])
            temp = self.__waternode(lon, lat, depth, url)
            nodes['lon'].extend(temp['lon'])
            nodes['lat'].extend(temp['lat'])
        else:
            nodes = dict(lon=[lon],lat=[lat])
            for i in url:
                temp = self.__waternode(nodes['lon'][-1], nodes['lat'][-1], depth, i)
                nodes['lat'].extend(temp['lat'])
                nodes['lon'].extend(temp['lon'])
        return nodes
    def __waternode(self, lon, lat, depth, url):
        '''
        start, end: indices of some period
        data: a dict that has 'u' and 'v'
        '''
        data = self.get_data(url)
        lonc, latc = data['lonc'][:], data['latc'][:]
        lonv, latv = data['lon'][:], data['lat'][:]
        h = data['h'][:]
        siglay = data['siglay'][:]
        if lon>90:
            lat, lon = dm2dd(lat, lon)
        nodes = dict(lon=[], lat=[])
        kf,distanceF = self.nearest_point_index(lon,lat,lonc,latc,num=1)
        kv,distanceV = self.nearest_point_index(lon,lat,lonv,latv,num=1)
        print('kf', kf)
        if h[kv] < 0:
            sys.exit('Sorry, your position is on land, please try another point')
        depth_total = siglay[:,kv]*h[kv]
        ###############layer###########################
        layer = np.argmin(abs(depth_total-depth))
        # for i in range(len(data['u'])):
        for i in range(self.hours):
            # u_t = np.array(data['u'])[i,layer,kf]
            # v_t = np.array(data['v'])[i,layer,kf]
            u_t = data['u'][i, layer, kf[0][0]]
            v_t = data['v'][i, layer, kf[0][0]]
            print('u_t, v_t, i', u_t, v_t, i)
            dx = 60*60*u_t
            dy = 60*60*v_t
            lon = lon + (dx/(111111*np.cos(lat*np.pi/180)))
            lat = lat + dy/111111
            nodes['lon'].append(lon)
            nodes['lat'].append(lat)
            kf, distanceF = self.nearest_point_index(lon, lat, lonc, latc,num=1)
            kv, distanceV = self.nearest_point_index(lon, lat, lonv, latv,num=1)
            # depth_total = siglay[:][kv]*h[kv]
            if distanceV>=.3:
                if i==start:
                    print('Sorry, your start position is NOT in the model domain')
                    break
        return nodes
    def watertemp(self, lon, lat, depth, time, url):
        '''
        Using __decideTimeRange() which return a string:
        'v0': 1978.01.01 ~ 2011.01.01
        'v1': 2011.01.01 ~ 2011.11.11
        'v2': 2011.11.11 ~ 2013.05.09
        'v3': 2013.05.09 ~ 2013.12.01
        not finished
        '''
        url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?lat[0:1:48450],lon[0:1:48450],siglay[0:1:44][0:1:48450],h[0:1:48450]'
        data = netCDF4.Dataset(url)
        lons = data.variables['lon']
        lats = data.variables['lat']
        siglay = data.variables['siglay']
        h = data.variables['h']
        t = []
        for i in range(len(time)):
            v = self.__decideTimeRange(time[i])
            index, = self.nearest_point_index(lon, lat, lons, lats, num=1)
            #layer = closest_num(depth[i], h[index[i]]*siglay[])
            if v == 'v0':
                indTime = time[i] - datetime(1978, 1, 1)/3600
                url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?temp[indTime:1:indTime+1][layer:1:layer+1][index[0]:1:index[0]+1]'
                d = netCDF4.Dataset(url).variables['temp'][0,0,0]
            elif v == 'v1':
                indTime = time[i] - datetime()
def angle_conversion(a):
    a = np.array(a)
    return a/180*np.pi
def dist(lon1, lat1, lon2, lat2):
    R = 6371.004
    lon1, lat1 = angle_conversion(lon1), angle_conversion(lat1)
    lon2, lat2 = angle_conversion(lon2), angle_conversion(lat2)
    l = R*np.arccos(np.cos(lat1)*np.cos(lat2)*np.cos(lon1-lon2)+\
                        np.sin(lat1)*np.sin(lat2))
    return l
def mon_alpha2num(m):
    month = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
    if m in month:
        n = month.index(m)
    else:
        raise Exception('Wrong month abbreviation')
    return n+1
def np_datetime(m):
    if type(m) is str:
        year = int(m[5:9])
        month = mon_alpha2num(m[2:5])
        day =  int(m[0:2])
        hour = int(m[10:12])
        minute = int(m[13:15])
        second = int(m[-2:])
        dt = datetime(year,month,day,hour=hour,minute=minute,second=second)
    else:
        dt = []
        for i in m:
            year = int(i[5:9])
            month = mon_alpha2num(i[2:5])
            day =  int(i[0:2])
            hour = int(i[10:12])
            minute = int(i[13:15])
            second = int(i[-2:])
            temp = datetime(year,month,day,hour=hour,minute=minute,second=second)
            dt.append(temp)
            dt = np.array(dt)
    return dt
def get_nc_data(url, *args):
    '''
    get specific dataset from url

    *args: dataset name, composed by strings
    ----------------------------------------
    example:
        url = 'http://www.nefsc.noaa.gov/drifter/drift_tcs_2013_1.dat'
        data = get_url_data(url, 'u', 'v')
    '''
    nc = netCDF4.Dataset(url)
    data = {}
    for arg in args:
        try:
            data[arg] = nc.variables[arg]
        except (IndexError, NameError, KeyError):
            print('Dataset {0} is not found'.format(arg))
    return data