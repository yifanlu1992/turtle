'''
If the deepest observation depth is >50m(or <50m, or all), draw the correlation of this observation and appriate model data.
'''
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import netCDF4
import matplotlib.pyplot as plt
from scipy import stats
from watertempModule import water_roms
import watertempModule as wtm
from turtleModule import mon_alpha2num, np_datetime, mean_value, bottom_value, index_by_depth,str2ndlist
def show2pic(x1, y1, fontsize):
    FONTSIZE = fontsize
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    x = np.arange(0.0, 30.0, 0.01)
    '''
    for i in range(10):
        # ax.plot(temp[index[i]], obstemp[index[i]], '.', color=colors[i], label='{0}'.format(i))
        ax.scatter(temp[index[i]], obsdata['temp'][index[i]], s=50, c=colors[i], label='{0}'.format(i))
    '''
    # ax.scatter(temp[index[0]], obsdata['temp'][index[0]], s=50, c='b', label='<45')
    # ax.scatter(temp[index[1]], obsdata['temp'][index[1]], s=50, c='r', label='>=45')
    ax1.scatter(x1, y1, s=50, c='b')
    ax1.plot(x, x, 'r-', linewidth=2)
    plt.axis([0, 30, 0, 30], fontsize=15)
    plt.xlabel('Model (degC)', fontsize=FONTSIZE)
    plt.ylabel('Observed (degC)', fontsize=FONTSIZE)
    i = x1[x1.isnull()==False].index
    fit = np.polyfit(x1[i], y1[i], 1)
    fit_fn = np.poly1d(fit)
    x2, y2 = x1[i], fit_fn(x1[i])
    plt.plot(x2, y2,'y-', linewidth=2)
    gradient, intercept, r_value, p_value, std_err = stats.linregress(y1[i], x1[i])
    r_squared = r_value**2
    # ax1.set_title('R-squard: %.4f' % r_squared, fontsize=FONTSIZE)
    plt.savefig('obsVSmodelDeepestBottom1.png',dpi=200)

    fig2 = plt.figure()
    ax2 =  fig2.add_subplot(111)
    nbins = 200
    H, xedges, yedges = np.histogram2d(x1[i], y1[i], bins=nbins)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0, H)
    xedges,yedges=np.meshgrid(xedges, yedges)
    plt.pcolormesh(xedges, yedges, Hmasked)
    plt.xlabel('Model (degC)', fontsize=FONTSIZE)
    plt.ylabel('Observed (degC)', fontsize=FONTSIZE)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Quantity', fontsize=FONTSIZE)
    cbar.ax.tick_params(labelsize=20)
    cbar.ax.set_yticks(fontsize=20)
    plt.axis([0, 30, 0, 30])
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.plot(x, x, 'r-', linewidth=2)
    plt.plot(x2, y2, 'y-', linewidth=2)
    # plt.title('R-squard: %.4f' % r_squared, fontsize=FONTSIZE)
    plt.savefig('obsVSmodelDeepestBottom2.png',dpi=200)
    return ax1, ax2, r_squared
#############################MAIN CODE###########################################
FONTSIZE = 14
# obs = pd.read_csv('ctd_extract_TF.csv')
obs = pd.read_csv('ctdWithModTempByDepth.csv')
obstime1 = np_datetime(obs['END_DATE'])
tf_index = np.where((obs['TF'].notnull()) & (obstime1>=datetime(2009,10,11,2,0)))[0]
obslat, obslon = obs['LAT'][tf_index].values, obs['LON'][tf_index].values
obstime = np_datetime(obs['END_DATE'][tf_index])
# obsdepth = mean_value(obs['TEMP_DBAR'][tf_index])
obsdepth = obs['MAX_DBAR'][tf_index].values
# obstemp = mean_value(obs['TEMP_VALS'][tf_index])
obstemp = bottom_value(obs['TEMP_VALS'][tf_index])
obsdata = pd.DataFrame({'depth':obsdepth, 'temp':obstemp, 'lon':obslon,
                        'lat':obslat, 'time':obstime}).sort_index(by='depth')

starttime = datetime(2009,10,11,2,0)# altough our starttime is(2009, 8, 24),but the roms time start on 2009,10,11
endtime = datetime(2013,12,13)
tempobj = wtm.water_roms()
url = tempobj.get_url(starttime, endtime)
# temp = tempobj.watertemp(obslon, obslat, obsdepth, obstime, url)
obsdata_time=[]# if we use the obsdata['time'].values in the below , we will get the datetime64,but we just want to get the datetime
for i in obsdata['time']:
    obsdata_time.append(i)
#temp = tempobj.watertemp(obsdata['lon'].values, obsdata['lat'].values,
                         #obsdata['depth'].values, obsdata_time, url)
#temp = pd.Series(temp, index = obsdata['temp'].index)
temp1 = pd.Series(str2ndlist(obs['modTempByDepth'][tf_index],bracket=True), index=tf_index)
temp1.ix[obsdata.index]
temp=[]
for i in temp1:
    temp.append(i[-1])
temp=pd.Series(temp,index=obsdata.index)

index = index_by_depth(obsdata['depth'], 50)
# colors = utilities.uniquecolors(10)
tp='all'
if tp == 'all':
    x1, y1 = temp, obsdata['temp']
    ax1, ax2, r_squared = show2pic(x1, y1, FONTSIZE)
    ax1.set_title('Bottom Temperature(R-squared=%.2f)' % r_squared, fontsize=FONTSIZE)
    ax2.set_title('Bottom Temperature(R-squared=%.2f)' % r_squared, fontsize=FONTSIZE)
elif tp == '<50':
    x1, y1 = temp[index[0]], obsdata['temp'][index[0]]
    ax1, ax2, r_squared = show2pic(x1, y1, FONTSIZE)
    ax1.set_title('%s, R-squared: %.4f' % (tp, r_squared), fontsize=FONTSIZE)
    ax2.set_title('%s, R-squared: %.4f' % (tp, r_squared), fontsize=FONTSIZE)
elif tp == '>50':
    x1, y1 = temp[index[1]], obsdata['temp'][index[1]]
    ax1, ax2, r_squared = show2pic(x1, y1, FONTSIZE)
    ax1.set_title('%s, R-squared: %.4f' % (tp, r_squared), fontsize=FONTSIZE)
    ax2.set_title('%s, R-squared: %.4f' % (tp, r_squared), fontsize=FONTSIZE)
plt.show()
'''
# Plot Deepest Data Quantity
fig = plt.figure()
ax = fig.add_subplot(111)
y = obsdata['depth'].values
x = np.arange(1, np.amax(y)+1)
bar = np.array([0]*np.amax(y))
for i in y:
    if i in x:
        bar[i-1] = bar[i-1]+1
plt.barh(x, bar)
plt.ylim([250, 0])
plt.ylabel('depth', fontsize=25)
plt.xlabel('Quantity', fontsize=25)
plt.title('Deepest data histogram', fontsize=25)
plt.show()
'''
