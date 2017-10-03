# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 13:51:36 2017
draw the correlation between turtle and models(fvcom,roms,hycom) and ship 
@author: yifan
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from scipy import stats
from turtleModule import str2ndlist, np_datetime, bottom_value, closest_num
def histogramPoints(x, y, bins):
    H, xedges, yedges = np.histogram2d(x, y, bins=bins)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0, H)
    return xedges, yedges, Hmasked
def get_max_color_value(Hmasked,bins):
    m=0
    for i in range(bins):
        for j in Hmasked[i]:
            if j>m:
               m=j
    return m
def show2pic(rng,tempobs, tempmod,fontsize,text,bins):
    n = np.arange(0, 30, 0.01)
    fig = plt.figure(figsize=[12,12])
    vmax=0
    for i in range(len(text)):
        tempObs=tempobs[i]
        tempMod=tempmod[i]
        x, y ,Hmasked = histogramPoints(tempObs, tempMod, bins) 
        m=get_max_color_value(Hmasked,bins)
        if m>vmax:
           vmax=m
    for j in range(len(text)):
        print(j)
        ax = fig.add_subplot(2,2,j+1)
        tempObs=tempobs[j]
        tempMod=tempmod[j]
        x, y ,Hmasked = histogramPoints(tempObs, tempMod, bins)
        c = ax.pcolormesh(x, y, Hmasked,vmin=0,vmax=vmax)#,vmin=number[i][0],vmax=number[i][1] 
        ax.plot(n, n, 'r-')
        fit = np.polyfit(tempObs, tempMod, 1)
        fit_fn = np.poly1d(fit)
        ax.plot(tempObs, fit_fn(tempObs), 'y-', linewidth=2)
        gradient, intercept, r_value, p_value, std_err = stats.linregress(tempObs, tempMod)
        r_squared=r_value**2
        ax.set_title('%s'%text[j],fontsize=15)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.ylim([-5,35])
        plt.text(x=1,y=31,s=r'$\mathregular{R^2}$='+str(round(r_squared,3)),fontsize=15) 
        #cbar = plt.colorbar(c)
        #cbar.ax.tick_params(labelsize=10)
        if j==1 or j==3: 
           plt.setp(ax.get_yticklabels() ,visible=False)
        if j==0 or j==1:
           plt.setp(ax.get_xticklabels() ,visible=False)
    cax = fig.add_axes([0.91, 0.12, 0.025, 0.78])
    fig.colorbar(c, cax=cax)#, cax=cax
    fig.text(0.5, 0.06, 'Observed temperature($^\circ$C)', ha='center', va='center', fontsize=FONTSIZE)
    fig.text(0.06, 0.5, 'Model temperature($^\circ$C)', ha='center', va='center', rotation='vertical',fontsize=FONTSIZE)
    fig.text(0.5, 0.94, '%s' %rng[0], ha='center', va='center', fontsize=FONTSIZE)
    fig.text(0.97, 0.5, 'Quantity', ha='center', va='center', rotation='vertical',fontsize=FONTSIZE)
    #fig.tight_layout()
    plt.savefig('correlation comparison_bottom.png', dpi=200)
#########################################MAIN CODE#####################################################################################################
FONTSIZE = 25
##################################ship#########################################
obsData=pd.read_csv('matched_turtleVSship_allship.csv')
obsdepth=pd.Series(str2ndlist(obsData['turtle_depth'], bracket=True))
obstemp=pd.Series(str2ndlist(obsData['turtle_temp'], bracket=True))
shipdepth=pd.Series(str2ndlist(obsData['ship_depth'], bracket=True))
shiptemp=pd.Series(str2ndlist(obsData['ship_temp'], bracket=True))
tempObs_ship=[] 
tempMod_ship=[]
for i in obsData.index:
    
    for k in range(len(shipdepth[i])):
        if obsdepth[i][-1]==shipdepth[i][k]:
            tempObs_ship.append(obstemp[i][-1])
            tempMod_ship.append(shiptemp[i][k])

##################################roms#########################################
obs = pd.read_csv('ctdWithModTempByDepth.csv') 
tf_index = np.where(obs['TF'].notnull())[0]    # Get the index of good data.
obsTemp = pd.Series(str2ndlist(obs['TEMP_VALS'][tf_index]), index=tf_index)
obsDepth = pd.Series(str2ndlist(obs['TEMP_DBAR'][tf_index]), index=tf_index)
tempMod = pd.Series(str2ndlist(obs['modTempByDepth'][tf_index],bracket=True), index=tf_index)
tempObs_roms = []
tempMod_roms = []
for i in range(len(obsTemp.values)):
    if tempMod.values[i][-1]>100:continue
    tempObs_roms.append(obsTemp.values[i][-1])
    tempMod_roms.append(tempMod.values[i][-1])
    '''for j in range(len(obsDepth.values[i])):
        d = obsDepth.values[i][j]
        if tempMod.values[i][j] > 100: continue
        
        if d>0:
            tempObs_roms.append(obsTemp.values[i][j])
            tempMod_roms.append(tempMod.values[i][j])'''
        
##################################fvcom########################################
obs = pd.read_csv('ctdWithdepthofbottom_fvcom.csv') 
tf_index = np.where(obs['in FVcom range'].notnull())[0]    # Get the index of good data.
obsTemp = pd.Series(str2ndlist(obs['TEMP_VALS'][tf_index]), index=tf_index)
obsDepth = pd.Series(str2ndlist(obs['TEMP_DBAR'][tf_index]), index=tf_index)
tempMod = pd.Series(str2ndlist(obs['modtempBYdepth'][tf_index],bracket=True), index=tf_index)
tempObs_fvcom = []
tempMod_fvcom = []
for i in range(len(obsTemp.values)):
    tempObs_fvcom.append(obsTemp.values[i][-1])
    tempMod_fvcom.append(tempMod.values[i][-1])
    '''for j in range(len(obsDepth.values[i])):
        d = obsDepth.values[i][j]
        if tempMod.values[i][j] > 100: continue
        
        if d>0:
            tempObs_fvcom.append(obsTemp.values[i][j])
            tempMod_fvcom.append(tempMod.values[i][j])'''
        
##################################hycom########################################
obs = pd.read_csv('ctd_withHYCOMtemp.csv') 
tf_index = np.where(obs['TF'].notnull())[0]    # Get the index of good data.
obsTemp = pd.Series(str2ndlist(obs['TEMP_VALS'][tf_index]), index=tf_index)
obsDepth = pd.Series(str2ndlist(obs['TEMP_DBAR'][tf_index]), index=tf_index)
tempMod = pd.Series(str2ndlist(obs['modtemp_HYCOM'][tf_index],bracket=True), index=tf_index)
tempObs_hycom = []
tempMod_hycom = []
for i in range(len(obsTemp.values)):
    if tempMod.values[i][-1]<-10:continue
    if tempMod.values[i][-1]>0: # here will just delete one value(-2.6550007)
        tempObs_hycom.append(obsTemp.values[i][-1])
        tempMod_hycom.append(tempMod.values[i][-1])
    '''for j in range(len(obsDepth.values[i])):
        d = obsDepth.values[i][j]
        if tempMod.values[i][j] < -10: continue
        
        if d>0:
            tempObs_hycom.append(obsTemp.values[i][j])
            tempMod_hycom.append(tempMod.values[i][j])'''
        
tempObs=[tempObs_ship ,tempObs_roms,tempObs_fvcom,tempObs_hycom]
tempMod=[tempMod_ship,tempMod_roms,tempMod_fvcom,tempMod_hycom]
#rng = ['25.0', '25.0~50.0', '50.0~75.0', '75.0','Using the entire profiles','<50']
rng=['bottomtemp_correlation comparison']
#number=[[0,200],[0,25],[0,12],[0,5],[0,200],[0,200]]
text=['SHIP','ROMS','FVCOM','HYCOM']
show2pic(rng,tempObs,tempMod,FONTSIZE,text,bins=150)
plt.show()
