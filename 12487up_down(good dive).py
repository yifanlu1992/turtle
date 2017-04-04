# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 13:48:49 2016
find out downcast and upcast for 12487 turtle that dive from 10% to 90% of ocean depth
@author: yifan
"""
import numpy as np
import pandas as pd
from gettopo import gettopo
def find_start(a,List):#find one dive behind of "a" in list
    start=len(List)-3
    for i in range(a,len(List)-2):
        if List[i]<List[i+1] and List[i]<2:
            if List[i+1]<List[i+2] and List[i+1]>=2 :
                start=i
                break                            #find farest on surface before diving
    for i in range(start,len(List)-2):
        if List[i]>=2 and List[i]>List[i+1]:
           if List[i+1]<2 and List[i+1]>List[i+2]:
               break                             #find nearest on surface after diving
    return [start,i+1]
    
def return_rawindex(oldindex,newlist):
    oldlist=[]
    for i in range(len(newlist)):
        va=[]
        for k in range(len(newlist[i])):
            for j in range(len(oldindex)):
                if j==newlist[i][k]:
                    va.append(oldindex[j])
        oldlist.append(va)      
    return oldlist            
#####################################
data=pd.read_csv('12487_location.csv')# /net/data5/jmanning/turtle/
inde=np.where(data['index'].notnull())[0] #delete the data without lan and lon
depth=pd.Series(data['depth'][inde],index=inde)

turtle_index=data['index'][inde]
lon=data['lon'][inde] 
lat=data['lat'][inde]
new_index=range(len(inde))
turtle_index.index=new_index
lon.index=new_index
lat.index=new_index
depth.index=new_index
'''turtle_index=pd.Series(data['index'][inde],index=new_index)
lon=pd.Series(data['lon'][inde],index=new_index) 
lat=pd.Series(data['lat'][inde],index=new_index)
 '''                        
water_depth=[]  # get water depth according to lat and lon
for i in new_index:
    if i%10000==0:   
        print(i)
    if i==0 or i==len(new_index)-1:
        wd1=-gettopo(lat[i],lon[i]) # this is a function to get water depth from a USGS database
    elif turtle_index[i-1]!=turtle_index[i+1]:
        wd1=-gettopo(lat[i],lon[i])
    water_depth.append(wd1)
wd=pd.Series(water_depth,index=new_index)

indx=[] #find the start of downcast
for n in new_index:
    if depth[n]<=wd[n]*0.1:      # use 2 beacuse some dives don`t go to air
            I=find_start(n,depth)
            indx.append(I)
    if n%1000==0:    
        print(I,n) 
Indx = [list(x) for x in set(tuple(x) for x in indx)]  #get rid of repeated data
Indx.sort()

INDX=[]  #find the deepest dives 
for i in range(len(Indx)):
    Max=max(depth[Indx[i][0]:Indx[i][1]+1])
    Max_indx=[]
    for j in range(Indx[i][0],Indx[i][1]+1):
        if int(depth[j])==int(Max):
            if  depth[j]>=wd[j]*0.9:
                Max_indx.append(j)    
    INDX.append(Max_indx)
    #print('i',i)
newIndx=return_rawindex(inde,Indx)
newINDX=return_rawindex(inde,INDX)   

va_index=[]
for i in range(len(newINDX)):
    if newINDX[i]!=[]:
       va_index.append(i)
downcast,upcast=[],[]
for i in va_index:
    downcast.append([newIndx[i][0],newINDX[i][0]])
    upcast.append([newINDX[i][-1],newIndx[i][-1]])
data=pd.DataFrame(range(len(upcast)))
data['down']=pd.Series(downcast)
data['up']=pd.Series(upcast)
data.to_csv('12487up_down(good position).csv')
