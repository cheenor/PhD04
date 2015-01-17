#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 20 21:03:22 2014

@author: jhchen
"""
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib.pyplot as plt
import calendar
import datetime 
import time
import math
dirin="D:/MyPaper/PhD04/CERES/"
dirout="D:/MyPaper/PhD04/Data/CERES/"
iyr=2010

rgns=["ETP","WTP","PRD","MLYR","NPC","NEC"]
imm =[  6,    6,    6,   6,     8,    7  ]
idd =[  4,    1,    1,   1,     2,    1  ]
ndds=[  30,   30,   30,  30,    30,   30 ]
lon1=[  90,   80,   110,  110,  112,  120 ]
lon2=[  100,  90,   118,  122,  120,  130 ]
lat1=[  27.5,  27.5,  27.5, 27, 34,   43  ]
lat2=[  37.5,  37.5,  35,  33,  42,   49  ]
nr=len(rgns)
filenm=["CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed3A_Subset_201004-201010_TOA_OBS.nc",
    "CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed3A_Subset_201004-201010_TOA_CMP.nc",
    "CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed3A_Subset_201004-201010_SRF_CMP.nc"]
ima=4
ida=1
ystr='%d'%iyr
mstr='%d'%ima
dstr='%d'%ida
text=ystr+"-"+mstr+"-"+dstr
d = datetime.datetime.strptime(text,'%Y-%m-%d')
daynumdata=d.timetuple().tm_yday
fpath=[dirin+filenm[0],
       dirin+filenm[1],
        dirin+filenm[2]]
for i in range(0,len(rgns)):  ### len(rgns)
    ims=imm[i]
    ids=idd[i]
    ystr='%d'%iyr
    mstr='%d'%ims
    dstr='%d'%ids
    text=ystr+"-"+mstr+"-"+dstr
    d = datetime.datetime.strptime(text,'%Y-%m-%d')
    text=ystr+"_"+mstr+"_"+dstr
    daynumt=d.timetuple().tm_yday
    rnst=daynumt-daynumdata
    ndd=ndds[i]
    nddstr='%d'%ndd
    lonw=lon1[i]
    lone=lon2[i]
    lats=lat1[i]
    latn=lat2[i]
    for ifl in range(0,len(fpath)):
        fout=dirout+rgns[i]+"_"+filenm[ifl][61:68]+"_"+text+"__"+nddstr+"d.txt"
        f=open(fout,'w')
        toaobs=Dataset(fpath[ifl],'a')
        lon= toaobs.variables['lon'][:]   # read 
        tmstp=toaobs.variables['time'][:]
        lat=toaobs.variables['lat'][:]
        varnms=[]
        ondim=[]
        for nm in toaobs.variables:
            if len(nm)>4:        
                varnms.append(nm)            
                tmpvar=toaobs.variables[nm][:] ## time lat lon
                nt=len(tmpvar[:,0,0])
                nlat=len(tmpvar[0,:,0])
                nlon=len(tmpvar[0,0,:])
                for it in range(rnst,rnst+ndd):
                    tmpavg=0.
                    ixy=0
                    for iy in range(0,nlat):
                        if lat[iy]>lats and lat[iy]<latn:
                            for ix in range(0,nlon):
                                if lon[ix]>lonw and lon[ix]<lone:
                                    if tmpvar[it,iy,ix]:  # skip the masked value
                                        tmpavg=tmpavg+tmpvar[it,iy,ix]
                                        ixy=ixy+1
                    if ixy==0:
                        tmpavg=-999.
                    if ixy !=0:
                        tmpavg=tmpavg/ixy
                    if tmpavg:
                        tmpavg=tmpavg
                    else: 
                        tmpavg=-999.
                    ondim.append(tmpavg)
        slen_var=len(ondim)
        lenv=len(varnms)
        for nmstr in varnms:
            head="%s "%nmstr
            f.write(head)
        f.write('\n')
        for ii in range(0,ndd):
            for ij in range(0,lenv):                  
                itme="%f "%ondim[ii+ij*ndd]
                f.write(itme)
            f.write('\n')
        f.close()
            
            
            
            
            
            
    