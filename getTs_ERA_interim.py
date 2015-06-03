#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 15:39:20 2015

@author: jhchen
"""
from netCDF4 import Dataset
import datetime 
import os
import calendar
os.system("cls")
#
dirin='X:/Data/ERA_interim/SRFX2.5/'
dirout='X:/Data/ERA_interim/SHLH/'
years =range(1979,2015,4)
rgnm=[]
slon=[]
elon=[]
slat=[]
elat=[]
#-----PRD----------
rgnm.append('PRD')
slon.append(107.5) ; elon.append(117.5)
slat.append(17.5)  ;  elat.append(27.5)
#-----MLYR----------
rgnm.append('MLYR')
slon.append(110.) ;  elon.append(120.)
slat.append(25.)  ;  elat.append(35.)
#-----NPC----------
rgnm.append('NPC')
slon.append(110.) ;  elon.append(120.)
slat.append(32.5)  ;  elat.append(42.5)
#-----NEC----------
rgnm.append('NEC')
slon.append(120.) ;  elon.append(130.)
slat.append(40.)  ;  elat.append(50.)
#-----WTP----------
rgnm.append('WTP')
slon.append(80.) ;  elon.append(90.)
slat.append(27.5)  ;  elat.append(37.5)
#-----ETP----------
rgnm.append('ETP')
slon.append(90.) ;  elon.append(100.)
slat.append(27.5)  ;  elat.append(37.5)
#
ngns=len(rgnm)
fpath=dirin+'1983-1986_ts_sst.nc'
f=Dataset(fpath,'a')
for a in f.variables:
    print a
lon= f.variables['longitude'][:]   # read 
tmstp=f.variables['time'][:]
lat=f.variables['latitude'][:]
skt=f.variables['skt'][:] # nt,ny,nx
nx=len(lon)
ny=len(lat)
for yy in years:
    ystr1="%d"%yy
    ystr2="%d"%(yy+3)
    filename=ystr1+"-"+ystr2+"_ts_sst.nc"
    fpath=dirin+filename
    f=Dataset(fpath,'a')
    skt=f.variables['skt'][:]
    sst=f.variables['sst'][:]
    itds=0
    itde=0
    for iy in range(0,4):
        iyr=iy+yy
        ndt=365*4
        if calendar.isleap(iyr):
            ndt=366*4
        itde=itds+ndt
        namestr="%d"%iyr
        for ig in range(0,ngns):
            fpath=dirout+namestr+rgnm[ig]+"_SST_SKT.txt"
            fout=open(fpath,"w")
            item="Time "
            fout.write(item)
            item="SST "
            fout.write(item)
            item="SKT"
            fout.write(item)            
            fout.write('\n')
            itcont=1
            for itt in range(itds,itde):
                cont1=0.
                cont2=0.
                tmp1=0.0
                tmp2=0.0
                item="%d "%itcont
                fout.write(item)
                itcont=itcont+1
                for ix in range(0,nx):
                    if lon[ix]>(slon[ig]) and lon[ix]<(elon[ig]):
                        for iy in range(0,ny):
                            if lat[iy]>slat[ig] and lat[iy]<elat[ig]:
                                if sst[itt,iy,ix]>0. :
                                    tmp1=tmp1+sst[itt,iy,ix]
                                    cont1=cont1+1.
                                if skt[itt,iy,ix]>0. :
                                    tmp2=tmp2+skt[itt,iy,ix]
                                    cont2=cont2+1.
                if cont1 > 0. :
                    tmp1=tmp1/cont1
                else:
                    tmp1=-99.
                if cont2 >0. :
                    tmp2=tmp2/cont2
                else:
                    tmp2=-99.
                item="%f "%tmp1
                fout.write(item)
                item="%f "%tmp2
                fout.write(item)
                fout.write('\n')
            fout.close()
        idts=itde
    f.close()
                
        