#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 20 21:03:22 2014

@author: jhchen
"""
from netCDF4 import Dataset
import datetime 
import os
os.system("cls")

dirin="D:/MyPaper/PhD04/CERES/"
dirout="D:/MyPaper/PhD04/Data/CERES/"
iyr=2010

rgns=["ETP","WTP","PRD","MLYR","NPC","NEC"]
iyrs=[2012, 2010, 2012, 2010,   2010, 2012]
imm =[  5,    7,    4,   6,     8,    7  ]
idd =[  20,    14,    1,   24,  2,    6  ]
ndds=[  30,   30,   30,  30,    30,   30 ]
lon1=[  90,   80,   110,  110,  112,  120 ]
lon2=[  100,  90,   118,  122,  120,  130 ]
lat1=[  27.5,  27.5,  27.5, 27, 34,   43  ]
lat2=[  37.5,  37.5,  35,  33,  42,   49  ]
filenm0=["CERES_SYN1deg-3H_Terra-Aqua-MODIS_Ed3A_Subset_20100401-20100930.nc",
    "CERES_SYN1deg-3H_Terra-Aqua-MODIS_Ed3A_Subset_20100401-20100930_OBS_TOA.nc"] #,
#    "CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed3A_Subset_201004-201010_SRF_CMP.nc"]
filenm1=["CERES_SYN1deg-3H_Terra-Aqua-MODIS_Ed3A_Subset_20120401-20120930.nc",
    "CERES_SYN1deg-3H_Terra-Aqua-MODIS_Ed3A_Subset_20120401-20120930_OBS_TOA.nc"] 
ima=4
ida=1
#,
#        dirin+filenm[2]]
for i in range(0,len(rgns)):  ### len(rgns)
    if iyrs[i] == 2010:
        filenm=filenm0
    elif iyrs[i] == 2012 :
        filenm=filenm1        
    iyr=iyrs[i]    
    ystr='%d'%iyr
    mstr='%d'%ima
    dstr='%d'%ida
    text=ystr+"-"+mstr+"-"+dstr
    d = datetime.datetime.strptime(text,'%Y-%m-%d')
    daynumdata=d.timetuple().tm_yday    
    fpath=[dirin+filenm[0],
           dirin+filenm[1]]
    ims=imm[i]
    ids=idd[i]
    ystr='%d'%iyr
    mstr='%d'%ims
    dstr='%d'%ids
    text=ystr+"-"+mstr+"-"+dstr
    d = datetime.datetime.strptime(text,'%Y-%m-%d')
    text=ystr+"_"+mstr+"_"+dstr
    daynumt=d.timetuple().tm_yday
    rnst=daynumt*8-daynumdata*8 ####!!!! daynumt is day, every 3 hours one data
    print rnst
    ndd=ndds[i]
    nddstr='%d'%ndd
    lonw=lon1[i]
    lone=lon2[i]
    lats=lat1[i]
    latn=lat2[i]
    for ifl in range(0,len(fpath)):
        if ifl ==0:
            fout=dirout+rgns[i]+"_"+text+"__"+nddstr+"d3h.txt"
        if ifl==1:
            fout=dirout+rgns[i]+"_"+text+"__"+nddstr+"d3h_TOA_OBS.txt"
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
                for it in range(rnst,rnst+ndd*8+1): #### dataset is every 3 hours
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
        for ii in range(0,ndd*8+1):            
            for ij in range(0,lenv):                  
                itme="%f "%ondim[ii+ij*(ndd*8+1)]
                f.write(itme)
            f.write('\n')
        f.close()
### get the upwoad lonfwave from ncep on gausee grid  4 times oneday  
filenm=["ulwrf.ntat.gauss.2010.nc"] #,
#    "CERES_SYN1deg-3H_Terra-Aqua-MODIS_Ed3A_Subset_20100501-20100930_OBS_TOA.nc"] #,
#    "CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed3A_Subset_201004-201010_SRF_CMP.nc"]
ima=1
ida=1
ystr='%d'%iyr
mstr='%d'%ima
dstr='%d'%ida
text=ystr+"-"+mstr+"-"+dstr
d = datetime.datetime.strptime(text,'%Y-%m-%d')
daynumdata=d.timetuple().tm_yday
fpath=[dirin+filenm[0]] #,
#       dirin+filenm[1]] #,
#        dirin+filenm[2]]
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
    rnst=daynumt*4-daynumdata*4 ####!!!! daynumt is day, every 3 hours one data
    print rnst
    ndd=ndds[i]
    nddstr='%d'%ndd
    lonw=lon1[i]
    lone=lon2[i]
    lats=lat1[i]
    latn=lat2[i]
    for ifl in range(0,len(fpath)):
        if ifl ==0:
            fout=dirout+rgns[i]+"_"+text+"__"+nddstr+"d6h_NCEP.txt"
        if ifl==1:
            fout=dirout+rgns[i]+"_"+text+"__"+nddstr+"d3h_TOA_OBS.txt"
        f=open(fout,'w')
        toaobs=Dataset(fpath[ifl],'a')
        lon= toaobs.variables['lon'][:]   # read 
        tmstp=toaobs.variables['time'][:]
        lat=toaobs.variables['lat'][:]
        varnms=[]
        ondim=[]
        for nm in toaobs.variables:
            if nm=='ulwrf':        
                varnms.append(nm)            
                tmpvar=toaobs.variables[nm][:] ## time lat lon
                nt=len(tmpvar[:,0,0])
                nlat=len(tmpvar[0,:,0])
                nlon=len(tmpvar[0,0,:])
                for it in range(rnst,rnst+ndd*4+1): #### dataset is every 3 hours
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
        for ii in range(0,ndd*4+1):            
            for ij in range(0,lenv):                  
                itme="%f "%ondim[ii+ij*(ndd*4+1)]
                f.write(itme)
            f.write('\n')
        f.close()        