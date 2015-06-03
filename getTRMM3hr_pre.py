#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 14:43:13 2015

@author: jhchen
"""
from pytrmm import TRMM3B40RTFile as TRM3B40
import os
import datetime
os.system("cls")
dirin='E:/Data/TRMM/3b40/'
dirout='D:/MyPaper/PhD04/Data/TRMM/3B40/'
regions=['ETP','WTP','PRD','MLYR','NPC','NEC']
slon =[90.  , 80. , 107.5 , 110. , 110. , 120. ]
elon =[100. , 90. , 117.5 , 120. , 120. , 130. ]
slat =[27.5 , 27.5, 17.5  , 25.  , 32.5 , 40.  ]
elat =[37.5 , 37.5, 27.5  , 35.  , 42.5 , 50.  ]
# the start date of ever region and how many days for every region.
ayear=[2012 ,2010 , 2012  , 2010 , 2010 , 2012 ]
amon =[5    ,  7  ,  4    ,  6   ,  8   , 7    ] 
aday =[20   ,  14 ,  1    ,  5   ,  1   , 6    ]
ndays=[31   ,  31 ,  31   , 31   ,  31  , 31   ]
nag=len(regions)
#get the basic information of the data TRMM 3B40
fpath="E:/Data/TRMM/3b40/2010/06/3B40RT.2010062718.7R2.bin.gz"
trmm_file = TRM3B40(fpath)
fpath='E:/Data/TRMM/3b40/3b40_FileHead.txt'
fout=open(fpath,'w')
print(trmm_file.header())
ax=trmm_file.header()
headnm=ax.keys()
headva=ax.values()
na=len(headnm)
for a in range(0,na) :
    item="%s: "%headnm[a]
    fout.write(item)
    item="%s "%headva[a]
    fout.write(item)
    fout.write('\n')
precip = trmm_file.precip()
print('Array dimensions:', precip.shape)
fout.write('Array dimensions:')
nxy=[]
for a in precip.shape:
    nxy.append(a)
    item="%d "%a
    fout.write(item)
del trmm_file,precip
ny=nxy[0]
nx=nxy[1]
lon=[]
lat=[]
for ix in range(0,nx):
    lon.append(0.125+ix*0.25)
for iy in range(0,ny):
    lat.append(89.875-iy*0.25)
#
def getdatestr(yy,mm,dd,nd):
    datestart=datetime.datetime(yy,mm,dd,0,0,0)
    det=datetime.timedelta(hours=3)            
    dateiso=[]
    nt=nd*8          
    for dt in range(0,nt):
        dateiso.append(datestart+dt*det)
    xdate=[]               
    for tm in dateiso:
        xdate.append(datetime.datetime.strftime(tm,"%b/%d %H:%M"))
    return xdate
for ig in range(0,nag):
    yearstr="%d"%ayear[ig]
    monstr="%2.2d"%amon[ig]
    daystr="%2.2d"%aday[ig]
    ndaystr="%2.2d"%ndays[ig]
    fpath=dirout+regions[ig]+"-"+yearstr+monstr+daystr+"-"+ndaystr+"d_TRMM3B40.txt"
    fout=open(fpath,"w")
    fout.write("Date ")
    fout.write("Precipitation")
    fout.write("\n")
    ndd=ndays[ig]
    xdate=getdatestr(ayear[ig],amon[ig],aday[ig],ndays[ig])
    idt=0
    for idd in range(0,ndd):
        for ih in range(0,24,3):
            fout.write(xdate[idt]+" ")
            idt=idt+1
            hourstr="%2.2d"%ih
            filename="3B40RT."+yearstr+monstr+daystr+hourstr+".7R2.bin.gz"
            fold=yearstr+"/"+monstr+"/"
            fpath=dirin+fold+filename
            if os.path.isfile(fpath):
                trmm_file=TRM3B40(fpath)
                precip = trmm_file.precip()
                cont=0.
                tmp=0.
                for ix in range(0,nx):
                    if lon[ix]>(slon[ig]) and lon[ix]<(elon[ig]):
                        for iy in range(0,ny):
                            if lat[iy]>slat[ig] and lat[iy]<elat[ig]:
                                if precip[iy,ix]>=0.0 :
                                    tmp=tmp+precip[iy,ix]
                                    cont=cont+1.
                item="%f "%(tmp/cont)
                fout.write(item)
                del trmm_file,precip
            else:
                fout.write("-99.0")  # there is no data file
            fout.write("\n")
    fout.close()
            
