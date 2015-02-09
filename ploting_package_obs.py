#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 19:53:14 2015

@author: jhchen
"""
import matplotlib
import numpy as np
import matplotlib.cm as cm
import datetime
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.dates as matdate
from matplotlib.dates import DateFormatter
import calendar
import string
import time 
casenm='ETP2D2'
strnm='ETP06'
dirin='D:/MyPaper/PhD04/Cases/ETP/20100604_0704/'
iy=2010
im=6
jd=4
nt=121
nz=52
nday=30
pic_out='D:/MyPaper/PhD04/Pics/'
fnm1=['OBS_Surface_input.txt','.43','.49','.dyn',
     '_lsforcing.37','_surface.39','_thetaqv_profile.41','_uv_profiles.35',
     '.99']
fnm=[]
for strs in fnm1:
     fnm.append(strnm+strs)    
cnname="ETP20100604_030.txt"
nv=[3,3,3,5,2,5,2,2,2]    # number variables of every file
ndim=[1,52,1,52,52,1,52,52,52]  # is the varlables has the vertical dimension
lev99=[-9,-6,-3,3,6,9]
color99=['g','g','g','r','r','r']
##################################################
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['contour.negative_linestyle'] = 'dashed'
matplotlib.rcParams['ytick.labelsize'] = 20
matplotlib.rcParams['xtick.labelsize'] = 20
plt.rc('lines', linewidth=4)
###################################################
ydat_r=[ -50.000 ,    50.000 ,   164.286,    307.143,    478.571  ,  678.571 ,
      907.143 ,  1164.286,   1450.000,   1764.286 ,  2107.143,   2478.572 ,
      2878.572,   3307.143,  3764.286,  4250.000,   4764.286,   5307.143, 
      5878.571,   6478.571,   7107.143,  7764.286,  8450.000,  9164.285,  
      9907.143,  10678.570,  11478.570,  12307.143,  13164.285,  14050.000,
      14964.285,  15907.143,  16878.572,  17878.572,  18907.145,  19964.285,
      21050.000,  22164.285,  23307.145,  24478.572,  25678.572,  26907.145,
      28164.285,  29450.000,  30764.285,  32107.145,  33478.570,  34878.570,
      36307.141,  37764.285,  39250.000,  40750.000]
ydat=[]
for yd in ydat_r:
    ydat.append(yd*0.001)
datestart=datetime.datetime(iy,im,jd,0,0,0)
det=datetime.timedelta(hours=6)            
dateiso=[] 
xdate=[]                
for dt in range(0,nt):
    dateiso.append(datestart+dt*det)           
for tm in dateiso:
        xdate.append(datetime.datetime.strftime(tm,"%b-%d"))
xxx=range(0,nt)
###############################################################################
#### forcing
starid=0
fpath=dirin+fnm[4]
onedim1=[]
linesplit=[]
f=open(fpath)
ff=f.readlines()  ## first line in obs file is legend 
for line in ff:
    line=string.lstrip(line)
    linesplit.append(line[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            onedim1.append(string.atof(strs))
fcq1=np.ndarray(shape=(nz,nt), dtype=float)
fcq2=np.ndarray(shape=(nz,nt), dtype=float)
scalefc=24*3600.
cp=1005.
hlat=2.5e6
for it in range(0,nt):
    for iz in range(0,nz):        
        k=it*(nz*2+1)+iz+1  ## the first record is time
        fcq1[iz,it]=onedim1[k]*scalefc
        fcq2[iz,it]=onedim1[k+nz]*scalefc*hlat/cp    ## convert to kg/kg to K per day   
fcq1obs_pf=np.ndarray(shape=(nz), dtype=float)
fcq2obs_pf=np.ndarray(shape=(nz), dtype=float)
for iz in range(0,nz):
    tmp1=0.0
    tmp2=0.0
    for it in range(starid,nt-1):
        tmp1=tmp1+fcq1[iz,it]/(nt-starid-1.)
        tmp2=tmp2+fcq2[iz,it]/(nt-starid-1.)
    fcq1obs_pf[iz]=tmp1
    fcq2obs_pf[iz]=tmp2
###############################################################################
########## obs q1 q2
fpath=dirin+fnm[len(fnm)-1]
onedim1=[]
linesplit=[]
f=open(fpath)
ff=f.readlines()  ## 
for line in ff:
    line=string.lstrip(line)
    linesplit.append(line[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            onedim1.append(string.atof(strs))
q1obs=np.ndarray(shape=(nz,nt), dtype=float)
q2obs=np.ndarray(shape=(nz,nt), dtype=float)
for it in range(0,nt):
    for iz in range(0,nz):        
        k=it*(nz*2+1)+iz+1  ## the first record is time
        q1obs[iz,it]=onedim1[k]
        q2obs[iz,it]=onedim1[k+nz]       
q1obs_pf=np.ndarray(shape=(nz), dtype=float)
q2obs_pf=np.ndarray(shape=(nz), dtype=float)
for iz in range(0,nz):
    tmp1=0.0
    tmp2=0.0
    for it in range(starid,nt-1):
        tmp1=tmp1+q1obs[iz,it]/(nt-starid-1.)
        tmp2=tmp2+q2obs[iz,it]/(nt-starid-1.)
    q1obs_pf[iz]=tmp1
    q2obs_pf[iz]=tmp2
### Time series
lev37=[-4,-2,-1,2,6,9]
color37=['g','g','g','r','r','r']
linetyp37=['dotted','dotted','dotted','solid','solid','solid'] 
charsize=20
font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 20,
        }    
fig,[axe1,axe2]=plt.subplots(nrows=2,ncols=1,figsize=(15,9))
plt.subplot(2,1,1)
axe1=plt.contour(xxx,ydat,q1obs,colors=color37,
linewidths=1.5,levels=lev37,linestyles=linetyp37)                           
#plt.title('Observation',fontsize=charsize)                       
plt.axis([0, nt, 0, 16])  ## x axis  y axis
plt.clabel(axe1,inline=1,fmt='%1d',fontsize=charsize-2)                                                 
axx=fig.add_subplot(2,1,1)                         
axx.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
axx.set_xticklabels(xticklabels, size=charsize)
text1=r"($a$)"
axx.text(1.5,14,text1,fontsize=charsize+4) 
text1=r"Observation  $Q_1$ ($K$ $d^{-1}$)"
axx.text(80,14,text1,fontsize=charsize)  
plt.ylabel('Height'+r' ($km$)', fontdict=font)
plt.show() 
###
plt.subplot(2,1,2)
axe2=plt.contour(xxx,ydat,q2obs,colors=color37,
linewidths=1.5,levels=lev37,linestyles=linetyp37)                           
#plt.title(casenm,fontsize=charsize)                        
plt.axis([0, nt, 0, 16])  ## x axis  y axis
plt.clabel(axe2,inline=1,fmt='%1d',fontsize=charsize-2)                       
axx=fig.add_subplot(2,1,2)                         
axx.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
axx.set_xticklabels(xticklabels, size=charsize)
text1=r"($b$)"
axx.text(1.5,14,text1,fontsize=charsize+4) 
text1=r"Observation  $Q_2$ ($K$ $d^{-1}$)"
axx.text(80,14,text1,fontsize=charsize) 
plt.ylabel('Height'+r' ($km$)', fontdict=font)
plt.show()                     
plt.savefig(pic_out+casenm+'_q1q2_input.pdf')          
plt.show()
plt.close()
f.close()
###forcing
lev37=[-6,-3,-1,1,3,6]
color37=['g','g','g','r','r','r']
linetyp37=['dotted','dotted','dotted','solid','solid','solid'] 
charsize=20
font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 20,
        } 
fig,[axe1,axe2]=plt.subplots(nrows=2,ncols=1,figsize=(15,9))
plt.subplot(2,1,1)
axe1=plt.contour(xxx,ydat,fcq1,colors=color37,
linewidths=1.5,levels=lev37,linestyles=linetyp37)                           
#plt.title('Observation',fontsize=charsize)                       
plt.axis([0, nt, 0, 16])  ## x axis  y axis
plt.clabel(axe1,inline=1,fmt='%1d',fontsize=charsize-2)                                                 
axx=fig.add_subplot(2,1,1)                         
axx.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
axx.set_xticklabels(xticklabels, size=charsize)
text1=r"($a$)"
axx.text(1.5,14,text1,fontsize=charsize+4) 
text1=r"Temperature forcing ($K$ $d^{-1}$)"
axx.text(80,14,text1,fontsize=charsize)  
plt.ylabel('Height'+r' ($km$)', fontdict=font)
plt.show() 
###
plt.subplot(2,1,2)
axe2=plt.contour(xxx,ydat,fcq2,colors=color37,
linewidths=1.5,levels=lev37,linestyles=linetyp37)                           
#plt.title(casenm,fontsize=charsize)                        
plt.axis([0, nt, 0, 16])  ## x axis  y axis
plt.clabel(axe2,inline=1,fmt='%1d',fontsize=charsize-2)                       
axx=fig.add_subplot(2,1,2)                         
axx.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
axx.set_xticklabels(xticklabels, size=charsize)
text1=r"($b$)"
axx.text(1.5,14,text1,fontsize=charsize+4) 
text1=r"Moisture forcing ($g$ $kg^{-1}$ $d^{-1}$)"
axx.text(80,14,text1,fontsize=charsize) 
plt.ylabel('Height'+r' ($km$)', fontdict=font)
plt.show()                     
plt.savefig(pic_out+casenm+'_lsforcing_input.pdf')          
plt.show()
plt.close()
f.close()
###############################################################################
#'OBS_Surface_input.txt' 0
fpath=dirin+fnm[0]
onedim1=[]
linesplit=[]
f=open(fpath)
ff=f.readlines()[1:]  ## skip the first line 
for line in ff:
    line=string.lstrip(line)
    linesplit.append(line[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            onedim1.append(string.atof(strs))
q1inted=np.ndarray(shape=(nt), dtype=float)
q2inted=np.ndarray(shape=(nt), dtype=float)
nceprain=np.ndarray(shape=(nt), dtype=float)
for it in range(0,nt):      
    k=it*3  ## the first record is time
    q1inted[it]=onedim1[k]
    q2inted[it]=onedim1[k+1]
    nceprain[it]=onedim1[k+2]
#  CN05 daily rain 
fpath=dirin+cnname
onedim1=[]
linesplit=[]
f=open(fpath)
ff=f.readlines()  ## skip the first line 
for line in ff:
    line=string.lstrip(line)
    linesplit.append(line[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            onedim1.append(string.atof(strs))
cnrain=np.ndarray(shape=(nday), dtype=float)
cntmp=np.ndarray(shape=(nday), dtype=float)
for it in range(0,nday):      
    k=it*3  ## the first record is time
    cntmp[it]=onedim1[k+1]
    cnrain[it]=onedim1[k+2]
###############################################################################
            