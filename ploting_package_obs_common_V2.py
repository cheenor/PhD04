#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 21:35:25 2015

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
import os
casenm='ETPNEW'
strnm='ETP06'
dirin='D:/MyPaper/PhD04/Cases/ETP/20100604_0704NEW/'
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
def readAscii(fpath,iskp):
    #iskp  the total line skipped of the file
    # fpath   the full path of the file
    # usage: onedim=readAscii(fpaht,iskp)
    onedim=[]
    linesplit=[]
    f=open(fpath)
    ff=f.readlines()[iskp:]  ## first line in obs file is legend 
    for line in ff:
        line=string.lstrip(line)
        linesplit.append(line[:-1].split(' '))
    for lnstrs in linesplit:
        for strs in lnstrs:
            if strs!='':
                onedim.append(string.atof(strs))
    del linesplit
    f.close()
    return onedim
#
def plotcontour(x,y,z,*args):
    # x  xdata 1 dim
    # y  ydata 1 dim
    # z  zdata 2 dims (ny*nx)
    # args selected parameters the order:
    nx=len(x)
    ny=len(y)
    dm1=len(z[:,0])
    dm2=len(z[0,:])
    if nx != dm1 or ny != dm2 :
        print 'the dimension of input data is wrong'
        os.system("pause") 
    nr=1
    for aa in args:
        if len(np.shape(aa)) ==2 :
            if len(aa[:,0])==nx and len(aa[0,:])==ny :
                nr=nr+1
    fig,[axe1,axe2]=plt.subplots(nrows=nr,ncols=1,figsize=(nr*7,9))
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
#### forcing
starid=0
fpath=dirin+fnm[4]
iskp=0
onedim1=readAscii(fpath,iskp)
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
del onedim1
###############################################################################
########## obs q1 q2
fpath=dirin+fnm[len(fnm)-1]
iskp=0
onedim1=readAscii(fpath,iskp)
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
del onedim1
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
###############################################################################
#'OBS_Surface_input.txt' 0
fpath=dirin+fnm[0]
iskp=1
onedim1=readAscii(fpath,iskp)
q1inted=np.ndarray(shape=(nt), dtype=float)
q2inted=np.ndarray(shape=(nt), dtype=float)
nceprain=np.ndarray(shape=(nt), dtype=float)
for it in range(0,nt):      
    k=it*4  ## the first record is time
    q1inted[it]=onedim1[k+1]
    q2inted[it]=onedim1[k+2]
    nceprain[it]=onedim1[k+3]
fig,(ax0,ax1) = plt.subplots(nrows=2,ncols=1,figsize=(12,6))
ax0.plot(xxx,q1inted,'g',label='Q1')#[0:lcc-3])
plt.axis([0, nt, -35, 35])  ## x axis  y axis
ax0.plot(xxx,q2inted,'b',label='Q2')
plt.axis([0, nt, -35, 35])  ## x axis  y axis
text1=r"($a$) Q1 and Q2"
ax0.text(2,33,text1,fontsize=16)
ax0.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
ax0.set_xticklabels(xticklabels, size=charsize)
#
ax1.plot(xxx,nceprain,'g',label='Q1')#[0:lcc-3])
plt.axis([0, nt, 0, 1])  ## x axis  y axis
text1=r"($a$) NCEP Rain"
ax1.text(2,0.9,text1,fontsize=16)
ax1.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
ax1.set_xticklabels(xticklabels, size=charsize)
plt.show()                     
plt.savefig(pic_out+casenm+'_intedQ1Q2_Rain.pdf')          
plt.show()
plt.close()    
del onedim1
#  CN05 daily rain 
fpath=dirin+cnname
iskp=0
onedim1=readAscii(fpath,iskp)
cnrain=np.ndarray(shape=(nday), dtype=float)
cntmp=np.ndarray(shape=(nday), dtype=float)
for it in range(0,nday):      
    k=it*3  ## the first record is time
    cntmp[it]=onedim1[k+1]
    cnrain[it]=onedim1[k+2]
###############################################################################
#File 39
fpath=dirin+fnm[5]
iskp=0
onedim1=readAscii(fpath,iskp)
ths=np.ndarray(shape=(nt), dtype=float)
qvss=np.ndarray(shape=(nt), dtype=float)
sst=np.ndarray(shape=(nt), dtype=float)
flh=np.ndarray(shape=(nt), dtype=float)
fsh=np.ndarray(shape=(nt), dtype=float)
for it in range(0,nt):
    k=it*6
    ths[it]= onedim1[k+1]
    qvss[it]= onedim1[k+2]
    sst[it]= onedim1[k+3]
    flh[it]= onedim1[k+4]
    fsh[it]= onedim1[k+5]
fig,(ax0,ax1) = plt.subplots(nrows=2,ncols=1,figsize=(12,6))
ax0.plot(xxx,flh,'g',label='Latent Heat')
text1=r"($a$) Latent Heat"
ax0.set_title(text1, loc='left',fontsize=16)
ax0.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
ax0.set_xticklabels(xticklabels, size=charsize) 
ax1.plot(xxx,fsh,'b',label='Sensible Heat')
text1=r"($b$) Sensible Heat"
ax1.set_title(text1, loc='left',fontsize=16)
ax1.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
ax1.set_xticklabels(xticklabels, size=charsize)
plt.show()                     
plt.savefig(pic_out+casenm+'_SHLH_input.pdf')          
plt.show()
plt.close() 
del onedim1 
###############################################################################
#File 35
fpath=dirin+fnm[7]
iskp=0
onedim1=readAscii(fpath,iskp)
uwnd=np.ndarray(shape=(nz,nt), dtype=float)
vwnd=np.ndarray(shape=(nz,nt), dtype=float)
for it in range(0,nt):
    for iz in range(0,nz):
        k=it*(2*nz+1)
        vwnd[iz,it]=onedim1[k+iz+1]*10
        uwnd[iz,it]=onedim1[k+iz+1+nz]*10
#        
#lev37=[-4,-2,-1,2,6,9]
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
axe1=plt.contour(xxx,ydat,vwnd,colors=color37,
linewidths=1.5,linestyles=linetyp37)                           
#plt.title('Observation',fontsize=charsize)                       
plt.axis([0, nt, 0, 16])  ## x axis  y axis
plt.clabel(axe1,inline=1,fmt='%1d',fontsize=charsize-2)                                                 
axx=fig.add_subplot(2,1,1)                         
axx.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
axx.set_xticklabels(xticklabels, size=charsize)
text1=r"($a$) V wind"
axx.set_title(text1, loc='left',fontsize=16) 
plt.ylabel('Height'+r' ($km$)', fontdict=font)
plt.show() 
###
plt.subplot(2,1,2)
axe2=plt.contour(xxx,ydat,uwnd,colors=color37,
linewidths=1.5,linestyles=linetyp37)                           
#plt.title(casenm,fontsize=charsize)                        
plt.axis([0, nt, 0, 16])  ## x axis  y axis
plt.clabel(axe2,inline=1,fmt='%1d',fontsize=charsize-2)                       
axx=fig.add_subplot(2,1,2)                         
axx.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
axx.set_xticklabels(xticklabels, size=charsize)
text1=r"($b$) U wind"
axx.set_title(text1, loc='left',fontsize=16)
plt.ylabel('Height'+r' ($km$)', fontdict=font)
plt.show()                     
plt.savefig(pic_out+casenm+'_wind_input.pdf')          
plt.show()
plt.close()
del onedim1
#
fpath=dirin+fnm[6]
iskp=0
onedim1=readAscii(fpath,iskp)
theta=np.ndarray(shape=(nz,nt), dtype=float)
qv=np.ndarray(shape=(nz,nt), dtype=float)
for it in range(0,nt):
    for iz in range(0,nz):
        k=it*(2*nz+1)
        theta[iz,it]=onedim1[k+iz+1]
        qv[iz,it]=onedim1[k+iz+1+nz]*1000.
fig,[axe1,axe2]=plt.subplots(nrows=2,ncols=1,figsize=(15,9))
plt.subplot(2,1,1)
axe1=plt.contour(xxx,ydat,theta,colors=color37,
linewidths=1.5,linestyles=linetyp37)                           
#plt.title('Observation',fontsize=charsize)                       
plt.axis([0, nt, 0, 16])  ## x axis  y axis
plt.clabel(axe1,inline=1,fmt='%1d',fontsize=charsize-2)                                                 
axx=fig.add_subplot(2,1,1)                         
axx.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
axx.set_xticklabels(xticklabels, size=charsize)
text1=r"($a$) Theta"
axx.set_title(text1, loc='left',fontsize=16) 
plt.ylabel('Height'+r' ($km$)', fontdict=font)
plt.show() 
###
plt.subplot(2,1,2)
axe2=plt.contour(xxx,ydat,qv,colors=color37,
linewidths=1.5,linestyles=linetyp37)                           
#plt.title(casenm,fontsize=charsize)                        
plt.axis([0, nt, 0, 16])  ## x axis  y axis
plt.clabel(axe2,inline=1,fmt='%1d',fontsize=charsize-2)                       
axx=fig.add_subplot(2,1,2)                         
axx.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
axx.set_xticklabels(xticklabels, size=charsize)
text1=r"($b$) qv"
axx.set_title(text1, loc='left',fontsize=16)
plt.ylabel('Height'+r' ($km$)', fontdict=font)
plt.show()                     
plt.savefig(pic_out+casenm+'_thetaqv_input.pdf')          
plt.show()
plt.close()
del onedim1