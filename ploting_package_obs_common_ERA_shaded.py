#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu May 07 23:43:21 2015

@author: jhchen
"""
import matplotlib as mpl
#mpl.use("TkAgg")  #Qt4Agg  ,   ps ,  TkAgg
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
#plt.rc('text', usetex=True)
#mpl.rcParams['ps.usedistiller'] = None
#mpl.rcParams['interactive'] = True
#mpl.rcParams['ps.fonttype'] = 42
#mpl.rcParams('text.latex', unicode = True)
#mpl.rcParams['ps.useafm'] = True
casenm='ETPCTR_ERA'
if casenm[0:3]=='ETP' :
    iy,im,jd=2010,6,3
    nt,nday=125,31
    namestr=casenm[0:3]
    marktr=r"($f$)"
    yearstr="%d"%iy
elif casenm[0:3]=='WTP' :
    iy,im,jd=2010,7,3
    nt,nday=125,31
    namestr=casenm[0:3]
    marktr=r"($e$)"
    yearstr="%d"%iy
elif casenm[0:3]=='PRD' :
    iy,im,jd=2012,4,1
    nt,nday=125,31
    namestr=casenm[0:3]
    marktr=r"($a$)"
    yearstr="%d"%iy
elif casenm[0:3]=='MLY' :
    iy,im,jd=2010,6,2
    nt,nday=125,31
    namestr=casenm[0:4]
    marktr=r"($b$)"
    yearstr="%d"%iy
elif casenm[0:3]=='NPC' :
    iy,im,jd=2010,8,2
    nt,nday=125,31
    namestr=casenm[0:3]
    marktr=r"($c$)"
    yearstr="%d"%iy
elif casenm[0:3]=='NEC' :
    iy,im,jd=2012,7,6
    nt,nday=125,31
    namestr=casenm[0:3]
    marktr=r"($d$)"
    yearstr="%d"%iy
topstr='' #'_250'
datestr="%4d"%iy+"%2.2d"%im+"%2.2d"%jd+"_031d"
strnm=namestr+'_'+datestr
dirin='D:/MyPaper/PhD04/Cases/ERA/FORCING/'+namestr+'/'
dircnrain='D:/MyPaper/PhD04/Data/RainCN05/'
cnname=namestr+datestr+".txt"
nz=52
pic_out='D:/MyPaper/PhD04/Pics/'
fnm1=['_Q1Q2_ERA.38','_ERA.43','_ERA.49','_lsforcing_ERA.37',
     '_surface_ERA.39','_thetaqv_profile_ERA.41','_uv_profiles_ERA.35',
     '_ERA.99','_SHLH_ERA.43','_OmegaComps_ERA.42']
fnm=[]
for strs in fnm1:
    l=len(strs)
    fnm.append(strnm+strs[0:l-3]+topstr+strs[l-3:])
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
#
det=datetime.timedelta(hours=3)          
dateiso2=[] 
xdate2=[]                
for dt in range(0,2*nt-1):
    dateiso2.append(datestart+dt*det)           
for tm in dateiso2:
        xdate2.append(datetime.datetime.strftime(tm,"%b-%d"))
xxx2=range(0,2*nt-1)
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
#### forcing
starid=0
fpath=dirin+fnm[3]
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
fpath=dirin+fnm[7]
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
lev37=[-15,-9,-6,-3,3,6,9,15]
color37=['limegreen','palegreen','lightsage','w','w',
         'w','lightsalmon','salmon','tomato','r']
levq1=[-15,-9,-6,-3,3,6,9,15]
colorq1=['r','r','r','r','r','r','r','r']
linetypq1=['dotted','dotted','dotted','dotted','solid','solid','solid','solid'] 
levq2=[-12,-9,-6,-3,3,6,9,12]
colorq2=['limegreen','palegreen','lightsage','w','w',
         'w','lightsalmon','salmon','tomato','r']   #['g','g','g','g','g','g','g','g']
linetypq2=['dotted','dotted','dotted','dotted','solid','solid','solid','solid']
charsize=20
font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 20,
        }    
#
fig,[axe1,axe2]=plt.subplots(nrows=2,ncols=1,figsize=(15,9))
plt.subplot(2,1,1)
#axe1=plt.contourf(xxx,ydat,q1obs,colors=color37,levels=lev37, extent='both')
#linewidths=1.5,linestyles=linetyp37)
axe1=plt.contourf(xxx,ydat,q2obs,colors=colorq2,levels=levq2) #, 
#    linewidths=1.5,linestyles=linetypq2) 
axe1=plt.contour(xxx,ydat,q1obs,colors=colorq1,levels=levq1, 
    linewidths=1.5,linestyles=linetypq1)                        
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
axe2=plt.contourf(xxx,ydat,q2obs,colors=color37,levels=lev37)
#linewidths=1.5,linestyles=linetyp37)                           
#plt.title(casenm,fontsize=charsize)                        
plt.axis([0, nt, 0, 16])  ## x axis  y axis
#plt.clabel(axe2,inline=1,fmt='%1d',fontsize=charsize-2)                       
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
levq1=[-14,-10,-6,-3,3,6,10,14]
#colorq1=['r','r','r','r','r','r','r','r']
colorq1=['b','dodgerblue','deepskyblue','aqua','w',
         'lightsalmon','salmon','tomato','r']
linetypq1=['dotted','dotted','dotted','dotted','solid','solid','solid','solid'] 
levq2=[-12,-8,-5,-2,2,5,8,12]
colorq2=['darkgreen','darkgreen','darkgreen','darkgreen',
        'darkgreen','darkgreen','darkgreen','darkgreen']
#colorq2=['darkred','darkred','darkred','darkred',
#         'darkred','darkred','darkred','darkred']
linetypq2=['dotted','dotted','dotted','dotted','solid','solid','solid','solid']
charsize=20
font = {'family' : 'serif',
#         'serif' : 'Bookman',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 20,
        } 
fig,axe1=plt.subplots(nrows=1,ncols=1,figsize=(15,6))
plt.subplot(1,1,1)
#plt.rc('text', usetex=False)
#mpl.rcParams['ps.fonttype'] = 3
axe1=plt.contourf(xxx,ydat,fcq1,colors=colorq1,
                  levels=levq1,extend='both')
plt.colorbar(orientation='horizontal',extend='both',
    extendfrac='auto',  spacing='uniform') 
axe1=plt.contour(xxx,ydat,fcq2,colors=colorq2,
linewidths=1,levels=levq2,linestyles=linetypq2)                           
#plt.title('Observation',fontsize=charsize)                       
plt.axis([0, nt, 0, 16])  ## x axis  y axis
plt.clabel(axe1,inline=1,fmt='%1d',fontsize=charsize-2)                                                 
axx=fig.add_subplot(1,1,1)                         
axx.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
axx.set_xticklabels(xticklabels, size=charsize)
text1=marktr
axx.text(1.5,16.3,text1,fontsize=charsize+2) 
text1=namestr+r" Temperature ($K$ $d^{-1}$) and miosture ($g$ $kg^{-1}$ $d^{-1}$) forcing"
#axx.text(80,14,text1,fontsize=charsize) 
plt.title(text1,fontsize=charsize+4) 
plt.ylabel('Height'+r' ($km$)', fontsize=20) #fontdict=font)
plt.xlabel(yearstr, fontdict=font)
plt.show()                     
plt.savefig(pic_out+casenm+'_lsforcing_input.png',dpi=300)        
plt.show()
plt.close()
#filename=pic_out+casenm+'_lsforcing_input.eps'
#os.system("pdfcrop %s %s" % (filename, filename)) 
###############################################################################
#'OBS_Surface_input.txt' 0
fpath=dirin+fnm[0]
iskp=0
onedim1=readAscii(fpath,iskp)
q1inted=np.ndarray(shape=(nt), dtype=float)
q2inted=np.ndarray(shape=(nt), dtype=float)
for it in range(0,nt):      
    k=it*2  ## the first record is time
    q1inted[it]=onedim1[k]
    q2inted[it]=onedim1[k+1]
del onedim1
#
fpath=dirin+fnm[8]
iskp=0
onedim1=readAscii(fpath,iskp)
train=np.ndarray(shape=(nt*2-1), dtype=float)
crain=np.ndarray(shape=(nt*2-1), dtype=float)
fsh=np.ndarray(shape=(nt*2-1), dtype=float)
flh=np.ndarray(shape=(nt*2-1), dtype=float)
for it in range(0,nt*2-1):
    k=it*4
    flh[it]=onedim1[k]
    fsh[it]=onedim1[k+1]
    train[it]=onedim1[k+2]
    crain[it]=onedim1[k+3]
fig,(ax0,ax1) = plt.subplots(nrows=2,ncols=1,figsize=(15,6))
ax0=plt.subplot(2,1,1)
ax0.plot(xxx,q1inted,'g',label='Q1')#[0:lcc-3])
plt.axis([0, nt, -3, 5])  ## x axis  y axis
ax0.plot(xxx,q2inted,'b',label='Q2')
plt.axis([0, nt, -3, 5])  ## x axis  y axis
text1=r"($a$) Q1 and Q2"
ax0.text(2,5.5,text1,fontsize=16)
ax0.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
ax0.set_xticklabels(xticklabels, size=charsize)
#
ax1=plt.subplot(2,1,2)
ax1.plot(xxx2,train,'g',label='rain')#[0:lcc-3])
plt.axis([0, nt*2-1, 0, 1])  ## x axis  y axis
text1=r"($b$) ERA Rain"
ax1.text(2,0.9,text1,fontsize=16)
ax1.set_xticks(range(0,nt*2-1,32))
xticklabels = [xdate2[nn] for nn in range(0,nt*2-1,32)] 
ax1.set_xticklabels(xticklabels, size=charsize)
plt.show()                     
plt.savefig(pic_out+casenm+'_intedQ1Q2_Rain.pdf')          
plt.show()
plt.close()    
del onedim1
#  CN05 daily rain 
#fpath=dircnrain+cnname
#iskp=0
#onedim1=readAscii(fpath,iskp)
cnrain=np.ndarray(shape=(nday), dtype=float)
cntmp=np.ndarray(shape=(nday), dtype=float)
#for it in range(0,nday):      
#    k=it*3  ## the first record is time
#    cntmp[it]=onedim1[k+1]
#    cnrain[it]=onedim1[k+2]
###############################################################################
#File 39
fpath=dirin+fnm[4]
iskp=0
onedim1=readAscii(fpath,iskp)
ths=np.ndarray(shape=(nt), dtype=float)
qvss=np.ndarray(shape=(nt), dtype=float)
sst=np.ndarray(shape=(nt), dtype=float)
for it in range(0,nt):
    k=it*4
    ths[it]= onedim1[k+1]
    qvss[it]= onedim1[k+2]
    sst[it]= onedim1[k+3]
fig,(ax0,ax1) = plt.subplots(nrows=2,ncols=1,figsize=(12,6))
ax0.plot(xxx2,flh,'g',label='Latent Heat')
text1=r"($a$) Latent Heat"
ax0.set_title(text1, loc='left',fontsize=16)
ax0.set_xticks(range(0,nt*2-1,32))
xticklabels = [xdate2[nn] for nn in range(0,nt*2-1,32)] 
ax0.set_xticklabels(xticklabels, size=charsize) 
ax1.plot(xxx2,fsh,'b',label='Sensible Heat')
text1=r"($b$) Sensible Heat"
ax1.set_title(text1, loc='left',fontsize=16)
ax1.set_xticks(range(0,nt*2-1,32))
xticklabels = [xdate2[nn] for nn in range(0,nt*2-1,32)] 
ax1.set_xticklabels(xticklabels, size=charsize)
plt.show()                     
plt.savefig(pic_out+casenm+'_SHLH_input.pdf')          
plt.show()
plt.close() 
del onedim1 
###############################################################################
#File 35
fpath=dirin+fnm[6]
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
fpath=dirin+fnm[5]
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
# OMEGA
lev37=[-6,-4,-2,-1,0,1,2,4,6]
color37=['g','g','g','g','b','r','r','r','r']
linetyp37=['dotted','dotted','dotted','dotted','solid','solid','solid','solid','solid'] 
fpath=dirin+fnm[9]
iskp=0
onedim1=readAscii(fpath,iskp)
omegaa=np.ndarray(shape=(nz,nt), dtype=float)
omegao=np.ndarray(shape=(nz,nt), dtype=float)
q1a=np.ndarray(shape=(nz,nt), dtype=float)
q1v=np.ndarray(shape=(nz,nt), dtype=float)
q2a=np.ndarray(shape=(nz,nt), dtype=float)
q2v=np.ndarray(shape=(nz,nt), dtype=float)
for it in range(0,nt):
    for iz in range(0,nz):
        k=it*(6*nz+1)
        omegaa[iz,it]=onedim1[k+iz+1]*36
        omegao[iz,it]=onedim1[k+iz+1+nz]*36
        q1a[iz,it]=onedim1[k+iz+1+nz*2]
        q1v[iz,it]=onedim1[k+iz+1+nz*3]
        q2a[iz,it]=onedim1[k+iz+1+nz*4]
        q2v[iz,it]=onedim1[k+iz+1+nz*5]
fig,[axe1,axe2]=plt.subplots(nrows=2,ncols=1,figsize=(15,9))
plt.subplot(2,1,1)
axe1=plt.contour(xxx,ydat,omegaa,colors=color37,levels=lev37,
    linewidths=1.5,linestyles=linetyp37)                           
#plt.title('Observation',fontsize=charsize)                       
plt.axis([0, nt, 0, 16])  ## x axis  y axis
plt.clabel(axe1,inline=1,fmt='%1d',fontsize=charsize-2)                                                 
axx=fig.add_subplot(2,1,1)                         
axx.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
axx.set_xticklabels(xticklabels, size=charsize)
text1=r"($a$) Omega adjust"
axx.set_title(text1, loc='left',fontsize=16) 
plt.ylabel('Height'+r' ($km$)', fontdict=font)
plt.show() 
###
plt.subplot(2,1,2)
axe2=plt.contour(xxx,ydat,omegao,colors=color37,levels=lev37,
    linewidths=1.5,linestyles=linetyp37)                           
#plt.title(casenm,fontsize=charsize)                        
plt.axis([0, nt, 0, 16])  ## x axis  y axis
plt.clabel(axe2,inline=1,fmt='%1d',fontsize=charsize-2)                       
axx=fig.add_subplot(2,1,2)                         
axx.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
axx.set_xticklabels(xticklabels, size=charsize)
text1=r"($b$) Omega origin"
axx.set_title(text1, loc='left',fontsize=16)
plt.ylabel('Height'+r' ($km$)', fontdict=font)
plt.show()                     
plt.savefig(pic_out+casenm+'_omega_compare.pdf')          
plt.show()
plt.close()
del onedim1
#############################################################################
fig,[axe1,axe2]=plt.subplots(nrows=2,ncols=1,figsize=(15,9))
plt.subplot(2,1,1)
axe1=plt.contour(xxx,ydat,q1a,colors=color37,levels=lev37,
    linewidths=1.5,linestyles=linetyp37)                           
#plt.title('Observation',fontsize=charsize)                       
plt.axis([0, nt, 0, 16])  ## x axis  y axis
plt.clabel(axe1,inline=1,fmt='%1d',fontsize=charsize-2)                                                 
axx=fig.add_subplot(2,1,1)                         
axx.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
axx.set_xticklabels(xticklabels, size=charsize)
text1=r"($a$) Q1 advection"
axx.set_title(text1, loc='left',fontsize=16) 
plt.ylabel('Height'+r' ($km$)', fontdict=font)
plt.show() 
###
plt.subplot(2,1,2)
axe2=plt.contour(xxx,ydat,q1v,colors=color37,levels=lev37,
    linewidths=1.5,linestyles=linetyp37)                           
#plt.title(casenm,fontsize=charsize)                        
plt.axis([0, nt, 0, 16])  ## x axis  y axis
plt.clabel(axe2,inline=1,fmt='%1d',fontsize=charsize-2)                       
axx=fig.add_subplot(2,1,2)                         
axx.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
axx.set_xticklabels(xticklabels, size=charsize)
text1=r"($b$) Q1 vertical"
axx.set_title(text1, loc='left',fontsize=16)
plt.ylabel('Height'+r' ($km$)', fontdict=font)
plt.show()                     
plt.savefig(pic_out+casenm+'_q1_comps.pdf')          
plt.show()
plt.close()
#############################################################################
fig,[axe1,axe2]=plt.subplots(nrows=2,ncols=1,figsize=(15,9))
plt.subplot(2,1,1)
axe1=plt.contour(xxx,ydat,q2a,colors=color37,levels=lev37,
    linewidths=1.5,linestyles=linetyp37)                           
#plt.title('Observation',fontsize=charsize)                       
plt.axis([0, nt, 0, 16])  ## x axis  y axis
plt.clabel(axe1,inline=1,fmt='%1d',fontsize=charsize-2)                                                 
axx=fig.add_subplot(2,1,1)                         
axx.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
axx.set_xticklabels(xticklabels, size=charsize)
text1=r"($a$) Q2 advection"
axx.set_title(text1, loc='left',fontsize=16) 
plt.ylabel('Height'+r' ($km$)', fontdict=font)
plt.show() 
###
plt.subplot(2,1,2)
axe2=plt.contour(xxx,ydat,q2v,colors=color37,levels=lev37,
    linewidths=1.5,linestyles=linetyp37)                           
#plt.title(casenm,fontsize=charsize)                        
plt.axis([0, nt, 0, 16])  ## x axis  y axis
plt.clabel(axe2,inline=1,fmt='%1d',fontsize=charsize-2)                       
axx=fig.add_subplot(2,1,2)                         
axx.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
axx.set_xticklabels(xticklabels, size=charsize)
text1=r"($b$) Q2 vertical"
axx.set_title(text1, loc='left',fontsize=16)
plt.ylabel('Height'+r' ($km$)', fontdict=font)
plt.show()                     
plt.savefig(pic_out+casenm+'_q2_comps.pdf')          
plt.show()
plt.close()