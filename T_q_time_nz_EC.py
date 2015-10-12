#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri May 08 11:16:43 2015

@author: jhchen
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import string
import numpy as np
import datetime
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20

casenm='WTPCTR_EC'
#casenm='WTP2D0'
#casenm='NPC2D2'
nt=121
nz=52
starid=4  # this parameter can discard the first day
dirs='D:/MyPaper/PhD04/Cases/'
diro='D:/MyPaper/PhD04/Cases/ERA/FORCING/'
if casenm[0:3]=='ETP':
    area=casenm[0:3]
#    folds='/20100604_0704/Simulated/'
#    datestr='20120520_031d'    
    iy,im,jd=2012,5,20
    iy,im,jd=2010,6,3
if casenm[0:3]=='WTP':
    area=casenm[0:3]   
    iy,im,jd=2010,7,14  
#    iy,im,jd=2010,7,3  
if casenm[0:3]=='NPC':
    area=casenm[0:3]
#    folds='/20100802/Simulated/'
#    datestr='20100802_031d'    
    iy,im,jd=2010,8,2
if casenm[0:3]=='PRD':
    area=casenm[0:3]
#    folds='/20100402/Simulated/'
#    datestr='20120401_031d'    
    iy,im,jd=2012,4,1 
if casenm[0:3]=='MLY':
    area=casenm[0:4]
    folds='/20100624/Simulated/'
#    datestr='20100605_031d' 
#    datestr='20100624_031d'
#    datestr='20100602_031d'
    iy,im,jd=2010,6,2
#    iy,im,jd=2010,6,24
if casenm[0:3]=='NEC':
    area=casenm[0:3]
#    folds='/20100801/Simulated/'
#    datestr='20120706_031d'    
    iy,im,jd=2012,7,6 
folds="/CTREC"+"%4d"%iy+"%2.2d"%im+"%2.2d"%jd+"/Simulation/"
datestr="%4d"%iy+"%2.2d"%im+"%2.2d"%jd+"_031d"
dirin=dirs+area+folds
dirobs=diro+area+'/'
f43=area+'_'+datestr+"_ERA.43"
nameforcing=area+'_'+datestr+"_LSFORCING_ERA.37"
dirpic='D:/MyPaper/PhD04/Pics/'
#fpath=dirin+'micro_202_ETP2D3'
#fff=np.fromfile(fpath,dtype=float)
#fc=open(fpath,"rb")
#fff=fc.read()
#print fff[0],fff[1],fff[12],fff[78]
#print len(fff)
datestart=datetime.datetime(iy,im,jd,0,0,0)
det=datetime.timedelta(hours=6)            
dateiso=[]            
for dt in range(0,nt):
    dateiso.append(datestart+dt*det)
xdate=[]    
xdat=range(0,121)            
for tm in dateiso:
        xdate.append(datetime.datetime.strftime(tm,"%b/%d")) 
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
ylevs=[1000, 925, 850, 700, 600, 500, 400, 300, 250, 
       200, 150, 100, 70, 50, 30, 20, 10]
for yd in ydat_r:
    ydat.append(yd*0.001)
del ydat_r
###############################################################################
fpath=dirin+casenm+'_All.txt'
onedim1=[]
linesplit=[]
f=open(fpath)
ff=f.readlines()  ## 
for line in ff:
    lines=string.lstrip(line)
    linesplit.append(lines[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            onedim1.append(string.atof(strs))
smt=np.ndarray(shape=(nz,nt),dtype=float)
smq=np.ndarray(shape=(nz,nt),dtype=float)
smu=np.ndarray(shape=(nz,nt),dtype=float)
smw=np.ndarray(shape=(nz,nt),dtype=float)
smqc=np.ndarray(shape=(nz,nt),dtype=float)
smqr=np.ndarray(shape=(nz,nt),dtype=float)
smqa=np.ndarray(shape=(nz,nt),dtype=float)
smqb=np.ndarray(shape=(nz,nt),dtype=float)
smrh=np.ndarray(shape=(nz,nt),dtype=float)
smfqv=np.ndarray(shape=(nz,nt),dtype=float)
smftc=np.ndarray(shape=(nz,nt),dtype=float)
smrat=np.ndarray(shape=(nz,nt),dtype=float)
smtc=np.ndarray(shape=(nz,nt),dtype=float)
smtco=np.ndarray(shape=(nz,nt),dtype=float)
iskp=13*nz
for it in range(0,nt):
    itts=it*iskp
    itte=itts+nz
    smt[0:nz,it]=onedim1[itts:itte]
    itts=itte
    itte=itts+nz
    smq[0:nz,it]=onedim1[itts:itte]
    itts=itte
    itte=itts+nz
    smu[0:nz,it]=onedim1[itts:itte]
    itts=itte
    itte=itts+nz
    smw[0:nz,it]=onedim1[itts:itte]
    itts=itte
    itte=itts+nz
    smqc[0:nz,it]=onedim1[itts:itte]
    itts=itte
    itte=itts+nz
    smqr[0:nz,it]=onedim1[itts:itte]
    itts=itte
    itte=itts+nz
    smqa[0:nz,it]=onedim1[itts:itte]
    itts=itte
    itte=itts+nz
    smqb[0:nz,it]=onedim1[itts:itte]
    itts=itte
    itte=itts+nz
    smrh[0:nz,it]=onedim1[itts:itte]
    itts=itte
    itte=itts+nz
    smfqv[0:nz,it]=onedim1[itts:itte]
    itts=itte
    itte=itts+nz
    smftc[0:nz,it]=onedim1[itts:itte]
    itts=itte
    itte=itts+nz
    smrat[0:nz,it]=onedim1[itts:itte]
    itts=itte
    itte=itts+nz
    smtc[0:nz,it]=onedim1[itts:itte]
###############################################################################
##  open obs  files
del onedim1,linesplit
fpath=dirobs+f43
onedim1=[]
linesplit=[]
f=open(fpath)
ff=f.readlines()  ## 
for line in ff:
    lines=string.lstrip(line)
    linesplit.append(lines[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            onedim1.append(string.atof(strs))
obstha=np.ndarray(shape=(nz,nt),dtype=float)
obsqv=np.ndarray(shape=(nz,nt),dtype=float)
obstmp=np.ndarray(shape=(nz,nt),dtype=float)
#obsrh=np.ndarray(shape=(nz,nt),dtype=float)
obsu=np.ndarray(shape=(nz,nt),dtype=float)
obsv=np.ndarray(shape=(nz,nt),dtype=float)
obsw=np.ndarray(shape=(nz,nt),dtype=float)           
iskp=6*nz+1
for it in range(0,nt):
    for iz in range(0,nz):
        k=iskp*it+1+iz            
        obstha[iz,it]=onedim1[k]
        k=iskp*it+1+iz+nz*1            
        obsqv[iz,it]=onedim1[k]*1000. #convert kg/kg to g/kg  
        k=iskp*it+1+iz+nz*2         
        obstmp[iz,it]=onedim1[k]  
        k=iskp*it+1+iz+nz*3          
        obsu[iz,it]=onedim1[k]  
        k=iskp*it+1+iz+nz*4           
        obsv[iz,it]=onedim1[k]  
        k=iskp*it+1+iz+nz*5          
        obsw[iz,it]=onedim1[k]
        smtco[iz,it]=(smt[iz,it]-obstha[iz,it])/(obstha[iz,it]/obstmp[iz,it])
###############################################################################
levs1=[-15,-12,-6,-3,4,8,12,15]
colors1=['g','g','g','g','r','r','r','r']
linetype1=['dotted','dotted','dotted','dotted','solid','solid','solid','solid'] 
levs2=[-3,-1,1,2,3,4]
colors2=['g','g','r','r','r','r']
linetype2=['dotted','dotted','solid','solid','solid','solid']
titlename=[r"Temperature Bias ($K$)",r"Water Vapor Mixing Ratio Bias ($g$ $kg^{-1}$) "]
font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 16,
        }     
fig,(axe1,axe2)=plt.subplots(nrows=2,ncols=1,figsize=(15,8))
plt.subplot(2,1,1)
zdat=smtc+273.16-obstmp
#zdat=smtco   #-obstmp
##zdat[0,:]=0.0   ## the first level is below surface ground
axe1=plt.contour(xdat,ydat,zdat,colors=colors1,
    linewidths=1.5,levels=levs1,linestyles=linetype1)                           
plt.title(titlename[0],fontsize=16)                          
plt.axis([0, 121, 0, 16])
plt.clabel(axe1,inline=1,fmt='%1.0f',fontsize=12)
axx=fig.add_subplot(2,1,1) 
text1=r"($a$)"
axx.text(2,16.5,text1,fontsize=18)                        
axx.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
axx.set_xticklabels(xticklabels, size=16)
plt.ylabel(r'Height ($km$)', fontdict=font)
plt.show()                     
#
plt.subplot(2,1,2)
zdat=smq-obsqv
zdat[0,:]=0.0   ## the first level is below surface ground
axe2=plt.contour(xdat,ydat,zdat,colors=colors2,
    linewidths=1.5,levels=levs2,linestyles=linetype2)                           
plt.title(titlename[1],fontsize=16)                          
plt.axis([0, 121, 0, 16])
plt.clabel(axe2,inline=1,fmt='%1.0f',fontsize=12) 
axx=fig.add_subplot(2,1,2)   
text1=r"($b$)"
axx.text(2,16.5,text1,fontsize=18)                       
axx.set_xticks(range(0,nt,16))
xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
axx.set_xticklabels(xticklabels, size=16)
plt.ylabel(r'Height ($km$)', fontdict=font)
plt.show()
plt.savefig(dirpic+casenm+datestr+"_T&qv_Bias2.pdf")          
plt.show()
plt.close()
#--------------PROFILES -------------------------------------------------------
#obstmp[iz,it]=onedim1[k] obsqv[iz,it]
#smtco[iz,it]   smq
meanobs=np.ndarray(shape=(nz,2),dtype=float)
obsstd=np.ndarray(shape=(nz,2),dtype=float)
meansim=np.ndarray(shape=(nz,2),dtype=float)
simstd=np.ndarray(shape=(nz,2),dtype=float)
for k in range(0,nz):
    tmp1=0.
    tmp2=0.
    tmp3=0.
    tmp4=0.
    tmp=np.ndarray(shape=(nt,4),dtype=float)
    for it in range(0,nt):
        tmp1=tmp1+obstha[k,it]
        tmp2=tmp2+obsqv[k,it]
        tmp3=tmp3+smt[k,it]
        tmp4=tmp4+smq[k,it]
        tmp[it,0]=obstha[k,it]
        tmp[it,1]=obsqv[k,it]
        tmp[it,2]=smt[k,it]
        tmp[it,3]=smq[k,it]
    meanobs[k,0]=tmp1/nt
    meanobs[k,1]=tmp2/nt
    meansim[k,0]=tmp3/nt
    meansim[k,1]=tmp4/nt
    obsstd[k,0]=tmp[:,0].std()
    obsstd[k,1]=tmp[:,1].std()
    simstd[k,0]=tmp[:,2].std()
    simstd[k,1]=tmp[:,3].std()
#
colors=['g','r']
width=[1.5,1.5]
fig,(ax0,ax1)=plt.subplots(nrows=1,ncols=2,figsize=(10,12))
ax0=plt.subplot(1,2,1)
ax0.errorbar(meanobs[1:32,0],ydat[1:32],label="OBS",
    c=colors[0],lw=width[0],xerr=obsstd[1:32,0]) #qc
ax0.errorbar(meansim[1:32,0],ydat[1:32],label="SIM",
    c=colors[1],lw=width[1],xerr=simstd[1:32,0]) #qc
plt.ylabel(r'Height ($km$)', fontdict=font)
#axx=fig.add_subplot(1,2,1) 
#text1=r"($a$)"
#axx.text(285,16.5,text1,fontsize=14)
plt.xlim(300,500)                        
ax0.set_xticks(range(300,500,50))
#ax0.yaxis.limit_range_for_scale(0,16)
ax0.set_title(r'$(a)$  Potential Tempe. ($K$)')
ax1=plt.subplot(1,2,2)
ax1.errorbar(meanobs[1:32,1],ydat[1:32],label="OBS",
    c=colors[0],lw=width[0],xerr=obsstd[1:32,1]) #qc
ax1.errorbar(meansim[1:32,1],ydat[1:32],label="SIM",
    c=colors[1],lw=width[1],xerr=simstd[1:32,1]) #qc
#plt.ylabel(r'Height ($km$)', fontdict=font)
ax1.set_title(r'($b$) Vapor Mixing Ratio ($g$ $kg^{-1}$)')
ax1.set_xticks(range(0,25,4))
#ax1.yaxis.limit_range_for_scale(0,16)
#axx=fig.add_subplot(1,2,2) 
#text1=r"($b$)"
#axx.text(1,16.5,text1,fontsize=14)                        
#axx.set_xticks(range(0,35,5))
#axx.set_yticks(range(0,16,3))
plt.show()
plt.savefig(dirpic+casenm+datestr+"_obsVSsim_t_q2.pdf")          
plt.show()
plt.close()
#-------------------profiles of basis
bias=np.ndarray(shape=(nz,2),dtype=float)
biasstd=np.ndarray(shape=(it,2),dtype=float)
for k in range(0,nz):
    tmp1=0.
    tmp2=0.
    tmp=np.ndarray(shape=(nt,2),dtype=float)
    for it in range(0,nt):
        tmp1=tmp1+smtc[k,it]+273.16-obstmp[k,it]
        tmp2=tmp2+smq[k,it]-obsqv[k,it]
        tmp[it,0]=smtc[k,it]+273.16-obstmp[k,it]
        tmp[it,1]=smq[k,it]-obsqv[k,it]
    bias[k,0]=tmp1/nt
    bias[k,1]=tmp2/nt
    biasstd[k,0]=tmp[:,0].std()
    biasstd[k,1]=tmp[:,1].std()
#
colors=['k','lightgrey']
colors=['g','r']
width=[2,2]
fig,ax0=plt.subplots(nrows=1,ncols=1,figsize=(6,10))
#ax0=plt.subplot(1,2,1)
ax0.errorbar(bias[1:32,0],ydat[1:32],label="Temperature",
    c=colors[1],lw=width[0],xerr=biasstd[1:32,0]) #qc
ax0.errorbar(bias[1:32,1],ydat[1:32],label="Vapor",
    c=colors[0],lw=width[0],xerr=biasstd[1:32,1]) #qc
ax0.plot([0,0],[0,16],c='k',lw=1.5)
#ax0.errorbar(meansim[1:32,0],ydat[1:32],label="SIM",
#    c=colors[1],lw=width[1],xerr=simstd[1:32,0]) #qc
plt.ylabel(r'Height ($km$)', fontdict=font)
#axx=fig.add_subplot(1,2,1) 
#text1=r"($a$)"
#axx.text(285,16.5,text1,fontsize=14)
#plt.xlim(-4,4)                        
#ax0.set_xticks(range(4,4,2))
#ax0.yaxis.limit_range_for_scale(0,16)
ax0.set_title(area)
#ax1=plt.subplot(1,2,2)

#ax1.errorbar(meansim[1:32,1],ydat[1:32],label="SIM",
#    c=colors[1],lw=width[1],xerr=simstd[1:32,1]) #qc
#plt.ylabel(r'Height ($km$)', fontdict=font)
#ax1.set_title(r'($b$) Vapor Mixing Ratio ($g$ $kg^{-1}$)')
#ax1.set_xticks(range(0,25,4))
#ax1.yaxis.limit_range_for_scale(0,16)
#axx=fig.add_subplot(1,2,2) 
#text1=r"($b$)"
#axx.text(1,16.5,text1,fontsize=14)                        
#axx.set_xticks(range(0,35,5))
#axx.set_yticks(range(0,16,3))
plt.show()
plt.savefig(dirpic+casenm+datestr+"_BiasProfile_t_q2.pdf")          
plt.show()
plt.close()