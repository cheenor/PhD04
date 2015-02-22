#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 16 13:05:23 2015

@author: jhchen
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import string
import numpy as np
import datetime
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20

casenm='ETP2D3'
#casenm='WTP2D0'
#casenm='NPC2D2'
nt=121
nz=52
starid=4  # this parameter can discard the dirst day
if casenm[0:3]=='ETP':
    dirin='D:/MyPaper/PhD04/Cases/ETP/20100604_0704/Simulated/'
    dirobs='D:/MyPaper/PhD04/Cases/ETP/20100604_0704/'
    f43="ETP06.43"
    nameforcing="ETP06_lsforcing.37"
    iy=2010
    im=6
    jd=4
if casenm[0:3]=='WTP':    
    dirin='D:/MyPaper/PhD04/Cases/WTP/20100624_0723/Simulated/'
    dirobs='D:/MyPaper/PhD04/Cases/WTP/20100624_0723/'
    f43="WTP06.43"
    nameforcing="WTP06_lsforcing.37"
    iy=2010
    im=6
    jd=24
if casenm[0:3]=='NPC':
    dirin='D:/MyPaper/PhD04/Cases/NPC/20100802/Simulated/'
    dirobs='D:/MyPaper/PhD04/Cases/NPC/20100802/'
    f43="NPC.43"
    nameforcing="NPC_lsforcing.37"
    iy=2010
    im=8
    jd=2

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
smqv=np.ndarray(shape=(nz,nt),dtype=float)
smtc=np.ndarray(shape=(nz,nt),dtype=float)
smrat=np.ndarray(shape=(nz,nt),dtype=float)
iskp=12*nz
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
    smqv[0:nz,it]=onedim1[itts:itte]
    itts=itte
    itte=itts+nz
    smtc[0:nz,it]=onedim1[itts:itte]
    itts=itte
    itte=itts+nz
    smrat[0:nz,it]=onedim1[itts:itte]
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
obsrh=np.ndarray(shape=(nz,nt),dtype=float)
obsu=np.ndarray(shape=(nz,nt),dtype=float)
obsv=np.ndarray(shape=(nz,nt),dtype=float)
obsw=np.ndarray(shape=(nz,nt),dtype=float)           
iskp=7*nz+1
for it in range(0,nt):
    itts=it*iskp+1
    itte=itts+nz            
    obstha[0:nz,it]=onedim1[itts:itte]      
    itts=itte
    itte=itts+nz            
    obsqv[0:nz,it]=onedim1[itts:itte]    
    itts=itte
    itte=itts+nz            
    obstmp[0:nz,it]=onedim1[itts:itte]       
    itts=itte
    itte=itts+nz            
    obsrh[0:nz,it]=onedim1[itts:itte]       
    itts=itte
    itte=itts+nz            
    obsu[0:nz,it]=onedim1[itts:itte]       
    itts=itte
    itte=itts+nz            
    obsv[0:nz,it]=onedim1[itts:itte]       
    itts=itte
    itte=itts+nz            
    obsw[0:nz,it]=onedim1[itts:itte]       
obsqv=obsqv*1000. ## kg/kg to g/kg
###############################################################################
levs1=[-6,-3,4,8,12,15]
colors1=['g','g','r','r','r','r']
linetype1=['dotted','dotted','solid','solid','solid','solid'] 
levs2=[-3,-1,1,2,3,4]
colors2=['g','g','r','r','r','r']
linetype2=['dotted','dotted','solid','solid','solid','solid']
titlename=[r"Potential Temperature Bias ($K$)",r"Water Vapor Mixing Ratio Bias ($g$ $kg^{-1}$) "]
font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 16,
        }     
fig,(axe1,axe2)=plt.subplots(nrows=2,ncols=1,figsize=(15,8))
plt.subplot(2,1,1)
zdat=smt-obstha
zdat[0,:]=0.0   ## the first level is below surface ground
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
plt.savefig(dirpic+casenm+"_T&qv_Bias.pdf")          
plt.show()
plt.close()