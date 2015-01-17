#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 31 13:27:10 2014

@author: jhchen
"""
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
import calendar
import string
import numpy as np
nt=2880
nx=200
nz=52
casenm='ETP2D1'
iy=2010
im=6
jd=4
dirin='D:/MyPaper/PhD04/Cases/ETP/20100604_0704/Simulated/'
dirout='D:/MyPaper/PhD04/Cases/ETP/20100604_0704/Simulated/animate/'
filenm=casenm+'_Raw_qcqaqbqr.txt'
nlines=14976000  # the total lines of the file
reslu=3  ### reslution is 3km
####
mpl.rcParams['ytick.labelsize'] = 24
mpl.rcParams['xtick.labelsize'] = 24
mpl.rcParams[  'savefig.dpi'  ] = 150
###
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
xdat=[]
xdsp=[]
xlabs=[]
for i in range(0,nx):
    xdsp.append(i*reslu)
    xdat.append(i*reslu)
    xstr="%d"%(i*reslu)
    xlabs.append(xstr)
xlabs[nx-1]=xlabs[nx-1]+r' $Km$'
xlabs[nx-2]=xlabs[nx-2]+r' $Km$'  
xlabs[nx-3]=xlabs[nx-3]+r' $Km$'                  
#  
daymean=np.ndarray(shape=(4,96,nz,nx), dtype=float)      
filepath=dirin+filenm
for itd in range(96,nt,96): #nt):96 skip the first day, 96 every day
    for itm in range(0,96):
        it=itd+itm
        timstr="%04d"%(it) 
        nls=it*5200
        nle=nls+5200  #### there is 5200 lines in every time,52*200*4/8 , this depens on the output format of fortran in file mean_tquwql.f
        onedim1=[]
        linesplit=[]
        f=open(filepath)
        ff=f.readlines()[nls:nle] ## we just read the data of one timesetp erevy time to save the memory
        for line in ff:
            line=string.lstrip(line)
            linesplit.append(line[:-1].split(' '))
        for lnstrs in linesplit:
            for strs in lnstrs:
                if strs!='':
                    onedim1.append(string.atof(strs))
        nl=len(onedim1)
        del linesplit
        del ff
#    qc=np.ndarray(shape=(nz,nx), dtype=float)
#    qa=np.ndarray(shape=(nz,nx), dtype=float)
#    qb=np.ndarray(shape=(nz,nx), dtype=float)
#    qr=np.ndarray(shape=(nz,nx), dtype=float)
#    qds=np.ndarray(shape=(nz,nx), dtype=float)
        for ix in range(0,nx):
            for iz in range(0,nz):
                itxz=ix*nz*4+iz
                daymean[0,itm,iz,ix]=onedim1[itxz]/29.0+daymean[0,itm,iz,ix]
#            qc[iz,ix]=onedim1[itxz]
                itxz=ix*nz*4+iz+nz
                daymean[1,itm,iz,ix]=onedim1[itxz]/29.0+daymean[1,itm,iz,ix]
#            qa[iz,ix]=onedim1[itxz]
                itxz=ix*nz*4+iz+nz*2
                daymean[2,itm,iz,ix]=onedim1[itxz]/29.0+daymean[2,itm,iz,ix]
#            qb[iz,ix]=onedim1[itxz]
                itxz=+ix*nz*4+iz+nz*3
                daymean[3,itm,iz,ix]=onedim1[itxz]/29.0+daymean[3,itm,iz,ix]
#            qr[iz,ix]=onedim1[itxz]  
        del onedim1
        f.close()        
varnm=["qc","qa","qb","qr","Total_Cloud_Water"]
filepath=dirin+'qabcr_daily.txt'
f=open(filepath,'w')  ## daymean=np.ndarray(shape=(4,96,nz,nx), dtype=float)
for iv in range(0,4):                  
    itme="%s "%varnm[iv]
    f.write(itme)
f.write('\n')         
for it in range(0,96):
    for iz in range(0,nz):
        for ix in range(0,nx):
            for iv in range(0,4):                  
                itme="%f "%daymean[iv,it,iz,ix]
                f.write(itme)
            f.write('\n')
f.close()        
        
############ end reading file  ##################
monstr="%02d"%(im)  ### number to string, 1 to 01, 10 to 10
dnm=calendar.monthrange(iy,im)[1]  #calendar.monthrange(1997,7) 
#                              #reture two index, sencond is the day number
datestart=datetime.datetime(iy,im,jd,0,0,0)
det=datetime.timedelta(minutes=15)            
dateiso=[] 
xdate=[]                
for dt in range(0,96):
    dateiso.append(datestart+dt*det)           
for tm in dateiso:
        xdate.append(datetime.datetime.strftime(tm,"%H:%M")) 
# some parameters for plotting
varnm=["qc","qa","qb","qr","Total_Cloud_Water"]
levqc=[0.0,0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7]
levqa=[0.0,0.02,0.04,0.1,0.3,0.5,0.7,0.9,1.2,1.5,1.8,2.4]
levqb=[0.0,0.02,0.04,0.1,0.3,0.5,0.7,0.9,1.2,1.5,1.8,2.4]
levqr=[0.0,0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7]
levsum=[0.0,0.1,0.4,0.8,1.2,1.5,1.8,2.1,2.4,2.7,3.0]
colorr=['skyblue','1.0','0.95','0.9','0.85','0.8','0.75','0.7','0.65','0.6','0.5','0.4','0.3','0.2']
skycolor=['royalblue','royalblue','royalblue','royalblue',# 00
    'royalblue','royalblue','royalblue','royalblue',# 01
    'deepskyblue','deepskyblue','deepskyblue','deepskyblue',# 02
    'deepskyblue','deepskyblue','deepskyblue','deepskyblue',# 03
    'skyblue','skyblue','skyblue','skyblue',# 04
    'deepskyblue','deepskyblue','deepskyblue','deepskyblue',# 05
    'deepskyblue','deepskyblue','deepskyblue','deepskyblue',# 06
    'royalblue','royalblue','royalblue','royalblue',# 07
    'royalblue','royalblue','royalblue','royalblue',# 08
    'b','b','b','b',# 09
    'b','b','b','b',# 10
    'darkblue','darkblue','darkblue','darkblue',# 11
    'darkblue','darkblue','darkblue','darkblue',# 12
    'navy','navy','navy','navy',# 13
    'navy','navy','navy','navy',# 14
    'midnightblue','midnightblue','midnightblue','midnightblue',# 15
    'midnightblue','midnightblue','midnightblue','midnightblue',# 16
    'midnightblue','midnightblue','midnightblue','midnightblue',# 17
    'navy','navy','navy','navy',# 18
    'navy','navy','navy','navy',# 19
    'darkblue','darkblue','darkblue','darkblue',# 20
    'darkblue','darkblue','darkblue','darkblue',# 21
    'b','b','b','b',# 22
    'b','b','b','b' # 23
    ]
    
#####skyblue 
for it in range(0,1): 
    timstr="%02d"%(it)
    color=colorr
    color[0]=skycolor[it]
    for i in range(0,5):
        qds=np.ndarray(shape=(nz,nx), dtype=float)
        if i== 0 :
            qds[:,:]=daymean[i,it,:,:]
            levs=levqc
        if i== 1 :
            qds[:,:]=daymean[i,it,:,:]
            levs=levqa
        if i== 2 :
            qds[:,:]=daymean[i,it,:,:]  
            levs=levqb
        if i== 3 :
            qds[:,:]=daymean[i,it,:,:]
            levs=levqr
        if i== 4 :
            qds[:,:]=daymean[0,it,:,:]+daymean[1,it,:,:]+daymean[2,it,:,:]+daymean[3,it,:,:]
            levs=levsum
        qds[0,:]=0.0 
        if np.max(qds) == np.min(qds):
            qds[1,0]=0.005
        if np.max(qds) != np.min(qds):
            fig,ax0 = plt.subplots(nrows=1,ncols=1,figsize=(18,6))
            ax0=plt.contourf(xdsp,ydat[0:32],qds[0:32,:],colors=color,
                          levels=levs,extend='both')
            varnm[4]='Total Cloud Water'
            plt.text(600.,-1.,r'$km$',fontsize=28)
            titlstr=varnm[i]+"   UTC "+xdate[it] 
            plt.title(titlstr,fontsize=28) 
            ylabs='Height'+r' ($km$)'
            plt.ylabel(ylabs,fontsize=28)
            plt.xticks(range(0,nx*3,50))
#            l,b,w,h = plt.gca().get_position().bounds
#            cbar=plt.colorbar(ax0,shrink=0.7)
#            ll,bb,ww,hh = cbar.ax.get_position().bounds
#            cbar.ax.set_position([ll, b+0.1*h, ww, h*0.8])
#            plt.setp(cbar.ax.get_yticklabels(), visible=False)
#            axx=fig.add_subplot(111)
#            axx.set_ylabel(ylabs,fontsize=28)
#            axx.set_xticks(range(0,nx,20))
#            xticklabels = [xlabs for nn in range(0,nx,20)] 
#            axx.set_xticklabels(xticklabels, size=16)    
#            plt.show()
            varnm[4]='Total_Cloud_Water'                    
            plt.savefig(dirout+varnm[i]+'_'+casenm+'_fordaily_'+timstr+'.png')          
            plt.show()
            plt.close()
        del qds
#    print np.max(qds),np.min(qds)