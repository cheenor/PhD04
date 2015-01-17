#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 31 13:27:10 2014
This script handles the output data of CRM 
output: 96 figures of a day
every figure is the mean of 29 days. we don't wan
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
del ydat_r 
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
monstr="%02d"%(im)  ### number to string, 1 to 01, 10 to 10
dnm=calendar.monthrange(iy,im)[1]  #calendar.monthrange(1997,7) 
#                              #reture two index, sencond is the day number
datestart=datetime.datetime(iy,im,jd,0,0,0)
det=datetime.timedelta(minutes=15)            
dateiso=[] 
xdate=[]                
for dt in range(0,dnm*4*24):
    dateiso.append(datestart+dt*det)           
for tm in dateiso:
        xdate.append(datetime.datetime.strftime(tm,"%b-%d %H:%M")) 
# some parameters for plotting
varnm=["qc","qa","qb","qr","Total_Cloud_Water"]
levqc=[0.0,0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7]
levqa=[0.0,0.02,0.04,0.1,0.3,0.5,0.7,0.9,1.2,1.5,1.8,2.4]
levqb=[0.0,0.02,0.04,0.1,0.3,0.5,0.7,0.9,1.2,1.5,1.8,2.4]
levqr=[0.0,0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7]
levsum=[0.0,0.1,0.4,0.8,1.2,1.5,1.8,2.1,2.4,2.7,3.0]
color=['skyblue','1.0','0.95','0.9','0.85','0.8','0.75','0.7','0.65','0.6','0.5','0.4','0.3','0.2']
#####skyblue                
#        
filepath=dirin+filenm
for it in range(1394,1500): #nt): 1643
    print it
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
    del linesplit,line,lnstrs,strs
    del ff
    qc=np.ndarray(shape=(nz,nx), dtype=float)
    qa=np.ndarray(shape=(nz,nx), dtype=float)
    qb=np.ndarray(shape=(nz,nx), dtype=float)
    qr=np.ndarray(shape=(nz,nx), dtype=float)
    qds=np.ndarray(shape=(nz,nx), dtype=float)
    for ix in range(0,nx):
        for iz in range(0,nz):
            itxz=ix*nz*4+iz
            qc[iz,ix]=onedim1[itxz]
            itxz=ix*nz*4+iz+nz
            qa[iz,ix]=onedim1[itxz]
            itxz=ix*nz*4+iz+nz*2
            qb[iz,ix]=onedim1[itxz]
            itxz=+ix*nz*4+iz+nz*3
            qr[iz,ix]=onedim1[itxz]  
    del onedim1
    qc[0,:]=0.0 
    qa[0,:]=0.0 
    qb[0,:]=0.0 
    qc[1,:]=0.0 
    qa[1,:]=0.0 
    qb[1,:]=0.0 
    for i in range(0,5): 
        if i== 0 :
            qds[:,:]=qc[:,:]
            levs=levqc
        if i== 1 :
            qds[:,:]=qa[:,:]
            levs=levqa
        if i== 2 :
            qds[:,:]=qb[:,:]  
            levs=levqb
        if i== 3 :
            qds[:,:]=qr[:,:]
            levs=levqr
        if i== 4 :
            qds[:,:]=qc[:,:]+qa[:,:]+qb[:,:]+qr[:,:]
            levs=levsum
        qds[0,:]=0.0        
        if np.max(qds) == np.min(qds):
            qds[2,0]=0.005
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
#            axx.set_ylabel(ylabs,fontsize=28)
#            axx.set_xticks(range(0,nx,20))
#            xticklabels = [xlabs for nn in range(0,nx,20)] 
#            axx.set_xticklabels([xlabs[nn] for nn in range(0,nx,20)], size=16)    
#            plt.show()
            varnm[4]='Total_Cloud_Water'                    
            plt.savefig(dirout+varnm[i]+'_'+casenm+'_'+timstr+'.png')          
            plt.show()
            plt.close()
    del qa,qb,qc,qr,qds
    f.close()
#    print np.max(qds),np.min(qds)