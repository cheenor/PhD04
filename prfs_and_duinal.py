#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 08 10:12:57 2015
read and plot data that are generated by getdataplot.f90
@author: jhchen
"""
import matplotlib as mpl
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
#
def readAscii(fpath,iskp,nrl):
    #iskp  the total line skipped of the file
    # fpath   the full path of the file
    # usage: onedim=readAscii(fpaht,iskp)
    onedim=[]
    linesplit=[]
    f=open(fpath)
    print iskp,nrl
    ff=f.readlines()[iskp:nrl]  ## first line in obs file is legend 
    for line in ff:
        line=string.lstrip(line)
        linesplit.append(line[:-1].split(' '))
    for lnstrs in linesplit:
        for strs in lnstrs:
            if strs!='':
                onedim.append(string.atof(strs))
    del linesplit,ff
    f.close()
    print len(onedim)
    return onedim
dt=15
ndt=24*60/dt
nz=34
nbin=99
dirpic="D:/MyPaper/PhD04/Pics/"
casenm=['PRDCTR_EC','MLYRCTR_EC', 'NPCCTR_EC',
           'NECCTR_EC','WTPCTR_EC' , 'ETPCTR_EC']   
astr=[r'$(a)$',r'$(b)$', r'$(c)$',r'$(d)$',r'$(e)$',r'$(f)$']
dirin="D:/MyPaper/PhD04/Cases/postdata/CTREC/"
cloudtype=["DEEPCONVECTION","ALLCELLS","STRATIFORM","CIRRUS"]
nty=len(cloudtype)
nr=len(casenm)
datestart=datetime.datetime(1993,10,16,0,0,0)
det=datetime.timedelta(minutes=15)            
dateiso=[] 
xtime=[]                
for dt in range(0,ndt):
    dateiso.append(datestart+dt*det)           
for tm in dateiso:
        xtime.append(datetime.datetime.strftime(tm,"%H"))
xdat=range(0,ndt)
#
zdat=[  0.0500000, 0.1643000, 0.3071000, 0.4786000
        , 0.6786000, 0.9071000, 1.1640000, 1.4500000, 1.7640001
        , 2.1070001, 2.4790001, 2.8789999, 3.3069999, 3.7639999
        , 4.2500000, 4.7639999, 5.3070002, 5.8790002, 6.4790001
        , 7.1069999, 7.7639999, 8.4499998, 9.1639996, 9.9069996
        ,10.6800003,11.4799995,12.3100004,13.1599998,14.0500002
        ,14.9600000,15.9099998,16.8799992,17.8799992,18.9099998]
FCDY=np.ndarray(shape=(nz,ndt,nty,nr),dtype=float)
DUQRL=np.ndarray(shape=(nz,ndt,nty,nr),dtype=float)
DUQRS=np.ndarray(shape=(nz,ndt,nty,nr),dtype=float)
DUQL=np.ndarray(shape=(nz,ndt,nty,nr),dtype=float)
DUQI=np.ndarray(shape=(nz,ndt,nty,nr),dtype=float)
DUOMG=np.ndarray(shape=(nz,ndt,nty,nr),dtype=float)
MPQRL=np.ndarray(shape=(nz,nty,nr),dtype=float)
MPQRS=np.ndarray(shape=(nz,nty,nr),dtype=float)
MPQL=np.ndarray(shape=(nz,nty,nr),dtype=float)
MPQI=np.ndarray(shape=(nz,nty,nr),dtype=float)
MPOMG=np.ndarray(shape=(nz,nty,nr),dtype=float)
BINMAX=np.ndarray(shape=(nz,nbin,nty,nr),dtype=float)
#  could water bin must be same as that in getdataforplot.f90
WATERBIN=np.ndarray(shape=(nbin),dtype=float)
WATERBIN[0]=0.005  # MIN
C1=WATERBIN[0]
IC1=0
C3=0.01
IC2=-1
IC3=-1
IC4=-1
#xbin=range(0,nbin)
for I in range(1, nbin):
    C2=C1+(I-IC1)*C3
    if C2>=0.005 and C2 <0.5 and IC2< 0 :
        C1=0.005 ; IC1=I ; C3=0.05
        IC2=1    #! MAKE SURE THE IF BLOCK JUST CALLED ONCE, SO IC1 IS RIGHT
    elif C2>=0.5 and C2<1. and IC3< 0 :
        C1=0.5 ; IC1=I ; C3=0.05
        IC3=1
    elif C2>=10 and IC4< 0 :
        C1=10. ; IC1=I ; C3=1.
        IC4=1
    WATERBIN[I]=C2
xbin=WATERBIN
xbinstr=[]
for xx in WATERBIN:
    xbinstr.append("%e"%xx)
#
#------------------------------------------------------------------------------
for rg in range(0,nr):
    for ty in range(0,nty):
        fpath=dirin+casenm[rg]+"_"+cloudtype[ty]+"_GETPLOTF90.TXT"
        iskp=0 ; nrl=6*ndt+iskp
        onedim=readAscii(fpath,iskp,nrl)
        for it in range(0,ndt):
            for k in range(0,nz):
                kk=it*6*nz+k+nz*0
                FCDY[k,it,ty,rg]=onedim[kk]
                kk=it*6*nz+k+nz*1
                DUQRL[k,it,ty,rg]=onedim[kk]
                kk=it*6*nz+k+nz*2
                DUQRS[k,it,ty,rg]=onedim[kk]
                kk=it*6*nz+k+nz*3
                DUQL[k,it,ty,rg]=onedim[kk]
                kk=it*6*nz+k+nz*4
                DUQI[k,it,ty,rg]=onedim[kk]
                kk=it*6*nz+k+nz*5
                DUOMG[k,it,ty,rg]=onedim[kk]
        del onedim
        iskp=nrl; nrl=iskp+5
        onedim=readAscii(fpath,iskp,nrl)
        for k in range(0,nz):
            kk=k+nz*0
            MPQRL[k,ty,rg]=onedim[kk]
            kk=k+nz*1
            MPQRS[k,ty,rg]=onedim[kk]
            kk=k+nz*2
            MPQL[k,ty,rg]=onedim[kk]
            kk=k+nz*3
            MPQI[k,ty,rg]=onedim[kk]
            kk=k+nz*4
            MPOMG[k,ty,rg]=onedim[kk]
        del onedim
        iskp=nrl ; nrl=iskp+nbin
        onedim=readAscii(fpath,iskp,nrl)
        for ib in range(0,nbin):
            for k in range(0,nz):
                kk=ib*nz+k
                BINMAX[k,ib,ty,rg]=onedim[kk]
        del onedim
#----------- end of reading data--------------------------------------------                
cloudlevs=[2,5,10,15,20,30,40,50,60,70,80,90,100,110]
cloudclors=['w','lightgray','plum','darkorchid','b','dodgerblue','skyblue','aqua',
            'lime','greenyellow','yellow','salmon','pink','orangered','r','darkred']         
# 
#FCDY[it,k,ty,rg]=onedim[kk]
fig,ax=plt.subplots(nrows=nr,ncols=nty,figsize=(12,18))
titlename=r"Frequency of all cloud cells ($10^{-2}%$)"
font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 14,
        }        
ir=0
jc=0
ij=1  
for i in range(0,nr):
    if casenm[i][0:3] == "MLY" :
        regioname= casenm[i][0:4]
    else:
        regioname= casenm[i][0:3]
    if jc==nty:
        jc=0
        ir=ir+1
    for j in range(0,nty):        
        plt.subplot(nr,nty,ij)   # plot ax[i,j]
        ax[ir,jc]=plt.contourf(xdat,zdat,FCDY[:,:,j,i]*100.,colors=cloudclors, levels=cloudlevs,extend='both')
#        plt.colorbar(orientation='horizontal',extend='both',
#                     extendfrac='auto',  spacing='uniform')                           
        marknm=regioname+cloudtype[j]
        plt.title(marknm,fontsize=12)                          
#        plt.axis([0, 16, 0, ndt])
        if ij in range(1,22,4):
            plt.ylabel(r'Cloud Top Height ($km$)', fontdict=font)     
        axx=fig.add_subplot(nr,nty,ij)                         
        axx.set_xticks(range(0,ndt,12))
        if ij in range(21,25):   
            xticklabels = [xtime[nn] for nn in range(0,ndt,12)] 
            axx.set_xticklabels(xticklabels, size=14)        
        plt.show()
        ij=ij+1
        jc=jc+1    
cax = fig.add_axes([0.2, 0.08, 0.6, 0.04])
fig.colorbar(ax[0,0], cax,extend='both',
             spacing='uniform', orientation='horizontal')
plt.show()
plt.savefig(dirpic+'AllCASES_cloudtopduinalcycle.png',dpi=300)        
plt.show()
plt.close() 
#-------------------deep convection -------------------------------------------
cloudlevs=[2,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,110]
cloudclors=['w','lightgray','plum','darkorchid','darkviolet','b','dodgerblue','skyblue','aqua',
            'greenyellow','lime','limegreen','yellow','darkorange','tomato','r']       
# 
#FCDY[it,k,ty,rg]=onedim[kk]
fig,ax=plt.subplots(nrows=2,ncols=3,figsize=(12,8))
titlename=r"Frequency of all cloud cells ($10^{-2}%$)"
font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 14,
        }        
ir=0
jc=0
ij=1  
for i in range(0,nr):
    if casenm[i][0:3] == "MLY" :
        regioname= casenm[i][0:4]
    else:
        regioname= casenm[i][0:3]
    if jc==3:
        jc=0
        ir=ir+1
    j=0 # for deep convection        
    plt.subplot(2,3,ij)   # plot ax[i,j]
    ax[ir,jc]=plt.contourf(xdat,zdat,FCDY[:,:,j,i]*100.,colors=cloudclors, levels=cloudlevs,extend='both')
#        plt.colorbar(orientation='horizontal',extend='both',
#                     extendfrac='auto',  spacing='uniform')                           
    marknm=astr[i]+' '+ regioname
    plt.title(marknm,fontsize=14)
    axx=axx=fig.add_subplot(2,3,ij)
    plt.axis([0, 96, 5, 18])                           
    if jc==0:
        plt.ylabel(r'Cloud Top Height ($km$)', fontdict=font)                        
    axx.set_xticks(range(0,ndt,12))  
    xticklabels = [xtime[nn] for nn in range(0,ndt,12)] 
    axx.set_xticklabels(xticklabels, size=14,rotation=90)
    if ir==1:
        plt.xlabel(r'UTC', fontdict=font)        
    if i in(1,2,4,5)  :
        for tick in axx.yaxis.get_major_ticks():
            tick.label1On = False
    if ir==0  :
        for tick in axx.xaxis.get_major_ticks():
            tick.label1On = False    
    plt.show()
    ij=ij+1
    jc=jc+1    
plt.subplots_adjust(left = 0.1, wspace = 0.1, hspace = 0.2, \
    bottom = 0.25, top = 0.90)
cax = fig.add_axes([0.1, 0.08, 0.8, 0.04])
fig.colorbar(ax[0,0], cax,extend='both',
             spacing='uniform', orientation='horizontal')                                                     
titlename=r"Frequency of deep convection top ($10^{-2}%$)"
plt.title(titlename,fontsize=14)
plt.show()
plt.savefig(dirpic+'deepconvection_cloudtopduinalcycle.png',dpi=300)        
plt.show()
plt.close()            
#---------------profiles ------------------------------------------------------
colors=["b","b","b"]
sty=["solid",'dotted','solid']
width=[1.2,1.2,2]
mker=["","",""]
fig,ax=plt.subplots(nrows=nr,ncols=nty,figsize=(21,8))
ir=0
jc=0
ij=1  
for i in range(0,nr):
    if casenm[0:3] == "MLY" :
        regioname= casenm[i][0:4]
    else:
        regioname= casenm[i][0:3]
    if jc==nty:
        jc=0
        ir=ir+1
    for j in range(0,nty):          
        plt.subplot(nr,nty,ij)   # plot ax[i,j]
        plt.ylim(0,16)
        #plt.xlim(0,0.08)
        ax[ir,jc].plot(MPQRL[:,j,i],zdat,label="Longwave heating rate",
                 c=colors[0],ls=sty[0],marker=mker[0],lw=width[0],)
        ax[ir,jc].plot(MPQRS[:,j,i],zdat,label="Shortwave heating rate",
                 c=colors[1],ls=sty[1],marker=mker[1],lw=width[1],) 
        marknm=regioname
        plt.title(marknm,fontsize=12)  
        if ij in range(1,22,4):
            plt.ylabel(r'Height ($km$)', fontdict=font)            
        plt.show()
        ij=ij+1
        jc=jc+1   
plt.show()
plt.savefig(dirpic+'AllCases_heatingrate.png',dpi=300)        
plt.show()
plt.close()
#        
fig,ax=plt.subplots(nrows=nr,ncols=nty,figsize=(21,8))
ij=1  
ir=0
jc=0 
for i in range(0,nr):
    if casenm[0:3] == "MLY" :
        regioname= casenm[i][0:4]
    else:
        regioname= casenm[i][0:3]
    if jc==nty:
        jc=0
        ir=ir+1
    for j in range(0,nty):          
        plt.subplot(nr,nty,ij)   # plot ax[i,j]
        plt.ylim(0,16)
        #plt.xlim(0,0.08)
        ax[ir,jc].plot(MPQL[:,j,i],zdat,label="Liquid",
                 c=colors[0],ls=sty[0],marker=mker[0],lw=width[0],)
        ax[ir,jc].plot(MPQI[:,j,i],zdat,label="Ice",
                 c=colors[1],ls=sty[1],marker=mker[1],lw=width[1],) 
        ax[ir,jc].plot(MPQI[:,j,i]+MPQL[:,j,i],zdat,label="Ice",
                 c=colors[2],ls=sty[2],marker=mker[2],lw=width[2],) 
        marknm=regioname
        plt.title(marknm,fontsize=12)  
        if ij in range(1,22,4):
            plt.ylabel(r'Height ($km$)', fontdict=font)            
        plt.show()
        ij=ij+1 
        jc=jc+1
plt.show()
plt.savefig(dirpic+'AllCases_meanprfs.png',dpi=300)        
plt.show()
plt.close()          
#
fig,ax=plt.subplots(nrows=nr,ncols=nty,figsize=(21,8))
ij=1
ir=0
jc=0 
for i in range(0,nr):
    if casenm[0:3] == "MLY" :
        regioname= casenm[i][0:4]
    else:
        regioname= casenm[i][0:3]
    if jc==nty:
        jc=0
        ir=ir+1
    for j in range(0,nty):          
        plt.subplot(nr,nty,ij)   # plot ax[i,j]
        plt.ylim(0,16)
        #plt.xlim(0,0.08)
        ax[ir,jc].plot(MPOMG[:,j,i],zdat,label="Vertical velocity",
                 c=colors[0],ls=sty[0],marker=mker[0],lw=width[0],)
        marknm=regioname
        plt.title(marknm,fontsize=12)  
        if ij in range(1,22,4):
            plt.ylabel(r'Height ($km$)', fontdict=font)            
        plt.show()
        ij=ij+1
        jc=jc+1    
plt.show()
plt.savefig(dirpic+'Allcases_meanomgprfs.png',dpi=300)        
plt.show()
plt.close()                     
#------------------------------------------------------------------------------            
fig,ax=plt.subplots(nrows=nr,ncols=nty,figsize=(21,8))
#maxlevs=[]
ij=1   
ir=0
jc=0
for i in range(0,nr):
    if casenm[0:3] == "MLY" :
        regioname= casenm[i][0:4]
    else:
        regioname= casenm[i][0:3]
    if jc==nty:
        jc=0
        ir=ir+1
    for j in range(0,nty):          
        plt.subplot(nr,nty,ij)   # plot ax[i,j]
        ax[ir,jc]=plt.contourf(WATERBIN,zdat,BINMAX[:,:,j,i],cmap=cm.YlOrRd, extend='both')
#        ax[ij]=plt.contourf(xdat,zdat,BINMAX[:,:,j,i],cmap=cm.Greys, levels=maxlevs,extend='both')  levels=maxlevs,          
        marknm=regioname+cloudtype[j]
        plt.title(marknm,fontsize=12)                          
#        plt.axis([0, 16, 0, nbin])
        if ij in range(1,22,4):
            plt.ylabel(r'Height ($km$)', fontdict=font)     
        axx=fig.add_subplot(nr,nty,ij)                         
        axx.set_xticks(range(0,ndt,2))
        if ij in range(21,25):   
            xticklabels = [xbinstr[nn] for nn in range(0,nbin,2)] 
            axx.set_xticklabels(xticklabels, rotation=90, size=12)        
        plt.show()
        ij=ij+1 
        jc=jc+1   
cax = fig.add_axes([0.2, 0.08, 0.6, 0.04])
fig.colorbar(ax[0,0], cax,extend='both',
              spacing='uniform', orientation='horizontal')
plt.show()
plt.savefig(dirpic+'AllCases_maxcwcvsheight.png',dpi=300)        
plt.show()
plt.close()             