#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 20:59:11 2015

@author: jhchen
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import string
import numpy as np
import os
from pylab import *
os.system("cls")
#################################
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
#################################
def readAscii(fpath,iskp,*nl):
    #iskp  the total line skipped of the file
    # fpath   the full path of the file
    # usage: onedim=readAscii(fpaht,iskp)
    onedim=[]
    linesplit=[]
    f=open(fpath)
    if nl:
        nrl=nl[0]
        ff=f.readlines()[iskp:nrl]  ## first line in obs file is legend 
    else:
        ff=f.readlines()[iskp:]
    for line in ff:
        line=string.lstrip(line)
        linesplit.append(line[:-1].split(' '))
    for lnstrs in linesplit:
        for strs in lnstrs:
#            if strs!='':
#                onedim.append(string.atof(strs))
            try:
                onedim.append(string.atof(strs))
            except:
                strs=""
#            else:
#                onedim.append(string.atof(strs))
    del linesplit,ff
    f.close()
    print len(onedim)
    return onedim
#########################-----------------------------------------------------
nt=121
nz=52
ntt=2880
CASE=['PRDCTR_EC','MLYRCTR_EC','NPCCTR_EC','NECCTR_EC','WTPCTR_EC','ETPCTR_EC']
astr=[r'$(a)$',r'$(b)$', r'$(c)$',r'$(d)$',r'$(e)$',r'$(f)$']
nx=200
nga=len(CASE)
dirs='D:/MyPaper/PhD04/Cases/'
dirout='D:/MyPaper/PhD04/Pics/'
ydat_r=[ -50.000 ,    50.000 ,   164.286,    307.143,    478.571  ,  678.571 ,
      907.143 ,  1164.286,   1450.000,   1764.286 ,  2107.143,   2478.572 ,
      2878.572,   3307.143,  3764.286,  4250.000,   4764.286,   5307.143, 
      5878.571,   6478.571,   7107.143,  7764.286,  8450.000,  9164.285,  
      9907.143,  10678.570,  11478.570,  12307.143,  13164.285,  14050.000,
      14964.285,  15907.143,  16878.572,  17878.572,  18907.145,  19964.285,
      21050.000,  22164.285,  23307.145,  24478.572,  25678.572,  26907.145,
      28164.285,  29450.000,  30764.285,  32107.145,  33478.570,  34878.570,
      36307.141,  37764.285,  39250.000,  40750.000]
zdat=[]
for yd in ydat_r:
    zdat.append(yd*0.001) 
upflux=np.ndarray(shape=(nga,nz,nt), dtype=float)
dnflux=np.ndarray(shape=(nga,nz,nt), dtype=float)
upflux_pf=np.ndarray(shape=(nga,nz), dtype=float)
dnflux_pf=np.ndarray(shape=(nga,nz), dtype=float)
alldatestr=[]
for iga in range(0,nga):
    casenm=CASE[iga]
    if casenm[0:3]=='ETP':
        area=casenm[0:3]   
        iy,im,jd=2010,6,3
        atr=r'$(f)$'
    if casenm[0:3]=='WTP':
        area=casenm[0:3]    
        iy,im,jd=2010,7,3
        atr=r'$(e)$'    
    if casenm[0:3]=='NPC':
        area=casenm[0:3]   
        iy,im,jd=2010,8,2
        atr=r'$(c)$'
    if casenm[0:3]=='PRD':
        area=casenm[0:3]   
        iy,im,jd=2012,4,1
        atr=r'$(a)$'
    if casenm[0:3]=='MLY':
        area=casenm[0:4]
        iy,im,jd=2010,6,2 
        atr=r'$(b)$'
    if casenm[0:3]=='NEC':
        area=casenm[0:3]   
        iy,im,jd=2012,7,6 
        atr=r'$(d)$'
    folds="/CTREC"+"%4d"%iy+"%2.2d"%im+"%2.2d"%jd+"/Simulation/"
    datestr="%4d"%iy+"%2.2d"%im+"%2.2d"%jd+'_031d'
    dirpic='D:/MyPaper/PhD04/Pics/'
    dirin=dirs+area+folds
    dirq12=dirin
    nameq1q2=area+'_'+"%4d"%iy+"%2.2d"%im+"%2.2d"%jd +"_031d_regrid_ERA.99"
    nameforcing=area+'_'+"%4d"%iy+"%2.2d"%im+"%2.2d"%jd +"_031d_lsforcing_regrid_ERA.37"
    ydat=[]
    for i in range(0,nz):
        ydat.append(i*0.5) 
####### file 1
    filenm=casenm+'_meanmassflux_updown.txt'
    fpath=dirin+filenm
    iskp=0
    onedim=readAscii(fpath, iskp)
    for it in range(0,nt):
        for iz in range(0,nz):
            k=iz+2*it*nz
            upflux[iga,iz,it]=onedim[k]
            kk=iz+2*it*nz+nz
            dnflux[iga,iz,it]=onedim[kk]
    for iz in range(0,nz):
        upflux_pf[iga,iz]=upflux[iga,iz,:].mean()
        dnflux_pf[iga,iz]=dnflux[iga,iz,:].mean()
    datestart=datetime.datetime(iy,im,jd,0,0,0)
    det=datetime.timedelta(hours=6)            
    dateiso=[]            
    for dt in range(0,nt):
        dateiso.append(datestart+dt*det)
    xdate=[]    
    xdat=range(0,nt)            
    for tm in dateiso:
        xdate.append(datetime.datetime.strftime(tm,"%b/%d")) 
    alldatestr.append(xdate)    
    #
xdat=range(0,nt) 
###########################################################################
font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 14,
        }  
fig,ax=plt.subplots(nrows=2,ncols=3,figsize=(21,8))
color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
wd=[2,2,2,2,2]
ir=0
jc=0
ij=1
for i in range(0,nga):
    if CASE[i][0:3] == "MLY" :
        regioname= CASE[i][0:4]
    else:
        regioname= CASE[i][0:3]
    if jc==3:
        jc=0
        ir=ir+1
    datestr=alldatestr[i]    
    plt.subplot(2,3,ij)   # plot ax[i,j]
    plt.ylim(0,16)
    ax[ir,jc]=plt.contour(xdat,zdat,dnflux[i,:,:],colors='cyan', extend='both')
    plt.ylim(0,16)
    ax[ir,jc]=plt.contourf(xdat,zdat,upflux[i,:,:],cmap=cm.OrRd)
#        plt.colorbar(orientation='horizontal',extend='both',
#                     extendfrac='auto',  spacing='uniform')                           
    marknm=astr[i]+' '+ regioname
    plt.title(marknm,fontsize=14)
    axx=fig.add_subplot(2,3,ij)
    ymajorLocator   = MultipleLocator(4)
    axx.yaxis.set_major_locator(ymajorLocator)                          
    if jc==0:
        plt.ylabel(r'Height ($km$)', fontdict=font)
    if ir==1:
        plt.xlabel(r'Mass Flux ($km$)', fontdict=font)                        
    axx.set_xticks(range(0,nt,16))  
    xticklabels = [datestr[nn] for nn in range(0,nt,16)] 
    axx.set_xticklabels(xticklabels, size=14,rotation=90)       
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
    bottom = 0.3, top = 0.90)
cax = fig.add_axes([0.1, 0.08, 0.8, 0.04])
fig.colorbar(ax[0,0], cax,extend='both',
             spacing='uniform', orientation='horizontal')                                                     
titlename=r"Mass Flux ($kg$ $m^{-2}$ $hr^{-1}$)"
plt.title(titlename,fontsize=14)
plt.show()
plt.savefig(dirpic+'MassFlux_VS_date.png',dpi=300)        
plt.show()
plt.close()
#
fig,ax=plt.subplots(nrows=2,ncols=3,figsize=(10,12))
color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
wd=[2,2,2,2,2]
ir=0
jc=0
ij=1
for i in range(0,nga):
    if CASE[i][0:3] == "MLY" :
        regioname= CASE[i][0:4]
    else:
        regioname= CASE[i][0:3]
    if jc==3:
        jc=0
        ir=ir+1
    datestr=alldatestr[i] 
    ax[ir,jc].plot(dnflux_pf[i,:],zdat,c=color_cycle[1], lw=3,label=r'Down')
    ax[ir,jc].plot(upflux_pf[i,:],zdat,c=color_cycle[0], lw=3,label=r'Up')
#        plt.colorbar(orientation='horizontal',extend='both',
#                     extendfrac='auto',  spacing='uniform')                           
    marknm=astr[i]+' '+ regioname
    ax[ir,jc].set_title(marknm,fontsize=14)
    ax[ir,jc].set_ylim(0,16)
    ax[ir,jc].set_xlim(-5,9)
    ymajorLocator   = MultipleLocator(4)
    ax[ir,jc].yaxis.set_major_locator(ymajorLocator) 
    xmajorLocator   = MultipleLocator(4)
    ax[ir,jc].xaxis.set_major_locator(xmajorLocator)                          
    if jc==0:
        ax[ir,jc].set_ylabel(r'Height ($km$)', fontdict=font)
    if ir==1:
        ax[ir,jc].set_xlabel(r'Mass Flux', fontdict=font)                              
    if i in(1,2,4,5)  :
        for tick in ax[ir,jc].yaxis.get_major_ticks():
            tick.label1On = False
    if ir==0  :
        for tick in ax[ir,jc].xaxis.get_major_ticks():
            tick.label1On = False    
    if ij==6:
        legend(loc=1,frameon=False)
    plt.show()
    ij=ij+1
    jc=jc+1    
plt.show()
plt.savefig(dirpic+'MassFlux_profiles.png',dpi=300)        
plt.show()
plt.close() 




