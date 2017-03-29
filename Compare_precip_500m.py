#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 21:22:34 2015

@author: jhchen
"""
import string
import matplotlib.pyplot as plt
import numpy as np
import datetime
from pylab import *
import scipy.stats as scista
mpl.rcParams['ytick.labelsize'] = 18
#
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
#
CASE=['PRDCTR_H','MLYRCTR_H','NPCCTR_H','NECCTR_H','WTPCTR_H','ETPCTR_H',]
nt=2880
nx=1202
nga=len(CASE)
dirs='Z:/CRM/500m/'
dirobs='D:/MyPaper/PhD04/Data/TRMM/3B40/'
direc='D:/MyPaper/PhD04/Cases/ERA/FORCING/'
dirout='D:/MyPaper/PhD04/Pics/'
fig,ax=plt.subplots(nrows=2,ncols=3,figsize=(18,9))
#fig,axs=plt.subplots(nrows=2,ncols=3,figsize=(12,12))
color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
wd=[2,2,2,2,2]
#color_cycle=['lightgrey', 'k', 'grey', 'y','indigo', 'cyan']
#wd=[4,1.5,1.5,2,2]
ir=0
ic=0
for iga in range(0,nga):
    if ic==3:
        ic=0
        ir=ir+1
    print ir,ic
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
    dirin=dirs+area+'/run/postdata/'
#------------------------------------------------------------------------------
    fpath=dirin+'preci_'+casenm+'.txt'
    iskp=0
    onedim=readAscii(fpath,iskp)
    preci=np.ndarray(shape=(nt),dtype=float)
    pbl=np.ndarray(shape=(nt),dtype=float)
    fsh=np.ndarray(shape=(nt),dtype=float)
    flh=np.ndarray(shape=(nt),dtype=float)
    tmp=np.ndarray(shape=(4),dtype=float)
    for it in range(0,nt):
        tmp[0]=0.0
        tmp[1]=0.0
        tmp[2]=0.0
        tmp[3]=0.0
        for ix in range(0,nx): # for 2-201 200 grids
            if ix!=0 or ix!=nx-1 :
                k=it*(nx*4)+ix+nx*0
                tmp[0]=tmp[0]+onedim[k]*1e3*3600. # convert m/s to mm/hr
                k=it*(nx*4)+ix+nx*1
                tmp[1]=tmp[1]+onedim[k]
                k=it*(nx*4)+ix+nx*2
                tmp[2]=tmp[2]+onedim[k]
                k=it*(nx*4)+ix+nx*3
                tmp[3]=tmp[3]+onedim[k]
            preci[it]=tmp[0]/(nx*1.0)
            pbl[it]=tmp[1]/(nx*1.)
            fsh[it]=tmp[2]/(nx*1.)
            flh[it]=tmp[3]/(nx*1.)
    del onedim
    ntt=nt/(4*3)
    rainsim=np.ndarray(shape=(ntt),dtype=float)
    for it in range(0,nt/(4*3)):
        a=0.
        for i in range(0,12):
            k=it*12+i
            a=a+preci[k]
        rainsim[it]=a/12.
#------------------------------------------------------------------------------
    fpath=dirobs+area+'-'+datestr+'_TRMM3B40.txt'
    iskp=1
    rainobs=readAscii(fpath,iskp)
    for i in range(0,len(rainsim)):
        if rainobs[i]<0:
            rainobs[i]=rainsim[i]
    nl=len(rainobs)
    fpath=direc+area+'/'+area+'_'+datestr+'_SHLH_ERA.43'
    iskp=0
    onedim=readAscii(fpath,iskp)
    rainect=np.ndarray(shape=(nl),dtype=float)
    rainecc=np.ndarray(shape=(nl),dtype=float)
    for it in range(0,nl):
        k=it*4+2
        rainect[it]=onedim[k]
        rainecc[it]=onedim[k+1]
    for i in range(0,len(rainect)):
        if rainect[i]<0:
            rainect[i]=0.    
    datestart=datetime.datetime(iy,im,jd,0,0,0)
    det=datetime.timedelta(hours=3)            
    dateiso=[]            
    for dt in range(0,nl):
        dateiso.append(datestart+dt*det)
    xdate=[]    
    xxx=range(0,nl)            
    for tm in dateiso:
        xdate.append(datetime.datetime.strftime(tm,"%d/%b")) 
#
    ns=min(nl,ntt)
#    ax[ir,ic]=plt.subplot(1,1,1)    
    ax[ir,ic].plot(xxx[0:ns],rainsim[0:ns],color=color_cycle[0],lw=wd[0],label="CRM") 
    ax[ir,ic].plot(xxx[0:ns],rainobs[0:ns],color=color_cycle[1],lw=wd[1],label="TRMM")
    ax[ir,ic].plot(xxx[0:ns],rainect[0:ns],color=color_cycle[2],lw=wd[2],label="ERA")
    tilstr=atr+' '+area
    ax[ir, ic].set_title(tilstr, fontsize=18)
    xx1=[]
    yy1=[]
    yy2=[]
    for irn in range(0,len(rainsim)):
        if rainsim[irn]>0.0:
            xx1.append(rainsim[irn])
            yy1.append(rainobs[irn])
            yy2.append(rainect[irn])
    r1,p1=scista.pearsonr(xx1,yy1)
    r1str="r1= "+"%.2f"%r1
    r2,p2=scista.pearsonr(xx1,yy2)
    r2str="r2= "+"%.2f"%r2
#    plt.ylabel(r'Precipitation ($mm$ $hr^{-1}$)', size=18)
    if iga==5 :
        ax[ir,ic].legend(loc=(0.65,0.65),frameon=False)
    if iga>3 :
        ax[ir,ic].set_yticks(range(0,3))
        ax[ir,ic].text(4,1.8,r1str,fontsize=18)
        ax[ir,ic].text(4,1.6,r2str,fontsize=18)
    else:
        ax[ir,ic].set_yticks(range(0,7))
        ax[ir,ic].text(4,5.4,r1str,fontsize=18)
        ax[ir,ic].text(4,4.8,r2str,fontsize=18)
    ax[ir,ic].set_xticks(range(0,ns,72))
    xticklabels = [xdate[nn] for nn in range(0,ns,72)] 
    ax[ir,ic].set_xticklabels(xticklabels, size=18)
    ic=ic+1
fig.text(0.03, 0.7, r'Precipitation ($mm$ $hr^{-1}$)', ha = 'left',fontsize=24,rotation=90)
fig.subplots_adjust(left=0.1,bottom=0.12,right=1-0.05,top=1-0.1,hspace=0.4)
plt.show()                     
plt.savefig(dirout+'ALLCASE_rainfall_Color_500m.png',dpi=300)          
plt.show()
plt.close()
"""
#---- plot scatter --------------------------------------------------
    fig,axe1=plt.subplots(nrows=1,ncols=1,figsize=(6,6))
    axe1.set_ylim(0,8)
    axe1.set_xlim(0,8)  
    plt.scatter(rainobs[0:ns],rainsim[0:ns],c=color_cycle[0],marker="+",
            label='SIM')
    plt.scatter(rainobs[0:ns],rainect[0:ns],c=color_cycle[4],marker="o",
            linewidths=0.0,label='ERA',alpha=0.75)
    plt.ylabel(r'SIM and ERA', size=20)
    plt.xlabel(r'TRMM(3B40)', size=20)
    plt.show()                     
    plt.savefig(dirout+casenm+'_rainfall_scatter.pdf')          
    plt.show()
    plt.close()
"""