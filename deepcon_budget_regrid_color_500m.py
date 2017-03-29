#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 10:12:17 2015

@author: chenjh
"""
import matplotlib as mpl
import numpy as np
import matplotlib.cm as cm
import datetime
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.dates as matdate
import calendar
import string
import time
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['contour.negative_linestyle'] = 'dashed'
matplotlib.rcParams['ytick.labelsize'] = 16
matplotlib.rcParams['xtick.labelsize'] = 16
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
#--------------------------------------·---------------
add_suffix='_V3'
if add_suffix[2]=='3':
    a,b=30,50
if add_suffix[2]=='2':
    a,b=3,10
if add_suffix[2]=='1':
    a,b=3,40
nz=33
km=nz
CASENM=["PRDCTR_H","MLYRCTR_H","NPCCTR_H","NECCTR_H","WTPCTR_H","ETPCTR_H"]
datestr=['20120401' , '20100602' , '20100802' ,
           '20120706' , '20100703' , '20100603' ]
orderstr=[r'($a$)',r'($b$)',r'($c$)',r'($d$)',r'($e$)',r'($f$)',]
nga=len(CASENM)
dirin="Z:/CRM/500m/"
dirpic="D:/MyPaper/PhD04/Pics/"
#
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
for yd in range (0,33):#ydat_r:
    ydat.append(yd*0.5)
ql=np.ndarray(shape=(km,nga),dtype=float)
qi=np.ndarray(shape=(km,nga),dtype=float)
qrl=np.ndarray(shape=(km,nga),dtype=float)
qrs=np.ndarray(shape=(km,nga),dtype=float)
envdc=np.ndarray(shape=(3,km,nga),dtype=float)
condc=np.ndarray(shape=(5,km,nga),dtype=float)
q1q2com=np.ndarray(shape=(12,km,nga),dtype=float)
for iga in range(0,nga):
    casenm=CASENM[iga]
    if casenm[0:3]=='MLY' :
        areastr=casenm[0:4]
    else:
        areastr=casenm[0:3]
    filename=casenm+'_DeepCons_budgetandprof_regrid'+add_suffix+'.txt'
    fpath=dirin+areastr+'/run/'+filename
    iskp=0
    onedim=readAscii(fpath, iskp)
    alldat=np.ndarray(shape=(24,km),dtype=float)
    for j in range(0,24):
        for i in range(0,km):
            k=j*km+i
            alldat[j,i]=onedim[k]
    ql[:,iga]=alldat[0,:]
    qi[:,iga]=alldat[1,:]
    qrl[:,iga]=alldat[2,:]
    qrs[:,iga]=alldat[3,:]
    for i in range(0,3):
    	j=4+i
    	envdc[i,:,iga]=alldat[j,:]
    for i in range(0,5):
    	j=4+i+3
    	condc[i,:,iga]=alldat[j,:]
    for i in range(0,12):
    	j=4+3+5+i
    	q1q2com[i,:,iga]=alldat[j,:]
    del onedim,alldat
envdc[:,0,:]=0.
condc[:,0,:]=0.
q1q2com[:,0,:]=0.
#----- Plotting ---------------------------------
font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 16,
        } 
fig,ax=plt.subplots(nrows=2,ncols=3,figsize=(9,20))
jr=0
jc=0
ij=1
for iga in range(0,nga):
    casenm=CASENM[iga]
    if casenm[0:3]=='MLY' :
        areastr=casenm[0:4]
    else:
        areastr=casenm[0:3]
    if jc==3:
        jc=0
        jr=jr+1
#    plt.subplot(2,3,ij)
    ax[jr,jc].plot(qrl[:,iga],ydat,label='Longwave',c='k',lw=3)
    ax[jr,jc].plot(qrs[:,iga],ydat,label='Shortwave',c='0.65',lw=4)
    marknm=orderstr[iga]+' '+areastr
    ax[jr,jc].set_title(marknm,fontsize=14)    
    ax[jr,jc].set_ylim(0,16)
    ax[jr,jc].set_xlim(-15,15)
    xmajorLocator   = MultipleLocator(6) #将x或y轴主刻度标签设置为6的倍数  
    ax[jr,jc].xaxis.set_major_locator(xmajorLocator)
    if jr==1:
        ax[jr,jc].set_xlabel(r'Heating rate ($k$ $s^{-1}$)', fontdict=font)
    if jc==0:
        ax[jr,jc].set_ylabel(r'Height ($km$)', fontdict=font)
    if jr in(0,1,2) and jc in(1,2,3):
            setp(ax[jr,jc].get_yticklabels(), visible=False) #    
    jc=jc+1
#    ij=ij+1
plt.show()
plt.savefig(dirpic+'DCC_heatingprofiles_fortran'+add_suffix+'_regrid_Color_500m.png',dpi=300)          
plt.show()
plt.close() 
#
color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
lnstygrey=['-','--','-','--','-','--','-','--','-','--','-','--']
lnwd=[2,2.5,2,2.5,2,2.5,2,2.5,2,2.5,2,2.5,2,2.5]
fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(8,12))
ij=0
for iga in range(0,nga):
    casenm=CASENM[iga]
    if casenm[0:3]=='MLY' :
        areastr=casenm[0:4]
    else:
        areastr=casenm[0:3]
    if jc==3:
        jc=0
        jr=jr+1
#    plt.subplot(2,3,ij)
    ccolor=color_cycle[iga]
    marknm=areastr+' liquid'
    sty=lnstygrey[ij]
    ax.plot(ql[:,iga],ydat,label=marknm,color=ccolor,ls=sty,lw=lnwd[ij]+2)
    ij=ij+1
    marknm=areastr+' ice'
    sty=lnstygrey[ij]
    ax.plot(qi[:,iga],ydat,label=marknm,color=ccolor,ls=sty,lw=lnwd[ij]+2)
    ij=ij+1   
    #plt.ylim(0, 16)
    ax.set_ylim(0,16) 
    ax.set_xlabel(r'Water content ($g$ $kg^{-1}$)', fontdict=font)
    ax.set_ylabel(r'Height ($km$)', fontdict=font)
    jc=jc+1
#    ij=ij+1
plt.legend()
plt.show()
plt.savefig(dirpic+'DCC_qlqi_fortran'+add_suffix+'_regrid_Gray_500m.png',dpi=300)          
plt.show()
plt.close() 
#
matplotlib.rcParams['ytick.labelsize'] = 18
matplotlib.rcParams['xtick.labelsize'] = 18
fig,ax=plt.subplots(nrows=1,ncols=2,figsize=(12,12))
color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
wd=[2,2,2,2,2]
jr=0
jc=0
for iga in range(0,nga):
    if jc==2:
        jc=0
        jr=jr+1
    casenm=CASENM[iga]
    if casenm[0:3]=='MLY':
        area=casenm[0:4]
    else:
        area=casenm[0:3]   
    atr=orderstr[iga]
    lnstycolor=['-','-','-','-','-','-','-','-']
    lncolor=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
    lncolor=['r','g','b','orange','blueviolet','aqua'] #lncolor=
    lnmkcolor=['None','None','None','None','None','None'] 
    lnwidcolor=[4.0,4.0,4.0,4.0,4.0,4.0]  
    lnstygrey=['-',':','-','--','-','-.']
    lngrey=['k','k','gray','gray','lightgrey','lightgrey']
    lnmkgrey='None','None','None','None','None','None'
    lnwidgrey=[3.5,3.5,4.0,4.5,5,5.5]   
    colors=lncolor
    sty=lnstygrey
    mker=lnmkgrey
    width=lnwidcolor 
    q1m=np.ndarray(shape=(km),dtype=float)
    q2m=np.ndarray(shape=(km),dtype=float)
    q2m[:]=(condc[0,:,iga]+condc[1,:,iga]+condc[2,:,iga]+condc[3,:,iga])*2.5e10/2.834e10 \
        +q1q2com[7,:,iga]+q1q2com[8,:,iga]+q1q2com[9,:,iga]
    q1m[:]=condc[0,:,iga]+condc[1,:,iga]+condc[2,:,iga]+condc[3,:,iga]+condc[4,:,iga] \
        + q1q2com[0,:,iga]+q1q2com[1,:,iga]+q1q2com[2,:,iga] +q1q2com[3,:,iga]                          \
        + q1q2com[4,:,iga]+q1q2com[5,:,iga]  
    ax[0].set_ylim(0,16)           
    ax[0].plot(q1m,ydat,label=area,
        color=colors[iga],ls=sty[iga],marker=mker[iga],lw=width[iga],) #
    ax[0].set_title(r'($a$) $Q_1$',fontsize=20)
    #allvar_mean[5,0]=0.
    ax[1].set_ylim(0,16) 
    ax[1].plot(q2m,ydat,label=area , #r'$Q_2$',
        color=colors[iga],ls=sty[iga],marker=mker[iga],lw=width[iga],)  #
    ax[1].set_title(r'($b$) $Q_2$',fontsize=20)
    plt.legend()
#plt.subplots_adjust(left = 0.1, wspace = 0.2, hspace = 0.3, \
#    bottom = 0.1, top = 0.90)
plt.show()                     
plt.savefig(dirpic+'ALLCASE_DCC_Q1Q2'+add_suffix+'_regrid_Color_500m.png',dpi=300)          
plt.show()
plt.close()
matplotlib.rcParams['ytick.labelsize'] = 16
matplotlib.rcParams['xtick.labelsize'] = 16
#------------------------------------------------------------------------------
fig,ax=plt.subplots(nrows=3,ncols=4,figsize=(12,35))
color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
wd=[2,2,2,2,2]
jr=0
jc=0
for iga in range(0,nga):
    if jc==4:
        jc=0
        jr=jr+1
    casenm=CASENM[iga]
    if casenm[0:3]=='MLY':
        area=casenm[0:4]
    else:
        area=casenm[0:3]   
    atr=orderstr[iga]
    lnstycolor=['-','-','-','-']
    lncolor=['orangered','orangered','yellowgreen','yellowgreen']
    lncolor=['r','darkgoldenrod','g','b','darkorchid']
    lncolor=['deeppink','b','green','y']
    lncolor=['r','g','b','orange','blueviolet','aqua']
    lnmkcolor=['None','None','None','None','None'] 
    lnwidcolor=[4.0,4.0,4.0,4.0,4.0]  
    lnstygrey=['-',':','-','--','-']
    lngrey=['k','k','gray','gray','lightgrey']
    lnmkgrey=['None','None','None','None','None'] 
    lnwidgrey=[4.0,4.0,4.0,4.0,4.0]   
    colors=lncolor
    sty=lnstygrey
    mker=lnmkgrey
    width=lnwidcolor
    size_title=18     
    ax[jr,jc].set_ylim(0,16)
    ax[jr,jc].set_xlim(-50,90)            
    ax[jr,jc].plot(q1q2com[0,:,iga],ydat,label=r'$Q_1$e',
        color=colors[0],ls=sty[0],marker=mker[0],lw=width[0],) #
    #allvar_mean[5,0]=0.
    ax[jr,jc].plot(q1q2com[1,:,iga]+q1q2com[2,:,iga],ydat,label=r'$Q_1$d'+'\n'+'+ $Q_1$s',
        color=colors[1],ls=sty[1],marker=mker[1],lw=width[1],)  #
    q1cm=np.ndarray(shape=(km),dtype=float)
    q1cm[:]=condc[0,:,iga]+condc[1,:,iga]+condc[2,:,iga]+condc[3,:,iga]+condc[4,:,iga]
    ax[jr,jc].plot(q1cm,ydat,label=r'$Q_1$c',
        color=colors[2],ls=sty[2],marker=mker[2],lw=width[2],)  #
    radsim=np.ndarray(shape=(km),dtype=float)
    radsim[:]=q1q2com[4,:,iga]+q1q2com[5,:,iga]
    ax[jr,jc].plot(radsim,ydat,label=r'$Q_R$',
        color=colors[3],ls=sty[3],marker=mker[3],lw=width[3],)   #q
    #ax[ir,ic].set_title('Case '+casenm+r'   $Q_1$ and $Q_2$'+ r' ($K$ $d^{-1}$)',fontsize=size_title)
    titlestr=atr+" "+area+r' $Q_1$'# ($K$ $day^{-1}$)'
    ax[jr,jc].set_title(titlestr,fontsize=size_title)
    xmajorLocator   = MultipleLocator(a) #将y轴主刻度标签设置为2的倍数  
#    ymajorFormatter = FormatStrFormatter('%1.1f') #设置y轴标签文本的格式 
    ax[jr,jc].xaxis.set_major_locator(xmajorLocator) 
    ymajorLocator   = MultipleLocator(4) 
    ax[jr,jc].yaxis.set_major_locator(ymajorLocator)
    if jr==1 and jc==3 :
        ax[jr,jc].legend(loc=(0.97,0.2),frameon=False)
    if jr==0 and jc==2 :
        ax[jr,jc].legend(loc=(2.17,0.1),frameon=False)
    if jc==0:
        ylabs='Height'+r' ($km$)'
        ax[jr,jc].set_ylabel(ylabs,fontsize=size_title)
    if jr in(0,1,2) and jc in(1,2,3):
            setp(ax[jr,jc].get_yticklabels(), visible=False) #
    # Q2
    jc=jc+1
    lnstycolor=['-','-','-','-']
    lncolor=['r','g','b','orange','blueviolet','aqua']
    lnmkcolor=['None','None','None','None','None'] 
    lnwidcolor=[4.0,4.0,4.0,4.0,4.0]  
    lnstygrey=['-',':','-','--','-']
    lngrey=['k','k','gray','gray','lightgrey']
    lnmkgrey=['None','None','None','None','None'] 
    lnwidgrey=[4.0,4.0,4.0,4.0,4.0]   
    colors=lncolor
    sty=lnstygrey
    mker=lnmkgrey
    width=lnwidcolor 
    ax[jr,jc].set_ylim(0,16) 
    ax[jr,jc].set_xlim(-50,90)          
    ax[jr,jc].plot(q1q2com[7,:,iga],ydat,label=r'$Q_2$e',
        color=colors[0],ls=sty[0],marker=mker[0],lw=width[0],) #
    #allvar_mean[5,0]=0.
    ax[jr,jc].plot(q1q2com[8,:,iga]+q1q2com[9,:,iga],ydat,label=r'$Q_2$d'+'\n'+'+ $Q_2$s',
        color=colors[1],ls=sty[1],marker=mker[1],lw=width[1],)  #
    q2cm=np.ndarray(shape=(km),dtype=float)
    q2cm[:]=(condc[0,:,iga]+condc[1,:,iga]+condc[2,:,iga]+condc[3,:,iga])*2.5e10/2.834e10
    ax[jr,jc].plot(q2cm,ydat,label=r'$Q_2$c',
        color=colors[2],ls=sty[2],marker=mker[2],lw=width[2],)  #
    #ax[ir,ic].set_title('Case '+casenm+r'   $Q_1$ and $Q_2$'+ r' ($K$ $d^{-1}$)',fontsize=size_title)
    titlestr=atr+" "+area+r' $Q_2$' # ($K$ $day^{-1}$)'
    ax[jr,jc].set_title(titlestr,fontsize=size_title)
    xmajorLocator   = MultipleLocator(a) #将y轴主刻度标签设置为2的倍数  
#    ymajorFormatter = FormatStrFormatter('%1.1f') #设置y轴标签文本的格式 
    ax[jr,jc].xaxis.set_major_locator(xmajorLocator) 
    ymajorLocator   = MultipleLocator(4) 
    ax[jr,jc].yaxis.set_major_locator(ymajorLocator) 
    if jr==1 and jc==3 :
        ax[jr,jc].legend(loc=(0.97,0.2),frameon=False)
    if jr==0 and jc==2 :
        ax[jr,jc].legend(loc=(2.17,0.1),frameon=False)
    if jc==0:
        ylabs='Height'+r' ($km$)'
        ax[jr,jc].set_ylabel(ylabs,fontsize=size_title)
    if jr in(0,1,2) and jc in(1,2,3):
            setp(ax[jr,jc].get_yticklabels(), visible=False) #
    jc=jc+1
plt.subplots_adjust(left = 0.1, wspace = 0.2, hspace = 0.3, \
    bottom = 0.1, top = 0.90)
plt.show()                     
plt.savefig(dirpic+'ALLCASE_DCC_Q1Q2Comps'+add_suffix+'_regrid_Color_500m.png',dpi=300)          
plt.show()
plt.close()
#
fig,ax=plt.subplots(nrows=2,ncols=3,figsize=(10,18))
color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
wd=[2,2,2,2,2]
iro=0
ic=0
for iga in range(0,nga):
    if ic==3:
        ic=0
        iro=iro+1
    casenm=CASENM[iga]
    if casenm[0:3]=='MLY':
        area=casenm[0:4]
    else:
        area=casenm[0:3]   
    atr=orderstr[iga]  
##############################
# Q1 and Q2
###################
    ####plot setup 
    lnstycolor=['-','-','-','-','-']
    lncolor=['orangered','orangered','yellowgreen','yellowgreen']
    lncolor=['r','darkgoldenrod','g','b','darkorchid']
    lncolor=['lime','b','green','y','magenta']
    lncolor=['r','g','b','orange','blueviolet','aqua']
    lnmkcolor=['None','None','None','None','None'] 
    lnwidcolor=[4.0,4.0,4.0,4.0,4.0,4.0,4.0]  
    lnstygrey=['-',':','-','--','-']
    lngrey=['k','k','gray','gray','lightgrey']
    lnmkgrey=['None','None','None','None','None'] 
    lnwidgrey=[4.0,4.0,4.0,4.0,4.0]   
    colors=lncolor
    sty=lnstygrey
    mker=lnmkgrey
    width=lnwidcolor 
    micro_com=np.ndarray(shape=(5,km),dtype=float)
    micro_com[:,:]=condc[:,:,iga]
    size_title=18     
    ax[iro,ic].set_ylim(0,16)
    ax[iro,ic].set_xlim(-80,130)            
    ax[iro,ic].plot(micro_com[0,:],ydat,label=r'$Cond$',
        color=colors[0],ls=sty[0],marker=mker[0],lw=width[0]) #Condensation
    #allvar_mean[5,0]=0.
    ax[iro,ic].plot(micro_com[1,:],ydat,label=r'$Evap$',
        color=colors[1],ls=sty[1],marker=mker[1],lw=width[1])  #Evaporation
    ax[iro,ic].plot(micro_com[2,:],ydat,label=r'$Dep+Frez$',
        color=colors[2],ls=sty[2],marker=mker[2],lw=width[2])   #Deposition
    ax[iro,ic].plot(micro_com[3,:],ydat,label=r'$Sub$',
        color=colors[3],ls=sty[3],marker=mker[3],lw=width[3])  #Sublimation
    ax[iro,ic].plot(micro_com[4,:],ydat,label=r'$Fus$',
        color=colors[4],ls=sty[4],marker=mker[4],lw=width[4])  #Fustion
    #ax[ir,ic].set_title('Case '+casenm+r'   $Q_1$ and $Q_2$'+ r' ($K$ $d^{-1}$)',fontsize=size_title)
    titlestr=atr+" "+area+r' ($K$ $day^{-1}$)'
    ax[iro,ic].set_title(titlestr,fontsize=size_title)
    xmajorLocator   = MultipleLocator(b) #将y轴主刻度标签设置为3的倍数  
    ax[iro,ic].xaxis.set_major_locator(xmajorLocator) 
    if iro==0 and ic==2 :
        ax[iro,ic].legend(loc=(1.0,0.5),frameon=False)
    if ic==0:    
        ylabs='Height'+r' ($km$)'
        ax[iro,ic].set_ylabel(ylabs,fontsize=size_title)
    if iro in(0,1,2) and ic in(1,2,3):
            setp(ax[iro,ic].get_yticklabels(), visible=False) #
    # Q2
    ic=ic+1
plt.subplots_adjust(left = 0.1,right=0.85, wspace = 0.2, hspace = 0.25, \
    bottom = 0.1, top = 0.90)
plt.show()                     
plt.savefig(dirpic+'ALLCASE_DCC_Q1Q2CondComps'+add_suffix+'_regrid_Color_500m.png',dpi=300)          
plt.show()
#plt.close()
