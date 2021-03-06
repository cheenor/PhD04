#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 08:21:53 2015

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
##################################################
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['contour.negative_linestyle'] = 'dashed'
matplotlib.rcParams['ytick.labelsize'] = 16
matplotlib.rcParams['xtick.labelsize'] = 16
plt.rc('lines', linewidth=4)
###################################################
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
nt=121
nz=52
CASE=['PRDCTR_EC','MLYRCTR_EC','NPCCTR_EC','NECCTR_EC','WTPCTR_EC','ETPCTR_EC']
nx=202
nga=len(CASE)
dirin1='D:/MyPaper/PhD04/Cases/ERA/FORCING/'
dircnrain='D:/MyPaper/PhD04/Data/RainCN05/'
pic_out='D:/MyPaper/PhD04/Pics/'
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
for yd in ydat_r:
    ydat.append(yd*0.001)
fig,ax=plt.subplots(nrows=3,ncols=2,figsize=(18,12))
#fig,axs=plt.subplots(nrows=2,ncols=3,figsize=(12,12))
color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
wd=[2,2,2,2,2]
jc=0
jr=0
ij=1
for iga in range(0,nga):
    if jc==2:
        jc=0
        jr=jr+1
    casenm=CASE[iga]
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
    dirin=dirin1+namestr+'/'
    cnname=namestr+datestr+".txt"
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
    charsize=18
    font = {'family' : 'serif',
#         'serif' : 'Bookman',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 18,
        } 
    plt.subplot(3,2,ij)
#    plt.colorbar(orientation='horizontal',extend='both',
#                 extendfrac='auto',  spacing='uniform') 
    ax[jr,jc]=plt.contour(xxx,ydat,fcq2,colors=colorq2,
        linewidths=1,levels=levq2,linestyles=linetypq2)                           
#plt.title('Observation',fontsize=charsize)                       
    plt.axis([0, nt, 0, 16])  ## x axis  y axis
    plt.clabel(ax[jr,jc],inline=1,fmt='%1d',fontsize=charsize-2)
    ax[jr,jc]=plt.contourf(xxx,ydat,fcq1,colors=colorq1,
                  levels=levq1,extend='both')                                                   
    axx=fig.add_subplot(3,2,ij)                      
    axx.set_xticks(range(0,nt,20))
    xticklabels = [xdate[nn] for nn in range(0,nt,20)] 
    axx.set_xticklabels(xticklabels, size=charsize)
    ymajorLocator   = MultipleLocator(4) 
    axx.yaxis.set_major_locator(ymajorLocator)    
    text1=marktr+' '+namestr+' ('+yearstr+')'
    axx.text(1.5,16.5,text1,fontsize=charsize+2) 
#    text1=namestr+' ('+yearstr+')'#+r" Temperature ($K$ $d^{-1}$) and miosture ($g$ $kg^{-1}$ $d^{-1}$) forcing"
    #axx.text(80,14,text1,fontsize=charsize) 
#    plt.title(text1,fontsize=charsize) 
    plt.ylabel('Height'+r' ($km$)', fontsize=20) #fontdict=font)
#    plt.xlabel(yearstr, fontdict=font)
    plt.show()
    ij=ij+1
    jc=jc+1
cax = fig.add_axes([0.2, 0.045, 0.6, 0.03])
fig.colorbar(ax[0,0], cax,extend='both',
             spacing='uniform', orientation='horizontal')
fig.subplots_adjust(left=0.1,bottom=0.15,right=1-0.05,top=1-0.1,hspace=0.4)
plt.show()                     
plt.savefig(pic_out+'ALL_lsforcing_input.png',dpi=300)        
plt.show()
plt.close()
