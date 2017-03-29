#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 08 07:56:28 2015

@author: jhchen
"""
import matplotlib.pyplot as plt
import matplotlib.axes as mplaxes
import numpy as np
import string
from pylab import *
import matplotlib as mpl
nz=34
ng=5
CASENMSTR=['PRDCTR_H','MLYRCTR_H', 'NPCCTR_H',
           'NECCTR_H','WTPCTR_H' , 'ETPCTR_H']   
astr=[r'$(a)$',r'$(b)$', r'$(c)$',r'$(d)$',r'$(e)$',r'$(f)$'] 
DATESTR  =['20120401' , '20100602' , '20100802' ,
           '20120706' , '20100703' , '20100603' ]
dirin="Z:/CRM/500m//Postdata/"
dirpic="D:/MyPaper/PhD04/Pics/"
#-----------------------------------------------------------------------
zdat=[              0.0500000, 0.1643000, 0.3071000, 0.4786000
    , 0.6786000, 0.9071000, 1.1640000, 1.4500000, 1.7640001
    , 2.1070001, 2.4790001, 2.8789999, 3.3069999, 3.7639999
    , 4.2500000, 4.7639999, 5.3070002, 5.8790002, 6.4790001
    , 7.1069999, 7.7639999, 8.4499998, 9.1639996, 9.9069996
    ,10.6800003,11.4799995,12.3100004,13.1599998,14.0500002
    ,14.9600000,15.9099998,16.8799992,17.8799992,18.9099998]
#
def readAscii(fpath,iskp,nrl):
    #iskp  the total line skipped of the file
    # fpath   the full path of the file
    # usage: onedim=readAscii(fpaht,iskp)
    onedim=[]
    linesplit=[]
    f=open(fpath)
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
#--------------------read data from cloudcell.f ---------------------------
cdict = {'red': ((0., 1, 1),
                 (0.15, 0, 0),
                 (0.66, 1, 1),
                 (0.89, 1, 1),
                 (1, 0.5, 0.5)),
         'green': ((0., 1, 1),
                   (0.15, 0, 0),
                   (0.38, 1, 1),
                   (0.64, 1, 1),
                   (0.91, 0, 0),
                   (1, 0, 0)),
         'blue': ((0., 1, 1),
                  (0.15, 1, 1),
                  (0.35, 1, 1),
                  (0.7, 0, 0),
                  (1, 0, 0))}
my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
nga=len(CASENMSTR)
cloudlevs=[3,5,10,15,20,25,30,35,40,50,60,70,80,90,100,110,130]
cloudclors=['w','lightgray','plum','darkorchid','darkviolet','b','dodgerblue','skyblue','aqua',
            'greenyellow','lime','limegreen','yellow','darkorange','tomato','r']
fig,ax=plt.subplots(nrows=2,ncols=3,figsize=(12,12))
ir=0
jc=0
ij=1
filestring='ALLCLOUDCELSS'
for iga in range(0,nga):
    casenm=CASENMSTR[iga]
    if jc==3:
        jc=0
        ir=ir+1
    print ir,jc
    if casenm[0:3]=='MLY':
        area=casenm[0:4]
    else:
        area=casenm[0:3]
    ft=np.ndarray(shape=(nz,nz,ng),dtype=float)
    fpath=dirin+casenm+'_'+filestring+'_FREQUENCY_f90.TXT'
    for i in range(0,ng):
        iskp=i*(nz+2)+1
        nrl=nz+iskp
        onedim=readAscii(fpath,iskp,nrl)
        for ke in range(0,nz):
            for kb in range(0,nz):
                k=ke*(nz+1)+kb+1
                ft[kb,ke,i]=onedim[k]            
#plot  ----------  ft  # xdat,ydat,zdat
    font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 16,
        }  
#zdat[0,:]=0.0   ## the first level is below surface ground
    ft0=np.ndarray(shape=(nz,nz), dtype=float) #(km,km)  For exchange the dims
    for i1 in range(0,nz):
        i10=nz-i1-1
        for i2 in range(0,nz):
            ft0[i10,i2]=ft[i2,i1,4]
    plt.subplot(2,3,ij)        
    ax[ir,jc]=plt.contourf(zdat,zdat,ft0,cmap=my_cmap, levels=cloudlevs,extend='both')
    plt.axis([0, 16, 0, 16])
    tilstr=astr[iga]+' '+area
    plt.title(tilstr, fontsize=16)
    if jc==0:
        plt.ylabel(r'Cloud Top Height ($km$)', fontdict=font)
    if ir==1:            
        plt.xlabel(r'Cloud Base Height ($km$)', fontdict=font)
    if iga in(1,2,4,5)  :
        axx= plt.subplot(2,3,ij)
        for tick in axx.yaxis.get_major_ticks():
            tick.label1On = False
    if ir==0  :
        axx= plt.subplot(2,3,ij)
        for tick in axx.xaxis.get_major_ticks():
            tick.label1On = False
#        plt.yticks() #
    jc=jc+1
    ij=ij+1
plt.subplots_adjust(left = 0.1, wspace = 0.1, hspace = 0.2, \
    bottom = 0.25, top = 0.90)
cax = fig.add_axes([0.1, 0.08, 0.8, 0.04])
fig.colorbar(ax[0,0], cax,extend='both',
             spacing='uniform', orientation='horizontal')                                                     
titlename=r"Frequency of all cloud cells ($10^{-2}%$)"
plt.title(titlename,fontsize=16)
plt.show()
plt.savefig(dirpic+"AllCases_CloudCellsTopBase_fortran_color_500m.png",dpi=300)          
plt.show()
plt.close()
###############################################################################
nga=len(CASENMSTR)
cloudlevs=[3,5,10,15,20,25,30,35,40,50,60,70,80,90,100,110,130]
cloudclors=['w','lightgray','plum','darkorchid','darkviolet','b','dodgerblue','skyblue','aqua',
            'greenyellow','lime','limegreen','yellow','darkorange','tomato','r']
fig,ax=plt.subplots(nrows=2,ncols=3,figsize=(12,12))
ir=0
jc=0
ij=1
filestring='ALLCLOUDCELSS'
for iga in range(0,nga):
    casenm=CASENMSTR[iga]
    if jc==3:
        jc=0
        ir=ir+1
    print ir,jc
    if casenm[0:3]=='MLY':
        area=casenm[0:4]
    else:
        area=casenm[0:3]
    ft=np.ndarray(shape=(nz,nz,ng),dtype=float)
    fpath=dirin+casenm+'_'+filestring+'_FREQUENCY_f90.TXT'
    for i in range(0,ng):
        iskp=i*(nz+2)+1
        nrl=nz+iskp
        onedim=readAscii(fpath,iskp,nrl)
        for ke in range(0,nz):
            for kb in range(0,nz):
                k=ke*(nz+1)+kb+1
                ft[kb,ke,i]=onedim[k]            
#plot  ----------  ft  # xdat,ydat,zdat
    font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 16,
        }  
#zdat[0,:]=0.0   ## the first level is below surface ground
    ft0=np.ndarray(shape=(nz,nz), dtype=float) #(km,km)  For exchange the dims
    for i1 in range(0,nz):
        i10=nz-i1-1
        for i2 in range(0,nz):
            ft0[i10,i2]=ft[i2,i1,4]
    plt.subplot(2,3,ij)        
    ax[ir,jc]=plt.contourf(zdat,zdat,ft0,cmap=cm.binary, levels=cloudlevs,extend='both')
    plt.axis([0, 16, 0, 16])
    tilstr=astr[iga]+' '+area
    plt.title(tilstr, fontsize=16)
    if jc==0:
        plt.ylabel(r'Cloud Top Height ($km$)', fontdict=font)
    if ir==1:            
        plt.xlabel(r'Cloud Base Height ($km$)', fontdict=font)
    if iga in(1,2,4,5)  :
        axx= plt.subplot(2,3,ij)
        for tick in axx.yaxis.get_major_ticks():
            tick.label1On = False
    if ir==0  :
        axx= plt.subplot(2,3,ij)
        for tick in axx.xaxis.get_major_ticks():
            tick.label1On = False
#        plt.yticks() #
    jc=jc+1
    ij=ij+1
plt.subplots_adjust(left = 0.1, wspace = 0.1, hspace = 0.2, \
    bottom = 0.25, top = 0.90)
cax = fig.add_axes([0.1, 0.08, 0.8, 0.04])
fig.colorbar(ax[0,0], cax,extend='both',
             spacing='uniform', orientation='horizontal')                                                     
titlename=r"Frequency of all cloud cells ($10^{-2}%$)"
plt.title(titlename,fontsize=16)
plt.show()
plt.savefig(dirpic+"AllCases_CloudCellsTopBase_fortran_Gray_500m.png",dpi=300)          
plt.show()
plt.close()   