#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 02 20:02:42 2015

@author: chenjh
"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import string
from matplotlib.font_manager import FontProperties
nz=34
ng=5
CASENMSTR=['PRDCTR_EC','MLYRCTR_EC', 'NPCCTR_EC',
           'NECCTR_EC', 'WTPCTR_EC','ETPCTR_EC']   
orderstr=[r'($a$)',r'($b$)',r'($c$)',r'($d$)',r'($e$)',r'($f$)']
DATESTR  =['20120401' , '20100602' , '20100802' ,
           '20120706' , '20100703' , '20100603' ]
nga=len(CASENMSTR)
dirin="D:/MyPaper/PhD04/Cases/postdata/"
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
#--------------------read data from cloudcell.f ---------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
#plot  ----------  ft  # xdat,ydat,zdat
font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 16,
        }  
cloudlevs=[5,10,15,20,30,40,50,60,70,80,90,100,110,120,130]
cloudclors=['w','lightgray','plum','darkorchid','darkviolet','b','dodgerblue','skyblue','aqua',
            'greenyellow','lime','yellow','darkorange','chocolate','tomato','r']
fig,ax=plt.subplots(nrows=2,ncols=3,figsize=(12,12))

jr=0
jc=0
ij=1
for iga in range(0,nga):
    casenm=CASENMSTR[iga]
    if casenm[0:3]=='MLY' :
        areastr=casenm[0:4]
    else:
        areastr=casenm[0:3]
    if jc==3:
        jc=0
        jr=jr+1
    ft=np.ndarray(shape=(nz,nz,ng),dtype=float)
    fpath=dirin+casenm+'_ALLCLOUDCELSS_FREQUENCY_f90.TXT'
    for i in range(0,ng):
        iskp=i*(nz+2)+1
        nrl=nz+iskp
        onedim=readAscii(fpath,iskp,nrl)
        for ke in range(0,nz):
            for kb in range(0,nz):
                k=ke*(nz+1)+kb+1
                ft[kb,ke,i]=onedim[k]
    plt.subplot(2,3,ij)
#zdat[0,:]=0.0   ## the first level is below surface ground
    ft0=np.ndarray(shape=(nz,nz), dtype=float) #(km,km)  For exchange the dims
    for i1 in range(0,nz):
       i10=nz-i1-1
       for i2 in range(0,nz):
           ft0[i10,i2]=ft[i2,i1,4]
    ax[jr,jc]=plt.contourf(zdat,zdat,ft0,cmap=cm.binary, levels=cloudlevs,extend='both')
    marknm=orderstr[iga]+' '+areastr
    plt.title(marknm,fontsize=12)    
    plt.axis([0, 16, 0, 16])
    if jr==1:
        plt.xlabel(r'Cloud Base Height ($km$)', fontdict=font)
    if jc==0:
        plt.ylabel(r'Cloud Top Height ($km$)', fontdict=font)
    jc=jc+1
    ij=ij+1
figtitle = r"Frequency of all cloud cells ($10^{-2}$ %)"
fig.text(0.5, 0.95, figtitle,
    horizontalalignment='center',
    fontproperties=FontProperties(size=18))
#cax = fig.add_axes([0.2, 0.025, 0.6, 0.02])
cax = fig.add_axes([0.91, 0.2, 0.03, 0.6])
fig.colorbar(ax[0,0],cax,orientation='vertical',extend='both', # 'horizontal'
    extendfrac='auto',  spacing='uniform')
plt.show()
plt.savefig(dirpic+"AllCase_CloudCellsTopBase_fortran_grey.png",dpi=300)          
plt.show()
plt.close()
###############################################################################
font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 16,
        }  
cloudlevs=[5,10,15,20,30,40,50,60,70,80,90,100,110,120,130]
cloudclors=['w','lightgray','plum','darkorchid','darkviolet','b','dodgerblue','skyblue','aqua',
            'greenyellow','lime','yellow','darkorange','chocolate','tomato','r']
fig,ax=plt.subplots(nrows=2,ncols=3,figsize=(12,12))
jr=0
jc=0
ij=1
for iga in range(0,nga):
    casenm=CASENMSTR[iga]
    if casenm[0:3]=='MLY' :
        areastr=casenm[0:4]
    else:
        areastr=casenm[0:3]
    if jc==3:
        jc=0
        jr=jr+1
    ft=np.ndarray(shape=(nz,nz,ng),dtype=float)
    fpath=dirin+casenm+'_ALLCLOUDCELSS_FREQUENCY_f90.TXT'
    for i in range(0,ng):
        iskp=i*(nz+2)+1
        nrl=nz+iskp
        onedim=readAscii(fpath,iskp,nrl)
        for ke in range(0,nz):
            for kb in range(0,nz):
                k=ke*(nz+1)+kb+1
                ft[kb,ke,i]=onedim[k]
    plt.subplot(2,3,ij)
#zdat[0,:]=0.0   ## the first level is below surface ground
    ft0=np.ndarray(shape=(nz,nz), dtype=float) #(km,km)  For exchange the dims
    for i1 in range(0,nz):
       i10=nz-i1-1
       for i2 in range(0,nz):
           ft0[i10,i2]=ft[i2,i1,4]
    ax[jr,jc]=plt.contourf(zdat,zdat,ft0,cmap=cm.binary, levels=cloudlevs,extend='both')
    marknm=orderstr[iga]+' '+areastr
    plt.title(marknm,fontsize=12)    
    plt.axis([0, 16, 0, 16])
    if jr==1:
        plt.xlabel(r'Cloud Base Height ($km$)', fontdict=font)
    if jc==0:
        plt.ylabel(r'Cloud Top Height ($km$)', fontdict=font)
    jc=jc+1
    ij=ij+1
figtitle = r"Frequency of all cloud cells ($10^{-2}$ %)"
fig.text(0.5, 0.95, figtitle,
    horizontalalignment='center',
    fontproperties=FontProperties(size=18))
#cax = fig.add_axes([0.2, 0.025, 0.6, 0.02])
cax = fig.add_axes([0.91, 0.2, 0.03, 0.6])
fig.colorbar(ax[0,0],cax,orientation='vertical',extend='both', # 'horizontal'
    extendfrac='auto',  spacing='uniform')
plt.show()
plt.savefig(dirpic+"AllCase_CloudCellsTopBase_fortran_grey.png",dpi=300)          
plt.show()
plt.close()  