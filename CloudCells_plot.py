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
CASENMSTR=['PRDCTR_EC','MLYRCTR_EC', 'NPCCTR_EC',
           'NECCTR_EC','WTPCTR_EC' , 'ETPCTR_EC']   
astr=[r'$(a)$',r'$(b)$', r'$(c)$',r'$(d)$',r'$(e)$',r'$(f)$'] 
DATESTR  =['20120401' , '20100602' , '20100802' ,
           '20120706' , '20100703' , '20100603' ]
dirin="D:/MyPaper/PhD04/Cases/postdata/"
dirpic="D:/MyPaper/PhD04/Pics/"
#
diro='D:/MyPaper/PhD04/Cases/ERA/FORCING/'
dis_pressure_ea=[750.,550.,400.,250.,150.]
dis_pressure_tp=[450.,300.,200.,100.]
nzz=52
ntx=121
f52
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
def readAscii2(fpath,iskp):
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
def pressure2heigh(pres,tmp4prs,ydat,prelevel):
    nt=len(tmp4prs[0,:])
    nz=len(tmp4prs[:,0])
    menatmp=np.zeros(shape=(nz),dtype=float)
    for iz in range(0,nz):
        for it in range(0,nt):
            menatmp[iz]=menatmp[iz]+tmp4prs[iz,it]/nt    
    for iz in range(1,nz):
        if pres<prelevel[iz-1] and pres>=prelevel[iz]:
            z1=ydat[iz-1]*1000. # km to m
            p1=prelevel[iz-1]
            at=((menatmp[iz]+menatmp[iz-1])/2.-273.15)*1./273.
            z2=z1+18400*(1+at)*math.log10(p1/pres)
            return z2
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
fig,ax=plt.subplots(nrows=2,ncols=3,figsize=(15,12))
ir=0
jc=0
ij=1
filestring='ALLCLOUDCELSS'
for iga in range(0,nga):
    casenm=CASENMSTR[iga]
    dis_pressure=dis_pressure_ea
    if iga ==4 or iga==5:
        dis_pressure=dis_pressure_tp
    if jc==3:
        jc=0
        ir=ir+1
    print ir,jc
    if casenm[0:3]=='MLY':
        area=casenm[0:4]
    else:
        area=casenm[0:3]
    f52=area+'_'+DATESTR[iga]+"_031d_ERA_52pressure.52"
    dirobs=diro+area+'/'
    fpath=dirobs+f52
    iskp=0
    prelevel=readAscii2(fpath, iskp)
    tmp4prs=np.zeros(shape=(nzz,ntx),dtype=float)       
    fpath=dirobs+'temperature.txt'
    onedim1=readAscii2(fpath, 0)   
    for it in range(0,ntx):
        for iz in range(0,nzz):
            k=it*nzz+iz
            #print k,len(onedim1)              
            tmp4prs[iz,it]=onedim1[k]
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
        'size'   : 18,
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
    plt.title(tilstr, fontsize=20)
    if jc==0:
        plt.ylabel(r'Cloud Top Height ($km$)', fontdict=font)
    if ir==1:            
        plt.xlabel(r'Cloud Base Height ($km$)', fontdict=font)
    axx= plt.subplot(2,3,ij)
    ymajorLocator   = MultipleLocator(4) 
    axx.yaxis.set_major_locator(ymajorLocator)
    xmajorLocator   = MultipleLocator(4) 
    axx.xaxis.set_major_locator(xmajorLocator)
    #if iga in(1,2,4,5)  :
    #    for tick in axx.yaxis.get_major_ticks():
    #        tick.label1On = False
    #if ir==0  :
    #    for tick in axx.xaxis.get_major_ticks():
    #        tick.label1On = False
    presstr='%d'%prelevel[0]
    axx.text(16.5,0,presstr, fontdict=font)
    #axx.text(0,15,presstr, fontdict=font,rotation=-90)
    for pres in dis_pressure:    
        hp=pressure2heigh(pres,tmp4prs,ydat,prelevel)
        print hp
        hp=hp/1000.
        presstr='%d'%pres
        axx.text(16.5,hp,presstr, fontdict=font)
        #axx.text(hp,15,presstr, fontdict=font,rotation=-90)
    #axx.text(18,12,'Pressure '+r'($hPa$)', fontdict=font,rotation=-90)
    if ij in(3,6):
        axx.text(19.5,11,'Pressure '+r'($hPa$)', fontdict=font,rotation=-90)
#        plt.yticks() #
    jc=jc+1
    ij=ij+1
plt.subplots_adjust(left = 0.1, wspace = 0.4, hspace = 0.25, \
    bottom = 0.20, top = 0.90)
cax = fig.add_axes([0.1, 0.06, 0.8, 0.04])
fig.colorbar(ax[0,0], cax,extend='both',
             spacing='uniform', orientation='horizontal')                                                     
titlename=r"Frequency of all cloud cells ($10^{-2}%$)"
plt.title(titlename,fontsize=16)
plt.show()
plt.savefig(dirpic+"AllCases_CloudCellsTopBase_fortran_color_p.png",dpi=300)          
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
    greys=np.linspace(0.0,0.8,17)
    ncl=len(greys)
    mycolor=['1.0']
    for icl in range(0,ncl):
        mycolor.append('%.2f'%greys[ncl-icl-1])       
    ax[ir,jc]=plt.contourf(zdat,zdat,ft0,colors=mycolor, levels=cloudlevs,extend='both')
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
plt.savefig(dirpic+"AllCases_CloudCellsTopBase_fortran_Gray_New.png",dpi=300)          
plt.show()
plt.close()
##############################################################################
fig,ax=plt.subplots(nrows=2,ncols=3,figsize=(15,12))
ir=0
jc=0
ij=1
filestring='ALLCLOUDCELSS'
for iga in range(0,nga):
    casenm=CASENMSTR[iga]
    dis_pressure=dis_pressure_ea
    if iga ==4 or iga==5:
        dis_pressure=dis_pressure_tp
    if jc==3:
        jc=0
        ir=ir+1
    print ir,jc
    if casenm[0:3]=='MLY':
        area=casenm[0:4]
    else:
        area=casenm[0:3]
    f52=area+'_'+DATESTR[iga]+"_031d_ERA_52pressure.52"
    dirobs=diro+area+'/'
    fpath=dirobs+f52
    iskp=0
    prelevel=readAscii2(fpath, iskp)
    tmp4prs=np.zeros(shape=(nzz,ntx),dtype=float)       
    fpath=dirobs+'temperature.txt'
    onedim1=readAscii2(fpath, 0)   
    for it in range(0,ntx):
        for iz in range(0,nzz):
            k=it*nzz+iz
            #print k,len(onedim1)              
            tmp4prs[iz,it]=onedim1[k]
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
        'size'   : 18,
        }  
#zdat[0,:]=0.0   ## the first level is below surface ground
    ft0=np.ndarray(shape=(nz,nz), dtype=float) #(km,km)  For exchange the dims
    for i1 in range(0,nz):
        i10=nz-i1-1
        for i2 in range(0,nz):
            ft0[i10,i2]=ft[i2,i1,4]
    greys=np.linspace(0.0,0.86,17)
    ncl=len(greys)
    mycolor=['1.0']
    for icl in range(0,ncl):
        mycolor.append('%.2f'%greys[ncl-icl-1]) 
    plt.subplot(2,3,ij)        
    ax[ir,jc]=plt.contourf(zdat,zdat,ft0,colors=mycolor, levels=cloudlevs,extend='both')
    plt.axis([0, 16, 0, 16])
    tilstr=astr[iga]+' '+area
    plt.title(tilstr, fontsize=20)
    if jc==0:
        plt.ylabel(r'Cloud Top Height ($km$)', fontdict=font)
    if ir==1:            
        plt.xlabel(r'Cloud Base Height ($km$)', fontdict=font)
    axx= plt.subplot(2,3,ij)
    ymajorLocator   = MultipleLocator(4) 
    axx.yaxis.set_major_locator(ymajorLocator)
    xmajorLocator   = MultipleLocator(4) 
    axx.xaxis.set_major_locator(xmajorLocator)
    #if iga in(1,2,4,5)  :
    #    for tick in axx.yaxis.get_major_ticks():
    #        tick.label1On = False
    #if ir==0  :
    #    for tick in axx.xaxis.get_major_ticks():
    #        tick.label1On = False
    presstr='%d'%prelevel[0]
    axx.text(16.5,0,presstr, fontdict=font)
    #axx.text(0,15,presstr, fontdict=font,rotation=-90)
    for pres in dis_pressure:    
        hp=pressure2heigh(pres,tmp4prs,ydat,prelevel)
        print hp
        hp=hp/1000.
        presstr='%d'%pres
        axx.text(16.5,hp,presstr, fontdict=font)
        #axx.text(hp,15,presstr, fontdict=font,rotation=-90)
    #axx.text(18,12,'Pressure '+r'($hPa$)', fontdict=font,rotation=-90)
    if ij in(3,6):
        axx.text(19.5,11,'Pressure '+r'($hPa$)', fontdict=font,rotation=-90)
#        plt.yticks() #
    jc=jc+1
    ij=ij+1
plt.subplots_adjust(left = 0.1, wspace = 0.4, hspace = 0.25, \
    bottom = 0.20, top = 0.90)
cax = fig.add_axes([0.1, 0.06, 0.8, 0.04])
fig.colorbar(ax[0,0], cax,extend='both',
             spacing='uniform', orientation='horizontal')                                                     
titlename=r"Frequency of all cloud cells ($10^{-2}%$)"
plt.title(titlename,fontsize=16)
plt.show()
plt.savefig(dirpic+"AllCases_CloudCellsTopBase_fortran_Grey_p.png",dpi=300)          
plt.show()
plt.close()   