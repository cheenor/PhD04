#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 02 20:02:42 2015

@author: chenjh
"""
import matplotlib.pyplot as plt
import numpy as np
import string
nz=34
ng=5
CASENMSTR=['PRDCTR_EC','MLYRCTR_EC', 'NPCCTR_EC',
           'NECCTR_EC','ETPCTR_EC' , 'WTPCTR_EC']   
DATESTR  =['20120401' , '20100624' , '20100802' ,
           '20120706' , '20120520' , '20100714' ]
dirin="D:/MyPaper/PhD04/Cases/postdata/"
dirpic="D:/MyPaper/PhD04/Pics/"
casenm=CASENMSTR[4]
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
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
#plot  ----------  ft  # xdat,ydat,zdat
font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 18,
        }  
cloudlevs=[2,5,10,15,20,30,40,50,60,70,80,90,100,110,130]
cloudclors=['w','lightgray','plum','darkorchid','darkviolet','b','dodgerblue','skyblue','aqua',
            'greenyellow','lime','yellow','darkorange','chocolate','tomato','r']
fig,axe1=plt.subplots(nrows=1,ncols=1,figsize=(6,6))
plt.subplot(1,1,1)
#zdat[0,:]=0.0   ## the first level is below surface ground
ft0=np.ndarray(shape=(nz,nz), dtype=float) #(km,km)  For exchange the dims
for i1 in range(0,nz):
    i10=nz-i1-1
    for i2 in range(0,nz):
        ft0[i10,i2]=ft[i2,i1,4]
titlename=r"Frequency of all cloud cells ($10^{-2}%$)"
axe1=plt.contourf(zdat,zdat,ft0,colors=cloudclors, levels=cloudlevs,extend='both')
plt.colorbar(orientation='horizontal',extend='both',
    extendfrac='auto',  spacing='uniform')                           
plt.title(titlename,fontsize=16)                          
plt.axis([0, 16, 0, 16])
plt.xlabel(r'Cloud Base Height ($km$)', fontdict=font)
plt.ylabel(r'Cloud Top Height ($km$)', fontdict=font)
plt.show()
plt.savefig(dirpic+casenm+"_CloudCellsTopBase_fortran.pdf")          
plt.show()
plt.close()   