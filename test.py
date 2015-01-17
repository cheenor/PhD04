# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 15:32:04 2014

@author: Chenjh
"""
import datetime
from matplotlib.dates import drange
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import calendar
from time import *
import string

a=10
b="%02d"%(a)
print b
b=1998
c=str(b)
print c[2:]
c=datetime.datetime(1990,3,1,0,0,0)
d=datetime.datetime(1990,3,31,0,0,0)
e=datetime.timedelta(hours=6)
dates=drange(c,d,e)
print c+2*e
c=calendar.monthrange(1997,7)
print calendar.monthrange(1997,7)[1]
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
dirp='Z:/DATA/LargeScale/TP/NcepR2_Pre/000101-001231/input/'
f=open(dirp+'ETP01.43')
ff=f.readlines()
linesplit=[]
for line in ff:
    line=string.lstrip(line)   ### delete the left ''
    linesplit.append(line[:-1].split(' '))
ff=f.readlines()
onedim=[]
for lnstrs in linesplit:
        for strs in lnstrs:
            if strs!='':
               onedim.append(string.atof(strs))  
ndim=52
dnm=31
nv=3
xdat=range(0,dnm*4)
dat=np.ndarray(shape=(ndim,dnm*4,nv), dtype=float) 
ll=0                                           
for it in range(0,dnm*4):
    for inv in range(0,nv):
        for ik in range(0,ndim):                        
            ll=it*(nv*ndim+1)+inv*ndim+ik+1
            dat[ik,it,inv]=onedim[ll] 
fig,axes=plt.subplots(nrows=nv,ncols=1,figsize=(15,4*nv))
#fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(6,6))
nv=3
for i in range(0,nv):                 
    zdat=dat[:,:,i]
    axes[i]=plt.subplot(nv,1,i)
#    matplotlib.scale.ScaleBase.limit_range_for_scale(0,100,0)
#    matplotlib.rcParams['contour.negative_colors'] = 'green'     
    axes[i]=plt.contour(xdat,ydat,zdat,10,cmap='seismic',
        linestyles='solid')
    plt.axis([0, dnm*4, 3.5, 20])
    plt.clabel(axes[i],inline=1,fontsize=10)
#    matplotlib.scale.ScaleBase.limit_range_for_scale(-100,0,0)
    axes[i]=plt.contour(xdat,ydat,zdat,10,
        cmap='seismic',linestyles='dotted')
    plt.clabel(axes[i],inline=1,fontsize=10)
    plt.axis([0, dnm*4, 3.5, 20])
    plt.show()                  
fig.subplots_adjust(hspace=0.4)
plt.savefig('D:/MyPaper/PhD04/Pics/test.png')
plt.close()


            