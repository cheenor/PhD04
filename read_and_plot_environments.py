#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
read and plot environments of the CRM input data

@author: Chenjh
"""
import matplotlib
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
rgns=['ETP','WTP']
fnm=['OBS_Surface_input.txt','Q1Q2_profiles_input.txt','.43','.49','.dyn',
     '_lsforcing.37','_surface.39','_thetaqv_profile.41','_uv_profiles.35',
     '.99']
nv=[3,2,3,3,5,2,5,2,2,2]    # number variables of every file
ndim=[1,17,52,1,52,52,1,52,52,52]  # is the varlables has the vertical dimension
dirin='Z:/DATA/LargeScale/TP/NcepR2_Pre/'
fold2='/input/'
pic_out='D:/MyPaper/PhD04/RawPics_MJJAS/'
det=datetime.timedelta(hours=6)
lev99=[-8,-4,-2,2,6,9]
color99=['g','g','g','r','r','r']
##################################################
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['contour.negative_linestyle'] = 'dashed'
matplotlib.rcParams['savefig.dpi'] = 100
plt.rc('lines', linewidth=4)
###################################################
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
ylevs=[1000, 925, 850, 700, 600, 500, 400, 300, 250, 
       200, 150, 100, 70, 50, 30, 20, 10]
for yd in ydat_r:
    ydat.append(yd*0.001)
del ydat_r    
###########################################################################
for iy in range(1990,2013):
    year=str(iy)  # number to string
    fold1=year[2:]+'0101'+'-'+year[2:]+'1231'  #### connect the string
    dirpath=dirin+fold1+fold2
    
    for rg in rgns:
        for im in range(4,9):
            monstr="%02d"%(im+1)  ### number to string, 1 to 01, 10 to 10
            dnm=calendar.monthrange(iy,im+1)[1]  #calendar.monthrange(1997,7) 
#                              #reture two index, sencond is the day number
            datestart=datetime.datetime(iy,im+1,1,0,0,0)            
            dateiso=[]            
            for dt in range(0,dnm*4):
                dateiso.append(datestart+dt*det)
            xdate=[]    
            xdat=range(0,dnm*4)            
            for tm in dateiso:
                xdate.append(datetime.datetime.strftime( tm,"%b/%d"))  
            del dateiso    
            nf=len(fnm)
            for n in range(nf-1,nf): #range(0,nf):
                onedim=[]
                linesplit=[]
                nsrts=0                
                if n<2:
                  nstrs=1  
                  filepath=dirpath+year+rg+'_'+monstr+fnm[n]
                else:
                  filepath=dirpath+rg+monstr+fnm[n]  
                f=open(filepath)
                if n>1:
                    nstrs=0
                ff=f.readlines()[nstrs:]
                for line in ff:
                    line=string.lstrip(line)
                    linesplit.append(line[:-1].split(' '))
                for lnstrs in linesplit:
                    for strs in lnstrs:
                        if strs!='':
                            onedim.append(string.atof(strs))
#                linesplit 
                dat=np.ndarray(shape=(ndim[n],dnm*4,nv[n]), dtype=float) 
                ll=0                                           
                for it in range(0,dnm*4):
                    for inv in range(0,nv[n]):
                        for ik in range(0,ndim[n]):                        
#                            ll=it*nv[n]+inv*ndim[n]+1+ik
                            ll=it*(nv[n]*ndim[n]+1)+inv*ndim[n]+ik+1
                            dat[ik,it,inv]=onedim[ll]
###################################################################### 
                if ndim[n]<2:                    
                    fig=plt.figure(figsize=(8,6))  ### all in one panle
                    xlend=[]
                else:
                    fig,axes=plt.subplots(nrows=nv[n],ncols=1,figsize=(15,4*nv[n]))
#                    fig = plt.subplots(nrows=nv[n])                            
                for i in range(0,nv[n]):
                    if ndim[n]<2:                    
                        zdat=dat[0,:,i]                        
                    else:
                        zdat=dat[:,:,i]
                    if ndim[n] <2: 
                        plt.plot(xdat, zdat)
                        axe=fig.add_subplot(111)
                        axe.set_xticks(range(0,dnm*4,16))
                        xticklabels = [xdate[nn] for nn in range(0,dnm*4,16)]
                        axe.set_xticklabels(xticklabels, size=10)
                    else:                                                    
                        axes[i]=plt.subplot(nv[n],1,i)
                        if ndim[n]==17:
                            ydatx=ylevs
                        else:
                            ydatx=ydat    
                        if zdat.any>0:
                            matplotlib.rcParams['contour.negative_linestyle'] = 'dashed'
                            axes[i]=plt.contour(xdat,ydatx,zdat,colors=color99,
                                    linewidths=1,levels=lev99)
                            if len(ydatx)==17:
                                plt.axis([0, dnm*4, 600,50])
                            else:
                                plt.axis([0, dnm*4, 3.5, 18])
                            plt.clabel(axes[i],inline=1,fontsize=10) 
                            plt.show()                        
#                        if zdat.any<0:
#                            axes[i]=plt.contour(xdat,ydatx,zdat,linewidths=1,
#                                    colors='green',linestyles='dotted')
#                            if len(ydatx)==17:
#                                plt.axis([0, dnm*4,600,50])
#                            else:
#                                plt.axis([0, dnm*4, 3.5, 18])                                  
#                            plt.clabel(axes[i],inline=1,fontsize=10) 
#                            plt.show()
                        axx=fig.add_subplot(nv[n],1,i)                         
                        axx.set_xticks(range(0,dnm*4,16))
                        xticklabels = [xdate[nn] for nn in range(0,dnm*4,16)] 
                        axx.set_xticklabels(xticklabels, size=14)
                        plt.show()                     
                plt.savefig(pic_out+rg+year+monstr+fnm[n]+'.png')          
                plt.show()
                plt.close()
                del zdat,dat,onedim,linesplit
                f.close()