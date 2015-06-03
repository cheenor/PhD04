#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 11:06:37 2015

@author: jhchen
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import string
import numpy as np
import datetime
import os
import calendar
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20

dirint="D:/MyPaper/PhD04/Cases/AllMonths/"
dirdata="X:/Data/ERA_interim/"
rgname=["ETP","WTP","PRD","MLYR","NPC","NEC"]
nas=len(rgname)


###############################################################################
varnames=["Temp. Forc", "Moisture Forcing", "U wind"         ,"V wind", "Vapor",
          "Height"    , "Temperature"     , "Adjusted Omega" , "Theta",
          "Q1"        , "Q2"              , "HADQ"           , "VADQ" ,
          "TCHQ"      , "HADT"            , "VADT"           , "TCHT" ,
          "Original Omega"]
varunits=[r"$K$ $day^{-1}$", r"$K$ $day^{-1}$", r"$m$ $s^{-1}$"   ,   r"$m$ $s^{-1}$", r"$g$ $kg^{-1}$",
          r"$m$"           ,         r"$K$"   , r"$hPa$ $s^{-1}$" ,   r"$K$"         ,
          r"$K$ $day^{-1}$", r"$K$ $day^{-1}$", r"$K$ $day^{-1}$" , r"$K$ $day^{-1}$" ,
          r"$K$ $day^{-1}$", r"$K$ $day^{-1}$", r"$K$ $day^{-1}$", r"$K$ $day^{-1}$" ,
          r"$hPa$ $s^{-1}$"]
nvs=len(varnames)
iwant=[0,1]
wantednm=[]
for i in iwant:
    wantednm.append(iwant[i])
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
ydat_p=[1000,	975,	950,	925,	900,	875,	850,	825,	800,	775,	750,	700,	650,
      600,	550,	500,	450,	400,	350,	300,	250,	225,	200,	175,	150,	125,
      100,	70,	50,	30,	20,	10,	7,	5,	3,	2,	1]
for yd in ydat_p:
    ydat.append(yd*1.0)
nz=len(ydat)
###############################################################################
def readAscii(fpath,iskp,nvs,nt,nz):
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
    raw=np.ndarray(shape=(nvs,nz,nt),dtype=float)
    print len(onedim),nt*nz*(nvs+1)
    for it in range(0,nt):
        for iz in range(0,nz):
            for iv in range(0,nvs):
                k=it*nz*(nvs+1)+iz*(nvs+1)+iv+1
#                print iv,iz,it,k
                raw[iv,iz,it]=onedim[k]
    del onedim    
    return raw
def getdatestr(iyr,im,jd,nt):
    datestart=datetime.datetime(iyr,im+1,jd,0,0,0)
    det=datetime.timedelta(hours=6)            
    dateiso=[]            
    for dt in range(0,nt):
        dateiso.append(datestart+dt*det)
    xdate=[]    
    xxx=range(0,nt)            
    for tm in dateiso:
        xdate.append(datetime.datetime.strftime(tm,"%b/%d")) 
    return xxx, xdate
def plotting(raw,its,ite,xxx,xdate,ndt,ydat,iwant,varnames,varunits,pic_out):
    levs1=[-15,-9,-6,-3,3,6,9,15]
    colors1=['g','g','g','g','r','r','r','r']
    linetype1=['dotted','dotted','dotted','dotted','solid','solid','solid','solid'] 
    levs2=[-15,-9,-6,-3,3,6,9,15]
    colors2=['g','g','g','g','r','r','r','r']
    linetype2=['dotted','dotted','dotted','dotted','solid','solid','solid','solid']
    nl=len(pic_out)
    titlename=[]
    for i in iwant:
        titlename.append(varnames[i]+" ("+varunits[i]+")  "+pic_out[nl-6:] )
    font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 16,
        }
    charsize=18 
    if len(iwant)>1 :
        fig,[axe1,axe2]=plt.subplots(nrows=2,ncols=1,figsize=(15,9))
        plt.subplot(2,1,1)
        i=iwant[0]
        axe1=plt.contour(xxx,ydat,raw[i,:,its:ite],colors=colors1,
            linewidths=1.5,levels=levs1,linestyles=linetype1)                           
#plt.title('Observation',fontsize=charsize)                       
        plt.axis([0, ndt, 1000, 50])  ## x axis  y axis
        plt.clabel(axe1,inline=1,fmt='%1d',fontsize=charsize-2)                                                 
        axx=fig.add_subplot(2,1,1)                         
        axx.set_xticks(range(0,ndt,16))
        xticklabels = [xdate[nn] for nn in range(0,ndt,16)] 
        axx.set_xticklabels(xticklabels, size=charsize)
        text1=r"($a$)"
        axx.text(1.5,14,text1,fontsize=charsize+4) 
        text1=titlename[i]
        axx.text(80,14,text1,fontsize=charsize)  
        plt.ylabel('Pressure'+r' ($hPa$)', fontdict=font)
        plt.show() 
###
        plt.subplot(2,1,2)
        i=iwant[1]
        axe2=plt.contour(xxx,ydat,raw[i,:,its:ite],colors=colors2,
            linewidths=1.5,levels=levs2,linestyles=linetype2)                           
#plt.title(casenm,fontsize=charsize)                        
        plt.axis([0, ndt, 1000, 50])  ## x axis  y axis
        plt.clabel(axe2,inline=1,fmt='%1d',fontsize=charsize-2)                       
        axx=fig.add_subplot(2,1,2)                         
        axx.set_xticks(range(0,ndt,16))
        xticklabels = [xdate[nn] for nn in range(0,ndt,16)] 
        axx.set_xticklabels(xticklabels, size=charsize)
        text1=r"($b$)"
        axx.text(1.5,14,text1,fontsize=charsize+4) 
        text1=titlename[i]
        axx.text(80,14,text1,fontsize=charsize) 
        plt.ylabel('Pressure'+r' ($hPa$)', fontdict=font)
        plt.show()                     
        plt.savefig(pic_out+'_lsforcing_input.png')          
        plt.show()
        plt.close()
############################################################################### 
ndays=[31,28,31,30,31,30,31,31,30,31,30,31]
for iyr in range(1979,2013):
    yearstr="%d"%iyr
    leny=len(yearstr)
    fold1=yearstr[leny-2:]+'0101-'+yearstr[leny-2:]+'1231/'
    nt=365*4
    ndays[1]=28
    if calendar.isleap(iyr):
        nt=366*4
        ndays[1]=29
    for ig in range(0,nas):
        if ig <2 :
            fold0='ERA_PRE_O/'
        else:
            fold0='ERA_EA/'        
        filename=yearstr+rgname[ig]+'_RAW.txt' 
        fpath=dirdata+fold0+fold1+filename
        iskp=1
        raw=readAscii(fpath,iskp,nvs,nt,nz)
        for im in range(3,9):
            monstr="%2.2d"%(im+1)
            its=0
            for imm in range(0,im):
                its=its+ndays[imm]*4
            ndd=ndays[im] # day number of m month of year iy
            ite=its+ndd*4
            ndt=ndd*4
            jd=1
            xxx,xdate=getdatestr(iyr,im,jd,ndt)
            #
            dirout=dirint+rgname[ig]
            if os.path.exists(dirout):
                pic_out=dirout+"/"+rgname[ig]+"-"+yearstr+monstr
            else:
                os.makedirs(dirout)
                pic_out=dirout+"/"+rgname[ig]+"-"+yearstr+monstr
            plotting(raw,its,ite,xxx,xdate,ndt,ydat,iwant,varnames,varunits,pic_out)
            #
#
#
#