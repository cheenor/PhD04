#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 02 00:10:35 2015

@author: jhchen
"""
import string
import numpy as np
dirin=""
dirout=""
nt=2880
nx=200
nz=52
bin0=1.
abin=[]
for i in range(1,1000):
    abin.append(i*0.1+bin0)
nbin=len(abin)
rgname=["ETP","WTP","PRD","MLYR","NPC","NEC"]
nrg=len(rgname)

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
###############################################################################
def cdf(rain,abin,train,nbin,nt,nx):
    cdfr=np.ndarray(shape=(nbin),dtype=float)
    for i in range(0,nbin):
        bin0=abin[i]
        cont=0.0
        for it in range(0,nt):
            for ix in range(0,nx):
                if rain[ix,it] >= bin0:
                    cont=cont+1.0
        cdfr[i]=cont/train
    return cdfr
###############################################################################
cdf_reu=np.ndarray(shape=(nbin,nrg),dtype=float)
for i in range(0,nrg):    
    fpath=dirin+rgname[i]
    iskp=0
    onedim=readAscii(fpath,iskp)
    rain=np.ndarray(shape=(nx,nt),dtype=float)
    train_bin0=0.
    for it in range(0,nt):
        for ix in range(0,nx):
            k=it*nx+ix
            rain[ix,it]=onedim[k]
            if rain[ix,it] >= bin0:
                train_bin0=train_bin0+1.0
    cdf_reu[:,i]=cdf(rain,abin,train_bin0,nbin,nt,nx) 
###############################################################################
#  Plotting for cdf_reu 
