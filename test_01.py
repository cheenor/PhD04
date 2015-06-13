# -*- coding: utf-8 -*-
"""
Created on Wed May 27 13:17:30 2015

@author: jhchen
"""
from pytrmm import TRMM3B40RTFile as TRM3B40
import os
import numpy as np
import matplotlib.pyplot as plt
"""
os.system("cls")
dirin='E:/Data/TRMM/3b40/2010/06/'
fpath=dirin+"3B40RT.2010062718.7R2.bin.gz"
trmm_file = TRM3B40(fpath)
fpath='E:/Data/TRMM/3b40/3b40_FileHead.txt'
fout=open(fpath,'w')
print(trmm_file.header())
ax=trmm_file.header()
headnm=ax.keys()
headva=ax.values()
na=len(headnm)
for a in range(0,na) :
    item="%s: "%headnm[a]
    fout.write(item)
    item="%s "%headva[a]
    fout.write(item)
    fout.write('\n')
precip = trmm_file.precip()
print('Array dimensions:', precip.shape)
fout.write('Array dimensions:')
for a in precip.shape:
    item="%d "%a
    fout.write(item)
fout.close()
del trmm_file
nbin=99
WATERBIN=np.ndarray(shape=(nbin),dtype=float)
WATERBIN[0]=0.01  # MIN
C1=WATERBIN[0]
IC1=0
C3=0.01
IC2=-1
IC3=-1
IC4=-1
for I in range(1, nbin):
    C2=C1+(I-IC1)*C3
    if C2>=0.1 and C2 <1.1 and IC2< 0 :
        C1=0.1 ; IC1=I ; C3=0.1
        IC2=1    #! MAKE SURE THE IF BLOCK JUST CALLED ONCE, SO IC1 IS RIGHT
    elif C2>=1. and C2<10. and IC3< 0 :
        C1=1. ; IC1=I ; C3=0.5
        IC3=1
    elif C2>=10 and IC4< 0 :
        C1=10. ; IC1=I ; C3=1.
        IC4=1
    WATERBIN[I]=C2
print WATERBIN
fig,ax0=plt.subplots(nrows=6,ncols=4,figsize=(12,8))
"""
for i in range(0,12):
    print i
    i+=4
