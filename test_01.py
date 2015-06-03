# -*- coding: utf-8 -*-
"""
Created on Wed May 27 13:17:30 2015

@author: jhchen
"""
from pytrmm import TRMM3B40RTFile as TRM3B40
import os
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