#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 11:00:56 2015

@author: jhchen
"""
from netCDF4 import Dataset
import os
import datetime
import calendar
os.system("cls")
dirin='X:/Data/TRMM/3b40/'
dirout='D:/MyPaper/PhD04/Data/TRMM/3B40/'
regions=['ETP','WTP','PRD','MLYR','NPC','NEC']
slon =[90.  , 80. , 107.5 , 110. , 110. , 120. ]
elon =[100. , 90. , 117.5 , 120. , 120. , 130. ]
slat =[27.5 , 27.5, 17.5  , 25.  , 32.5 , 40.  ]
elat =[37.5 , 37.5, 27.5  , 35.  , 42.5 , 50.  ]
# the start date of ever region and how many days for every region.
ayear=[2012 ,2010 , 2012  , 2010 , 2010 , 2012 ]
amon =[5    ,  7  ,  4    ,  6   ,  8   , 7    ] 
aday =[20   ,  14 ,  1    ,  2   ,  1   , 6    ]
ndays=[31   ,  31 ,  31   , 31   ,  31  , 31   ]
nag=len(regions)
ima=1
ida=1
for i in range(0,nag):
    iya=ayear[i]
    ystr='%d'%iya
    fpath=dirin+ystr+'rainfall.nc'
    fnc=Dataset(fpath,'a')
    