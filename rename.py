#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 10 11:31:23 2015

@author: jhchen
"""

import os,sys
dirin = "D:/MyPaper/PhD04/Cases/ETP/20100604_0704/Simulated/Dayliy/TCWC/"
filenames=os.listdir(dirin)
print filenames[2]
l=len(filenames[2])
print filenames[2][l-6:l-4]
for a in range(0,len(filenames)):
    l=len(filenames[a])
    os.rename(dirin+filenames[a],dirin+filenames[a][l-6:l-4]+'.png')