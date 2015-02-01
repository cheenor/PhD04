#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 10 11:31:23 2015

@author: jhchen
"""

import os,sys
os.system("cls")
#dirin = "D:/MyPaper/PhD04/Cases/ETP/20100604_0704/Simulated/Dayliy/TCWC/"
#dirin = "D:/MyPaper/PhD04/Cases/NPC/20100802/Simulated/animate/"
dirin = "D:/MyPaper/PhD04/Cases/NPC/20100802/Simulated/npc2d0/"
filenames=os.listdir(dirin)
print filenames[2]
l=len(filenames[2])
print filenames[2][l-8:l-4]
i=1
for a in range(0,len(filenames)):
    l=len(filenames[a])
    timestr="%04d"%(i)
    if i==1682 :    
        print i
        print timestr
        print a
        filenames[a]
    os.rename(dirin+filenames[a],dirin+timestr+'.png')
    i+=1