#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 19 11:52:19 2014

@author: jhchen
"""
import datetime
import matplotlib.pyplot as plt
import calendar
import string
casenm='NPC2D1'
dirobs='D:/MyPaper/PhD04/Cases/NPC/20100802/'
dirsim=dirobs+'Simulated/'
dircn='D:/MyPaper/PhD04/Data/RainCN05/'
cn05nm=casenm[0:3]+'20100601_030.txt'
dirout='D:/MyPaper/PhD04/Pics/'
simnm='rain_'+casenm+'.txt'
obsnm='2010'+casenm[0:3]+'_OBS_Surface_input.txt'
dignos=casenm[0:3]+'_diagnosed_Rain.txt'
nt=121
#
#;;; set for x axis labels
iy=2010
im=6
jd=1
monstr="%02d"%(im)  ### number to string, 1 to 01, 10 to 10
dnm=calendar.monthrange(iy,im)[1]  #calendar.monthrange(1997,7) 
#                              #reture two index, sencond is the day number
datestart=datetime.datetime(iy,im,jd,0,0,0)
det=datetime.timedelta(hours=6)            
dateiso=[]            
for dt in range(0,dnm*4):
    dateiso.append(datestart+dt*det)
    xdate=[]    
    xdat=range(0,dnm*4)            
    for tm in dateiso:
        xdate.append(datetime.datetime.strftime(tm,"%b/%d")) 
###### fro daily
det=datetime.timedelta(hours=24)            
dateiso=[]            
for dt in range(0,dnm):
    dateiso.append(datestart+dt*det)
    xdated=[]    
    xdatd=range(0,dnm)            
    for tm in dateiso:
        xdated.append(datetime.datetime.strftime(tm,"%b/%d"))        
########## reading files
filepath=dirobs+obsnm
onedim=[]
linesplit=[]
f=open(filepath)
ff=f.readlines()[1:]  ## first line in obs file is legend 
for line in ff:
    line=string.lstrip(line)
    linesplit.append(line[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            onedim.append(string.atof(strs))
nl=len(onedim)
timestep=[]
q1=[]
q2=[]
rainobs=[]
for i in range(0,nl-1,4):
    timestep.append(onedim[i])
    q1.append(onedim[i+1])
    q2.append(onedim[i+2])
    rainobs.append(onedim[i+3])
del onedim
del linesplit
###
filepath=dirsim+simnm
rainsim=[]
linesplit=[]
f=open(filepath)
ff=f.readlines()  ## first line in obs file is 
for line in ff:
    line=string.lstrip(line)
    linesplit.append(line[:-1].split(' '))
for lnstrs in linesplit:
        for strs in lnstrs:
            if strs!='':
                rainsim.append(string.atof(strs))    
### rainsim every six hour and unite is mm/hr
lsim=len(rainsim)
dayrainsim=[]
for i in range(0,lsim-1,4):
    temp=0.
    temp=(rainsim[i]*24.+rainsim[i+1]*24.+rainsim[i+2]*24.+rainsim[i+3]*24.)/4.
    dayrainsim.append(temp)                    

### read cn
#filepath=dircn+cn05nm
#onedimcn=[]
#linesplit=[]
#f=open(filepath)
#ffcn=f.readlines()[0:]  
#for line in ffcn:
#    line=string.lstrip(line)
#    linesplit.append(line[:-1].split(' '))
#for lnstrs in linesplit:
#    for strs in lnstrs:
#        if strs!='':
#            onedimcn.append(string.atof(strs))
#nl=len(onedimcn)
#raincn=[]
#dayid=[]
#temp=[]
#for i in range(0,nl-1,3):
#    dayid.append(onedimcn[i])
#    temp.append(onedimcn[i+1])
#    raincn.append(onedimcn[i+2])  
### read diagnosed
filepath=dirobs+dignos
onedim=[]
linesplit=[]
f=open(filepath)
ff=f.readlines()[0:]  
for line in ff:
    line=string.lstrip(line)
    linesplit.append(line[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            onedim.append(string.atof(strs))
nl=len(onedim)
tmid=[]
tmcn=[]
q1x=[]
q2x=[]
q1rain=[]
q2rain=[]
for i in range(0,nl-1,5):
    tmid.append(onedim[i])
    q1x.append(onedim[i+1])
    q1rain.append(onedim[i+2])
    q2x.append(onedim[i+3])
    q2rain.append(onedim[i+4])    
#####################
bias=[]
for i in range(0,nt-1):
    bias.append(rainsim[i]-rainobs[i])
######
#####　　　　plotting
fig,(ax0,ax1) = plt.subplots(nrows=2,ncols=1,figsize=(12,6))
line1,=ax0.plot(xdat[0:nt-1],rainobs[0:nt-1],'g',label='NCEP')
line2,=ax0.plot(xdat[0:nt-1],rainsim[0:nt-1],'b',label='Simulated')
line3,=ax0.plot(xdat[0:nt-1],q1rain[0:nt-1],'k',label='Q1toRain')
line4,=ax0.plot(xdat[0:nt-1],q2rain[0:nt-1],'c',label='Q2roRain')
ax0.set_title('Observation and Simulated Rainfall (mm/hr)',fontsize=16)
ax0.set_xticks(range(0,dnm*4,12))
xticklabels = [xdate[nn] for nn in range(0,dnm*4,12)] 
ax0.set_xticklabels(xticklabels, size=16)
#
line5,=ax1.plot(xdat[0:nt-1],bias[0:nt-1],'r',label='Bias')
ax1.set_title('Simulated-Observation Rainfall (mm/hr)',fontsize=16)
ax1.set_xticks(range(0,nt-1,12))
xticklabels = [xdate[nn] for nn in range(0,nt-1,12)] 
ax1.set_xticklabels(xticklabels, size=16)
plt.legend(ax0,handles=[line1,line2,line3,line4,line5])
plt.show()                     
plt.savefig(dirout+casenm+'_rain.pdf')          
plt.show()
plt.close()
################### cn05
#biasday=[]
#nd=min(len(dayrainsim),len(raincn))
#for i in range(0,nd-1):
#    biasday.append(dayrainsim[i]-raincn[i])
#####　　　　plotting
#fig,(ax0,ax1) = plt.subplots(nrows=2,ncols=1,figsize=(12,6))
#ax0.plot(xdatd[0:nd-1],raincn[0:nd-1],'g')
#ax0.plot(xdatd[0:nd-1],dayrainsim[0:nd-1],'b')
#ax0.set_title('Observation and Simulated Rainfall (mm/d)',fontsize=16)
#ax0.set_xticks(range(0,nd-1,3))
#xticklabels = [xdated[nn] for nn in range(0,nd-1,3)] 
#ax0.set_xticklabels(xticklabels, size=16)
#
#ax1.plot(xdatd[0:nd-1],biasday[0:nd-1],'r')
#ax1.set_title('Simulated-Observation Rainfall (mm/d)',fontsize=16)
#ax1.set_xticks(range(0,nd-1,3))
#xticklabels = [xdated[nn] for nn in range(0,nd-1,3)] 
#ax1.set_xticklabels(xticklabels, size=16)
#plt.show()                     
#plt.savefig(dirout+casenm+'_rain_cn05.pdf')          
#plt.show()
#plt.close()
