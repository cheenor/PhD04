#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 12:01:30 2015

@author: jhchen
"""
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
import calendar
import string
import numpy as np
from pylab import *
casenm='WTPCTR_EC'
#casenm='WTP2D0'
#casenm='NPC2D2'
if casenm[0:3]=='ETP':
    dirceres='D:/MyPaper/PhD04/Data/CERES/'
    dirsim='D:/MyPaper/PhD04/Cases/ETP/20100604_0704/Simulated/'
    simnm='rad_3hr_'+casenm
    ceresnm="ETP_2010_6_4__30d3h.txt"
    ncepnm="ETP_2010_6_4__30d6h_NCEP.txt"
    iy,im,jd=2010,6,4
if casenm[0:3]=='WTP':    
    dirceres='D:/MyPaper/PhD04/Data/CERES/'
    dirsim='D:/MyPaper/PhD04/Cases/WTP/20100624_0723/Simulated/'
    simnm='rad_3hr_'+casenm
    ceresnm="WTP_2010_6_24__30d3h.txt"
    ncepnm="WTP_2010_6_24__30d6h_NCEP.txt"
    iy,im,jd=2010,6,24  
if casenm[0:3]=='NPC':
    dirceres='D:/MyPaper/PhD04/Data/CERES/'
    dirsim='D:/MyPaper/PhD04/Cases/NPC/20100725/Simulated/'
    simnm='rad_3hr_'+casenm
    ceresnm="NPC_2010_7_25__30d3h.txt"
    ncepnm="NPC_2010_7_25__30d6h_NCEP.txt"
    iy,im,jd=2010,7,25
dirout='D:/MyPaper/PhD04/Pics/'
nt=241
#
#;;; set for x axis labels
monstr="%02d"%(im)  ### number to string, 1 to 01, 10 to 10
dnm=calendar.monthrange(iy,im)[1]  #calendar.monthrange(1997,7) 
#                              #reture two index, sencond is the day number
mpl.rcParams['ytick.labelsize'] = 24
mpl.rcParams['xtick.labelsize'] = 24
datestart=datetime.datetime(iy,im,jd,0,0,0)
det=datetime.timedelta(hours=3)            
dateiso=[]            
for dt in range(0,nt):
    dateiso.append(datestart+dt*det)
xdate=[]    
xdat=range(0,nt)            
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
#################
det=datetime.timedelta(hours=6)            
dateiso=[]            
for dt in range(0,nt):
    dateiso.append(datestart+dt*det)
xdatencep=[]    
xdatncep=range(0,121)            
for tm in dateiso:
        xdatencep.append(datetime.datetime.strftime(tm,"%b/%d")) 
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
filepath=dirceres+ceresnm
onedim1=[]
linesplit=[]
f=open(filepath)
ff=f.readlines()[1:]  ## first line in obs file is legend 
for line in ff:
    lines=string.lstrip(line)
    linesplit.append(lines[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            onedim1.append(string.atof(strs))
nl=len(onedim1)
timestep=[]
ceresvrnm=["toa_sw_all_3h", "toa_sw_clr_3h", "toa_lw_all_3h",  # 0,1,2
          "toa_lw_clr_3h", "toa_wn_all_3h", "toa_wn_clr_3h",    #3,4,5
          "toa_comp_sw-up_pri_3h", "toa_comp_sw-up_clr_3h", "toa_comp_sw-up_naer_3h",  #6,7,8
          "toa_comp_sw-up_all_3h", "toa_comp_lw-up_pri_3h", "toa_comp_lw-up_clr_3h",   #9,10,11
          "toa_comp_lw-up_naer_3h", "toa_comp_lw-up_all_3h", "toa_comp_wn-up_pri_3h",   #12,13,14
          "toa_comp_wn-up_clr_3h", "toa_comp_wn-up_naer_3h", "toa_comp_wn-up_all_3h",    #15,16,17
          "toa_comp_sw-down_all_3h", "sfc_comp_sw-up_pri_3h", "sfc_comp_sw-up_clr_3h",   #18,19,20
          "sfc_comp_sw-up_naer_3h", "sfc_comp_sw-up_all_3h", "sfc_comp_sw-down_pri_3h",  #21,22,23
          "sfc_comp_sw-down_clr_3h", "sfc_comp_sw-down_naer_3h", "sfc_comp_sw-down_all_3h",  #24,25,26
          "sfc_comp_lw-up_pri_3h", "sfc_comp_lw-up_clr_3h", "sfc_comp_lw-up_naer_3h",   #27,28,29
          "sfc_comp_lw-up_all_3h", "sfc_comp_lw-down_pri_3h", "sfc_comp_lw-down_clr_3h",  #30,31,32
          "sfc_comp_lw-down_naer_3h", "sfc_comp_lw-down_all_3h", "sfc_comp_wn-up_pri_3h",  #33,34,35
          "sfc_comp_wn-up_clr_3h", "sfc_comp_wn-up_naer_3h", "sfc_comp_wn-up_all_3h",     #36,37,38
          "sfc_comp_wn-down_pri_3h", "sfc_comp_wn-down_clr_3h", "sfc_comp_wn-down_naer_3h",  #39,40,41
          "sfc_comp_wn-down_all_3h"]                                                          #42
nvar=len(ceresvrnm)        
ceresvar=np.ndarray(shape=(nvar,nt), dtype=float) 
for i in range(0,nt):
    for j in range(0,nvar):
        k=i*nvar+j
        ceresvar[j,i]=onedim1[k]
        if onedim1[k] ==-999.0:
            ceresvar[j,i]=0.0            
ceresalb=[]
for i in range(0,nt):
    if ceresvar[9,i]!=0.0:
        ceresalb.append(ceresvar[18,i]/ceresvar[9,i])
    if ceresvar[9,i]==0.0:
        ceresalb.append(0.0) 
############################################################################
### ncep
filepath=dirceres+ncepnm
ncepvar=["ulwrf"] #,"adnst","auplt","adnlt",
#         "aupss","adbss","aupls","adnls"]
nvar=len(ncepvar)
nsncep=121
ulwrfncep=np.ndarray(shape=(nvar,nsncep), dtype=float)          
linesplit=[]
onedim=[]
f=open(filepath)
ff=f.readlines()[1:] ## first line in obs file is 
for line in ff:
    lines=string.lstrip(line)
    linesplit.append(lines[:-1].split(' '))
for lnstrs in linesplit:
        for strs in lnstrs:
            if strs!='':
                onedim.append(string.atof(strs))
### radiaiton every six hours and unite is W/m**2
for i in range(0,nsncep,1):
    for j in range(0,nvar):
        k=i*nvar+j
        ulwrfncep[j,i]= onedim[k]                                
####################################      
filepath=dirsim+simnm
simvars=["aupst","adnst","auplt","adnlt",
         "aupss","adnss","aupls","adnls"]
nvar=len(simvars)
nst=241
simvar=np.ndarray(shape=(nvar,nst), dtype=float)          
linesplit=[]
onedim=[]
f=open(filepath)
ff=f.readlines()  ## first line in obs file is 
for line in ff:
    lines=string.lstrip(line)
    linesplit.append(lines[:-1].split(' '))
for lnstrs in linesplit:
        for strs in lnstrs:
            if strs!='':
                onedim.append(string.atof(strs))    
### radiaiton every three hours and unite is W/m**2
for i in range(0,nst,1):
    for j in range(0,nvar):
        k=i*nvar+j
        simvar[j,i]= onedim[k]                     
alb_sim=[]
for i in range(0,nst,1):
    if simvar[1,i]!=0.0:
        alb_sim.append(simvar[0,i]/simvar[1,i])
    if simvar[1,i]==0:
        alb_sim.append(0.0)
############### for compare with ncep data
ulwrfsim=simvar[6,:]
ulwrfsim6h=np.ndarray(shape=(nsncep), dtype=float) 
for i in range(0,nsncep):
    if i==0:
        ulwrfsim6h[i]=(ulwrfsim[i]+ulwrfsim[i+1])/2.
    elif i==nsncep-1:
        ulwrfsim6h[i]=(ulwrfsim[i*2]+ulwrfsim[i*2-1])/2.
    else:    
        ulwrfsim6h[i]=(ulwrfsim[i*2]+ulwrfsim[i*2-1]+ulwrfsim[i*2+1])/3.        
####################
nd=30
meansimvar=np.ndarray(shape=(nvar,nd), dtype=float)  
meansimalb=np.ndarray(shape=(nd), dtype=float)        
idd=0
for i in range(0,nst-8,8): 
    for k in range(0,nvar): 
        temp=0.0
        for j in range(0,8):
            kk=i+j
            temp=temp+simvar[k,kk]            
        meansimvar[k,idd]=temp/8.0
    idd=idd+1
idd=0
for i in range(0,nst-8,8):
    ik=0
    temp=0.0
    for j in range(0,8):
        kk=i+j
        if alb_sim[kk] !=0.0:
            temp=temp+alb_sim[kk]
            ik=ik+1
    meansimalb[idd]=temp/ik
    idd=idd+1
##########################################################3        
################### toa sw
nd=30
avg1=0.
avg2=0.
for i in range(0,nt):
    avg1=avg1+ceresvar[9,i]
    avg2=avg2+simvar[0,i]
avg1=avg1/nt
avg2=avg2/nt 
avgstr1="%.2f"%avg1
avgstr2="%.2f"%avg2
ndstr="%d"%nd  
ndstr=ndstr+"days"
#####　　　　plotting
fig,(ax0,ax1) = plt.subplots(nrows=2,ncols=1,figsize=(12,6))
lcc=len(ceresvrnm[9])
ax0.set_ylim(0,700)
ax0.plot(xdat[0:nt-1],ceresvar[9,0:nt-1],'g',label=ceresvrnm[9])#[0:lcc-3])
ax0.plot(xdat[0:nt-1],simvar[0,0:nt-1],'b',label=simvars[0])
#ax0.set_title('Upward shorwave radiaiton at TOA '+ r' $W m^{-2}$',fontsize=16)
text1=r"($a$) Upward shorwave radiaiton at TOA"
ax0.text(2,610,text1,fontsize=16)
text1=ndstr+r" averaged "+"CERES: "+avgstr1+ r' $W$ $m^{-2}$'
ax0.text(150,610,text1)
text1=ndstr+r" averaged "+"Simulated: "+avgstr2+ r' $W$ $m^{-2}$'
ax0.text(150,550,text1)
ax0.set_yticks(range(0,700,150))
ax0.set_xticks(range(0,nt,24))
xticklabels = [xdate[nn] for nn in range(0,nt,24)] 
ax0.set_xticklabels(xticklabels, size=16)
#
avg1=0.
avg2=0.
for i in range(0,nt):
    avg1=avg1+ceresvar[18,i]
    avg2=avg2+simvar[1,i]
avg1=avg1/nt
avg2=avg2/nt 
avgstr1="%.2f"%avg1
avgstr2="%.2f"%avg2
lcc=len(ceresvrnm[18])
ax1.plot(xdat[0:nt-1],ceresvar[18,0:nt-1],'g',label=ceresvrnm[18])#[0:lcc-3])
ax1.plot(xdat[0:nt-1],simvar[1,0:nt-1],'b',label=simvars[1])
#ax1.set_title('Downward shorwave radiaiton at TOA '+ r' $W m^{-2}$',fontsize=16)
#text1="CERES "+ndstr+" averaged is "+avgstr1+ r' $W m^{-2}$'
#ax1.text(1,1300,text1)
#text1="Simulated "+ndstr+" averaged is "+avgstr2+ r' $W m^{-2}$'
#ax1.text(1,1100,text1)
ax1.set_ylim(0,1500,300)
text1=r"($b$) Downward shorwave radiaiton at TOA"
ax1.text(2,1450,text1,fontsize=16)
text1=ndstr+r" averaged "+"CERES: "+avgstr1+ r' $W$ $m^{-2}$'
ax1.text(150,1450,text1)
text1=ndstr+r" averaged "+"Simulated: "+avgstr2+ r' $W$ $m^{-2}$'
ax1.text(150,1300,text1)
ax1.set_yticks(range(0,1600,300))
ax1.set_xticks(range(0,nt,24))
xticklabels = [xdate[nn] for nn in range(0,nt,24)] 
ax1.set_xticklabels(xticklabels, size=16)
plt.show()                     
plt.savefig(dirout+casenm+'_TOA_SW-3h.pdf')          
plt.show()
plt.close()
#  TOA LW
avg1=0.
avg2=0.
for i in range(0,nt):
    avg1=avg1+ceresvar[13,i]
    avg2=avg2+simvar[2,i]
avg1=avg1/nt
avg2=avg2/nt 
avgstr1="%.2f"%avg1
avgstr2="%.2f"%avg2
#fig,(ax0,ax1) = plt.subplots(nrows=2,ncols=1,figsize=(12,6))
fig,(ax0) = plt.subplots(nrows=1,ncols=1,figsize=(12,6))
lcc=len(ceresvrnm[13])
ax0.plot(xdat[0:nt-1],ceresvar[13,0:nt-1],'g',label=ceresvrnm[13])#[0:lcc-3])
ax0.plot(xdat[0:nt-1],simvar[2,0:nt-1],'b',label=simvars[0])
#ax0.set_title('Upward longwave radiaiton at TOA '+ r' $W m^{-2}$',fontsize=16)
#text1="CERES "+ndstr+" averaged is "+avgstr1+ r' $W m^{-2}$'
#ax0.text(1,248,text1)
#text1="Simulated "+ndstr+" averaged is "+avgstr2+ r' $W m^{-2}$'
#ax0.text(1,238,text1)
ax0.set_ylim(150,280)
text1=r"($a$) Upward longwave radiaiton at TOA"
ax0.text(2,265,text1,fontsize=16)
text1=ndstr+r" averaged "+"CERES: "+avgstr1+ r' $W$ $m^{-2}$'
ax0.text(150,265,text1)
text1=ndstr+r" averaged "+"Simulated: "+avgstr2+ r' $W$ $m^{-2}$'
ax0.text(150,255,text1)
ax0.set_yticks(range(150,280,30))
ax0.set_xticks(range(0,nt,24))
xticklabels = [xdate[nn] for nn in range(0,nt,24)] 
ax0.set_xticklabels(xticklabels, size=16)
#
fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(10,8))
ax.set_ylim(170,250)
ax.set_xlim(170,250)
plot(ceresvar[13,0:nt-1],simvar[2,0:nt-1],'bo')
plot([170,250],[170,250],'b-',lw=2.6)
text1=r"($a$) Upward Longwave Radiaiton at TOA"
ax.text(170,251,text1,fontsize=16)
text1=ndstr+r" averaged "+"CERES: "+avgstr1+ r' $W$ $m^{-2}$'
ax.text(172,246,text1)
text1=ndstr+r" averaged "+"Simulated: "+avgstr2+ r' $W$ $m^{-2}$'
ax.text(172,243,text1)
ax.set_yticks(range(170,250,20))
ax.set_xticks(range(170,250,20))
ax.set_xlabel("CERES",fontsize=24)
ax.set_ylabel("Simulated",fontsize=24)
plt.show()                     
plt.savefig(dirout+casenm+'_TOA_UPLW-3h_scatter.pdf')          
plt.show()
plt.close()
"""
avg1=0.
avg2=0.
for i in range(0,nt):
    avg1=avg1+ceresvar[7,i]
    avg2=avg2+simvar[3,i]
avg1=avg1/nt
avg2=avg2/nt 
avgstr1="%.2f"%avg1
avgstr2="%.2f"%avg2
#ax1.plot(xdat[0:nt-1],toavar[12,nt-1],'g',label=toanms[12])
ax1.plot(xdat[0:nt-1],simvar[3,0:nt-1],'b',label=simvars[3])
#ax1.set_title('Downward longwave radiaiton at TOA '+ r' $W m^{-2}$',fontsize=16)
#text1="Simulated "+ndstr+" averaged is "+avgstr2+ r' $W m^{-2}$'
text1=r"($b$) Downward longwave radiaiton at TOA"
ax1.text(1,0.40,text1)
text1=ndstr+r" averaged "+"Simulated: "+avgstr2+ r' $W$ $m^{-2}$'
ax1.text(150,0.40,text1)
#ax1.set_yticks(range(0.3,0.5))
ax1.set_xticks(range(0,nt,24))
xticklabels = [xdate[nn] for nn in range(0,nt,24)] 
ax1.set_xticklabels(xticklabels, size=16)
"""
plt.show()                     
plt.savefig(dirout+casenm+'_TOA_LW-3h.pdf')          
plt.show()
plt.close()
###############################################################################
#####　　　　plotting
avg1=0.
avg2=0.
for i in range(0,nt):
    avg1=avg1+ceresvar[22,i]
    avg2=avg2+simvar[4,i]
avg1=avg1/nt
avg2=avg2/nt 
avgstr1="%.2f"%avg1
avgstr2="%.2f"%avg2
###
fig,(ax0,ax1) = plt.subplots(nrows=2,ncols=1,figsize=(12,6))
lcc=len(ceresvrnm[22])
ax0.plot(xdat[0:nt-1],ceresvar[22,0:nt-1],'g',label=ceresvrnm[22])#[0:lcc-3])
ax0.plot(xdat[0:nt-1],simvar[4,0:nt-1],'b',label=simvars[4])
#ax0.set_title('Upward shorwave radiaiton at SRF '+ r' $W m^{-2}$',fontsize=16)
#text1="CERES "+ndstr+" averaged is "+avgstr1+ r' $W m^{-2}$'
#ax0.text(1,175,text1)
#text1="Simulated "+ndstr+" averaged is "+avgstr2+ r' $W m^{-2}$'
#ax0.text(1,150,text1)
ax0.set_ylim(0,210)
text1=r"($a$) Upward shorwave radiaiton at SRF"
ax0.text(2,185,text1,fontsize=16)
text1=ndstr+r" averaged "+"CERES: "+avgstr1+ r' $W$ $m^{-2}$'
ax0.text(150,185,text1)
text1=ndstr+r" averaged "+"Simulated: "+avgstr2+ r' $W$ $m^{-2}$'
ax0.text(150,165,text1)
ax0.set_yticks(range(0,210,40))
ax0.set_xticks(range(0,nt,24))
xticklabels = [xdate[nn] for nn in range(0,nt,24)] 
ax0.set_xticklabels(xticklabels, size=16)
#
avg1=0.
avg2=0.
for i in range(0,nt):
    avg1=avg1+ceresvar[26,i]
    avg2=avg2+simvar[5,i]
avg1=avg1/nt
avg2=avg2/nt 
avgstr1="%.2f"%avg1
avgstr2="%.2f"%avg2
lcc=len(ceresvrnm[26])
ax1.plot(xdat[0:nt-1],ceresvar[26,0:nt-1],'g',label=ceresvrnm[26])#[0:lcc-3])
ax1.plot(xdat[0:nt-1],simvar[5,0:nt-1],'b',label=simvars[5])
#ax1.set_title('Downward shorwave radiaiton at SRF '+ r' $W m^{-2}$',fontsize=16)
#text1="CERES "+ndstr+" averaged is "+avgstr1+ r' $W m^{-2}$'
#ax1.text(1,900,text1)
#text1="Simulated "+ndstr+" averaged is "+avgstr2+ r' $W m^{-2}$'
#ax1.text(1,750,text1)
ax1.set_ylim(0,1000)
text1=r"($b$) Downward shorwave radiaiton at SRF"
ax1.text(2,900,text1,fontsize=16)
text1=ndstr+r" averaged "+"CERES: "+avgstr1+ r' $W$ $m^{-2}$'
ax1.text(150,900,text1)
text1=ndstr+r" averaged "+"Simulated: "+avgstr2+ r' $W$ $m^{-2}$'
ax1.text(150,800,text1)
ax1.set_yticks(range(0,1000,200))
ax1.set_xticks(range(0,nt,24))
xticklabels = [xdate[nn] for nn in range(0,nt,24)] 
ax1.set_xticklabels(xticklabels, size=16)
plt.show()                     
plt.savefig(dirout+casenm+'_SRF_SW-3h.pdf')          
plt.show()
plt.close()
#  TOA LW
avg1=0.
avg2=0.
for i in range(0,nt):
    avg1=avg1+ceresvar[30,i]
    avg2=avg2+simvar[6,i]
avg1=avg1/nt
avg2=avg2/nt 
avgstr1="%.2f"%avg1
avgstr2="%.2f"%avg2
fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(10,8))
ax.set_ylim(260,430)
ax.set_xlim(260,430)
plot(ceresvar[30,0:nt-1],simvar[6,0:nt-1],'bo')
plot([250,450],[250,450],'b-',lw=2.6)
text1=r"($a$) Upward Longwave Radiaiton at SRF"
ax.text(260,432,text1,fontsize=16)
text1=ndstr+r" averaged "+"CERES: "+avgstr1+ r' $W$ $m^{-2}$'
ax.text(262,420,text1)
text1=ndstr+r" averaged "+"Simulated: "+avgstr2+ r' $W$ $m^{-2}$'
ax.text(262,410,text1)
ax.set_yticks(range(260,430,30))
ax.set_xticks(range(260,430,30))
ax.set_xlabel("CERES",fontsize=24)
ax.set_ylabel("Simulated",fontsize=24)
plt.show()                     
plt.savefig(dirout+casenm+'_SRF_UPLW-3h_scatter.pdf')          
plt.show()
plt.close()
lcc=len(ceresvrnm[30])
fig,(ax0,ax1) = plt.subplots(nrows=2,ncols=1,figsize=(12,6))
ax0.plot(xdat[0:nt-1],ceresvar[30,0:nt-1],'g',label=ceresvrnm[30])#[0:lcc-3])
ax0.plot(xdat[0:nt-1],simvar[6,0:nt-1],'b',label=simvars[6])
#ax0.set_title('Upward longwave radiaiton at SRF '+ r' $W m^{-2}$',fontsize=16)
#text1="CERES "+ndstr+" averaged is "+avgstr1+ r' $W m^{-2}$'
#ax0.text(1,400,text1)
#text1="Simulated "+ndstr+" averaged is "+avgstr2+ r' $W m^{-2}$'
#ax0.text(1,380,text1)
ax0.set_ylim(250,450)
text1=r"($a$) Upward longwave radiaiton at SRF"
ax0.text(2,420,text1,fontsize=16)
text1=ndstr+r" averaged "+"CERES: "+avgstr1+ r' $W$ $m^{-2}$'
ax0.text(150,420,text1)
text1=ndstr+r" averaged "+"Simulated: "+avgstr2+ r' $W$ $m^{-2}$'
ax0.text(150,390,text1)
ax0.set_yticks(range(250,450,50))
ax0.set_xticks(range(0,nt,24))
xticklabels = [xdate[nn] for nn in range(0,nt,24)] 
ax0.set_xticklabels(xticklabels, size=16)
#
avg1=0.
avg2=0.
for i in range(0,nt):
    avg1=avg1+ceresvar[34,i]
    avg2=avg2+simvar[7,i]
avg1=avg1/nt
avg2=avg2/nt 
avgstr1="%.2f"%avg1
avgstr2="%.2f"%avg2
lcc=len(ceresvrnm[34])
ax1.plot(xdat[0:nt-1],ceresvar[34,0:nt-1],'g',label=ceresvrnm[34])#[0:lcc-3])
ax1.plot(xdat[0:nt-1],simvar[7,0:nt-1],'b',label=simvars[7])
#ax1.set_title('Downward longwave radiaiton at SRF '+ r' $W m^{-2}$',fontsize=16)
#text1="CERES "+ndstr+" averaged is "+avgstr1+ r' $W m^{-2}$'
#ax1.text(1,340,text1)
#text1="Simulated "+ndstr+" averaged is "+avgstr2+ r' $W m^{-2}$'
#ax1.text(1,320,text1)
ax1.set_ylim(180,380)
text1=r"($b$) Downward longwave radiaiton at SRF"
ax1.text(2,360,text1,fontsize=16)
text1=ndstr+r" averaged "+"CERES: "+avgstr1+ r' $W$ $m^{-2}$'
ax1.text(150,360,text1)
text1=ndstr+r" averaged "+"Simulated: "+avgstr2+ r' $W$ $m^{-2}$'
ax1.text(150,340,text1)
ax1.set_yticks(range(180,380,40))
ax1.set_xticks(range(0,nt,24))
xticklabels = [xdate[nn] for nn in range(0,nt,24)] 
ax1.set_xticklabels(xticklabels, size=16)
plt.show()                     
plt.savefig(dirout+casenm+'_SRF_LW-3h.pdf')          
plt.show()
plt.close()
##########################
fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(10,8))
ax.set_ylim(180,380)
ax.set_xlim(180,380)
plot(ceresvar[34,0:nt-1],simvar[7,0:nt-1],'bo')
plot([180,380],[180,380],'b-',lw=2.6)
text1=r"($b$)  Downward Longwave Radiaiton at SRF"
ax.text(180,382,text1,fontsize=16)
text1=ndstr+r" averaged "+"CERES: "+avgstr1+ r' $W$ $m^{-2}$'
ax.text(182,370,text1)
text1=ndstr+r" averaged "+"Simulated: "+avgstr2+ r' $W$ $m^{-2}$'
ax.text(182,360,text1)
ax.set_yticks(range(180,380,30))
ax.set_xticks(range(180,380,30))
ax.set_xlabel("CERES",fontsize=24)
ax.set_ylabel("Simulated",fontsize=24)
plt.show()                     
plt.savefig(dirout+casenm+'_SRF_DNLW-3h_scatter.pdf')          
plt.show()
plt.close()
############
#  TOA LW
avg1=0.
avg2=0.
for i in range(0,nt):
    avg1=avg1+ceresalb[i]
    avg2=avg2+alb_sim[i]
avg1=avg1/nt
avg2=avg2/nt 
avgstr1="%.3f"%avg1
avgstr2="%.3f"%avg2
fig,ax0 = plt.subplots(nrows=1,ncols=1,figsize=(12,8))
ax0.plot(xdat[0:nt-1],ceresalb[0:nt-1],'g',label='Observation')
ax0.plot(xdat[0:nt-1],alb_sim[0:nt-1],'b',label='Simulated')
ax0.set_title('All Sky Albedo',fontsize=16)
text1="CERES "+ndstr+" averaged is "+avgstr1
ax0.text(1,0.54,text1)
text1="Simulated "+ndstr+" averaged is "+avgstr2
ax0.text(1,0.52,text1)
ax0.set_xticks(range(0,nt,24))
xticklabels = [xdate[nn] for nn in range(0,nt,24)] 
ax0.set_xticklabels(xticklabels, size=16)
#
plt.show()                     
plt.savefig(dirout+casenm+'_all_sky_albedo-3h.pdf')          
plt.show()
plt.close()
###############################################################################
for i in range(0,nsncep):
    avg1=avg1+ulwrfsim6h[i]
    avg2=avg2+ulwrfncep[0,i]
avg1=avg1/nsncep
avg2=avg2/nsncep 
avgstr1="%.3f"%avg1
avgstr2="%.3f"%avg2
fig,ax0 = plt.subplots(nrows=1,ncols=1,figsize=(12,8))
ax0.plot(xdatncep[0:nsncep-1],ulwrfncep[0,0:nsncep-1],'g',label='Observation')
ax0.plot(xdatncep[0:nsncep-1],ulwrfsim6h[0:nsncep-1],'b',label='Simulated')
ax0.set_title('Upward Longwave at Surface',fontsize=16)
text1="NCEP "+ndstr+" averaged is "+avgstr2
ax0.text(1,390,text1)
text1="Simulated "+ndstr+" averaged is "+avgstr1
ax0.text(1,370,text1)
ax0.set_xticks(range(0,nsncep,12))
xticklabels = [xdatencep[nn] for nn in range(0,nsncep,12)] 
ax0.set_xticklabels(xticklabels, size=16)
#
plt.show()                     
plt.savefig(dirout+casenm+'_Upward_Longwave_Surface_NCEP.pdf')          
plt.show()
plt.close()
