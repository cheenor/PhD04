#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 22 11:30:22 2014

@author: jhchen
"""
import datetime
import matplotlib.pyplot as plt
import calendar
import string
import numpy as np
casenm='ETP'
dirceres='D:/MyPaper/PhD04/Data/CERES/'
dirsim='D:/MyPaper/PhD04/Cases/ETP/201006/Simulated/'
dirout='D:/MyPaper/PhD04/Pics/'
simnm='rad_3hr_ETP2D0'
ceres_toa=casenm+"_TOA_CMP_2010_6_4__30d.txt"
ceres_srf=casenm+"_SRF_CMP_2010_6_4__30d.txt"
ceres_toa_obs=casenm+"_TOA_OBS_2010_6_4__30d.txt"
nt=121
#
#;;; set for x axis labels
iy=2010
im=6
jd=4
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
filepath=dirceres+ceres_toa
onedim1=[]
linesplit=[]
f=open(filepath)
ff=f.readlines()[1:]  ## first line in obs file is legend 
for line in ff:
    line=string.lstrip(line)
    linesplit.append(line[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            onedim1.append(string.atof(strs))
nl=len(onedim1)
timestep=[]
toanms=["toa_sw-up_pri", "toa_sw-up_clr","toa_sw-up_naer","toa_sw-up_all",
        "toa_lw-up_pri","toa_lw-up_clr","toa_lw-up_naer","toa_lw-up_all",
        "toa_wn-up_pri","toa_wn-up_clr","toa_wn-up_naer","toa_wn-up_all",
        "toa_sw-down_all"]
nvar=len(toanms)        
toavar=np.ndarray(shape=(nvar,30), dtype=float) 
nd=30
for i in range(0,nd):
    for j in range(0,nvar):
        k=i*nvar+j        
        toavar[j,i]=onedim1[k]
############################################################################
filepath=dirceres+ceres_toa_obs
onedim1=[]
linesplit=[]
f=open(filepath)
ff=f.readlines()[1:]  ## first line in obs file is legend 
for line in ff:
    line=string.lstrip(line)
    linesplit.append(line[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            onedim1.append(string.atof(strs))
nl=len(onedim1)
timestep=[]
toanms_obs=["toa_sw_all_daily", "toa_sw_clr_daily", "toa_lw_all_daily",
            "toa_lw_clr_daily", "toa_wn_all_daily", "toa_wn_clr_daily",
            "toa_net_all_daily", "toa_net_clr_daily","toa_alb_all_daily",
            "toa_alb_clr_daily", "toa_solar_all_daily"]
nvar=len(toanms_obs)        
toavarobs=np.ndarray(shape=(nvar,30), dtype=float) 
nd=30
for i in range(0,nd):
    for j in range(0,nvar):
        k=i*nvar+j        
        toavarobs[j,i]=onedim1[k]        
############################################################################        
filepath=dirceres+ceres_srf
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
srfnms=["sfc_sw-up_pri", "sfc_sw-up_clr", "sfc_sw-up_naer", "sfc_sw-up_all",
        "sfc_sw-down_pri","sfc_sw-down_clr","sfc_sw-down_naer","sfc_sw-down_all",
        "sfc_lw-up_pri","sfc_lw-up_clr","sfc_lw-up_naer", "sfc_lw-up_all",
        "sfc_lw-down_pri","sfc_lw-down_clr","sfc_lw-down_naer","sfc_lw-down_all",
        "sfc_wn-up_pri","sfc_wn-up_clr","sfc_wn-up_naer","sfc_wn-up_all", 
        "sfc_wn-down_pri","sfc_wn-down_clr","sfc_wn-down_naer", "sfc_wn-down_all"]
nvar=len(srfnms)        
srfvar=np.ndarray(shape=(nvar,30), dtype=float) 
nd=30
for i in range(0,nd):
    for j in range(0,nvar):
        k=i*nvar+j        
        srfvar[j,i]=onedim[k]
####################################      
filepath=dirsim+simnm
simvars=["aupst","adnst","auplt","adnlt",
         "aupss","adbss","aupls","adnls"]
nvar=len(simvars)
nst=241
simvar=np.ndarray(shape=(nvar,nst), dtype=float)          
linesplit=[]
onedim=[]
f=open(filepath)
ff=f.readlines()  ## first line in obs file is 
for line in ff:
    line=string.lstrip(line)
    linesplit.append(line[:-1].split(' '))
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
for i in range(0,nd):
    avg1=avg1+toavar[3,i]
    avg2=avg2+meansimvar[0,i]
avg1=avg1/nd
avg2=avg2/nd 
avgstr1="%.2f"%avg1
avgstr2="%.2f"%avg2
ndstr="%d"%nd  
ndstr=ndstr+"days"
#####　　　　plotting
fig,(ax0,ax1) = plt.subplots(nrows=2,ncols=1,figsize=(12,6))
ax0.plot(xdatd[0:nd-1],toavar[3,0:nd-1],'g',label=toanms[3])
ax0.plot(xdatd[0:nd-1],meansimvar[0,0:nd-1],'b',label=simvars[0])
ax0.set_title('Upward shorwave radiaiton at TOA '+ r' $W m^{-2}$',fontsize=16)
text1="CERES "+ndstr+" averaged is "+avgstr1+ r' $W m^{-2}$'
ax0.text(1,225,text1)
text1="Simulated "+ndstr+" averaged is "+avgstr2+ r' $W m^{-2}$'
ax0.text(1,210,text1)
ax0.set_xticks(range(0,nd,3))
xticklabels = [xdated[nn] for nn in range(0,nd,3)] 
ax0.set_xticklabels(xticklabels, size=16)
#
avg1=0.
avg2=0.
for i in range(0,nd):
    avg1=avg1+toavar[12,i]
    avg2=avg2+meansimvar[1,i]
avg1=avg1/nd
avg2=avg2/nd 
avgstr1="%.2f"%avg1
avgstr2="%.2f"%avg2
ax1.plot(xdatd[0:nd-1],toavar[12,0:nd-1],'g',label=toanms[12])
ax1.plot(xdatd[0:nd-1],meansimvar[1,0:nd-1],'b',label=simvars[1])
ax1.set_title('Downward shorwave radiaiton at TOA '+ r' $W m^{-2}$',fontsize=16)
text1="CERES "+ndstr+" averaged is "+avgstr1+ r' $W m^{-2}$'
ax1.text(1,487.5,text1)
text1="Simulated "+ndstr+" averaged is "+avgstr2+ r' $W m^{-2}$'
ax1.text(1,485,text1)
ax1.set_xticks(range(0,nd,3))
xticklabels = [xdated[nn] for nn in range(0,nd,3)] 
ax1.set_xticklabels(xticklabels, size=16)
plt.show()                     
plt.savefig(dirout+casenm+'_TOA_SW.pdf')          
plt.show()
plt.close()
#  TOA LW
avg1=0.
avg2=0.
for i in range(0,nd):
    avg1=avg1+toavar[7,i]
    avg2=avg2+meansimvar[2,i]
avg1=avg1/nd
avg2=avg2/nd 
avgstr1="%.2f"%avg1
avgstr2="%.2f"%avg2
fig,(ax0,ax1) = plt.subplots(nrows=2,ncols=1,figsize=(12,6))
ax0.plot(xdatd[0:nd-1],toavar[7,0:nd-1],'g',label=toanms[7])
ax0.plot(xdatd[0:nd-1],meansimvar[2,0:nd-1],'b',label=simvars[0])
ax0.set_title('Upward longwave radiaiton at TOA '+ r' $W m^{-2}$',fontsize=16)
text1="CERES "+ndstr+" averaged is "+avgstr1+ r' $W m^{-2}$'
ax0.text(1,240,text1)
text1="Simulated "+ndstr+" averaged is "+avgstr2+ r' $W m^{-2}$'
ax0.text(1,228,text1)
ax0.set_xticks(range(0,nd,3))
xticklabels = [xdated[nn] for nn in range(0,nd,3)] 
ax0.set_xticklabels(xticklabels, size=16)
#
avg1=0.
avg2=0.
for i in range(0,nd):
    avg1=avg1+toavar[7,i]
    avg2=avg2+meansimvar[3,i]
avg1=avg1/nd
avg2=avg2/nd 
avgstr1="%.2f"%avg1
avgstr2="%.2f"%avg2
#ax1.plot(xdatd[0:nd-1],toavar[12,nd-1],'g',label=toanms[12])
ax1.plot(xdatd[0:nd-1],meansimvar[3,0:nd-1],'b',label=simvars[3])
ax1.set_title('Downward longwave radiaiton at TOA '+ r' $W m^{-2}$',fontsize=16)
text1="Simulated "+ndstr+" averaged is "+avgstr2+ r' $W m^{-2}$'
ax1.text(1,0.36,text1)
ax1.set_xticks(range(0,nd,3))
xticklabels = [xdated[nn] for nn in range(0,nd,3)] 
ax1.set_xticklabels(xticklabels, size=16)
plt.show()                     
plt.savefig(dirout+casenm+'_TOA_LW.pdf')          
plt.show()
plt.close()
###############################################################################
#####　　　　plotting
avg1=0.
avg2=0.
for i in range(0,nd):
    avg1=avg1+srfvar[3,i]
    avg2=avg2+meansimvar[4,i]
avg1=avg1/nd
avg2=avg2/nd 
avgstr1="%.2f"%avg1
avgstr2="%.2f"%avg2
fig,(ax0,ax1) = plt.subplots(nrows=2,ncols=1,figsize=(12,6))
ax0.plot(xdatd[0:nd-1],srfvar[3,0:nd-1],'g',label=srfnms[3])
ax0.plot(xdatd[0:nd-1],meansimvar[4,0:nd-1],'b',label=simvars[4])
ax0.set_title('Upward shorwave radiaiton at SRF '+ r' $W m^{-2}$',fontsize=16)
text1="CERES "+ndstr+" averaged is "+avgstr1+ r' $W m^{-2}$'
ax0.text(1,65,text1)
text1="Simulated "+ndstr+" averaged is "+avgstr2+ r' $W m^{-2}$'
ax0.text(1,40,text1)
ax0.set_xticks(range(0,nd,3))
xticklabels = [xdated[nn] for nn in range(0,nd,3)] 
ax0.set_xticklabels(xticklabels, size=16)
#
avg1=0.
avg2=0.
for i in range(0,nd):
    avg1=avg1+srfvar[7,i]
    avg2=avg2+meansimvar[5,i]
avg1=avg1/nd
avg2=avg2/nd 
avgstr1="%.2f"%avg1
avgstr2="%.2f"%avg2
ax1.plot(xdatd[0:nd-1],srfvar[7,0:nd-1],'g',label=srfnms[7])
ax1.plot(xdatd[0:nd-1],meansimvar[5,0:nd-1],'b',label=simvars[5])
ax1.set_title('Downward shorwave radiaiton at SRF '+ r' $W m^{-2}$',fontsize=16)
text1="CERES "+ndstr+" averaged is "+avgstr1+ r' $W m^{-2}$'
ax1.text(1,320,text1)
text1="Simulated "+ndstr+" averaged is "+avgstr2+ r' $W m^{-2}$'
ax1.text(1,200,text1)
ax1.set_xticks(range(0,nd,3))
xticklabels = [xdated[nn] for nn in range(0,nd,3)] 
ax1.set_xticklabels(xticklabels, size=16)
plt.show()                     
plt.savefig(dirout+casenm+'_SRF_SW.pdf')          
plt.show()
plt.close()
#  TOA LW
avg1=0.
avg2=0.
for i in range(0,nd):
    avg1=avg1+srfvar[11,i]
    avg2=avg2+meansimvar[6,i]
avg1=avg1/nd
avg2=avg2/nd 
avgstr1="%.2f"%avg1
avgstr2="%.2f"%avg2
fig,(ax0,ax1) = plt.subplots(nrows=2,ncols=1,figsize=(12,6))
ax0.plot(xdatd[0:nd-1],srfvar[11,0:nd-1],'g',label=srfnms[11])
ax0.plot(xdatd[0:nd-1],meansimvar[6,0:nd-1],'b',label=simvars[6])
ax0.set_title('Upward longwave radiaiton at SRF '+ r' $W m^{-2}$',fontsize=16)
text1="CERES "+ndstr+" averaged is "+avgstr1+ r' $W m^{-2}$'
ax0.text(1,360,text1)
text1="Simulated "+ndstr+" averaged is "+avgstr2+ r' $W m^{-2}$'
ax0.text(1,350,text1)
ax0.set_xticks(range(0,nd,3))
xticklabels = [xdated[nn] for nn in range(0,nd,3)] 
ax0.set_xticklabels(xticklabels, size=16)
#
avg1=0.
avg2=0.
for i in range(0,nd):
    avg1=avg1+srfvar[15,i]
    avg2=avg2+meansimvar[7,i]
avg1=avg1/nd
avg2=avg2/nd 
avgstr1="%.2f"%avg1
avgstr2="%.2f"%avg2
ax1.plot(xdatd[0:nd-1],srfvar[15,0:nd-1],'g',label=srfnms[15])
ax1.plot(xdatd[0:nd-1],meansimvar[7,0:nd-1],'b',label=simvars[7])
ax1.set_title('Downward longwave radiaiton at SRF '+ r' $W m^{-2}$',fontsize=16)
text1="CERES "+ndstr+" averaged is "+avgstr1+ r' $W m^{-2}$'
ax1.text(1,320,text1)
text1="Simulated "+ndstr+" averaged is "+avgstr2+ r' $W m^{-2}$'
ax1.text(1,300,text1)
ax1.set_xticks(range(0,nd,3))
xticklabels = [xdated[nn] for nn in range(0,nd,3)] 
ax1.set_xticklabels(xticklabels, size=16)
plt.show()                     
plt.savefig(dirout+casenm+'_SRF_LW.pdf')          
plt.show()
plt.close()
############
#  TOA LW
avg1=0.
avg2=0.
for i in range(0,nd):
    avg1=avg1+toavarobs[8,i]
    avg2=avg2+meansimalb[i]
avg1=avg1/nd
avg2=avg2/nd 
avgstr1="%.3f"%avg1
avgstr2="%.3f"%avg2
fig,ax0 = plt.subplots(nrows=1,ncols=1,figsize=(12,8))
ax0.plot(xdatd[0:nd-1],toavarobs[8,0:nd-1],'g',label='Observation')
ax0.plot(xdatd[0:nd-1],meansimalb[0:nd-1],'b',label='Simulated')
ax0.set_title('All Sky Albedo',fontsize=16)
text1="CERES "+ndstr+" averaged is "+avgstr1
ax0.text(1,0.54,text1)
text1="Simulated "+ndstr+" averaged is "+avgstr2
ax0.text(1,0.52,text1)
ax0.set_xticks(range(0,nd,3))
xticklabels = [xdated[nn] for nn in range(0,nd,3)] 
ax0.set_xticklabels(xticklabels, size=16)
#
plt.show()                     
plt.savefig(dirout+casenm+'_all_sky_albedo.pdf')          
plt.show()
plt.close()