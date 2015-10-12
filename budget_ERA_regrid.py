#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 08:45:15 2015

@author: jhchen
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import string
import numpy as np
import os
os.system("cls")
casenm='WTPCTR_EC'
nt=121
nz=33
dirs='D:/MyPaper/PhD04/Cases/'
diro='D:/MyPaper/PhD04/Cases/ERA/FORCING/'
if casenm[0:3]=='ETP':
    area=casenm[0:3]   
    iy,im,jd=2012,5,20
if casenm[0:3]=='WTP':
    area=casenm[0:3]    
    iy,im,jd=2010,7,3    
if casenm[0:3]=='NPC':
    area=casenm[0:3]  
    iy,im,jd=2010,8,2
if casenm[0:3]=='PRD':
    area=casenm[0:3]   
    iy,im,jd=2012,4,1 
if casenm[0:3]=='MLY':
    area=casenm[0:4]    
    iy,im,jd=2010,6,24 
if casenm[0:3]=='NEC':
    area=casenm[0:3]    
    iy,im,jd=2012,7,6 
folds="/CTREC"+"%4d"%iy+"%2.2d"%im+"%2.2d"%jd+"/Simulation/"
datestr="%4d"%iy+"%2.2d"%im+"%2.2d"%jd+"_031d"
dirpic='D:/MyPaper/PhD04/Pics/'
dirin=dirs+area+folds
dirq12=dirin
nameq1q2=area+'_'+"%4d"%iy+"%2.2d"%im+"%2.2d"%jd +"_031d_regrid_ERA.99"
nameforcing=area+'_'+"%4d"%iy+"%2.2d"%im+"%2.2d"%jd +"_031d_lsforcing_regrid_ERA.37"
#################################
mpl.rcParams['ytick.labelsize'] = 24
mpl.rcParams['xtick.labelsize'] = 24
#################################
ydat_r=[ -50.000 ,    50.000 ,   164.286,    307.143,    478.571  ,  678.571 ,
      907.143 ,  1164.286,   1450.000,   1764.286 ,  2107.143,   2478.572 ,
      2878.572,   3307.143,  3764.286,  4250.000,   4764.286,   5307.143, 
      5878.571,   6478.571,   7107.143,  7764.286,  8450.000,  9164.285,  
      9907.143,  10678.570,  11478.570,  12307.143,  13164.285,  14050.000,
      14964.285,  15907.143,  16878.572,  17878.572,  18907.145,  19964.285,
      21050.000,  22164.285,  23307.145,  24478.572,  25678.572,  26907.145,
      28164.285,  29450.000,  30764.285,  32107.145,  33478.570,  34878.570,
      36307.141,  37764.285,  39250.000,  40750.000]
ydato=[]
for yd in ydat_r:
    ydato.append(yd*0.001) 
ydat=[]
"""
for i in range(0,nz):
    ydat.append(ydato[i]) 
"""
for i in range(0,nz):
    ydat.append(i*0.5) 
####### file 1
filenm=casenm+'_regrid_qabcr_All.txt'
fpath=dirin+filenm
onedim1=[]
linesplit=[]
f=open(fpath)
ff=f.readlines()  ## first line in obs file is legend 
for line in ff:
    line=string.lstrip(line)
    linesplit.append(line[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            onedim1.append(string.atof(strs))
nl=len(onedim1)
allnms=["Temperature", "Water vapor mixing ratio","Horizontal velocity","Vertical velocity",
        "Cloud water mixing ratio","Rain water mixing ratio","Type A ice water mixing ratio",
        "Type B ice water mixing ratio","Relative humidity","fqv","ftc","rat","tc"]
allunts=[r'K', r'$g kg^{-1}$',r'$m s^{-1}$',r'$m s^{-1}$',
        r'$g kg^{-1}$',r'$g kg^{-1}$',r'$g kg^{-1}$',
        r'$g kg^{-1}$',r'$%$','fqv',"ftc","rat","C"]        
nvar=len(allnms)        
allvar=np.ndarray(shape=(nvar,nt,nz), dtype=float) # nvar 0-nvar t,q,u,w,qc,qr,qa,qb,rh,fqv,fqc,rat
for it in range(0,nt):
    for iz in range(0,nz):
        for ir in range(0,nvar):       
            k=it*nz*nvar+iz+nz*ir
            allvar[ir,it,iz]=onedim1[k]
#            allvar[ir,it,0]=0
del onedim1, linesplit
allvar_mean=np.ndarray(shape=(nvar,nz), dtype=float)
for iz in range(0,nz):
    for ir in range(0,nvar): 
        temp=0.0
        for it in range(1,nt-1):
            temp=temp+allvar[ir,it,iz]/(nt-2.) 
        allvar_mean[ir,iz]=temp   
del allvar
############### simulation Q1 Q2 phase change
filenm=casenm+'_micro_regrid_6hour.txt'
fpath=dirin+filenm
onedim2=[]
linesplit=[]
f=open(fpath)
ff=f.readlines()  ## first line in obs file is legend 
for line in ff:
    line=string.lstrip(line)
    linesplit.append(line[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            onedim2.append(string.atof(strs))
nl=len(onedim2)
micro=np.ndarray(shape=(5,nt,nz), dtype=float) # 1 condensation 2 evaporation 
#3 deposition 4 sublimation 5 fus freezing and melting
for it in range(0,nt):
    for iv in range(0,5):
        for iz in range(0,nz):
            izv=it*nz*5+iv*nz+iz
            micro[iv,it,iz]=onedim2[izv]
q1c=np.ndarray(shape=(nt,nz), dtype=float) # 1 condensation 2 evaporation
q2c=np.ndarray(shape=(nt,nz), dtype=float) # 1 condensation 2 evaporation
lv=2.5e6
lf=2.835e6
ls=2.835e6
cp=1850.
scq1=lv/cp
scq2=lf/cp
scq3=lf/cp
scq4=lv/cp
scq1=1.
scq2=1.
scq3=1.
scq4=1.
"""
scq1=0.25
scq2=0.25
scq3=0.25
scq4=0.25
"""
for it in range(0,nt):
    for iz in range(0,nz):           
        q1c[it,iz]=(micro[0,it,iz]+micro[1,it,iz])*scq1      \
                    +  micro[4,it,iz]*scq2                   \
                    + (micro[2,it,iz]+micro[3,it,iz])*scq3  
        q2c[it,iz]=( micro[0,it,iz]+micro[1,it,iz]      \
                    + (micro[2,it,iz]+micro[3,it,iz])*2.5e10/2.834e10 )*scq4 
    q1c[it,0]=0.0
    q2c[it,0]=0.0
q1cm=np.ndarray(shape=(nz), dtype=float) # 1 condensation 2 evaporation
micro_com=np.ndarray(shape=(5,nz), dtype=float) 
q2cm=np.ndarray(shape=(nz), dtype=float) # 1 condensation 2 evaporation
for iz in range(0,nz):
    temp1=0.0
    temp2=0.0
    for it in range(1,nt-1): #### abandon the firt and the last timestep 
        temp1=temp1+q1c[it,iz]/(nt-2.) 
        temp2=temp2+q2c[it,iz]/(nt-2.)
    q1cm[iz]=temp1
    q2cm[iz]=temp2
for iz in range(0,nz):    
    for iv in range(0,5):
        tmp3=0.0
        for it in range(1,nt-1):
            tmp3=tmp3+micro[iv,it,iz]/(nt-2)
        micro_com[iv,iz]= tmp3
q1cm2=np.ndarray(shape=(nt), dtype=float) # 1 condensation 2 evaporation
q2cm2=np.ndarray(shape=(nt), dtype=float) # 1 condensation 2 evaporation 
for it in range(0,nt):
    temp1=0.0
    temp2=0.0
    for iz in range(1,nz): 
        temp1=temp1+q1c[it,iz]
        temp2=temp2+q2c[it,iz]
    q1cm2[it]=temp1
    q2cm2[it]=temp2    
q1cm[0]=0.0
q2cm[0]=0.0
########## obs q1 q2
#########################################################            
####### file 2
filenm=casenm+'_regrid_eddydiffrad_All.txt'
fpath=dirin+filenm
onedim3=[]
linesplit=[]
f=open(fpath)
ff=f.readlines()  ## first line in obs file is legend 
for line in ff:
    line=string.lstrip(line)
    linesplit.append(line[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            onedim3.append(string.atof(strs))
nl=len(onedim3)
eddynms=["q1e", "q2e","q1d","q2d",
        "q1s","q2s","rlw",
        "rsw","fx","e1flux","e2flux","eflux"]
eddyunts=[ r'$K$ $d^{-1}$', r'$K$ $d^{-1}$',r'$K$ $d^{-1}$',r'$K$ $d^{-1}$',
        r'$K$ $d^{-1}$',r'$K$ $d^{-1}$',r'$K$ $d^{-1}$',
        r'$K$ $d^{-1}$',r'$m$ $m^{-1}$ $d^{-1}$',
        r'$K$ $kg^{-1}$ $m^2$ $s^{-1}$',
        r'$g$ $kg^{-1}$ $m^2$ $s^{-1}$',
        r'$kg$ $m^{-1}$ $s^2$']        
nvar=len(eddynms)        
eddyvar=np.ndarray(shape=(nvar,nt,nz), dtype=float) # nvar 0-nvar t,q,u,w,qc,qr,qa,qb,rh,fqv,fqc,rat
for it in range(0,nt):
    for iz in range(0,nz):
        for ir in range(0,nvar):       
            k=it*nz*nvar+iz+nz*ir
            eddyvar[ir,it,iz]=onedim3[k]
#del onedim1, linesplit
eddyvar_mean=np.ndarray(shape=(nvar,nz), dtype=float)
for iz in range(0,nz):
    for ir in range(0,nvar): 
        temp=0.0
        for it in range(1,nt-1): #### abandon the firt and the last timestep 
            temp=temp+eddyvar[ir,it,iz]/(nt-2.) 
        eddyvar_mean[ir,iz]=temp   
eddyvar_mean[:,0]=0.0
#eddyvar_mean[:,1]=0.0
#eddyvar_mean[:,2]=0.0
#eddyvar_mean[:,3]=0.0
#eddyvar_mean[0,1]=0.0
#eddyvar_mean[1,1]=0.0
#eddyvar_mean[0,2]=0.0
#eddyvar_mean[1,2]=0.0
#eddyvar_mean[0,3]=0.0
#eddyvar_mean[1,3]=0.0
#############################################################################3#
### sim q1q2
fpath=dirin+casenm+"_q1q2_regrid.txt"
onedim1=[]
linesplit=[]
f=open(fpath)
ff=f.readlines()  ## first line in obs file is legend 
for line in ff:
    line=string.lstrip(line)
    linesplit.append(line[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            onedim1.append(string.atof(strs))
q1simre=np.ndarray(shape=(nt,nz), dtype=float)
q2simre=np.ndarray(shape=(nt,nz), dtype=float)
for it in range(0,nt):
    for iz in range(0,nz):        
        k=it*(nz*2)+iz  ## the first record is time
        q1simre[it,iz]=onedim1[k]
        q2simre[it,iz]=onedim1[k+nz]       
q1simre_pf=np.ndarray(shape=(nz), dtype=float)
q2simre_pf=np.ndarray(shape=(nz), dtype=float)
for iz in range(0,nz):
    tmp1=0.0
    tmp2=0.0
    for it in range(1,nt-1):
        tmp1=tmp1+q1simre[it,iz] 
        tmp2=tmp2+q2simre[it,iz]
    q1simre_pf[iz]=tmp1/(nt-2)
    q2simre_pf[iz]=tmp2/(nt-2.)
q1simre_pf[0]=0.0
q2simre_pf[0]=0.0    
########## obs q1 q2
fpath=dirq12+nameq1q2
onedim1=[]
linesplit=[]
f=open(fpath)
ff=f.readlines()  ## first line in obs file is legend 
for line in ff:
    line=string.lstrip(line)
    linesplit.append(line[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            onedim1.append(string.atof(strs))
q1obs=np.ndarray(shape=(nt,nz), dtype=float)
q2obs=np.ndarray(shape=(nt,nz), dtype=float)
for it in range(0,nt):
    for iz in range(0,nz):        
        k=it*(nz*2)+iz  ## the first record is time
        q1obs[it,iz]=onedim1[k]
        q2obs[it,iz]=onedim1[k+nz]       
q1obs_pf=np.ndarray(shape=(nz), dtype=float)
q2obs_pf=np.ndarray(shape=(nz), dtype=float)
for iz in range(0,nz):
    tmp1=0.0
    tmp2=0.0
    for it in range(1,nt-1):
        tmp1=tmp1+q1obs[it,iz]/(nt-2.)
        tmp2=tmp2+q2obs[it,iz]/(nt-2.)
    q1obs_pf[iz]=tmp1
    q2obs_pf[iz]=tmp2
#########################################################
#### forcing
fpath=dirq12+nameforcing
onedim1=[]
linesplit=[]
f=open(fpath)
ff=f.readlines()  ## first line in obs file is legend 
for line in ff:
    line=string.lstrip(line)
    linesplit.append(line[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            onedim1.append(string.atof(strs))
fcq1=np.ndarray(shape=(nt,nz), dtype=float)
fcq2=np.ndarray(shape=(nt,nz), dtype=float)
scalefc=24*3600.
cp=1005.
hlat=2.5e6
for it in range(0,nt):
    for iz in range(0,nz):        
        k=it*(nz*2)+iz  ## the first record is time
        fcq1[it,iz]=onedim1[k]*scalefc
        fcq2[it,iz]=onedim1[k+nz]*scalefc*hlat/cp    ## convert to kg/kg to K per day   
fcq1obs_pf=np.ndarray(shape=(nz), dtype=float)
fcq2obs_pf=np.ndarray(shape=(nz), dtype=float)
for iz in range(0,nz):
    tmp1=0.0
    tmp2=0.0
    for it in range(1,nt-1):
        tmp1=tmp1+fcq1[it,iz]/(nt-2.)
        tmp2=tmp2+fcq2[it,iz]/(nt-2.)
    fcq1obs_pf[iz]=tmp1
    fcq2obs_pf[iz]=tmp2            
##### Plotting set up   
lnstycolor=['-','-','-','-','-']
lncolor=['orangered','darkgoldenrod','yellowgreen','deepskyblue','darkorchid']
lnmkcolor=['None','None','None','None','None'] 
lnwidcolor=[3.0,3.0,3.0,3.0,3.0]  
lnstygrey=['-','--',':','-',':']
lngrey=['silver','gray','darkgray','gainsboro','k']
lnmkgrey=['o','v','x','+','*']
lnwidgrey=[4.0,4.0,4.0,4.0,4.0]   
colors=lncolor
sty=lnstycolor
mker=lnmkcolor
width=lnwidcolor  
fig,ax0 = plt.subplots(nrows=1,ncols=1,figsize=(8,15))
#allvar_mean[4,0]=0.
allvar_mean[4,0]=0.
ax0.plot(allvar_mean[4,:],ydat,label=allnms[4],
    c=colors[0],ls=sty[0],marker=mker[0],lw=width[0],) #qc
#allvar_mean[5,0]=0.
ax0.plot(allvar_mean[5,:],ydat,label=allnms[5],
    c=colors[1],ls=sty[1],marker=mker[1],lw=width[1],)  #qr
#allvar_mean[6,1]=0.
allvar_mean[6,0]=0.
ax0.plot(allvar_mean[6,:],ydat,label=allnms[6],
    c=colors[2],ls=sty[2],marker=mker[2],lw=width[2],)  #qa

#allvar_mean[7,1]=0.
allvar_mean[7,0]=0.
ax0.plot(allvar_mean[7,:],ydat,label=allnms[7],
    c=colors[3],ls=sty[3],marker=mker[3],lw=width[3],)   #qb
totalwater=allvar_mean[4,:]+allvar_mean[5,:]+allvar_mean[6,:]+allvar_mean[7,:]
ax0.plot(totalwater,ydat,label='Total water content',
    c=colors[4],ls=sty[4],marker=mker[4],lw=width[4],)   #total water
ax0.set_title('Case '+casenm+' Water content profiles'+ r' $g$ $kg^{-1}$',fontsize=24)
ylabs='Height'+r' $km$'
ax0.set_ylabel(ylabs,fontsize=30)
ax0.legend()
#
plt.ylim(0,16)
plt.xlim(0,0.08)
xlb=range(0,8,2)
xlbs=[]
for xx in xlb:
    xlbs.append(xx/100.)
ax0.set_xticks(xlbs)
plt.show()                     
plt.savefig(dirpic+casenm+'_qaqbqcqr_profiles_regrid.pdf')          
#plt.show()
plt.close() 
##############################
# Q1 and Q2
###################
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['xtick.labelsize'] = 16
# q1q2
#eddyvar_mean[1,0]=0.0
#eddyvar_mean[1,1]=0.0
#eddyvar_mean[1,2]=0.0
q1sim_pf=np.ndarray(shape=(nz), dtype=float)
q2sim_pf=np.ndarray(shape=(nz), dtype=float)
q1sim_pf[:]=eddyvar_mean[0,:]+eddyvar_mean[2,:]+  q1cm[:] \
            +eddyvar_mean[6,:]+eddyvar_mean[7,:] +eddyvar_mean[4,:]
q2sim_pf[:]=eddyvar_mean[1,:]+eddyvar_mean[3,:]+ q2cm[:] +eddyvar_mean[5,:] 
q1sim_pf[0]=0.0
#q1sim_pf[0]=0.0
#q1sim_pf[1]=0.0
#q1sim_pf[1]=0.0
q2sim_pf[0]=0.0
#q2sim_pf[0]=0.0
#q2sim_pf[1]=0.0
#q2sim_pf[1]=0.0
####plot setup 
lnstycolor=['-',':','-',':']
lncolor=['orangered','orangered','yellowgreen','yellowgreen']
lncolor=['r','darkgoldenrod','g','b','darkorchid']
lncolor=['red','red','yellowgreen','yellowgreen']
lnmkcolor=['None','None','None','None','None'] 
lnwidcolor=[3.0,3.0,3.0,3.0,3.0]  
lnstygrey=['-',':','-',':']
lngrey=['silver','silver','darkgray','darkgray']
lnmkgrey=['o','x','o','x']
lnwidgrey=[4.0,4.0,4.0,4.0,4.0]   
colors=lncolor
sty=lnstycolor
mker=lnmkcolor
width=lnwidcolor 
size_title=16     
fig,([ax0,ax1,ax2],[ax3,ax4,ax5]) = plt.subplots(nrows=2,ncols=3,figsize=(9,9)) 
ax0.set_ylim(0,16)           
ax0.plot(q1sim_pf,ydat,label=r'Simulated $Q_1$',
    c=colors[0],ls=sty[0],marker=mker[0],lw=width[0],) #q1sim
#allvar_mean[5,0]=0.
ax0.plot(q2sim_pf,ydat,label=r'Simulated $Q_2$',
    c=colors[1],ls=sty[1],marker=mker[1],lw=width[1],)  #q1sim
ax0.plot(q1obs_pf,ydat,label=r'Obsevation $Q_1$',
    c=colors[2],ls=sty[2],marker=mker[2],lw=width[2],)  #q1obs
ax0.plot(q2obs_pf,ydat,label=r'Obsevation $Q_2$',
    c=colors[3],ls=sty[3],marker=mker[3],lw=width[3],)   #q2obs
#ax0.set_title('Case '+casenm+r'   $Q_1$ and $Q_2$'+ r' ($K$ $d^{-1}$)',fontsize=size_title)
ax0.set_title(r'$Q_1$ and $Q_2$'+ r' ($K$ $d^{-1}$)',fontsize=size_title)
ylabs='Height'+r' ($km$)'
ax0.set_ylabel(ylabs,fontsize=size_title)
#ax0.legend()
############forcing
ax1.set_ylim(0,16)           
ax1.plot(q1simre_pf,ydat,label=r'Simulated $Q_1$',
    c=colors[0],ls=sty[0],marker=mker[0],lw=width[0],) #q1sim
#allvar_mean[5,0]=0.
ax1.plot(q2simre_pf,ydat,label=r'Simulated $Q_2$',
    c=colors[1],ls=sty[1],marker=mker[1],lw=width[1],)  #q1sim
ax1.plot(q1obs_pf,ydat,label=r'Obsevation $Q_1$',
    c=colors[2],ls=sty[2],marker=mker[2],lw=width[2],)  #q1obs
ax1.plot(q2obs_pf,ydat,label=r'Obsevation $Q_2$',
    c=colors[3],ls=sty[3],marker=mker[3],lw=width[3],)   #q2obs
#ax1.set_title('Case '+casenm+r'   $Q_1$ and $Q_2$'+ r' ($K$ $d^{-1}$)',fontsize=size_title)
ax1.set_title(r'$Q_1$ and $Q_2$'+ r' ($K$ $d^{-1}$)',fontsize=size_title)
ylabs='Height'+r' ($km$)'
ax1.set_ylabel(ylabs,fontsize=size_title)
#ax1.legend()
###
lnstycolor=['-','-','-','-','-']
lncolor=['orangered','darkgoldenrod','yellowgreen','deepskyblue','darkorchid']
#lncolor=['r','darkgoldenrod','g','b','darkorchid']
#lncolor=['r','darkgoldenrod','g','b','darkorchid']
lnmkcolor=['None','None','None','None','None'] 
lnwidcolor=[3.0,3.0,3.0,3.0,3.0]  
lnstygrey=['-','--',':','-',':']
lngrey=['silver','gray','darkgray','gainsboro','k']
lnmkgrey=['o','v','x','+','*']
lnwidgrey=[4.0,4.0,4.0,4.0,4.0]     
colors=lncolor
sty=lnstycolor
mker=lnmkcolor
width=lnwidcolor 
ax2.set_ylim(0,16)           
ax2.plot(eddyvar_mean[0,:],ydat,label=r'Simulated $Q_1$e',
    c=colors[0],ls=sty[0],marker=mker[0],lw=width[0],) #q1sim
#allvar_mean[5,0]=0.
ax2.plot(eddyvar_mean[2,:]+eddyvar_mean[4,:],ydat,label=r'Simulated $Q_1$d and $Q_1$s',
    c=colors[1],ls=sty[1],marker=mker[1],lw=width[1],)  #q1sim
ax2.plot(q1cm,ydat,label=r'Simulated $Q_1$c',
    c=colors[2],ls=sty[2],marker=mker[2],lw=width[2],)  #q1obs
radsim=eddyvar_mean[6,:]+eddyvar_mean[7,:]    
ax2.plot(radsim,ydat,label=r'Simulated radiation',
    c=colors[3],ls=sty[3],marker=mker[3],lw=width[3],)   #q2obs
#ax2.set_title('Case '+casenm+r'  Terms of $Q_1$'+ r' ($K$ $d^{-1}$)',fontsize=size_title)
ax2.set_title(r'Terms of $Q_1$'+ r' ($K$ $d^{-1}$)',fontsize=size_title)
ylabs='Height'+r' ($km$)'
ax2.set_ylabel(ylabs,fontsize=size_title)
ax2.legend()
#
ax3.set_ylim(0,16)           
ax3.plot(eddyvar_mean[1,:],ydat,label=r'Simulated $Q_2$e',
    c=colors[0],ls=sty[0],marker=mker[0],lw=width[0],) #q1sim
#allvar_mean[5,0]=0.
ax3.plot(eddyvar_mean[3,:]+eddyvar_mean[5,:],ydat,label=r'Simulated $Q_2$d and $Q_2$s ',
    c=colors[1],ls=sty[1],marker=mker[1],lw=width[1],)  #q1sim
ax3.plot(q2cm,ydat,label=r'Simulated $Q_2$c',
    c=colors[2],ls=sty[2],marker=mker[2],lw=width[2],)  #q1obs
#ax3.set_title('Case '+casenm+r'  Terms of $Q_2$'+ r' ($K$ $d^{-1}$)',fontsize=size_title)
ax3.set_title(r'Terms of $Q_2$'+ r' ($K$ $d^{-1}$)',fontsize=size_title)
ylabs='Height'+r' ($km$)'
ax3.set_ylabel(ylabs,fontsize=size_title)
ax3.legend()
#
ax4.set_ylim(0,16)           
ax4.plot(eddyvar_mean[8,:],ydat,label=r'Simulated fx',
    c=colors[0],ls=sty[0],marker=mker[0],lw=width[0],) #fx
#allvar_mean[5,0]=0.
ax4.plot(eddyvar_mean[9,:],ydat,label=r'Simulated e1flux',
    c=colors[1],ls=sty[1],marker=mker[1],lw=width[1],)  #q1sim
ax4.plot(eddyvar_mean[10,:],ydat,label=r'Simulated e2flux',
    c=colors[2],ls=sty[2],marker=mker[2],lw=width[2],)  #q1obs
ax4.plot(eddyvar_mean[11,:],ydat,label=r'Simulated eflux',
    c=colors[3],ls=sty[2],marker=mker[2],lw=width[2],)  #q1obs
ax4.set_title('Case '+casenm+r'  simulated   ',fontsize=size_title)
ylabs='Height'+r' ($km$)'
ax4.set_ylabel(ylabs,fontsize=size_title)
ax4.legend()
#######
ax5.axis('off')
#####
fig.subplots_adjust(hspace=0.4)
plt.title(casenm,fontsize=size_title)
plt.show()                   
plt.savefig(dirpic+casenm+'_Q1Q2_profiles_regrid.pdf')          
#plt.show()
plt.close() 
################### plotting every simulated copms of Qc 
sty2=['-','-','-','-','-','-']
colors2=['yellowgreen','orangered','cyan','darkgoldenrod','deepskyblue','darkorchid']
colors2=['g','r','cyan','darkgoldenrod','b','darkorchid']
mker2=['None','None','None','None','None','None'] 
width2=[3.0,3.0,3.0,3.0,3.0,3.0]   
fig,([ax0,ax1,ax2],[ax3,ax4,ax5]) = plt.subplots(nrows=2,ncols=3,figsize=(9,9)) 
ax0.set_ylim(0,16)           
ax0.plot(micro_com[0,:],ydat,   #label=r'Simulated $Q_1$',
    c=colors2[0],ls=sty2[0],marker=mker2[0],lw=width2[0],) #q1sim
ax0.set_title('Case '+casenm+r' Condensation'+ r' ($K$ $d^{-1}$)',fontsize=size_title)
ylabs='Height'+r' ($km$)'
ax0.set_ylabel(ylabs,fontsize=size_title)
##############################
ax1.set_ylim(0,16)           
ax1.plot(micro_com[1,:],ydat,   #label=r'Simulated $Q_1$',
    c=colors2[1],ls=sty2[1],marker=mker2[1],lw=width2[1],) #q1sim
ax1.set_title('Case '+casenm+r' Evaporation'+ r' ($K$ $d^{-1}$)',fontsize=size_title)
ylabs='Height'+r' ($km$)'
ax1.set_ylabel(ylabs,fontsize=size_title)
##############################
ax2.set_ylim(0,16)           
ax2.plot(micro_com[2,:],ydat,   #label=r'Simulated $Q_1$',
    c=colors2[2],ls=sty2[2],marker=mker2[2],lw=width2[2],) #q1sim
ax2.set_title('Case '+casenm+r' Deposition'+ r' ($K$ $d^{-1}$)',fontsize=size_title)
ylabs='Height'+r' ($km$)'
ax2.set_ylabel(ylabs,fontsize=size_title)
##############################
ax3.set_ylim(0,16)           
ax3.plot(micro_com[3,:],ydat,   #label=r'Simulated $Q_1$',
    c=colors2[3],ls=sty2[3],marker=mker2[3],lw=width2[3],) #q1sim
ax3.set_title('Case '+casenm+r' Sublimation'+ r' ($K$ $d^{-1}$)',fontsize=size_title)
ylabs='Height'+r' ($km$)'
ax3.set_ylabel(ylabs,fontsize=size_title)
##############################
ax4.set_ylim(0,16)           
ax4.plot(micro_com[4,:],ydat,   #label=r'Simulated $Q_1$',
    c=colors2[4],ls=sty2[4],marker=mker2[4],lw=width2[4],) #q1sim
ax4.set_title('Case '+casenm+r' Fustion'+ r' ($K$ $d^{-1}$)',fontsize=size_title)
ylabs='Height'+r' ($km$)'
ax4.set_ylabel(ylabs,fontsize=size_title)
ax5.axis('off')
fig.subplots_adjust(hspace=0.4)
plt.title(casenm,fontsize=size_title)
plt.show()                 
plt.savefig(dirpic+casenm+'_Q1Q2_profiles_Composition_regrid.pdf')          
#plt.show()
plt.close() 
         