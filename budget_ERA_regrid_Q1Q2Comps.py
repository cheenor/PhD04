#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 20:59:11 2015

@author: jhchen
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import string
import numpy as np
from pylab import *
import os
os.system("cls")
#################################
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
#################################
def readAscii(fpath,iskp,*nl):
    #iskp  the total line skipped of the file
    # fpath   the full path of the file
    # usage: onedim=readAscii(fpaht,iskp)
    onedim=[]
    linesplit=[]
    f=open(fpath)
    if nl:
        nrl=nl[0]
        ff=f.readlines()[iskp:nrl]  ## first line in obs file is legend 
    else:
        ff=f.readlines()[iskp:]
    for line in ff:
        line=string.lstrip(line)
        linesplit.append(line[:-1].split(' '))
    for lnstrs in linesplit:
        for strs in lnstrs:
#            if strs!='':
#                onedim.append(string.atof(strs))
            try:
                onedim.append(string.atof(strs))
            except:
                strs=""
#            else:
#                onedim.append(string.atof(strs))
    del linesplit,ff
    f.close()
    print len(onedim)
    return onedim
#########################-----------------------------------------------------
nt=121
nz=33
CASE=['PRDCTR_EC','MLYRCTR_EC','NPCCTR_EC','NECCTR_EC','WTPCTR_EC','ETPCTR_EC']
nx=202
nga=len(CASE)
dirs='D:/MyPaper/PhD04/Cases/'
diro='D:/MyPaper/PhD04/Cases/ERA/FORCING/'
dirout='D:/MyPaper/PhD04/Pics/'
fig,ax=plt.subplots(nrows=3,ncols=4,figsize=(12,35))
#fig,axs=plt.subplots(nrows=2,ncols=3,figsize=(12,12))
color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
wd=[2,2,2,2,2]
iro=0
ic=0
for iga in range(0,nga):
    if ic==4:
        ic=0
        iro=iro+1
    print iro,ic
    casenm=CASE[iga]
    if casenm[0:3]=='ETP':
        area=casenm[0:3]   
        iy,im,jd=2012,5,20
        atr=r'$(f)$'
    if casenm[0:3]=='WTP':
        area=casenm[0:3]    
        iy,im,jd=2010,7,14
        atr=r'$(e)$'    
    if casenm[0:3]=='NPC':
        area=casenm[0:3]   
        iy,im,jd=2010,8,2
        atr=r'$(c)$'
    if casenm[0:3]=='PRD':
        area=casenm[0:3]   
        iy,im,jd=2012,4,1
        atr=r'$(a)$'
    if casenm[0:3]=='MLY':
        area=casenm[0:4]
        iy,im,jd=2010,6,24 
        atr=r'$(b)$'
    if casenm[0:3]=='NEC':
        area=casenm[0:3]   
        iy,im,jd=2012,7,6 
        atr=r'$(d)$'
    folds="/CTREC"+"%4d"%iy+"%2.2d"%im+"%2.2d"%jd+"/Simulation/"
    datestr="%4d"%iy+"%2.2d"%im+"%2.2d"%jd+'_031d'
    dirpic='D:/MyPaper/PhD04/Pics/'
    dirin=dirs+area+folds
    dirq12=dirin
    nameq1q2=area+'_'+"%4d"%iy+"%2.2d"%im+"%2.2d"%jd +"_031d_regrid_ERA.99"
    nameforcing=area+'_'+"%4d"%iy+"%2.2d"%im+"%2.2d"%jd +"_031d_lsforcing_regrid_ERA.37"
    ydat=[]
    for i in range(0,nz):
        ydat.append(i*0.5) 
####### file 1
    filenm=casenm+'_regrid_qabcr_All.txt'
    fpath=dirin+filenm
    iskp=0
    onedim1=readAscii(fpath, iskp)
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
    del onedim1
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
    iskp=0
    onedim2=readAscii(fpath, iskp)
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
    iskp=0
    onedim3=readAscii(fpath, iskp)
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
    iskp=0
    onedim1=readAscii(fpath,iskp)
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
    iskp=0
    onedim1=readAscii(fpath, iskp)
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
    iskp=0
    onedim1=readAscii(fpath, iskp)    
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
##############################
# Q1 and Q2
###################
    q1sim_pf=np.ndarray(shape=(nz), dtype=float)
    q2sim_pf=np.ndarray(shape=(nz), dtype=float)
    q1sim_pf[:]=eddyvar_mean[0,:]+eddyvar_mean[2,:]+  q1cm[:] \
            +eddyvar_mean[6,:]+eddyvar_mean[7,:] +eddyvar_mean[4,:]
    q2sim_pf[:]=eddyvar_mean[1,:]+eddyvar_mean[3,:]+ q2cm[:] +eddyvar_mean[5,:] 
    q1sim_pf[0]=0.0
    q2sim_pf[0]=0.0
    ####plot setup 
    lnstycolor=['-','-','-','-']
    lncolor=['orangered','orangered','yellowgreen','yellowgreen']
    lncolor=['r','darkgoldenrod','g','b','darkorchid']
<<<<<<< HEAD
    lncolor=['deeppink','b','green','y']
    lnmkcolor=['None','None','None','None','None'] 
    lnwidcolor=[3.0,3.0,3.0,3.0,3.0]  
=======
    lncolor=['r','g','b','orange','blueviolet','aqua']
    lnmkcolor=['None','None','None','None','None'] 
    lnwidcolor=[4.0,4.0,4.0,4.0,4.0]  
<<<<<<< HEAD
>>>>>>> e7f6294ce64f9ff8e82dba507be001724e7f2df1
    lnstygrey=['-','-','-','-']
    lngrey=['silver','silver','darkgray','darkgray']
=======
    lnstygrey=['-',':','-','--'] 
    lngrey=['k','k','darkgray','darkgray']
>>>>>>> updated 20151022
    lnmkgrey=['o','x','o','x']
    lnwidgrey=[3.0,3.0,4.0,4.0,4.0]   
    colors=lngrey
    sty=lnstygrey
    mker=lnmkcolor
    width=lnwidgrey
    size_title=18     
<<<<<<< HEAD
    ax[iro,ic].set_ylim(0,16)           
=======
    ax[iro,ic].set_ylim(0,16)
    ax[iro,ic].set_xlim(-6,6)           
>>>>>>> e7f6294ce64f9ff8e82dba507be001724e7f2df1
    ax[iro,ic].plot(eddyvar_mean[0,:],ydat,label=r'$Q_1$e',
        color=colors[0],ls=sty[0],marker=mker[0],lw=width[0],) #
    #allvar_mean[5,0]=0.
    ax[iro,ic].plot(eddyvar_mean[2,:]+eddyvar_mean[4,:],ydat,label=r'$Q_1$d'+'\n'+'+ $Q_1$s',
        color=colors[1],ls=sty[1],marker=mker[1],lw=width[1],)  #
    ax[iro,ic].plot(q1cm,ydat,label=r'$Q_1$c',
        color=colors[2],ls=sty[2],marker=mker[2],lw=width[2],)  #
    radsim=eddyvar_mean[6,:]+eddyvar_mean[7,:]
    ax[iro,ic].plot(radsim,ydat,label=r'$Q_R$',
        color=colors[3],ls=sty[3],marker=mker[3],lw=width[3],)   #q
    #ax[ir,ic].set_title('Case '+casenm+r'   $Q_1$ and $Q_2$'+ r' ($K$ $d^{-1}$)',fontsize=size_title)
    titlestr=atr+" "+area+r' $Q_1$'# ($K$ $day^{-1}$)'
    ax[iro,ic].set_title(titlestr,fontsize=size_title)
<<<<<<< HEAD
    xmajorLocator   = MultipleLocator(2) #将y轴主刻度标签设置为2的倍数  
=======
    xmajorLocator   = MultipleLocator(3) #将y轴主刻度标签设置为2的倍数  
>>>>>>> e7f6294ce64f9ff8e82dba507be001724e7f2df1
#    ymajorFormatter = FormatStrFormatter('%1.1f') #设置y轴标签文本的格式 
    ax[iro,ic].xaxis.set_major_locator(xmajorLocator) 
    ymajorLocator   = MultipleLocator(4) 
    ax[iro,ic].yaxis.set_major_locator(ymajorLocator)
    if iro==1 and ic==3 :
        ax[iro,ic].legend(loc=(0.97,0.2),frameon=False)
    if iro==0 and ic==2 :
        ax[iro,ic].legend(loc=(2.17,0.1),frameon=False)
    if ic==0:
        ylabs='Height'+r' ($km$)'
        ax[iro,ic].set_ylabel(ylabs,fontsize=size_title)
    if iro in(0,1,2) and ic in(1,2,3):
            setp(ax[iro,ic].get_yticklabels(), visible=False) #
    # Q2
    ic=ic+1
    lnstycolor=['-','-','-','-']
<<<<<<< HEAD
    lncolor=['r','orange','lime','y']
    lnmkcolor=['None','None','None','None','None'] 
    lnwidcolor=[3.0,3.0,3.0,3.0,3.0]  
=======
    lncolor=['r','g','b','k']
    lnmkcolor=['None','None','None','None','None'] 
    lnwidcolor=[4.0,4.0,4.0,4.0,4.0]  
<<<<<<< HEAD
>>>>>>> e7f6294ce64f9ff8e82dba507be001724e7f2df1
    lnstygrey=['-','-','-','-']
    lngrey=['silver','silver','darkgray','darkgray']
=======
    lnstygrey=['-',':','-',':'] 
    lngrey=['k','k','darkgray','darkgray']
>>>>>>> updated 20151022
    lnmkgrey=['o','x','o','x']
    lnwidgrey=[4.0,4.0,4.0,4.0,4.0]   
    colors=lngrey
    sty=lnstygrey
    mker=lnmkcolor
<<<<<<< HEAD
    width=lnwidcolor 
<<<<<<< HEAD
    ax[iro,ic].set_ylim(0,16)           
=======
=======
    width=lnwidgrey
>>>>>>> updated 20151022
    ax[iro,ic].set_ylim(0,16) 
    ax[iro,ic].set_xlim(-6,6)           
>>>>>>> e7f6294ce64f9ff8e82dba507be001724e7f2df1
    ax[iro,ic].plot(eddyvar_mean[1,:],ydat,label=r'$Q_2$e',
        color=colors[0],ls=sty[0],marker=mker[0],lw=width[0],) #
    #allvar_mean[5,0]=0.
    ax[iro,ic].plot(eddyvar_mean[3,:]+eddyvar_mean[5,:],ydat,label=r'$Q_2$d'+'\n'+'+ $Q_2$s',
        color=colors[1],ls=sty[1],marker=mker[1],lw=width[1],)  #
    ax[iro,ic].plot(q2cm,ydat,label=r'$Q_2$c',
        color=colors[2],ls=sty[2],marker=mker[2],lw=width[2],)  #
    #ax[ir,ic].set_title('Case '+casenm+r'   $Q_1$ and $Q_2$'+ r' ($K$ $d^{-1}$)',fontsize=size_title)
    titlestr=atr+" "+area+r' $Q_2$' # ($K$ $day^{-1}$)'
    ax[iro,ic].set_title(titlestr,fontsize=size_title)
<<<<<<< HEAD
    xmajorLocator   = MultipleLocator(2) #将y轴主刻度标签设置为2的倍数  
=======
    xmajorLocator   = MultipleLocator(3) #将y轴主刻度标签设置为2的倍数  
>>>>>>> e7f6294ce64f9ff8e82dba507be001724e7f2df1
#    ymajorFormatter = FormatStrFormatter('%1.1f') #设置y轴标签文本的格式 
    ax[iro,ic].xaxis.set_major_locator(xmajorLocator) 
    ymajorLocator   = MultipleLocator(4) 
    ax[iro,ic].yaxis.set_major_locator(ymajorLocator) 
    if iro==1 and ic==3 :
        ax[iro,ic].legend(loc=(0.97,0.2),frameon=False)
    if iro==0 and ic==2 :
        ax[iro,ic].legend(loc=(2.17,0.1),frameon=False)
    if ic==0:
        ylabs='Height'+r' ($km$)'
        ax[iro,ic].set_ylabel(ylabs,fontsize=size_title)
    if iro in(0,1,2) and ic in(1,2,3):
            setp(ax[iro,ic].get_yticklabels(), visible=False) #
    ic=ic+1
plt.subplots_adjust(left = 0.1, wspace = 0.2, hspace = 0.3, \
    bottom = 0.1, top = 0.90)
plt.show()                     
plt.savefig(dirout+'ALLCASE_Q1Q2Comps_Gray.png',dpi=300)          
plt.show()
plt.close()
    #
    