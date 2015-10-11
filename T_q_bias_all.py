#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wen Aug 05 11:16:43 2015

@author: jhchen
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import string
import numpy as np
import datetime
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['contour.negative_linestyle'] = 'dashed'
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['xtick.labelsize'] = 16
plt.rc('lines', linewidth=4)
def readAscii(fpath,iskp):
    #iskp  the total line skipped of the file
    # fpath   the full path of the file
    # usage: onedim=readAscii(fpaht,iskp)
    onedim=[]
    linesplit=[]
    f=open(fpath)
    ff=f.readlines()[iskp:]  ## first line in obs file is legend 
    for line in ff:
        line=string.lstrip(line)
        linesplit.append(line[:-1].split(' '))
    for lnstrs in linesplit:
        for strs in lnstrs:
            if strs!='':
                onedim.append(string.atof(strs))
    del linesplit
    f.close()
    return onedim
#
nt=121
nday=31
nz=52
CASE=['PRDCTR_EC','MLYRCTR_EC','NPCCTR_EC','NECCTR_EC','WTPCTR_EC','ETPCTR_EC']
nx=202
nga=len(CASE)
diro='D:/MyPaper/PhD04/Cases/ERA/FORCING/'
dirin1='D:/MyPaper/PhD04/Cases/'
dirs=dirin1
dircnrain='D:/MyPaper/PhD04/Data/RainCN05/'
pic_out='D:/MyPaper/PhD04/Pics/'
dirpic=pic_out
#
ydat_r=[ -50.000 ,    50.000 ,   164.286,    307.143,    478.571  ,  678.571 ,
      907.143 ,  1164.286,   1450.000,   1764.286 ,  2107.143,   2478.572 ,
      2878.572,   3307.143,  3764.286,  4250.000,   4764.286,   5307.143, 
      5878.571,   6478.571,   7107.143,  7764.286,  8450.000,  9164.285,  
      9907.143,  10678.570,  11478.570,  12307.143,  13164.285,  14050.000,
      14964.285,  15907.143,  16878.572,  17878.572,  18907.145,  19964.285,
      21050.000,  22164.285,  23307.145,  24478.572,  25678.572,  26907.145,
      28164.285,  29450.000,  30764.285,  32107.145,  33478.570,  34878.570,
      36307.141,  37764.285,  39250.000,  40750.000]
ydat=[]
for yd in ydat_r:
    ydat.append(yd*0.001)
smt=np.ndarray(shape=(nz,nt,nga),dtype=float)
smq=np.ndarray(shape=(nz,nt,nga),dtype=float)
smu=np.ndarray(shape=(nz,nt,nga),dtype=float)
smw=np.ndarray(shape=(nz,nt,nga),dtype=float)
smqc=np.ndarray(shape=(nz,nt,nga),dtype=float)
smqr=np.ndarray(shape=(nz,nt,nga),dtype=float)
smqa=np.ndarray(shape=(nz,nt,nga),dtype=float)
smqb=np.ndarray(shape=(nz,nt,nga),dtype=float)
smrh=np.ndarray(shape=(nz,nt,nga),dtype=float)
smfqv=np.ndarray(shape=(nz,nt,nga),dtype=float)
smftc=np.ndarray(shape=(nz,nt,nga),dtype=float)
smrat=np.ndarray(shape=(nz,nt,nga),dtype=float)
smtc=np.ndarray(shape=(nz,nt,nga),dtype=float)
smtco=np.ndarray(shape=(nz,nt,nga),dtype=float)
#
obstha=np.ndarray(shape=(nz,nt,nga),dtype=float)
obsqv=np.ndarray(shape=(nz,nt,nga),dtype=float)
obstmp=np.ndarray(shape=(nz,nt,nga),dtype=float)
#obsrh=np.ndarray(shape=(nz,nt,nga),dtype=float)
obsu=np.ndarray(shape=(nz,nt,nga),dtype=float)
obsv=np.ndarray(shape=(nz,nt,nga),dtype=float)
obsw=np.ndarray(shape=(nz,nt,nga),dtype=float)  
xdatet=[]
titlestr=[]
for iga in range(0,nga):
    casenm=CASE[iga]
    if casenm[0:3]=='ETP' :
        iy,im,jd=2010,6,3
#        nt,nday=125,31
        namestr=casenm[0:3]
        marktr=r"($f$)"
        yearstr="%d"%iy
    elif casenm[0:3]=='WTP' :
        iy,im,jd=2010,7,3
#        nt,nday=125,31
        namestr=casenm[0:3]
        marktr=r"($e$)"
        yearstr="%d"%iy
    elif casenm[0:3]=='PRD' :
        iy,im,jd=2012,4,1
#        nt,nday=125,31
        namestr=casenm[0:3]
        marktr=r"($a$)"
        yearstr="%d"%iy
    elif casenm[0:3]=='MLY' :
        iy,im,jd=2010,6,2
#        nt,nday=125,31
        namestr=casenm[0:4]
        marktr=r"($b$)"
        yearstr="%d"%iy
    elif casenm[0:3]=='NPC' :
        iy,im,jd=2010,8,2
#        nt,nday=125,31
        namestr=casenm[0:3]
        marktr=r"($c$)"
        yearstr="%d"%iy
    elif casenm[0:3]=='NEC' :
        iy,im,jd=2012,7,6
#        nt,nday=125,31
        namestr=casenm[0:3]
        marktr=r"($d$)"
        yearstr="%d"%iy
    topstr='' #'_250'
    datestr="%4d"%iy+"%2.2d"%im+"%2.2d"%jd+"_031d"
    strnm=namestr+'_'+datestr
    dirin=dirin1+namestr+'/'
    folds="/CTREC"+"%4d"%iy+"%2.2d"%im+"%2.2d"%jd+"/Simulation/"
    datestr="%4d"%iy+"%2.2d"%im+"%2.2d"%jd+"_031d"
    dirin=dirs+namestr+folds
    dirobs=diro+namestr+'/'
    f43=namestr+'_'+datestr+"_ERA.43"
    nameforcing=namestr+'_'+datestr+"_LSFORCING_ERA.37"
    titlestr.append(marktr+' '+namestr+' ('+yearstr+')')
#fpath=dirin+'micro_202_ETP2D3'
#fff=np.fromfile(fpath,dtype=float)
#fc=open(fpath,"rb")
#fff=fc.read()
#print fff[0],fff[1],fff[12],fff[78]
#print len(fff)
    datestart=datetime.datetime(iy,im,jd,0,0,0)
    det=datetime.timedelta(hours=6)            
    dateiso=[]            
    for dt in range(0,nt):
        dateiso.append(datestart+dt*det)
    xdate=[]    
    xdat=range(0,nt)            
    for tm in dateiso:
        xdate.append(datetime.datetime.strftime(tm,"%b/%d")) 
    xdatet.append(xdate)
###################################################################
    fpath=dirin+casenm+'_All.txt'
    iskp=0
    onedim1=readAscii(fpath, iskp)
    iskp=13*nz
    ikkk=0
    for it in range(0,nt):
        for iz in range(0,nz):
            ikkk=it*iskp+nz*0+iz
            print ikkk,it,iskp,nz,iz
            smt[iz,it,iga]=onedim1[ikkk]
            ikkk=it*iskp+nz*1+iz
            smq[iz,it,iga]=onedim1[ikkk]
            ikkk=it*iskp+nz*2+iz
            smu[iz,it,iga]=onedim1[ikkk]
            ikkk=it*iskp+nz*3+iz
            smw[iz,it,iga]=onedim1[ikkk]
            ikkk=it*iskp+nz*4+iz
            smqc[iz,it,iga]=onedim1[ikkk]
            ikkk=it*iskp+nz*5+iz
            smqr[iz,it,iga]=onedim1[ikkk]
            ikkk=it*iskp+nz*6+iz
            smqa[iz,it,iga]=onedim1[ikkk]
            ikkk=it*iskp+nz*7+iz
            smqb[iz,it,iga]=onedim1[ikkk]
            ikkk=it*iskp+nz*8+iz
            smrh[iz,it,iga]=onedim1[ikkk]
            ikkk=it*iskp+nz*9+iz
            smfqv[iz,it,iga]=onedim1[ikkk]
            ikkk=it*iskp+nz*10+iz
            smftc[iz,it,iga]=onedim1[ikkk]
            ikkk=it*iskp+nz*11+iz
            smrat[iz,it,iga]=onedim1[ikkk]
            ikkk=it*iskp+nz*12+iz
            smtc[iz,it,iga]=onedim1[ikkk]
###############################################################################
##  open obs  files
    del onedim1
    fpath=dirobs+f43
    iskp=0
    onedim1=readAscii(fpath, iskp)         
    iskp=6*nz+1
    for it in range(0,nt):
        for iz in range(0,nz):
            k=iskp*it+1+iz            
            obstha[iz,it,iga]=onedim1[k]
            k=iskp*it+1+iz+nz*1            
            obsqv[iz,it,iga]=onedim1[k]*1000. #convert kg/kg to g/kg  
            k=iskp*it+1+iz+nz*2         
            obstmp[iz,it,iga]=onedim1[k]  
            k=iskp*it+1+iz+nz*3          
            obsu[iz,it,iga]=onedim1[k]  
            k=iskp*it+1+iz+nz*4           
            obsv[iz,it,iga]=onedim1[k]  
            k=iskp*it+1+iz+nz*5          
            obsw[iz,it,iga]=onedim1[k]*3600./100. # convert pa/s to hpa/hr
            smtco[iz,it,iga]=(smt[iz,it,iga]-obstha[iz,it,iga])/\
                (obstha[iz,it,iga]/obstmp[iz,it,iga])
###############################################################################
fig,ax=plt.subplots(nrows=3,ncols=2,figsize=(18,12))
#fig,axs=plt.subplots(nrows=2,ncols=3,figsize=(12,12))
color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
wd=[2,2,2,2,2]
jc=0
jr=0
ij=1
for iga in range(0,nga):
    if jc==2:
        jc=0
        jr=jr+1
    mker=titlestr[iga]
    xdate=xdatet[iga]
    levs1=[-15,-12,-6,-3,4,8,12,15]
    colors1=['g','g','g','g','r','r','r','r']
    linetype1=['dotted','dotted','dotted','dotted','solid','solid','solid','solid'] 
    levs2=[-3,-1,1,2,3,4]
    colors2=['g','g','r','r','r','r']
    linetype2=['dotted','dotted','solid','solid','solid','solid']
    titlename=[r"Temperature Bias ($K$)",r"Water Vapor Mixing Ratio Bias ($g$ $kg^{-1}$) "]
    font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 16,
        }     
    plt.subplot(3,2,ij)
    zdat=smtc[:,:,iga]+273.16-obstmp[:,:,iga]
    #zdat=smtco   #-obstmp
    ##zdat[0,:]=0.0   ## the first level is below surface ground
    ax[jr,jc]=plt.contour(xdat,ydat,zdat,colors='r',
        linewidths=1.5,levels=levs1,linestyles=linetype1)                           
#    plt.title(titlename[0],fontsize=16)                          
    plt.axis([0, 121, 0, 16])
    plt.clabel(ax[jr,jc],inline=1,fmt='%1d',fontsize=12)
    zdat=smq[:,:,iga]-obsqv[:,:,iga]
    zdat[0,:]=0.0   ## the first level is below surface ground
    ax[jr,jc]=plt.contour(xdat,ydat,zdat,colors='g',
        linewidths=1.5,levels=levs2,linestyles=linetype2)  
    plt.axis([0, 121, 0, 16])
    plt.clabel(ax[jr,jc],inline=1,fmt='%1d',fontsize=12)
    axx=fig.add_subplot(3,2,ij)
    text1=mker  #r"($a$)"
    axx.text(1.5,16.5,text1,fontsize=18)                        
    axx.set_xticks(range(0,nt,16))
    xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
    axx.set_xticklabels(xticklabels, size=16)
    plt.ylabel(r'Height ($km$)', fontdict=font)
    jc=jc+1
    ij=ij+1                
plt.show()
fig.subplots_adjust(left=0.1,bottom=0.1,right=1-0.1,top=1-0.1,hspace=0.4)
plt.savefig(dirpic+"ALLCASE_T&qv_Bias.png",dpi=300)          
plt.show()
plt.close()
#--------------PROFILES -------------------------------------------------------
fig,ax=plt.subplots(nrows=2,ncols=6,figsize=(15,50))
#fig,axs=plt.subplots(nrows=2,ncols=3,figsize=(12,12))
#color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
#color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
#wd=[2,2,2,2,2]
jc=0
jr=0
ij=1
for iga in range(0,nga):
#    if jc==6:
#        jc=0
#        jr=jr+1
    jc=iga
    ij=jc+1
    mker=titlestr[iga]
    xdate=xdatet[iga]
#obstmp[iz,it]=onedim1[k] obsqv[iz,it]
#smtco[iz,it]   smq
    meanobs=np.ndarray(shape=(nz,2),dtype=float)
    obsstd=np.ndarray(shape=(nz,2),dtype=float)
    meansim=np.ndarray(shape=(nz,2),dtype=float)
    simstd=np.ndarray(shape=(nz,2),dtype=float)
    for k in range(0,nz):
        tmp1=0.
        tmp2=0.
        tmp3=0.
        tmp4=0.
        tmp=np.ndarray(shape=(nt,4),dtype=float)
        for it in range(0,nt):
            tmp1=tmp1+obstha[k,it,iga]
            tmp2=tmp2+obsqv[k,it,iga]
            tmp3=tmp3+smt[k,it,iga]
            tmp4=tmp4+smq[k,it,iga]
            tmp[it,0]=obstha[k,it,iga]
            tmp[it,1]=obsqv[k,it,iga]
            tmp[it,2]=smt[k,it,iga]
            tmp[it,3]=smq[k,it,iga]
        meanobs[k,0]=tmp1/nt
        meanobs[k,1]=tmp2/nt
        meansim[k,0]=tmp3/nt
        meansim[k,1]=tmp4/nt
        obsstd[k,0]=tmp[:,0].std()
        obsstd[k,1]=tmp[:,1].std()
        simstd[k,0]=tmp[:,2].std()
        simstd[k,1]=tmp[:,3].std()
#
    colors=['g','r']
    width=[1.5,1.5]
    ax0=plt.subplot(2,6,ij)
    ax0.errorbar(meanobs[1:32,0],ydat[1:32],label=r"$O$",
        c=colors[0],lw=width[0],xerr=obsstd[1:32,0]) #qc
    ax0.errorbar(meansim[1:32,0],ydat[1:32],label=r"$S$",
        c=colors[1],lw=width[1],xerr=simstd[1:32,0]) #qc
#    plt.ylabel(r'Height ($km$)', fontdict=font)
    #axx=fig.add_subplot(1,2,1) 
    #text1=r"($a$)"
    #axx.text(285,16.5,text1,fontsize=14)
    plt.xlim(300,500)                        
    ax0.set_xticks(range(300,500,80))
    ax0.set_title(mker,fontsize=18)
    #ax0.yaxis.limit_range_for_scale(0,16)
#    ax0.set_title(r'$(a)$  Potential Tempe. ($K$)')
    leftc,widthc=.25,.5
    bottomc,heightc=.25,.5
    rightc=leftc+widthc
    topc=bottomc+heightc
#    ax0.text(rightc,topc,r'$\theta$',horizontalalignment='right',
#        verticalalignment='top',fontsize=18)
    ax0.text(460,14.5,r'$\theta$',fontsize=20)
    axx=fig.add_subplot(2,6,ij)
    if ij==1 or ij==7:    
        plt.ylabel(r'Height ($km$)', fontdict=font)
        ax0.legend(loc=(0.4,0.7),frameon=False)
    else:
        for tick in axx.yaxis.get_major_ticks():
            tick.label1On = False    
    ij=ij+6
    ax1=plt.subplot(2,6,ij)
    ax1.errorbar(meanobs[1:32,1],ydat[1:32],label=r"$O$",
        c=colors[0],lw=width[0],xerr=obsstd[1:32,1]) #
    ax1.errorbar(meansim[1:32,1],ydat[1:32],label=r"$S$",
        c=colors[1],lw=width[1],xerr=simstd[1:32,1]) #qc
#plt.ylabel(r'Height ($km$)', fontdict=font)
#    ax1.set_title(r'($b$) Vapor Mixing Ratio ($g$ $kg^{-1}$)')
#    ax1.text(leftc,topc,r'$q$',horizontalalignment='left',
#        verticalalignment='top',fontsize=18)
    ax1.text(420.*23/500.,14.5,r'$q$',fontsize=20)
    ax1.set_title(mker,fontsize=18)
    plt.xlim(0,23)
    ax1.set_xticks(range(0,23,6))
    axx=fig.add_subplot(2,6,ij) 
    if ij==1 or ij==7:    
        plt.ylabel(r'Height ($km$)', fontdict=font)
        ax1.legend(loc=(0.4,0.7),frameon=False)
        
    else:
        for tick in axx.yaxis.get_major_ticks():
            tick.label1On = False 
#    ij=ij+1
#    jc=jc+1
plt.show()
plt.savefig(dirpic+"ALLCASE_obsVSsim_t_q_profs.png",dpi=300)          
plt.show()
plt.close()
#-------------------profiles of basis
fig,ax=plt.subplots(nrows=2,ncols=3,figsize=(12,21))
#fig,axs=plt.subplots(nrows=2,ncols=3,figsize=(12,12))
color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
wd=[2,2,2,2,2]
jc=0
jr=0
ij=1
for iga in range(0,nga):
    if jc==3:
        jc=0
        jr=jr+1
    bias=np.ndarray(shape=(nz,2),dtype=float)
    biasstd=np.ndarray(shape=(it,2),dtype=float)
    for k in range(0,nz):
        tmp1=0.
        tmp2=0.
        tmp=np.ndarray(shape=(nt,2),dtype=float)
        for it in range(0,nt):
            tmp1=tmp1+smtc[k,it,iga]+273.16-obstmp[k,it,iga]
            tmp2=tmp2+smq[k,it,iga]-obsqv[k,it,iga]
            tmp[it,0]=smtc[k,it,iga]+273.16-obstmp[k,it,iga]
            tmp[it,1]=smq[k,it,iga]-obsqv[k,it,iga]
        bias[k,0]=tmp1/nt
        bias[k,1]=tmp2/nt
        biasstd[k,0]=tmp[:,0].std()
        biasstd[k,1]=tmp[:,1].std()
#
    colors=['k','lightgrey']
    colors=['g','r']
    width=[2,2]
    #ax0=plt.subplot(1,2,1)
    ax[jr,jc].errorbar(bias[1:32,0],ydat[1:32],label=r"$T$",
        c=colors[1],lw=width[0],xerr=biasstd[1:32,0]) #qc
    ax[jr,jc].errorbar(bias[1:32,1],ydat[1:32],label=r"$q$",
        c=colors[0],lw=width[0],xerr=biasstd[1:32,1]) #qc  
    ax[jr,jc].plot([0,0],[0,16],c='k',lw=1.5)
    #ax[jr,jc].errorbar(meansim[1:32,0],ydat[1:32],label="SIM",
    #    c=colors[1],lw=width[1],xerr=simstd[1:32,0]) #qc
    strs=titlestr[iga]    
    ax[jr,jc].set_xlim(-12,10) 
    ax[jr,jc].set_title(strs,fontsize=18)
    if jr==1 and jc==2 :
        ax[jr,jc].legend(loc=(0.05,0.8),frameon=False)
    if jc==0:    
        ax[jr,jc].set_ylabel(r'Height ($km$)', fontdict=font)
#    ax[jr,jc].set_title(area)
    ij=ij+1
    jc=jc+1
plt.show()
plt.savefig(dirpic+"ALLCASE_BiasProfile_t_q.png",dpi=300)          
plt.show()
plt.close()
#----------the origin q and tmeperature
###############################################################################
fig,ax=plt.subplots(nrows=3,ncols=2,figsize=(18,12))
#fig,axs=plt.subplots(nrows=2,ncols=3,figsize=(12,12))
color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
wd=[2,2,2,2,2]
jc=0
jr=0
ij=1
for iga in range(0,nga):
    if jc==2:
        jc=0
        jr=jr+1
    mker=titlestr[iga]
    xdate=xdatet[iga]
    levs1=[-12,-9,-6,-3,3,6,9,12]
    colors1=['g','g','g','g','r','r','r','r']
    linetype1=['solid','solid','solid','solid','dotted','dotted','dotted','dotted'] 
    levs2=[3,5,9,12,15,18,21]
    colors2=['g','g','r','r','r','r']
    linetype2=['dotted','dotted','solid','solid','solid','solid']
    titlename=[r"Temperature Bias ($K$)",r"Water Vapor Mixing Ratio Bias ($g$ $kg^{-1}$) "]
    font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 16,
        }     
    plt.subplot(3,2,ij)
    zdat=obsw[:,:,iga]
    #zdat=smtco   #-obstmp
    ##zdat[0,:]=0.0   ## the first level is below surface ground
    ax[jr,jc]=plt.contour(xdat,ydat,zdat,colors='r',
        linewidths=1.5,levels=levs1,linestyles=linetype1)                           
#    plt.title(titlename[0],fontsize=16)                          
    plt.axis([0, 121, 0, 16])
    plt.clabel(ax[jr,jc],inline=1,fmt='%1d',fontsize=12)
    zdat=obsqv[:,:,iga]
    zdat[0,:]=0.0   ## the first level is below surface ground
    ax[jr,jc]=plt.contour(xdat,ydat,zdat,colors='g',
        linewidths=1.5,levels=levs2)#,linestyles=linetype2)  
    plt.axis([0, 121, 0, 16])
    plt.clabel(ax[jr,jc],inline=1,fmt='%1d',fontsize=12)
    axx=fig.add_subplot(3,2,ij)
    text1=mker  #r"($a$)"
    axx.text(1.5,16.5,text1,fontsize=18)                        
    axx.set_xticks(range(0,nt,16))
    xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
    axx.set_xticklabels(xticklabels, size=16)
    plt.ylabel(r'Height ($km$)', fontdict=font)
    jc=jc+1
    ij=ij+1                
plt.show()
fig.subplots_adjust(left=0.1,bottom=0.1,right=1-0.1,top=1-0.1,hspace=0.4)
plt.savefig(dirpic+"ALLCASE_EC_T&qv.png",dpi=300)          
plt.show()
plt.close()