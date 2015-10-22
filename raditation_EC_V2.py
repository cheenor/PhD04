#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 09:22:47 2015

@author: jhchen
"""
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
import calendar
import string
import numpy as np
from pylab import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter 
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
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
    return onedim
#
CASES=['PRDCTR_EC','MLYRCTR_EC','NPCCTR_EC','NECCTR_EC','WTPCTR_EC','ETPCTR_EC']
nt=241
dirs='D:/MyPaper/PhD04/Cases/'
diro='D:/MyPaper/PhD04/Cases/ERA/FORCING/'
<<<<<<< HEAD
dirpic='D:/MyPaper/PhD04/Pics/'
dirceres='D:/MyPaper/PhD04/Data/CERES/'
nga=len(CASES)
fig,axs = plt.subplots(nrows=2,ncols=3,figsize=(18,8))
color_cycle=['deeppink','lime','y','indigo','b','darkorange','cyan','k','grey']
=======
dirout='D:/MyPaper/PhD04/Pics/'
dirceres='D:/MyPaper/PhD04/Data/CERES/'
nga=len(CASES)
fig,axs = plt.subplots(nrows=2,ncols=3,figsize=(18,8))
color_cycle=['deeppink','lime','y','indigo','b','r','cyan','k','grey']
<<<<<<< HEAD
>>>>>>> e7f6294ce64f9ff8e82dba507be001724e7f2df1
=======
color_cycle=['deeppink','lime','y','indigo','k','0.65','cyan','k','0.75']
>>>>>>> updated 20151022
wds=[2,2,2,2,2,2,2,2,2]
linests=['-','--','-','--','-','--','-','--','-','--']
linests=['-','-','-','-','-','-','-','-','-','-']
ir=0
jc=0
for iga in range(0,nga):
    if jc==3:
        jc=0
        ir=ir+1
    casenm=CASES[iga]
    if casenm[0:3]=='ETP':
        area=casenm[0:3]    
        iy,im,jd=2010,6,3
        astr=r'$(f)$'
    if casenm[0:3]=='WTP':
        area=casenm[0:3]   
        iy,im,jd=2010,7,3 
        astr=r'$(e)$'  
    if casenm[0:3]=='NPC':
        area=casenm[0:3]    
        iy,im,jd=2010,8,2
        astr=r'$(c)$'
    if casenm[0:3]=='PRD':
        area=casenm[0:3]   
        iy,im,jd=2012,4,1 
        astr=r'$(a)$'
    if casenm[0:3]=='MLY':
        area=casenm[0:4]   
        iy,im,jd=2010,6,2 
        astr=r'$(b)$'
    if casenm[0:3]=='NEC':
        area=casenm[0:3]   
        iy,im,jd=2012,7,6 
        astr=r'$(d)$'
    folds="/CTREC"+"%4d"%iy+"%2.2d"%im+"%2.2d"%jd+"/Simulation/"
    datestr="%4d"%iy+"%2.2d"%im+"%2.2d"%jd+'_031d'  
    dirin=dirs+area+folds
    dirsim=dirin
    simnm='rad_3hr_'+casenm
    ceresnm=area+"_"+"%4d"%iy+'_'+"%2.2d"%im+'_'+"%2.2d"%jd+"__30d3h.txt"
    ceresnmtoa=area+"_"+"%4d"%iy+'_'+"%2.2d"%im+'_'+"%2.2d"%jd+"__30d3h_TOA_OBS.txt"
    #;;; set for x axis labels
    monstr="%02d"%(im)  ### number to string, 1 to 01, 10 to 10
    dnm=calendar.monthrange(iy,im)[1]  #calendar.monthrange(1997,7) 
    #                              #reture two index, sencond is the day number
    datestart=datetime.datetime(iy,im,jd,0,0,0)
    det=datetime.timedelta(hours=3)            
    dateiso=[]            
    for dt in range(0,nt):
        dateiso.append(datestart+dt*det)
    xdate=[]    
    xdat=range(0,nt)            
    for tm in dateiso:
<<<<<<< HEAD
        xdate.append(datetime.datetime.strftime(tm,"%b/%d")) 
=======
        xdate.append(datetime.datetime.strftime(tm,"%d/%b")) 
>>>>>>> e7f6294ce64f9ff8e82dba507be001724e7f2df1
###### fro daily
    det=datetime.timedelta(hours=24)            
    dateiso=[]            
    for dt in range(0,dnm):
        dateiso.append(datestart+dt*det)
    xdated=[]    
    xdatd=range(0,dnm)            
    for tm in dateiso:
<<<<<<< HEAD
        xdated.append(datetime.datetime.strftime(tm,"%b/%d"))
=======
        xdated.append(datetime.datetime.strftime(tm,"%d/%b"))
>>>>>>> e7f6294ce64f9ff8e82dba507be001724e7f2df1
#################
    det=datetime.timedelta(hours=6)            
    dateiso=[]            
    for dt in range(0,nt):
        dateiso.append(datestart+dt*det)
    xdatencep=[]    
    xdatncep=range(0,121)            
    for tm in dateiso:
<<<<<<< HEAD
        xdatencep.append(datetime.datetime.strftime(tm,"%b/%d")) 
=======
        xdatencep.append(datetime.datetime.strftime(tm,"%d/%b")) 
>>>>>>> e7f6294ce64f9ff8e82dba507be001724e7f2df1
###### fro daily
    det=datetime.timedelta(hours=24)            
    dateiso=[]            
    for dt in range(0,dnm):
        dateiso.append(datestart+dt*det)
    xdated=[]    
    xdatd=range(0,dnm)            
    for tm in dateiso:
<<<<<<< HEAD
        xdated.append(datetime.datetime.strftime(tm,"%b/%d"))         
=======
        xdated.append(datetime.datetime.strftime(tm,"%d/%b"))         
>>>>>>> e7f6294ce64f9ff8e82dba507be001724e7f2df1
########## reading files
    fpath=dirceres+ceresnm
    iskp=1
    onedim1=readAscii(fpath, iskp)
    nl=len(onedim1)
    timestep=[]
    ceresvrnm=['toa_comp_sw-up_pri_3h','toa_comp_sw-up_clr_3h', 'toa_comp_sw-up_naer_3h', # 0 1,2
           'toa_comp_sw-up_all_3h','toa_comp_lw-up_pri_3h', 'toa_comp_lw-up_clr_3h',  # 3 4 5
           'toa_comp_lw-up_naer_3h','toa_comp_lw-up_all_3h','toa_comp_wn-up_pri_3h',  #6 7 8 
           'toa_comp_wn-up_clr_3h', 'toa_comp_wn-up_naer_3h','toa_comp_wn-up_all_3h',  # 9 10 11
           'toa_comp_sw-down_all_3h','sfc_comp_sw-up_pri_3h', 'sfc_comp_sw-up_clr_3h',   # 12 13 14
           'sfc_comp_sw-up_naer_3h', 'sfc_comp_sw-up_all_3h', 'sfc_comp_sw-down_pri_3h',  # 15 16 17
           'sfc_comp_sw-down_clr_3h', 'sfc_comp_sw-down_naer_3h', 'sfc_comp_sw-down_all_3h', # 18 19 20
           'sfc_comp_lw-up_pri_3h','sfc_comp_lw-up_clr_3h', 'sfc_comp_lw-up_naer_3h',      # 21 22 23
           'sfc_comp_lw-up_all_3h', 'sfc_comp_lw-down_pri_3h', 'sfc_comp_lw-down_clr_3h',  # 24 25 26
           'sfc_comp_lw-down_naer_3h', 'sfc_comp_lw-down_all_3h', 'sfc_comp_wn-up_pri_3h',  # 27 28 29
           'sfc_comp_wn-up_clr_3h', 'sfc_comp_wn-up_naer_3h', 'sfc_comp_wn-up_all_3h',      # 30 31 32
           'sfc_comp_wn-down_pri_3h', 'sfc_comp_wn-down_clr_3h', 'sfc_comp_wn-down_naer_3h', # 33 34 35
           'sfc_comp_wn-down_all_3h' ]                                                          #36
    nvar=len(ceresvrnm)        
    ceresvar=np.ndarray(shape=(nvar,nt), dtype=float) 
    for i in range(0,nt):
        for j in range(0,nvar):
            k=i*nvar+j
            ceresvar[j,i]=onedim1[k]          
    del onedim1
    rsrflw=np.ndarray(shape=(nt), dtype=float) 
    rsrfsw=np.ndarray(shape=(nt), dtype=float)
    rtoalw=np.ndarray(shape=(nt), dtype=float) 
    rtoasw=np.ndarray(shape=(nt), dtype=float)
    for i in range(0,nt):
        if ceresvar[24,i]>0.0 and ceresvar[28,i]!= -999. :
            rsrflw[i]=ceresvar[28,i]/ceresvar[24,i]
        else:
            rsrflw[i]=NaN
        if ceresvar[20,i]>0.0 and ceresvar[16,i]!= -999. :
            rsrfsw[i]=ceresvar[16,i]/ceresvar[20,i]
        else:
            rsrfsw[i]=NaN
        if ceresvar[12,i]>0.0 and ceresvar[3,i]!= -999. :
            rtoasw[i]=ceresvar[3,i]/ceresvar[12,i]
        else:
            rtoasw[i]=NaN
    cerestoavar=['toa_sw_all_3h', 'toa_sw_clr_3h', 'toa_lw_all_3h', 'toa_lw_clr_3h',
             'toa_wn_all_3h', 'toa_wn_clr_3h' ]
    ntoa=len(cerestoavar)
    cerestoa=np.ndarray(shape=(ntoa,nt), dtype=float) 
    fpath=dirceres+ceresnmtoa
    iskp=1
    onedim1=readAscii(fpath, iskp)
    for i in range(0,nt):
        for j in range(0,ntoa):
            k=i*ntoa+j
            cerestoa[j,i]=onedim1[k]
            if onedim1[k] ==-999.0:
                cerestoa[j,i]=NaN    
    ceresalb=[]
    for i in range(0,nt):
        if ceresvar[9,i]!=0.0:
            ceresalb.append(ceresvar[18,i]/ceresvar[9,i])
        if ceresvar[9,i]==0.0:
            ceresalb.append(0.0) 
############################################################################
#                             
####################################      
    fpath=dirsim+simnm
    simvars=["aupst","adnst","auplt","adnlt",  # 0 1 2 3
         "aupss","adnss","aupls","adnls"]        # 4 5 6  7
    nvar=len(simvars)
    nst=241
    simvar=np.ndarray(shape=(nvar,nst), dtype=float)          
    linesplit=[]
    iskp=0
    onedim=readAscii(fpath, iskp)    
    ### radiaiton every three hours and unite is W/m**2
    for i in range(0,nst,1):
        for j in range(0,nvar):
            k=i*nvar+j
            simvar[j,i]= onedim[k]                     
    alb_sim=[]
    rssrflw=np.ndarray(shape=(nt), dtype=float) 
    rssrfsw=np.ndarray(shape=(nt), dtype=float)
    rstoalw=np.ndarray(shape=(nt), dtype=float) 
    rstoasw=np.ndarray(shape=(nt), dtype=float)
    for i in range(0,nt):
        if simvar[6,i]>0.0 :
            rssrflw[i]=simvar[7,i]/simvar[6,i]
        else:
            rssrflw[i]=NaN
        if simvar[5,i]>0.0 :
            rssrfsw[i]=simvar[4,i]/simvar[5,i]
        else:
            rssrfsw[i]-NaN
        if simvar[1,i]>0.0 :
            rstoasw[i]=simvar[0,i]/simvar[1,i]
        else:
            rstoasw[i]=NaN   
    nd=30
    avg1=0.
    avg2=0.
    for i in range(0,nt):
        avg1=avg1+ceresvar[3,i]
        avg2=avg2+simvar[0,i]
    avg1=avg1/nt
    avg2=avg2/nt 
    avgstr1="%.2f"%avg1
    avgstr2="%.2f"%avg2
    ndstr="%d"%nd  
    ndstr=ndstr+"days"
#####------- plotting------------------------------------------------
    #color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
    #color=color_cycle[0],lw=wd[0]
<<<<<<< HEAD
    line1=axs[ir,jc].plot(xdat[0:nt-1],rsrflw[0:nt-1],color=color_cycle[0],lw=wds[0],
        linestyle=linests[0],label=r'CERES $r_{sl}$')#[0:lcc-3])
    line2=axs[ir,jc].plot(xdat[0:nt-1],rssrflw[0:nt-1],color=color_cycle[1],
        lw=wds[1],linestyle=linests[1],label=r'CRM $r_{sl}$')
=======
    """
    line1=axs[ir,jc].plot(xdat[0:nt-1],rsrflw[0:nt-1],color=color_cycle[0],lw=wds[0],
        linestyle=linests[0],label=r'CERES $r_{sl}$')#[0:lcc-3])
    line2=axs[ir,jc].plot(xdat[0:nt-1],rssrflw[0:nt-1],color=color_cycle[1],
        lw=wds[1],linestyle=linests[1],label
        =r'CRM $r_{sl}$')
    """
>>>>>>> e7f6294ce64f9ff8e82dba507be001724e7f2df1
#    axs[ir,jc].plot(xdat[0:nt-1],rsrfsw[0:nt-1],color=color_cycle[2],lw=wds[2],
#        linestyle=linests[0],label='CERRES SRF Shortwave')#[0:lcc-3])
#    axs[ir,jc].plot(xdat[0:nt-1],rssrfsw[0:nt-1],color=color_cycle[3],lw=wds[3],
#        linestyle=linests[3],label='SIM SRF Shortwave')
    line3=axs[ir,jc].plot(xdat[0:nt-1],rtoasw[0:nt-1],color=color_cycle[4],lw=wds[4],
<<<<<<< HEAD
        linestyle=linests[4],label=r'CERES $r_{ts}$')#[0:lcc-3])
    line4=axs[ir,jc].plot(xdat[0:nt-1],rstoasw[0:nt-1],color=color_cycle[5],lw=wds[5],
        linestyle=linests[5],label=r'CRM $r_{ts}$')
#ax0.set_title('Upward shorwave radiaiton at TOA '+ r' $W m^{-2}$',fontsize=16)
    tilstr=astr+' '+area
#    xmajorLocator   = MultipleLocator(20) #将x主刻度标签设置为20的倍数  
#    xmajorFormatter = FormatStrFormatter('%5.1f') #设置x轴标签文本的格式  
#    xminorLocator   = MultipleLocator(5) #将x轴次刻度标签设置为5的倍数    
    ymajorLocator   = MultipleLocator(0.3) #将y轴主刻度标签设置为0.5的倍数  
#    ymajorFormatter = FormatStrFormatter('%1.1f') #设置y轴标签文本的格式 
    axs[ir,jc].yaxis.set_major_locator(ymajorLocator) 
    axs[ir,jc].set_title(tilstr, fontsize=18)
    axs[ir,jc].set_xticks(range(0,nt-1,24*4))
    xticklabels = [xdate[nn] for nn in range(0,nt-1,24*4)] 
=======
        linestyle=linests[4],label=r'CERES') # $r_{ts}$')#[0:lcc-3])
    line4=axs[ir,jc].plot(xdat[0:nt-1],rstoasw[0:nt-1],color=color_cycle[5],lw=wds[5],
        linestyle=linests[5],label=r'CRM') # $r_{ts}$')
#ax0.set_title('Upward shorwave radiaiton at TOA '+ r' $W m^{-2}$',fontsize=16)
    tilstr=astr+' '+area+' ('+"%4d"%iy+')'
#    xmajorLocator   = MultipleLocator(20) #将x主刻度标签设置为20的倍数  
#    xmajorFormatter = FormatStrFormatter('%5.1f') #设置x轴标签文本的格式  
#    xminorLocator   = MultipleLocator(5) #将x轴次刻度标签设置为5的倍数    
    ymajorLocator   = MultipleLocator(0.2) #将y轴主刻度标签设置为0.5的倍数  
#    ymajorFormatter = FormatStrFormatter('%1.1f') #设置y轴标签文本的格式
    axs[ir,jc].set_ylim(0.1,0.8)
    axs[ir,jc].yaxis.set_major_locator(ymajorLocator) 
    axs[ir,jc].set_title(tilstr, fontsize=18)
    axs[ir,jc].set_xticks(range(0,nt-1,24*3))
    xticklabels = [xdate[nn] for nn in range(0,nt-1,24*3)] 
>>>>>>> e7f6294ce64f9ff8e82dba507be001724e7f2df1
    axs[ir,jc].set_xticklabels(xticklabels, size=18)
    print iga
    if iga==2 :
        axs[ir,jc].legend(loc=(1.0,0.5),frameon=False)
    if iga in(1,2,4,5)  :
        setp(axs[ir,jc].get_yticklabels(), visible=False) #
    jc=jc+1
#
<<<<<<< HEAD
plt.show()                     
plt.savefig(dirout+'ALLCASE_radiation.png',dpi=300)          
=======
plt.show()
fig.subplots_adjust(left=0.1,bottom=0.1,right=1-0.1,top=1-0.1,hspace=0.4)                  
<<<<<<< HEAD
plt.savefig(dirout+'ALLCASE_radiation_albedo.png',dpi=300)          
>>>>>>> e7f6294ce64f9ff8e82dba507be001724e7f2df1
=======
plt.savefig(dirout+'ALLCASE_radiation_albedo_grey.png',dpi=300)          
>>>>>>> updated 20151022
plt.show()
plt.close()

###############################################################################