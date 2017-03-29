# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 13:20:18 2015

@author: chenjh
"""
import matplotlib.pyplot as plt
#plt.rc('text', usetex=True)
import matplotlib.cm as cm
import matplotlib as mpl
from pylab import *
import string
import numpy as np
import datetime
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['contour.negative_linestyle'] = 'dashed'
mpl.rcParams['ytick.labelsize'] = 18
mpl.rcParams['xtick.labelsize'] = 18
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
dirin='Z:/CRM/500m/Postdata/'
#dirin='D:/MyPaper/PhD04/Cases/alldata/'
dirpic='D:/MyPaper/PhD04/Pics/'
casenm=['PRDCTR_H','MLYRCTR_H', 'NPCCTR_H',
           'NECCTR_H','WTPCTR_H' , 'ETPCTR_H']   
astr=[r'$(a)$',r'$(b)$', r'$(c)$',r'$(d)$',r'$(e)$',r'$(f)$']
allyear=[2012,2010,2010,2012,2010,2010]
allmon=[4,6,8,7,7,6]
allday=[1,2,2,6,3,3]
nt=2880
nz=34
nga=len(casenm)
ydat_r=[50.000 ,   164.286,    307.143,    478.571  ,  678.571 ,
      907.143 ,  1164.286,   1450.000,   1764.286 ,  2107.143,   2478.572 ,
      2878.572,   3307.143,  3764.286,  4250.000,   4764.286,   5307.143, 
      5878.571,   6478.571,   7107.143,  7764.286,  8450.000,  9164.285,  
      9907.143,  10678.570,  11478.570,  12307.143,  13164.285,  14050.000,
      14964.285,  15907.143,  16878.572,  17878.572,  18907.145]
ydat=[]
for yd in ydat_r:
    ydat.append(yd*0.001)
rawdata=np.zeros(shape=(nga,2,nz,nt),dtype=float)
rawdata2=np.zeros(shape=(nga,2,nz,nt),dtype=float)
rawdata3=np.zeros(shape=(nga,2,nz,nt),dtype=float)
cloudpath=np.zeros(shape=(nga,5,2,nt),dtype=float)
ntt=nt/12
rawdata3h=np.zeros(shape=(nga,2,nz,ntt),dtype=float)
cloudpath3h=np.zeros(shape=(nga,5,2,ntt),dtype=float)
alldatestr=[]
maker=[]
xdat=range(0,nt) 
alldatestr3h=[]
xdat3h=range(0,ntt) 
nt1=nz+nz+5+5
print nt1
for iga in range(0,nga):
    fpath=dirin+casenm[iga]+'_ALLTYPES_CLOUDWATER_PATH.TXT'
    iskp=0
    iks=0
    ondedim=readAscii(fpath, iskp)
    itt=0
    for it in range(0,nt):
        iks=it*nt1 
        print iks
        for iz in range(0,nz):
            ikk=iks+iz               
            rawdata[iga,0,iz,it]=ondedim[ikk]               
            rawdata[iga,1,iz,it]=ondedim[ikk+nz]
        iks=iks+2*nz         
        for itk in range(0,5):
            ikk=iks+itk
            cloudpath[iga,itk,0,it]=ondedim[ikk]
            if ondedim[ikk]<0:
                cloudpath[iga,itk,0,it]=0.
        iks=iks+5
        for itk in range(0,5):
            ikk=iks+itk
            cloudpath[iga,itk,1,it]=ondedim[ikk]
            if ondedim[ikk]<0:
                cloudpath[iga,itk,1,it]=0.
        for itk in range(0,2):
            for iz in range(0,nz):
                if rawdata[iga,itk,iz,it]<0:
                    rawdata[iga,0,iz,it]=None
        print iks
    for itt in range(0,ntt):
        nnt1=itt*12-6
        nnt2=itt*12+6
        if itt==0:
            nnt1=0
            nnt2=7
        if itt==ntt-1:
            nnt1=nt-6
            nnt2=nt
        dnn=float(nnt2-nnt1)
        for iz in range(0,nz):
            for ik in range(0,2):
                rawdata3h[iga,ik,iz,itt]=0.
                for it1 in range(nnt1,nnt2):                    
                    rawdata3h[iga,ik,iz,itt]=rawdata3h[iga,ik,iz,itt]+rawdata[iga,ik,iz,it1]/dnn                    
        for itk in range(0,5):
            for ik in range(0,2):
                cloudpath3h[iga,itk,ik,itt]=0.
                for it1 in range(nnt1,nnt2):                   
                    cloudpath3h[iga,itk,ik,itt]=cloudpath3h[iga,itk,ik,itt]+cloudpath[iga,itk,ik,it1]/dnn              
    del ondedim
    fpath=dirin+casenm[iga]+'_TEST01_CHEN.TXT'
    iskp=0
    ondedim=readAscii(fpath, iskp)
    for it in range(0,nt):
        for iz in range(0,nz):
            iks=it*nz+it*nz+iz*2
            rawdata2[iga,0,iz,it]=ondedim[iks]
            rawdata2[iga,1,iz,it]=ondedim[iks+1]            
    fpath=dirin+casenm[iga]+'_TEST02_CHEN.TXT'
    iskp=0
    ondedim=readAscii(fpath, iskp)
    for it in range(0,nt):
        for iz in range(0,nz):
            iks=it*nz+it*nz++iz
            rawdata3[iga,0,iz,it]=ondedim[iks]
            iks=it*nz+it*nz++iz+nz
            rawdata3[iga,1,iz,it]=ondedim[iks]      
    # end reading and generate the data strings
    year,mon,day=allyear[iga],allmon[iga],allday[iga]
    datestart=datetime.datetime(year,mon,day,0,0,0)
    det=datetime.timedelta(minutes=15)            
    xdate=[]             
    for dt in range(0,nt):
        tempdate=datestart+dt*det
        xdate.append(datetime.datetime.strftime(tempdate,"%d/%b"))
    alldatestr.append(xdate)
    del xdate
    det=datetime.timedelta(hours=3)            
    xdate=[]             
    for dt in range(0,ntt):
        tempdate=datestart+dt*det
        xdate.append(datetime.datetime.strftime(tempdate,"%d/%b"))
    alldatestr3h.append(xdate)
    if casenm[iga][0:3]=='MLY' :
        maker.append(astr[iga]+casenm[iga][0:4])
    else:
        maker.append(astr[iga]+casenm[iga][0:3])
#############################################################################

#############################################################################
fig,ax=plt.subplots(nrows=3,ncols=2,figsize=(18,12))
#fig,axs=plt.subplots(nrows=2,ncols=3,figsize=(12,12))
color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
wd=[2,2,2,2,2]
jc=0
jr=0
ij=1
mycolor=cm.YlOrRd
cdict = {'red': ((0., 1, 1),
                 (0.11, 0, 0),
                 (0.66, 1, 1),
                 (0.89, 1, 1),
                 (1, 0.5, 0.5)),
         'green': ((0., 1, 1),
                   (0.11, 0, 0),
                   (0.375, 1, 1),
                   (0.64, 1, 1),
                   (0.91, 0, 0),
                   (1, 0, 0)),
         'blue': ((0., 1, 1),
                  (0.11, 1, 1),
                  (0.34, 1, 1),
                  (0.65, 0, 0),
                  (1, 0, 0))}

my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
for iga in range(0,nga):
    if jc==2:
        jc=0
        jr=jr+1
    mker=maker[iga]
    xdate=alldatestr[iga]
    levs1=[0.2,0.4,0.6,0.8,1.0,1.2,1.5,2]
    levs2=[0.2,0.4,0.6,0.8,1.0,1.2,1.5,2]
    levs1=[0.0,0.01,0.03,0.04,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.3]
    levs2=[0.0,0.01,0.03,0.04,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.3]
    font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'bold',
        'size'   : 16,
        }     
    plt.subplot(3,2,ij)
#    zdat=rawdata[iga,0,:,:]+rawdata[iga,1,:,:]
    #zdat=smtco   #-obstmp
    ##zdat[0,:]=0.0   ## the first level is below surface ground
#    ax[jr,jc]=plt.contour(xdat,ydat,zdat,colors='k',
#        linewidths=1.5,levels=levs1)                           
#    plt.title(titlename[0],fontsize=16)                          
#    plt.axis([0, nt+1, 0, 16])
#    plt.clabel(ax[jr,jc],inline=1,fmt='%2.2f',fontsize=12)   
#    zdat=rawdata[iga,1,:,:]+rawdata[iga,0,:,:]
    zdat=rawdata[iga,1,:,:]
    ax[jr,jc]=plt.contourf(xdat,ydat,zdat,cmap=my_cmap ,extend='both', #colors='grey',
        linewidths=1.5,levels=levs2)  
#    ax[jr,jc]=plt.contourf(xdat,ydat,zdat, cmap=cm.jet,extend='both',
#        levels=levs2)#,linestyles=linetype2)  #linewidths=1.5,
    plt.axis([0, nt+1, 0, 16])
#    plt.clabel(ax[jr,jc],inline=1,fmt='%2.1f',fontsize=12)
    axx=fig.add_subplot(3,2,ij)
    ymajorLocator   = MultipleLocator(4)
    axx.yaxis.set_major_locator(ymajorLocator) 
    text1=mker  #r"($a$)"
    axx.text(1.5,16.5,text1,fontsize=18)                        
    axx.set_xticks(range(0,nt,96*4))
    xticklabels = [xdate[nn] for nn in range(0,nt,96*4)] 
    axx.set_xticklabels(xticklabels, size=16)
    plt.ylabel(r'Height ($km$)', fontdict=font)
    jc=jc+1
    ij=ij+1                
plt.show()
fig.subplots_adjust(left=0.1,bottom=0.1,right=1-0.1,top=1-0.1,hspace=0.4)
cax = fig.add_axes([0.2, 0.03, 0.6, 0.03])
fig.colorbar(ax[0,0], cax,extend='both',
              spacing='uniform', orientation='horizontal')
plt.savefig(dirpic+"ALLCASE_IceCloudWater_500m.png",dpi=300)          
plt.show()
plt.close()
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
    mker=maker[iga]
    xdate=alldatestr[iga]
    levs1=[0.2,0.4,0.6,0.8,1.0,1.2,1.5,2]
    levs2=[0.2,0.4,0.6,0.8,1.0,1.2,1.5,2]
    levs1=[0.0,0.01,0.03,0.04,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.3]
    levs2=[0.0,0.01,0.03,0.04,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.3]
    font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 16,
        }     
    plt.subplot(3,2,ij)
    zdat=rawdata[iga,0,:,:]+rawdata[iga,1,:,:]
    #zdat=smtco   #-obstmp
    ##zdat[0,:]=0.0   ## the first level is below surface ground
    ax[jr,jc]=plt.contourf(xdat,ydat,zdat,cmap=my_cmap ,extend='both',
        linewidths=1.5,levels=levs1)                           
#    plt.title(titlename[0],fontsize=16)                          
    plt.axis([0, nt+1, 0, 16])
#    plt.clabel(ax[jr,jc],inline=1,fmt='%2.2f',fontsize=12)   
#    zdat=rawdata[iga,1,:,:]+rawdata[iga,0,:,:]
#    zdat=rawdata[iga,1,:,:]
#    ax[jr,jc]=plt.contourf(xdat,ydat,zdat,cmap=cm.jet , #colors='grey',
#        linewidths=1.5,levels=levs2)  
#    ax[jr,jc]=plt.contourf(xdat,ydat,zdat, cmap=cm.jet,extend='both',
#        levels=levs2)#,linestyles=linetype2)  #linewidths=1.5,
#    plt.axis([0, nt+1, 0, 16])
#    plt.clabel(ax[jr,jc],inline=1,fmt='%2.1f',fontsize=12)
    axx=fig.add_subplot(3,2,ij)
    ymajorLocator   = MultipleLocator(4)
    axx.yaxis.set_major_locator(ymajorLocator) 
    text1=mker  #r"($a$)"
    axx.text(1.5,16.5,text1,fontsize=18)                        
    axx.set_xticks(range(0,nt,96*4))
    xticklabels = [xdate[nn] for nn in range(0,nt,96*4)] 
    axx.set_xticklabels(xticklabels, size=16)
    plt.ylabel(r'Height ($km$)', fontdict=font)
    jc=jc+1
    ij=ij+1                
plt.show()
fig.subplots_adjust(left=0.1,bottom=0.1,right=1-0.1,top=1-0.1,hspace=0.4)
cax = fig.add_axes([0.2, 0.03, 0.6, 0.03])
fig.colorbar(ax[0,0], cax,extend='both',
              spacing='uniform', orientation='horizontal')
plt.savefig(dirpic+"ALLCASE_TotalCloudWater_500m.png",dpi=300)          
plt.show()
plt.close()        
#############################################################################
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
    mker=maker[iga]
    xdate=alldatestr[iga]
    font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 16,
        } 
    width=0.
#    zdat=cloudpath[iga,3,0,:]+cloudpath[iga,3,1,:]
#    zdat=[]
    plt.subplot(3,2,ij)
    zdat=cloudpath[iga,1,0,:]+cloudpath[iga,1,1,:]
#    zdat.append(cloudpath[iga,1,0,:])
#    zdat.append(cloudpath[iga,1,1,:])
    #zdat=smtco   #-obstmp
    ##zdat[0,:]=0.0   ## the first level is below surface ground
#    ax[jr,jc].plot(xdat,zdat,color='k',linewidth=3)
    p1=plt.bar(xdat,zdat, width, color='b',edgecolor='b')
#    ax[jr,jc].hist( zdat,len(xdat), histtype='step', stacked=True,facecolor='g', fill=True)                           
#    plt.title(titlename[0],fontsize=16)                          
    plt.xlim([0, nt+1])
#    plt.ylim([0, 601])
    zdat=cloudpath[iga,1,1,:] 
    p2=plt.bar(xdat,zdat, width, color='cyan',edgecolor='cyan')
#    ax[jr,jc].hist( zdat,xdat, histtype='step', stacked=True,facecolor='g', fill=True) 
#    ax[jr,jc].plot(xdat,zdat,color='grey',linewidth=3) 
#    plt.xlim([0, nt+1])
    axx=fig.add_subplot(3,2,ij)   
    if ij<5: 
        ymax=2801
        ymajorLocator   = MultipleLocator(450)
        axx.yaxis.set_major_locator(ymajorLocator) 
    else:
        ymax=401
        ymajorLocator   = MultipleLocator(150)
        axx.yaxis.set_major_locator(ymajorLocator) 
    plt.ylim([0, ymax])    
    text1=mker  #r"($a$)"
    axx.text(1.5,ymax+5,text1,fontsize=18)                        
    axx.set_xticks(range(0,nt+1,96*4))
    xticklabels = [xdate[nn] for nn in range(0,nt+1,96*4)] 
    axx.set_xticklabels(xticklabels, size=16)
    plt.grid(ls=':',lw=0.5)
#    plt.ylabel(r'Height ($km$)', fontdict=font)
    jc=jc+1
    ij=ij+1                
plt.show()
fig.subplots_adjust(left=0.1,bottom=0.1,right=1-0.1,top=1-0.1,hspace=0.4)
#cax = fig.add_axes([0.2, 0.03, 0.6, 0.03])
#fig.colorbar(ax[0,0], cax,extend='both',
#              spacing='uniform', orientation='horizontal')
plt.savefig(dirpic+"ALLCASE_CloudWaterPath_500m.png",dpi=300)          
plt.show()
plt.close()
###############################################################################
###   for 3 hours mean 
###############################################################################
fig,ax=plt.subplots(nrows=3,ncols=2,figsize=(18,12))
#fig,axs=plt.subplots(nrows=2,ncols=3,figsize=(12,12))
color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
wd=[2,2,2,2,2]
jc=0
jr=0
ij=1
mycolor=cm.YlOrRd
cdict = {'red': ((0., 1, 1),
                 (0.11, 0, 0),
                 (0.66, 1, 1),
                 (0.89, 1, 1),
                 (1, 0.5, 0.5)),
         'green': ((0., 1, 1),
                   (0.11, 0, 0),
                   (0.375, 1, 1),
                   (0.64, 1, 1),
                   (0.91, 0, 0),
                   (1, 0, 0)),
         'blue': ((0., 1, 1),
                  (0.11, 1, 1),
                  (0.34, 1, 1),
                  (0.65, 0, 0),
                  (1, 0, 0))}

my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
my_cmap=cm.binary
for iga in range(0,nga):
    if jc==2:
        jc=0
        jr=jr+1
    mker=maker[iga]
    xdate=alldatestr3h[iga]
    levs1=[0.2,0.4,0.6,0.8,1.0,1.2,1.5,2]
    levs2=[0.2,0.4,0.6,0.8,1.0,1.2,1.5,2]
    levs1=[0.0,0.01,0.03,0.04,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.3]
    levs2=[0.0,0.01,0.03,0.04,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.3]
    font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 16,
        }     
    plt.subplot(3,2,ij)
#    zdat=rawdata[iga,0,:,:]+rawdata[iga,1,:,:]
    #zdat=smtco   #-obstmp
    ##zdat[0,:]=0.0   ## the first level is below surface ground
#    ax[jr,jc]=plt.contour(xdat,ydat,zdat,colors='k',
#        linewidths=1.5,levels=levs1)                           
#    plt.title(titlename[0],fontsize=16)                          
#    plt.axis([0, nt+1, 0, 16])
#    plt.clabel(ax[jr,jc],inline=1,fmt='%2.2f',fontsize=12)   
#    zdat=rawdata[iga,1,:,:]+rawdata[iga,0,:,:]
    zdat=rawdata3h[iga,1,:,:]
    ax[jr,jc]=plt.contourf(xdat3h,ydat,zdat,cmap=my_cmap ,extend='both', #colors='grey',
        linewidths=1.5,levels=levs2)  
#    ax[jr,jc]=plt.contourf(xdat,ydat,zdat, cmap=cm.jet,extend='both',
#        levels=levs2)#,linestyles=linetype2)  #linewidths=1.5,
    plt.axis([0, ntt+1, 0, 16])
#    plt.clabel(ax[jr,jc],inline=1,fmt='%2.1f',fontsize=12)
    axx=fig.add_subplot(3,2,ij)
    ymajorLocator   = MultipleLocator(4)
    axx.yaxis.set_major_locator(ymajorLocator) 
    text1=mker  #r"($a$)"
    axx.text(1.5,16.5,text1,fontsize=18)                        
    axx.set_xticks(range(0,ntt,8*4))
    xticklabels = [xdate[nn] for nn in range(0,ntt,8*4)] 
    axx.set_xticklabels(xticklabels, size=16)
    plt.ylabel(r'Height ($km$)', fontdict=font)
    jc=jc+1
    ij=ij+1                
plt.show()
fig.subplots_adjust(left=0.1,bottom=0.1,right=1-0.1,top=1-0.1,hspace=0.4)
cax = fig.add_axes([0.2, 0.03, 0.6, 0.03])
fig.colorbar(ax[0,0], cax,extend='both',
              spacing='uniform', orientation='horizontal')
plt.savefig(dirpic+"ALLCASE_IceCloudWater_3h_500m.png",dpi=300)          
plt.show()
plt.close()
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
    mker=maker[iga]
    xdate=alldatestr3h[iga]
    levs1=[0.2,0.4,0.6,0.8,1.0,1.2,1.5,2]
    levs2=[0.2,0.4,0.6,0.8,1.0,1.2,1.5,2]
    levs1=[0.0,0.01,0.03,0.04,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.3]
    levs2=[0.0,0.01,0.03,0.04,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.3]
    font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 16,
        }     
    plt.subplot(3,2,ij)
    zdat=rawdata3h[iga,0,:,:]+rawdata3h[iga,1,:,:]
    #zdat=smtco   #-obstmp
    ##zdat[0,:]=0.0   ## the first level is below surface ground
    ax[jr,jc]=plt.contourf(xdat3h,ydat,zdat,cmap=my_cmap ,extend='both',
        linewidths=1.5,levels=levs1)                           
#    plt.title(titlename[0],fontsize=16)                          
    plt.axis([0, ntt+1, 0, 16])
#    plt.clabel(ax[jr,jc],inline=1,fmt='%2.2f',fontsize=12)   
#    zdat=rawdata[iga,1,:,:]+rawdata[iga,0,:,:]
#    zdat=rawdata[iga,1,:,:]
#    ax[jr,jc]=plt.contourf(xdat,ydat,zdat,cmap=cm.jet , #colors='grey',
#        linewidths=1.5,levels=levs2)  
#    ax[jr,jc]=plt.contourf(xdat,ydat,zdat, cmap=cm.jet,extend='both',
#        levels=levs2)#,linestyles=linetype2)  #linewidths=1.5,
#    plt.axis([0, nt+1, 0, 16])
#    plt.clabel(ax[jr,jc],inline=1,fmt='%2.1f',fontsize=12)
    axx=fig.add_subplot(3,2,ij)
    ymajorLocator   = MultipleLocator(4)
    axx.yaxis.set_major_locator(ymajorLocator) 
    text1=mker  #r"($a$)"
    axx.text(1.5,16.5,text1,fontsize=18)                        
    axx.set_xticks(range(0,ntt,8*4))
    xticklabels = [xdate[nn] for nn in range(0,ntt,8*4)] 
    axx.set_xticklabels(xticklabels, size=16)
    plt.ylabel(r'Height ($km$)', fontdict=font)
    jc=jc+1
    ij=ij+1                
plt.show()
fig.subplots_adjust(left=0.1,bottom=0.1,right=1-0.1,top=1-0.1,hspace=0.4)
cax = fig.add_axes([0.2, 0.03, 0.6, 0.03])
fig.colorbar(ax[0,0], cax,extend='both',
              spacing='uniform', orientation='horizontal')
plt.savefig(dirpic+"ALLCASE_TotalCloudWater_3h_gray_500m.png",dpi=300)          
plt.show()
plt.close()        
#############################################################################
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
    mker=maker[iga]
    xdate=alldatestr3h[iga]
    font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 16,
        } 
    width=1
#    zdat=cloudpath[iga,3,0,:]+cloudpath[iga,3,1,:]
#    zdat=[]
    plt.subplot(3,2,ij)
    zdat=cloudpath3h[iga,1,0,:]+cloudpath3h[iga,1,1,:]
#    zdat.append(cloudpath[iga,1,0,:])
#    zdat.append(cloudpath[iga,1,1,:])
    #zdat=smtco   #-obstmp
    ##zdat[0,:]=0.0   ## the first level is below surface ground
#    ax[jr,jc].plot(xdat,zdat,color='k',linewidth=3)
    p1=plt.bar(xdat3h,zdat/1000., width, color='k',edgecolor='k')
#    ax[jr,jc].hist( zdat,len(xdat), histtype='step', stacked=True,facecolor='g', fill=True)                           
#    plt.title(titlename[0],fontsize=16)                          
    plt.xlim([0, ntt+1])
#    plt.ylim([0, 601])
    zdat=cloudpath3h[iga,1,1,:] 
    p2=plt.bar(xdat3h,zdat/1000., width, color='gray',edgecolor='gray')
#    ax[jr,jc].hist( zdat,xdat, histtype='step', stacked=True,facecolor='g', fill=True) 
#    ax[jr,jc].plot(xdat,zdat,color='grey',linewidth=3) 
#    plt.xlim([0, nt+1])
    axx=fig.add_subplot(3,2,ij)   
    if ij<5: 
        ymax=1.9
        ymajorLocator   = MultipleLocator(0.4)
        axx.yaxis.set_major_locator(ymajorLocator) 
    else:
        ymax=0.4
        ymajorLocator   = MultipleLocator(0.2)
        axx.yaxis.set_major_locator(ymajorLocator) 
    plt.ylim([0, ymax])    
    text1=mker  #r"($a$)"
    axx.text(1.5,ymax+5/20,text1,fontsize=18)                        
    axx.set_xticks(range(0,ntt+1,8*4))
    xticklabels = [xdate[nn] for nn in range(0,ntt+1,8*4)] 
    axx.set_xticklabels(xticklabels, size=16)
    plt.grid(ls=':',lw=0.5)
#    plt.ylabel(r'Height ($km$)', fontdict=font)
    jc=jc+1
    ij=ij+1                
plt.show()
fig.subplots_adjust(left=0.1,bottom=0.1,right=1-0.1,top=1-0.1,hspace=0.4)
#cax = fig.add_axes([0.2, 0.03, 0.6, 0.03])
#fig.colorbar(ax[0,0], cax,extend='both',
#              spacing='uniform', orientation='horizontal')
plt.savefig(dirpic+"ALLCASE_CloudWaterPath_3h_gray_500m.png",dpi=300)          
plt.show()
plt.close()