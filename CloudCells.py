#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 09:56:43 2015

@author: jhchen
"""
import string
import numpy as np
import matplotlib.pyplot as plt

def readAscii(fpath,iskp,nrl):
    #iskp  the total line skipped of the file
    # fpath   the full path of the file
    # usage: onedim=readAscii(fpaht,iskp)
    onedim=[]
    linesplit=[]
    f=open(fpath)
    print iskp,nrl
    ff=f.readlines()[iskp:nrl]  ## first line in obs file is legend 
    for line in ff:
        line=string.lstrip(line)
        linesplit.append(line[:-1].split(' '))
    for lnstrs in linesplit:
        for strs in lnstrs:
            if strs!='':
                onedim.append(string.atof(strs))
    del linesplit,ff
    f.close()
    print len(onedim)
    return onedim
def infcld(c,km):
    na=0
    kb=np.ndarray(shape=(99), dtype=int)
    ke=np.ndarray(shape=(99), dtype=int)
    cm=np.ndarray(shape=(99), dtype=float)
    if c[0]>0. :
        k1=0
        k2=0
        aa=c[0]
    else:
        aa=0.0
    for k in range(1,km):
        if c[k-1] <=0. and c[k] > 0. :
            k1=k
            k2=k
            aa=max(aa,c[k])
        elif c[k-1] >0. and c[k]>0. :
            k2=k
            aa=max(aa,c[k])
        elif c[k-1] > 0.  and c[k]<= 0.0 :           
            kb[na]=k1
            ke[na]=k2
            cm[na]=aa
            aa=0.
            na=na+1
    if c[km-1] > 0.0 :
        kb[na]=k1
        ke[na]=k2
        cm[na]=aa
        na=na+1
    return kb,ke,cm,na 
def infcldx(c,km):
    na=0
    kb=np.ndarray(shape=(299), dtype=int)
    ke=np.ndarray(shape=(299), dtype=int)
    if c[0]>0. :
        k1=0
        k2=0
        aa=c[0]
    else:
        aa=0.0
    for k in range(1,km):
        if c[k-1] <=0. and c[k] > 0. :
            k1=k
            k2=k
            aa=max(aa,c[k])
        elif c[k-1] >0. and c[k]>0. :
            k2=k
            aa=max(aa,c[k])
        elif c[k-1] > 0.  and c[k]<= 0.0 :           
            kb[na]=k1
            ke[na]=k2
#            cm(na)=aa
            aa=0.
            na=na+1
    if c[km-1] > 0.0 :
        kb[na]=k1
        ke[na]=k2
#        cm(na)=aa
        na=na+1
    return kb,ke,na      
def deepcc(qw,im,km):
    c=np.ndarray(shape=(im,km), dtype=float)
    kb=np.ndarray(shape=(im,99), dtype=int)
    ke=np.ndarray(shape=(im,99), dtype=int)
    cm=np.ndarray(shape=(im,99), dtype=float)
    na=np.ndarray(shape=(im), dtype=int)
    for i in range(0,im):
        for k in range(0,km):
            if qw[i,k] >= 1.0e-3:
                c[i,k]=1. # cloud
            else:
                c[i,k]=0.
        kb[i,:],ke[i,:],cm[i,:],na[i]=infcld(c[i,:],km)
    return kb,ke,c,cm,na       
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
casenm="ETPCTR_EC"
if casenm == "MLY" :
    dirin ="D:/MyPaper/PhD04/Cases/"+casenm[0:4]+"/"
    regionname=casenm[0:4]
else:
    dirin ="D:/MyPaper/PhD04/Cases/"+casenm[0:3]+"/"
    regionname=casenm[0:3]
dirout="D:/MyPaper/PhD04/Cases/Postdata/"
dirpic="D:/MyPaper/PhD04/Pics/"
fold="20100604_0704/Simulated/"
fold="CTREC20120520/Simulation/"
nd=30 # days
dt=15 # output timestep in min
nx=200
nz=34
nt=nd*24*60/dt
namenp=["Total Coloud Water Content","Cloud Water Path",""]
z=[              0.0500000, 0.1643000, 0.3071000, 0.4786000
    , 0.6786000, 0.9071000, 1.1640000, 1.4500000, 1.7640001
    , 2.1070001, 2.4790001, 2.8789999, 3.3069999, 3.7639999
    , 4.2500000, 4.7639999, 5.3070002, 5.8790002, 6.4790001
    , 7.1069999, 7.7639999, 8.4499998, 9.1639996, 9.9069996
    ,10.6800003,11.4799995,12.3100004,13.1599998,14.0500002
    ,14.9600000,15.9099998,16.8799992,17.8799992,18.9099998]
qcc,qce,cwp=1.0e-2,1.0e-4,0.2   #! low limit for convection cloud (g/kg)
                                # low limit for cloud ensemble (g/kg)
                                #cloud water path threshold (g/m2)
                                #--if cwp>0 then qc0=cwp/(d dz) else =qcc
                                # --qc0 is used for cloud identification
ztd,zbd=10.0,7.0 # deep convection top, base heights (km)
ztp,zbp=14.0,2.0  #PBL convection top, base heights (km)
ksc=4             #shallow convection layer thickness in grid
zbmin=0.5  #minimum base height to count cloud (km)
no=3  # over lap cloud layer
npp=3
ne=1
spv=1.0e20 # default value
fout=dirout+regionname+"_deepcon_basetop.txt"
fdc=open(fout,'w')
fout=dirout+regionname+"_cloudcells.txt"
fcells=open(fout,'w')
fout=dirout+"_Plotting_1d.txt"
fo1=open(fout,'w')
fout=dirout+"_Plotting_2d.txt"
fo2=open(fout,'w')
#
qz=np.ndarray(shape=(nz), dtype=float)
wz=np.ndarray(shape=(nz), dtype=float)
qd=np.ndarray(shape=(nz,nz,npp), dtype=float)
td=np.ndarray(shape=(nz,nz), dtype=float)
wd=np.ndarray(shape=(nz,nz), dtype=float)
#
qb=np.ndarray(shape=(nz), dtype=float)
wb=np.ndarray(shape=(nz), dtype=float)
qp=np.ndarray(shape=(nz,nz,npp), dtype=float)
tp=np.ndarray(shape=(nz,nz), dtype=float)
wp=np.ndarray(shape=(nz,nz), dtype=float)
#
qs=np.ndarray(shape=(nz,npp), dtype=float)
ts=np.ndarray(shape=(nz), dtype=float)
ws=np.ndarray(shape=(nz), dtype=float)
#
fa=np.ndarray(shape=(no+2,ne), dtype=float)
fb=np.ndarray(shape=(nz,nz,no,ne), dtype=float)
sa=np.ndarray(shape=(ne), dtype=float)
cels=[0.00,1.00]
#
ft=np.ndarray(shape=(nz,nz), dtype=float) #(km,km)
tr=np.ndarray(shape=(nz), dtype=float) #(km)
#
#   
fpath=dirin+fold+casenm+"_qlqit.txt"
iskp=0
sg=0.
for it in range(0,nt):
    iskp=it*76*nz   # (200+200+200+1)/8
    nrlines=76*nz+iskp
    onedim=readAscii(fpath,iskp,nrlines)
    ql=np.ndarray(shape=(nx,nz), dtype=float)  # liquid
    qi=np.ndarray(shape=(nx,nz), dtype=float)  # ice
    tc=np.ndarray(shape=(nx,nz), dtype=float)  # temp
    den=np.ndarray(shape=(nz), dtype=float)  # density
    wk=np.ndarray(shape=(nx,nz), dtype=float)  # total
    lcb=np.ndarray(shape=(nz), dtype=int)  # density
    lce=np.ndarray(shape=(nz), dtype=int)  # density
    for iz in range(0,nz):
        k=iz*601
        den[iz]=onedim[k]
        for ix in range(0,nx):
            k=iz*601+ix+1
            ql[ix,iz]=onedim[k]
            k=iz*601+ix+200+1
            qi[ix,iz]=onedim[k]
            k=iz*601+ix+400+1
            tc[ix,iz]=onedim[k]
            wk[ix,iz]=ql[ix,iz]+qi[ix,iz]
    del onedim
    zd=np.ndarray(shape=(nz), dtype=float)
    qc0=np.ndarray(shape=(nz), dtype=float)
    zd[0]=0.5*(z[1]-z[0])*den[0]
    for k in range(1,nz-1):
        zd[k]=0.5*(z[k+1]-z[k-1])*den[k]
    zd[nz-1]=0.5*(z[nz-1]-z[nz-2])*den[nz-1]
    if cwp > 0.0 :
        for k in range(0,nz):
            qc0[k]=1e-3*cwp/zd[k]
    elif qcc > 0.0 :
        for k in range(0,nz):
            qc0[k]=qcc
    else:
        print "xxxx wrong cloud id threshold"
        exit()
    szd=0.0
    for iz in range(0,nz):
        szd=szd+zd[iz]
    kdb,kde,c,cm,na =deepcc(wk,nx,nz)
    for ix in range(0,nx):
        for i in range(0,na[ix]):
            if kdb[ix,i] <=5 and kde[ix,i] > 24 :
                for iz in range(0,nz):
                    itme="%f "%(c[ix,iz]+0.1)
                    fdc.write(itme)
                itme="%d "%kdb[ix,i]
                fdc.write(itme)
                itme="%d "%kde[ix,i]
                fdc.write(itme)
        fdc.write('\n')  # Please note the file structure when the file is read.
#
#
    frc=1./(nx*nz)
    frg=1./nx
    frr=1./(nx*nt)
    for k in range(0,nz):
        for ix in range(0,nx):
            tr[k]=tr[k]+frr*tc[ix,k]
    sumf=0.0
    ezd=0.0
    qmx=0.0
    for i in range(0,nx):   # statistics 200 grids from 0 to nx
        for k in range(0,nz):
            sumf=sumf+frc*wk[i,k]
            ezd=ezd+frg*wk[i,k]*zd[k]
            qmx=max(qmx,wk[i,k])
    ezd=ezd/szd
    #--------get pdf for cloud ensemble
    if sumf > qce :
        sg=sg+1.
        wa=np.ndarray(shape=(no+2), dtype=float)
        ww=np.ndarray(shape=(nz,nz,no), dtype=float)
        kcb=np.ndarray(shape=(nx,99), dtype=int)
        kce=np.ndarray(shape=(nx,99), dtype=int)
        cm=np.ndarray(shape=(nx,99), dtype=float)
        ns=np.ndarray(shape=(nx), dtype=int)
        cc=np.ndarray(shape=(nx,nz), dtype=float)
        ce=0.
        for ix in range(0,nx):
            for iz in range(0,nz):
                cc[ix,iz]=wk[ix,iz]
                if cc[ix,iz] < qc0[iz] :
                    cc[ix,iz]=0.
            kcb[ix,:],kce[ix,:],cm[ix,:],ns[ix]=infcld(cc[ix,:],nz)       
            for ii in range(0,ns[ix]):
                kb=kcb[ix,ii]
                ke=kce[ix,ii]
                ft[kb,ke]=ft[kb,ke]+frg  # for ft[kb,ke],the value of ft is cloud cover, 
                                         # z[kb] is base,z[ke]is cloud top
#c...............vertical profiles for deep convections                
            for ii in range(0,ns[ix]):
                kb=kcb[ix,ii]
                ke=kce[ix,ii]
                if z[ke]>ztd and z[kb] <= zbd :
                    czd=0.
                    qzd=0.
                    for k in range(kb,ke):
                        czd=czd+zd[k]
                        qzd=qzd+zd[k]*cc[ix,k]
                    qzd=qzd/czd
                    qz[kb]=qz[kb]+qzd
                    wz[kb]=wz[kb]+1.
                    for k in range(kb,ke):
                        qd[k,kb,0]=qd[k,kb,0]+cc[iz,k]  # total cloud water
                        qd[k,kb,1]=qd[k,kb,1]+cc[iz,k]/qzd
                        qd[k,kb,2]=qd[k,kb,2]+cc[iz,k]/qzd*den[k]
                        td[k,kb]=td[k,kb]+tc[ix,k]
                        wd[k,kb]=wd[k,kb]+1.
#c...............vertical profiles for PBL convections                
            for ii in range(0,ns[ix]):
                kb=kcb[ix,ii]
                ke=kce[ix,ii]
                if z[ke]>ztp and z[kb] <= zbp :
                    czd=0.
                    qzd=0.
                    for k in range(kb,ke):
                        czd=czd+zd[k]
                        qzd=qzd+zd[k]*cc[ix,k]
                    qzd=qzd/czd
                    qb[kb]=qb[kb]+qzd
                    wb[kb]=wb[kb]+1.
                    for k in range(kb,ke):
                        qp[k,ke,0]=qp[k,ke,0]+cc[iz,k]
                        qp[k,ke,1]=qp[k,ke,1]+cc[iz,k]/qzd
                        qp[k,ke,2]=qp[k,ke,2]+cc[iz,k]/qzd*den[k]
                        tp[k,ke]=tp[k,ke]+tc(ix,k)
                        wp[k,ke]=wp[k,ke]+1.
#c...............vertical profiles for shallow convections                
            if ns[ix]==1:
                ii=0
                kb=kcb[ix,ii]
                ke=kce[ix,ii]
                if ke-kb < ksc and ke > kb:
                    czd=0.
                    qzd=0.
                    tzd=0.
                    for k in range(kb,ke):
                        czd=czd+zd[k]
                        qzd=qzd+zd[k]*cc[ix,k]
                        tzd=tzd+zd[k]*tc[ix,k]
                    qzd=qzd/czd
                    tzd=tzd/czd
                    qs[kb,0]=qs[kb,0]+qzd
                    qs[kb,1]=qs[kb,1]+ezd
                    qs[kb,2]=qs[kb,2]+qzd/ezd*den[k]
                    ts[kb]=ts[kb]+tzd
                    ws[kb]=ws[kb]+1.
#...............for overlap cloud cells
#               ----ignore thin cells
            lc=0
            lccc=-1
            for ii in range(0,ns[ix]):
                if (kce[ix,ii]-kcb[ix,ii]) >= 2 :
                    lce[lc]=kce[ix,ii]
                    lcb[lc]=kcb[ix,ii]
                    lc=lc+1
                    lccc=9
            #----ignore small gaps
            if lccc >0 : # more than one cloud layer
                for ii in range(1,lc):
                    if lcb[ii]-lce[ii-1]<2 :
                        lce[ii-1]=lce[ii]
                        lc=lc-1
                        for ic in range(ii,lc):
                            lcb[ic]=lcb[ic+1]
                            lce[ic]=lce[ic+1]
            #----get frequency according to overlap cells
            if lccc >0 :
                ce=ce+frg
                io=min(lc,no+1)
                wa[io]=wa[io]+frg
            if lc == 1 :
                kb=lcb[0] # bas of the cell
                ke=lce[0] # top of the cell                
            elif lc == 2 :
                kb=lce[0]  # top of low cell
                ke=lcb[1]  # bas of upp cell                        
            elif lc == 3 :
                kb=lcb[1]  # top of mid cell
                ke=lce[1]  # bas of mid cell  
            if lc>=1 and lc<=3 :
                lcc=lc-1
                ww[kb,ke,lcc]=ww[kb,ke,lcc]+frg
#        if ce ==0. :
#            ce=0.0001
        wa[0]=1.0-ce # clear sky
        ich=-1
        for ii in range(0,ne):
            if cels[ii] < ce and ce <=cels[ii+1] :
                ic=ii
                ich=1
        if ich>0 :
            sa[ic]=sa[ic]+1.
            for io in range(0,no+2):
                fa[io,ic]=fa[io,ic]+wa[io]
            for io in range(0,no):
                for ke in range(0,nz):
                    for kb in range(0,nz):
                        fb[kb,ke,io,ic]=fb[kb,ke,io,ic]+ww[kb,ke,io]
# end loop of it
fdc.close()
#
if sg > 0. :
    for kb in range(0,nz):
        for kb in range(0,nz):
            ft[kb,ke]=  ft[kb,ke]*100./sg    #  # for ft[kb,ke],the value of ft is cloud cover, 
                                         # z[kb] is base,z[ke]is cloud top                         
# deep convection
for ip in range(0,npp):
    for k in range(0,nz):
        for kb in range(0,nz):
            if wd[k,kb]>0 :
                qd[k,kb,ip]=qd[k,kb,ip]/wd[k,kb]
for kb in range(0,nz):
    if wz[kb]>0 :
        qz[kb]=qz[kb]/wz[kb]
# PBL convection
for ip in range(0,npp):
    for k in range(0,nz):
        for ke in range(0,nz):
            if wp[k,ke]>0 :
                qp[k,ke,ip]=qp[k,ke,ip]/wp[k,ke]
for ke in range(0,nz):
    if wz[ke]>0 :
        qb[ke]=qb[ke]/wz[ke]
# shallow convection     
for ip in range(0,npp):
    for k in range(0,nz):
        for kb in range(0,nz):
            if ws[kb]>0 :
                qs[kb,ip]=qs[kb,ip]/ws[kb]                 
for kb in range(0,nz):
    qs[kb,npp-1]=qs[kb,npp-1]*qs[kb,npp-2]
#                                              
for k in range(0,nz):
    for kb in range(0,nz):
        if wd[k,kb] > 0. :
            td[k,kb]=td[k,kb]/wd[k,kb] 
for kb in range(0,nz):
    if ws[kb] > 0. :
        ts[kb]=ts[kb]/ws[kb]
#
for ic in range(0,ic):
    if sa[ic]>0 :
        for io in range(0,no+2):
            fa[io,ic]=fa[io,ic]*100./sa[ic]
        for io in range(0,no):
            for ke in range(0,nz):
                for kb in range(0,nz):
                    fb[kb,ke,io,ic]=fb[kb,ke,io,ic]*100./sa[ic]
# output
itme="%d "%(sg+0.01)
fcells.write(itme)                     
itme="%s "%"<sg; Frequency of all cloud cells"
fcells.write(itme)
fcells.write('\n')
for ke in range(0,nz):
    itme="%f "%(z[ke]*10.)
    fcells.write(itme)
    for kb in range(0,nz):
        itme="%f "%(ft[kb,ke]*100.) # for ft[kb,ke],the value of ft is cloud cover, 
                                         # z[kb] is base,z[ke]is cloud top
        fcells.write(itme)
    fcells.write('\n')
for kb in range(0,nz):
        itme="%d "%(z[kb]*100.)
        fcells.write(itme)  
fcells.write('\n')
#
for lc in range(0,ne):
    for io in range(0,no):
        itme="%d "%(lc)
        fcells.write(itme) 
        itme="%f "%(sa[lc]+0.01)
        fcells.write(itme) 
        itme="%s "%' Frequency of overlap cloud cells'
        fcells.write(itme)
        itme="%s "%'....... na= '
        fcells.write(itme)
        itme="%d "%io
        fcells.write(itme)
        itme="%s "%'area= '
        fcells.write(itme)
        fcells.write('\n')
        for ic in range(0,no+2):
            itme="%f "%(fa[ic,lc])
            fcells.write(itme)
        fcells.write('\n')
        for ke in range(0,nz):  # sur to top
            for kb in range(0,nz):
                itme="%f "%(fb[kb,ke,io,lc]*100.)
                fcells.write(itme) 
            fcells.write('\n')    
        for kb in range(0,nz):
            itme="%f "%(z[kb]*10.)
            fcells.write(itme)
        fcells.write('\n') 
#
itme="%s "%' ======) CEM deep convection ... samples'
fcells.write(itme)
fcells.write('\n')
for k in range(0,nz):
    itme="%f "%(z[k]*10.)
    fcells.write(itme)
    for kb in range(0,20):
        itme="%f "%(wd[k,kb]*0.1+0.01)
        fcells.write(itme)
    fcells.write('\n')     
for kb in range(0,20):
    itme="%f "%(z[kb]*10.)
    fcells.write(itme)
fcells.write('\n')
itme="%s "%' ======) CEM deep convection ... pot temperature'
fcells.write(itme)
fcells.write('\n')        
for k in range(0,nz):
    itme="%f "%(z[k]*10.)
    fcells.write(itme)
    for kb in range(0,20):
        itme="%f "%(td[k,kb])
        fcells.write(itme)
    fcells.write('\n')     
for kb in range(0,20):
    itme="%f "%(z[kb]*10.)
    fcells.write(itme)
for ip in range(0,npp):
    itme="%s "%' ======) CEM deep convection ... profile, ip='
    fcells.write(itme)
    itme="%d "%ip
    fcells.write(itme)
    fcells.write('\n')         
    for k in range(0,nz):
        itme="%f "%(z[k]*10.)
        fcells.write(itme)
        for kb in range(0,20):
            itme="%f "%(qd[k,kb,ip]*100)
            fcells.write(itme)
        fcells.write('\n')     
    for kb in range(0,20):
        itme="%f "%(z[kb]*10.)
        fcells.write(itme)
    fcells.write('\n')
#
itme="%s "%' ======) CEM PBL convection ... samples'
fcells.write(itme)
fcells.write('\n')
for k in range(0,nz):
    itme="%f"%(z[k]*10.)
    fcells.write(itme)
    for kb in range(0,30):
        itme="%f "%(wp[k,kb]*0.1+0.01)
        fcells.write(itme)
    fcells.write('\n')    
for kb in range(0,30):
    itme="%f "%(z[kb]*10.)
    fcells.write(itme)
fcells.write('\n') 
for ip in range(0,npp):
    itme="%s "%' ======) CEM PBL convection ... profile, ip='
    fcells.write(itme)
    itme="%d "%ip
    fcells.write(itme)
    fcells.write('\n')         
    for k in range(0,nz):
        itme="%f "%(z[k]*10.)
        fcells.write(itme)
        for ke in range(0,30):
            itme="%f "%(qp[k,ke,ip]*100)
            fcells.write(itme)
        fcells.write('\n')     
    for ke in range(0,30):
        itme="%f "%(z[ke]*10.)
        fcells.write(itme)
    fcells.write('\n')
#
itme="%s "%' ======) CEM shallow convection ... profile'
fcells.write(itme)
fcells.write('\n')
itme="%s "%' z base  sample  profiles'
fcells.write(itme)
fcells.write('\n')
for kb in range(0,nz):
    item="%f "%(z[kb]*10.)
    fcells.write(item)
    item="%f "%(ws[kb]+0.01)
    fcells.write(item)
    item="%f "%(tr[kb])
    fcells.write(item)
    item="%f "%(ts[kb])
    fcells.write(item)    
    for ip in range(0,npp):
        item="%f "%(qs[kb,ip]*10000.)
        fcells.write(item)
    fcells.write('\n')
fcells.close() 
#------  1-d plot
#   deep convection 
for ip in range(0,npp):
    for kb in range(0,nz):
        qzd=0.1*(kb-1.)
        for k in range(0,nz):
            ww[k,kb,0]=(qd[k,kb,ip]+qzd)*100.
            if wd[k,kb] <= 0 :
                ww[k,kb,0]=spv
#
#Plotting








# pbl convection
for ip in range(0,npp):
    for ke in range(0,nz):
        qzd=0.1*(ke-1.)
        for k in range(0,nz):
            ww[k,kb,0]=(qp[k,ke,ip]+qzd)*100.
            if wd[k,ke] <= 0 :
                ww[k,ke,0]=spv

#


# shallow convection
for ip in range(0,npp):
    for kb in range(0,nz):
        cm[kb]=qs[kb,ip]
        if ws[kb] <= 0 :
            cm[kb]=spv


#
tmpout=np.ndarray(shape=(nz), dtype=float)
for kb in range(0,nz):
    tmpout[kb]=ts[kb]
    if ws[kb] <= 0. :
        tmpout[kb]=spv
#
        
        
#
for kb in range(0,nz):
    tmpout[kb]=qz[kb]
    if wz[kb] <= 0. :
        tmpout[kb]=spv
#
        
        
#       
for kb in range(0,nz):
    tmpout[kb]=qb[kb]
    if wb[kb] <= 0. :
        tmpout[kb]=spv



#########################################################################
if zbmin >0. :
    for kb in range(0,nz):
        if z[kb] < zbmin :
            for ke in range(0,nz):
                ft[kb,ke]=0
                for lc in range(0,ne):
                    for io in range(0,no):
                        fb[ke,ke,io,lc]=0.
nrec=0
#
#plot  ----------  ft  # xdat,ydat,zdat
font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 18,
        }  
cloudlevs=[2,5,10,15,20,30,40,50,60,70,80,90,100,110]
cloudclors=['w','lightgray','plum','darkorchid','b','dodgerblue','skyblue','aqua',
            'lime','greenyellow','yellow','salmon','pink','orangered','r','darkred']
fig,axe1=plt.subplots(nrows=1,ncols=1,figsize=(6,6))
plt.subplot(1,1,1)
#zdat[0,:]=0.0   ## the first level is below surface ground
ft0=np.ndarray(shape=(nz,nz), dtype=float) #(km,km)  For exchange the dims
for i1 in range(0,nz):
    for i2 in range(0,nz):
        ft0[i1,i2]=ft[i2,i1]
titlename=r"Frequency of all cloud cells ($10^{-2}%$)"
axe1=plt.contourf(z,z,ft0,colors=cloudclors, levels=cloudlevs,extend='both')
plt.colorbar(orientation='horizontal',extend='both',
    extendfrac='auto',  spacing='uniform')                           
plt.title(titlename,fontsize=16)                          
plt.axis([0, 16, 0, 16])
plt.xlabel(r'Cloud Base Height ($km$)', fontdict=font)
plt.ylabel(r'Cloud Top Height ($km$)', fontdict=font)
plt.show()
plt.savefig(dirpic+casenm+"_CloudCellsTopBase_20120520.pdf")          
plt.show()
plt.close()   
#fb
#
for ip in range(0,npp):
    for kb in range(0,nz):
        if wd[kb,kb] <= 0. :
            for k in range(0,nz):
                ww[k,kb,0]=spv
        else:
            for k in range(0,nz):
                ww[k,kb,0]=qd[k,kb,ip] #np 
#
#               
#
for ip in range(0,npp):
    for kb in range(0,nz):
        if wd[kb,kb] <= 0. :
            for k in range(0,nz):
                ww[k,kb,1]=spv
        else:
            for k in range(0,nz):
                ww[k,kb,1]=qd[k,kb,ip]
#
#
#
for kb in range(0,nz):
    for k in range(0,nz):
        ww[k,kb,1]=td[k,kb]
        if wd[k,kb]<=0. :
            ww[k,kb,1]=spv
            
##
#
#
for ip in range(0,npp):
    for ke in range(0,nz):
        if wp[kb,kb] <= 0. :
            for k in range(0,nz):
                ww[k,ke,1]=spv
        else:
            for k in range(0,nz):
                ww[k,ke,1]=qp[k,kb,ip]