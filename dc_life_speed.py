#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat May 02 23:01:24 2015

@author: jhchen
"""
import string
import numpy as np

def readAscii(fpath,iskp,nrl):
    #iskp  the total line skipped of the file
    # fpath   the full path of the file
    # usage: onedim=readAscii(fpaht,iskp)
    onedim=[]
    linesplit=[]
    f=open(fpath)
    ff=f.readlines()[iskp:nrl]  ## first line in obs file is legend 
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
            kb(na)=k1
            ke(na)=k2
            cm(na)=aa
            aa=0.
            na=na+1
    if c[km-1] > 0.0 :
        kb(na)=k1
        ke(na)=k2
        cm(na)=aa
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
            kb(na)=k1
            ke(na)=k2
#            cm(na)=aa
            aa=0.
            na=na+1
    if c[km-1] > 0.0 :
        kb(na)=k1
        ke(na)=k2
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
def shapecc(qw,im,km):
    c=np.ndarray(shape=(im,km), dtype=float)
    kzb=np.ndarray(shape=(im,99), dtype=int)
    kze=np.ndarray(shape=(im,99), dtype=int)
    """ kxl=np.ndarray(shape=(km,im+1), dtype=int) #from west  left
    kxr=np.ndarray(shape=(km,im+1), dtype=int) # right
    naz=np.ndarray(shape=(im), dtype=int)
    nax=np.ndarray(shape=(im), dtype=int)
    for i in range(0,im):
        for k in range(0,km):
            if qw[i,k] >= 1.0e-3:
                c[i,k]=1. # cloud
            else:
                c[i,k]=0.
        kzb[i,:],kze[i,:],cm,na[i]=infcld(c[i,:],km)
    for k in range(0,km):
        kxl[i,:],kxr[i,:],nax[k]=infcld(c[:,k],im) 
    """
    dcell=np.ndarray(shape=(2,im+1), dtype=int) # position number
    qmx=np.ndarray(shape=(im+1), dtype=float)
    qpmx=np.ndarray(shape=(2,im+1), dtype=int)
    ncl=0
    for ilp in range(0,nx):
        if ix==0 :
            ixs=ilp
        x1,k1,k2=fc1(kzb,kze,ixs)
        qmax,pqmax=gmax(qw,x1,k1,k2)
        if x1>0 and x1 <= nx :
           qmax,pqmax,ixl,ixr=fc2(kzb,kze,qw,qmax,pqmax,k1,k2,x1)
           dcell[0,ncl]=ixl
           dcell[1,ncl]=ixr
           qpmx[0,ncl]=pqmax[0]
           qpmx[1,ncl]=pqmax[1]
           qmx=[ncl]=qmax
           ncl=ncl+1
           ixs=ixr+1
           if ixs >=nx :
               break                
        elif x1<0:
            ixs=ixs+1
    return dcell,qmx,qpmx,ncl
def fc1(kzb,kze,ix):
    x1=-1
    k1=-1
    k2=-1
    for k in range(0,99):
        if kzb[ix,k] <=5 and kze[ix,k] > 18 :  # deep convection origin:24 9.9km,  19 6.5km Luo et al 2010
            x1=ix
            k1=kzb[ix,k] 
            k2=kze[ix,k]
            break
    return x1,k1,k2
def gmax(qw,x1,k1,k2):
    nx=len(qw[:,0])
    nz=len(qw[0,:])
    pqmax=np.ndarray(shape=(2), dtype=int)  # postion of the maximun
    pqmax[0]=x1
    pqmax[1]=k2
    qmax=qw[x1,k2]
    for k in range(k1,k2):
        if qw[x1,k] > qmax:
            pqmax[1]=k
            qmax=qw[x1,k]
    return qmax,pqmax
def fc2(kzb,kze,qw,qmax,pqmax,k1,k2,ix):
    nx=len(qw[:,0])
    xl=ix
    xr=ix
    xl0=ix
    xr0=ix+1 
    x1=-1
    ke=k2
    kb=k1
    deep=ke-kb
    fco=0.5  # the vertical coincident part of cloud between two nearby profiles
    for ixx in range(x1-1,-1,-1): # left
        icount=0
        for k in range(0,99):
            kb0 = kzb[ixx,k]
            ke0 = kze[ixx,k]
            deep0=ke0-kb0
            if kb0 != ke0 :
                if kb0 <=5 and ke0 > 18 :  # deep convection
                    qmax0,pqmax0=gmax(qw,ixx,kb0,ke0)
                    if qmax0 > qmax :
                        pqmax[0]=pqmax0[0]
                        pqmax[1]=pqmax0[1]
                        qmax=qmax0
                    kb=kb0
                    ke=ke0
                    deep=ke-kb
                    icount=icount+1
                elif deep0 >= int(0.8*deep) and kb0 < ke :
                    if (kb0 <= kb and (ke0 -kb)>= deep0*fco) or kb0>kb and (ke0-ke)<deep0*(1-fc0):
                        qmax0,pqmax0=gmax(qw,ixx,kb0,ke0)
                        if qmax0 > qmax :
                            pqmax[0]=pqmax0[0]
                            pqmax[1]=pqmax0[1]
                            qmax=qmax0
                        kb=kb0
                        ke=ke0
                        deep=ke-kb
                        icount=icount+1
        if not(icount>0) :
            xl=ixx
            break
#
    for ixx in range(x1+1,nx): # right
        icount=0
        for k in range(0,99):
            kb0 = kzb[ixx,k]
            ke0 = kze[ixx,k]
            deep0=ke0-kb0
            if kb0 != ke0 :
                if kb0 <=5 and ke0 > 18 :  # deep convection
                    qmax0,pqmax0=gmax(qw,ixx,kb0,ke0)
                    if qmax0 > qmax :
                        pqmax[0]=pqmax0[0]
                        pqmax[1]=pqmax0[1]
                        qmax=qmax0
                    kb=kb0
                    ke=ke0
                    deep=ke-kb
                    icount=icount+1
                elif deep0 >= int(0.8*deep) and kb0 < ke :
                    if (kb0 <= kb and (ke0 -kb)>= deep0*fco) or kb0>kb and (ke0-ke)<deep0*(1-fc0):
                        qmax0,pqmax0=gmax(qw,ixx,kb0,ke0)
                        if qmax0 > qmax :
                            pqmax[0]=pqmax0[0]
                            pqmax[1]=pqmax0[1]
                            qmax=qmax0
                        kb=kb0
                        ke=ke0
                        deep=ke-kb
                        icount=icount+1
        if not(icount>0) :
            xr=ixx
            break
    return qmax,pqmax,xl,xr
def dcspeed(dcell,cloud,qmx,qpmx,ncl):  # dcell horizontal 
    nt=len(dcell[0,0,:])
    nclmx=len(dcell[0,:,0])
    nx=nclmx-1
    maxspeed= 60  # the convection system maximum speed
    maxlife=4   # max life hours
    dt=15       # mins  output time
    delt=maxlife*60/dt  # seach period
    cellsp=np.ndarray(shape=(2,nt,nt), dtype=int) #
    ptm=[]   # record the time of erever cell
    idx1=0
    count=0
    dctopbs=[18,5] # deep convection top and base
    markcell=np.ndarry(shape=(nx,nt),dtype=int) # MARKED THE cells that has been readed
    fout=f.open('dccells.txt')
    item="%s "%"CellCount TimeFromBorned left right qmax qmax_x qmax_z ts"
    fout.write(item)
    fout.write('\n')
    for it in range(0,nt):
        if ncl[it]>0 : # the first snapshot that has dc cells
            nc0=ncl[it]
            its=it
            ite=it+delt            
            for i0 in range(0,nc0):
                markcell[i0,it]=1
                count=count+1
                contm=1
                pxmax=cellsp[0,i0,it]
                pzmaz=cellsp[1,i0,it]
                xlf=dcell[0,i0,it]  #left
                xrf=dcell[1,i0,it]  # right
                qmax=qmx[i0,it]
                cell0=[qmax,pxmax,pzmaz,xlf,xrf]
                outstring=getoutstring(count,contm,qmax,pxmax,pzmaz,xlf,xrf,it)
                fout.write(outstring)
                fout.write('\n')                

def gettc1(cloud,it,ix,bstop): # cloud nx,nz,nt
    bs=bstop[0]
    top=bstop[1]
    a=0
    for iz in range(bs,top):
        a=a+cloud[ix,iz,it]
    return a            
def getoutstring(count,contm,qmax,pxmax,pzmaz,xlf,xrf,it):
    item="%s "%count+"%s "%contm+"%s "%xlf+"%s "%xrf+"%s "%qmax+"%s "%pxmax+"%s "%pzmaz+"%s "%it
    return item
def dected(its,ite,ncl,cloud,cellsp,dcell,qmx,dctopbs,cell0,markcell):
    qmax0=cell0[0], pxmax0=cell0[1], pzmaz0=cell0[2]
    xlf0=cell[3]  , xrf0=cell[4]
    itt=its+1
    maxspeed= 60  # the convection system maximum speed  km per hour
    dt=15       # mins  output time
    res=3 # reslution km
    deltx=dt/60.*maxspeed/3  # the maxmuim x grids ervey step
    pxmax1=cellsp[0,0,itt]
    pzmaz1=cellsp[1,0,itt]
    xlf1=dcell[0,0,itt]  #left
    xrf1=dcell[1,0,itt]  # right
    qmax1=qmx[0,itt]
    nclmx=len(dcell[0,:,0])
    nx=nclmx-1
    a0=gettc1(cloud,its,pxmax0,dctopbs) #
    tmcont=0
    outstring=[]
    icn=1
    for itt in range(its+1,ite):
        nc1=ncl[itt]
        if nc1>0 and icn >0: # there is somedc cells
            icn=-1
            a1=gettc1(cloud,itt,pxmax1,dctopbs)
            if pxmax0<nx-deltx and pxmax1>deltx : #the leftest dc cell of nt+1 dont reach the left boundary
                tmc,xmax,zmax,qmax,xlft,xrgh,ixx=normal(nc1,itt,cloud,dctopbs,dcell,cellsp,qmx,qmax0,pxmax0,pzmaz0,xlf0,xrf0)
                if ixx>=0 :                
                    tmcont=tmcont+tmc
                    outstrs=getoutstring(count,tmcont,qmax,xmax,zmaz,xlft,xrgh,itt)
                    outstring.append(outstrs)
                    qmax0=qmax, pxmax0=xmax, pzmaz0=zmaz
                    xlf0=xlft , xrf0=xrgh
                    pxmax1=cellsp[0,0,itt+1]
                    pzmaz1=cellsp[1,0,itt+1]
                    xlf1=dcell[0,0,itt+1]  #left
                    xrf1=dcell[1,0,itt+1]  # right
                    qmax1=qmx[0,itt+1]
                    a0=gettc1(cloud,itt,pxmax0,dctopbs) #
                    markcell[ixx,itt]=1
                    icn=1
            else:
                tmc,ixx,itx=boundy(nc1,itt,cloud,dctopbs,qmx,qmax0,pxmax0,pzmaz0,xlf0,xrf0,  \
                    qmax1,qmax1,pxmax1,pzmaz1,xlf1,xrf1,deltx,markcell)
                if ixx==1:
                    tmcont=tmcont+tmc
                    outstrs=getoutstring(count,tmcont,qmax,xmax,zmaz,xlft,xrgh,itt)
                    outstring.append(outstrs)
                    qmax0=qmax, pxmax0=xmax, pzmaz0=zmaz
                    xlf0=xlft , xrf0=xrgh
                    pxmax1=cellsp[0,0,itt+1]
                    pzmaz1=cellsp[1,0,itt+1]
                    xlf1=dcell[0,0,itt+1]  #left
                    xrf1=dcell[1,0,itt+1]  # right
                    qmax1=qmx[0,itt+1]
                    a0=gettc1(cloud,itt,pxmax0,dctopbs) #
                    markcell[ixx,itt]=1
                    icn=1

        else:        # there is no dc cell  at this time
            tmcont=0
            break

    return                   
def boundy(nc1,itt,cloud,dctopbs,qmx,qmax0,pxmax0,pzmaz0,xlf0,xrf0,   \
    qmax1,qmax1,pxmax1,pzmaz1,xlf1,xrf1,deltx,markcell):
    #
    tmcn=0
    ixx=-1
    deltx=dt/60.*maxspeed/3  # the maxmuim x grids ervey step
    if (qxmax0 >nx-deltx and qxmax1<deltx) and markcell[0,itt+1] <1 and markcell[0,itt] <1 :
        if abs((xrf0-xlf0)-(xrf1-xlf1))< 0.2*(xrf0-xlf0)  and \
           (nx-pxmax0+pxmax1+2)<deltx                    and \
           abs(qmax1-qmax0)<0.2*qmax0   :
           tmcn=1
           ixx=1
           itx=itt  
    elif xrf0 <nx and xlf1<1 and xrf0>(nx-deltx) :
            
    elif xrf0==nx and xlf1<1 :   
        ix=0
        a1=gettc1(cloud,itt+1,ix,dctopbs) # 
        ix=nx
        a0=gettc1(cloud,itt,ix,dctopbs) # 
        if abs(a1-a0)<a0*0.3 and abs(qmax1-qmax0) <qmax0*0.3 and (nx-pxmax0+pxmax1+2)<deltx :
           tmcn=1
           ixx=2
           itx=itt
    return tmcn,ixx,itx
def normal(nc1,itt,cloud,dctopbs,dcell,cellsp,qmx,qmax0,pxmax0,pzmaz0,xlf0,xrf0):
    tmcont=0
    ixx=-1
    for i1 in range(0,nc1):
        x1=cellsp[0,i1,itt]
        pzmaz1=cellsp[1,i1,itt]
        xlf1=dcell[0,i1,itt]  #left
        xrf1=dcell[1,i1,itt]  # right
        qmax1=qmx[i1,itt]
        a1=gettc1(cloud,itt,x1,dctopbs) #
        if xrf0-xlf1>0 and abs(a1-a0)<a0*0.3 and abs(qmax1-qmax0) <qmax0*0.3 :
            tmcont=tmcont+1
            xmax=x1,zmax=pamaz1
            xlft=xlf1,xrgh=xrf1
            qmax=qmax1
            ixx=i1
            break
    return tmcont,xmax,zmax,qmax,xlft,xrgh,ixx
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
casenm="ETPCTR_EC"
dirin =""
dirout=""
nd=30 # days
dt=15 # output timestep in min
nx=200
nz=34
nt=nd*24*60/dt
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
np=3
ne=1
fout=dirout+""
fdc=open(fout,'w')
fout=dirout+""
fcells=open(fout,'w')
fout=dirout+""
fo1=open(fout,'w')
fout=dirout+""
fo2=open(fout,'w')
#
qz=np.ndarray(shape=(nz), dtype=float)
wz=np.ndarray(shape=(nz), dtype=float)
qd=np.ndarray(shape=(nz,nz,np), dtype=float)
td=np.ndarray(shape=(nz,nz), dtype=float)
wd=np.ndarray(shape=(nz,nz), dtype=float)
#
qb=np.ndarray(shape=(nz), dtype=float)
wb=np.ndarray(shape=(nz), dtype=float)
qp=np.ndarray(shape=(nz,nz,np), dtype=float)
tp=np.ndarray(shape=(nz,nz), dtype=float)
wp=np.ndarray(shape=(nz,nz), dtype=float)
#
qs=np.ndarray(shape=(nz,np), dtype=float)
ts=np.ndarray(shape=(nz), dtype=float)
ws=np.ndarray(shape=(nz), dtype=float)
#
fa=np.ndarray(shape=(no+2,ne), dtype=float)
fb=np.ndarray(shape=(nz,nz,no,ne), dtype=float)
sa=np.ndarray(shape=(ne), dtype=float)
cels=[0.00,1.00]
#
ft=np.ndarray(shape=(nz,nz), dtype=float) #(km,km,ng)
tr=np.ndarray(shape=(nz), dtype=float) #(km,km,ng)
#
dcell=np.ndarray(shape=(2,nx+1,nt), dtype=int) # left and right boundarr  
qmx=np.ndarray(shape=(nx+1,nt), dtype=float)   # nx+1 the possible  maximum cells of every snapshot 
pqmx=np.ndarray(shape=(2,nx+1,nt), dtype=int) # max position
ncl=np.ndarray(shape=(nt), dtype=float)
cloud=np.ndarray(shape=(nx,nz,nt), dtype=int) # left and right boundarr  
#       
fpath=dirin+""
iskp=0
sg=0.
for it in range(0,nt):
    iskp=it*76   # (200+200+200+1)/8
    nrlines=76
    onedim=readAscii(fpath,iskp,nrlines)
    ql=np.ndarray(shape=(nx,nz), dtype=float)  # liquid
    qi=np.ndarray(shape=(nx,nz), dtype=float)  # ice
    tc=np.ndarray(shape=(nx,nz), dtype=float)  # temp
    den=np.ndarray(shape=(nz), dtype=float)  # density
    qw=np.ndarray(shape=(nx,nz), dtype=float)  # total
    lcb=np.ndarray(shape=(nz), dtype=int)  # density
    lce=np.ndarray(shape=(nz), dtype=int)  # density
    for iz in range(0,nz):
        k=iz*601
        den[iz]=onedim[k]
        for ix in range(0,nx):
            k=iz*601+ix
            ql[ix,iz]=onedim[k]
            k=iz*601+ix+200
            qi[ix,iz]=onedim[k]
            k=iz*601+ix+400
            tc[ix,iz]=onedim[k]
            qw[ix,iz]=ql[ix,iz]+qi[ix,iz]
    del onedim
    zd=np.ndarray(shape=(nz), dtype=float)
    qc0=np.ndarray(shape=(nz), dtype=float)
    zd[0]=0.5*(z[1]-z[0])*den[0]
    for k in range(1,nz-1):
        zd[k]=0.5*(z[k+1]-z[k-1])*den[k]
    zd[nz-1]=0.5(z[nz-1]-z[nz-2])*den[nz-1]
    if cwp > 0.0 :
        for k in range(0,nz):
            qc0[k]=1e-3*cwp/zd(k)
    elif qcc > 0.0 :
        for k in range(0,nz):
            qc0[k]=qcc
    else:
        print "xxxx wrong cloud id threshold"
        exit()
    szd=0.0
    for iz in range(0,nz):
        szd=szd+zd[iz]
    kdb,kde,c,cm,na =deepcc(qw,nx,nz)
    cloud[:,:,it]=c[:,:]
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
                fdc.write('\n')
#
    dcell[:,:,it],qmx[:,it],ncl[it] =shapecc(qw,nx,nz)
#
    frc=1./(nx*nz)
    frg=1./nx
    frr=1./nx*nt
    for k in range(0,nz):
        for ix in range(0,nx):
            tr[k]=tr[k]+frr*tc[ix,k]
    sumf=0.0
    ezd=0.0
    qmx=0.0
    for i in range(0,nx):
        for k in range(0,nz):
            sumf=sumf+frc*qw[i,k]
            ezd=ezd+frg*qw[i,k]*zd[k]
            qmx=max(qmx,qw[i,k])
    ezd=ezd/szd
    #--------get pdf for cloud ensemble
    if sumf > qce :
        sg=sg+1.
        wa=np.ndarray(shape=(no+1), dtype=float)
        ww=np.ndarray(shape=(nz,nz,no), dtype=float)
        kcb=np.ndarray(shape=(im,99), dtype=int)
        kce=np.ndarray(shape=(im,99), dtype=int)
        cm=np.ndarray(shape=(im,99), dtype=float)
        ns=np.ndarray(shape=(im), dtype=float)
        cc=np.ndarray(shape=(nx,nz), dtype=float)
        ce=0.
        for ix in range(0,nx):
            for iz in range(0,nz):
                cc[ix,iz]=qw[ix,iz]
                if cc[ix,iz] < qc0[iz] :
                    cc[ix,iz]=0.
            kcb[ix,:],kce[ix,:],ecm[ix,:],ns[ix]=infcld(cc[ix,:],nz)       
            for ii in range(0,ns[ix]):
                kb=kcb[ix,ii]
                ke=kce[ix,ii]
                ft[kb,ke]=ft[kb,ke]+frg
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
                        qd[k,kb,1]=qd[k,kb,1]+cc[iz,k]
                        qd[k,kb,2]=qd[k,kb,2]+cc[iz,k]/qzd
                        qd[k,kb,3]=qd[k,kb,3]+cc[iz,k]/qzd*den[k]
                        td[k,kb]=td[k,kb]+tc(ix,k)
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
                        qp[k,ke,1]=qp[k,ke,1]+cc[iz,k]
                        qp[k,ke,2]=qp[k,ke,2]+cc[iz,k]/qzd
                        qp[k,ke,3]=qp[k,ke,3]+cc[iz,k]/qzd*den[k]
                        tp[k,ke]=tp[k,ke]+tc(ix,k)
                        wp[k,ke]=wp[k,ke]+1.
#c...............vertical profiles for shallow convections                
            if ns[ix]==1:
                ii=0
                kb=kcb[ix,ii]
                ke=kce[ix,ii]
                if ke-kb < ksc :
                    czd=0.
                    qzd=0.
                    tzd=0.
                    for k in range(kb,ke):
                        czd=czd+zd[k]
                        qzd=qzd+zd[k]*cc[ix,k]
                        tzd=tzd+zd[k]*tc[ix,k]
                    qzd=qzd/czd
                    tzd=tzd/czd
                    qs[kb,1]=qs[kb,1]+qzd
                    qs[kb,2]=qs[kb,2]+ezd
                    qs[kb,3]=qs[kb,3]+qzd/ezd*den[k]
                    ts[kb]=ts[kb]+tzd
                    ws[kb]=ws[kb]+1.
#...............for overlap cloud cells
#               ----ignore thin cells
            lc=0
            for ii in range(0,ns[ix]):
                if (kce[ix,ii]-kcb[ix,ii]) >= 2 :
                    lce[lc]=kce[ii]
                    lcb[lc]=kcb[ii]
                    lc=lc+1
            #----ignore small gaps
            if lc >1 :
                for ii in range(1,lc):
                    if lcb[ii]-lce[ii-1]<2 :
                        lce[ii-1]=lce[ii]
                        lc=lc-1
                        for ic in range(ii,lc):
                            lcb[ic]=lcb[ic+1]
                            lce[ic]=lce[ic+1]
            #----get frequency according to overlap cells
            if lc >0 :
                ce=ce+frg
                io=min(lc,no+1)
                wa[io]=wa[io]+frg
            if lc == 1 :
                kb=lcb[1] # bas of the cell
                ke=lce[1] # top of the cell                
            elif lc == 2 :
                kb=lce[1]  # top of low cell
                ke=lcb[2]  # bas of upp cell                        
            elif lc == 3 :
                kb=lcb[2]  # top of mid cell
                ke=lce[2]  # bas of mid cell  
            if lc>=1 and lc<=3 :
                ww[kb,ke,lc]=ww[kb,ke,lc]+frg
        wa[0]=1.0-ce # clear sky
        for ii in range(1,ne+1):
            if cels[ii-1] < ce and ce <=cels[ii] :
                ic=ii
        sa[ic]=sa[ic]+1.
        for io in range(0,no+2):
            fa[io,ic]=fa[io,ic]+wa[io]
        for io in range(0,no):
            for ke in range(0,nz):
                for kb in range(0,nz):
                    fb[kb,ke,io,ic]=fb[kb,ke,io,ic]+wwfb[kb,ke,io]
# end loop of it

if sg > 0. :
    for kb in range(0,nz):
        for kb in range(0,nz):
            ft[kb,ke,ig]=  ft[kb,ke,ig]*100./sg                              
# deep convection
for ip in range(0,np):
    for k in range(0,nz):
        for kb in range(0,nz):
            if wd[k,kb]>0 :
                qd[k,kb,ip]=qd[k,kb,ip]/wd[k,kb]
for kb in range(0,nz):
    if w[kb]>0 :
        qz[kb]=qz[kb]/wz[kb]
# PBL convection
for ip in range(0,np):
    for k in range(0,nz):
        for ke in range(0,nz):
            if wp[k,ke]>0 :
                qp[k,ke,ip]=qp[k,ke,ip]/wp[k,kp]
for ke in range(0,nz):
    if w[ke]>0 :
        qb[ke]=qb[ke]/wz[ke]
# shallow convection     
for ip in range(0,np):
    for k in range(0,nz):
        for kb in range(0,nz):
            if ws[wb]>0 :
                qs[kb,ip]=qd[kb,ip]/ws[kb]                 
for kb in range(0,nz):
    qs[kb,np]=qs[kb,np]*qs[kb,np-1]
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
    itme="%d "%(z[ke]*10.)
    fcells.write(itme)
    for kb in range(0,nz):
        itme="%d "%(ft[kb,ke,ig]*100.)
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
        itme="%d "%(sa[lc]+0.01)
        fcells.write(itme) 
        itme="%s "%' Frequency of overlap cloud cells'
        fcells.write(itme)
        fcells.write('\n')
        itme="%s "%'....... na= '
        fcells.write(itme)
        itme="%d "%io
        fcells.write(itme)
        itme="%s "%'area= '
        fcells.write(itme)
        for ic in range(0,no+2):
            itme="%f "%(fa[ic,lc])
            fcells.write(itme)
        fcells.write('\n')
        for ke in range(0,nz):  # sur to top
            for kb in range(0,nz):
                itme="%d "%(fb[kb,ke,io,lc]*100.)
                fcells.write(itme) 
            fcells.write('\n')    
        for kb in range(0,nz):
            itme="%d "%(z[kb]*10.)
            fcells.write(itme)
        fcells.write('\n') 
#
itme="%s "%' ======) CEM deep convection ... samples'
fcells.write(itme)
fcells.write('\n')
for k in range(0,nz):
    itme="%d "%(z[k]*10.)
    fcells.write(itme)
    for kb in range(0,20):
        itme="%d "%(wd[k,kb]+0.01)
        fcells.write(itme)
    fcells.write('\n')     
for kb in range(0,20):
    itme="%d "%(z[kb]*10.)
    fcells.write(itme)
fcells.write('\n')
itme="%s "%' ======) CEM deep convection ... pot temperature'
fcells.write(itme)
fcells.write('\n')        
for k in range(0,nz):
    itme="%d "%(z[k]*10.)
    fcells.write(itme)
    for kb in range(0,20):
        itme="%d "%(td[k,kb])
        fcells.write(itme)
    fcells.write('\n')     
for kb in range(0,20):
    itme="%d "%(z[kb]*10.)
    fcells.write(itme)
for ip in range(0,np):
    itme="%s "%' ======) CEM deep convection ... profile, ip='
    fcells.write(itme)
    itme="%d "%ip
    fcells.write(itme)
    fcells.write('\n')         
    for k in range(0,nz):
        itme="%d "%(z[k]*10.)
        fcells.write(itme)
        for kb in range(0,20):
            itme="%d "%(qd[k,kb,ip]*100)
            fcells.write(itme)
        fcells.write('\n')     
    for kb in range(0,20):
        itme="%d "%(z[kb]*10.)
        fcells.write(itme)
    fcells.write('\n')
#
itme="%s "%' ======) CEM PBL convection ... samples'
fcells.write(itme)
fcells.write('\n')
for k in range(0,nz):
    itme="%d "%(z[k]*10.)
    fcells.write(itme)
    for kb in range(0,30):
        itme="%d "%(wp[k,kb,ip]+0.01)
        fcells.write(itme)
    fcells.write('\n')    
for kb in range(0,30):
    itme="%d "%(z[kb]*10.)
    fcells.write(itme)
fcells.write('\n') 
for ip in range(0,np):
    itme="%s "%' ======) CEM PBL convection ... profile, ip='
    fcells.write(itme)
    itme="%d "%ip
    fcells.write(itme)
    fcells.write('\n')         
    for k in range(0,nz):
        itme="%d "%(z[k]*10.)
        fcells.write(itme)
        for ke in range(0,30):
            itme="%d "%(qp[k,ke,ip]*100)
            fcells.write(itme)
        fcells.write('\n')     
    for ke in range(0,30):
        itme="%d "%(z[ke]*10.)
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
    item="%d "%(z[kb]*10.)
    fcells.write(itme)
    item="%d "%(ws[kb]+0.01)
    fcells.write(itme)
    item="%d "%(tr[kb])
    fcells.write(itme)
    item="%d "%(ts[kb])
    fcells.write(itme)    
    for ip in range(0,np):
        item="%d "%(qs[kb,ip]*10000.)
        fcells.write(itme)
    fcells.write('\n')
fcells.close()
fdc.close()    
    