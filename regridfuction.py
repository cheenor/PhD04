#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 09 13:00:39 2015

@author: jhchen
"""
import numpy as np
def gridint(ain,nx,nz,xx,zz,nx1,nz1,xx1,zz1,x0,z0):
#cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#cc this subroutine performs interpolation of the data given in the
#cc input array ain(nx,nz) on the grid defined by xx(nx) and zz(nz)
#cc into the grid given by xx1(nx1) and zz1(nz1). NIETHER OF THE
#cc GRIDS HAS TO BE REGULAR. Data is returned in the ain(nx1,nz1)
#cc part of the input array.
#cc 
#cc    levels in the input array are given in zz(nz), 
#cc    levels in the output array are given in zz1(nz1)
#cc      x-coordinate in the input array are in xx(nx)
#cc      x-coordinate in the output array are in xx1(nx1)
#cc        x0(nx,nz) and z0(nx,nz) are working arrays
#cc
#cc NOTE that nx1 (or nz1) must be smaller than nx (or nz) and xx1 (zz1)
#cc  must be a subdomain of xx (zz)
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#      dimension ain(nx,nz),zz(nz),xx(nx)
#      dimension zz1(nz1),xx1(nx1) 
#      dimension x0(nx,nz),z0(nx,nz)
#
#c      check consistency of the input data
    ier=0
#cc nx1,nz1 not larger than nx,nz
    if nx1 > nx:
        ier=1
    if nz1 > nz :
        ier=2
#cc limits (zz1(nz1).le.zz(nz) ?) 
    if zz1[0] < zz[0] :
        ier=3
    if zz1[nz1]>zz[nz] :
        ier=4
#cc limits (xx1(nx1).le.xx(nx) ?) 
    if xx1[1] < xx[1]:
        ier=5
    if xx1[nx1] > xx[nx] :
        ier=6
    if ier>0 :
        print 'problems with input data. will stop.',' ier = ',ier
    nxz=nx*nz
    for  i in range(0,nxz):
        z0[i,1]=1.
        x0[i,1]=1.
#cc  map vertical grid positions:
    for k1 in range(0,nz1):
        zzh=zz1[k1]
        for k in range(0,nz):
            kk=k
            if zz[k] > zzh :
                break
        kkm=max0(1,kk-1)
        z0[1,k1]=float(kkm)+(zzh-zz[kkm])/(zz[kk]-zz[kkm]+1.e-6)
    do 3 i1=2,nx1
      do 3 k1=1,nz1
  3   z0(i1,k1)=z0(1,k1)
c
cc  map horizontal grid positions:
      do 11 i1=1,nx1
      xxh=xx1(i1)
      do 12 i=1,nx
      ii=i
      if(xx(i).ge.xxh) go to 16
 12   continue
 16   iim=max0(1,ii-1)
      x0(i1,1)=float(iim)+(xxh-xx(iim))/(xx(ii)-xx(iim)+1.e-6)
 11   continue
      do 13 i1=1,nx1
      do 13 k1=2,nz1
 13   x0(i1,k1)=x0(i1,1)
cc
cc  call Piotr's interpolation routine
      call inter2(ain,x0,z0,nx,nz)
cc
      return
      end
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE INTER2(XF,XD1,XD2,NX,NZ)
C IOR=ORDER OF ACCURACY/2; ONLY EVEN ORDER TRMBACK SCHEMES ARE CONSIDERED
      PARAMETER(IOR=2)
      PARAMETER(LINER=0)
      DIMENSION XF(*),XD1(*),XD2(*)
CC  N1 - HORIZONTAL INDEX, N2 - VERTICAL INDEX
      PARAMETER (N1=57, N2=41,NN=N1*N2)
      DIMENSION Z(NN,-IOR:IOR)
      DATA  EP/ 1.E-10/
      PARAMETER(NONOS=1)
      REAL      MX,MN
      PARAMETER(IBC=0)
      COMMON // IG0(NN),JG0(NN),X(-IOR+1:N1+IOR,-IOR+1:N2+IOR)
C  next is for shavano: 
C      DONOR(Y1,Y2,A)=CVMGM(Y2,Y1,A)*A
C  next is for workstation:
      DONOR(Y1,Y2,A)=AMAX1(0.,A)*Y1 + AMIN1(0.,A)*Y2
      TR2(Y1,Y2,A)=A*.5*(Y1+Y2)-A**2*.5*(Y2-Y1)
      TR4(YM1,Y0,YP1,YP2,A)=A/12.*(7.*(YP1+Y0)-(YP2+YM1))
     1 -A**2/24.*(15.*(YP1-Y0)-(YP2-YM1))-A**3/12.*((YP1+Y0)
     2 -(YP2+YM1))+A**4/24.*(3.*(YP1-Y0)-(YP2-YM1))
      TR6(YM2,YM1,Y0,YP1,YP2,YP3,A)=-A/60.*(-YM2+8.*YM1-37.*Y0
     1                                     -37.*YP1+8.*YP2-YP3)
     2-A**2/360.*(-2.*YM2+25.*YM1-245.*Y0+245.*YP1-25.*YP2+2.*YP3)
     3-A**3/48.*(YM2-7.*YM1+6.*Y0+6.*YP1-7.*YP2+YP3)
     4-A**4/144.*(YM2-11.*YM1+28.*Y0-28.*YP1+11.*YP2-YP3)
     5-A**5/240.*(-YM2+3.*YM1-2.*Y0-2.*YP1+3.*YP2-YP3)
     6-A**6/720.*(-YM2+5.*YM1-10.*Y0+10.*YP1-5.*YP2+YP3)
      PP(XI)=AMAX1(0.,XI)
      PN(XI)=AMIN1(0.,XI)
C
CC CHECK COSISTENCY OF THE DATA:
      IF(NX.NE.N1.OR.NZ.NE.N2) THEN
      PRINT 777
 777  FORMAT(2X,'!!! CALLS TO INTER2 WITH NON-MATCHING DIMENSIONS.'
     1  ,' STOP.')
      STOP
      ENDIF
CC
      DO 1 K=1,NN
      IG0(K)=NINT(XD1(K))
    1 JG0(K)=NINT(XD2(K))

C  GRID EXTENSION FOR BC REMOVAL 
      DO 508 I=1,N1      
      DO 509 J=1,N2
      II=(J-1)*N1+I
  509 X(I,J)=XF(II) 
      DO 5091 IS=1-IOR,0
C     II=(1-1)*N1+I
 5091 X(I,IS)=XF(I)
      DO 5092 IS=1,IOR
      II=(N2-1)*N1+I
 5092 X(I,N2+IS)=XF(II)
  508 CONTINUE
      DO 507 J=-IOR+1,N2+IOR
      DO 5071 IS=-IOR+1,1
 5071 X(IS,J)=X(1,J)*(1-IBC)+IBC*X(N1+IS-1,J)
      DO 5072 IS=0,IOR
 5072 X(N1+IS,J)=X(N1,J)*(1-IBC)+IBC*X(1+IS,J)
  507 CONTINUE
C  END OF GRID EXTENSION
C                     
C
C  HERE STARTS REZIDUAL ADVECTION
C
                     DO 50 J=-IOR,IOR
C
      IF(LINER.EQ.1) THEN
      DO 211 II=1,NN
      U=IG0(II)-XD1(II)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      FL0=DONOR(YM1, Y0,U)
      FL1=DONOR(Y0 ,YP1,U)
  211 Z(II,J)=Y0-(FL1-FL0) 
      GO TO 50
      ENDIF
C
      IF(IOR.EQ.1) THEN
        IF(NONOS.EQ.1) THEN
      DO 311 II=1,NN
      U=IG0(II)-XD1(II)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      F0=TR2(YM1, Y0,U)
      F1=TR2(Y0 ,YP1,U)
      FL0=DONOR(YM1, Y0,U)
      FL1=DONOR(Y0 ,YP1,U)
      W=Y0-(FL1-FL0) 
      MX=AMAX1(YM1,Y0,YP1,W)
      MN=AMIN1(YM1,Y0,YP1,W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
  311 Z(II,J)=W-(F1-F0) 
        ELSE
      DO 321 II=1,NN
      U=IG0(II)-XD1(II)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      F0=TR2(YM1, Y0,U)
      F1=TR2(Y0 ,YP1,U)
  321 Z(II,J)=Y0-(F1-F0) 
        ENDIF
      ENDIF
C
      IF(IOR.EQ.2) THEN
        IF(NONOS.EQ.1) THEN
      DO 312 II=1,NN
      U=IG0(II)-XD1(II)
      YM2=X(IG0(II)-2,JG0(II)+J)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      YP2=X(IG0(II)+2,JG0(II)+J)
      F0=TR4(YM2,YM1,Y0 ,YP1,U)
      F1=TR4(YM1,Y0 ,YP1,YP2,U)
      FL0=DONOR(YM1, Y0,U)
      FL1=DONOR(Y0 ,YP1,U)
      W=Y0-(FL1-FL0) 
      MX=AMAX1(YM1,Y0,YP1,W)
      MN=AMIN1(YM1,Y0,YP1,W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
  312 Z(II,J)=W-(F1-F0) 
        ELSE
      DO 322 II=1,NN
      U=IG0(II)-XD1(II)
      YM2=X(IG0(II)-2,JG0(II)+J)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      YP2=X(IG0(II)+2,JG0(II)+J)
      F0=TR4(YM2,YM1,Y0 ,YP1,U)
      F1=TR4(YM1,Y0 ,YP1,YP2,U)
  322 Z(II,J)=Y0-(F1-F0) 
        ENDIF
      ENDIF
C
      IF(IOR.EQ.3) THEN
        IF(NONOS.EQ.1) THEN
      DO 313 II=1,NN
      U=IG0(II)-XD1(II)
      YM3=X(IG0(II)-3,JG0(II)+J)
      YM2=X(IG0(II)-2,JG0(II)+J)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      YP2=X(IG0(II)+2,JG0(II)+J)
      YP3=X(IG0(II)+2,JG0(II)+J)
      F0=TR6(YM3,YM2,YM1,Y0 ,YP1,YP2,U)
      F1=TR6(YM2,YM1,Y0 ,YP1,YP2,YP3,U)
      FL0=DONOR(YM1, Y0,U)
      FL1=DONOR(Y0 ,YP1,U)
      W=Y0-(FL1-FL0) 
      MX=AMAX1(YM1,Y0,YP1,W)
      MN=AMIN1(YM1,Y0,YP1,W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
  313 Z(II,J)=W-(F1-F0) 
        ELSE
      DO 323 II=1,NN
      U=IG0(II)-XD1(II)
      YM3=X(IG0(II)-3,JG0(II)+J)
      YM2=X(IG0(II)-2,JG0(II)+J)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      YP2=X(IG0(II)+2,JG0(II)+J)
      YP3=X(IG0(II)+2,JG0(II)+J)
      F0=TR6(YM3,YM2,YM1,Y0 ,YP1,YP2,U)
      F1=TR6(YM2,YM1,Y0 ,YP1,YP2,YP3,U)
  323 Z(II,J)=Y0-(F1-F0) 
        ENDIF
      ENDIF
C
C
   50 CONTINUE
C  
      IF(LINER.EQ.1) THEN
      DO 212 II=1,NN
      U=JG0(II)-XD2(II)
      FL0=DONOR(Z(II,-1),Z(II,0),U)
      FL1=DONOR(Z(II, 0),Z(II,1),U)
  212 XF(II)=Z(II,0)-(FL1-FL0) 
      RETURN
      ENDIF
C
      IF(IOR.EQ.1) THEN
        IF(NONOS.EQ.1) THEN
      DO 411 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR2(Z(II,-1),Z(II,0),U)
      F1=TR2(Z(II, 0),Z(II,1),U)
      FL0=DONOR(Z(II,-1),Z(II,0),U)
      FL1=DONOR(Z(II, 0),Z(II,1),U)
      W=Z(II,0)-(FL1-FL0) 
      MX=AMAX1(Z(II,-1),Z(II,0),Z(II,1),W)
      MN=AMIN1(Z(II,-1),Z(II,0),Z(II,1),W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
      XF(II)=W-(F1-F0) 
  411 CONTINUE
        ELSE
      DO 421 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR2(Z(II,-1),Z(II,0),U)
      F1=TR2(Z(II, 0),Z(II,1),U)
  421 XF(II)=Z(II,0)-(F1-F0) 
        ENDIF
      ENDIF

      IF(IOR.EQ.2) THEN
        IF(NONOS.EQ.1) THEN
      DO 412 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR4(Z(II,-2),Z(II,-1),Z(II,0),Z(II,1),U)
      F1=TR4(Z(II,-1),Z(II, 0),Z(II,1),Z(II,2),U)
      FL0=DONOR(Z(II,-1),Z(II,0),U)
      FL1=DONOR(Z(II, 0),Z(II,1),U)
      W=Z(II,0)-(FL1-FL0) 
      MX=AMAX1(Z(II,-1),Z(II,0),Z(II,1),W)
      MN=AMIN1(Z(II,-1),Z(II,0),Z(II,1),W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
      XF(II)=W-(F1-F0) 
  412 CONTINUE
        ELSE
      DO 422 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR4(Z(II,-2),Z(II,-1),Z(II,0),Z(II,1),U)
      F1=TR4(Z(II,-1),Z(II, 0),Z(II,1),Z(II,2),U)
  422 XF(II)=Z(II,0)-(F1-F0) 
        ENDIF
      ENDIF

      IF(IOR.EQ.3) THEN
        IF(NONOS.EQ.1) THEN
      DO 413 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR6(Z(II,-3),Z(II,-2),Z(II,-1),Z(II,0),
     1                     Z(II, 1),Z(II, 2),U)
      F1=TR6(Z(II,-2),Z(II,-1),Z(II, 0),Z(II,1),
     1                     Z(II, 2),Z(II, 3),U)
      FL0=DONOR(Z(II,-1),Z(II,0),U)
      FL1=DONOR(Z(II, 0),Z(II,1),U)
      W=Z(II,0)-(FL1-FL0) 
      MX=AMAX1(Z(II,-1),Z(II,0),Z(II,1),W)
      MN=AMIN1(Z(II,-1),Z(II,0),Z(II,1),W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
      XF(II)=W-(F1-F0) 
  413 CONTINUE
        ELSE
      DO 423 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR6(Z(II,-3),Z(II,-2),Z(II,-1),Z(II,0),
     1                     Z(II, 1),Z(II, 2),U)
      F1=TR6(Z(II,-2),Z(II,-1),Z(II, 0),Z(II,1),
     1                     Z(II, 2),Z(II, 3),U)
  423 XF(II)=Z(II,0)-(F1-F0) 
        ENDIF
      ENDIF
      RETURN
      END  