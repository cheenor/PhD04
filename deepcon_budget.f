      PARAMETER (im=202, km=52, inn=480, itt=2880, ii=121, ifi=6)
      PARAMETER (mx=im,nxg=im,nz=km-1,nzg=33,nx=im)
      dimension XI(km)
      dimension x(nx),z(nz),zm(nz),xg(nxg),zg(nzg)
      dimension work1(nx,nz),work2(nx,nz)
c      dimension qc(im,km),qr(im,km),qa(im,km),qb(im,km),rho(km)
      dimension q1e(im,km),q2e(im,km),fx(im,km)
      dimension q1d(im,km),q2d(im,km)
      dimension q1s(im,km),q2s(im,km),rlw(im,km),rsw(im,km)
      dimension e1flux(im,km),e2flux(im,km),eflux(im,km)
      dimension sub(im,km),con(im,km),eva(im,km)
      dimension fus(im,km),dep(im,km)
!
      dimension qc(im,km),qr(im,km),qa(im,km),qb(im,km) 
      dimension th(im,km),q(im,km),rh(im,km),w(im,km)
      dimension te(km),ps(km),potf(km),pote(km)
      dimension qrl(im,km),qrs(im,km),u(im,km)
      dimension tq(im,km)
!
      dimension condc(5,km),q1q2com(12,km),envdc(3,km)
      dimension qlp(km),qip(km),qrldc(km),qrsdc(km)
      integer KB(99),KE(99),NA
c      dimension u(im,km),w(im,km),t(im,km),q(im,km),rh(im,km)

      character*100 chenm, path ! case name ,must be recorded 
      character casenm*20,fold*20
      casenm='ETPCTR_EC'
      fold='runctr2'
      if(casenm(1:3)=='MLY')then
        path='/home/jhchen/jhchen/ERA_Interim/'//casenm(1:4)//
     +   '/'//trim(fold)
      else
        path='/home/jhchen/jhchen/ERA_Interim/'//casenm(1:3)//
     +   '/'//trim(fold)
      endif
!      path='/home/jhchen/jhchen/ERA_Interim/ETP/'//trim(fold)
      chenm=trim(path)//'/eddydiffradcon_'//trim(casenm)    !'whereisthe_data_tgn2d1_(1-6)'      
       open(81,file=trim(chenm)//'_1'
     * ,form='UNFORMATTED',status='OLD',convert='big_endian') 
      open(82,file=trim(chenm)//'_2'
     *,form='unformatted',status='old',convert='big_endian') 
      open(83,file=trim(chenm)//'_3'
     *,form='unformatted',status='old',convert='big_endian') 
      open(84,file=trim(chenm)//'_4'
     *,form='unformatted',status='old',convert='big_endian') 
      open(85,file=trim(chenm)//'_5'
     *,form='unformatted',status='old',convert='big_endian') 
      open(86,file=trim(chenm)//'_6'
     *,form='unformatted',status='old',convert='big_endian') 
!
      chenm=trim(path)//'/'//trim(casenm)    !'whereisthe_data_tgn2d1_(1-6)'      
       open(91,file=trim(chenm)//'_1'
     * ,form='UNFORMATTED',status='OLD',convert='big_endian') 
      open(92,file=trim(chenm)//'_2'
     *,form='unformatted',status='old',convert='big_endian') 
      open(93,file=trim(chenm)//'_3'
     *,form='unformatted',status='old',convert='big_endian') 
      open(94,file=trim(chenm)//'_4'
     *,form='unformatted',status='old',convert='big_endian') 
      open(95,file=trim(chenm)//'_5'
     *,form='unformatted',status='old',convert='big_endian') 
      open(96,file=trim(chenm)//'_6'
     *,form='unformatted',status='old',convert='big_endian') 
c      open(51,file=trim(chenm)//'_rho.txt')
      open(52,file=trim(chenm)//'_DeepCons_budgetandprofs.txt')
c      open(53,file=trim(chenm)//'_Raw_q12seflux.txt')
c      open(54,file=trim(chenm)//'_Raw_rlwrsw.txt')
      rat=15.
      XI(2)=0.
      DO 162 K=2,NZM
      RATZ=RAT
      DEL=100.
      nzm1=nz
      k1=k
      XI(K+1)=XI(K)+((RATZ-1.)/FLOAT(NZM1-2)*FLOAT(K1-2)+1.)*DEL
 162  CONTINUE
      do 153 k=1,nz
 153  z(k)=XI(K+1)/1000.
      zm(1)=0.
      do 154 k=2,nz
 154  zm(k)=0.5*(z(k-1)+z(k))
      do 150 i=1,nx
 150  x(i)=float(i-1)
      do 252 k=1,nzg
 252  zg(k)=float(k-1)/2.
      do 250 i=1,nxg
 250  xg(i)=float(i-1)
!-----------------------------------------------------------
      qlp=0.
      qip=0.
      qrsdc=0
      qrldc=0.
      condc=0.
      q1q2com=0.
      icount=0

      sec=10.           !!!Time Step
      it=0
      do 999 if=1,ifi
      IH=80+if
      IH2=90+IF
      do 100 in=1,inn
      it=it+1
      read(IH) q1e,q2e,fx 
      read(IH) q1d,q2d,q1s,q2s
      read(IH) rlw,rsw
      read(IH) e1flux,e2flux,eflux
      read(IH) con
      read(IH) eva
      read(IH) dep
      read(IH) sub
      read(IH) fus
!      
      read(IH2) qc,qr 
      read(IH2) qa,qb 
      read(IH2) th,q
      read(IH2) rh,te,ps
      read(IH2) u,w
      read(IH2) 
      read(IH2) rho,potf,pote
      read(IH2) qrl, qrs
C
CC  the following loop output the data of every time stemp and grid. 
c      do i=2,im-1
c      write(52,99)q1e(i,:),q2e(i,:),fx(i,:)
c     +            ,q1d(i,:),q2d(i,:)
c      write(53,99)q1s(i,:),q2s(i,:)
c     +            ,e1flux(i,:),e2flux(i,:)
c     +            ,eflux(i,:)
c      write(54,99)rlw(i,:),rsw(i,:) 
c      end do
C
     
c      write(51,*)k,rho(k)
C      call gridint(q1e,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(q2e,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(fx,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(q1d,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(q2d,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(q1s,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(q2s,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(rlw,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(rsw,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(e1flux,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(con,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(eva,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(dep,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(sub,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(fus,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c !          
c      call gridint(qc,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(qr,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(qa,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(qb,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(th,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(q,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(rh,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(u,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(w,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(qrl,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
c      call gridint(qrs,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
!
      do 30 i=2,im-1
      do 40 k=1,km !nzg
        tq(i,k)=qc(i,k)*1000.+ qr(i,k)*1000.
     &              + qa(i,k)*1000.+qb(i,k)*1000.           
  40  continue
      call DEEPCC(tq,im,km,kb,ke,na,i)
      IDC=0
      if(na>0)print*,'NA = ',na
      DO L = 1,NA
          IF (KB(L).LE.5 .AND. KE(L).GE.24) THEN  !! 5 AND 24 ARE THE VERTICAL LEVELS
          IDC=IDC+1
          ENDIF
      ENDDO
      if (idc>0) then ! there is a deep convection in the columm
        do kk=1,km          
          qlp(kk)=qlp(kk)+qc(i,kk)*1000.+ qr(i,kk)*1000.
          qip(kk)=qip(kk)+qa(i,kk)*1000.+qb(i,kk)*1000.
          qrldc(kk)=qrldc(kk)+qrl(i,kk)
          qrsdc(kk)=qrsdc(kk)+qrs(i,kk)
          tmptc=th(i,kk)*te(kk)/potf(kk)/(1.+pote(kk))-273.16 
          envdc(1,kk)=envdc(1,kk)+tmptc
          envdc(2,kk)=envdc(2,kk)+q(i,kk)
          envdc(3,kk)=envdc(3,kk)+w(i,kk)
          !con    eva   dep   sub    fus      
          condc(1,kk)=condc(1,kk)+con(i,kk)
          condc(2,kk)=condc(2,kk)+eva(i,kk)
          condc(3,kk)=condc(3,kk)+dep(i,kk)
          condc(4,kk)=condc(4,kk)+sub(i,kk)
          condc(5,kk)=condc(5,kk)+fus(i,kk)
          !   q1e,q2e,fx, q1d,q2d,q1s,q2s, rlw,rsw,e1flux,e2flux,eflux
          q1q2com(1,kk)=q1q2com(1,kk)+q1e(i,kk)
          q1q2com(2,kk)=q1q2com(2,kk)+q1d(i,kk)
          q1q2com(3,kk)=q1q2com(3,kk)+q1s(i,kk)
          q1q2com(4,kk)=q1q2com(4,kk)+fx(i,kk)
          q1q2com(5,kk)=q1q2com(5,kk)+rlw(i,kk)
          q1q2com(6,kk)=q1q2com(6,kk)+rsw(i,kk)
          q1q2com(7,kk)=q1q2com(7,kk)+e1flux(i,kk)
          q1q2com(8,kk)=q1q2com(8,kk)+q2e(i,kk)
          q1q2com(9,kk)=q1q2com(9,kk)+q2d(i,kk)
          q1q2com(10,kk)=q1q2com(10,kk)+q2s(i,kk)
          q1q2com(11,kk)=q1q2com(11,kk)+e2flux(i,kk)
          q1q2com(12,kk)=q1q2com(12,kk)+eflux(i,kk)
        enddo
        icount=icount+1
      endif
  30  continue

 100  continue   
 999  continue   
      if (icount>0)then
        do 200 kk=1,km
          qlp(kk)=qlp(kk)/icount
          qip(kk)=qip(kk)/icount
          qrldc(kk)=qrldc(kk)/icount
          qrsdc(kk)=qrsdc(kk)/icount
          do ik =1,3 
            envdc(ik,kk)=envdc(ik,kk)/icount
          enddo
          do ik =1,5
            condc(ik,kk)=condc(ik,kk)/icount
          enddo
          do ik =1,12
            q1q2com(ik,kk)=q1q2com(ik,kk)/icount
          enddo
 200    continue
      endif
      write(52,99)(qlp(kk),kk=1,km)
      write(52,99)(qip(kk),kk=1,km)
      write(52,99)(qrldc(kk),kk=1,km)
      write(52,99)(qrsdc(kk),kk=1,km)
      do ik=1,3
        write(52,99)(envdc(ik,kk),kk=1,km)
      enddo
      do ik=1,5
        write(52,99)(condc(ik,kk),kk=1,km)
      enddo
      do ik=1,12
        write(52,99)(q1q2com(ik,kk),kk=1,km)
      enddo
99    format(1X,52(1X,E12.4))
c  units:              
c     all the variables are in K/day
      stop
      end
C     ------------------------------
      SUBROUTINE DEEPCC(QW,IM,KM,KB,KE,NA,imm) ! QW TOTAL CLOUD WATER,IM: HORIZONTAL GRID; KM : VERTICAL GRID
C     ------------------------------
C     ------------------------------
C
C     OUTPUT DEEP CONVECTION
C
      INTEGER IM,KM,OU,I,K,L,KB(99),KE(99),NA
      REAL QW(IM,KM),C(99),CM(99)
      INTEGER IDC,imm

!      WRITE(OU,'(A3,1X,I4,1X,\)')'ITT', IT !  \ NOT CHANGE LINE
!      DO I = 1,IM
        DO K = 1,KM
          IF (QW(imm,K).GE.1.0E-3) C(K) = 1. ! THIS GRIG IS COVERED BY CLOUD
          IF (QW(imm,K).LT.1.0E-3) C(K) = 0.
        ENDDO
        CALL INFCLD(C,KM,KB,KE,CM,NA)      
!      ENDDO
      RETURN
      END
C     -----------------------------------
C     -----------------------------------
      SUBROUTINE INFCLD(C,NL,KB,KE,CM,NA) !KB BASE; KE TOP, NA: LAYERS OF CLOUDS
C     -----------------------------------
C     -----------------------------------
C
C     FIND INFORMATION ABOUT ADJACENT CLOUD LAYERS
C
      IMPLICIT NONE
C
      INTEGER NL,KB(*),KE(*),NA,K1,K2,K
      REAL C(NL),CM(*),AA
C
      NA = 0
      IF (C(1).GT.0.) THEN
          K1 = 1
          K2 = 1
          AA = C(1)
      ELSE
          AA = 0.
      ENDIF
      DO K = 2,NL
         IF (C(K-1).LE.0..AND.C(K).GT.0.) THEN
             K1 = K
             K2 = K
             AA = MAX(AA,C(K))
         ELSEIF (C(K-1).GT.0..AND.C(K).GT.0.) THEN
             K2 = K
             AA = MAX(AA,C(K))
         ELSEIF (C(K-1).GT.0..AND.C(K).LE.0.) THEN
             NA = NA + 1
             KB(NA) = K1
             KE(NA) = K2
             CM(NA) = AA
             AA = 0.
         ENDIF
      ENDDO
      IF (C(NL).GT.0.) THEN  !TOP
          NA = NA + 1
          KB(NA) = K1
          KE(NA) = K2
          CM(NA) = AA
      ENDIF
      RETURN
      END
!
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine gridint
     1 (ain,nx,nz,xx,zz,nx1,nz1,xx1,zz1,x0,z0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc this subroutine performs interpolation of the data given in the
cc input array ain(nx,nz) on the grid defined by xx(nx) and zz(nz)
cc into the grid given by xx1(nx1) and zz1(nz1). NIETHER OF THE
cc GRIDS HAS TO BE REGULAR. Data is returned in the ain(nx1,nz1)
cc part of the input array.
cc 
cc    levels in the input array are given in zz(nz), 
cc    levels in the output array are given in zz1(nz1)
cc      x-coordinate in the input array are in xx(nx)
cc      x-coordinate in the output array are in xx1(nx1)
cc        x0(nx,nz) and z0(nx,nz) are working arrays
cc
cc NOTE that nx1 (or nz1) must be smaller than nx (or nz) and xx1 (zz1)
cc  must be a subdomain of xx (zz)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension ain(nx,nz),zz(nz),xx(nx)
      dimension zz1(nz1),xx1(nx1) 
      dimension x0(nx,nz),z0(nx,nz)

c      check consistency of the input data
      ier=0
cc nx1,nz1 not larger than nx,nz
      if(nx1.gt.nx) ier=1
      if(nz1.gt.nz) ier=2
cc limits (zz1(nz1).le.zz(nz) ?) 
      if(zz1(1).lt.zz(1)) ier=3
      if(zz1(nz1).gt.zz(nz)) ier=4
cc limits (xx1(nx1).le.xx(nx) ?) 
      if(xx1(1).lt.xx(1)) ier=5
      if(xx1(nx1).gt.xx(nx)) ier=6
      if(ier.ne.0) then
      print 999,ier
 999  format(2x,' ** problems with input data. will stop.'/
     1 ' ier = ',i3,'. see code why stoped.')
      stop
      endif
cc
      nxz=nx*nz
      do 99 i=1,nxz
      z0(i,1)=1.
  99  x0(i,1)=1.
cc  map vertical grid positions:
      do 1 k1=1,nz1
      zzh=zz1(k1)
      do 2 k=1,nz
      kk=k
      if(zz(k).ge.zzh) go to 6
  2   continue
  6   kkm=max0(1,kk-1)
      z0(1,k1)=float(kkm)+(zzh-zz(kkm))/(zz(kk)-zz(kkm)+1.e-6)
  1   continue
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
c      print*,x0(:,1)
c      print*,z0(1,:)
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
      PARAMETER (N1=121, N2=51,NN=N1*N2)  !!!!!!! N1 is im in, N2 is km-1
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
