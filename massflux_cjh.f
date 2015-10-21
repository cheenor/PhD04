      PARAMETER (im=202, km=52, inn=480, itt=2880, ii=121, ifi=6)
      dimension qc(im,km),qr(im,km),qa(im,km),qb(im,km)
      dimension xi(km+1)
      dimension w(im,km),rho(km),uf(im,km),df(im,km)
      dimension smuf(itt,km),smdf(itt,km)
      dimension tmuf(ii,km),tmdf(ii,km)

      character*100 chenm,path
      character casenm*20, fold*20
      casenm='casename'
      fold='runfold'
      if(casenm(1:3)=='MLY')then
        path='/nuist/scratch/chenjh/ERA_Interim/'//casenm(1:4)//
     +   '/'//trim(fold)
      else
        path='/nuist/scratch/chenjh/ERA_Interim/'//casenm(1:3)//
     +   '/'//trim(fold)
      endif
      chenm=trim(path)//'/'//trim(casenm)
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
!----------- output ---------------------------------------
      open(50,file=trim(chenm)//'_massflux_up.txt')
      open(51,file=trim(chenm)//'_massflux_down.txt')
      open(52,file=trim(chenm)//'_meanmassflux_updown.txt')
!
      rat=15.
      XI(2)=0.
      DO 152 K=2,km-1
      RATZ=RAT
      DEL=100.
      nzm1=km-1
      k1=k
      XI(K+1)=XI(K)+((RATZ-1.)/FLOAT(NZM1-2)*FLOAT(K1-2)+1.)*DEL
  152 CONTINUE
      xi(1)=2.*xi(2)-xi(3)
      xi(km+1)=2.*xi(km)-xi(km-1)

      it=0
      do 999 if=1,ifi
      IH=80+if
      do 100 in=1,inn
      it=it+1
      read(IH) qc,qr 
      read(IH) qa,qb 
      read(IH) 
      read(IH) 
      read(IH) w,w
      read(IH) 
      read(IH) rho
      read(IH)
      do  k=1,km
      do  j=1,im
      uf(j,k)=0.
      df(j,k)=0.
      enddo
      enddo
      do 2103 k=2,km-1
      do 2103 j=1,im-1
      upmass=0.5*(0.5*(rho(k-1)+rho(k))*w(j,k-1)
     *       +0.5*(rho(k)+rho(k+1))*w(j,k))*10.
      tc=(qc(j,k)+qr(j,k)+qa(j,k)+qb(j,k))*1000.
      if(tc.ge.0.1) then
      uf(j,k)=uf(j,k)+cvmgm(0.,upmass,upmass)*9.8/100.*3600.
      df(j,k)=df(j,k)+cvmgm(upmass,0.,upmass)*9.8/100.*3600.
      endif
c       cwp=rho(k)*(qc(j,k)+qa(j,k))*(xi(k+1)-xi(k))*1.e3  ! g/m**3
c       if(cwp.gt.0.2) then  
c       uf(j,k)=uf(j,k)+cvmgm(0.,upmass,upmass)*9.8/100.*3600.  !mb/hour
c       df(j,k)=df(j,k)+cvmgm(upmass,0.,upmass)*9.8/100.*3600.
c       endif
 2103 continue

      do 2103 j=1,im-1
        write(50,99)(uf(j,k),k=1,km)
        write(51,99)(df(j,k),k=1,km)
      enddo
      do 40 k=1,km
      smuf(it,k)=0.
      smdf(it,k)=0.
      do 30 i=2,im-1
      smuf(it,k)=smuf(it,k)+uf(i,k)/float(im-2)
      smdf(it,k)=smdf(it,k)+df(i,k)/float(im-2)
  30  continue
  40  continue
 100  continue
 999  continue
      do 200 k=1,km
      do 50 mt=1,ii
      n1=(mt-1)*24-11
      n2=(mt-1)*24+12
      if(mt.eq.1) then
      n1=1
      n2=12
      endif
      if(mt.eq.ii) then
      n1=itt-11
      n2=itt
      endif
      tmuf(mt,k)=0.
      tmdf(mt,k)=0.
      do 11 n=n1,n2
      tmuf(mt,k)=tmuf(mt,k)+smuf(n,k)/float(n2-n1+1)
      tmdf(mt,k)=tmdf(mt,k)+smdf(n,k)/float(n2-n1+1)
  11  continue
  50  continue
 200  continue
      do mt=1,ii
      write(52,98) (tmuf(mt,k),k=1,km),(tmdf(mt,k),k=1,km)
      enddo
99    format(1X,52(1X,E12.4))
98    format(8e12,4)
c  units:
c  updraft mass flux (tmuf   mb/hr)
c  downdraft mass flux (tmdf   mb/hr)
      stop
      end
       REAL FUNCTION CVMGM(A,B,C)
       REAL A,B,C
       IF(C.GE.0.0) THEN
       CVMGM=B
       ELSE
       CVMGM=A
       ENDIF
       RETURN
       END 

