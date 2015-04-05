      PROGRAM ERA_interim
!------------this program reads the EAR_interim dataset download from
!                      http://apps.ecmwf.int/datasets/    (includeing normal reanalyses and surface variables)
!--- run in the WinXP system and the netcdf library is needed for the fortran compiler to decode the netcdf files
!--- you can change the data start date and end date as you need
!-- NOTE: when more than year data is handled, the momery maybe is unsunfficiency.
!----------------------------------------------------------------------------------
!------------ 3-D variables are u v omega sh Temp. and Hgt with 17 levels--------
!------the surface variables include latent heat/sensible heat  
!  t2m td2m u10m,v10m ps
!-------------------------------last modified 26/2/2015  Jinghua Chen
      include 'netcdf.inc'
      integer,parameter :: nrec=1460
      integer,parameter :: nrecr=1464
      integer,parameter :: nt=5   ! it-2,it-1,it,it+1,it+1
      integer,parameter :: nz=37 !! add surface
      integer,parameter :: n3d=6 !1u 1v 3t  4oemga 5RH 6HGT
      integer,parameter :: ndx=10
      integer,parameter :: ndy=10
      integer,parameter :: nrg=2
      integer,parameter :: nday=366
      COMMON/D3D/ air(ndx,ndy,nz,nt),hgt(ndx,ndy,nz,nt),
     +            qv(ndx,ndy,nz,nt),uwnd(ndx,ndy,nz,nt), !1u 1v 3t  4oemga 5RH 6HGT
     +            vwnd(ndx,ndy,nz,nt),omega(ndx,ndy,nz,nt), ! !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top 11cpr 12DLR 13ULR,14t2m
     +            theta(ndx,ndy,nz,nt)
      COMMON/SRF/ t2m(ndx,ndy,nt),td2m(ndx,ndy,nt),
     +            u10m(ndx,ndy,nt),v10m(ndx,ndy,nt),
     +            psf(ndx,ndy,nt),qvs(ndx,ndy,nt),
     +            Geo(ndx,ndy)
    
!      COMMON/TQ/ tvs(ndx,ndy,nt),qvs(ndx,ndy,nt),
!     +      tv(ndx,ndy,nz,nt),qv(ndx,ndy,nz,nt)
      COMMON/Q12/ Q1(ndx,ndy,NZ),Q2(ndx,ndy,NZ),
     +  tls(ndx,ndy,NZ),qls(ndx,ndy,NZ)
      COMMON/VQ12/VQ1(ndx,ndy),VQ2(ndx,ndy)
!      COMMON/OMGOUT/OMGOO(ndx,ndy),DIVP(ndx,ndy)
      COMMON/AVEG/HAD_Q(ndx,ndy,NZ),VAD_Q(ndx,ndy,NZ),
     +  TCH_Q(ndx,ndy,NZ),HADT(ndx,ndy,NZ),VADT(ndx,ndy,NZ),
     +  TCHT(ndx,ndy,NZ) 
      REAL DAT3D(ndx,ndy,nz,nday,4,6)
!!---------- output data ----------------------------------------------------------------------
      COMMON/DXY/ DX(NDY),DY
!!!!!!!!!!!!
      real OMGS(ndx,ndy)

!!!!!!!!!!!!!!!!!!!!!!
      character(len=150) :: dirin,fpath,fold,filename,dirout
      character(len=30) :: name3d(n3d),name2d(5),area(nrg),foldout,date
!      
      character(len=4) :: years,yeare,year
      character(len=2) months,monthe,dys,dye,
     +               hour(4),hrstr
      integer IXX(nrg),IYY(nrg),IXG(nrg),IYG(nrg),ntec
      integer X1(nrg),Y1(nrg),X2(nrg),Y2(nrg),  
     +        XG1(nrg),YG1(nrg),XG2(nrg),YG2(nrg) 
      real plv(nz),dp(nz), bslatG(nrg),bslat(nrg),bslon(nrg)
      integer iys,iye,ims,ime,ids,ide,days(12)
!      integer ktop(ndx,ndy) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----variables with the suffix of '_xy' are average of the grid box and for output 
      real tls_xy(nz),
     +     qls_xy(nz),
     +     u_xy(nz),
     +     the_xy(nz),
     +     v_xy(nz),
     +     rh_xy(nz),
     +     gz_xy(nz),
     +     Q1_xy(nz),
     +     air_xy(nz),
     +     Q2_xy(nz),
     +     qv_xy(nz),
     +     omg_xy(nz),
     +     omgo_xy(nz),
     +     psrf_xy,
     +     t2m_xy,
     +     u10m_xy, 
     +     v10m_xy, 
     +     qv2m_xy,
     +     omgco_xy(nz),
     +     DIV_xy(nz)
!!!-----------------------------------------------------                
      real DIVO(ndx,ndy,nz),OMGO(ndx,ndy,nz)
     *          ,DIV(ndx,ndy,nz),OMG(ndx,ndy,nz)
      real XY_OUT(18),XYN_OUT(nz,18)
!----------------------------------------------------------------------------------------------
      real omegas(ndx,ndy,nt) !WSFC(ndx,ndy,nt),
      real qlss(ndx,ndy,nt),tlss(ndx,ndy,nt)
      real Q1temps(ndx,ndy,nt),Q2temps(ndx,ndy,nt)
      real TS_XY(nt),US_XY(nt),VS_XY(nt),GZS_XY,QVS_XY(nt)
      real omegaS_XY(nt),Q1S_XY(nt),Q2S_XY(nt),dyn(nz,nt,4) ! 1-3 for vorticity 4 for divergence
      real TLSS_XY(nt),QLSS_XY(nt), rhs_xy(nt),thes_xy(nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real tmp3d(ndx,ndy,nz,nday),tmp2d(5,ndx,ndy,4,366)
     +	,tmp2d0(ndx,ndy,366)
!
      DIMENSION plw(nz)
!
      DATA rd,cp,g,rv,hlat /287.,1005.,9.81,461.,2.5e6/
      DATA plv/1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.
     +         ,750.,700.,650.,600.,550.,500.,450.,400.,350.,300.
     +         ,250.,225.,200.,175.,150.,125.,100.,70.,50.,30.,20.
     +         ,10.,7.,5.,3.,2.,1./
!      data plw/962.5, 887.5, 775., 650., 550., 450., 350., 275., 225., 
!     +      175., 125., 85., 60., 40., 25., 15., 10./   !plw(K)=(plv(k)+plv(k+1))/2 plw(17)=plv(16)/2
!      data dp/75.,75.,150.,100.,100.,100.,100.,100.,50.,
!     +        50.,50.,50.,30.,20.,20.,10.,10./   ! top pressure is zero donw-up
      plw=0.0
      do iz=1,nz-1
        plw(iz)=plv(iz)/2.+plv(iz+1)/2.
        dp(iz)=plv(iz)-plv(iz+1)
      enddo
      plw(nz)=plv(nz)
      dp(nz)=plv(nz)-0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------setting ---------------------------------------------------
!!!!-------------------------set ----------------------------------------------------
! for NCEP : variables that named ending with character 'G' stand for the Gauss Grid 
! all ERA dataset  on the same grid system
! the following variables are for determining different regions 
! area stores the region's name, X1(:) Y1(:)store the start  gird ID,
!  X2(:) Y2(:) is end grid ID, IXX(:) IYY(:)  the number grid of different regions. 
! bslat bslon the star lat. and lon. of different regions 
! target area 3 grid X 3 grid, but when Q1 and Q2 are calculated, the loop is:
!    2 to IXX-1 and 2 to IYY-1, so 5 grid X 5 grid is needed
!-------WTP=Weast TP   77.5-90   27.5 37.5 
      area(1)='WTP'
      X1(1)  = 32;  X2(1)  = 37;  Y1(1)= 26;  Y2(1)=22
      IXX(1) = 6 ;  IYY(1) = 5 ;  bslat(1)=27.5; bslon(1)=77.5  
!      XG1(1) = 42;  XG2(1) = 49;  YG1(1)=33 
!      YG2(1) = 27;  IYG(1) =  7;  bslatG(1)=27.6186
!      IXG(1) =  8;  
!    East TP  90-102.5
      area(2)='ETP'  
      X1(2)  = 37;  X2(2) = 42;  Y1(2) =  26;    Y2(2)   = 22
      IXX(2) =  6;  IYY(2)=  5;  bslat(2)= 27.5; bslon(2)= 90.
!      XG1(2)= 49;  XG2(2)= 56;  YG1(2)= 33;  YG2(2)=27
!      IYG(2)=  7;  IXG(2)=  8;  bslatG(2)=27.6186
!     bay of bengal 85-95 15 -25
!      area(3)='BOB'; X1(3)=31; X2(3)=35;Y1(3)=31; Y2(3)=27; IXX(3)=5
!      IYY(3)=5; bslat(3)=15.;XG1(3)=46;XG2(3)=52;YG1(3)=40;YG2(3)=34
!      IYG(3)=7; bslatG(3)=14.3;IXG(3)=7;bslon(3)=85. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      name3d(1)='hgt' 
      name3d(2)='sh' 
      name3d(3)='air' 
      name3d(4)='uwnd' 
      name3d(5)='vwnd' 
      name3d(6)='omega'  
      name2d(1)='ps' 
      name2d(2)='t2m' 
      name2d(3)='td2m' 
      name2d(4)='u10m' 
      name2d(5)='v10m'
      do i=1,12
        days(i)=31
      enddo
      days(2)=28
      days(6)=30
      days(4)=30
      days(9)=30
      days(11)=30
      hour(1)='00'
      hour(2)='06'
      hour(3)='12'
      hour(4)='18'
      ims=1      ;   ime=12           !
      idss=1      ;   ide=31
      dirout='X:\Data\ERA_interim\ERA_Pre_O\'
      dirin='X:\Data\ERA_interim\'
      open(404,file='404.txt')
      do 1111 iyear=1979,2012
	  iys=iyear
	  iye=iyear
	  print*,iyear
        write(years,'(I4.4)')iyear
        write(yeare,'(I4.4)')iyear
        write(months,'(I2.2)')ims
        write(monthe,'(I2.2)')ime
        write(dys,'(I2.2)')idss
        write(dye,'(I2.2)')ide
        foldout=years(3:4)//months//dys//'-'//yeare(3:4)//monthe//dye
        istatus1=CHDIR(trim(dirout))
        istatus2=system("Md "//trim(foldout))

        days(2)=28
        ntt=nrec
        ndd=365
        if(mod(iyear,4)==0.and.mod(iyear,100)/=0)then
          days(2)=29 ! nrec=1464
          ntt=nrecr
          ndd=366
        elseif(mod(iyear,400)==0)then
          days(2)=29
          ntt=nrecr
          ndd=366
        endif
        do 1110 ig=1,nrg
	     print*,area(ig)
!------------open output files---------------------------------------
!-----------------surface heat flux and rainfall water
!          fpath=trim(dirout)//'/'//trim(foldout)//'/'//trim(area(ig))//
!     +     '_HF.txt'
!          open(10,file=trim(fpath))
!          write(10,905)'From',iys,ims,ids,'to',iye,ime,ide,
!     +     'Every six hour one data'
!          write(10,*)'TimeID latentHeat(W/m^2) SensibleHeat(W/m^2) Precipipitable_water_for_entire_atmosphere(kg/m^2)'
!-----------------sounding ------------------------------
          fpath=trim(dirout)//'/'//trim(foldout)//'/'//trim(area(ig))//
     +      '_sounding.txt'
          open(11,file=trim(fpath))
          write(11,905)'From',iys,ims,ids,'to',iye,ime,ide,
     +          'Every six hour one data'
          write(11,'(33A,1X,44A)')'Press(hPa) heigh(m) U(m/s) V(m/s)'
     +               ,  'omega(pa/s) Temp(K) Theta(K) Qv(kg/kg) RH(%)'  
!-----------large scale forcing -------------------
          fpath=trim(dirout)//'/'//trim(foldout)//'/'//trim(area(ig))//
     +         '_Forcing.txt'
          open(12,file=trim(fpath))
          write(12,905)'From',iys,ims,ids,'to',iye,ime,ide,
     +      'Every six hour one data'
          write(12,'(33A,1X,30A)')'TimeID 37_levels_T_forcing(K/day)'
     +             ,  '37_levels_qv_forcing(g/kg/day)'
!---------------Q1 and Q2---------------------------------------------    
          fpath=trim(dirout)//'/'//trim(foldout)//'/'//trim(area(ig))//
     +      '_q1q2.txt'
          open(13,file=trim(fpath))
          write(13,905)'From',iys,ims,ids,'to',iye,ime,ide,
     +     'Every six hour one data'
          write(13,*)'TimeID Q1(K/day)  Q2(K/day)'
!---------------------------------------------------------------------
          fpath=trim(dirout)//'/'//trim(foldout)//'/'//trim(area(ig))//
     +      '_SURFACE22.txt'
          open(14,file=trim(fpath))
          write(14,905)'From',iys,ims,ids,'to',iye,ime,ide,
     +        'Every six hour one data'
          write(14,'(36A,1X,38A)')'ID Press(hPa) heigh(m) U(m/s) V(m/s)'
     +           ,  'omega(pa/s) Temp(K) Theta(K) Qv(kg/kg) VQ1 VQ2'
!
          fpath=trim(dirout)//'/'//trim(foldout)//'/'//years//
     +          trim(area(ig))//'_lowlevel.txt'
          open(15,file=trim(fpath))
          fpath=trim(dirout)//'/'//trim(foldout)//'/'//years//
     +          trim(area(ig))//'_Middlelevel.txt'
          open(16,file=trim(fpath))
          fpath=trim(dirout)//'/'//trim(foldout)//'/'//years//
     +          trim(area(ig))//'_Hightlevel.txt'
          open(17,file=trim(fpath))
          fpath=trim(dirout)//'/'//trim(foldout)//'/'//years//
     +          trim(area(ig))//'_Abovelevel.txt'
          open(18,file=trim(fpath))
          fpath=trim(dirout)//'/'//trim(foldout)//'/'//years//
     +      trim(area(ig))//'_Tropsphere.txt'
          open(19,file=trim(fpath))
          fpath=trim(dirout)//'/'//trim(foldout)//'/'//years//
     +      trim(area(ig))//'_AllLevels.txt'
          open(20,file=trim(fpath))
          fpath=trim(dirout)//'/'//trim(foldout)//'/'//years//
     +      trim(area(ig))//'_RAW.txt'
          open(21,file=trim(fpath))
          do jj=15,21
          write(jj,908)'ID','T_ls(k/day)','Q_ls(k/day)','U(m/s)',
     +      'V(m/s)','moisture(kg/kg)','HGT(m)',
     +      'AIR(K)','Adj_omega(pa/s)',
     +      'Theta(K)','Q1(k/day)','Q2(k/day)','HADQ(K/day)',
     +      'VADQ(K/day)','TCHQ(K/day)','HADT(K/day)','VADT(K/day)',
     +      'TCHT(K/day)','Ori_omega(pa/s)'
          enddo
         fpath=trim(dirout)//'/'//trim(foldout)//'/'//years//
     +    trim(area(ig))//'_OMGCOMPSOUT.txt'
          open(22,file=trim(fpath))
          write(22,*)'OMGO DIV'
!     +      'prate','cprate','DLR(W/m^2)','ULR(W/m^2)'
!
!          fpath=trim(dirout)//'/'//trim(foldout)//'/'//years//
!     +    trim(area(ig))//'_dyn.txt'
!          open(23,file=trim(fpath))
!          write(23,909)'DATE','HOUR','VorX','VorY',
!     +      'VorZ','DIV'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ----------   GET THE GEO --------------
          fpath=trim(dirin)//'SRFX2.5/surface_geo.nc'
          GEO=0.0
          Call ERAGEO(fpath,X1(ig), X2(ig), Y1(ig), Y2(ig),
     +                         IXX(ig), IYY(ig), GEO)
!----------------------------------------
!                get the surface dateset
          print*,'GEO',GEO(3,3)  
          ids=1
          do i=1979,iyear-1
            nddd=365
            if(mod(iyear,4)==0.and.mod(iyear,100)/=0)then
              nddd=366
            elseif(mod(iyear,400)==0)then
              nddd=366
            endif
            ids=ids+nddd
          enddo
          tmp2d0=0.0
          do iv=1,5
            do ih=1,4 
              fpath=trim(dirin)//'SRFX2.5\'//'1979-2014_'//
     +             trim(name2d(iv))//'_'//hour(ih)//'.nc'
              Call ReadSRF(fpath, ih, ids, ndd, X1(ig), X2(ig), 
     +                   Y1(ig), Y2(ig),  IXX(ig), IYY(ig), tmp2d0)
              do i1=1,IXX(ig)
                do i2=1,IYY(ig)
                  do i3=1,ndd
                    tmp2d(iv,i1,i2,ih,i3)=tmp2d0(i1,i2,i3)
                  enddo
                enddo
              enddo               
            enddo 
          enddo
          dat3d=0.0
          do iv=1,6
            do ih=1,4
              fpath=trim(dirin)//'X2.5\'//years//'_'//
     +               trim(name3d(iv))//'_'//hour(ih)//'.nc'
!			write(66,*)ndd,fpath
!			pause
              call Read3d(fpath, ndd, X1(ig), X2(ig),
     +           Y1(ig), Y2(ig), IXX(ig), IYY(ig), tmp3d) !REAL DAT3D(6,ndx,ndy,nz,4,nday)
              DAT3D(:,:,:,:,IH,IV)=TMP3D(:,:,:,:) !tmp3d(ndx,ndy,nz,nday)  !DAT3D(ndx,ndy,nz,nday,4,6)
!	        print*,'AAA'
            enddo
          enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ibound=0
          do 1109 it=1,ndd*4 
            hgt  = 0.0
            qv   = 0.0
            air  = 0.0
            uwnd = 0.0
            vwnd = 0.0
            omega= 0.0
!	          print*,it
            if (it.eq.1)then  !!!! time boundary condition 1
              id=1
              ibound=-2
!!!!!!!!!!!!!   SRF  --------------------------------
              do itt=it,it+2
                ih=mod((itt-1),4)+1
			  write(404,*)'2D',it,id,ih
                do iv=1,5
                  do i1=1,IXX(ig)
                    do i2=1,IYY(ig)                     
                      if(iv.eq.1)then
                        psf(i1,i2,itt)=tmp2d(iv,i1,i2,ih,id)
                      elseif(iv.eq.2)then
                        t2m(i1,i2,itt)=tmp2d(iv,i1,i2,ih,id)
                      elseif(iv.eq.3)then
                        td2m(i1,i2,itt)=tmp2d(iv,i1,i2,ih,id)
                        if(tmp2d(3,i1,i2,ih,id).LT.0.+273.16)then
                          tc=tmp2d(3,i1,i2,ih,id)-273.15
                          es=6.112*exp(17.67*tc/(tc+243.5))    !!! Bolton
                        else
                          tk=tmp2d(3,i1,i2,ih,id)
                          es=6.1078*exp( 
     +                           17.2693882*(tk-273.16)/(tk-35.86) ) !!Tetens 
                        endif
                        qvs(i1,i2,itt)=0.622*es*100.
     +                          /tmp2d(1,i1,i2,ih,id)  !
                      elseif(iv.eq.4)then
                        u10m(i1,i2,itt)=tmp2d(iv,i1,i2,ih,id)
                      elseif(iv.eq.5)then
                        v10m(i1,i2,itt)=tmp2d(iv,i1,i2,ih,id)
                      endif 
                    enddo
                  enddo
                enddo
              enddo
!!!---------------3D ---------------------------------              
              do itt=it,it+2
                ih=mod((itt-1),4)+1
                write(404,*)'3D',it,id,ih,'AAA'
                do iv=1,n3d
                  if(iv.eq.1)then
                    hgt(:,:,:,itt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.2)then
                    qv(:,:,:,itt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.3)then
                    air(:,:,:,itt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.4)then
                    uwnd(:,:,:,itt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.5)then
                    vwnd(:,:,:,itt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.6)then
                    omega(:,:,:,itt)=DAT3D(:,:,:,ID,IH,IV)
                  endif 
                enddo                         
              enddo
!!!! ---------------------time boundary condition 2
            elseif (it.eq.2)then  !!!! time boundary condition 2
              id=1
              ibound=-1
!!!!!!!!!!!!!   SRF  --------------------------------
              do itt=it-1,it+1
                ih=mod((itt-1),4)+1
                do iv=1,5
                  do i1=1,IXX(ig)
                    do i2=1,IYY(ig)
                      if(iv.eq.1)then
                        psf(i1,i2,itt)=tmp2d(iv,i1,i2,ih,id)
                      elseif(iv.eq.2)then
                        t2m(i1,i2,itt)=tmp2d(iv,i1,i2,ih,id)
                      elseif(iv.eq.3)then
                        td2m(i1,i2,itt)=tmp2d(iv,i1,i2,ih,id)
                        if(tmp2d(3,i1,i2,ih,id).LT.0.+273.16)then
                          tc=tmp2d(3,i1,i2,ih,id)-273.15
                          es=6.112*exp(17.67*tc/(tc+243.5))    !!! Bolton
                        else
                          tk=tmp2d(3,i1,i2,ih,id)
                          es=6.1078*exp( 
     +                         17.2693882*(tk-273.16)/(tk-35.86) ) !!Tetens 
                        endif
                        qvs(i1,i2,itt)=0.622*es*100.
     +                                /tmp2d(1,i1,i2,ih,id)  !
                      elseif(iv.eq.4)then
                        u10m(i1,i2,itt)=tmp2d(iv,i1,i2,ih,id)
                      elseif(iv.eq.5)then
                        v10m(i1,i2,itt)=tmp2d(iv,i1,i2,ih,id)
                      endif 
                    enddo
                  enddo
                enddo
              enddo
!!!---------------3D ---------------------------------              
              do itt=it-1,it+1
                ih=mod((itt-1),4)+1
                do iv=1,n3d
!                  fpath=trim(dirin)//'X2.5\'//years//'_'//
!     +                    trim(name3d(iv))//'_'//hour(ih)//'.nc'
!                  tmp3d=0.             
!                  call Read3d(fpath, id,  ndd, X1(ig), X2(ig),
!     +                       Y1(ig), Y2(ig), IXX(ig), IYY(ig), tmp3d)
                  if(iv.eq.1)then
                    hgt(:,:,:,itt)=DAT3D(:,:,:,ID,IH,IV)!DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.2)then
                    qv(:,:,:,itt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.3)then
                    air(:,:,:,itt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.4)then
                    uwnd(:,:,:,itt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.5)then
                    vwnd(:,:,:,itt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.6)then
                    omega(:,:,:,itt)=DAT3D(:,:,:,ID,IH,IV)
                  endif 
                enddo                         
              enddo
!!!! -----------------time boundary condition 3
            elseif(it.eq.(ndd*4))then
              ibound=2
              id=ndd
!!!!!!!!!!!!!   SRF  --------------------------------
              ittt=1
              do itt=it-2,it
                ih=mod((itt-1),4)+1
                do iv=1,5
                  do i1=1,IXX(ig)
                    do i2=1,IYY(ig)
                      if(iv.eq.1)then
                        psf(i1,i2,ittt)=tmp2d(iv,i1,i2,ih,id)
                      elseif(iv.eq.2)then
                        t2m(i1,i2,ittt)=tmp2d(iv,i1,i2,ih,id)
                      elseif(iv.eq.3)then
                        td2m(i1,i2,ittt)=tmp2d(iv,i1,i2,ih,id)
                        if(tmp2d(3,i1,i2,ih,id).LT.0.+273.16)then
                          tc=tmp2d(3,i1,i2,ih,id)-273.15
                          es=6.112*exp(17.67*tc/(tc+243.5))    !!! Bolton
                        else
                          tk=tmp2d(3,i1,i2,ih,id)
                          es=6.1078*exp( 
     +                        17.2693882*(tk-273.16)/(tk-35.86) ) !!Tetens 
                        endif
                        qvs(i1,i2,ittt)=0.622*es*100.
     +                                   /tmp2d(1,i1,i2,ih,id)  !
                      elseif(iv.eq.4)then
                        u10m(i1,i2,ittt)=tmp2d(iv,i1,i2,ih,id)
                      elseif(iv.eq.5)then
                        v10m(i1,i2,ittt)=tmp2d(iv,i1,i2,ih,id)
                      endif 
                    enddo
                  enddo
                enddo
                ittt=ittt+1
              enddo              
!!!!!-----------------  3D ------------------------------------------------
              ittt=1
              do itt=it-2,it
                ih=mod((itt-1),4)+1
                do iv=1,n3d
!                  fpath=trim(dirin)//'X2.5/'//years//'_'//
!     +                    trim(name3d(iv))//'_'//hour(ih)//'.nc'
!                  tmp3d=0             
!                  call Read3d(fpath, id,  ndd, X1(ig), X2(ig),
!     +                       Y1(ig), Y2(ig), IXX(ig), IYY(ig), tmp3d)
                  if(iv.eq.1)then
                    hgt(:,:,:,ittt)=DAT3D(:,:,:,ID,IH,IV) !DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.2)then
                    qv(:,:,:,ittt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.3)then
                    air(:,:,:,ittt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.4)then
                    uwnd(:,:,:,ittt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.5)then
                    vwnd(:,:,:,ittt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.6)then
                    omega(:,:,:,ittt)=DAT3D(:,:,:,ID,IH,IV)
                  endif 
                enddo
                ittt=ittt+1                         
              enddo
!!!! -----------------time boundary condition 4
            elseif(it.eq.(ndd*4-1))then
              ibound=1
              id=ndd
!!!!!!!!!!!!!   SRF  --------------------------------
              ittt=1
              do itt=it-1,it+1
                ih=mod((itt-1),4)+1
                do iv=1,5
                  do i1=1,IXX(ig)
                    do i2=1,IYY(ig)
                      if(iv.eq.1)then
                        psf(i1,i2,ittt)=tmp2d(iv,i1,i2,ih,id)
                      elseif(iv.eq.2)then
                        t2m(i1,i2,ittt)=tmp2d(iv,i1,i2,ih,id)
                      elseif(iv.eq.3)then
                        td2m(i1,i2,ittt)=tmp2d(iv,i1,i2,ih,id)
                        if(tmp2d(3,i1,i2,ih,id).LT.0.+273.16)then
                          tc=tmp2d(3,i1,i2,ih,id)-273.15
                          es=6.112*exp(17.67*tc/(tc+243.5))    !!! Bolton
                        else
                          tk=tmp2d(3,i1,i2,ih,id)
                          es=6.1078*exp( 
     +                         17.2693882*(tk-273.16)/(tk-35.86) ) !!Tetens 
                        endif
                        qvs(i1,i2,ittt)=0.622*es*100.
     +                                /tmp2d(1,i1,i2,ih,id)  !
                      elseif(iv.eq.4)then
                        u10m(i1,i2,ittt)=tmp2d(iv,i1,i2,ih,id)
                      elseif(iv.eq.5)then
                        v10m(i1,i2,ittt)=tmp2d(iv,i1,i2,ih,id)
                      endif 
                    enddo
                  enddo
                enddo
                ittt=ittt+1
              enddo              
!!!!!-----------------  3D ------------------------------------------------
              ittt=1
              do itt=it-1,it+1
                ih=mod((itt-1),4)+1
                do iv=1,n3d
!                  fpath=trim(dirin)//'X2.5/'//years//'_'//
!     +                   trim(name3d(iv))//'_'//hour(ih)//'.nc'
!                  tmp3d=0             
!                  call Read3d(fpath, id,  ndd, X1(ig), X2(ig),
!     +                       Y1(ig), Y2(ig), IXX(ig), IYY(ig), tmp3d)
                  if(iv.eq.1)then
                    hgt(:,:,:,ittt)=DAT3D(:,:,:,ID,IH,IV) !DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.2)then
                    qv(:,:,:,ittt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.3)then
                    air(:,:,:,ittt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.4)then
                    uwnd(:,:,:,ittt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.5)then
                    vwnd(:,:,:,ittt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.6)then
                    omega(:,:,:,ittt)=DAT3D(:,:,:,ID,IH,IV)
                  endif 
                enddo
                ittt=ittt+1                         
              enddo
!!!! 
            else
              ibound=0  !!! normal 
              ih=mod((it-1),4)+1
              if(ih.eq.1)then   !!! it is hour 00
                id3=it/4+1
                id4=id3
                id5=id3
                id1=id3-1
                id2=id3-1
              elseif(ih.eq.2)then !!! it is hour 06
                id3=it/4+1
                id4=id3
                id5=id3
                id1=id3-1
                id2=id3
              elseif(ih.eq.3)then !!! it is hour 12
                id3=it/4+1
                id4=id3
                id5=id3+1
                id1=id3
                id2=id3
              elseif(ih.eq.4)then !!! it is hour 18
                id3=it/4
                id4=id3+1
                id5=id3+1
                id1=id3
                id2=id3
              endif
!!!!!!!!!!!!!   SRF  --------------------------------
              ittt=1
              do itt=it-2,it+2
                ih=mod((itt-1),4)+1
                if(ittt.eq.1)id=id1
                if(ittt.eq.2)id=id2
                if(ittt.eq.3)id=id3
                if(ittt.eq.4)id=id4
                if(ittt.eq.5)id=id5
                write(404,*)'2D, normal',it,id,ih
                do iv=1,5
                  do i1=1,IXX(ig)
                    do i2=1,IYY(ig)
                      if(iv.eq.1)then
                        psf(i1,i2,ittt)=tmp2d(iv,i1,i2,ih,id)
                      elseif(iv.eq.2)then
                        t2m(i1,i2,ittt)=tmp2d(iv,i1,i2,ih,id)
                      elseif(iv.eq.3)then
                        td2m(i1,i2,ittt)=tmp2d(iv,i1,i2,ih,id)
                        if(tmp2d(3,i1,i2,ih,id).LT.0.+273.16)then
                          tc=tmp2d(3,i1,i2,ih,id)-273.15
                          es=6.112*exp(17.67*tc/(tc+243.5))    !!! Bolton
                        else
                          tk=tmp2d(3,i1,i2,ih,id)
                          es=6.1078*exp( 
     +                         17.2693882*(tk-273.16)/(tk-35.86) ) !!Tetens 
                        endif
                        qvs(i1,i2,ittt)=0.622*es*100.
     +                                /tmp2d(1,i1,i2,ih,id)  ! !
                      elseif(iv.eq.4)then
                        u10m(i1,i2,ittt)=tmp2d(iv,i1,i2,ih,id)
                      elseif(iv.eq.5)then
                        v10m(i1,i2,ittt)=tmp2d(iv,i1,i2,ih,id)
                      endif 
                    enddo
                  enddo
                enddo
                ittt=ittt+1
              enddo
!!!!  ---------   3D ------------------------                            
              ittt=1
              do itt=it-2,it+2
                ih=mod((itt-1),4)+1
                if(ittt.eq.1)id=id1
                if(ittt.eq.2)id=id2
                if(ittt.eq.3)id=id3
                if(ittt.eq.4)id=id4
                if(ittt.eq.5)id=id5
                write(404,*)'3D, normal',it,id,ih,'BBB'
!                pause
                do iv=1,n3d
!                  fpath=trim(dirin)//'X2.5/'//years//'_'//
!     +                    trim(name3d(iv))//'_'//hour(ih)//'.nc'
!                  tmp3d=0             
!                  call Read3d(fpath, id, ndd,  X1(ig), X2(ig), 
!     +                       Y1(ig), Y2(ig),  IXX(ig), IYY(ig), tmp3d)
                  if(iv.eq.1)then
                    hgt(:,:,:,ittt)=DAT3D(:,:,:,ID,IH,IV) !DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.2)then
                    qv(:,:,:,ittt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.3)then
                    air(:,:,:,ittt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.4)then
                    uwnd(:,:,:,ittt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.5)then
                    vwnd(:,:,:,ittt)=DAT3D(:,:,:,ID,IH,IV)
                  elseif(iv.eq.6)then
                    omega(:,:,:,ittt)=DAT3D(:,:,:,ID,IH,IV)
                  endif 
                enddo
                ittt=ittt+1                         
              enddo
            endif  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(ibound.eq.0)then
              NTT=5
            else
              NTT=3
            endif
            do i1=1,IXX(ig)
              do i2=1, IYY(ig)
                do iz=1,nz  !!!! get the theta
                  do itt=1,NTT
                    ptmp=plv(iz)
                    p00=1000.
                    TK=air(i1,i2,iz,itt)
                    Theta(i1,i2,iz,itt)=TK*(p00/ptmp)**0.286
                  enddo
                enddo
              enddo
            enddo
            DY=2.5*111.17E3
		        PIE=3.141592657
            DO J=1,IYY(IG)
              Y=((IYY(IG)-J)*2.5+bslat(IG))/180.
              DX(J)=COS(Y*PIE)*2.5*111.17E3
            ENDDO
            DIVO=0.0
            OMGO=0.0
            DIV=0.0
            OMGS=0.0  !!! surface omega
            OMG=0.0
! if(lo==4) ikk=1
            ikt=27  !!!! pressure is 100hpa
            iks=14   !!!! pressure is 600 hpa
            call DIVVOR(PLV,    PLW      , OMGS ,  DIVO, IXX(ig), 
     +                  IYY(ig),bslat(ig), ibound, ikt, iks) !DIVVOR(PM,WSFC,DIVO,NXDIR,NYDIR,YMIN,IT)
            call VERVEL(PLV    , PLW   , OMGS, DIVO, OMGO, IXX(ig),
     +                  IYY(ig), ibound, ikt, iks)
            call ADJDV(PLV,PLW,DIVO,OMGO,DIV,IXX(ig),
     +                 IYY(ig),ibound, ikt, iks)   !ADJDV(PM,PW,DIVO,OMGO,DIV,NXDIR,NYDIR,ibound,ikk)
            call VERVEL(PLV    ,  PLW   , OMGS,  DIV, OMG, IXX(ig),
     +                  IYY(ig),  ibound, ikt, iks)
            call Q1Q2(PLV    , PLW   , OMG, IXX(ig),
     +                IYY(ig), ibound, ikt, iks)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            tls_xy = 0.0
            qls_xy = 0.0
            u_xy   = 0.0
            the_xy = 0.0
            v_xy   = 0.0
            rh_xy  = 0.0
            gz_xy  = 0.0
            Q1_xy  = 0.0
            air_xy = 0.0
            Q2_xy  = 0.0
            qv_xy  = 0.0
            omg_xy = 0.0
            psrf_xy= 0.0
            t2m_xy = 0.0
            u10m_xy= 0.0 
            v10m_xy= 0.0 
            qv2m_xy= 0.0
            tvs_xy = 0.0
            geo_xy = 0.0
            omgs_xy = 0.0
            vq1_xy = 0.0
            vq2_xy = 0.0
            omgco_xy= 0.0
            DIV_xy   =  0.0
            if(ibound.eq.0)then
              ICT=3
            elseif(ibound.eq.-1)then
              ICT=2
            elseif(ibound.eq.-2)then
              ICT=1
            elseif(ibound.eq.1)then
              ICT=2
            elseif(ibound.eq.2)then
              ICT=3
            endif
            Tscale=3600.*24.  !!! convert to K/day   Q1 and temperature forcing
            Qscale=3600.*24.*hlat/cp   !!! convert to k/day   Q2 and moisture forcing
            do i1=2,IXX(ig)-1
              do i2=2,IYY(ig)-1
                geo_xy = geo_xy + geo(i1,i2)
                psrf_xy=psrf_xy + psf(i1,i2,ICT)
                t2m_xy =t2m_xy  + t2m(i1,i2,ICT)
                u10m_xy=u10m_xy + u10m(i1,i2,ICT)
                v10m_xy=v10m_xy + v10m(i1,i2,ICT)
                qv2m_xy=qv2m_xy + qvs(i1,i2,ICT)
                omgs_xy=omgs_xy + OMGS(i1,i2)
                vq1_xy = vq1_xy + VQ1(i1,i2)
                vq2_xy = vq2_xy + VQ2(i1,i2)
                ptmp=psf(i1,i2,ICT)
                p00=1000.
                TK=t2m(i1,i2,ICT)
                theta_srf=TK*(p00/ptmp)**0.286
                tvs_xy = tvs_xy + theta_srf
                do iz=1,nz
                  tls_xy(iz) = tls(i1,i2,iz)*Tscale+tls_xy(iz)
                  qls_xy(iz) = qls(i1,i2,iz)*Qscale+qls_xy(iz)
                  u_xy(iz)   = uwnd(i1,i2,iz,ICT)+u_xy(iz)
                  the_xy(iz) = theta(i1,i2,iz,ICT)+the_xy(iz)
                  v_xy(iz)   = vwnd(i1,i2,iz,ICT)+v_xy(iz)
                  Q1_xy(iz)  = q1(i1,i2,iz)*Tscale+q1_xy(iz)
                  air_xy(iz) = air(i1,i2,iz,ICT)+air_xy(iz)
                  Q2_xy(iz)  = q2(i1,i2,iz)*Qscale+q2_xy(iz)
                  qv_xy(iz)  = qv(i1,i2,iz,ICT)+qv_xy(iz)
                  omg_xy(iz) = OMG(i1,i2,iz)+omg_xy(iz)
                  omgo_xy(iz)= omega(i1,i2,iz,ICT)+omgo_xy(iz)
                  gz_xy(iz)  = hgt(i1,i2,iz,ICT)+gz_xy(iz)
                  omgco_xy(iz) = OMGO(i1,i2,iz) +omgco_xy(iz)
                  DIV_xy(iz)   =   DIV(i1,i2,iz) + div_xy(iz)
!				IF(Q1_XY(IZ).GT. 0.)THEN
!					PRINT*,Q1(I1,I2,IZ)
!					PRINT*,TSCALE,QSCALE
!					PAUSE
!				ENDIF
                enddo
              enddo
            enddo
            XYC=(IXX(ig)-2.)*(IYY(ig)-2.)
            psrf_xy=psrf_xy/XYC
            t2m_xy =t2m_xy/XYC
            u10m_xy=u10m_xy/XYC
            v10m_xy=v10m_xy/XYC
            qv2m_xy=qv2m_xy/XYC
            tvs_xy = tvs_xy/XYC
            geo_xy = geo_xy/XYC
            vq1_xy =  VQ1_XY/XYC
            vq2_xy =  VQ2_XY/XYC
            write(14,904)      it, psrf_xy, geo_xy, u10m_xy
     +                  , v10m_xy, omgs_xy, t2m_xy, tvs_xy       
     +                  , qv2m_xy, vq1_xy , vq2_xy
            write(11,*)it
            write(13,*)it
            do iz=1,nz
!			IF(Q1_XY(IZ).GT. 0.)THEN
!					PRINT*,Q1_XY(IZ),IZ
!					PRINT*,XYC
!					PAUSE
!			ENDIF
              tls_xy(iz) = tls_xy(iz)/XYC
              qls_xy(iz) = qls_xy(iz)/XYC
              u_xy(iz)   = u_xy(iz)/XYC
              the_xy(iz) = the_xy(iz)/XYC
              v_xy(iz)   = v_xy(iz)/XYC
              Q1_xy(iz)  = q1_xy(iz)/XYC
              air_xy(iz) = air_xy(iz)/XYC
              Q2_xy(iz)  = q2_xy(iz)/XYC
              qv_xy(iz)  = qv_xy(iz)/XYC
              omg_xy(iz) = omg_xy(iz)/XYC   ! adjust
              omgo_xy(iz)= omgo_xy(iz)/XYC ! origional
              gz_xy(iz)  = gz_xy(iz)/XYC
              omgco_xy(iz)=omgco_xy(iz)/XYC
              DIV_xy(iz) = DIV_xy(iz)/XYC
              write(11,901) plv(iz)  , gz_xy(iz) , u_xy(iz)  ,v_xy(iz)
     +                   , omg_xy(iz), air_xy(iz), the_xy(iz),qv_xy(iz)
              write(13,903) plv(iz)  , q1_xy(iz) , q2_xy(iz)
              write(22,'(2e12.4)') omgco_xy(iz) , DIV_xy(iz)    
            enddo
            write(12,902) it, (tls_xy(ikk),ikk=1,nz),
     +                          (qls_xy(ikk),ikk=1,nz)  
            iks=1  ; ike=12    ! 1000 to 700
            XY_OUT=0.0
            call ERA_SAS(   PLV  , DP    , OMG , IXX(ig),
     +                   IYY(ig), ibound, iks , ike    , XY_OUT)
            write(15,906)IT,(XY_OUT(kk),kk=1,18)
            iks=13  ; ike=20    ! 650 to 300
            XY_OUT=0.0
            call ERA_SAS(   PLV, DP    , OMG , IXX(ig),
     +                   IYY(ig), ibound, iks , ike    , XY_OUT)
            write(16,906)IT,(XY_OUT(kk),kk=1,18)
            iks=21  ; ike=27    ! 300to 100
            XY_OUT=0.0
            call ERA_SAS(   PLV, DP    , OMG , IXX(ig),
     +                   IYY(ig), ibound, iks , ike    , XY_OUT)
            write(17,906)IT,(XY_OUT(kk),kk=1,18)
            iks=28  ; ike=37    ! 100 to 1
            XY_OUT=0.0
            call ERA_SAS(   PLV, DP    , OMG , IXX(ig),
     +                   IYY(ig), ibound, iks , ike    , XY_OUT)
            write(18,906)IT,(XY_OUT(kk),kk=1,18)
            iks=14 ; ike=27    ! 600 to 100
            XY_OUT=0.0
            call ERA_SAS(   PLV, DP    , OMG , IXX(ig),
     +                   IYY(ig), ibound, iks , ike    , XY_OUT)
            write(19,906)IT,(XY_OUT(kk),kk=1,18)
            iks=1  ; ike=37    ! 1000 to 1
            XY_OUT=0.0
            call ERA_SAS(   PLV, DP    , OMG , IXX(ig),
     +                   IYY(ig), ibound, iks , ike    , XY_OUT)
            write(20,906)IT,(XY_OUT(kk),kk=1,18)
!------------ no levels ave
            do iz=1,nz
              iks=iz ; ike=iz
              XY_OUT = 0.0
            call ERA_SAS(   PLV, DP    , OMG , IXX(ig),
     +                   IYY(ig), ibound, iks , ike    , XY_OUT)
              write(21,906)IT,(XY_OUT(kk),kk=1,18)
            enddo
1109      continue ! it
          do ifl=11,21
            close(ifl)
          enddo
1110    continue !  ig
1111  continue ! iyear
!---------------------------------------------------------------------
901   format(1X,8(1X,e12.4))
902   format(1X,I4,74(1X,e12.4))  ! 37*2
903   format(1X,F9.3,2(1X,e12.4))
904   format(1X,I4,10(1X,e12.4))
905   format(1X,A4,1X,I4,2I2.2,1X,A2,1X,I4,2I2.2,1X,A14)
906   format(1X,I4,18(1X,e12.4))
907   format(1X,A16,1X,6(1X,e12.4))
908   format(1X,19(1X,A11))
!909   format(1X,A10,1X,A5,6(1X,A11))
!910   format(1X,A16,1X,4(1X,e12.4))

      END PROGRAM
!!!    
!------------------  the following code copy from hmbudps.f-----------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  DIVVOR COMPUTES DERIVED DIVERGENCE AND VORTICITY
CC  Note the dx and dy are vector, when the HAD is handle with, carefully note the data store squence
CC  and the array  subscript.
CC  NCEP: EAST-WEAT grids, the longitude array(x+1,y) is great than array(x,y)
CC        South North: the latitude array(x,y+1) is smaller than array(x,y)
      SUBROUTINE DIVVOR(PM,PW,WSFC,DIVO,NXDIR,NYDIR,YMIN,ibound,ikt,iks)
      PARAMETER (GRID=2.50)
      integer,parameter :: nz=37 !! add surface
      integer,parameter :: n3d=6 !1u 1v 3t  4oemga 5RH 6HGT
      integer,parameter :: ndx=10
      integer,parameter :: ndy=10
      integer,parameter :: nt=5
      integer NXDIR,NYDIR,it,ktop
      DIMENSION WSFC(10,10),DIVO(10,10,nz),PM(nz),PW(nz)
      COMMON/D3D/ air(ndx,ndy,nz,nt),hgt(ndx,ndy,nz,nt),
     +            qv(ndx,ndy,nz,nt),uwnd(ndx,ndy,nz,nt), !1u 2v 3t  4oemga 5RH 6HGT
     +            vwnd(ndx,ndy,nz,nt),omega(ndx,ndy,nz,nt), ! !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top 11cpr 12DLR 13ULR,14t2m
     +            theta(ndx,ndy,nz,nt)
      COMMON/SRF/ t2m(ndx,ndy,nt),td2m(ndx,ndy,nt),
     +            u10m(ndx,ndy,nt),v10m(ndx,ndy,nt),
     +            psf(ndx,ndy,nt),qvs(ndx,ndy,nt),
     +            Geo(ndx,ndy)
      DATA PIE/0.0174532925/   !! pie/180 mean every dreeg to 
CC
      DY=1./(GRID*111.17E3)
      if(ibound.eq.0)then
        ITT=3
      elseif(ibound.eq.-1)then
        ITT=2
      elseif(ibound.eq.-2)then
        ITT=1
      elseif(ibound.eq.1)then
        ITT=2
      elseif(ibound.eq.2)then
        ITT=3
      endif
      DO 100 I=1,NXDIR
      DO 100 J=1,NYDIR
        Y=((NYDIR-J)*GRID+YMIN)/180.
        DX=1./(COS(Y*PIE)*GRID*111.17E3)
CC ******** COMPUTE TERRAIN INDUCED OMEGAS FOR LOWEST LAYER
        IF(I.EQ.1 .OR. I.EQ.NXDIR) GO TO 95
!        WX=data2d(i,j,it,4)*DX/2.*(hgt(I+1,J)-hgt(I-1,J))
        WX=u10m(I,J,ITT)*DX/2.*(Geo(I+1,J)-Geo(I-1,J))
   95   IF(I.EQ.1)WX=u10m(I,J,ITT)*DX*(Geo(I+1,J)-Geo(I,J)) ! WX=data2d(i,j,it,4)*DX*(hgt(I+1,J)-hgt(I,J))
        IF(I.EQ.NXDIR)WX=u10m(I,J,ITT)*DX*(Geo(I,J)-Geo(I-1,J))! WX=data2d(i,j,it,4)*DX*(hgt(I,J)-hgt(I-1,J))
        IF(J.EQ.1 .OR. J.EQ.NYDIR) GO TO 96
!      WY=data2d(i,j,it,5)*DY/2.*(hgt(I,J-1)-hgt(I,J+1))
        WY=v10m(I,J,ITT)*DY/2.*(Geo(I,J-1)-Geo(I,J+1))     ! MUST NOTE THE J, FROM NORTH TO SOUTH OR SOUTH TO NORTH
   96   IF(J.EQ.1)WY=v10m(I,J,ITT)*DY*(Geo(I,J)-Geo(I,J+1))! WY=data2d(i,j,it,5)*DY*(hgt(I,J)-hgt(I,J+1))
        IF(J.EQ.NYDIR)WY=v10m(I,J,ITT)*DY*(Geo(I,J-1)-Geo(I,J)) ! WY=data2d(i,j,it,5)*DY*(hgt(I,J-1)-hgt(I,J))
!      WSFC(I,J,ITT)=-data2d(i,j,it,9)/100*9.8/
!     +  (287.*data2d(i,j,it,6))*(WX+WY)
        WSFC(I,J)=-psf(I,J,ITT)/100.*9.8/
     +      (287.*t2m(I,J,ITT))*(WX+WY)
	  
CC UNIT: MB/S FOR WSFC
!      p_top=data2d(I,J,it,10)/100.
!       pres=data2d(I,J,it,9)/100.
        pres=PSF(I,J,ITT)/100.
!       ktop=10  ! 100hPa
!       do K=1,lev-1
!        if(pm(k)>p_top.and.pw(k)<P_top)then
!          ktop=k-1-ikk
!          goto 101
!        elseif(pw(k)>P_top.and.pm(k+1)<p_top)then
!            ktop=k-ikk
!             goto 101
!        endif
!        enddo
!  101 CONTINUE
        ktop=ikt
CC ******** COMPUTE DIVERGENCE AND VORTICITY
        DO 200 K=iks,ktop   !!!! 14 600hPa
!!!!
          IF(K.NE.1) GO TO 400
c if(lo<3.and.K<4) Go To 500
c      if(lo<3.and.k==4) Go To 400  ! tiblet surface is very high from   600
C%%%%%%%%%%%  SURFACE DIVERGENCE
CC BOUNDARY CALCULATIONS
          IF(I.EQ.1) THEN
!        DUX=(U(I+1,J,K)-U(I,J,K))*DX  
            DUX=(uwnd(I+1,J,K,ITT)-uwnd(I,J,K,ITT))*Dx
          ELSEIF(I.EQ.NXDIR) THEN
!        DUX=(U(I,J,K)-U(I-1,J,K))*DX
            DUX=(uwnd(I,J,K,ITT)-uwnd(I-1,J,K,ITT))*Dx
          ELSE
CC CENTER FINITE DIFFERENCE CALCULATIONS
!       DUX=(U(I+1,J,K)-U(I-1,J,K))*DX/2.
            DUX=(uwnd(I+1,J,K,ITT)-uwnd(I-1,J,K,ITT))*Dx/2.
          ENDIF
CC BOUNDARY CALCULATIONS
          IF(J.EQ.1) THEN
!        DVY=(V(I,J+1,K)-V(I,J,K))*DY
            DUY=-(vwnd(I,J+1,K,ITT)-vwnd(I,J,K,ITT))*DY
          ELSEIF(J.EQ.NYDIR) THEN
!        DVY=(V(I,J,K)-V(I,J-1,K))*DY
            DUY=-(vwnd(I,J,K,ITT)-vwnd(I,J-1,K,ITT))*DY
          ELSE
CC CENTER FINITE DIFFERENCE CALCULATIONS
!        DVY=(V(I,J+1,K)-V(I,J-1,K))*DY/2.
            DUY=-(vwnd(I,J+1,K,ITT)-vwnd(I,J-1,K,ITT))*DY/2.
          ENDIF
          GO TO 500
C%%%%%%%%%%%   UPPER LAYER DIVERGENCE
CC BOUNDARY CALCULATIONS
400       IF(I.EQ.1 .OR. (psf(I-1,J,ITT)/100.).LT.PM(K)) THEN
!        DUX=(U(I+1,J,K)-U(I,J,K))*DX
            DUX=(uwnd(I+1,J,K,ITT)-uwnd(I,J,K,ITT))*DX
          ELSEIF(I.EQ.NXDIR .OR. (psf(I+1,J,ITT)/100).LT.PM(K)) THEN
!        DUX=(U(I,J,K)-U(I-1,J,K))*DX
            DUX=(uwnd(I,J,K,ITT)-uwnd(I-1,J,K,ITT))*DX
          ELSE
CC CENTER FINITE DIFFERENCE CALCULATIONS
!        DUX=(U(I+1,J,K)-U(I-1,J,K))*DX/2.
            DUX=(uwnd(I+1,J,K,ITT)-uwnd(I-1,J,K,ITT))*DX/2.
          ENDIF
CC BOUNDARY CALCULATIONS
          IF(J.EQ.1 .OR. (psf(I,J-1,ITT)/100.).LT.PM(K)) THEN
!        DVY=(V(I,J+1,K)-V(I,J,K))*DY
            DUY=-(vwnd(I,J+1,K,ITT)-vwnd(I,J,K,ITT))*DY
          ELSEIF(J.EQ.NYDIR .OR. (psf(I,J+1,ITT)/100.).LT.PM(K)) THEN
!        DVY=(V(I,J,K)-V(I,J-1,K))*DY
            DUY=-(vwnd(I,J,K,ITT)-vwnd(I,J-1,K,ITT))*DY
          ELSE
CC CENTER FINITE DIFFERENCE CALCULATIONS
            DVY=-(vwnd(I,J+1,K,ITT)-vwnd(I,J-1,K,ITT))*DY/2.
!   DUY=(data3d(I,J+1,K,it,2)-data3d(I,J-1,K,it,2))*DY/2.
           ENDIF
C
C  COMPUTE DIVERGENCE
C
  500      DIVO(I,J,K)=DUX+DVY
c      ccc=-1.0/0
! print*,ccc
c      if(divo(I,j,k,it)==ccc)then
c        print*,DUX,DUY,'FFF'
c   stop
c endif
CC  UNITS:   DIVO(1/S)
  200 CONTINUE
  100 CONTINUE
      RETURN
      END SUBROUTINE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  VERVEL COMPUTES VERTICAL VELOCITY
      SUBROUTINE VERVEL(PM,PW,WSFC,DIVO,OMGO,NXDIR,NYDIR,ibound,ikt,iks)
      PARAMETER (GRID=2.50)
      integer,parameter :: nz=37 !! add surface
      integer,parameter :: n3d=6 !1u 1v 3t  4oemga 5RH 6HGT
      integer,parameter :: ndx=10
      integer,parameter :: ndy=10
      integer,parameter :: nt=5
      integer NXDIR,NYDIR,ibound
      DIMENSION WSFC(10,10),DIVO(10,10,nz),
     +          PM(nz),PW(nz), OMGO(10,10,nz)
      COMMON/D3D/ air(ndx,ndy,nz,nt),hgt(ndx,ndy,nz,nt),
     +            qv(ndx,ndy,nz,nt),uwnd(ndx,ndy,nz,nt), !1u 2v 3t  4oemga 5RH 6HGT
     +            vwnd(ndx,ndy,nz,nt),omega(ndx,ndy,nz,nt), ! !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top 11cpr 12DLR 13ULR,14t2m
     +            theta(ndx,ndy,nz,nt)
      COMMON/SRF/ t2m(ndx,ndy,nt),td2m(ndx,ndy,nt),
     +            u10m(ndx,ndy,nt),v10m(ndx,ndy,nt),
     +            psf(ndx,ndy,nt),qvs(ndx,ndy,nt),
     +            Geo(ndx,ndy)
!      COMMON/TQ/ tvs(10,10,nnr),qvs(10,10,nnr),
!     +  tv(10,10,lev,nnr),qv(10,10,lev,nnr)
      if(ibound.eq.0)then
        ITT=3
      elseif(ibound.eq.-1)then
        ITT=2
      elseif(ibound.eq.-2)then
        ITT=1
      elseif(ibound.eq.1)then
        ITT=2
      elseif(ibound.eq.2)then
        ITT=3
      endif
CC
!------ set the surface loopid  for tibet       
c  if(lo<3)then
c  i1=4
c  i2=5
c  i3=6
c  endif
       OMGO=0.0
       ktop=ikt  ! 100hPa
       DO 200 I=1,NXDIR
       DO 200 J=1,NYDIR
!!!-------------------the follow code get the nearest stand level of the tropopause of every gird!
!         p_top=data2d(I,J,it,10)/100.
!         pres=data2d(I,J,it,9)/100.
          pres=psf(I,J,ITT)/100.     
          do K=1,nz-1
            if(pm(k)>pres.and.pw(k)<pres)then
              isrf=k
              goto 101
            elseif(pw(k)>pres.and.pm(k+1)<pres)then
              isrf=k+1
             goto 101
            endif
          enddo
  101    CONTINUE
!          ktop=ikk
!       if(ktop(i,j)>11)print*,ktop(i,j)
!--------end of ensuring the level of tropopause        
c      OMGO(I,J,1,it)=WSFC(I,J,it)
c      OMGO(I,J,2,it)=OMGO(I,J,1,it)+DIVO(I,J,2,it)
        isrf=1
	      do is=1,iks-1
          OMGO(I,J,is)=0.0
        enddo 
        i1=iks
        i2=i1+1
        i3=i2+1
        OMGO(I,J,i1)=WSFC(I,J)
        OMGO(I,J,i2)=OMGO(I,J,i1)+DIVO(I,J,i2)*(pres-PW(i2))
        IF(pres.GT.PM(i2))
     *      OMGO(I,J,i2)=OMGO(I,J,i1)+DIVO(I,J,i1)*(pres-PW(i2))
        if(isnan(OMGO(I,J,i2)))then
          print*,OMGO(I,J,i2),DIVO(I,J,K)
          print*,pres,PW(i2),1
          stop
        endif
        DO 110 K=i3,ktop-1      !!!! when the k=lev, it is the tropopause, but actually it's not!!!! so need changed!!!
          IF(K.EQ.i3 .AND. pres.LT.PW(i2)) THEN
            OMGO(I,J,K)=OMGO(I,J,i1)+DIVO(I,J,K)*(pres-PW(i3))
            if(isnan(OMGO(I,J,K)))then
            print*,OMGO(I,J,i1),DIVO(I,J,K)
            print*,dp,k
            print*,pres,PW(i3),2
            stop
            endif
          ELSE
            dp=pm(K-1)-pm(k)
            OMGO(I,J,K)=OMGO(I,J,K-1)+DIVO(I,J,K)*dp   !!!50 d of pm
            if(isnan(OMGO(I,J,K)))then
              print*,OMGO(I,J,K-1),DIVO(I,J,K)
              print*,dp,k,3,i,j
              stop
            endif
          ENDIF
  110   CONTINUE
  !-----------above tropopause--------------------------------- 
        do k=ktop,nz-1
!          tv_t=tv(ix,iy, ik, it+1)- tv(ix,iy, ik, it)
               
          tv_p=0.5*(theta(I,J,K+1,ITT)-theta(I,J,K,ITT))
     +              + 0.5*(theta(I,J,K,ITT)-theta(I,J,K-1,ITT))
          dp=pm(K-1)-pm(k)
          OMGO(I,J,K)=-DIVO(I,J,K)
     +              *(tv_p/dp)  
          if(isnan(OMGO(I,J,K)))then
            print*,OMGO(I,J,K)
            print*,dp,k,4
            stop
          endif
        enddo
      IF(pres.LT.PW(i2)) THEN
        OMGO(I,J,i2)=(OMGO(I,J,i1)+OMGO(I,J,i3))*0.5
      endif
  200 CONTINUE
  100 CONTINUE
CC  UNITS:   OMGO(MB/S)
      RETURN
      END SUBROUTINE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  ADJDV ADJUSTS DIVERGENCE AND VERTICAL VELOCITY USING O'BRIEN (1970)
      SUBROUTINE ADJDV(PM,PW,DIVO,OMGO,DIV,NXDIR,NYDIR,ibound,ikt,iks)
      PARAMETER (GRID=2.50)
      integer,parameter :: nz=37 !! add surface
      integer,parameter :: n3d=6 !1u 1v 3t  4oemga 5RH 6HGT
      integer,parameter :: ndx=10
      integer,parameter :: ndy=10
      integer,parameter :: nt=5
      integer NXDIR,NYDIR,ibound
      DIMENSION OMGO(ndx,ndy,nz),DIVO(10,10,nz),DIV(10,10,nz),
     +          PM(nz),PW(nz) 
      COMMON/D3D/ air(ndx,ndy,nz,nt),hgt(ndx,ndy,nz,nt),
     +            qv(ndx,ndy,nz,nt),uwnd(ndx,ndy,nz,nt), !1u 2v 3t  4oemga 5RH 6HGT
     +            vwnd(ndx,ndy,nz,nt),omega(ndx,ndy,nz,nt), ! !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top 11cpr 12DLR 13ULR,14t2m
     +            theta(ndx,ndy,nz,nt)
      COMMON/SRF/ t2m(ndx,ndy,nt),td2m(ndx,ndy,nt),
     +            u10m(ndx,ndy,nt),v10m(ndx,ndy,nt),
     +            psf(ndx,ndy,nt),qvs(ndx,ndy,nt),
     +            Geo(ndx,ndy) ! !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top
CC
      if(ibound.eq.0)then
        ITT=3
      elseif(ibound.eq.-1)then
        ITT=2
      elseif(ibound.eq.-2)then
        ITT=1
      elseif(ibound.eq.1)then
        ITT=2
      elseif(ibound.eq.2)then
        ITT=3
      endif
!      
      DO 410 I=1,NXDIR
      DO 410 J=1,NYDIR
!!!-------------------the follow code get the nearest stand level of the tropopause of every gird!
!   p_top=data2d(I,J,it,10)/100.
!      pres=data2d(I,J,it,9)/100.
        pres=psf(I,J,ITT)/100.
        ktop=ikt  ! 100hPa
!     ktop=12
C     DC=(0.001-OMGO(I,J,18))/(PS(I,J)-PW(18))
        omg_top=OMGO(I,J,ktop)
        omg_top=0
        DC=(omg_top-OMGO(I,J,ktop-1))
     +  /(pres-PW(ktop-1))  !!!top omega
        isfl=1
        do K=1,nz-1
          if(pm(k)>pres.and.pw(k)<pres)then
            isrf=k
            goto 101
          elseif(pw(k)>pres.and.pm(k+1)<pres)then
            isrf=k+1
            goto 101
          endif
        enddo
  101   CONTINUE
c if(lo<3) isfl=4
        isfl=iks
        DO 500 K=isfl,ktop-1
          DIV(I,J,K)=DIVO(I,J,K)+DC
c ccc=-1.0/0
c      if(div(I,j,k,it)==ccc)then
c        print*,DIVO(I,J,K,it),DC,'DC',PW(ktop-1),ktop
c   stop
c endif
  500 CONTINUE
        do k=ktop,nz
          DIV(I,J,k)=DIVO(I,J,k)
        enddo
  410 CONTINUE
      RETURN
      END SUBROUTINE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  Q1Q2  COMPUTES Q1 AND Q2 USING S, Q, U, V, OMG
CC  Note the dx and dy are vector, when the HAD is handle with, carefully note the data store squence
CC  and the array  subscript.
CC  NCEP: EAST-WEAT grids, the longitude array(x+1,y) is great than array(x,y)
CC        South North: the latitude array(x,y+1) is smaller than array(x,y)
      SUBROUTINE Q1Q2(PM,PW,OMG,IX,IY,ibound,ikt,iks)
      PARAMETER (GRID=2.50)
      integer,parameter :: nz=37 !! add surface
      integer,parameter :: n3d=6 !1u 1v 3t  4oemga 5RH 6HGT
      integer,parameter :: ndx=10
      integer,parameter :: ndy=10
      integer,parameter :: nt=5           
      !1u 1v 3oemga 4t 5RH 6HGT
      DIMENSION PM(nz),PW(nz),OMG(10,10,nz)
      DIMENSION STOS(10,10,nz),HADS(10,10,nz),
     +          VADS(10,10,nz)      
      DIMENSION STOQ(10,10,nz),HADQ(10,10,nz),
     +          VADQ(10,10,nz)
      integer NXDIR,NYDIR,it,IX,IY,M,N,ibound
!      DIMENSION WSFC(10),DIVO(10,nz),
!     +          OMGO(10,10),DIV(10,10)
      COMMON/D3D/ air(ndx,ndy,nz,nt),hgt(ndx,ndy,nz,nt),
     +            qv(ndx,ndy,nz,nt),uwnd(ndx,ndy,nz,nt), !1u 2v 3t  4oemga 5RH 6HGT
     +            vwnd(ndx,ndy,nz,nt),omega(ndx,ndy,nz,nt), ! !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top 11cpr 12DLR 13ULR,14t2m
     +            theta(ndx,ndy,nz,nt)
      COMMON/SRF/ t2m(ndx,ndy,nt),td2m(ndx,ndy,nt),
     +            u10m(ndx,ndy,nt),v10m(ndx,ndy,nt),
     +            psf(ndx,ndy,nt),qvs(ndx,ndy,nt),
     +            Geo(ndx,ndy) ! !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top
      COMMON/DXY/ DX(10),DY
!      COMMON/TQ/ tvs(10),qvs(10),
!     +  tv(10,10),qv(10,10)
      COMMON/Q12/ Q1(10,10,NZ),Q2(10,10,NZ),
     +  tls(10,10,NZ),qls(10,10,NZ)
      COMMON/VQ12/VQ1(10,10),VQ2(10,10)
      COMMON/AVEG/HAD_Q(10,10,NZ),VAD_Q(10,10,NZ),
     +  TCH_Q(10,10,NZ),HADT(10,10,NZ),VADT(10,10,NZ),
     +  TCHT(10,10,NZ) 
    
      DT=6.*3600. 
      CP=0.24
      CL=600.
      C1=86400./CP/4.18684
      C2=-0.001*86400.*CL/CP
      TIMC=3600*24.
	Qscale=3600.*24.*2.5e6/1005.
      OMG=OMG*100.   !mb/s to pa/s
!	print*,'dx1=',dx(1)
      if(ibound.eq.0)then
        ITT=3
      elseif(ibound.eq.-1)then
        ITT=2
      elseif(ibound.eq.-2)then
        ITT=1
      elseif(ibound.eq.1)then
        ITT=2
      elseif(ibound.eq.2)then
        ITT=3
      endif
      DO 50 M=1,IX
      DO 50 N=1,IY
!      VQ1(M,N)=0.
!      VQ2(M,N)=0.
      DO 50 K=1,nz  !lev
      STOS(M,N,K)=0.
      HADS(M,N,K)=0.
      VADS(M,N,K)=0.
      STOQ(M,N,K)=0.
      HADQ(M,N,K)=0.
      VADQ(M,N,K)=0.
   50 CONTINUE
      DO 200 M=2,IX-1
      DO 200 N=2,IY-1
!      p_top=data2d(M,N,it,10)/100.
!     pres=data2d(M,N,it,9)/100.
        pres=psf(M,N,ITT)/100.
        ktop=ikt !10  ! 100hPa
!        do iK=1,nz-1
!          if(pm(ik)>p_top.and.pw(ik)<P_top)then
!            ktop=ik-1-ikk
!            goto 101
!          elseif(pw(ik)>P_top.and.pm(ik+1)<p_top)then
!            ktop=ik-ikk
!             goto 101
!          endif
!        enddo
!  101 CONTINUE
!
!      ktop=ikk
! the calculate above troposphause  ????
!  psf=data2d(M,N,IT,9)/100.
! psf=600.0  !!!!!!!!!!!!!!!!!!!!!!!!!!   for WTP
! if(psf>600.)print*,psf
        isfl=1
        do K=1,nz-1
          if(pm(k)>pres.and.pw(k)<pres)then
            isrf=k
            goto 101
          elseif(pw(k)>pres.and.pm(k+1)<pres)then
            isrf=k+1
            goto 101
          endif
        enddo
  101   CONTINUE
      DO 150 K=2,ktop-1
!    psf=data2d(M,N,IT,9)/100.
        STOS(M,N,K)=0.
        STOQ(M,N,K)=0.
        HADS(M,N,K)=0.
        VADS(M,N,K)=0.
        HADQ(M,N,K)=0.
        VADQ(M,N,K)=0.
!      IF(K.EQ.2 .AND. psf.LT.PM(K))then
!  print*,k,psf,pm(k)
!  GO TO 900   !!!find out why goto 900
!      endif
        IF(K.LT.14) GO TO 900 ! the second level convert to zero
!        IF(K.EQ.2 .AND. psf.LT.PM(K)) GO TO 900 ! the second level is not always zero
!        IF(K.EQ.3 .AND. psf.LT.PM(K)) GO TO 900
!        IF(K.EQ.4 .AND. psf.LT.PM(K)) GO TO 900
!        IF(K.EQ.5 .AND. psf.LT.PM(K)) GO TO 900
!        DO KK =1,isrf-1
!          IF(K.EQ.KK .AND. PRES .LT. PM(K))GO TO 900
!        ENDDO
!!!Time boundary----------------------------
!=================  local tendency=====================================
        IF(ibound.EQ.-2)then
!        STOS(M,N,K)=(data3d(M,N,K,IT+2,3)-data3d(M,N,K,IT,3))/(2.*DT) ! 
!        STOQ(M,N,K)=(QV(M,N,K,IT+2)-QV(M,N,K,IT))/(2.*DT)
          STOS(M,N,K)=(air(M,N,K,ITT+2)-air(M,N,K,ITT))/(2.*DT) ! 
          STOQ(M,N,K)=(QV(M,N,K,ITT+2)-QV(M,N,K,ITT))/(2.*DT)        
!!!!!!---------Time boundary
        elseif(ibound.eq.2)then
!            STOS(M,N,K)=(data3d(M,N,K,IT,3)-data3d(M,N,K,IT-2,3))/(2.*DT) ! SN= next time tv; SP=previous time
!            STOQ(M,N,K)=(QV(M,N,K,IT)-QV(M,N,K,IT-2))/(2.*DT)
          STOS(M,N,K)=(air(M,N,K,ITT)-air(M,N,K,ITT-2))/(2.*DT) ! 
          STOQ(M,N,K)=(QV(M,N,K,ITT)-QV(M,N,K,ITT-2))/(2.*DT) 
        else
!            STOS(M,N,K)=(data3d(M,N,K,3,3)-data3d(M,N,K,1,3))/(2.*DT) ! SN= next time tv; SP=previous time
!            STOQ(M,N,K)=(QV(M,N,K,3)-QV(M,N,K,1))/(2.*DT)
          STOS(M,N,K)=(air(M,N,K,ITT+1)-air(M,N,K,ITT-1))/(2.*DT) ! 
          STOQ(M,N,K)=(QV(M,N,K,ITT+1)-QV(M,N,K,ITT-1))/(2.*DT) 
        ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(M.EQ.2) THEN
C     IF(MM.EQ.1 .OR. PS(M-2,N).LT.PM(K)) THEN
          MM=M+2
          IF((M+2)>IX)then
            AXS=uwnd(M,N,K,ITT)*
     +          (air(M+1,N,K,ITT)-air(M,N,K,ITT))/(1.*DX(N))  !Temp.
            AXQ=uwnd(M,N,K,ITT)*
     +          (QV(M,N,K,ITT)-QV(M,N,K,ITT))/(1.*DX(N))
          ELSE
            AXS=0.5*UWND(M+1,N,K,ITT)*
     +          (AIR(M+2,N,K,ITT)-AIR(M,N,K,ITT))/(2.*DX(N))  !Temp.
            AXQ=0.5*UWND(M+1,N,K,ITT)*
     +          (QV(M+2,N,K,ITT)-QV(M,N,K,ITT))/(2.*DX(N))
          ENDIF
        ELSEIF(M.EQ.(IX-1)) THEN
C      ELSEIF(MM.EQ.7 .OR. PS(M+2,N).LT.PM(K)) THEN
!      AXS=0.5*U(M-1,N,K)*(S(M,N,K)-S(M-2,N,K))/(2.*DX(N))
          AXS=0.5*UWND(M-1,N,K,ITT)*
     +        (AIR(M,N,K,ITT)-AIR(M-2,N,K,ITT))/(2.*DX(N)) 
!      AXQ=0.5*U(M-1,N,K)*(Q(M,N,K)-Q(M-2,N,K))/(2.*DX(N))
          AXQ=0.5*UWND(M-1,N,K,ITT)*
     +        (QV(M,N,K,ITT)-QV(M-2,N,K,ITT))/(2.*DX(N))
    
        ELSE
          IF((M+2)>IX)then
            AXS=0.5*( UWND(M,N,K,ITT)*(AIR(M+1,N,K,ITT)
     *          -AIR(M,N,K,ITT)) +UWND(M-2,N,K,ITT)*
     +          (AIR(M,N,K,ITT)-AIR(M-1,N,K,ITT)) )/(DX(N))
            AXQ=0.5*(UWND(M,N,K,ITT)*(QV(M+1,N,K,ITT)-QV(M,N,K,ITT))
     *          +    UWND(M-2,N,K,ITT)*(QV(M,N,K,ITT)-QV(M-1,N,K,ITT))
     *               )/(DX(N))
          ELSE
            AXS=0.5*( UWND(M+1,N,K,ITT)*(AIR(M+2,N,K,ITT)
     *            -AIR(M,N,K,ITT)) +UWND(M-1,N,K,ITT)*
     +             (AIR(M,N,K,ITT)-AIR(M-2,N,K,ITT)))/(2.*DX(N))
            AXQ=0.5*(UWND(M+1,N,K,ITT)*(QV(M+2,N,K,ITT)-QV(M,N,K,ITT))
     *          +    UWND(M-1,N,K,ITT)*(QV(M,N,K,ITT)-QV(M-2,N,K,ITT))
     *               )/(2.*DX(N))
          ENDIF
        ENDIF
        IF(N.EQ.2) THEN
C     IF(NN.EQ.1 .OR. PS(M,N-2).LT.PM(K)) THEN
!      AYS=0.5*V(M,N+1,K)*(S(M,N+2,K)-S(M,N,K))/(2.*DY)
          NN=N+2
          DYY=0
          IF((N+2)>IY)then
            AYS= -VWND(M,N,K,ITT)*(AIR(M,N+1,K,ITT)  !!! why have a '-', for N+1 and N stand the actually
     +          -AIR(M,N,K,ITT))/(1.*DY)
            AYQ= -VWND(M,N,K,ITT)*(QV(M,N+1,K,ITT)
     +          -QV(M,N,K,ITT))/(1.*DY)
          ELSE
            AYS=-0.5*VWND(M,N+1,K,ITT)*(AIR(M,N+2,K,ITT)
     +          -AIR(M,N,K,ITT))/(2.*DY)
            AYQ=-0.5*VWND(M,N+1,K,ITT)*(QV(M,N+2,K,ITT)
     +          -QV(M,N,K,ITT))/(2.*DY)
          ENDIF
!      if(abs(AYS*TIMC)>30)then
! print*,'1XXXXX'
!      print*, AYS*TIMC
!      print*, data3d(M,N+2,K,IT,3)-data3d(M,N,K,IT,3)
! print*, data3d(M,N+1,K,IT,2)
! print*,M,N,K,'lo',lo
! stop
! endif
!      AYQ=0.5*V(M,N+1,K)*(Q(M,N+2,K)-Q(M,N,K))/(2.*DY)
        ELSEIF(N.EQ.(IY-1).and.(N-2)>0) THEN
C      ELSEIF(NN.EQ.7 .OR. PS(M,N+2).LT.PM(K)) THEN
!      AYS=0.5*V(M,N-1,K)*(S(M,N,K)-S(M,N-2,K))/(2.*DY)
          AYS=-0.5*VWND(M,N-1,K,ITT)*(AIR(M,N,K,ITT)
     +        -AIR(M,N-2,K,ITT))/(2.*DY)
! if(abs(AYS*TIMC)>30)then
! print*,'2XXXXX'
!      print*, AYS*TIMC
!      print*, data3d(M,N,K,IT,3)-data3d(M,N-2,K,IT,3)
! print*, data3d(M,N-1,K,IT,2)
! print*,M,N,K,'lo',lo
! stop
! endif
!      AYQ=0.5*V(M,N-1,K)*(Q(M,N,K)-Q(M,N-2,K))/(2.*DY)
          AYQ=-0.5*VWND(M,N-1,K,ITT)*(QV(M,N,K,ITT)
     +        -QV(M,N-2,K,ITT))/(2.*DY)
        ELSE
!      AYS=0.5*(V(M,N+1,K)*(S(M,N+2,K)-S(M,N,K))
!     *        +V(M,N-1,K)*(S(M,N,K)-S(M,N-2,K)))/(2.*DY)
!      AYQ=0.5*(V(M,N+1,K)*(Q(M,N+2,K)-Q(M,N,K))
!     *        +V(M,N-1,K)*(Q(M,N,K)-Q(M,N-2,K)))/(2.*DY)
          IF((N+2)>IY)then
            AYS=0.5*( VWND(M,N,K,ITT)*(-AIR(M,N+1,K,ITT) + !-
     +            AIR(M,N,K,ITT))+ VWND(M,N-1,K,ITT)*(-AIR(M,N,K,ITT)+   ! - 
     +            AIR(M,N-1,K,ITT)) )/(1.*DY)
            AYQ=0.5*( VWND(M,N,K,ITT)*(-QV(M,N+1,K,ITT)+QV(M,N,K,ITT))
     +            + VWND(M,N-1,K,ITT)*(-QV(M,N,K,ITT)+QV(M,N-1,K,ITT))
     +                )/(1.*DY)
          ELSE
            AYS=0.5*( VWND(M,N+1,K,ITT)*(-AIR(M,N+2,K,ITT) + !-
     +          AIR(M,N,K,ITT))+VWND(M,N-1,K,ITT)*(-AIR(M,N,K,ITT)+  !-
     +          AIR(M,N-2,K,ITT)) )/(2.*DY)
            AYQ=0.5*(VWND(M,N+1,K,ITT)*(-QV(M,N+2,K,ITT)+QV(M,N,K,ITT))
     +             + VWND(M,N-1,K,ITT)*(-QV(M,N,K,ITT)+QV(M,N-2,K,ITT))
     +               )/(2.*DY)
          ENDIF
        ENDIF
        HADS(M,N,K)=AXS+AYS
        if(isnan(HADS(M,N,K)))then
          print*,AXS,AYS
          print*,M,N,K
          stop
        endif
! if(abs(HADS(M,N,K)*TIMC)>30)then
!      print*, HADS(M,N,K)*TIMC,AXS*TIMC,AYS*TIMC
!      print*,M,N,K,'lo',lo
! print*,'*******************************************'
! stop
! endif
        HADQ(M,N,K)=AXQ+AYQ
        UKP1=.5*(THETA(M,N,K+1,ITT)-THETA(M,N,K,ITT))*OMG(M,N,K)
        UKM1=.5*(THETA(M,N,K,ITT)-THETA(M,N,K-1,ITT))*OMG(M,N,K-1)
        VKP1=.5*(QV(M,N,K+1,ITT)-QV(M,N,K,ITT))*OMG(M,N,K)
        VKM1=.5*(QV(M,N,K,ITT)-QV(M,N,K-1,ITT))*OMG(M,N,K-1)
        DPW=Pm(K)-Pm(K-1)
        kss=14
	  IF(ABS(VKP1)>7000.)THEN
	  PRINT*,VKP1,K
        PRINT*,QV(M,N,K+1,ITT),QV(M,N,K,ITT),OMG(M,N,K)
	  PRINT*,itt,QV(M,N,K,ITT+1)
	  PAUSE
	 ENDIF
!        if(lo<3)kss=4
       IF(K.GE.kss) GO TO 850   !!! k=5 in hmbudps is 800hPa this sentence is think that surface pressure great than 800hPa
!-------in my calculate, when k=15, the pressure=550, k=16,p=500, so changge to K GE 15     
!       IF(K.EQ.4 .AND. psf.GE.PM(3)) GO TO 850
!         IF(K.EQ.4 .AND. psf.LT.PM(3)) THEN
!          UKM1=.5*(S(M,N,K)-S(M,N,1))*OMG(M,N,K-1)
!          VKM1=.5*(Q(M,N,K)-Q(M,N,1))*OMG(M,N,K-1)
!         ENDIF
        IF(K.EQ.kss.AND. PRES.GE.PM(kss-1)) GO TO 850
        IF(K.EQ.kss.AND. PRES.LT.PW(kss-1)) THEN
          UKM1=.5*(THETA(M,N,K,ITT)-THETA(M,N,K-1,ITT))*OMG(M,N,1)
          VKM1=.5*(QV(M,N,K,ITT)-QV(M,N,K-1,ITT))*OMG(M,N,1)
          DPW=PW(K)-PRES
        ENDIF
        IF(K.EQ.kss.AND. PRES.GE.PW(kss-1)) THEN
          UKM1=.5*(THETA(M,N,K,ITT)-THETA(M,N,K-1,ITT))*OMG(M,N,K-1)
          VKM1=.5*(QV(M,N,K,ITT)-QV(M,N,K-1,ITT))*OMG(M,N,K-1)
        ENDIF
        IF(K.EQ.kss.AND. PRES.GE.PM(kss-1)) THEN
           if(kss==3)DPW=PW(K)-pres
           if(kss>3)DPW=PW(K)-PW(K-1)
        ENDIF
  850   VADS(M,N,K)=((pm(k)/1000)**0.286)*(UKM1+UKP1)/(DPW*100.) !(UKP1+UKM1)/DPW 
        VADQ(M,N,K)=(VKP1+VKM1)/(DPW*100.)
!      VADQ(M,N,K)=(VKP1+VKP1)/(DPW*100)
        IF(ABS(VADQ(M,N,K))>1)THEN
	        PRINT*,VADQ(M,N,K),pm(k)
          PRINT*,VKP1,VKM1,DPW
          PRINT*,DY,DX(1),DX(2),DX(3)
	       PAUSE
	      ENDIF
  900 CONTINUE
  150 CONTINUE
  200 CONTINUE
      open(999,file='F999.txt',status='unknown',position='APPEND')
      DO 11 I=1,IX
      DO 11 J=1,IY
      DO 11 K=1,NZ
      tls(I, J, K)=  -HADS(I,J,K) -VADS(I,J,K) 
      qls(I, J, K)=  -HADQ(I,J,K) -VADQ(I,J,K)
      Q1(I,J,K)=(STOS(I,J,K)+HADS(I,J,K)+VADS(I,J,K))
      Q2(I,J,K)=-(STOQ(I,J,K)+HADQ(I,J,K)+VADQ(I,J,K))
	IF(Q2(I,J,K)*Qscale.GT. 10000.)THEN
		PRINT*,Q2(I,J,K)
		PRINT*,STOQ(I,J,K), HADQ(I,J,K),VADQ(I,J,K)
	    PAUSE
	ENDIF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      HAD_Q(I,J,K)=-HADQ(I,J,K)
      VAD_Q(I,J,K)=-VADQ(I,J,K)
      TCH_Q(I,J,K)=-STOQ(I,J,K)
      HADT(I,J,K)=HADS(I,J,K)
      VADT(I,J,K)=VADS(I,J,K) 
      TCHT(I,J,K)=STOS(I,J,K)
C     Q1(I,J,K)=(HADS(I,J,K)+VADS(I,J,K))*C1
C     Q2(I,J,K)=(HADQ(I,J,K)+VADQ(I,J,K))*C2
   11 CONTINUE
CCC   UNIT:   DEG/DAY
      C1=86400./CP/4.18684
      C2=-0.001*86400.*CL/CP
	    VQ1=0.0
	    VQ2=0.0
      DO 20 I=1,IX
      DO 20 J=1,IY
      DO 30 K=2,ikt-1
       DPW=PW(K-1)-PW(K)
C     VQ1(I,J)=VQ1(I,J)+Q1(I,J,K)*DPW/(PW(1)-PW(18))
C     VQ2(I,J)=VQ2(I,J)+Q2(I,J,K)*DPW/(PW(1)-PW(18))
       pres=psf(I,J,ITT)/100.
       IF(K.EQ.2) DPW=Pres-PW(K)
C     VQ1(I,J)=VQ1(I,J)+Q1(I,J,K)*DPW*CP/CL/(9.8*24.)*100.
C     VQ2(I,J)=VQ2(I,J)+Q2(I,J,K)*DPW*CP/CL/(9.8*24.)*100.
C    UNIT:   MM/H
C     VQ1(I,J)=VQ1(I,J)+Q1(I,J,K)*DPW/(PS(I,J)-PW(18))
C     VQ2(I,J)=VQ2(I,J)+Q2(I,J,K)*DPW/(PS(I,J)-PW(18))
      Q1tmp=Q1(I,J,K)*24.*3600.
      Q2tmp=Q2(I,J,K)*24.*3600.*2.5e6/1005.
      VQ1(I,J)=VQ1(I,J)+Q1tmp*DPW*1004./9.8/86400.
      VQ2(I,J)=VQ2(I,J)+Q2tmp*DPW*1004./9.8/86400.
C    UNIT: 100W/M**2
   40 CONTINUE
   30 CONTINUE
   20 CONTINUE
      RETURN
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE ERA_SAS(PM,DP,OMG,IX,IY,ibound,iks,ike,XY_OUT)
      integer,parameter :: nz=37 !! add surface
      integer,parameter :: n3d=6 !1u 1v 3t  4oemga 5RH 6HGT
      integer,parameter :: ndx=10
      integer,parameter :: ndy=10
      integer,parameter :: nt=5           
      !1u 1v 3oemga 4t 5RH 6HGT
!      DIMENSION WSFC(10),DIVO(10,nz),
!     +          OMGO(10,10),DIV(10,10)
      COMMON/D3D/ air(ndx,ndy,nz,nt),hgt(ndx,ndy,nz,nt),
     +            qv(ndx,ndy,nz,nt),uwnd(ndx,ndy,nz,nt), !1u 2v 3t  4oemga 5RH 6HGT
     +            vwnd(ndx,ndy,nz,nt),omega(ndx,ndy,nz,nt), ! !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top 11cpr 12DLR 13ULR,14t2m
     +            theta(ndx,ndy,nz,nt)
      COMMON/SRF/ t2m(ndx,ndy,nt),td2m(ndx,ndy,nt),
     +            u10m(ndx,ndy,nt),v10m(ndx,ndy,nt),
     +            psf(ndx,ndy,nt),qvs(ndx,ndy,nt),
     +            Geo(ndx,ndy) ! !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top
      COMMON/AVEG/HAD_Q(10,10,nz),VAD_Q(10,10,nz),
     +  TCH_Q(10,10,nz),HADT(10,10,nz),VADT(10,10,nz),
     +  TCHT(10,10,nz) 
!      COMMON/TQ/ tvs(10),qvs(10),
!     +  tv(10,10),qv(10,10)
      COMMON/Q12/ Q1(10,10,NZ),Q2(10,10,NZ),
     +  tls(10,10,NZ),qls(10,10,NZ)
      DIMENSION PM(nz),DP(nz),OMG(10,10,nz)
      integer IX,IY,M,N,ibound
      real temp(20),XY_OUT(18)
!
!  
      if(ibound.eq.0)then
        ITT=3
      elseif(ibound.eq.-1)then
        ITT=2
      elseif(ibound.eq.-2)then
        ITT=1
      elseif(ibound.eq.1)then
        ITT=2
      elseif(ibound.eq.2)then
        ITT=3
      endif   
      temp=0.0 
      XY_OUT=0.0
      hlat=2.5e6
      cp=1005.
      Scale1=3600.*24.
      Scale2=3600.*24.*hlat/cp 
!---------low atomosphere
      SK=0.
      DO 198 ik=iks,ike
        DO 197 ix=2,IX-1
          DO 197 iy=2,IY-1
            temp(1)= temp(1)+tls(ix,iy, ik)*dp(ik)*Scale1 !day-1
            temp(2)= temp(2)+qls(ix,iy, ik)*dp(ik)*Scale2 !!! K/day day-1
            temp(3)= temp(3)+uwnd(ix,iy,ik,ITT)*dp(ik)
            temp(4)= temp(4)+vwnd(ix,iy,ik,ITT)*dp(ik)
            temp(5)= temp(5)+qv(ix,iy, ik, ITT)*dp(ik)
            temp(6)= temp(6)+hgt(ix,iy,ik,ITT)*dp(ik)/9.87  !HGT
            temp(7)= temp(7)+air(ix,iy,ik,ITT)*dp(ik)
            temp(8)= temp(8)+omg(ix,iy, ik)*dp(ik)
!            temp(9)= temp(9)+DATA3D(ix,iy,ik,ITT,5)  
!           temp(10)= temp(10)+qv(ix,iy, ik, it)
            temp(9)= temp(9)+theta(ix,iy, ik, ITT)*dp(ik)
            temp(10)= temp(10)+Q1(ix, iy, ik)*dp(ik)*Scale1   !Tmcl
            temp(11)= temp(11)+Q2(ix, iy, ik) *dp(ik)*Scale2   !*TMCL*hlat/cp     !k/day
            temp(12)= temp(12)+HAD_Q(ix,iy, ik)*dp(ik)*Scale2   ! day-1
            temp(13)= temp(13)+VAD_Q(ix,iy, ik)*dp(ik)*Scale2   !day-1
            temp(14)= temp(14)+TCH_Q(ix,iy, ik)*dp(ik)*Scale2  !!! K/day day-1
            temp(15)= temp(15)+HADT(ix,iy, ik)*dp(ik)*Scale1 !day-1
            temp(16)= temp(16)+VADT(ix,iy, ik)*dp(ik)*Scale1 !!! K/day day-1
            temp(17)= temp(17)+TCHT(ix,iy, ik)*dp(ik)*Scale1 !! K/day day-1
            temp(18)= temp(18)+omega(ix,iy,ik,ITT)*dp(ik)
197   continue
        SK=DP(ik)+SK
198   continue
      XY=(IX-2.0)*(IY-2.)
!      SK=PW(iks)-PW(ike)
      XYS=XY*SK
      do ii=1,18
        XY_OUT(ii)=temp(ii)/XYS
      enddo
!
      return
      end SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      SUBROUTINE Read3d(fpath,nd, X1, X2, Y1, Y2, NX, NY,tmp3d)
!                    call Read3d(fpath, ndd, X1(ig), X2(ig),
!     +           Y1(ig), Y2(ig), IXX(ig), IYY(ig), DAT3D,ih,iv) !REAL DAT3D(6,ndx,ndy,nz,4,nday)
      include 'netcdf.inc'
!
      integer,parameter :: nz=37 !! add surface
      integer,parameter :: ndx=10
      integer,parameter :: ndy=10
      integer,parameter :: nlon=144
      integer,parameter :: nlat=73
      integer,parameter :: nday=366
!      
      integer*4  ncid, status    ! file control
!-------------------------------------------------------------
 !     Below 4 variables is the data in netCDF file
      real*4       longitude( 144 )
      real*4       latitude( 73 )
	integer*4    level(37)
      integer*4, allocatable ::  time(:)
      integer*2, allocatable ::  raw( :, :, :,: )
!      real*8 , allocatable ::  time(:)
!      integer*2, allocatable  ::  f3d(:, :,:,: )
!     above4 variables is the data in netCDF file
!-------------------------------------------------------------
      integer*4   :: start(10)
      integer*4   :: count(10)
      integer*4   :: dimids(10)! allow up to 10 dimensions
      integer*4   :: dimid, xtype
      character(len=31) :: dummy
      real*8     ::  scale5(1), add5(1)
!----------------------------------------------------------
      character(len=150) :: fpath     
      integer X1, X2, Y1, Y2
      integer NX, NY
      real tmp3d(ndx,ndy,nz,nday)

!----------------------------------------------
      allocate(time(nd)) 
!  allocate(temp3d(144, 73, 17, nrec ))
      allocate(raw(nlon, nlat, nz, nd )) 
!
!      print*,nd
!	print*,trim(fpath)
      status=nf_open(trim(fpath),nf_nowrite,ncid)
      if ( status/=nf_noerr )	 write (*,*) nf_strerror(status)
!
	if (status/=nf_noerr)then
!	  write (*,*) nf_strerror(status)
	  print*,trim(fpath),'AAA'
	  stop
      endif
!   Retrieve data for Variable 'longitude'
!   Units of 'longitude' is 'degrees_east'
!   Long_name of 'longitude' is 'longitude'
      status=nf_inq_var(ncid,   1,dummy,xtype,ndim,dimids,natts)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
        status=nf_inq_dim(ncid,dimids(j),dummy,len)
        if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
        start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   1,start,count,longitude)
!   Retrieve data for Variable 'latitude'
!   Units of 'latitude' is 'degrees_north'
!   Long_name of 'latitude' is 'latitude'
      status=nf_inq_var(ncid,   2,dummy,xtype,ndim,dimids,natts)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   2,start,count,latitude)
!   Retrieve data for Variable 'level'
!   Units of 'level' is 'millibars'
!   Long_name of 'level' is 'pressure_level'
      status=nf_inq_var(ncid,   3,dummy,xtype,ndim,dimids,natts)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
        status=nf_inq_dim(ncid,dimids(j),dummy,len)
        if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
        start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_int(ncid,   3,start,count,level)
!   Retrieve data for Variable 'time'
!   Units of 'time' is 'hours since 1900-01-01 00:00:0.0'
!   Long_name of 'time' is 'time'
      status=nf_inq_var(ncid,   4,dummy,xtype,ndim,dimids,natts)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
        status=nf_inq_dim(ncid,dimids(j),dummy,len)
        if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
        start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_int(ncid,   4,start,count,time)
!   Retrieve data for Variable 
!   Units of 
!         
      status=nf_inq_var(ncid,   5,dummy,xtype,ndim,dimids,natts)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
        status=nf_inq_dim(ncid,dimids(j),dummy,len)
        if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
        start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_int2(ncid,   5,start,count,raw)
      scale5(1) =0.0 ; add5(1) =0.0
      status=nf_get_att_double(ncid,   5,'scale_factor',scale5(1))
      status=nf_get_att_double(ncid,   5,'add_offset',add5(1))
!------------------------------------------------------------------------
!-------  get data ------------------------------------------------------
!      open(901,file='901.txt')
!	write(901,*) nd,X1,X2,Y1,Y2,fpath
!      pause
      print*, latitude(1),latitude(73),'Checking latitude'
	do id=1,nd
        do ix=X1,X2 ! weat to east
          do iy=Y1,Y2,-1 ! south to north
            do ik=nz,1,-1
              inx=ix-X1+1
              iny=iy-Y2+1
              kk=nz-ik+1  !!!! in ERA dataset, when k=1,it's the top of atmossphere, and k=nz is the 1000hpa
              tmp3d(inx,iny,kk,id)=
     +            raw(ix,iy,ik,id )*scale5(1)+add5(1)
!	write(901,*)id,ix,iy,ik
!      write(901,*)inx,iny,kk
! 	write(901,*)tmp3d(inx,iny,kk,id)
!      pause
            enddo
          enddo
        enddo
      enddo  
!--------- close the netcdf file ------------------------------
      status=nf_close(ncid)
      return
      END SUBROUTINE
!#############################################################################################################
      SUBROUTINE ReadSRF(fpath,ih,ids,nd,X1,X2,Y1,Y2,NX,NY,tmpsrf)              
      include 'netcdf.inc'
!
      integer,parameter :: ndx=10
      integer,parameter :: ndy=10
      integer,parameter :: nlon=144
      integer,parameter :: nlat=73
!      
      integer*4  ncid, status    ! file control
!-------------------------------------------------------------
 !     Below 4 variables is the data in netCDF file
      real*4       longitude( 144 )
      real*4       latitude( 73 )
      integer*4    time(13149)
      integer*2    raw( 144, 73, 13149)
!      real*8 , allocatable ::  time(:)
!      integer*2, allocatable  ::  f3d(:, :,:,: )
!     above4 variables is the data in netCDF file
!-------------------------------------------------------------
      integer*4   :: start(10)
      integer*4   :: count(10)
      integer*4   :: dimids(10)! allow up to 10 dimensions
      integer*4   :: dimid, xtype
      character(len=31) :: dummy
      real*8     ::  scale4(1), add4(1)
!----------------------------------------------------------
      character(len=150) :: fpath     
      integer X1, X2, Y1, Y2, ih, ids, nd
      integer NX, NY
      real tmpsrf(ndx,ndy,366)
!----------------------------------------------
!
      status=nf_open(trim(fpath),nf_nowrite,ncid)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
	if ( ncid.lt.0 )then
!	  write (*,*) nf_strerror(status)
	  print*,trim(fpath)
	  stop
      endif
!   Retrieve data for Variable 'longitude'
!   Units of 'longitude' is 'degrees_east'
!   Long_name of 'longitude' is 'longitude'
      status=nf_inq_var(ncid,   1,dummy,xtype,ndim,dimids,natts)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
        status=nf_inq_dim(ncid,dimids(j),dummy,len)
        if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
        start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   1,start,count,longitude)
!   Retrieve data for Variable 'latitude'
!   Units of 'latitude' is 'degrees_north'
!   Long_name of 'latitude' is 'latitude'
      status=nf_inq_var(ncid,   2,dummy,xtype,ndim,dimids,natts)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   2,start,count,latitude)
!   Retrieve data for Variable 'time'
!   Units of 'time' is 'hours since 1900-01-01 00:00:0.0'
!   Long_name of 'time' is 'time'
      status=nf_inq_var(ncid,   3,dummy,xtype,ndim,dimids,natts)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
        status=nf_inq_dim(ncid,dimids(j),dummy,len)
        if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
        start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_int(ncid,   3,start,count,time)
!   Retrieve data for Variable 
!   Units of 
!         
      status=nf_inq_var(ncid,   4,dummy,xtype,ndim,dimids,natts)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
        status=nf_inq_dim(ncid,dimids(j),dummy,len)
        if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
        start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_int2(ncid,   4,start,count,raw)
      scale4(1) =0.0 ; add4(1) =0.0
      status=nf_get_att_double(ncid,   4,'scale_factor',scale4(1))
      status=nf_get_att_double(ncid,   4,'add_offset',add4(1))
!------------------------------------------------------------------------
!-------  get data ------------------------------------------------------
      itt=1
      do it=ids,ids+nd-1
        do ix=X1,X2 ! weat to east
          do iy=Y1,Y2,-1 ! south to north
!            do ik=nz,1,-1
              inx=ix-X1+1
              iny=iy-Y2+1
              tmpsrf(inx,iny,itt)=
     +            raw(ix,iy,it)*scale4(1)+add4(1)
!       if(i3d==4)then
!         print*,data3d(ix-X1(lo)+1,iy-Y2(lo)+1,ik,it+itt,i3d)
!               endif
!            enddo
          enddo
        enddo
        itt=itt+1
      enddo  
!--------- close the netcdf file ------------------------------
      status=nf_close(ncid)
      return
      END SUBROUTINE
!#############################################################################################################
      SUBROUTINE ERAGEO(fpath,X1, X2, Y1, Y2, NX, NY, tmpgeo)              
      include 'netcdf.inc'
!
      integer,parameter :: ndx=10
      integer,parameter :: ndy=10
      integer,parameter :: nlon=144
      integer,parameter :: nlat=73
!      
      integer*4  ncid, status    ! file control
!-------------------------------------------------------------
 !     Below 4 variables is the data in netCDF file
      real*4       longitude( 144 )
      real*4       latitude( 73 )
      integer*4    time(1)
      integer*2    raw( 144, 73, 1)
!      real*8 , allocatable ::  time(:)
!      integer*2, allocatable  ::  f3d(:, :,:,: )
!     above4 variables is the data in netCDF file
!-------------------------------------------------------------
      integer*4   :: start(10)
      integer*4   :: count(10)
      integer*4   :: dimids(10)! allow up to 10 dimensions
      integer*4   :: dimid, xtype
      character(len=31) :: dummy
      real*8     ::  scale4(1), add4(1)
!----------------------------------------------------------
      character(len=150) :: fpath     
      integer X1, X2, Y1, Y2
      integer NX, NY
      real tmpgeo(ndx,ndy)
!----------------------------------------------
!
      status=nf_open(trim(fpath),nf_nowrite,ncid)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
	if ( ncid.lt.0 )then
!	  write (*,*) nf_strerror(status)
	  print*,trim(fpath)
	  stop
      endif
!   Retrieve data for Variable 'longitude'
!   Units of 'longitude' is 'degrees_east'
!   Long_name of 'longitude' is 'longitude'
      status=nf_inq_var(ncid,   1,dummy,xtype,ndim,dimids,natts)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
        status=nf_inq_dim(ncid,dimids(j),dummy,len)
        if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
        start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   1,start,count,longitude)
!   Retrieve data for Variable 'latitude'
!   Units of 'latitude' is 'degrees_north'
!   Long_name of 'latitude' is 'latitude'
      status=nf_inq_var(ncid,   2,dummy,xtype,ndim,dimids,natts)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   2,start,count,latitude)
!   Retrieve data for Variable 'time'
!   Units of 'time' is 'hours since 1900-01-01 00:00:0.0'
!   Long_name of 'time' is 'time'
      status=nf_inq_var(ncid,   3,dummy,xtype,ndim,dimids,natts)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
        status=nf_inq_dim(ncid,dimids(j),dummy,len)
        if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
        start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_int(ncid,   3,start,count,time)
!   Retrieve data for Variable 
!   Units of 
!         
      status=nf_inq_var(ncid,   4,dummy,xtype,ndim,dimids,natts)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
        status=nf_inq_dim(ncid,dimids(j),dummy,len)
        if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
        start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_int2(ncid,   4,start,count,raw)
      scale4(1) =0.0 ; add4(1) =0.0
      status=nf_get_att_double(ncid,   4,'scale_factor',scale4(1))
      status=nf_get_att_double(ncid,   4,'add_offset',add4(1))
!------------------------------------------------------------------------
!-------  get data ------------------------------------------------------
        do ix=X1,X2 ! weat to east
          do iy=Y1,Y2,-1 ! south to north
!            do ik=nz,1,-1
              inx=ix-X1+1
              iny=iy-Y2+1
              tmpgeo(inx,iny)=
     +            (raw(ix,iy,1)*scale4(1)+add4(1))/9.8
!       if(i3d==4)then
!         print*,data3d(ix-X1(lo)+1,iy-Y2(lo)+1,ik,it+itt,i3d)
!               endif
!            enddo
          enddo
        enddo
!--------- close the netcdf file ------------------------------
      status=nf_close(ncid)
      return
      END SUBROUTINE