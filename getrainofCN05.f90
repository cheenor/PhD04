program GetRain
implicit none
integer, parameter:: NX=142, NY=82
real, parameter:: blat=14.75,blon=69.75,intr=0.5
character path*100,nm(2)*80,fpath*80
integer iy,im,id
real rain 
integer i,k,iyy,it,day(12),j,nd
real raw(2,50,12,31,NX,NY)
real tlon,tlat,tvalue
integer yy(6),mm(6),dd(6),ndd(6)            !!! 1-6 ETP,WTP,PRD,MLYR,NPC,NECl
character rgns(6)*4,dir*50
!------- path strings--------------
dir='Z:\DATA\CN05\CN05.2\'
nm(2)='2400_China_Pre_1961_2010_Full_daily_05x05.dat'
nm(1)='2400_China_Tm_1961_2010_Full_daily_05x05.dat'
!---------------------------------------------------------
!!!! ETP
yy(1)=2010; mm(1)=6; dd(1)=1; ndd(1)=30; rgns(1)='ETP'
!!!! WTP
yy(2)=2010; mm(2)=6; dd(2)=1; ndd(2)=30; rgns(2)='WTP'
!!!! PRD
yy(3)=2010; mm(3)=6; dd(3)=1; ndd(3)=30; rgns(3)='PRD'
!!!! MLYR
yy(4)=2010; mm(4)=6; dd(4)=1; ndd(4)=30; rgns(4)='MLYR'
!!!! NPC
yy(5)=2010; mm(5)=8; dd(5)=2; ndd(5)=30; rgns(5)='NPC'
!!!! NEC
yy(6)=2010; mm(6)=6; dd(6)=1; ndd(6)=30; rgns(6)='NEC'

do i=1,12
   day(i)=30
enddo
   day(1)=31;day(3)=31;day(5)=31;day(7)=31
   day(8)=31;day(10)=31;day(12)=31
do 333 k=1,2   
        open(12+k,file=trim(dir)//trim(nm(k)),form='binary') 
        do 333 iyy=1961,2010
                iy=iyy-1960
                day(2)=28
                if(mod(iyy,4)==0.and.mod(iyy,100)/=0)then
                day(2)=29
                elseif(mod(iyy,400)==0)then
                day(2)=29
                endif
                do 333 im=1,12
                  nd=day(im)
                    do 333 id=1,nd
                    do 333 j=1,NY
                    do 333 i=1,NX
                      read(12+k)raw(k,iy,im,id,i,j)
333 continue
!!Test  X Y
tlon=blon+(50-1)*intr
tlat=blat+(50-1)*intr
iy=1998-1960
im=7
id=9
open(53,file='test.txt')
write(53,*)'19980709,ix=  ,iy=  ,tlon= , tlat=  ,rainfall=',50,50,tlon,tlat,raw(2,iy,im,id,50,50)
tlon=blon+(51-1)*intr
tlat=blat+(50-1)*intr
write(53,*)'19980709,ix=  ,iy=  ,tlon= , tlat=  ,rainfall=',51,50,tlon,tlat,raw(2,iy,im,id,51,50)
tlon=blon+(50-1)*intr
tlat=blat+(51-1)*intr
write(53,*)'19980709,ix=  ,iy=  ,tlon= , tlat=  ,rainfall=',50,51,tlon,tlat,raw(2,iy,im,id,50,51)
close(53)
call output(raw,yy,mm,dd,ndd,rgns)

end program

  subroutine output(input,yy,mm,dd,ndd,rgns)
  implicit none
  integer,parameter :: NX=142
  integer,parameter :: NY=82
  real input(2,50,12,31,NX,NY)
  integer yy(6),mm(6),dd(6),ndd(6)            !!! 1-6 ETP,WTP,PRD,MLYR,NPC,NEC
  character rgns(6)*4
  real dout(2),lon,lat !
  real tlons,tlone,tlats,tlate
  real temp(2),XYD,kk(2)
  integer iys,iye,ims,ime,ids,ide
  real lons(6),lone(6),lats(6),late(6)
  integer day(12)
  integer iy,im,nd,id,ig,i,j,k,iyy,idd
  integer idss
  character syear*4,smonth*2,sday*2,snd*3
  character pathout*66, fpath*66

  pathout='D:\MyPaper\PhD04\Data\RainCN05\'
!-----------ETP--------------------------------------------------
  lons(1)=90;lone(1)=100
  lats(1)=27.5;late(1)=35
!-----------WTP--------------------------------------------------
  lons(2)=80;lone(2)=90
  lats(2)=27.5;late(2)=25
!-----------PRD--------------------------------------------------
  lons(3)=110;lone(3)=118
  lats(3)=21;late(3)=25
!-----------MLYR--------------------------------------------------
  lons(4)=110;lone(4)=122
  lats(4)=27;late(4)=33
!-----------NPC--------------------------------------------------
  lons(5)=112;lone(5)=120
  lats(5)=34;late(5)=42
!-----------NEC--------------------------------------------------
  lons(6)=120;lone(6)=130
  lats(6)=43;late(6)=49
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------------------------------------
  do i=1,12
    day(i)=30
  enddo
  day(1)=31;day(3)=31;day(5)=31;day(7)=31
  day(8)=31;day(10)=31;day(12)=31
  open(52,file='date.txt')
!!!!!------ the output data period
  do 330 ig=1,6  !!!! regions       
         tlons= lons(ig);tlone= lone(ig)
         tlats= lats(ig);tlate= late(ig)
         iyy=yy(ig)
         im =mm(ig)
         ids=dd(ig)
         write(syear,'(I4)') iyy
         write(smonth,'(I2.2)') im
         write(sday,'(I2.2)') ids
         write(snd,'(I3.3)') ndd(ig) 
         fpath=trim(pathout)//trim(rgns(ig))//syear//smonth//sday//'_'//snd//'.txt'
         open(101,file=trim(fpath))
         iy=iyy-1960
         day(2)=28
!-------------- 
        if(mod(iyy,4)==0.and.mod(iyy,100)/=0)then
          day(2)=29
        elseif(mod(iyy,400)==0)then
          day(2)=29
        endif
!----------------------------
        nd=day(im)
        idss=ids
        do 334 idd =1,ndd(ig)
           id=idss+idd-1
           if(id>nd)then
              id=1
              idss=1
              im=im+1
              if(im>12) then
                print*,'error'
                stop
              elseif(im<13)then
                nd=day(im)
              endif
            endif
            write(52,*)ig,iy+1960,im,id
           kk=0
           temp=0.
           dout=0.
          do 335 k=1,2  ! temp, rian 
          do 335 j=1,NY
             lat=14.75+0.5*(j-1)
          do 335  i=1,NX
             lon=69.75+0.5*(i-1)
             if(input(k,iy,im,id,i,j)/=-9999.0)then
             if(lon>=tlons.and.lon<=tlone)then 
             if(lat>=tlats.and.lat<=tlate)then
                temp(k)=temp(k)+input(k,iy,im,id,i,j)
                kk(k)=kk(k)+1.
             endif
             endif
             endif
335       continue 
           do k=1,2
            dout(k)=temp(k)/kk(k)
           end do 
          write(101,99) idd, dout(1),dout(2)
334 continue
        close(101)
330 continue
99    format(1X,I3,2(1X,f10.4))
  end subroutine