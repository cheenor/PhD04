program DownTRMM2A25
implicit none
integer i,j,k
character yrstr*4,timestr*20,dystr*3
character*200 path, filehead,ftpnm
integer day(12),ista
integer iyy,ims,ime,tmpday
integer days,daye

filehead='ftp://disc2.nascom.nasa.gov/'
ftpnm='ftp/data/s4pa//TRMM_L2/TRMM_2A25/' 
!ftp://disc2.nascom.nasa.gov/ftp/data/s4pa//TRMM_L2/TRMM_2A25/2008/006/2A25.080106.57786.6.HDF.Z/2008/006/2A25.080106.57786.6.HDF.Z
do i=1,12
	day(i)=31
enddo
day(4)=30
day(6)=30
day(9)=30
day(11)=30
day(2)=28
ims=4
ime=9
do iy=2006,2010
	day(2)=28
	if((mod(iy,4).eq.0 .and. mod(iy,100) .ne.0) .or. mod(iy,400).eq.0)then
		day(2)=29
	endif
	tmpday=1
	do i=1,ims-1
		tmpday=tmpday+day(i)
	enddo
	days=tmpday
	tmpday=0
	do i=1,ime
		tmpday=tmpday+day(i)
	enddo
	daye=tmpday
	write(yrstr,"(I4.4)")iy
	open(10,file=yrstr//'2A25.sh')
	do i=days,daye
		write(dystr,"(I3.3)")i
		timestr=yrstr//'/'//dystr//'/'
		path=trim(filehead)//trim(ftpnm)//trim(timestr)//'2A25*.HDF.Z'
		write(10,100)'wget',trim(path)
	enddo
	ista=system("mkdir "//yrstr)
	ista=system("chmod +x "//yrstr//"2A25.sh")
	ista=system("mv "//yrstr//"2A25.sh ./"//yrstr)
enddo

end program

!filehead='ftp://disc2.nascom.nasa.gov/'
!ftpnm='ftp/data/s4pa//TRMM_L2/TRMM_2A25' 
!ftp://disc2.nascom.nasa.gov/ftp/data/s4pa//TRMM_L2/TRMM_2A25.6/2008/006/2A25.080106.57790.6.HDF.Z