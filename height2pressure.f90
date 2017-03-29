implicit none
real heigh(52)
real press(52),pressout(52)
real tm(52),surface(6),p0
character*8 rgname(6)
integer i,j,k,it
real TIME,OUT1(52),OUT2(52),UE1(52),VE1(52),WE1(52),tmp
real p1,p2
character*30 filestring(6)
character*100 fpath
real dz
data heigh/-50.000 ,    50.000 ,   164.286,    307.143,    478.571  ,  678.571 , &
    &  907.143 ,  1164.286,   1450.000,   1764.286 ,  2107.143,   2478.572 , &
    &  2878.572,   3307.143,  3764.286,  4250.000,   4764.286,   5307.143, &
    &  5878.571,   6478.571,   7107.143,  7764.286,  8450.000,  9164.285,  &
    &  9907.143,  10678.570,  11478.570,  12307.143,  13164.285,  14050.000, &
    &  14964.285,  15907.143,  16878.572,  17878.572,  18907.145,  19964.285, &
    &  21050.000,  22164.285,  23307.145,  24478.572,  25678.572,  26907.145, &
    &  28164.285,  29450.000,  30764.285,  32107.145,  33478.570,  34878.570, &
    &  36307.141,  37764.285,  39250.000,  40750.000/
data surface/590.0514,566.6624,986.2728,970.1259,950.8006,975.0000/!ETP,WTP,PRD,MLYR,NPC
data rgname/'ETP','WTP','PRD','MLYR','NPC','NEC'/
data filestring/'20100603','20100703','20120401','20100602','20100802','20120706'/


do i=1,6 
	fpath='D:\MyPaper\PhD04\Cases\ERA\FORCING\'//trim(rgname(i))//'\'&
	&//trim(rgname(i))//'_'//trim(filestring(i))//'_031d_ERA.43'
	open(43,file=trim(fpath))
	fpath='D:\MyPaper\PhD04\Cases\ERA\FORCING\'//trim(rgname(i))//'\'&
	&//trim(rgname(i))//'_'//trim(filestring(i))//'_031d_ERA_52pressure.52'
	open(52,file=trim(fpath))
	press=0.0
	do it=1,124
    	p0=surface(i)
	    press(1)=p0
		read(43,804) TIME,OUT1,OUT2,tm,UE1,VE1,WE1
		do k =2,52
    		dz=heigh(k)-heigh(k-1)
    		tmp=dz/(18400*(1+(tm(k)-273.15)/273))           
            p2=p0/(10**tmp)
            press(k)=p2/124.+press(k)
            print*,tmp,p2
            p0=p2
        enddo
    enddo
    
    write(52,*)press
804  FORMAT(8E13.5)
    close(52)
    close(43)
enddo
pause
end