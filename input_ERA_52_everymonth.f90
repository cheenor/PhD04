PROGRAM SOUNDING
!-------------------------------------------------------------------
!      
!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CC
!!DATA PROVIDED BY XIN LIN (CSU, DICK JOHNSON'S STUDENT) 
!CC
!!WE ADD 3 LEVELS IN THE STRATOSPHERE (21, 31, 42 KM) ASSUMING
!!CONSTANT T TO GET DEEP SOUNDING.
!CC
!!NUMBER OF TIME LEVELS:
    PARAMETER(NTDAT=125) !!! 30DAYS 120+1, 31 DAYS 124+1
    DIMENSION TIDAT(NTDAT),PL1(NTDAT),PL2(NTDAT),PL3(NTDAT)
!!LEVEL OF TEST PLOT
    DATA KTEST /30/
!!BELOW ARE ARRAYS DEFINED AT CLARKS MODEL LEVELS:
    PARAMETER(L=52,LP=L+1)
    DIMENSION XI(LP),XIS(L),Z(L)
    DIMENSION TX(2),TZ(2)
    DATA RD,CP,G,RV,HLAT /287.,1005.,9.81,461.,2.5E6/
    DIMENSION TME1(L),THE1(L),QVE1(L),UE1(L),VE1(L),DTLS(L),DQLS(L)
    DIMENSION OUT1(L),OUT2(L),WE1(L),Q1LS(L),Q2LS(L),RHE1(L)
    DIMENSION OUTDY(L,4),OMG_ADJ(L)
!!NPIN IS NUMBER OF LEVELS
!     PARAMETER(NPIN0=39,NPIN=42)
    PARAMETER(NPIN0=37,NLTP=24,NPIN=26)
    DIMENSION PRESS(NPIN0),TEMP(NPIN0),THET(NPIN0),ZIN(NPIN0),  &
  &          VAP(NPIN0),UU(NPIN0),VV(NPIN0),WW(NPIN0),RH(NPIN0)
    DIMENSION TLSF(NPIN0),QLSF(NPIN0),XYN_OUT(NPIN0,19)
    REAL GZ(NPIN0),THE(NPIN0)
!CC
!!SST DATA FROM NMC ANALYSIS (TSST IS TIME IN DAYS, DSST IS IN DEG C)
    PARAMETER(NSST=9)
    CHARACTER*100 DIR,FILEPATH,PATH,FOLD,FILENAME,FILEPATH2
    CHARACTER*4 YEARSTR,AREA(4)
    CHARACTER MONSTR*2
    CHARACTER*16 DATE
    INTEGER IMS(4),IDS(4),IME(4),IDE(4),DAYS(12),TMID
    INTEGER IM,ID,KKK
    DIMENSION TSST(NSST),DSST(NSST)
!--------------------------------------------------------------------
    PARAMETER(NTM=125)
    REAL LHF(NTM),SHF(NTM)
!      REAL Q1TORAIN,Q2TORAIN
    DATA RD,CP,G,RV,HLAT /287.,1005.,9.81,461.,2.5E6/
!--------------------------------------------------------------------
    CHARACTER*16 CELORD,FOLD2
    REAL XYD(NPIN0,19),DYN(NPIN0,4)
!    DATA PRESS/1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.  &
!  &         ,750.,700.,650.,600.,550.,500.,450.,400.,350.,300.     &
!  &         ,250.,225.,200.,175.,150.,125.,100.,70.,50.,30.,20.    &
!  &         ,10.,7.,5.,3.,2.,1./
!---------------SET THE TIME ----------------------------------------
    DO I=1,12
      DAYS(I)=31
    ENDDO
! HOUR(1)='00:00'
! HOUR(2)='06:00'
! HOUR(3)='12:00'
! HOUR(4)='18:00'
    DAYS(2)=28
    DAYS(6)=30
    DAYS(4)=30
    DAYS(9)=30
    DAYS(11)=30   
    IYR=2010      
    IMS(1)=6  ;IME(1)=7
    IMS(2)=6  ;IME(2)=7
    IMS(3)=5  ;IME(3)=6
    IMS(4)=8  ;IME(4)=8
    IDS(1)=24  ;IDE(1)=23
    IDS(2)=24  ;IDE(2)=23
    IDS(3)=25  ;IDE(3)=24
    IDS(4)=1  ;IDE(4)=31
    AREA(1)='WTP'
    AREA(2)='ETP'
    AREA(3)='BOB'
    AREA(4)='NEC'  
    DO 1016 IYR=2010,2010
      FOLD2='INPUT'
      WRITE(YEARSTR,'(I4)')IYR
      FOLD=YEARSTR(3:4)//'0101-'//YEARSTR(3:4)//'1231\'
      DIR='X:\DATA\ERA_INTERIM\ERA_PRE_O\'
      ISTATUS1=CHDIR(TRIM(DIR)//TRIM(FOLD))
      ISTATUS2=SYSTEM("MD "//TRIM(FOLD2))
      FOLD2='INPUT/'
      DO 1015 IM=1,11
        DAYS(2)=28
        IF(MOD(IYR,4)==0.AND.MOD(IYR,100)/=0)THEN
          DAYS(2)=29 
        ELSEIF(MOD(IYR,400)==0)THEN
          DAYS(2)=29
        ENDIF
        IMS(1)=IM  ;IME(1)=IM
        IMS(2)=IM  ;IME(2)=IM
        IMS(3)=IM  ;IME(3)=IM
        IMS(4)=IM  ;IME(4)=IM
        IDS(1)=1  ;IDE(1)=DAYS(IM)
        IDS(2)=1  ;IDE(2)=DAYS(IM)
        IDS(3)=1  ;IDE(3)=DAYS(IM)
        IDS(4)=1  ;IDE(4)=DAYS(IM)
        WRITE(YEARSTR,'(I4)')IYR
        WRITE(MONSTR,'(I2.2)')IM
        FOLD=YEARSTR(3:4)//'0101-'//YEARSTR(3:4)//'1231\'
        DIR='Z:\DATA\LARGESCALE\TP\NCEPR2_PRE\'
        PATH=TRIM(DIR)//TRIM(FOLD)
        FILEPATH=TRIM(DIR)//TRIM(FOLD)//'DAYSTRT.TXT'
        OPEN(997,FILE=TRIM(FILEPATH))
        TX(1)=0.
        TX(2)=40.
        TZ(1)=0.
        TZ(2)=0.
!CCCCCCCCCCCCC!CODE BELOW TAKEN OUT FROM CLARKS MODEL SETUP:
        NZM=L-1 
!CC
!CC
        RAT=15.
        XI(2)=0.
        DO 152 K=2,NZM
          RATZ=RAT
          DEL=100.
          NZM1=L-1
          K1=K
          XI(K+1)=XI(K)+((RATZ-1.)/FLOAT(NZM1-2)*FLOAT(K1-2)+1.)*DEL
          PRINT 501,K+1,XI(K+1)
501       FORMAT(2X,'** GRID: K,XI:  ',I5,E16.8)
152     CONTINUE
        XI(1)=2.*XI(2)-XI(3)
        XI(LP)=2.*XI(L)-XI(L-1)
        DO K=1,L
          XIS(K)=.5*(XI(K)+XI(K+1))
          Z(K)=XIS(K)
!C PRINT*,Z(K)
        ENDDO
        OPEN(99,FILE='D:\MYPAPER\PHD04\FORTRANPROJECTS\INPUT_EVERYMONTH\Z-GEO.TXT')
        WRITE(99,9999)(Z(K),K=1,L)
9999    FORMAT(1X,52(1X,F10.3))
        CLOSE(99)  
!C!SST DATA: CONVERT TIME INTO HOURS
        DO III=1,NSST
          TSST(III)=TSST(III)*24.
        ENDDO
        NSTSST=2 ! STARTING INDEX FOR INTERPOLATION FROM SST DATASET
        IWRITE=0
!!!-----TOGA SOUNDING FILE OPEN--------------
!      OPEN(10,FILE=
!     *'/MNT/RAID50/HIBA/DATA_TOGA_OBS_DAT/IFA_DAT_NEW.SOUNDING'
!     *,STATUS='OLD')
!      OPEN(20,FILE=
!     *'/MNT/RAID50/HIBA/DATA_TOGA_OBS_DAT/IFA_DAT_NEW.FORCING'
!     *,STATUS='OLD')
!      OPEN(30,FILE=
!     *'/HOME/WUXQ/FORCING_TOGA/TOGA30/2DFORCING/MISC.IFA'
!     *,STATUS='OLD')
!      OPEN(40,FILE=
!     *'/MNT/RAID50/HIBA/DATA_TOGA_OBS_DAT/IFA_DAT.SURFACE'
!     *,STATUS='OLD')
!----------------------------------------------------------
        DO 1014 IP=1,2   ! AREA LOOPS   JUST FOR TP
          TEMPRESS=0
          TEMPTEMP=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IMTS=0
          IMTE=0 
          DO IMT=1,IMS(IP)-1
            IMTS=IMTS+DAYS(IMT)
          ENDDO
          IMTS=IMTS*4+IDS(IP)*4-4
          DO IMT=1,IME(IP)-1
            IMTE=IMTE+DAYS(IMT)
          ENDDO   
          IMTE=IMTE*4+IDE(IP)*4
          IMT=IMTE-IMTS    
          WRITE(997,*)IMTS ,IP
!-------------------------------------------------------------
          FILEPATH=TRIM(DIR)//TRIM(FOLD)//              &
     &        YEARSTR//TRIM(AREA(IP))//'_RAW.TXT'
          FILEPATH2=TRIM(DIR)//TRIM(FOLD)//             &
     &        YEARSTR//TRIM(AREA(IP))//'_SURFACE.TXT'    
          CALL OBS(IMTS,IP,FILEPATH,FILEPATH2,DAYS(IM)*4+1
     &                                       ,DAYS(IM)*4,MONSTR)
          OPEN(90,FILE=TRIM(FILEPATH))
          READ(90,*)
          DO I=1,IMTS*NPIN0
            READ(90,*)
          ENDDO
!----TARGET 1 2 3 .....IMTS IMTS+1 ......IMTE.......1460-----------------
!!!-------NCEP IMITATIVE SOUNDING FILES OPEN---------------
          ILEPATH=TRIM(DIR)//TRIM(FOLD)//TRIM(AREA(IP))//  &
     &                   '_SOUNDING.TXT'
!   PRESS(HPA) HEIGH(M) U(M/S) V(M/S)  OMEGA(PA/S) TEMP(K) THETA(K) QV(KG/KG) RH(%)
          OPEN(10,FILE=TRIM(FILEPATH))
          READ(10,*)
          READ(10,*)
          FILEPATH=TRIM(DIR)//TRIM(FOLD)//TRIM(AREA(IP))//  &
     &                   '_FORCING.TXT'
! TIMEID 17_LEVELS_T_FORCING(K/DAY) 17_LEVELS_QV_FORCING(K/DAY)
          OPEN(20,FILE=TRIM(FILEPATH))
          READ(20,*)
          READ(20,*)
!       FILEPATH=TRIM(DIR)//TRIM(FOLD)//'2010'//TRIM(AREA(IP))//'_RAW.TXT'
! 'DATE','HOUR','T_LS(K/DAY)','Q_LS(K/DAY)','U(M/S)',
!     +'V(M/S)','MOISTURE(KG/KG)','HGT(M)','AIR(K)','ADJ_OMEGA(PA/S)',
!     +'RH(%)','THETA(K)','Q1(K/DAY)','Q2(K/DAY)','HADQ(K/DAY)',
!     +'VADQ(K/DAY)','TCHQ(K/DAY)','HADT(K/DAY)','HADT(K/DAY)',
!     +'TCHT(K/DAY)','ORI_OMEGA(PA/S)'
!      OPEN(20,FILE=TRIM(FILEPATH))
! READ(20,*)

          FILEPATH=TRIM(DIR)//TRIM(FOLD)//TRIM(AREA(IP))//'_SURFACE22.TXT'
! PRESS(HPA) HEIGH(M) U(M/S) V(M/S)  OMEGA(PA/S) TEMP(K) THETA(K) QV(KG/KG) RH(%) FOR SURFACE
          OPEN(40,FILE=TRIM(FILEPATH))
          READ(40,*)
          READ(40,*)
          FILEPATH=TRIM(DIR)//TRIM(FOLD)//TRIM(AREA(IP))//'_HF.TXT'
! TIMEID LATENTHEAT(W/M^2) SENSIBLEHEAT(W/M^2)  PRECIPIPITABLE_WATER_FOR_ENTIRE_ATMOSPHERE(KG/M^2)
          OPEN(30,FILE=TRIM(FILEPATH))
          READ(30,*)
          READ(30,*)
          READ(30,*)
!          FILEPATH=TRIM(DIR)//TRIM(FOLD)//YEARSTR//         &
!    &             TRIM(AREA(IP))//'_DYN.TXT'
! TIMEID LATENTHEAT(W/M^2) SENSIBLEHEAT(W/M^2)  PRECIPIPITABLE_WATER_FOR_ENTIRE_ATMOSPHERE(KG/M^2)
!          OPEN(50,FILE=TRIM(FILEPATH))
!          READ(50,*)
!------------SKIP ----------------------------------------
          DO I=1,IMTS
            READ(30,*)
          ENDDO
          DO I=1,IMTS*(NPIN0+1)
            READ(10,*)
          ENDDO
          DO I=1,IMTS 
            READ(20,*)
          ENDDO
          DO I=1,IMTS
            READ(40,*)
          ENDDO
!          DO I=1,IMTS*17
!            READ(50,*)
!          ENDDO
!---------------------------------------------------
          OPEN(111,FILE='D:\MYPAPER\PHD04\FORTRANPROJECTS\'  &
     &           //'INPUT_EVERYMONTH\INPUT_EVERYMONTHZINGZ.TXT')
          DO 999 ITIM=1,DAYS(IM)*4+1   !!! WHY +1`, ONE MORE TIMESTEP THAN ONE MONTH  
          READ(30,301) TMID,FLH,FSH,PEWR,CPEWR
!                         PRINT*, TMID,'30'
301       FORMAT(1X,I4,1X,F8.2,1X,F8.2,1X,E12.4,1X,E12.4)
!C 876  FORMAT(4I5,4F8.2)
!READ SOUNDING DATA FOR THIS TIME LEVEL:
          READ(10,*)TMID
!                            PRINT*,TMID,'10',DAYS(IM)*4
          DO K=1,NPIN0
            READ(10,101) PRESS(K),GZ(K),UU(K),VV(K),WW(K), &
     &                  TEMP(K),THE(K),VAP(K)
          ENDDO
! PRINT*,PRESS(K),GZ(K),UU(K),VV(K),WW(K),TEMP(K),
!     + THE(K),VAP(K),RH(K)
!C 176  FORMAT(1X,7E18.8)
!READ FIRST LEVEL (SURFACE)
          ISFL=15  !!!! FOR TP, THE 14TH LEVLE (600HPA) IS THE SURFACE
          READ(40,102) ITTSP,PRESS(ISFL),GZ(ISFL),UU(ISFL),VV(ISFL),  &   !!! 4 = SURFACE
     &       WW(ISFL),TEMP(ISFL), THE(ISFL),VAP(ISFL),VQ1,VQ2
! PRINT*, PRESS(1),GZ(1),UU(1),VV(1),WW(1),TEMP(1),
!     + THE(1),VAP(1),RH(1)
          TEMPRESS=PRESS(ISFL)+TEMPRESS
          TEMPTEMP=TEMP(ISFL)+TEMPTEMP
101   FORMAT(1X,8(1X,e12.4))
102   FORMAT(1X,I4,10(1X,e12.4))
!             PRESS(1)=1008.
          WW(ISFL)=0.
! EXTRAPOLATE FIRST LEVEL (SURFACE)
!   PRESS(1)=1008.
!   COE2=(PRESS(1)-PRESS(2))/(PRESS(3)-PRESS(2))
!   TEMP(1)=COE2*TEMP(3) + (1.-COE2)*TEMP(2)
!   VAP(1)=COE2*VAP(3) + (1.-COE2)*VAP(2)
!   UU(1)=COE2*UU(3) + (1.-COE2)*UU(2)
!   VV(1)=COE2*VV(3) + (1.-COE2)*VV(2)

! READ LS FORCING DATA FOR THIS TIME LEVEL:
          TLSF(ISFL)=0.   ! 4=SURDACE
          QLSF(ISFL)=0.
          READ(20,201)TMID,(TLSF(K),K=1,NPIN0),     &
     &                     (QLSF(K),K=1,NPIN0)
201       FORMAT(1X,I4,74(1X,E12.4))

C      READ(20,776) (QLSF(K),K=2,NPIN0)
          TLSF(ISFL)=0.
          QLSF(ISFL)=0.
C 776  FORMAT(1X,5E20.8)
!--------------READ THE Q1 AND Q2 WHCIH ARE WRITTEN BY CHENJINGHUA 
          DO IK=1,NPIN0
            READ(90,906)TMID,(XYD(IK,KK),KK=1,18)
          ENDDO
906       FORMAT(1X,I4,18(1X,E12.4))
!          DO IK=1,NPIN0
!            READ(50,517)DATE,(DYN(IK,IR),IR=1,4)
!          ENDDO  
!517       FORMAT(1X,A16,1X,4(1X,E12.4))
          XYD(ISFL,10)=0.  !!!Q1 
          XYD(ISFL,11)=0.   !!!Q2
          XYD(ISFL,8) =0.   !!OMG_ADJ
!!!!!! THE FOLLOWING CODES ARE FOR TP
          DO K=1,NLTP
            IK=K+ISFL-1
            TLSF(K)=TLSF(IK)
            QLSF(K)=QLSF(IK)
            PRESS(K)=PRESS(IK)
            GZ(K)=GZ(IK)
            UU(K)=UU(IK)
            VV(K)=VV(IK)                          
            WW(K)=WW(IK)
            TEMP(K)=TEMP(IK)
            THE(K)=THE(IK)
            VAP(K)=VAP(IK)
            RH(K)=RH(IK)
            XYD(K,11)=XYD(IK,11)  !!!Q1 
            XYD(K,12)=XYD(IK,12)   !!!Q2
            XYD(K,8) =XYD(IK,8)   !!OMG_ADJ
!            DO IR=1,3
!              DYN(K,IR)=DYN(IK,IR)
!            END DO
          END DO
          DO K=NLTP+1,NPIN
            TLSF(K)=0.
            QLSF(K)=0.
          ENDDO
!CONVERT FROM TEMPERATURE (DEG C OR K) INTO POTENTIAL TEMPERATURE
          DO K=1,NLTP  !NPIN0   !!!!  MODEL LEVEL  NPIN0
            TEMP(K)=TEMP(K) ! +273.16
            THET(K)=TEMP(K)*(1.E3/PRESS(K))**(RD/CP)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            TLSF(K)=TLSF(K)*THET(K)/TEMP(K) ! THETA FORCING IN K/DAY
C      TLSF(K)=TLSF(K) ! TEMPERATURE FORCING IN K/DAY
            QLSF(K)=QLSF(K)*CP/HLAT           ! QV FORCING IN KG/KG/DAY
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ENDDO
!COMPUTE APPROXIMATED HEIGHT OF PRESSURE LEVELS:
          ZIN(1)=0.
          DO K=2,NLTP   !NPIN0  !!!!!!!
            KM=K-1
            TEMPK =TEMP(K ) * (1.+.6E-3*VAP(K ))   !!!!VAP  VAPOR MIXING  KG/KG
            TEMPKM=TEMP(KM) * (1.+.6E-3*VAP(KM))
            DELT=TEMPK-TEMPKM
            IF(DELT.GT.1.E-4) THEN
                TAVI=ALOG(TEMPK/TEMPKM)/DELT
            ELSE
                TAVI=1./TEMPK
            ENDIF
            DELTZ=-RD/(TAVI*G) * ALOG(PRESS(K)/PRESS(KM))
            ZIN(K)=ZIN(KM)+DELTZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ZIN(K)=GZ(K)/9.8
            WRITE(111,*)K, '  ZIN',ZIN(K),'   GZ',GZ(K)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ENDDO
! EXTRAPOLATE STRATOSPHERE:
!      ZIN(NPIN0+1)=21.E3
!      ZIN(NPIN0+2)=31.E3
!      ZIN(NPIN0+3)=42.5E3
!         PRINT*,ZIN(16),ZIN(17),ZIN(18)
! ZIN(16)=21.E3
!          ZIN(14)=35.E3
!          ZIN(15)=42.5E3
! STOP
          TEMP00=TEMP(NLTP)
          Z00=ZIN(NLTP)
          P00=PRESS(NLTP)
          DEN00=P00*1.E2/(RD*TEMP00)
          COE=G/(RD*TEMP00) 
! GET RH AT K=NPIN0:
          RH00=RH(NLTP)/100.
!          PRINT*,' RH AT K=NPIN0: ',RH00
!          DO K=NLTP-3,NLTP   !NPIN0-3,NPIN0
!            DEN=PRESS(K)*1.E2/(RD*TEMP(K))
!            ESAT=611.*               &
!     &           EXP(HLAT/RV * (1./273.16 - 1./TEMP(K)))
!            QVS00=ESAT/(RV*DEN*TEMP(K))
!            RH(K)=0.1*RH(K-1)
!            VAP(K)=RH(K)/100.*QVS00*1.E3
!            VAP(K)=VAP(K)*1E-3  !!!!! UNIT  MUST FOLLOW THE INPUT DATA
!          ENDDO
          DO K=NLTP+1,NPIN
            PRESS(K)=P00*EXP(-COE*(ZIN(K)-Z00))
!                            IF(K==NPIN) PRINT*,PRESS(K),P00,ZIN(K),Z00,K
            TEMP(K)=TEMP(NLTP)
            THET(K)=TEMP(K)*(1.E3/PRESS(K))**(RD/CP)
            DEN=PRESS(K)*1.E2/(RD*TEMP(K))
            ESAT=611.*                    &
     +          EXP(HLAT/RV * (1./273.16 - 1./TEMP(K)))
            QVS00=ESAT/(RV*DEN*TEMP(K))
C     VAP(K)=0.1*RH00*QVS00*1.E3
C     RH(K)=0.1*RH00*100.
            RH(K)=0.1*RH(K-1)
            VAP(K)=RH(K)/100.*QVS00*1.E3
            VAP(K)=VAP(K)*1E-3  !!!!! UNIT  MUST FOLLOW THE INPUT DATA
            UU(K)=UU(NLTP)
            VV(K)=VV(NLTP)
            WW(K)=WW(NLTP)
            XYD(K,11)= XYD(NLTP,11)
            XYD(K,12)= XYD(NLTP,12)
            XYD(K,8)= XYD(NLTP,8)
          ENDDO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FROM 1000 OR SURFACE
C WRITE INITIAL SOUNDING TO FORT.33 IF ITIM IS THE STARTING TIME:
          FILENAME=TRIM(DIR)//TRIM(FOLD)//TRIM(FOLD2)//  &
     &             TRIM(AREA(IP))//MONSTR//'.33'
          OPEN(33,FILE=TRIM(FILENAME))
!          IF(ITIM.EQ.17) THEN
          IF(ITIM.EQ.1) THEN   !!! THE FIRST STEP OF SOUNDING
!            MAKE SURE THE PRESSURE IS DECLINING ALONG WITH THE ALTITUDE  
            NPIN00=NPIN
            KKK=1
            WRITE(33,701) (PRESS(K),K=KKK,NPIN00)
701         FORMAT(5X,16H  DATA PRESS  / /  &
     &        5X,3H1  ,7(F7.2,1H,)/   &
     &        5X,3H1  ,7(F7.2,1H,)/   &
     &        5X,3H1  ,7(F7.2,1H,)/   &
!     &        5X,3H1  ,7(F7.2,1H,)/   &
!     &        5X,3H1  ,7(F7.2,1H,)/   &
     &        5X,3H1  ,2(F7.2,1H,),F7.2,1H/) 
              WRITE(33,702) (TEMP(K)-273.16,K=KKK,NPIN00)
702           FORMAT(5X,15H  DATA TEMP  / /   &
     &        5X,3H1  ,7(F7.2,1H,)/   &
     &        5X,3H1  ,7(F7.2,1H,)/   &
     &        5X,3H1  ,7(F7.2,1H,)/   &
!     &        5X,3H1  ,7(F7.2,1H,)/   &
!     &        5X,3H1  ,7(F7.2,1H,)/   &
     &        5X,3H1  ,2(F7.2,1H,),F7.2,1H/)  
              WRITE(33,703) (VAP(K)*1000,K=KKK,NPIN00)
703           FORMAT(5X,14H  DATA VAP  / / &
     &        5X,3H1  ,5(E10.3,1H,)/ &
     &        5X,3H1  ,5(E10.3,1H,)/ &
     &        5X,3H1  ,5(E10.3,1H,)/ &
     &        5X,3H1  ,5(E10.3,1H,)/ &
!     &        5X,3H1  ,5(E10.3,1H,)/ &
!     &        5X,3H1  ,5(E10.3,1H,)/ &
!     &        5X,3H1  ,5(E10.3,1H,)/ &
!     &        5X,3H1  ,5(E10.3,1H,)/ &
     &        5X,3H1  ,3(E10.3,1H,),E10.3,1H/) 
              WRITE(33,704) (UU(K),K=KKK,NPIN00)
704           FORMAT(5X,12H  DATA U  / / &
     &        5X,3H1  ,7(F7.2,1H,)/   &
     &        5X,3H1  ,7(F7.2,1H,)/   &
     &        5X,3H1  ,7(F7.2,1H,)/   &
!     &        5X,3H1  ,7(F7.2,1H,)/   &
!     &        5X,3H1  ,7(F7.2,1H,)/   &
     &        5X,3H1  ,2(F7.2,1H,),F7.2,1H/) 
              WRITE(33,705) (VV(K),K=KKK,NPIN00)
705           FORMAT(5X,12H  DATA V  / / &
     &        5X,3H1  ,7(F7.2,1H,)/   &
     &        5X,3H1  ,7(F7.2,1H,)/   &
     &        5X,3H1  ,7(F7.2,1H,)/   &
!     &        5X,3H1  ,7(F7.2,1H,)/   &
!     &        5X,3H1  ,7(F7.2,1H,)/   &
     &        5X,3H1  ,2(F7.2,1H,),F7.2,1H/) 
! FLAG TO WRITE OTHER DATA
              ITIMS=ITIM
              IWRITE=1
            ENDIF  !END OF FIRTST STEP
!COMPUTE ENVIRONMENTAL PROFILES FROM SOUNDING ASSUMING NO TOPOGRAPHY:
            IISN=1
            RHE1(IISN)=RH(IISN)
            THE1(IISN)=THET(IISN)
            TME1(IISN)=TEMP(IISN)
            QVE1(IISN)=VAP(IISN)      !*1.E-3
            UE1(IISN)=UU(IISN)
            VE1(IISN)=VV(IISN)
            WE1(IISN)=WW(IISN)
            PRESST=PRESS(IISN)
C      Q1LS(IISN)=XYD(IISN,11)
C   Q2LS(IISN)=XYD(IISN,12)
!INTEGRATE UPWARDS:
            FILENAME=TRIM(DIR)//TRIM(FOLD)//TRIM(FOLD2)//  &    
    &                   TRIM(AREA(IP))//MONSTR//'_UV_PROFILES.35'
            OPEN(35,FILE=TRIM(FILENAME))  
            FILENAME=TRIM(DIR)//TRIM(FOLD)//TRIM(FOLD2)//  & 
    &                 TRIM(AREA(IP))//MONSTR//'_LSFORCING.37'
            OPEN(37,FILE=TRIM(FILENAME))
            FILENAME=TRIM(DIR)//TRIM(FOLD)//TRIM(FOLD2)//  &
    &                 TRIM(AREA(IP))//MONSTR//'_SURFACE.39' 
            OPEN(39,FILE=TRIM(FILENAME))
            FILENAME=TRIM(DIR)//TRIM(FOLD)//TRIM(FOLD2)//  &
    &                 TRIM(AREA(IP))//MONSTR//'.49'
            OPEN(49,FILE=TRIM(FILENAME))
            FILENAME=TRIM(DIR)//TRIM(FOLD)//TRIM(FOLD2)//  &
    &                     TRIM(AREA(IP))//MONSTR//'_THETAQV_PROFILE.41'   !!!!!! UNIT THETA(K)  QV G/KG 
            OPEN(41,FILE=TRIM(FILENAME))
            FILENAME=TRIM(DIR)//TRIM(FOLD)//TRIM(FOLD2)//  &
    &               TRIM(AREA(IP))//MONSTR//'.43'
            OPEN(43,FILE=TRIM(FILENAME))
            FILENAME=TRIM(DIR)//TRIM(FOLD)//TRIM(FOLD2)//  &
    &                 TRIM(AREA(IP))//MONSTR//'.99'
            OPEN(99,FILE=TRIM(FILENAME))
            FILENAME=TRIM(DIR)//TRIM(FOLD)//TRIM(FOLD2)//  &
    &                 TRIM(AREA(IP))//MONSTR//'.DYN'
            OPEN(999,FILE=TRIM(FILENAME))
            FILENAME=TRIM(DIR)//TRIM(FOLD)//TRIM(FOLD2)//  &
    &                 TRIM(AREA(IP))//MONSTR//'_DIAGNOSED_RAIN.TXT'
            OPEN(44,FILE=TRIM(FILENAME))
!                   WRITE(44,*)'TIMEID Q1TORAIN Q2TORAIN'
            DO 64 K=2,L
              DO KK=2,NPIN  !!!! FROM 2 OR
                IISN=KK-1
!     PRINT*,ZIN(KK),Z(K),K,KK
                IF(ZIN(KK).GE.Z(K)) GO TO 665
              ENDDO
              PRINT*,' *** INPUT SOUNDING DOES NOT GO HIGH ENOUGH. STOP.'
              STOP 'SOUNDING'
 665          CONTINUE 
              COE2=(Z(K)-ZIN(IISN))/(ZIN(IISN+1)-ZIN(IISN))
C PRINT*,COE2,Z(K),ZIN(IISN)
              THE1(K)=COE2*THET(IISN+1) + (1.-COE2)*THET(IISN)
              TME1(K)=COE2*TEMP(IISN+1) + (1.-COE2)*TEMP(IISN)
              QVE1(K)=(COE2*VAP(IISN+1) + (1.-COE2)*VAP(IISN))!*1.E-3
              UE1(K)=COE2*UU(IISN+1) + (1.-COE2)*UU(IISN)
              VE1(K)=COE2*VV(IISN+1) + (1.-COE2)*VV(IISN)
              WE1(K)=COE2*WW(IISN+1) + (1.-COE2)*WW(IISN)
              DTLS(K)=COE2*TLSF(IISN+1) + (1.-COE2)*TLSF(IISN)
              DQLS(K)=COE2*QLSF(IISN+1) + (1.-COE2)*QLSF(IISN)
              Q1LS(K)=COE2*XYD(IISN+1,11) + (1.-COE2)*XYD(IISN,11)
              Q2LS(K)=COE2*XYD(IISN+1,12) + (1.-COE2)*XYD(IISN,12)
              OMG_ADJ(K)=COE2*XYD(IISN+1,8) + (1.-COE2)*XYD(IISN,8) 
              RHE1(K)=COE2*RH(IISN+1) + (1.-COE2)*RH(IISN)
!              DO IDY=1,4
!                OUTDY(K,IDY)=COE2*DYN(IISN+1,IDY) + 
!    &                       (1.-COE2)*DYN(IISN,IDY)
!              ENDDO  
64            CONTINUE

!SCALE AND WRITE TO FILES:
              IF(IWRITE.EQ.1) THEN
                ITIMS=1
                TIME = FLOAT(ITIM-ITIMS)*6.
                TIDAT(ITIM-ITIMS+1)=(ITIM-ITIMS)*6./24.
!  VELOCITY PROFILES FOR SELECTED PERIOD (FORT.35); NOTE ROTATION
                SVEL=10.    ! VELOCITY SCALE IN CLARKS MODEL
                DO K=1,L
                  OUT1(K)=-VE1(K)/SVEL
                OUT2(K)= UE1(K)/SVEL
                ENDDO
                WRITE(35,801) TIME,OUT1,OUT2
801             FORMAT(10F8.3)
!  PROFILES FOR L-S FORCING TERMS (FORT.37)
                DAY=24.*3600.
                OUT1(1)=0.
                OUT2(1)=0.
                DO K=2,L
                  OUT1(K)=DTLS(K)/DAY    ! NOW IN K/SEC
                  OUT2(K)=DQLS(K)/DAY    ! NOW IN KG/KG/SEC
!
!  SET FORCING TO ZERO ABOVE 17KM (17.19KM AT LEVEL 33)
C         IF(K.GE.33) THEN
C         OUT1(K)=0.
C         OUT2(K)=0.
C         ENDIF
!CCCCCCC
                ENDDO
                WRITE(37,802) TIME,OUT1,OUT2
                WRITE(99,802) TIME,Q1LS,Q2LS
                WRITE(999,802)TIME,OUTDY(:,1),OUTDY(:,2),  &
     &                OUTDY(:,3),OUTDY(:,4),OMG_ADJ
802             FORMAT(8E12.4)
                PL1(ITIM-ITIMS+1)=OUT1(KTEST)*DAY
                PL2(ITIM-ITIMS+1)=OUT2(KTEST)*1.E3*DAY

!  TIME SERIES OF OCEAN SURFACE THETA AND QV (FORT.39)
C     TIME2=TIME+ITIMS*6.
C     DO IIII=1,100
C     IF(TIME2.GT.TSST(NSTSST))   NSTSST=NSTSST+1 
C     ENDDO
CCCCC!INTERPOLATE OCEAN TEMP:
C      COE2=(TIME2-TSST(NSTSST-1))/(TSST(NSTSST)-TSST(NSTSST-1))
C        SST=COE2*DSST(NSTSST) + (1.-COE2)*DSST(NSTSST-1)
                SST=TEMP(1)
                THS=SST*(THE1(1)+THE1(2))/(TME1(1)+TME1(2))
                DEN=PRESS(1)*1.E2/(RD*TEMP(1))
                ESAT=611.*EXP(HLAT/RV * (1./273.16 - 1./SST))
                QVSS=ESAT/(RV*DEN*SST)
                WRITE(39,803) TIME,THS,QVSS,SST  !!! WHAT ARE THE UNITS?
                WRITE(49,903) TIME,SST,THS,QVSS
                TMPQ1=0.0
                TMPQ2=0.0
                DO IC=1,NPIN0
                  TMPQ1=TMPQ1+XYD(IC,11)
                  TMPQ2=TMPQ2+XYD(IC,12) !  DATA RD,CP,G,RV,HLAT /287.,1005.,9.81,461.,2.5E6/
                ENDDO
                Q2TORAIN=(((TMPQ2/24./3600.)*CP)/(HLAT*G*1000.)     &   !!! *1000 WATER DENITY
     &                               +FLH/(HLAT*1000.))*3600*1000.     !!!! M/S CONVERT TO MM/HR

                Q1TORAIN=(((TMPQ1/24./3600.)*CP/G)-FSH/CP)         &
                                        /(HLAT)                 !  *1000 WATER DENITY
     &                             *3600*1000   !!! CONVERT TO MM/HR
                WRITE(44,803)TIME,TMPQ1,Q1TORAIN,TMPQ2,Q2TORAIN

803             FORMAT(6E16.5)
903             FORMAT(4E16.5)
                PL3(ITIM-ITIMS+1)=SST
!  THETA AND QV PROFILES FOR SELECTED PERIOD (FORT.41)
                DO K=1,L
                  OUT1(K)=THE1(K)
                  OUT2(K)=QVE1(K)
                ENDDO
                WRITE(41,804) TIME,OUT1,OUT2
                WRITE(43,804) TIME,OUT1,OUT2,TME1,RHE1,UE1,VE1,WE1
804             FORMAT(8E13.5)
              ENDIF

999         CONTINUE
         
            WRITE(997,*)TEMPRESS/(DAYS(IM)*4+1.),IP,DAYS(IM),IM
            WRITE(997,*)TEMPTEMP/(DAYS(IM)*4+1.), IP,DAYS(IM),IM
1014      CONTINUE       
          CLOSE(997)
          CLOSE(99)
          CLOSE(90)
          CLOSE(10)
          CLOSE(20)
          CLOSE(30)
          CLOSE(40)
          CLOSE(50)
          CLOSE(33)
          CLOSE(35)
          CLOSE(37)
          CLOSE(39)
          CLOSE(49)
          CLOSE(43)
          CLOSE(41)
          CLOSE(999)
1015    CONTINUE
1016  CONTINUE
      STOP
    END PROGRAM


  SUBROUTINE OBS(ITT,IP,FRAW,FSURF,NTDAT,NTDYS,MONSTR)

  CHARACTER*100 FRAW,FSURF,FOUTS
  CHARACTER*16 CELORD(NTDAT) , CELORS(NTDAT)
  REAL XYD(NTDAT,37,19)
  REAL XYS(NTDAT,10)
  REAL Q1Q2(NTDAT,2)
  INTEGER NTDYS
  CHARACTER MONSTR*2
! 'DATE','HOUR','1 T_LS(K/DAY)','2 Q_LS(K/DAY)','3 U(M/S)',
!     +'4 V(M/S)','5 MOISTURE(KG/KG)','6 HGT(M)','7 AIR(K)','8 ADJ_OMEGA(PA/S)',
!     +'9 RH(%)','10 THETA(K)','11 Q1(K/DAY)','12 Q2(K/DAY)','13 HADQ(K/DAY)',
!     +'14 VADQ(K/DAY)','15 TCHQ(K/DAY)','16 HADT(K/DAY)','17 HADT(K/DAY)',
!     +'18 TCHT(K/DAY)','19 ORI_OMEGA(PA/S)'
  OPEN(20,FILE=TRIM(FRAW))
  READ(20,*)
  DO I=1,ITT*37
      READ(20,*)
  ENDDO
  OPEN(30,FILE=TRIM(FSURF))
  READ(30,*)
  DO I=1,ITT
      READ(30,*)
  ENDDO
  DO IT=1,NTDYS
  Q1T=0.
  Q2T=0.
    DO IK=1,37
      READ(20,906)ISKP,(XYD(IT,IK,KK),KK=1,18)
      Q1T=Q1T+XYD(IT,IK,11)
      Q2T=Q2T+XYD(IT,IK,12)
    ENDDO
    Q1Q2(IT,1)=Q1T
    Q1Q2(IT,2)=Q2T
    READ(30,907)ISKP,(XYS(IT,KK),KK=1,10) !LH(W/M^2)   SH(W/M^2)       PRATE      CPRATE  DLR(W/M^2)  ULR(W/M^2)
  ENDDO
!---------------OUTPUT-----------------------------------------
  ILEN=LEN_TRIM(FRAW)
  FOUTS=FRAW(1:ILEN-15)//'INPUT/'//FRAW(ILEN-15:ILEN-7)// & 
&        MONSTR//'OBS_SURFACE_INPUT.TXT'
  OPEN(40,FILE=TRIM(FOUTS))
  WRITE(40,*)'TIMESTEP  [Q1](K/D)   [Q2](K/D)  RAINFALL(MM/HR)'
  FOUTS=FRAW(1:ILEN-15)//'INPUT'//FRAW(ILEN-15:ILEN-7)//  &
&        MONSTR//'Q1Q2_PROFILES_INPUT.TXT'
  OPEN(50,FILE=TRIM(FOUTS))
  WRITE(50,*)'TIMESTEP Q1_LEVEL5_TO_17 Q2_LEVEL5_TO_17'  !!!! THE FIRST LEVEL IS NOT ZERO 
!
  DO IT=1,NTDAT
      WRITE(40,1014)IT*1.0,XYS(IT,9),XYS(IT,10),XYS(IT,3)
      WRITE(50,417)IT*1.0,(XYD(IT,KK,10),KK=14,37),   &
&        (XYD(IT,KK,11),KK=14,37)
  ENDDO 


906   FORMAT(1X,I4,18(1X,E12.4))
907   FORMAT(1X,I4,10(1X,E12.4))
1014  FORMAT(4E12.4)
417   FORMAT(8E12.4)
  CLOSE(20)
  CLOSE(30)
  RETURN
  END SUBROUTINE