PROGRAM GETDATAFORPLOTTING
IMPLICIT NONE
!--------PARAMETERS OF MODLE (INPUT DATA)----------------
INTEGER, PARAMETER :: NT=2880,NZ=34,NX=200,DT=15  ! DT OUTPUT INTERVAL TIMESTEP
INTEGER, PARAMETER :: ND=30,NDT=24*60/DT
INTEGER, PARAMETER :: NBIN=99  ! BINS FOR CLOUD WATER CONTENT
REAL QL(NZ,NX),QI(NZ,NX),TL(NZ)
REAL DEN(NZ),QRL(NZ,NX),QRS(NZ,NX)
REAL OMG(NZ,NX),MAXOMG(NX,2)
REAL PRECI(NX+2),PBL(NX+2),FSH(NX+2),FLH(NX+2)
REAL RAINRATE

! THRESHOLD VALUES SETTINGS  -----------------------------
REAL QCC , QCE , CWP ! LOW LIMIT FOR CONVECTION CLOUD (G/KG)! LOW LIMIT FOR CLOUD ENSEMBLE (G/KG)! CLOUD WATER PATH THRESHOLD (G/M2)                   
                                ! --IF CWP>0 THEN QC0=CWP/(D DZ) ELSE =QCC
                                ! --QC0 IS USED FOR CLOUD IDENTIFICATION
DATA QCC,QCE,CWP/1.0E-2,1.0E-4,0.2/

REAL WATERBIN(NBIN)  ! IN G/KG
!-----------------------------------------------------------
REAL RAINBIN(35)
REAL DUMXRAINPDF(NDT,35),DUMXRAINDC(NDT,NZ)
INTEGER IRAIN,KK,IDUMXRAINDC(NDT)
!----------------------------------------------------------
!------- RESULT ARRAYS  ----------------------------------
! MEAN PROFILES 
REAL  MPDCQRL(NZ), MPDCQRS(NZ), MPDCQL(NZ),MPDCQI(NZ)  & ! DEEP CONVECTION
&    ,MPDCOMG(NZ)
REAL  AMPDCQRL(NZ), AMPDCQRS(NZ), AMPDCQL(NZ),AMPDCQI(NZ)  & ! DEEP CONVECTION 
&    ,AMPDCOMG(NZ)
REAL  MPSTQRL(NZ), MPSTQRS(NZ), MPSTQL(NZ),MPSTQI(NZ)  & !  STRATIFORM
&    ,MPSTOMG(NZ)
REAL  MPCRQRL(NZ), MPCRQRS(NZ), MPCRQL(NZ),MPCRQI(NZ)  & ! CIRRUS
&    ,MPCROMG(NZ)
REAL  MPACQRL(NZ), MPACQRS(NZ), MPACQL(NZ),MPACQI(NZ)  & ! ALLCELLS
&    ,MPACOMG(NZ)
! MEAN DIUNAL VALUE
REAL  DUDCQRL(NDT,NZ), DUDCQRS(NDT,NZ), DUDCQL(NDT,NZ),DUDCQI(NDT,NZ) & ! DEEP CONVECTION
&    ,DUDCOMG(NDT,NZ)
REAL  ADUDCQRL(NDT,NZ), ADUDCQRS(NDT,NZ), ADUDCQL(NDT,NZ),ADUDCQI(NDT,NZ) & ! DEEP CONVECTION
&    ,ADUDCOMG(NDT,NZ)
REAL  DUSTQRL(NDT,NZ), DUSTQRS(NDT,NZ), DUSTQL(NDT,NZ),DUSTQI(NDT,NZ) & !  STRATIFORM
&    ,DUSTOMG(NDT,NZ)
REAL  DUCRQRL(NDT,NZ), DUCRQRS(NDT,NZ), DUCRQL(NDT,NZ),DUCRQI(NDT,NZ) & ! CIRRUS
&    ,DUCROMG(NDT,NZ)
REAL  DUACQRL(NDT,NZ), DUACQRS(NDT,NZ), DUACQL(NDT,NZ),DUACQI(NDT,NZ) & ! ALLCELLS
&    ,DUACOMG(NDT,NZ)
! THE PDF OF THE MAXIMUM CLOUD WATER CONTENT WITH HEIGH
REAL  BINACMAX(NBIN,NZ), BINDCMAX(NBIN,NZ), BINSTMAX(NBIN,NZ) &
&    ,BINCRMAX(NBIN,NZ),ABINDCMAX(NBIN,NZ) 
REAL  DEEPCON(NX/2,NZ,2)
REAL FCDCDY(NDT,NZ),FCSTDY(NDT,NZ),FCCRDY(NDT,NZ),FCACDY(NDT,NZ) &
&   ,AFCDCDY(NDT,NZ)
REAL FCDDC(NDT,NZ),FCDST(NDT,NZ),FCDCR(NDT,NZ),AFCDDC(NDT,NZ)
!--------- RAINFALL ANALYS--------------------------------------------
REAL CON_PR(NDT,3), STR_PR(NDT,3), TAL_PR(NDT,3)
!---------------------------------------------------------
INTEGER I,J,K,IK,IZ,L,IX
INTEGER K1,K2,IDT,IKK
INTEGER KB(99),KE(99),NA
INTEGER ITDC,ITST,ITCR, ITADC              ! CONT FOR DEEP CONVECTION, STRATIFORM,CIRRUS
INTEGER IDDC,IDST,IDCR,IDADC,IDTAL                 ! CONT FOR DEEP CONVECTION, STRATIFORM,CIRRUS, DIUNAL CYCLE
INTEGER IDCZ(NZ),ISTZ(NZ),ICRZ(NZ),IADCZ(NZ),IACZ(NZ)
INTEGER IUDC(NDT),IUDCK(NDT,NZ)
INTEGER IUADC(NDT),IUADCK(NDT,NZ)
INTEGER IUST(NDT),IUSTK(NDT,NZ)
INTEGER IUCR(NDT),IUCRK(NDT,NZ)
INTEGER IUAC(NDT),IUACK(NDT,NZ)
INTEGER IMAXDC,IMAXST,IMAXAC,IMAXCR,IMAXADC
INTEGER ICONP,ISTRP,ITALP
INTEGER IUCONP(NDT),IUSTRP(NDT),IUTALP(NDT)
INTEGER IUPREMX
!=========================================================================
REAL FRC,FRR
REAL CM(99)
INTEGER IT
REAL C1,C2,C3,C4,C5,C6,C7,C8,C9    ! TEMPORARY VARABLES 
INTEGER IC1,IC2,IC3,IC4,IC5
CHARACTER*100 FPATH,DIRIN,DIROUT
CHARACTER FOLD*30,CASENM(6)*20,REGNM*20
CHARACTER DATESTR(6)*8
REAL RMXTP
INTEGER IPRE(200),IPREMX
!

CASENM(1)="ETPCTR_EC"  ; DATESTR(1)='20100603'
CASENM(2)="WTPCTR_EC"  ; DATESTR(2)='20100703'
CASENM(3)="NPCCTR_EC"  ; DATESTR(3)='20100802'
CASENM(4)="NECCTR_EC"  ; DATESTR(4)='20120706'
CASENM(5)="MLYRCTR_EC" ; DATESTR(5)='20100602'
CASENM(6)="PRDCTR_EC"  ; DATESTR(6)='20120401'
DIRIN="D:\MyPaper\PhD04\Cases\"
DIROUT="D:\MyPaper\PhD04\Cases\postdata\CTREC\"
!
WATERBIN(1)=0.005  ! MIN
C1=WATERBIN(1)
IC1=1
C3=0.005
IC2=-1
IC3=-1
IC4=-1
DO I =2, NBIN
    C2=C1+(I-IC1)*C3
    IF (C2>=0.005 .AND. C2 <0.5 .AND. IC2< 0)THEN
        C1=0.005 ; IC1=I ; C3=0.05
        IC2=1    ! MAKE SURE THE IF BLOCK JUST CALLED ONCE, SO IC1 IS RIGHT
    ELSEIF(C2>=0.5 .AND. C2<1. .AND. IC3< 0)THEN
        C1=0.5 ; IC1=I ; C3=0.05
        IC3=1
    ELSEIF(C2>=1 .AND. IC4< 0)THEN
        C1=1. ; IC1=I ; C3=0.15
        IC4=1
    ENDIF
    WATERBIN(I)=C2
ENDDO
DO I=1,35
    RAINBIN(I)=I*1.0
ENDDO
PRINT*,WATERBIN
DO I =1,6
	IF (CASENM(I)(1:3)=="MLY") THEN
		REGNM=CASENM(I)(1:4)
	ELSE
		REGNM=CASENM(I)(1:3)
	ENDIF
	FPATH=TRIM(DIRIN)//TRIM(REGNM)//'/CTREC'//DATESTR(I)  &
   &   //'/Simulation/'//TRIM(CASENM(I))//'_OMGQRLQRS.TXT'
	OPEN(20,FILE=TRIM(FPATH))
    FPATH=TRIM(DIRIN)//TRIM(REGNM)//'/CTREC'//DATESTR(I)  &
   &   //'/Simulation/'//'PRECI_'//TRIM(CASENM(I))//'.TXT'
    OPEN(30,FILE=TRIM(FPATH))
!
! -------- INTINITAL ARRAYS ----------------------	
	ITDC = 0 ; ITST = 0 ; ITCR = 0 ; ITADC = 0
	IDDC = 0 ; IDST = 0 ; IDCR = 0 ; IDADC = 0
    MPDCQRL = 0. ; MPDCQRS = 0. ; MPDCQL = 0. ; MPDCQI = 0.   ! DEEP CONVECTION
    MPDCOMG = 0.
    AMPDCQRL = 0. ; AMPDCQRS = 0. ; AMPDCQL = 0. ; AMPDCQI = 0.   ! DEEP CONVECTION
    AMPDCOMG = 0.
    MPSTQRL = 0. ; MPSTQRS = 0. ; MPSTQL = 0. ; MPSTQI = 0.   !  STRATIFORM
    MPSTOMG = 0.
    MPCRQRL = 0. ; MPCRQRS = 0. ; MPCRQL = 0. ; MPCRQI = 0.   ! CIRRUS
    MPCROMG = 0. 
    DUDCQRL = 0. ; DUDCQRS = 0. ; DUDCQL = 0. ; DUDCQI = 0.   ! DEEP CONVECTION
    DUDCOMG = 0.
    ADUDCQRL = 0. ; ADUDCQRS = 0. ; ADUDCQL = 0. ; ADUDCQI = 0.   ! DEEP CONVECTION
    ADUDCOMG = 0.
    DUSTQRL = 0. ; DUSTQRS = 0. ; DUSTQL = 0. ; DUSTQI = 0.   !   
    DUSTOMG = 0.
    DUCRQRL = 0. ; DUCRQRS = 0. ; DUCRQL = 0. ; DUCRQI = 0.   ! CIRRUS
    DUCROMG = 0. 
	FCDCDY  = 0. ; FCSTDY  = 0. ; FCCRDY = 0. ; FCACDY = 0.
    AFCDCDY  = 0.
    CON_PR  = 0. ; STR_PR  = 0. ; TAL_PR  = 0. 
    BINACMAX =0. ; BINDCMAX=0.  ; BINSTMAX=0. ; BINCRMAX =0.
    ABINDCMAX=0.
    DUMXRAINPDF=0.
    IMAXDC  = 0  ; IMAXST  = 0  ; IMAXAC = 0. ; IMAXCR=0 ; IMAXADC= 0
    IDCZ = 0  ; ISTZ = 0 ; ICRZ = 0  ; IADCZ = 0 ; IACZ=0
    IUDC =0   ; IUDCK= 0 ; IUSTK= 0  ; IUADCK= 0
    IUST = 0  ; IUCRK= 0 ; IUCR = 0  ; IUADC = 0
    ICONP =0  ; ISTRP= 0 ; ITALP =0  ; IPREMX= 0 ; IUPREMX=0
    IUCONP =0  ; IUSTRP= 0 ; IUTALP =0 ;IDUMXRAINDC=0
	IDT=1
	DO IT =1,NT
		IF (IDT < NDT)THEN
			IDT=IDT+1
		ELSE
			IDT=1
		ENDIF
        READ(30,'(8E12.4)')PRECI(:),PBL(:),FSH(:),FLH(:)
        RMXTP=PRECI(2)*1000.*3600.
        IPREMX=1
        IPRE=0
        DO IX=1,NX
            RAINRATE=PRECI(IX+1)*1000.*3600. ! CONVERT M/S TO MM/HR
            IF (RMXTP<RAINRATE) THEN
                RMXTP= RAINRATE
                IPREMX=IX
            ENDIF
        ENDDO
		DO IX=1,NX
            MAXOMG(IX,1)= -999.
            MAXOMG(IX,2)=  999.
            RAINRATE=PRECI(IX+1)*1000.*3600. ! CONVERT M/S TO MM/HR          
			DO IK=1,NZ
				READ(20,'(8E12.4)')QL(IK,IX),QI(IK,IX),OMG(IK,IX),   &
		   &	  QRS(IK,IX),QRL(IK,IX)
				TL(IK)=QL(IK,IX)+QI(IK,IX) ! TOTAL WATER CONTENT
                IF(TL(IK)<1.0E-3) TL(IK)=0.0   ! IF THERE IS NO CLOUD TL=0
                IF (OMG(IK,IX)>MAXOMG(IK,1))THEN
                    MAXOMG(IX,1)=OMG(IK,IX)*10. ! MAX POSITIVE VALUE
               !     IF (MAXOMG(IX,1)>10.) PRINT*,MAXOMG(IX,1),'AAA'
                ELSEIF(OMG(IK,IX)<MAXOMG(IK,2))THEN
                    MAXOMG(IX,2)=OMG(IK,IX)*10. ! MAX NEGITIVE VALUE
                ENDIF              
			ENDDO
!------------------ PRECIPITATION ------------------------------------
            IF (RAINRATE>25. .OR. MAXOMG(IX,1)>10.)THEN  ! CONVETIVE RAIN
                CON_PR(IDT,1)=CON_PR(IDT,1)+1
                ICONP=ICONP+1
                CON_PR(IDT,2)=CON_PR(IDT,2)+RAINRATE
                IUCONP(IDT)=IUCONP(IDT)+1
            ELSEIF(RAINRATE>0.001)THEN   !
                STR_PR(IDT,1)=STR_PR(IDT,1)+1
                ISTRP=ISTRP+1
                STR_PR(IDT,2)=STR_PR(IDT,2)+RAINRATE
                IUSTRP(IDT)=IUSTRP(IDT)+1
            ENDIF
            IF (RAINRATE>0.001)THEN
                TAL_PR(IDT,1)=TAL_PR(IDT,1)+1
                ITALP=ITALP+1
                TAL_PR(IDT,2)=TAL_PR(IDT,2)+RAINRATE
                IUTALP(IDT)=IUTALP(IDT)+1
            ENDIF
!-----------FOR DEEP CONVECTIONS BY CLOUD TOP AND CLOUD BASE---------------------------------------
			CALL DEEPCC(TL,NZ,KB,KE,NA)
			DO L = 1, NA
    			IF (KB(L).LE.5 .AND. KE(L).GE.21 .AND. RAINRATE>1. ) THEN  !! 5 AND 24 ARE THE VERTICAL LEVELS
!----------------------- FOR PDF -----------------------------------------------
     				K1 = KB(L)  ; K2 = KE(L)
     				FCDCDY(IDT,K2)=FCDCDY(IDT,K2)+1.  ! RECORD THE SAMPLE NUMBER WITH K2 TOP
                    IUDC(IDT)=IUDC(IDT)+1 ! RECORD ALL THE SAMPLES AT THE TIME OF IDT
                    C1=TL(K1)
                    IC1=K1
                    IF(IX==IPREMX)IUPREMX=K2
    				DO IKK =K1,K2 !LOOP FOR MEAN PROFILES
						MPDCQRL(IKK)=MPDCQRL(IKK)+QRL(IKK,IX)
						MPDCQRS(IKK)=MPDCQRS(IKK)+QRS(IKK,IX)
						MPDCQL(IKK)=MPDCQL(IKK)+QL(IKK,IX)
						MPDCQI(IKK)=MPDCQI(IKK)+QI(IKK,IX)
						MPDCOMG(IKK)=MPDCOMG(IKK)+OMG(IKK,IX)
						IDCZ(IKK)=IDCZ(IKK)+1
						DUDCQRL(IDT,IKK)=DUDCQRL(IDT,IKK)+QRL(IKK,IX)
						DUDCQRS(IDT,IKK)=DUDCQRS(IDT,IKK)+QRS(IKK,IX)
						DUDCQL(IDT,IKK)=DUDCQL(IDT,IKK)+QL(IKK,IX)
						DUDCQI(IDT,IKK)=DUDCQI(IDT,IKK)+QI(IKK,IX)
						DUDCOMG(IDT,IKK)=DUDCOMG(IDT,IKK)+OMG(IKK,IX)
                        IUDCK(IDT,IKK)=IUDCK(IDT,IKK)+1
                        IF (TL(IKK)>C1) THEN
                            C1=TL(IKK)
                            IC1=IKK
                        ENDIF
    				ENDDO
                    CALL GETBIN(WATERBIN,NBIN,C1,IC2)
                    BINDCMAX(IC2,IC1)=BINDCMAX(IC2,IC1)+1.
                    IMAXDC=IMAXDC+1
    			ENDIF
			ENDDO
            CALL GETBIN(RAINBIN,35,RMXTP,IRAIN)
            DUMXRAINPDF(IDT,IRAIN)=DUMXRAINPDF(IDT,IRAIN)+1.
            IF (IUPREMX>0)THEN
                KK=IUPREMX
                DUMXRAINDC(IDT,KK)=DUMXRAINDC(IDT,KK)+1
                IDUMXRAINDC(IDT)=IDUMXRAINDC(IDT)+1
            ENDIF
! FOR ALL CLOUD CELLS AND OTHER CLOUD TYPES INCLUDING THE DEEPCONVECTION BUT BY PRECIPITATION ADN MAX VELOCITY
            NA=0
            CALL INFCLD(TL,NZ,KB,KE,CM,NA) !KB BASE; KE TOP, NA: LAYERS OF CLOUDS
            DO L =1, NA
                K1 = KB(L)  ; K2 = KE(L)
                FCACDY(IDT,K2)=FCACDY(IDT,K2)+1.  ! RECORD THE SAMPLE NUMBER WITH K2 TOP
                IUAC(IDT)=IUAC(IDT)+1 ! RECORD ALL THE SAMPLES AT THE TIME OF IDT
                C1=TL(K1)
                IC1=K1
                DO IKK =K1,K2 !LOOP FOR MEAN PROFILES
                    MPACQRL(IKK)=MPACQRL(IKK)+QRL(IKK,IX)
                    MPACQRS(IKK)=MPACQRS(IKK)+QRS(IKK,IX)
                    MPACQL(IKK)=MPACQL(IKK)+QL(IKK,IX)
                    MPACQI(IKK)=MPACQI(IKK)+QI(IKK,IX)
                    MPACOMG(IKK)=MPACOMG(IKK)+OMG(IKK,IX)
                    IACZ(IKK)=IACZ(IKK)+1
                    DUACQRL(IDT,IKK)=DUACQRL(IDT,IKK)+QRL(IKK,IX)
                    DUACQRS(IDT,IKK)=DUACQRS(IDT,IKK)+QRS(IKK,IX)
                    DUACQL(IDT,IKK)=DUACQL(IDT,IKK)+QL(IKK,IX)
                    DUACQI(IDT,IKK)=DUACQI(IDT,IKK)+QI(IKK,IX)
                    DUACOMG(IDT,IKK)=DUACOMG(IDT,IKK)+OMG(IKK,IX)
                    IUACK(IDT,IKK)=IUACK(IDT,IKK)+1
                    IF (TL(IKK)>C1) THEN
                        C1=TL(IKK)
                        IC1=IKK
                    ENDIF
                ENDDO
                CALL GETBIN(WATERBIN,NBIN,C1,IC2)
                BINACMAX(IC2,IC1)=BINACMAX(IC2,IC1)+1.
                IMAXAC=IMAXAC+1
                !  FOLLOWING ARE FOR DIFFERENT CLOUD TYPE
                IF ((K1.LE.5 .AND. K2.GE.21) .OR. MAXOMG(IX,1)>10. &    !!!! k2 GE 24 OR 20 
               &     .OR. RAINRATE>25.) THEN
                    AFCDCDY(IDT,K2)=AFCDCDY(IDT,K2)+1.  ! RECORD THE SAMPLE NUMBER WITH K2 TOP
                    IUADC(IDT)=IUADC(IDT)+1 ! RECORD ALL THE SAMPLES AT THE TIME OF IDT
                    C1=TL(K1)
                    IC1=K1
                    DO IKK =K1,K2 !LOOP FOR MEAN PROFILES
                        AMPDCQRL(IKK)=AMPDCQRL(IKK)+QRL(IKK,IX)
                        AMPDCQRS(IKK)=AMPDCQRS(IKK)+QRS(IKK,IX)
                        AMPDCQL(IKK)=AMPDCQL(IKK)+QL(IKK,IX)
                        AMPDCQI(IKK)=AMPDCQI(IKK)+QI(IKK,IX)
                        AMPDCOMG(IKK)=AMPDCOMG(IKK)+OMG(IKK,IX)
                        IADCZ(IKK)=IADCZ(IKK)+1
                        ADUDCQRL(IDT,IKK)=ADUDCQRL(IDT,IKK)+QRL(IKK,IX)
                        ADUDCQRS(IDT,IKK)=ADUDCQRS(IDT,IKK)+QRS(IKK,IX)
                        ADUDCQL(IDT,IKK)=ADUDCQL(IDT,IKK)+QL(IKK,IX)
                        ADUDCQI(IDT,IKK)=ADUDCQI(IDT,IKK)+QI(IKK,IX)
                        ADUDCOMG(IDT,IKK)=ADUDCOMG(IDT,IKK)+OMG(IKK,IX)
                        IUADCK(IDT,IKK)=IUADCK(IDT,IKK)+1
                        IF (TL(IKK)>C1) THEN
                            C1=TL(IKK)
                            IC1=IKK
                        ENDIF
                    ENDDO
                    CALL GETBIN(WATERBIN,NBIN,C1,IC2)
                    ABINDCMAX(IC2,IC1)=ABINDCMAX(IC2,IC1)+1.
                    IMAXADC=IMAXADC+1
                ELSEIF (K2.LT.21) THEN ! STRATIFORM 
                    FCSTDY(IDT,K2)=FCSTDY(IDT,K2)+1.  ! RECORD THE SAMPLE NUMBER WITH K2 TOP
                    IUST(IDT)=IUST(IDT)+1 ! RECORD ALL THE SAMPLES AT THE TIME OF IDT
                    C1=TL(K1)
                    IC1=K1
                    DO IKK =K1,K2 !LOOP FOR MEAN PROFILES
                        MPSTQRL(IKK)=MPSTQRL(IKK)+QRL(IKK,IX)
                        MPSTQRS(IKK)=MPSTQRS(IKK)+QRS(IKK,IX)
                        MPSTQL(IKK)=MPSTQL(IKK)+QL(IKK,IX)
                        MPSTQI(IKK)=MPSTQI(IKK)+QI(IKK,IX)
                        MPSTOMG(IKK)=MPSTOMG(IKK)+OMG(IKK,IX)
                        ISTZ(IKK)=ISTZ(IKK)+1
                        DUSTQRL(IDT,IKK)=DUSTQRL(IDT,IKK)+QRL(IKK,IX)
                        DUSTQRS(IDT,IKK)=DUSTQRS(IDT,IKK)+QRS(IKK,IX)
                        DUSTQL(IDT,IKK)=DUSTQL(IDT,IKK)+QL(IKK,IX)
                        DUSTQI(IDT,IKK)=DUSTQI(IDT,IKK)+QI(IKK,IX)
                        DUSTOMG(IDT,IKK)=DUSTOMG(IDT,IKK)+OMG(IKK,IX)
                        IUSTK(IDT,IKK)=IUSTK(IDT,IKK)+1
                        IF (TL(IKK)>C1) THEN
                            C1=TL(IKK)
                            IC1=IKK
                        ENDIF
                    ENDDO
                    CALL GETBIN(WATERBIN,NBIN,C1,IC2)
                    BINSTMAX(IC2,IC1)=BINSTMAX(IC2,IC1)+1.
                    IMAXST=IMAXST+1
                ELSEIF(K1.GE.24 ) THEN  ! CIRRUS 
                    FCCRDY(IDT,K2)=FCCRDY(IDT,K2)+1.  ! RECORD THE SAMPLE NUMBER WITH K2 TOP
                    IUCR(IDT)=IUCR(IDT)+1 ! RECORD ALL THE SAMPLES AT THE TIME OF IDT
                    C1=TL(K1)
                    IC1=K1
                    DO IKK =K1,K2 !LOOP FOR MEAN PROFILES
                        MPCRQRL(IKK)=MPCRQRL(IKK)+QRL(IKK,IX)
                        MPCRQRS(IKK)=MPCRQRS(IKK)+QRS(IKK,IX)
                        MPCRQL(IKK)=MPCRQL(IKK)+QL(IKK,IX)
                        MPCRQI(IKK)=MPCRQI(IKK)+QI(IKK,IX)
                        MPCROMG(IKK)=MPCROMG(IKK)+OMG(IKK,IX)
                        ICRZ(IKK)=ICRZ(IKK)+1
                        DUCRQRL(IDT,IKK)=DUCRQRL(IDT,IKK)+QRL(IKK,IX)
                        DUCRQRS(IDT,IKK)=DUCRQRS(IDT,IKK)+QRS(IKK,IX)
                        DUCRQL(IDT,IKK)=DUCRQL(IDT,IKK)+QL(IKK,IX)
                        DUCRQI(IDT,IKK)=DUCRQI(IDT,IKK)+QI(IKK,IX)
                        DUCROMG(IDT,IKK)=DUCROMG(IDT,IKK)+OMG(IKK,IX)
                        IUCRK(IDT,IKK)=IUCRK(IDT,IKK)+1
                        IF (TL(IKK)>C1) THEN
                            C1=TL(IKK)
                            IC1=IKK
                        ENDIF
                    ENDDO
                    CALL GETBIN(WATERBIN,NBIN,C1,IC2)
                    BINCRMAX(IC2,IC1)=BINCRMAX(IC2,IC1)+1.
                    IMAXCR=IMAXCR+1
                ENDIF
            ENDDO ! NA
		ENDDO ! NX
    ENDDO ! NT
! ------- DOING THE AVERAGED  AND OUTPUT --------------------------------------
!-----------FOR DEEP CONVECTION 
    FPATH=TRIM(DIROUT)//TRIM(CASENM(I))//"_DEEPCONVECTION_GETPLOTF90.TXT"
    CALL MEANOUTPUT(NDT,NZ,NBIN,FCDCDY,IUDC,DUDCQRL,DUDCQRS,DUDCQL,&
        &    DUDCQI,  DUDCOMG,IUDCK, MPDCQRL,MPDCQRS,MPDCQL,MPDCQI, &
        &    MPDCOMG,  IDCZ, BINDCMAX, IMAXDC,FPATH)
    !-----------FOR DEEP CONVECTION 
    FPATH=TRIM(DIROUT)//TRIM(CASENM(I))//"_DEEPCONVECTION_A_GETPLOTF90.TXT"
    CALL MEANOUTPUT(NDT,NZ,NBIN,AFCDCDY,IUADC,ADUDCQRL,ADUDCQRS,ADUDCQL,&
        &    ADUDCQI,  ADUDCOMG,IUADCK, AMPDCQRL,AMPDCQRS,AMPDCQL,AMPDCQI, &
        &    AMPDCOMG,  IADCZ, ABINDCMAX, IMAXADC,FPATH)
    FPATH=TRIM(DIROUT)//TRIM(CASENM(I))//"_ALLCELLS_GETPLOTF90.TXT"
    CALL MEANOUTPUT(NDT,NZ,NBIN,FCACDY,IUAC,DUACQRL,DUACQRS,DUACQL,&
        &    DUACQI,  DUACOMG,IUACK, MPACQRL,MPACQRS,MPACQL,MPACQI, &
        &    MPACOMG,  IACZ, BINACMAX, IMAXAC,FPATH)
    FPATH=TRIM(DIROUT)//TRIM(CASENM(I))//"_STRATIFORM_GETPLOTF90.TXT"
    CALL MEANOUTPUT(NDT,NZ,NBIN,FCSTDY,IUST,DUSTQRL,DUSTQRS,DUSTQL,&
        &    DUSTQI,  DUSTOMG,IUSTK, MPSTQRL,MPSTQRS,MPSTQL,MPSTQI, &
        &    MPSTOMG,  ISTZ, BINSTMAX, IMAXST,FPATH)
    FPATH=TRIM(DIROUT)//TRIM(CASENM(I))//"_CIRRUS_GETPLOTF90.TXT"
    CALL MEANOUTPUT(NDT,NZ,NBIN,FCCRDY,IUCR,DUCRQRL,DUCRQRS,DUCRQL,&
        &    DUCRQI,  DUCROMG,IUCRK, MPCRQRL,MPCRQRS,MPCRQL,MPCRQI, &
        &    MPCROMG,  ICRZ, BINCRMAX, IMAXCR,FPATH)
!----------- FOR STRATIFORM     
ENDDO ! REGIONS

END PROGRAM
!
SUBROUTINE DEEPCC(QW,KM,KB,KE,NA) ! QW TOTAL CLOUD WATER,IM: HORIZONTAL GRID; KM : VERTICAL GRID
!     ------------------------------
!     ------------------------------
!
!     OUTPUT DEEP CONVECTION
!
INTEGER IM,KM,OU,I,K,L,KB(99),KE(99),NA
REAL QW(KM),C(99),CM(99)
INTEGER IDC
!      WRITE(OU,'(A3,1X,I4,1X,\)')'ITT', IT !  \ NOT CHANGE LINE
DO K = 1,KM
    IF (QW(K).GE.1.0E-3) C(K) = 1. ! THIS GRIG IS COVERED BY DEEP CONVECTION CLOUD
    IF (QW(K).LT.1.0E-3) C(K) = 0.
ENDDO
CALL INFCLD(C,KM,KB,KE,CM,NA)
RETURN
END SUBROUTINE
!#
SUBROUTINE INFCLD(C,NL,KB,KE,CM,NA) !KB BASE; KE TOP, NA: LAYERS OF CLOUDS
!     -----------------------------------
!     -----------------------------------
!
!     FIND INFORMATION ABOUT ADJACENT CLOUD LAYERS
!
IMPLICIT NONE
!
INTEGER NL,KB(*),KE(*),NA,K1,K2,K
REAL C(NL),CM(*),AA
!
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
END SUBROUTINE
SUBROUTINE GETBIN(BIN,NB,C1,IC)
IMPLICIT NONE
INTEGER NB,IC
REAL BIN(NB),C1
INTEGER I,K
IC=0
!print*,C1
IF (BIN(NB)<=C1) THEN 
   IC=NB
   GOTO 100
ENDIF
IF (BIN(1)>=C1) THEN 
   IC=1
   GOTO 100
ENDIF
DO I =2, NB
    IF (BIN(I-1)<=C1 .AND. BIN(I)>C1) THEN
        IC=I
        GOTO 100
    ENDIF
ENDDO
100 CONTINUE
if (ic==0)then
print*,c1,BIN(NB),BIN(1)
endif
RETURN
END SUBROUTINE
SUBROUTINE MEANOUTPUT(NDT , NZ  , NBIN , FCDY , IU , DUQRL, &
        &          DUQRS , DUQL, DUQI ,DUOMG ,IUK , MPQRL, MPQRS,  &
        &          MPQL  , MPQI, MPOMG, IZ, BINMAX, IMAX , FPATH  )
IMPLICIT NONE
INTEGER NDT,NZ,NBIN
INTEGER IU(NDT),IUK(NDT,NZ),IZ(NZ),IMAX
REAL    FCDY(NDT,NZ),DUQRL(NDT,NZ),DUQRS(NDT,NZ),DUQL(NDT,NZ), &
    &   DUQI(NDT,NZ),DUOMG(NDT,NZ), MPQRL(NZ), MPQRS(NZ), MPQL(NZ), &
    &   MPQI(NZ), MPOMG(NZ), BINMAX(NBIN,NZ)
CHARACTER FPATH*100
INTEGER I,J,K
REAL TMP,TMPU(NZ)
!
OPEN(10,FILE=TRIM(FPATH))
print*,NDT
TMP=0.
TMPU=0.
DO I=1,NDT
    TMP=TMP+IU(I)
ENDDO
    DO I =1,NDT
        IF(TMP>0)THEN
            DO K =1,NZ
                FCDY(I,K)=FCDY(I,K)*100./TMP  ! FREQUENCY
            ENDDO
        ELSE
            FCDY(I,:)=-999         !THERE IS NO DEEP CLOUD THIS TIME
        ENDIF
        WRITE(10,99)(FCDY(I,K),K=1,NZ)
        DO K=1,NZ
            IF (IUK(I,K)>0)THEN 
                DUQRL(I,K)=DUQRL(I,K)/IUK(I,K)
                DUQRS(I,K)=DUQRS(I,K)/IUK(I,K)
                DUQL(I,K)=DUQL(I,K)/IUK(I,K)
                DUQI(I,K)=DUQI(I,K)/IUK(I,K)
                DUOMG(I,K)=DUOMG(I,K)/IUK(I,K)
            ELSE
                DUQRL(I,K)=-999
                DUQRS(I,K)=-999
                DUQL(I,K) =-999
                DUQI(I,K) =-999
                DUOMG(I,K)=-999
            ENDIF
        ENDDO
        WRITE(10,99)(DUQRL(I,K),K=1,NZ)
        WRITE(10,99)(DUQRS(I,K), K=1,NZ)
        WRITE(10,99)(DUQL(I,K) , K=1,NZ)
        WRITE(10,99)(DUQI(I,K) , K=1,NZ)
        WRITE(10,99)(DUOMG(I,K), K=1,NZ)
    ENDDO
    DO K =1,NZ 
        IF (IZ(K)>0)THEN 
            MPQRL(K)=MPQRL(K)/IZ(K)
            MPQRS(K)=MPQRS(K)/IZ(K)
            MPQL(K)=MPQL(K)/IZ(K)
            MPQI(K)=MPQI(K)/IZ(K)
            MPOMG(K)=MPOMG(K)/IZ(K)
        ELSE
            MPQRL(K)=-999
            MPQRS(K)=-999
            MPQL(K) =-999
            MPQI(K) =-999
            MPOMG(K)=-999
        ENDIF                      
    ENDDO
    WRITE(10,99)(MPQRL(K),K=1,NZ)
    WRITE(10,99)(MPQRS(K), K=1,NZ)
    WRITE(10,99)(MPQL(K) , K=1,NZ)
    WRITE(10,99)(MPQI(K) , K=1,NZ)
    WRITE(10,99)(MPOMG(K), K=1,NZ)
    DO I=1,NBIN
        DO K =1, NZ
            IF(IMAX>1)THEN 
                BINMAX(I,K)=BINMAX(I,K)/IMAX
            ELSE
                BINMAX(I,K)=-999
            ENDIF
        ENDDO
        WRITE(10,99)(BINMAX(I,K), K=1,NZ)
    ENDDO
CLOSE(10)
99 FORMAT(1X,34(1X,E12.4))
RETURN
END SUBROUTINE