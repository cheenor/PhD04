      PROGRAM DETECHING_CLOUD_CELLS
C      C     CLOUD COVER OVER ALL SAMPLES
C           DATA RESOLUTION: DT=15(MIN), DX=3(KM), DZ VARIES
C                            0015UTC JUN.22 1997 TO 2400UTC JUL.17 1997 (MD=26)
C
      IMPLICIT NONE
      INTEGER IM,KM,DT,MT,MD,NF,NTT
      PARAMETER (IM=200,KM=34,DT=15,MT=(24*60)/DT,MD=30,NF=7)
      PARAMETER (NTT=2880)
      INTEGER AT,MA,NT
      PARAMETER (AT=6*60,NT=AT/DT,MA=(MD*MT)/NT)
      REAL QL(IM,KM)            ! CLOUD WATER (G/KG)
     $    ,QI(IM,KM)            ! CLOUD ICE (G/KG)
C     $    ,QV(IM,KM)            ! WATER VAPOR MIXING RATIO (G/KG)
C     $    ,RH(IM,KM)            ! RELATIVE HUMIDITY (%)
C     $    , W(IM,KM)            ! VERTICAL VELOCITY (M/S)
     $    , TC(IM,KM)            ! TEMPERATURE (K)
C     $    , U(IM,KM)            ! WIND (M/S)
     $    ,    Z(KM)            ! VERTICAL COORDINATE (KM)
     $    ,    P(KM)            ! PRESSURE (MB)
     $    ,    D(KM)            ! AIR DENSITY (KG/M3)
     $    ,    DEN(KM)
C
      DATA Z/            0.0500000, 0.1643000, 0.3071000, 0.4786000
     $      , 0.6786000, 0.9071000, 1.1640000, 1.4500000, 1.7640001
     $      , 2.1070001, 2.4790001, 2.8789999, 3.3069999, 3.7639999
     $      , 4.2500000, 4.7639999, 5.3070002, 5.8790002, 6.4790001
     $      , 7.1069999, 7.7639999, 8.4499998, 9.1639996, 9.9069996
     $      ,10.6800003,11.4799995,12.3100004,13.1599998,14.0500002
     $      ,14.9600000,15.9099998,16.8799992,17.8799992,18.9099998/

      DATA P/             964.6000,  952.3000,  937.0000,  918.9000
     $      ,  898.1000,  874.7000,  849.0000,  821.1000,  791.4000
     $      ,  760.0000,  727.2000,  693.1000,  658.1000,  622.3000
     $      ,  586.0000,  549.3000,  512.6000,  476.1000,  440.0000
     $      ,  404.6000,  369.9000,  336.3000,  303.9000,  272.8000
     $      ,  243.2000,  215.2000,  188.9000,  164.7000,  142.5000
     $      ,  122.6000,  104.9000,   89.3700,   75.8500,   64.1600/

      DATA D/             1.101309,  1.086979,  1.069333,  1.048541
     $      ,  1.024801,  0.998338,  0.969394,  0.938233,  0.905128
     $      ,  0.870366,  0.834235,  0.797026,  0.759030,  0.720527
     $      ,  0.681790,  0.643081,  0.604644,  0.566707,  0.529478
     $      ,  0.493144,  0.457871,  0.423803,  0.391061,  0.359743
     $      ,  0.329927,  0.301668,  0.275005,  0.249953,  0.226517
     $      ,  0.204680,  0.184416,  0.165686,  0.148440,  0.132622/

C     ------------------------------------------------------------------
      INTEGER NG                ! NO OF GCM GRIDS IN THE X-DOMAIN
      PARAMETER (NG=5)
      INTEGER IGRID(2,NG),NGSS
      DATA IGRID/001,100, 101,200
     $          ,051,150, 151,050
     $          ,001,200/
      DATA NGSS/1/            !START LOOP FOT NG
      REAL QCC                  ! LOW LIMIT FOR CONVECTION CLOUD (G/KG)
     $    ,QCE                  ! LOW LIMIT FOR CLOUD ENSEMBLE (G/KG)
     $    ,CWP                  ! CLOUD WATER PATH THRESHOLD (G/M2)
                                ! --IF CWP>0 THEN QC0=CWP/(D DZ) ELSE =QCC
                                ! --QC0 IS USED FOR CLOUD IDENTIFICATION
      DATA QCC,QCE,CWP/1.0E-2,1.0E-4,0.0/
C      DATA QCC,QCE,CWP/1.0E-2,1.0E-4,0.2/
C     DATA QCC,QCE,CWP/1.0E-4,1.0E-6,0.2/

      REAL QC0(KM)              ! CLOUD WATER THRESHOLD (G/KG)
C
      INTEGER NE                ! NO CLOUD ENSEMBLE CLASSES
      PARAMETER (NE=1)
      REAL CELS(0:NE)           ! BOUNDARIES OF CE CLASSES
CX    DATA CELS/0.00,0.20,0.40,0.60,0.80,1.00/
      DATA CELS/0.00,1.00/
C
      INTEGER NO                ! NO OF VERTICAL OVERLAP CLOUD CELLS
      PARAMETER (NO=3)
C
      INTEGER NP                ! NO OF VERTICAL PROFILES
      PARAMETER (NP=3)
      REAL ZTD,ZBD              ! DEEP CONVECTION TOP, BASE HEIGHTS (KM)
      DATA ZTD,ZBD/10.0,7.0/  !   ?????
      REAL ZTP,ZBP              ! PBL CONVECTION TOP, BASE HEIGHTS (KM)
      DATA ZTP,ZBP/14.0,2.0/
      INTEGER KSC               ! SHALLOW CONVECTION LAYER THICKNESS
      DATA KSC/4/

      REAL ZBMIN                ! MINIMUM BASE HEIGHT TO COUNT CLOUD (KM)
      DATA ZBMIN/0.5/
C     ------------------------------------------------------------------
      REAL AL(MA,KM,NG),AI(MA,KM,NG)
C
      REAL ZD(KM),SZD
      REAL FT(KM,KM,NG),SG(NG)
      REAL FA(0:NO+1,NE),FB(KM,KM,NO,NE),SA(NE)
      REAL QD(KM,KM,NP),WD(KM,KM),QZ(KM),WZ(KM),TD(KM,KM)
      REAL QP(KM,KM,NP),WP(KM,KM),QB(KM),WB(KM)
      REAL QS(KM,NP),TS(KM),WS(KM),TR(KM)
      REAL FRG,FRC,FRR,SUM,QMX,SPV,QZD,TZD,EZD,CZD
     $    ,WK(IM,KM),C(KM),CM(KM),CE,WA(0:NO+1),WW(KM,KM,NO)
      INTEGER IT,ID,I,K,IG,IB,IE,II,NI,IS,NS,KB,KE,KCB(KM),KCE(KM),IO,IC
     $       ,LCB(KM),LCE(KM),LC,IP,IGN
      DATA SPV/1.0E20/
C
      INTEGER IU,OU,O1,O2,O3,IREC,NREC,LENCHR
      DATA IU,OU,O1,O2,O3/10,20,31,32,33/
      CHARACTER DIR*100,FINAME*120,FONAME*120
      CHARACTER PATH*100,FOLD*30,CASENMSTR(6)*20,DIRIN*100
     *          , FILENU*3,FPATH*150,DIROUT*100
      CHARACTER DATESTR(6)*30,CASENM*30
      INTEGER IH,IA,IK
!
      CASENMSTR(1)='PRDCTR_EC'   ; DATESTR(1)='20120401'
      CASENMSTR(2)='MLYRCTR_EC'  ; DATESTR(2)='20100602'
      CASENMSTR(3)='NPCCTR_EC'   ; DATESTR(3)='20100802'
      CASENMSTR(4)='NECCTR_EC'   ; DATESTR(4)='20120706'
      CASENMSTR(5)='ETPCTR_EC'   ; DATESTR(5)='20100603'
      CASENMSTR(6)='WTPCTR_EC'   ; DATESTR(6)='20100703'
      DO IGN =1,6
      CASENM=CASENMSTR(IGN)
      IF (CASENM(1:3) == "MLY") THEN
        FOLD=CASENM(1:4)
      ELSE
        FOLD=CASENM(1:3)
      ENDIF
      DIRIN='D:\MyPaper\PhD04\Cases\'
      DIROUT='D:\MyPaper\PhD04\Cases\'
      DIRIN=TRIM(DIRIN)//TRIM(FOLD)//'\CTREC'//
     + TRIM(DATESTR(IGN))//'\Simulation\'
      DIROUT=TRIM(DIROUT)//'postdata\'
      IH=81  
      FPATH=TRIM(DIRIN)//'/'//TRIM(CASENM)//'_QLQIT.TXT'
      OPEN(IH,FILE=TRIM(FPATH))
!     * ,FORM='UNFORMATTED',STATUS='OLD')
C
      FPATH = DIROUT(1:LENCHR(DIROUT))//TRIM(CASENM)//
     +  '_QL_QI_QT_f90.TXT'   
      OPEN(OU,FILE=FPATH(1:LENCHR(FPATH)))   
C-----INITIALIZATION

      DO IG = 1,NG
         DO K = 1,KM
            DO IA = 1,MA
               AL(IA,K,IG) = 0.
               AI(IA,K,IG) = 0.
            ENDDO
         ENDDO
      ENDDO
C-----LOOP FOR ALL TIME
      IA=0
      DO 111 IT=1,NTT
        DO IK=1,KM 
          READ(IH,99)DEN(IK),QL(:,IK),QI(:,IK),TC(:,IK)
!          print*, QL
        ENDDO
C-----LOOP FOR EACH GCM GRID
        IF (IA < NTT/30)THEN
            IA=IA+1
        ELSE
            IA=1         
        ENDIF  !IA = IT   !((ID-1)*MT + IT-1)/NT + 1
        DO IG = NGSS,NG ! WE WANT GRID FROM 1 TO 200 IN X     
      ! 1-200 GRIDS DIVIDED TO 1-100,101-200,51-150,151-50,1-200
          IB = IGRID(1,IG) ! GRID BEGIN  
          IE = IGRID(2,IG) ! GRID END
          NI = IE-IB
          IF (IE.LT.IB) NI = IM-IB+IE
          NI = NI+1
C........GET GCM GRID MEAN PROFILE
          FRC = NI
          FRC = 1./FRC
          DO I = 1,NI    ! NI
            II = I + IB-1
            IF (II.GT.IM) II = II - IM
            IF (II.LT.01) II = II + IM
            DO K = 1,KM
              AL(IA,K,IG) = AL(IA,K,IG) + FRC*QL(I,K)
              AI(IA,K,IG) = AI(IA,K,IG) + FRC*QI(I,K)
            ENDDO
         ENDDO
      ENDDO
  111 CONTINUE
      CLOSE(IH)
C-----OUTPUT
      DO IG = 1,NG
        WRITE(OU,99) ((AL(IA,K,IG),K=1,KM),IA=1,MA)   ! CLOUD LIQUID WATER
        WRITE(OU,99) ((AI(IA,K,IG),K=1,KM),IA=1,MA)   ! CLOUD ICE WATER
        WRITE(OU,99) ((AL(IA,K,IG)
     $                         +AI(IA,K,IG),K=1,KM),IA=1,MA)  ! TATOL CLOUD WATER
      ENDDO
      CLOSE(OU)    
C     --------------------
C     --------------------
!      SUBROUTINE COVER_ALL(CASENM,FOLD,DIRIN,DIROUT)
C     --------------------
C     --------------------
      IH=81  
      FPATH=TRIM(DIRIN)//'/'//TRIM(CASENM)//'_QLQIT.TXT'
      OPEN(IH,FILE=TRIM(FPATH))
C
C      FINAME = DIR(1:LENCHR(DIR))//'CRM_SGP.1997JUN22_JUL17'
C      OPEN(IU,FILE=FINAME(1:LENCHR(FINAME)),STATUS='OLD'
C-CRY$       ,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=((IM*KM+1)/2)*8)
C     $       ,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=(IM*KM)*4)
      OU=82
      FPATH =DIROUT(1:LENCHR(DIROUT))//TRIM(CASENM)//
     + '_OUTPUT_DH_f90.TXT'      ! ASCII
      OPEN(OU,FILE=FPATH(1:LENCHR(FPATH)) )   
C-----INITIALIZATION
!      ZD(1) = 0.5*(Z(2)-Z(1))*D(1)
!      DO K = 2,KM-1
!         ZD(K) = 0.5*(Z(K+1)-Z(K-1))*D(K)
!~      ENDDO
!     ZD(KM) = 0.5*(Z(KM)-Z(KM-1))*D(KM)
!
!      IF (CWP.GT.0.0) THEN
!         DO K = 1,KM
!            QC0(K) = 1.E-3*CWP/ZD(K)
!         ENDDO
!      ELSEIF (QCC.GT.0.0) THEN
!         DO K = 1,KM
!            QC0(K) = QCC
!         ENDDO
!      ELSE
!         STOP 'XXXX WRONG CLOUD ID THRESHOLD'
!      ENDIF
!
!      SZD = 0.
!      DO K = 1,KM
!         SZD = SZD + ZD(K)
!      ENDDO
      DO IP = 1,NP
         DO K = 1,KM
            DO KB = 1,KM
               QD(K,KB,IP) = 0.
            ENDDO
         ENDDO
         DO KB = 1,KM
            QS(KB,IP) = 0.
         ENDDO
      ENDDO
      DO K = 1,KM
         DO KB = 1,KM
            WD(K,KB) = 0.
            TD(K,KB) = 0.
         ENDDO
      ENDDO
      DO KB = 1,KM
         TS(KB) = 0.
         WS(KB) = 0.
         QZ(KB) = 0.
         WZ(KB) = 0.
      ENDDO
      DO IP = 1,NP
         DO K = 1,KM
            DO KE = 1,KM
               QP(K,KE,IP) = 0.
            ENDDO
         ENDDO
      ENDDO
      DO K = 1,KM
         DO KE = 1,KM
            WP(K,KE) = 0.
         ENDDO
      ENDDO
      DO KE = 1,KM
         QB(KE) = 0.
         WB(KE) = 0.
      ENDDO
      DO IG = 1,NG
         SG(IG) = 0.
         DO KE = 1,KM
            DO KB = 1,KM
               FT(KB,KE,IG) = 0.
            ENDDO
         ENDDO
      ENDDO
      DO IC = 1,NE
         SA(IC) = 0.
         DO IO = 0,NO+1
            FA(IO,IC) = 0.
         ENDDO
         DO IO = 1,NO
            DO KE = 1,KM
               DO KB = 1,KM
                  FB(KB,KE,IO,IC) = 0.
               ENDDO
            ENDDO
         ENDDO
      ENDDO

C-----INITIALIZE REFERENCE TEMPERATURE

C      FRR = MD*MT*IM !MD DAYS, MT TIMESTEPS PER DAY
      FRR = NTT*IM ! TOTAL TIME STEPS
      FRR = 1./FRR
      DO K = 1,KM
         TR(K) = 0.
      ENDDO
C-----LOOP FOR ALL TIMESTEPS
      QL=0.
      QI=0.
      TC=0.
      DO 100 IT=1,NTT
        DO IK=1,KM 
          READ(IH,99)DEN(IK),QL(:,IK),QI(:,IK),TC(:,IK)
        ENDDO
C      READ(IU,REC=IREC+1) QL
C      READ(IU,REC=IREC+2) QI
C      READ(IU,REC=IREC+6) T
        D(:)=DEN(:)
        ZD(1) = 0.5*(Z(2)-Z(1))*D(1) ! DH*DEN
        DO K = 2,KM-1
          ZD(K) = 0.5*(Z(K+1)-Z(K-1))*D(K)
        ENDDO
        ZD(KM) = 0.5*(Z(KM)-Z(KM-1))*D(KM)
        IF (CWP.GT.0.0) THEN
          DO K = 1,KM
            QC0(K) = 1.E-3*CWP/ZD(K)
          ENDDO
        ELSEIF (QCC.GT.0.0) THEN
          DO K = 1,KM
            QC0(K) = QCC
          ENDDO
        ELSE
          STOP 'XXXX WRONG CLOUD ID THRESHOLD'
        ENDIF
        SZD = 0.
        DO K = 1,KM
          SZD = SZD + ZD(K)
        ENDDO
C-----GET TOTAL CONDENSATE
        DO K = 1,KM
          DO I = 1,IM
            WK(I,K) = QL(I,K) + QI(I,K)
          ENDDO
        ENDDO
        CALL DEEPCC(WK,IM,KM,OU,IT) ! WK TOTAL CLOUD WATER
C-----GET REFERENCE TEMPERATURE = AVERAGE OVER TIME AND DOMAIN
        DO K = 1,KM
          DO I = 1,IM
            TR(K) = TR(K) + FRR*TC(I,K)
          ENDDO
        ENDDO
C-----LOOP FOR EACH GCM GRID
        DO IG = NGSS,NG
          IB = IGRID(1,IG)
          IE = IGRID(2,IG)
          NI = IE-IB
          IF (IE.LT.IB) NI = IM-IB+IE
          NI = NI+1
C........GET GCM GRID MEAN
          FRG = NI  ! GRID NUMBER OF EVERY EXPANDED GRID 
          FRG = 1./FRG
          FRC = KM*NI
          FRC = 1./FRC
          SUM = 0.
          EZD = 0.
          QMX = 0.
          DO I = 1,NI
            II = I + IB-1
            IF (II.GT.IM) II = II - IM
            IF (II.LT.01) II = II + IM
            DO K = 1,KM
              SUM = SUM + FRC*WK(II,K)
              EZD = EZD + FRG*WK(II,K)*ZD(K) !!!! WATER PATH?
              QMX = MAX(QMX , WK(II,K))
            ENDDO
          ENDDO
          EZD = EZD/SZD
!CX        WRITE(OU,'(A,3I4,2(1PE12.5))') 'SUM>',ID,IT,IG,SUM,QMX
C--------GET PDF FOR CLOUD ENSEMBLE
          IF (SUM.GE.QCE) THEN
            SG(IG) = SG(IG) + 1.
            IF (IG.LE.3) THEN
              DO IO = 1,NO+1   !NO: NO OF VERTICAL OVERLAP CLOUD CELLS
                WA(IO) = 0.
              ENDDO
              DO IO = 1,NO
                DO KE = 1,KM
                  DO KB = 1,KM
                    WW(KB,KE,IO) = 0.
                  ENDDO
                ENDDO
              ENDDO
              CE = 0.
            ENDIF
            DO I = 1,NI
              II = I + IB-1
              IF (II.GT.IM) II = II - IM
              IF (II.LT.01) II = II + IM
              DO K = 1,KM
                C(K) = WK(II,K)
                IF (C(K).LT.QC0(K)) C(K) = 0.
              ENDDO
              CALL INFCLD(C,KM,KCB,KCE,CM,NS) ! C WK TOTAL WATER; CM MAX WATER 
C...............FOR ALL CLOUD CELLS
              DO IS = 1,NS ! NS: LAYERS OF CLOUDS
                KB = KCB(IS)
                KE = KCE(IS)
                FT(KB,KE,IG) = FT(KB,KE,IG) + FRG !CLOUD COVER OF EVERY LAYER
              ENDDO
              IF (IG.LE.3) THEN ! 1-100,101-200,51-150
C...............VERTICAL PROFILES FOR DEEP CONVECTIONS
                DO IS = 1,NS  ! LAYER OF CLOUDS
                  KB = KCB(IS)
                  KE = KCE(IS)
                  IF (Z(KE).GE.ZTD.AND.Z(KB).LE.ZBD) THEN ! ZTD ZTB
                    CZD = 0.
                    QZD = 0.
                    DO K = KB,KE
                      CZD = CZD + ZD(K)  !!! TOTAL MASS FORM CLOUD BASE TO CLOUD TOP
                      QZD = QZD + ZD(K)*C(K) ! ! TOTAL WATER MASS FORM CLOUD BASE TO CLOUD TOP
                    ENDDO
                    QZD = QZD/CZD
                    QZ(KB) = QZ(KB) + QZD
                    WZ(KB) = WZ(KB) + 1.
                    DO K = KB,KE
                      QD(K,KB,1) = QD(K,KB,1) + C(K) ! TOTAL WATER CONTENT 
                      QD(K,KB,2) = QD(K,KB,2) + C(K)/QZD  !
                      QD(K,KB,3) = QD(K,KB,3) + C(K)/QZD*D(K) !
CX                    QD(K,KB,3) = QD(K,KB,3) + C(K)/EZD
                      TD(K,KB)   = TD(K,KB)   + TC(II,K)
                      WD(K,KB)   = WD(K,KB)   + 1.
                    ENDDO
                  ENDIF
                ENDDO
C...............VERTICAL PROFILES FOR PBL CONVECTIONS
                DO IS = 1,NS
                  KB = KCB(IS)
                  KE = KCE(IS)
                  IF (Z(KE).LE.ZTP.AND.Z(KB).LE.ZBP) THEN
                    CZD = 0.
                    QZD = 0.
                    DO K = KB,KE
                      CZD = CZD + ZD(K)
                      QZD = QZD + ZD(K)*C(K)
                    ENDDO
                    QZD = QZD/CZD
                    QB(KE) = QB(KE) + QZD
                    WB(KE) = WB(KE) + 1.
                    DO K = KB,KE
                      QP(K,KE,1) = QP(K,KE,1) + C(K)
                      QP(K,KE,2) = QP(K,KE,2) + C(K)/QZD
                      QP(K,KE,3) = QP(K,KE,3) + C(K)/EZD*D(K)
CX                     QP(K,KE,3) = QP(K,KE,3) + C(K)/EZD
                      WP(K,KE) = WP(K,KE) + 1.
                    ENDDO
                  ENDIF
                ENDDO
C..............VERTICAL PROFILES FOR SHALLOW CONVECTIONS
                IF (NS.EQ.1) THEN
                  KB = KCB(1)
                  KE = KCE(1)
                  IF (KE-KB.LT.KSC) THEN
                    CZD = 0.
                    QZD = 0.
                    TZD = 0.
                    DO K = KB,KE
                      CZD = CZD + ZD(K)
                      QZD = QZD + ZD(K)*C(K)
                      TZD = TZD + ZD(K)*TC(II,K)
                    ENDDO
                    QZD = QZD/CZD
                    TZD = TZD/CZD
                    QS(KB,1) = QS(KB,1) + QZD
                    QS(KB,2) = QS(KB,2) + EZD
                    QS(KB,3) = QS(KB,3) + QZD/EZD*D(KB)
CX                  QS(KB,3) = QS(KB,3) + QZD/EZD
                    TS(KB) = TS(KB) + TZD
                    WS(KB) = WS(KB) + 1.
                  ENDIF
                ENDIF
C...............FOR OVERLAP CLOUD CELLS
C               ----IGNORE THIN CELLS
                LC = 0
                DO IS = 1,NS
                  IF (KCE(IS)-KCB(IS).GE.2) THEN
                    LC = LC + 1
                    LCB(LC) = KCB(IS)
                    LCE(LC) = KCE(IS)
                  ENDIF
                ENDDO
C               ----IGNORE SMALL GAPS
                IF (LC.GT.1) THEN
                  DO IS = 2,LC
                    IF (LCB(IS)-LCE(IS-1).LT.2) THEN
                      LCE(IS-1) = LCE(IS)
                      LC = LC - 1
                      DO IC = IS,LC
                        LCB(IC) = LCB(IC+1)
                        LCE(IC) = LCE(IC+1)
                      ENDDO
                    ENDIF
                  ENDDO
                ENDIF
C              ----GET FREQUENCY ACCORDING TO OVERLAP CELLS
                IF (LC.GT.0) THEN 
                  CE = CE + FRG
                  IO = MIN(LC,NO+1)
                  WA(IO) = WA(IO) + FRG
                ENDIF
                IF (LC.EQ.1) THEN
                  KB = LCB(1)    ! BAS OF THE CELL
                  KE = LCE(1)    ! TOP OF THE CELL
                ELSEIF (LC.EQ.2) THEN
                  KB = LCE(1)    ! TOP OF LOW CELL
                  KE = LCB(2)    ! BAS OF UPP CELL
                ELSEIF (LC.EQ.3) THEN
                  KB = LCB(2)    ! BAS OF MID CELL
                  KE = LCE(2)    ! TOP OF MID CELL
                ENDIF
                IF (LC.GE.1.AND.LC.LE.3)
     &              WW(KB,KE,LC) = WW(KB,KE,LC) + FRG
              ENDIF
            ENDDO
            IF (IG.LE.3) THEN
              WA(0) = 1.0 - CE         ! CLEAR SKY 
              IC=0
              DO IS = 1,NE
                IF (CELS(IS-1).LT.CE.AND.CE.LE.CELS(IS)) IC = IS
              ENDDO
              IF(IC>0)THEN
                SA(IC) = SA(IC) + 1.
                DO IO = 0,NO+1
                    FA(IO,IC) = FA(IO,IC) + WA(IO)
                ENDDO
                DO IO = 1,NO
                    DO KE = 1,KM
                        DO KB = 1,KM
                          FB(KB,KE,IO,IC)=FB(KB,KE,IO,IC)+ WW(KB,KE,IO)
                        ENDDO
                    ENDDO
                ENDDO
              ENDIF
            ENDIF
          ENDIF
        ENDDO ! IG
C  200 CONTINUE
  100 CONTINUE  
C      CLOSE(IU)
      DO K = 1,KM
         TR(K) = TR(K) - 273.15
      ENDDO
C-----DO AVERAGE
      DO IG = 1,NG
        IF (SG(IG).GT.0.) THEN
          DO KB = 1,KM
            IF (Z(KB)<0.5) THEN   ! THE MINMUM CLOUD BASE IS 0.5KM
              DO KE = 1,KM
                FT(KB,KE,IG) = 0.0
              ENDDO
            ELSE
              DO KE = 1,KM
                FT(KB,KE,IG) = FT(KB,KE,IG)/SG(IG) * 100.
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO
!
      DO IP = 1,NP
        DO K = 1,KM
          DO KB = 1,KM
            IF (WD(K,KB).GT.0.) QD(K,KB,IP) = QD(K,KB,IP)/WD(K,KB)
          ENDDO
        ENDDO
      ENDDO
!
      DO KB = 1,KM
        IF (WZ(KB).GT.0.) QZ(KB) = QZ(KB)/WZ(KB)
      ENDDO
      DO IP = 1,NP
        DO K = 1,KM
          DO KE = 1,KM
            IF (WP(K,KE).GT.0.) QP(K,KE,IP) = QP(K,KE,IP)/WP(K,KE)
          ENDDO
        ENDDO
      ENDDO
!
      DO KE = 1,KM
        IF (WB(KE).GT.0.) QB(KE) = QB(KE)/WB(KE)
      ENDDO

      DO IP = 1,NP
        DO KB = 1,KM
          IF (WS(KB).GT.0.) QS(KB,IP) = QS(KB,IP)/WS(KB)
        ENDDO
      ENDDO
      DO KB = 1,KM
        QS(KB,NP) = QS(KB,NP)*QS(KB,NP-1)
      ENDDO
!
      DO K = 1,KM
        DO KB = 1,KM
          IF (WD(K,KB).GT.0.) TD(K,KB) = TD(K,KB)/WD(K,KB)
        ENDDO
      ENDDO
      DO KB = 1,KM
        IF (WS(KB).GT.0.) TS(KB) = TS(KB)/WS(KB)
      ENDDO
!
      DO IC = 1,NE
        IF (SA(IC).GT.0.) THEN
          DO IO = 0,NO+1
            FA(IO,IC) = FA(IO,IC)/SA(IC) * 100.
          ENDDO
          DO IO = 1,NO
            DO KE = 1,KM
              DO KB = 1,KM
                FB(KB,KE,IO,IC) = FB(KB,KE,IO,IC)/SA(IC) * 100.
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
!
C-----OUTPUT ASCII INFORMATION
      CLOSE(OU)
      FPATH =DIROUT(1:LENCHR(DIROUT))//TRIM(CASENM)//
     &       '_ALLCLOUDCELSS_FREQUENCY_f90.TXT'      ! ASCII
      OPEN(OU,FILE=FPATH(1:LENCHR(FPATH)) )
      DO IG = 1,NG
        WRITE(OU,'(A,I2,I6,A)')
     $         '======) GCM GRID ',IG,(SG(IG)+.01)
     $        ,' FREQUENCY OF ALL CLOUD CELLS'
        DO KE = KM,1,-1                          ! TOP TO SFC
          WRITE(OU,'(37E12.4)')
     $            (Z(KE)*10.),((FT(KB,KE,IG)*100.),KB=1,KM) !ft[kb,ke],the value of ft is cloud cover, 
!                                         # z[kb] is base,z[ke]is cloud top
        ENDDO
        WRITE(OU,'(A,34E12.4)') 'Z+>', ((Z(KB)*10.),KB=1,KM)
      ENDDO
      CLOSE(OU)
      FPATH =DIROUT(1:LENCHR(DIROUT))//TRIM(CASENM)//
     &       '_OVERLAPCLOUDCELSS_FREQUENCY_f90.TXT'      ! ASCII
      OPEN(OU,FILE=FPATH(1:LENCHR(FPATH)) )
      DO LC = 1,NE
        DO IO = 1,NO
          WRITE(OU,'(A,I2,I6,A)')
     $         '======) CEM CLASS',LC,(SA(LC)+.01)
     $        ,' FREQUENCY OF OVERLAP CLOUD CELLS'
          WRITE(OU,'(A,I2,A,10G12.5)') '....... NA=',IO
     $        ,' AREA= ',(FA(IC,LC),IC=0,NO+1)
          DO KE = KM,1,-1                          ! TOP TO SFC
            WRITE(OU,'(37E12.4)')
     $            (Z(KE)*10.),((FB(KB,KE,IO,LC)*100.),KB=1,KM)
          ENDDO
          WRITE(OU,'(A,34E12.4)') 'Z+>', ((Z(KB)*10.),KB=1,KM)
        ENDDO
      ENDDO
      CLOSE(OU)
      FPATH =DIROUT(1:LENCHR(DIROUT))//TRIM(CASENM)//
     &       '_DEEPCONVECTION_INFO_f90.TXT'      ! ASCII
      OPEN(OU,FILE=FPATH(1:LENCHR(FPATH)) )
      WRITE(OU,'(A)')
     $         '======) CEM DEEP CONVECTION ... SAMPLES'
      DO K = KM,1,-1                           ! TOP TO SFC
        WRITE(OU,'(21E12.4)')
     $            (Z(K)*10.), ((WD(K,KB)*.1+.01),KB=1,20)
      ENDDO
      WRITE(OU,'(A,20E12.4)') 'Z+> ',((Z(KB)*10.),KB=1,20)
      WRITE(OU,'(A)')
     $         '======) CEM DEEP CONVECTION ... POT TEMPERATURE'
      DO K = KM,1,-1                           ! TOP TO SFC
        WRITE(OU,'(21E12.4)')
     $            (Z(K)*10.), ((TD(K,KB)),KB=1,20)
      ENDDO
      WRITE(OU,'(A,20E12.4)') 'Z+> ',((Z(KB)*10.),KB=1,20)
      DO IP = 1,NP
        WRITE(OU,'(A,I2)')
     $         '======) CEM DEEP CONVECTION ... PROFILE ',IP
        DO K = KM,1,-1                           ! TOP TO SFC
          WRITE(OU,'(21E12.4)')
     $            (Z(K)*10.),((QD(K,KB,IP)*100.),KB=1,20)
        ENDDO
        WRITE(OU,'(A,20E12.4)') 'Z+> ',((Z(KB)*10.),KB=1,20)
      ENDDO
      CLOSE(OU)
      FPATH =DIROUT(1:LENCHR(DIROUT))//TRIM(CASENM)//
     &       '_PBLDEEPCONVECTION_INFO.TXT'      ! ASCII
      WRITE(OU,'(A)')
     $         '======) CEM PBL CONVECTION ... SAMPLES'
      DO K = KM,1,-1                           ! TOP TO SFC
        WRITE(OU,'(31E12.4)')
     $            (Z(K)*10.), ((WP(K,KE)*.1+.01),KE=1,30)
      ENDDO
      WRITE(OU,'(A,30E12.4)') 'Z+> ',((Z(KE)*10.),KE=1,30)
      DO IP = 1,NP
        WRITE(OU,'(A,I2)')
     $         '======) CEM PBL CONVECTION ... PROFILE ',IP
        DO K = KM,1,-1                           ! TOP TO SFC
          WRITE(OU,'(31E12.4)')
     $            (Z(K)*10.),((QP(K,KE,IP)*100.),KE=1,30)
        ENDDO
        WRITE(OU,'(A,30E12.4)') 'Z+> ',((Z(KE)*10.),KE=1,30)
      ENDDO
      CLOSE(OU)
      FPATH =DIROUT(1:LENCHR(DIROUT))//TRIM(CASENM)//
     &       '_SHALLOWCONVECTION_INFO_f90.TXT'      ! ASCII
      WRITE(OU,'(A)')
     $         '======) CEM SHALLOW CONVECTION ... PROFILES'
      WRITE(OU,'(A)')
     $         ' Z BASE  SAMPLE  PROFILES'
      DO KB = KM,1,-1                             ! TOP TO SFC
        WRITE(OU,'(10E12.4)')
     $         (Z(KB)*10.),(WS(KB)*.1+.01)
     $       , (TR(KB)),(TS(KB))
     $       ,((QS(KB,IP)*10000.),IP=1,NP)
      ENDDO
      CLOSE(OU)
C-----OUTPUT DATA FOR PLOTTING
      FPATH =DIROUT(1:LENCHR(DIROUT))//TRIM(CASENM)//'_PLOT-1D_f90.TXT'     ! ASCII 1D
      OPEN(O1,FILE=FPATH(1:LENCHR(FPATH)) )
      FPATH =DIROUT(1:LENCHR(DIROUT))//TRIM(CASENM)//'_PLOT-2D_f90.TXT'     ! ASCII 2D
      OPEN(O2,FILE=FPATH(1:LENCHR(FPATH)) )
      FPATH =DIROUT(1:LENCHR(DIROUT))//TRIM(CASENM)//'_R1DCCM3_f90.TXT'     ! ASCII FOR CCM3 R1D
      OPEN(O3,FILE=FPATH(1:LENCHR(FPATH)) )
!C-----OUTPUT BIN 1-D PLOTS
      DO IP = 1,NP
        DO KB = 1,KM
          QZD = 0.1*REAL(KB-1)
          DO K = 1,KM
            WW(K,KB,1) = (QD(K,KB,IP) + QZD)*100. ! QD IP=1 TOTAL WATER CONTENT
                                                  !    IP=2 
                                                  !    IP=3 
            IF (WD(K,KB).LE.0.) WW(K,KB,1) = SPV
          ENDDO
        ENDDO
        WRITE(O1,'(1X,4(1X,I4))') KM,KM*KM,KM*KM,KM
        WRITE(O1,99) ((WW(K,KB,1),K=1,KM),KB=1,KM)
     $            ,(( Z(K)     ,K=1,KM),KB=1,KM)
     $            , (               KM ,KB=1,KM)
      ENDDO

      DO KB = 1,KM
         TZD = 1.*REAL(KB-1)
         DO K = 1,KM
            WW(K,KB,1) = TD(K,KB) + TZD
            IF (WD(K,KB).LE.0.) WW(K,KB,1) = SPV
         ENDDO
      ENDDO
      WRITE(O1,'(1X,4(1X,I4))') KM,KM*KM,KM*KM,KM
      WRITE(O1,99) ((WW(K,KB,1),K=1,KM),KB=1,KM)
     $         ,(( Z(K)     ,K=1,KM),KB=1,KM)
     $         , (               KM ,KB=1,KM)

      DO IP = 1,NP
         DO KE = 1,KM
            QZD = 0.1*REAL(KE-1)
            DO K = 1,KM
               WW(K,KE,1) = (QP(K,KE,IP) + QZD)*100.
               IF (WP(K,KE).LE.0.) WW(K,KE,1) = SPV
            ENDDO
         ENDDO
         WRITE(O1,'(1X,4(1X,I4))') KM,KM*KM,KM*KM,KM
         WRITE(O1,99) ((WW(K,KB,1),K=1,KM),KB=1,KM)
     $            ,(( Z(K)     ,K=1,KM),KB=1,KM)
     $            , (               KM ,KB=1,KM)
      ENDDO

      DO IP = 1,NP
         DO KB = 1,KM
            CM(KB) = QS(KB,IP)
            IF (WS(KB).LE.0.) CM(KB) = SPV
         ENDDO
         WRITE(O1,'(1X,4(1X,I4))') 1,KM,KM,1
         WRITE(O1,99) (CM(KB),KB=1,KM),(Z(KB),KB=1,KM),KM
      ENDDO

      DO KB = 1,KM
         CM(KB) = TS(KB)
         IF (WS(KB).LE.0.) CM(KB) = SPV
      ENDDO
      WRITE(O1,'(1X,4(1X,I4))') 1,KM,KM,1
      WRITE(O1,99) (CM(KB),KB=1,KM),(Z(KB),KB=1,KM),KM

      WRITE(O1,'(1X,4(1X,I4))') 1,KM,KM,1
      WRITE(O1,99) (TR(K),K=1,KM),(Z(K),K=1,KM),KM

      DO KB = 1,KM
         CM(KB) = QZ(KB)
         IF (WZ(KB).LE.0.) CM(KB) = SPV
      ENDDO
      WRITE(O1,'(1X,4(1X,I4))') 1,KM,KM,1
      WRITE(O1,99) (CM(KB),KB=1,KM),(Z(KB),KB=1,KM),KM

      DO KE = 1,KM
         CM(KE) = QB(KE)
         IF (WB(KE).LE.0.) CM(KE) = SPV
      ENDDO
      WRITE(O1,'(1X,4(1X,I4))') 1,KM,KM,1
      WRITE(O1,99) (CM(KE),KE=1,KM),(Z(KE),KE=1,KM),KM

      CLOSE(O1)

C-----OUTPUT BIN 2-D PLOTS

      IF (ZBMIN.GT.0.0) THEN     ! DISCARD CLOUDD BELOW ZBMIN
         DO KB = 1,KM
            IF (Z(KB).LT.ZBMIN) THEN
               DO KE = 1,KM
                  DO IG = 1,NG
                     FT(KB,KE,IG) = 0.
                  ENDDO
                  DO LC = 1,NE
                     DO IO = 1,NO
                     FB(KB,KE,IO,LC) = 0.
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDIF

      NREC = 0
      DO IG = 1,NG
         NREC = NREC + 1
         WRITE(O2,99) ((FT(KB,KE,IG),KB=1,KM),KE=1,KM)
      ENDDO

      DO LC = 1,NE
         DO IO = 1,NO
            NREC = NREC + 1
            WRITE(O2,99) ((FB(KB,KE,IO,LC),KB=1,KM),KE=1,KM)
         ENDDO
      ENDDO

      DO IP = 1,NP
         DO KB = 1,KM
            IF (WD(KB,KB).LE.0.) THEN
                DO K = 1,KM
                   WW(K,KB,1) = SPV
                ENDDO
            ELSE
                DO K = 1,KM
                   WW(K,KB,1) = QD(K,KB,IP)
                ENDDO
            ENDIF
         ENDDO
         NREC = NREC + 1
         WRITE(O2,99) ((WW(K,KB,1),KB=1,KM),K=1,KM)
      ENDDO

      DO KB = 1,KM
         DO K = 1,KM
            WW(K,KB,1) = TD(K,KB)
            IF (WD(K,KB).LE.0.) WW(K,KB,1) = SPV
         ENDDO
      ENDDO
      NREC = NREC + 1
      WRITE(O2,99) ((WW(K,KB,1),KB=1,KM),K=1,KM)

      DO IP = 1,NP
         DO KE = 1,KM
            IF (WP(KE,KE).LE.0.) THEN
                DO K = 1,KM
                   WW(K,KE,1) = SPV
                ENDDO
            ELSE
                DO K = 1,KM
                   WW(K,KE,1) = QP(K,KE,IP)
                ENDDO
            ENDIF
         ENDDO
         NREC = NREC + 1
         WRITE(O2,99) ((WW(K,KE,1),KE=1,KM),K=1,KM)
      ENDDO

      CLOSE(O2)

C-----OUTPUT BIN CCM3 R1D CALCULATIONS
      DO KB = 1,KM
         CM(KB) = QZ(KB)
         IF (WZ(KB).LE.0.) CM(KB) = SPV
      ENDDO
      IP = 2
      DO KB = 1,KM
         DO K = 1,KM
            IF (WD(K,KB).LE.0.) THEN
                WW(K,KB,1) = SPV
            ELSE
                WW(K,KB,1) = QD(K,KB,IP)
            ENDIF
         ENDDO
      ENDDO
      WRITE(O3,99) ((WW(K,KB,1),KB=1,KM),K=1,KM)
     $          ,(CM(KB),KB=1,KM)
     $          ,( Z(K),K=1,KM)
     $          ,( D(K),K=1,KM)
     $          ,(TR(K),K=1,KM)
     $          , SPV

      CLOSE(O3)
      ENDDO
99    FORMAT(8E12.4)  

      END PROGRAM
C     ------------------------------
C     ------------------------------
      SUBROUTINE DEEPCC(QW,IM,KM,OU,IT) ! QW TOTAL CLOUD WATER,IM: HORIZONTAL GRID; KM : VERTICAL GRID
C     ------------------------------
C     ------------------------------
C
C     OUTPUT DEEP CONVECTION
C
      INTEGER IM,KM,OU,I,K,L,KB(99),KE(99),NA
      REAL QW(IM,KM),C(99),CM(99)
      INTEGER IDC

!      WRITE(OU,'(A3,1X,I4,1X,\)')'ITT', IT !  \ NOT CHANGE LINE
      DO I = 1,IM
        DO K = 1,KM
          IF (QW(I,K).GE.1.0E-3) C(K) = 1. ! THIS GRIG IS COVERED BY CLOUD
          IF (QW(I,K).LT.1.0E-3) C(K) = 0.
        ENDDO
        CALL INFCLD(C,KM,KB,KE,CM,NA)
        IDC=0
        DO L = 1,NA
          IF (KB(L).LE.5 .AND. KE(L).GE.24) THEN  !! 5 AND 24 ARE THE VERTICAL LEVELS
          IDC=IDC+1
          ENDIF
        ENDDO
        IF(IDC>0) THEN
          WRITE(OU,'(I2,1X,I2)') NA, IDC!
          IF (IDC>1) THEN
            PRINT*,IDC
          ENDIF 
        ELSE
          WRITE(OU,'(I2,1X,I2)') NA, -9 !  
        ENDIF
        DO L = 1,NA
          IF (KB(L).LE.5 .AND. KE(L).GE.24)  !! 5 AND 24 ARE THE VERTICAL LEVELS
     $      WRITE(OU,'(34e12.4,2I4)') ((C(K)+0.1),K=1,KM),KB(L),KE(L)  !!! DEEP CONVECTION CLOUD BASE AND TOP
          IDC=IDC+1
        ENDDO
      ENDDO
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
C     -------------------
C     -------------------
      FUNCTION LENCHR(CH)
C     -------------------
C     -------------------
C
C     RETURNS POSITION OF LAST NON-BLANK CHARACTER IN CH,
C     OR 1 IF CH IS ALL BLANKS.
C
      CHARACTER*(*) CH
      DO 10 I=LEN(CH),1,-1
         IF (CH(I:I).NE.' ') THEN
            LENCHR = I
            RETURN
         ENDIF
   10 CONTINUE
      LENCHR = 1
      RETURN
      END
