      PROGRAM SENN              ! VERSION 11.1

C                           REFERENCE                                  *
C  THE "SPATIALLY EXTENDED NONLINEAR NODE" (SENN) MODEL FOR THE        *
C  EXCITATION OF A MYELINATED NERVE IS ADAPTED FROM A PAPER ENTITLED   *
C  "ANALYSIS OF A MODEL FOR EXCITATION OF MYELINATED NERVE" BY DONALD  *
C  R. McNEAL, IEEE TRANSACTIONS ON BIOMEDICAL ENGINEERING, BME-23,     *
C  NO. 4, JULY 1976, pp 329-337.                                       *
C----------------------------------------------------------------------*
C  A DESCRIPTION OF THE SENN MODEL AND A COMPREHENSIVE USER'S GUIDE    *
C  ARE PRESENTED IN THE BOOK "ELECTROSTIMULATION: THEORY, APPLICATIONS *
C  AND COMPUTER MODEL" BY REILLY & DIAMANT, ARTECH HOUSE PUBLISHERS,   *
C  2011. THE LATEST MODEL SOURCE CODE, PC AND MAC EXECUTABLES, A       *
C  STARTUP GUIDE, AND SUPPORTING FILES ARE AVAILABLE FREE OF COPYRIGHT *
C  RESTRICTION ON THE PUBLISHER'S WEBSITE AT                           *
C     http://www.artechhouse.com/static/reslib/reilly/reilly.html      *
C----------------------------------------------------------------------*
C  *IMPORTANT*: FOR PROPER VIEWING AND PRINTING ALIGNMENT OF MODEL     *
C  SOURCE CODE AND OUTPUT FILES, USE THE COMPILER EDITOR OR A          *
C  MONOSPACED FONT SUCH AS COURIER NEW, COURIER PS, LUCIDA CONSOLE,    *
C  OR MONACO.                                                          *
C***********************************************************************
      include "COMMON.f"
C     WARNING
C          @@@@@@   PNAB=(GM/30.365)*PNAB    @@@@@@@
C                    @@@@@@@@@@@@@@@@@@@@
C   EPOTIN, SINEIN, XIN and YIN are filled from external EXCEL tab-delimited text files 
C       EPOTFILE.txt, SINEFILE.txt and XYFILE.txt. 
C   SINEIN2 is computed from SINEIN, substituting radian frequencies and phases.
C   NSINES is the number of sine wave parameter sets in SINEIN, currently up to 1100.

C   XCAL and YCAL are computed from XIN and YIN and are interpolated points for the 
C       Runge-Kutta subroutine RKGS when IWAVE = 13. The full interpolated Y array is 
C       YINTERP. The full interpolated (X,Y) array is in output file XYINTERP.
C   NTRP is the number of points to be interpolated between consecutive points in 
C       XIN,YIN for the formation of XCAL,YCAL and output file XYINTERP. 
C   Dimensioning is currently set for a maximum of 8000 input points and a maximum of 
C       299 interpolated points (between consecutive input points) in IWAVE = 13. 

      NAMELIST/FIBER/
     &   NNODES,NLIN1,NLIN2,NODE1,DIAM,GAP,CM,GM,RHOI,RHOE
      NAMELIST/STIMULUS/XC,YC,XA,YA,WIREL,
     &   IWAVE,UIO,XPD,UIO2,XPD2,DELAY,FREQ,PHASE,FREQ2,
     &   PHASE2,AMP2,NSINES,DCOFF,TAUS,VREF,NP,FS,S,NTRP 
         NAMELIST /CONTROL/
     &   TT,ITHR,VTH,NTHNODE,DELT,DELT2,DELT2M,FINAL,IPRNT,pltn
      NAMELIST/HILLOCK/WB,WH,DIAMB,DIAMH
        character*20 namlst

      OPEN(UNIT=66,FILE='data.out',FORM='FORMATTED',
     &STATUS = 'UNKNOWN')  ! L.S. 'NEW' changed to Absoft 'UNKNOWN' 12/4/01

C CONSTANTS FOR THE IONIC COMPONENT EQUATIONS OF TRANSMEMBRANE CURRENT
C AS GIVEN IN McNEAL (1976, p 336) 

      IRUN=0 !RUN COUNTER

      open(unit=7,file='inparam.txt',form='formatted',status='old',
     &access='sequential')  !!  6/6/07 

C      OPEN(UNIT=37,FILE='X&XMULT',STATUS='NEW') !! Diagnostic plot file, 
C                                   time series (X,XMULT).  See subr. FCT.
C      OPEN(UNIT=38,FILE='X&XMULT&I',STATUS='NEW') !! Diagnostic file, 
C                        (X,XMULT) & interp out array Index.  See subr. FCT.

3333    CONTINUE

      CROSS = .FALSE.  !! Resetting logical variables for run concatenation  10/16/06
      NFOUND = .FALSE.
      DO I = 1,6
        TFLAG(I) = .FALSE.
      ENDDO
C      Write(*,*) 'nfound = ',nfound, 'cross = ',cross, 'tflag = ',tflag  ! diagnostic 10/16/06
      
      R=8.3144   !gas constant
      T=295.16   !295.18K absolute temp McNeal
      F=96.487   !96514.0 C/g/mole Faraday's constant McNeal
      SODO=114.5 !(Na)o external sodium concentration
      SODI=13.74 !(Na)i internal sodium concentration
      POTO=2.5   !(K)o, external potassium concentration
      POTI=120.  !(K)i, internal potassium concentration
      PNAB=.008  !sodium permeability constant
      PKB=.0012  !potassium permeability constant
      PPB=.00054 !nonspecific permeability constant
      GL=30.3    !leak conductance
      VL=.0260430075 !.026mV, leak current equilibrium potential
      ER=-70.    !resting potential Vr
      TT=1.      !1= terminated axon, 2= cell body + Hillock
      wb=0.
      wh=0.
      diamb=0.
      diamh=0. !! Corrected HILLOCK typo; "diah" changed to "diamh" 3/27/2010

C MORE IONIC CURRENT PARAMETERS AS GIVEN IN McNEAL (1976, p 336) 
C THE PARAMETERS ARE USED IN THE ALPHA SUB X AND BETA SUB X EXPRESSIONS
C  (X=H,M,P,N) AS FOLLOWS:

C    ALPHA(X)=CA(X,1)(V-Y)(1-EXP((CA(X,2)-V)/CA(X,3)))**-1
C     SAME FOR BETA(X)   REPLACE CA BY BA
C      EQUATIONS FORMED IN FCT WITH A(X)=ALPHA(X)  B(X)=BETA(X)

C  IN ABOVE X=1,2,3,4 FOR H,M,P,N RESPECTIVELY
      CA(1,1)=0.1  !coeff alpha h
      CB(1,1)=4.5  !coeff beta h
      CA(2,1)=0.36 !coeff alpha m
      CB(2,1)=0.4  !coeff beta m
      CA(3,1)=0.006!coeff alpha p
      CB(3,1)=0.09 !coeff beta p
      CA(4,1)=0.02 !coeff alpha n
      CB(4,1)=0.05 !coeff beta n
      CA(1,2)=-10. !term alpha h
      CB(1,2)=45.  !term beta h
      CA(2,2)=22.  !term alpha m
      CB(2,2)=13.  !term beta m
      CA(3,2)=40.  !term alpha p
      CB(3,2)=-25. !term beta p
      CA(4,2)=35.  !term alpha n
      CB(4,2)=10.  !term beta n
      CA(1,3)=6.   !denom term alpha h
      CB(1,3)=10.  !denom term beta h
      CA(2,3)=3.   !denom term alpha m
      CB(2,3)=20.  !denom term beta m
      CA(3,3)=10.  !denom term alpha p
      CB(3,3)=20.  !denom term beta p
      CA(4,3)=10.  !denom term alpha n
      CB(4,3)=10.  !denom term beta n
      CM=2.        !membrane capacitance/unit area (Overridden in inparam.txt)
      RHOI=110.    !resistivity of axoplasm        (Overridden in inparam.txt)
      GM=30.365    !membrane conductance/unit area (Overridden in inparam.txt)
      GAP=.00025   !gap width                      (Overridden in inparam.txt)
      ELD=100.     !ratio of internodal space to fiber diameter L/DIAM
      SDD=0.7      !ratio of axon and fiber diameters
      NNGTT=0      !# nodes found gt threshold

C  HANDY PI EXPRESSIONS
      DATA PI/3.141593/
      TWOPI = 2.0*PI
      FOURPI = 4.0*PI
      PID180 = PI/180.0

C      READ INPUT VARIABLES

    1 CONTINUE
      READ(7,6666,END=5000)MES2
6666   FORMAT(20A4)
      READ(7,FIBER)
      READ(7,STIMULUS)

C   TRAPS FOR IMPROPER VALUES OF FS AND S
      IF(FS .GT. 3 .OR. FS .LT. 0)THEN
        WRITE(66,*) 'FS must be integer 0, 1, 2 or 3'
        WRITE(*,*) 'FS must be integer 0, 1, 2 or 3'
        STOP
      ENDIF
      IF(S .GT. 1 .OR. S .LT. 0)THEN
        WRITE(66,*) 'S Must be integer 0 or 1'
        WRITE(*,*) 'S Must be integer 0 or 1'
        STOP
      ENDIF
C   WARNING ABOUT USE OF FS=2
      IF(FS .EQ. 2 .AND. (IWAVE .NE. 1 .OR. NP .NE. 1)) THEN
        WRITE(66,*) 'FS = 2 VERIFIED ONLY FOR IWAVE = 1, SINGLE PULSE'
        WRITE(*,*) 'FS = 2 VERIFIED ONLY FOR IWAVE = 1, SINGLE PULSE'
        WRITE(*,*) 'Press Return to Continue, or Terminate Program'
        PAUSE 
      ENDIF 

      WRITE(*,600)
      WRITE(66,600)
  600 FORMAT(' SPATIALLY EXTENDED NON-LINEAR NODE (SENN) MODEL, 2010')
      IF(FS .EQ. 1 .OR. FS .EQ. 2)THEN
        WRITE(*,*)'IWAVE',IWAVE, '   FS',FS
        WRITE(66,*)'IWAVE',IWAVE, '   FS',FS
      ELSE  ! FS = 0 or 3
        WRITE(*,*)'IWAVE',IWAVE, '   FS',FS, '   S',S
        WRITE(66,*)'IWAVE',IWAVE, '   FS',FS, '   S',S
      ENDIF

C   PIMULT VALUES FOR ELECTRODE PROBE MODES
      IF(FS .EQ. 0 .AND. S .EQ. 0)THEN
        PIMULT=FOURPI    ! Point electrode in situ
      ELSEIF(FS .EQ. 0 .AND. S .EQ. 1)THEN
        PIMULT=TWOPI     ! Point electrode on surface
      ELSEIF(FS .EQ. 3 .AND. S .EQ. 0)THEN
        PIMULT=TWOPI     ! Line electrode in situ
      ELSEIF(FS .EQ. 3 .AND. S .EQ. 1)THEN
        PIMULT=PI        ! Line electrode on surface
      ENDIF             ! FS=0 or FS=3 setup of PIMULT

      IF(IWAVE .EQ. 13) THEN   !Input waveform
        WRITE(66,*) 'Note: In this implementation, X values for',
     &' IWAVE=13 are internally '
        WRITE(66,*) 'calculated, based on one space measurement',
     &' of the input data. '
        WRITE(66,*) 'FOR VALID RESULTS, THE INPUT DATA POINTS',
     &' MUST BE UNIFORMLY SPACED.'
        WRITE(*,*) 'Note: In this implementation, X values for',
     &' IWAVE=13 are internally '
        WRITE(*,*) 'calculated, based on one space measurement',
     &' of the input data. '
        WRITE(*,*) 'FOR VALID RESULTS, THE INPUT DATA POINTS',
     &' MUST BE UNIFORMLY SPACED.'
        IF (NTRP .GT. 299) THEN
          WRITE(66,*) '299 interpolated points is the current',
     &' maximum.'
          WRITE(*,*) '299 interpolated points is the current maximum.'
          STOP
        ENDIF
      ENDIF
	  
      IF(IWAVE .EQ. 8.OR. IWAVE .EQ. 9)THEN
        WRITE(*,*)'UIO2,FREQ, FOR IWAVE?',IWAVE
        READ(*,*)UIO2,FREQ
      ENDIF
      if(iwave .eq. 3)then
        write(*,*)'UIO,FREQ,Tp,DCOFF FOR IWAVE?',IWAVE
        read(*,*)uio,FREQ,Tp,DCOFF
      endif                           
      READ(7,CONTROL)
      DELT2 = DELT2M*DELT  !! 9/4/10 Inp param DELT2M to form DELT2
      IF(TT.EQ. 2)THEN
        READ(7,HILLOCK)
C       CLOSE(7)
      ENDIF

C     Import array(s) for certain modes. Presently import EPOT array for FS = 2, 
C     Sine Wave parameter array for IWAVE = 12, and X,Y arrays for IWAVE = 13.
      IF (FS .EQ. 2) THEN   ! Array of external potentials
        OPEN(UNIT=22,FILE='EPOTfile.txt',STATUS='OLD')
        READ(22,*,END=100) EPOTIN
  100   CLOSE(22)
        WRITE(*,*) (EPOTIN(I),I=1,NNODES)
      ENDIF

      IF (IWAVE .EQ. 12) THEN   ! Array of sine wave parameters
        OPEN(UNIT=22,FILE='SINEfile.txt',STATUS='OLD')
        READ(22,*,END=101) ((SINEIN(I,J),J=1,3),I=1,NSINES)
  101   CLOSE(22)
      ENDIF

      IF (IWAVE .EQ. 13) THEN   ! Array of waveform (x,y) values
        OPEN(UNIT=2,FILE='XYINTERP',STATUS='UNKNOWN') ! Interpolated (x,y) file,
C                   see subr. INTERP ! L.S. 'NEW' changed to Absoft 'UNKNOWN' 12/4/01
        OPEN(UNIT=22,FILE='XYfile.txt',STATUS='OLD')
        OPEN(UNIT=43,FILE='XYIN',STATUS='UNKNOWN')     !! Diagnostic
C                   see subr. INTERP ! L.S. 'NEW' changed to Absoft 'UNKNOWN' 5/30/02

        WRITE(66,*) 'External input file named XYFILE'
        WRITE(*,*) 'External input file named XYFILE'
        LENIN = 1
        DO WHILE(.TRUE.)
C		  WRITE (*,*) LENIN    !! Diagnostic
          READ(22,*,END=102) XIN(LENIN),YIN(LENIN)
          WRITE(43,*) XIN(LENIN),YIN(LENIN)    !! Diagnostic
          LENIN = LENIN+1
        END DO
  102   CLOSE(22)
        CLOSE(43)    !! Diagnostic
        LENIN = LENIN-1
C        WRITE(*,*) LENIN   !! Diagnostic
C        WRITE(*,*) XIN(27),YIN(27)   !! Diagnostic
        IF (LENIN .GT. 8001) THEN
          WRITE(66,*) '8001 input points is the current maximum.'
          WRITE(*,*) '8001 input points is the current maximum.'
          STOP
        ENDIF
        WRITE(66,*) LENIN, ' input (x,y) pairs'
        WRITE(*,*) LENIN, ' input (x,y) pairs'

C       NOTE: IN THIS IMPLEMENTATION, INPUT TEMPORAL SPACING MUST BE IN MILLISECONDS
        DELTIN = XIN(2)-XIN(1)    ! spacing in milliseconds
        WRITE(*,*) DELTIN,' ms input x-spacing'
        WRITE(*,*) NTRP, ' points interpolated'
        CALL INTERP(XIN,YIN,DELTIN,LENIN,NTRP,XCAL,YCAL,YINTERP,
     &DELTOT,LENOT)
        CLOSE(2)  !! close XYINTERP 1/29/10
        DELT = DELTOT*2.0   ! Overrides inp. param. list. Mult. by 2 for RKGS.
        DELT2 = DELT2M*DELT    !! 9/4/10
        WRITE(*,*) DELTOT, ' ms output (interpolated) spacing'
        WRITE(*,*) LENOT, ' output (interpolated) (x,y) pairs'
        WRITE(*,*) 'DELT = ',DELT, ' ms overriding inp. param.',
     &' list (interp spacing*2)'
        WRITE(*,*) 'DELT2 = ',DELT2, ' ms = DELT2M*DELT'
      ENDIF

      WRITE(*,*) 'NNODES', NNODES

      INNGTT=NTHNODE          !# of nodes needed to validate 
        NON = (NNODES-1)/2    !# of nodes to the rt or left of center node
      URATIO=1.
C      WRITE (*,*) 'URATIO = ',URATIO  !! DIAGNOSTIC 10/19/03

      IF(IWAVE .EQ. 6)THEN
        NP=2
      ENDIF
C     FIX TO ASSURE UIO2=0.0 WHEN NP=1
      IF(NP .EQ. 1.AND. IWAVE .NE. 8 .AND.IWAVE.NE.9)THEN
        UIO2=0.
        DELAY=0.
        XPD2=0.    !! Added 8/25/09
      ELSEIF(NP.GE. 2.AND. IWAVE .NE. 6)THEN
        UIO2=UIO
        XPD2=XPD
      ENDIF

C   Form radian angles and radian frequencies 

      IF(IWAVE .EQ. 2. OR. IWAVE .EQ. 3. .OR. IWAVE .EQ. 5
     & .OR. IWAVE .EQ. 10 .OR. IWAVE .EQ. 11) THEN
        ANGLE=PHASE*PID180
        PROD =TWOPI*FREQ
      ENDIF

      IF (IWAVE .EQ. 10 .OR. IWAVE .EQ. 11) THEN    ! SECOND SINUSOID
        ANGLE2=PHASE2*PID180
        PROD2 =TWOPI*FREQ2
      ENDIF

      IF (IWAVE .EQ. 12) THEN     ! MULTIPLE SINUSOIDS
        DO I = 1,NSINES
          SINEIN2(I,1) = SINEIN(I,1)         ! AMPLITUDE
          SINEIN2(I,2) = SINEIN(I,2)*TWOPI   ! RADIAN FREQUENCY
          SINEIN2(I,3) = SINEIN(I,3)*PID180  ! RADIAN PHASE
        END DO
C        WRITE(*,10) ((SINEIN2(I,J),J=1,3),I=1,NSINES)  !! Diagnostic
      ENDIF

      IF(IWAVE .EQ. 8)THEN
        PROD =TWOPI*FREQ
        angle = PID180*phase-TWOPI*freq*xpd
      ENDIF
      IF(IWAVE .EQ. 9)THEN
        PROD = TWOPI*FREQ
        ANGLE=PHASE*PID180
      ENDIF

      WRITE(66,6667)MES2  !! Run descriptor, top of inparam.txt
6667   FORMAT(20A4)      
      write(66,*) 
        namlst='Namelist : FIBER'
      write(66,*)namlst
      WRITE(66,*)
      write(66,*)'NNODES=',NNODES,
     &'         TOTAL NUMBER OF NODES (Odd)'
        WRITE(66,*)
     &'NLIN1= ',NLIN1,'          FIRST NONLINEAR NODE'
        WRITE(66,*)
     &'NLIN2= ',NLIN2,'         LAST NONLINEAR NODE'
       WRITE(66,*)
     &'NODE1= ',NODE1,'         FIRST PRINT NODE'
        WRITE(66,95)
     & ' DIAM=',DIAM,'         FIBER DIAMETER (cm)'
90      FORMAT(1A,F7.3,1A)
        WRITE(66,95)
     &  ' GAP= ',GAP,'         INTRANODAL GAP (cm)'
95      FORMAT(1A,F7.5,1A)
        WRITE(66,90)' CM=  ',CM,
     &'         MEMBRANE CAPACITY (uF/cm**2)'
        WRITE(66,90)
     &  ' GM=  ',GM,'         LIN. MEMBR. CONDUCTANCE/AREA (mS/cm**2)'
        WRITE(66,90)
     &' RHOI=',RHOI,'         AXOPLASM RESISTIVITY (ohm.cm)'
92      FORMAT(1A,F7.3,1A)
        WRITE(66,90)
     & ' RHOE=',RHOE,'         MEDIUM RESISTIVITY (ohm.cm)'
        namlst='Namelist : STIMULUS'
      WRITE(66,*)
      write(66,*)namlst
      WRITE(66,*)
      write(66,90)' XC=',XC,
     &'           X LOCUS CATHODE FOR FS=0 or FS=3 (cm)'
        WRITE(66,90)' YC=',YC,
     &'           Y LOCUS CATHODE FOR FS=0 or FS=3 (cm)'
        WRITE(66,90)
     &  ' XA=',XA,
     &'           X LOCUS ANODE FOR FS=0 or FS=3 (cm)'
        WRITE(66,90)' YA=',YA,
     &'           Y LOCUS ANODE FOR FS=0 or FS=3 (cm)'
        WRITE(66,*)'WIREL=',WIREL,
     &'    PROBE WIRE CONTACT LENGTH FOR FS=3 (cm)'

      IF(IWAVE .EQ. 1)THEN
        WRITE(66,*)'IWAVE=',IWAVE,'           MONOPHASIC'
        ELSEIF(IWAVE .EQ. 2)THEN
          WRITE(66,*)'IWAVE=',IWAVE,'          SINEWAVE'
        ELSEIF(IWAVE .EQ. 3)THEN
          WRITE(66,*)'IWAVE=',IWAVE,'          SINE + DC'
        ELSEIF(IWAVE .EQ. 4)THEN
          WRITE(66,*)'IWAVE=',IWAVE,'          EXPONENTIAL'
        ELSEIF(IWAVE .EQ. 5)THEN
          WRITE(66,*)'IWAVE=',IWAVE,'          SINE*EXPONENTIAL'
        ELSEIF(IWAVE .EQ. 6 .AND. UIO/UIO2 .GE. 0.0)THEN
          WRITE(66,*)'IWAVE=',IWAVE,'          MONOPHASIC DOUBLET'
        ELSEIF(IWAVE .EQ. 6 .AND. UIO/UIO2 .LT. 0.0)THEN
          WRITE(66,*)'IWAVE=',IWAVE,'          BIPHASIC DOUBLET'
        ELSEIF(IWAVE .EQ. 7)THEN
          WRITE(66,*)'IWAVE=',IWAVE,'          SPECIAL INPUT'
        ELSEIF(IWAVE .EQ. 8)THEN
          WRITE(66,*)'IWAVE=',IWAVE,
     &'          RECTANGULAR PULSE + SINUSOIDAL'
        ELSEIF(IWAVE .EQ. 9)THEN
          WRITE(66,*)'IWAVE=',IWAVE,
     &'          SINUSOIDAL + RECTANGULAR PULSE'
     
        ELSEIF(IWAVE .EQ. 10)THEN
          write(66,*)'IWAVE=',IWAVE,'          SUM OF TWO SINUSOIDS'
        ELSEIF(IWAVE .EQ. 11)THEN
          write(66,*)'IWAVE=',IWAVE,
     &'          COSINE AMPLITUDE MODULATION'
        ELSEIF(IWAVE .EQ. 12)THEN
          WRITE(66,*)'IWAVE=',IWAVE,
     &'          SUM OF ARRAY OF SINUSOIDS'
        ELSEIF(IWAVE .EQ. 13)THEN
          WRITE(66,*)'IWAVE=',IWAVE,'          INPUT WAVEFORM ARRAY'

      ELSE
        WRITE(66,*)'IWAVE=',IWAVE,'          UNDEFINED WAVEFORM_'
      STOP
      ENDIF	
        WRITE(66,90)
     &' UIO= ',UIO,
     &'         STIMULUS AMPLITUDE(mA for FS=0 or FS=3,',
     &'                     mA/cm**2 for FS=1)'
        WRITE(66,90)
     &' XPD= ',XPD,'         PULSE DURATION (ms)'
        WRITE(66,90)
     &' UIO2=',UIO2,'         2ND (DOUBLET)AMPLITUDE (mA or mA/cm**2)'
        WRITE(66,90)
     &' XPD2=',XPD2,'         2ND (DOUBLET) PULSE DURATION (ms)'
        WRITE(66,99)
     &' DELAY=',DELAY,'         DELAY BETWEEN PULSES,IWAVE=1;',
     &'                     or BETWEEN DOUBLET PHASES,IWAVE=6 (ms)'
        WRITE(66,201)
     &' FREQ= ',FREQ,'     SINEWAVE FREQUENCY,(kHz),',
     &'                     for IWAVE=2,3,5,8,9,10,11' 
201     FORMAT(1A,F10.3,1A) 
        WRITE(66,199)
     &' PHASE= ',PHASE,'        SINEWAVE PHASE,(deg),',
     &'                     for IWAVE=2,3,5,8,9,10,11'

      WRITE(66,201)
     &' FREQ2= ',FREQ2,'    SECOND SINEWAVE FREQUENCY,(kHz),',
     &'                     for IWAVE=10,11'
      WRITE(66,199)
     &' PHASE2= ',PHASE2,'       SECOND SINEWAVE PHASE,(deg),',
     &'                     for IWAVE=10,11'
      WRITE(66,90)
     &' AMP2= ',AMP2,'        SECOND SINEWAVE AMPLITUDE,',
     &'                     for IWAVE=10,11'
      WRITE(66,*)'NSINES=',NSINES,
     &'         NUMBER OF SINE WAVES FOR IWAVE=12'

199      FORMAT(1A,F6.1,1A)
        WRITE(66,99)
     &' DCOFF=',DCOFF,'         DC OFFSET, IWAVE=3'
       WRITE(66,90)
     &' TAUS=',TAUS,'         EXP TIME CONSTANT FOR IWAVE=4,5'
       WRITE(66,90)
     &' VREF=',VREF,'         FIRST NODE POTENTIAL FOR FS=1'
       WRITE(66,*)
     &'NP=',NP,'              NUMBER OF PULSE REPETITIONS, IWAVE=1'
       WRITE(66,*)
     &'FS=',FS,'              0=POINT ELECTRODE, 1=UNIFORM FIELD,'
       WRITE(66,*)
     &'                    2=IMPORT EPOT ARRAY, 3=WIRE ELECTRODE'
       WRITE(66,*)
     &'S=',S,'               0=ELECTRODE IN SITU,'
       WRITE(66,*)
     &'                    1=ELECTRODE ON SURFACE'
       WRITE(66,*)
     &'NTRP=',NTRP,'            NUMBER OF INTERPOLATED POINTS',
     &' FOR IWAVE=13'
       WRITE(66,*)
C       WRITE(66,*)
        namlst='Namelist : CONTROL'
       write(66,*)namlst
       WRITE(66,*)
       WRITE(66,*)'TT=',TT,'        1= TRUNCATED AXON,',
     &' 2= CELL BODY+HILLOCK'
       WRITE(66,*)
     &'ITHR=   ',ITHR,
     &'         0=SINGLE RUN,1=THRESHOLD SEEKING'
       WRITE(66,90)
     &' VTH=',VTH,
     &'          VOLTAGE CRITERIA FOR THRESHOLD(mV)'
       WRITE(66,*)
     &'NTHNODE=',NTHNODE,'         NO. NODES THRESHOLD FOR ITHR=1'
       WRITE(66,98)
     &' DELT=  ',DELT,'         TIME STEP DURING PULSE (ms)'
97      FORMAT(1A,F6.2,1A)   ! 9/10/10
98      FORMAT(1A,F5.4,1A)
      IF(DELT2. EQ. 0)DELT2=DELT
       WRITE(66,97)
     &' DELT2M= ',DELT2M,'       FACTOR FOR TIME STEP AFTER PULSE'
       WRITE(66,99)
     &' FINAL=',FINAL,'         MAX. SOLUTION TIME/RUN (ms)'
99      FORMAT(1A,F6.3,1A)
       WRITE(66,*)
     &'IPRNT=  ',IPRNT,'        PRINT INTERVAL'
      WRITE(66,*)
        namlst='Namelist : HILLOCK'
      write(66,*)namlst
       WRITE(66,98)
     &' WB= ',WB,'            WIDTH OF CELL BODY (cm)'
       WRITE(66,98)
     &' WH= ',WH,'            WIDTH OF CELL HILLOCK(cm)'
       WRITE(66,98)
     &' DIAMB= ',DIAMB,'         DIAMETER OF CELL BODY (cm)'
        write(66,98)
     &' DIAMH= ',DIAMH,'         DIAMETER OF CELL BODY (cm)'
      WRITE(66,*)
      WRITE(66,666)
666    FORMAT('1')
      IF(IWAVE .NE. 8.AND.IWAVE .NE. 9.AND.UIO.NE. 0.0) URATIO=UIO2/UIO
C       WRITE (*,*) 'URATIO = ',URATIO  !! DIAGNOSTIC 10/19/03

C COMPUTE TIMES OF PULSE LEAD AND PULSE TRAIL
C      IF(XPD2 .EQ. 0.0)THEN
      IF(IWAVE .NE. 6)THEN  !!Temp. fix for IWAVES 1 & 6  1/25/10
        DO I=1,NP 
          PT(I)= I*XPD  +(I-1)*DELAY
          PL(I)=PT(I)-XPD
          write(*,*)'pl',pl(i),' pt',pt(i),' dly',delay
        ENDDO
      ELSE       !! IWAVE=6 temp. fix 1/25/10
        DO I = 1,NP-1  !! 1/16/2010
          PT(I)= I*XPD+(I-1)*DELAY !time of trailing edge of stimuli +delay
          PL(I)=PT(I)-XPD !time of leading edge of stimulus
          write(*,*)'pl',pl(i),' pt',pt(i),' dly',delay
C          PT(I+1)=DELAY + (I+1)*XPD2  !! old
          PT(I+1)=DELAY + XPD+I*XPD2  !! 1/16/2010
          PL(I+1)=PT(I+1)-XPD2
          write(*,*)'pln',pl(i+1),' pt',pt(i+1),' dly',delay
        ENDDO
      ENDIF


C VARY IONIC NON LINEAR PARAMETERS DEPENDING ON GM

      XMFACT=GM/30.365
      PNAB=XMFACT*PNAB


C      WRITE(66,600)
C      WRITE(*,600)
C  600 FORMAT(' REVISED INPUT NEUROELECTRIC MODEL 24 NOV 92')


C      IF(UIO.LT.0.) INNGTT=2
C      UIOP=.80*UIO

C NORMALIZE YC TO SATISFIY Y0  BY YC=Y0*DIAM*ELD
C 2-5-91 INPUT CHANGED TO ACCEPT YC INSTEAD OF Y0 
C 2-5-91 NOW SOLVE FOR Y0 INSTEAD OF YC
C      YC=Y0*DIAM*ELD
      Y0=YC/(DIAM*ELD)     !! Not used!?  10/15/06

C      INITIALIZE VARIABLES AND CONSTANTS

      PRMT(1)=0.0
      PRMT(2)= FINAL
      TEND=FINAL	
      IF(ITHR.EQ.1) GO TO 408 !threshold seeking  !! 8/24/07 TEST CHANGE from 407 to 408
      IF(ITHR .EQ. 0) GO TO 408 !non-threshold seeking

C RUN OUT MODE modified by stimulus width under certain circumstances

  407 TEND=XPD +0.5
      IF(XPD .GE. 0.1 .AND. XPD .LE. 0.5) TEND= TEND + 0.47
      IF(XPD .GE. 1.0 .AND. XPD .LT. 1.5) TEND=XPD +0.3
      IF(XPD .GE. 2.0) TEND=XPD + 0.2   !!  Test (Old TEND=XPD)
  408 CONTINUE
      PRMT(3)=DELT
      PRMT(4)=100.
      UIOLD=UIO
      IT=0
      IF(IPRNT.EQ.0) IPRNT=.02/DELT+.5  !modifies printing interval
      NI=0
      ITA=0 
      ITB=0
      EL=ELD*DIAM
      NDIM=2*NON+4*(NLIN2-NLIN1)+4
      CGA=1000.*PI*DIAM*SDD**2/(4.*RHOI*ELD) !AXIAL INTERNODAL CONDUCTANCE
      AREA=PI*SDD*DIAM*GAP      !INTRANODAL SURFACE AREA
      WRITE(*,*)'AREA',AREA

C      WRITE(*,*)'IWAVE',IWAVE

      AB=PI*DIAMB*WB  !SURFACE AREA OF CELL BODY
      AH=PI*DIAMH*WH  !SURFACE AREA OF HILLOCK
      LBH=(WB+WH)/2. !WB=WIDTH OF CELL BODY WH=WIDTH OF HILLOCK
      LHN=(WH+GAP)/2.!GAP=NODAL GAP = WN
      AN= AREA       !INTRANODAL SURFACE AREA 
      CCM=CM*AN      !MEMBRANE CAPACITANCE 
      CGM=GM*AN      !MEMBRANE CONDUCTANCE
      IF(TT .EQ. 2)THEN
        RAH=RHOI*(WH/(2.*AH) + GAP/(2.*AN))
        GAH=1./RAH
        RAB=RHOI*(WB/(2.*AB)  + WH/(2.*AH))
        GAB=1./RAB
        CMH=CM*AH      !HILLOCK  CAPACITANCE
        GMH=GM*AH      !HILLOCK  CONDUCTANCE
        CMB=CM*AB      !CELL BOD CAPACITANCE
        GMB=GM*AB      !CELL BOD CONDUCTANCE
      ENDIF
	
      IF (NODE1.EQ.0) NODE1 = NON + 2
      J=NODE1
      NODEZ=NODE1-1
      IF(NODEZ.LE.0) NODEZ=1
C      IN(1)=NODEZ
      IN(1)=NODE1
C     MPN=10 ! MAX # OF PRINTABLE NODES STARTING WITH NODE1 ! added 1  12/26/08
      mpn=11  ! was already 10.  Added 1  12/28/08 for data.out V&I header print
      IF (NODE1 + 9 .GT. NNODES) mpn=11 -(node1+9 -nnodes)  ! added 1  12/26/08
      DO 409 I=2,MPN ! MAX NUMBER OF PRINTABLE NODES IS 10  ! added 1  12/26/08
        IN(I)=J +I-2
  409 CONTINUE

C WRITE OUT THE PARAMETERS USED FOR THIS RUN

  500 FORMAT(1H ,///,5X,'XC(CM)',1X,F10.4,13X,'DT(ms)',9X,F10.5,13X,
     1       'NON',4X,I3)
  501 FORMAT(1H ,4X,'YC    ',1X,F10.4,13X,'PD(ms)',9X,F10.4,13X,
     1'NLIN1,2',I3,','I3)
  502 FORMAT
     1(1H ,4X,'XA(CM)',1X,F10.4,13X,'TFinal(ms)',9X,F10.4,13X,'IPRNT',
     1       2X,I4)
  503 FORMAT(1H ,4X,'YA(CM)',1X,F10.4,13X,'RHOE(OHM-CM)',1X,F10.4,13X,
     1       'ITHR',3X,I3)
  504 FORMAT(1H ,4X,'DIAM(CM)',1X,F8.4,13X,'VTH(MV)',6X,F10.4,
     113X,'NODE1',2X,I3,/5X,'IWAVE',2X,I10,13X,'CM(UF/CM**2)',F10.2,
     214X,'RHOI(CM)',1X,F10.2,/,5X,'GAP(CM)',F10.6,
     313X,'GM(MHRO/CM**2)',1X,F10.3,11X,'NODEZ   ',I3,/,
     45X,'VREF',3X,F10.3,13X,'FS = 0 (POINT ELECTRODES)'//)
  104 FORMAT(1H ,4X,'DIAM(CM)',1X,F8.4,13X,'VTH(MV)',6X,F10.4,
     113X,'NODE1',2X,I3,/5X,'IWAVE',2X,I10,13X,'CM(UF/CM**2)',F10.2,
     214X,'RHOI(CM)',1X,F10.2,/,5X,'GAP(CM)',F10.6,
     313X,'GM(MHRO/CM**2)',1X,F10.3,11X,'NODEZ   ',I3,/,
     45X,'VREF',3X,F10.3,13X,'FS = 1 (UNIFORM FIELD)'//)
C     ENDIF

C       BASED ON OPTION SELECTED FOR STIMULUS CURRENT PRINT
C       IF OPTION SELECT NOT VALID, ASSUME SINGLE RECTANGULAR PULSE

      IA = 0
      IB = 0
      YMAX = -1.E38             !MAXIMUM VALUE BELOW THRESHOLD
      YMIN =  1.E38 		!MINIMUM VALUE ABOVE THRESHOLD

      GO TO (16,11,12,13,14,16,20,17,17,11,11,11,20),IWAVE

      GO TO 20
11    CONTINUE    ! Sinusoid(s)-only modes (IWAVE=2,10,11,12)
      IF (IWAVE .EQ. 2) THEN    ! Single sinusoid
        PER = 1.0/FREQ 
        WRITE(*,505) FREQ,PHASE,PER
      ELSEIF (IWAVE .EQ. 10 .OR. IWAVE .EQ. 11) THEN   ! Two sinusoids
        PER2 = 1.0/FREQ2
        PER = 1.0/FREQ
        WRITE(*,505) FREQ,PHASE,PER
        WRITE(*,555) FREQ2,PHASE2,PER2,AMP2
      ELSE       ! IWAVE=12, multiple sinusoids
        WRITE(*,*) NSINES,'  SINUSOIDS'
        WRITE(66,*) NSINES,'  SINUSOIDS'
        WRITE(*,*) '    AMP            FREQ        PHASE'
        WRITE(66,*) '    AMP            FREQ        PHASE'
        WRITE(*,10) ((SINEIN(I,J),J=1,3),I=1,NSINES)
        WRITE(66,10) ((SINEIN(I,J),J=1,3),I=1,NSINES)
        WRITE(*,*)
        WRITE(66,*)
C        WRITE(*,10) ((SINEIN2(I,J),J=1,3),I=1,NSINES)    !! Diagnostic
C        WRITE(66,10) ((SINEIN2(I,J),J=1,3),I=1,NSINES)   !! Diagnostic
        WRITE(*,*)
        WRITE(66,*)
10      FORMAT(1H ,F11.7,3X,F10.2,6X,F6.1)
      ENDIF
     
      GOTO 20

   12 WRITE(*,506) FREQ,PHASE,DCOFF !SINE OR COSINE ON PEDESTAL
      GO TO 20
   13 WRITE(*,507) TAUS
      GO TO 20
   14 WRITE(*,508) TAUS,FREQ,PHASE  !EXPONENTIAL SINUSOID
      GO TO 20
C   15 WRITE(*,509) UIO,XPD,DELAY,UIO2,XPD2 ! OLD IWAVE=6
C     GOTO 20
  16   CONTINUE
C   16 WRITE(*,1509) UIO,XPD,DELAY,UIO2,XPD2,NP ! IWAVE=1,6
1509  FORMAT(1H,4X,'INITIAL CURRENT(MA)',F12.5, 
     1 3X,'PULSE(ms)',F14.5,3X,'DELAY(ms)',F14.5,/,
     2 4X,'NTH CURRENT(MA)',F17.5,3X,'PULSE(ms)',F14.5,3X,'NPULSES 
     3  =',I3)
  505 FORMAT(2X,' SINUSOID',2X,'FREQUENCY',F10.3, 2X,'PHASE',F8.3,
     & 2X,'PERIOD',F10.3)

  555 FORMAT(' 2ND SINUSOID',2X,'FREQUENCY',F10.3,2X,'PHASE',F8.3,
     & 2X,'PERIOD',F10.3,2X,'AMP',F8.3)

  506 FORMAT(1H ,4X,'SINUSOID ON PEDESTAL',5X,'FREQUENCY',F12.5,5X,
     1 'PHASE',F12.5,'DC-OFFSET',F10.5)
  507 FORMAT(1H ,2X,'EXPONENTIAL',4X,'TAUS(ms)',F12.4)
  508 FORMAT(1H ,2X,'EXPONENTIAL SINUSOID',4X,'TAUS(ms)',F12.4,4X,
     1 'FREQUENCY',F10.4,4X,'PHASE',F10.4)
  509 FORMAT(1H ,4X,'DUAL RECTANGLE',3X,'INITIAL CURRENT(MA)',F12.5,
     1 3X,'PULSE(ms)',F14.5,3X,'DELAY(ms)',F14.5,/,
     2 22X,'SECOND  CURRENT(MA)',F12.5,3X,'PULSE(ms)',F14.5)

C       IF ITERATING WITH NEW VALUE OF IO, BEGIN HERE
C       INITIAL CONDITIONS

17    continue   ! special functions
C      WRITE(*,1509) UIO,XPD,UIO2,XPD2 ! IWAVE=8,IWAVE=9
      WRITE(*,*) ' UIO',UIO,' XPD',XPD,' UIO2',UIO2,' XPD2',XPD2       
   20 CONTINUE
        WRITE(17,*)5000,5000 
      write(30,*)5000,5000
C      IF (IWAVE.EQ.6 .AND. UIO.NE.UIOLD) UIO2=URATIO*UIO
      IF(UIO .NE. UIOLD)UIO2=URATIO*UIO
C      WRITE (*,*) 'UIO2 = ',UIO2  !! DIAGNOSTIC 10/19/03

      K=2*NON+1
202   CONTINUE
      DO 21 I = 1,K
        Y(I)=0.0
   21 CONTINUE
      JT = NLIN2-NLIN1
      IF (JT) 25,25,22
   22 JT = JT +1
      DO 23 I=1,JT
        L = K+4*I-3
        Y(L)  = .8249 !initial conditions h(0)
        Y(L+1)= .0005 !                   m(0)
        Y(L+2)= .0049 !                   p(0)
        Y(L+3)= .0268 !                   l(0)
   23 CONTINUE

C        CALCULATE EXTERNAL POTENTIALS

   25 CONTINUE
   
C      WRITE (*,*) 'EPOT CALC'   !! Diagnostic
   

      JT=2*NON+3
	  IF(TT .EQ. 2)THEN
	    IF(FS .EQ. 0)THEN
	      X3= (3-NON-1)*EL
	      X2= X3-LHN
	      X1= X3-LHN-LBH
	    ENDIF 
	  ENDIF

      DO 29 I=1,JT
        IF(FS .EQ. 0 .OR. FS .EQ. 3)THEN  ! point or wire electrode
C         FKT=I-NON-2 ! old index method
          FKT=I-NON-1
          XI=FKT*EL
            IF(TT .EQ. 2)THEN
              IF(I .EQ. 1)THEN
                WRITE(*,*)'FS',FS
                RC=SQRT((XC-X1)**2+YC**2)
                RA=SQRT((XA-X1)**2+YA**2)
                IF(FS .EQ. 0)THEN  ! point electrode
                  EPOT(1)=RHOE*UIO*(1./RA-1./RC)/PIMULT
                ELSE  ! FS=3, wire electrode
                  EPOT(1)=RHOE*UIO*(LOG(RC/RA))/(PIMULT*WIREL)
                ENDIF
                WRITE(*,*)'EPOT(1)',EPOT(1)
              ELSEIF(I .EQ. 2)THEN
                RC=SQRT((XC-X2)**2+YC**2)
                RA=SQRT((XA-X2)**2+YA**2)
                IF(FS .EQ. 0)THEN  ! point electrode
                  EPOT(2)=RHOE*UIO*(1./RA-1./RC)/PIMULT
                ELSE  ! FS=3, wire electrode
                  EPOT(2)=RHOE*UIO*(LOG(RC/RA))/(PIMULT*WIREL)
                ENDIF
              ELSEIF( I .EQ. 3)THEN
                RC=SQRT((XC-X3)**2+YC**2)
                RA=SQRT((XA-X3)**2+YA**2)
                IF(FS .EQ. 0)THEN  ! point electrode
                  EPOT(3)=RHOE*UIO*(1./RA-1./RC)/PIMULT
                ELSE  ! FS=3, wire electrode
                  EPOT(3)=RHOE*UIO*(LOG(RC/RA))/(PIMULT*WIREL)
                ENDIF
              ELSE         ! I GT 3
                RC =SQRT((XC-XI)**2+YC**2)
                RA =SQRT((XA-XI)**2+YA**2)
                IF(FS .EQ. 0)THEN  ! point electrode
                  EPOT(I)=RHOE*UIO*(1./RA-1./RC)/PIMULT
                ELSE  ! FS=3, wire electrode
                  EPOT(I)=RHOE*UIO*(LOG(RC/RA))/(PIMULT*WIREL)
                ENDIF
              ENDIF
            ELSEIF(TT .EQ. 1)THEN         
              RC =SQRT((XC-XI)**2+YC**2)
              RA =SQRT((XA-XI)**2+YA**2)
                IF(FS .EQ. 0)THEN  ! point electrode
                  EPOT(I)=RHOE*UIO*(1./RA-1./RC)/PIMULT
                ELSE  ! FS=3, wire electrode
                  EPOT(I)=RHOE*UIO*(LOG(RC/RA))/(PIMULT*WIREL)
                ENDIF
            ENDIF ! FS=0 or FS=3, point or wire electrode

        ELSE IF(FS .EQ. 1)THEN  !FS=1  uniform field, generated EPOTs
          IF(TT .EQ. 2)THEN    !CELL BODY + HILLOCK MODEL
            IF(I .EQ.1)THEN
              EPOT(1)=VREF
            ELSE IF(I .EQ.2)THEN
              EPOT(2)=EPOT(1) + UIO*(WB+WH)/2.*RHOE 
            ELSE IF(I .EQ. 3)THEN 
              EPOT(3)=EPOT(2) + UIO*(WH+GAP)/2.*RHOE
            ELSE ! ALL SUCCESSIVE NODES CELL BODY + HILLOCK
              EPOT(I)= EPOT(3) +(I-3)*UIO*RHOE*EL
            ENDIF !TT=2
          ELSE ! TT ASSUMED =1 AND FS=1
            EPOT(I) = VREF+(I-1)*UIO*RHOE*EL
          ENDIF
	    
        ELSE IF(FS .EQ. 2)THEN  !FS=2  uniform field, imported EPOTs
          EPOT(I) = EPOTIN(I) * UIO

        ENDIF
 29   CONTINUE
      WRITE(66,*)'CMB',CMB,' GMB',GMB,' CMH',CMH,' GMH',GMH
      WRITE(66,*)'CCM',CCM,' CGA',CGA,' AN',AN,' CGM',CGM
      WRITE(66,*)'AB',AB,' AH',AH,' GAH',GAH,' GAB',GAB
      WRITE(66,*)
      WRITE(66,*)'EPOT',(EPOT(MM),MM=1,NNODES)
C      WRITE(66,*)'EPT',(EPOT(MM),MM=1,NNODES) !! Diagnostic 11/9/03 [NOTE POSSIBLE 
C                                            !!  MISTAKE: Should be EPT(MM)? 12/12/06]
      DO 30 I = 1,NDIM
        DERY(I) =0.0
 30   CONTINUE
      DERY(NON+1)=1.0

C        INITIALIZE VARIABLE FOR A NEW RUN


      IPT=0
      TMAX=0.
      VMAX=0.0
C      WRITE(6,510)
  510 FORMAT(1H ,//,4X,'I(0)')
C      WRITE(6,511) UIO
       WRITE(66,*)
       WRITE(66,777)'    I=',UIO
	   write(*,777) '    I=',UIO
	   
 777   FORMAT (A6,F15.7)  !! Replaces * format for UIO, for
                          !! added precision under Absoft

C         WRITE THE CURRENT VALUES TO THE DISK FILE

  511 FORMAT(1H ,/,4X,F9.3,6X,F9.3,6X,F9.3,6X,F9.3)
      WRITE(66,512)(IN(I),I=2,MPN)
  512 FORMAT(1H ,/,5X,' TIME',10(7X,'V',I3,1x),6X,'    ')  ! 12/28/08, 9 to 10
      WRITE(66,513)(IN(I),I=2,MPN)
  513 FORMAT(1H ,9X,10(7X,'I',I3,1x),6X,'    ')  ! 12/28/08, 9 to 10


C       BEGIN INTEGRATION

      IF(NI .GT. 19)THEN
C       WRITE(5,915)
        WRITE(66,915)
        write(*,915)
 915    FORMAT(' ','MAXIMUM NUMBER OF ITERATIONS EXCEEDED, ...',
     &' PROGRAM TERMINATING') 
        GOTO 200
      ENDIF

      CALL RKGS(PRMT,Y,DERY,NDIM,IHLF,AUX,DELT2)

      NI=NI+1
      VM(NI)=VMAX
      NXGT(NI)=NNGTT         ! # nodes found exceeding threshold
      TM(NI)=TMAX
      NN(NI)=NODE
      UM(NI)=UIO
            IF(ITHR.LE.0) GO TO 200
C     SELECT NEW IO IF SEEKING THRESHOLD
      CROSS = (NXGT(NI) .GE. INNGTT)
      WRITE(*,69)'ITER#',NI,' VMAX',VMAX,'  UIO ',UIO, !!CROSS,removed 4/5/10
     &' #NODES',NXGT(NI),' TMAX',TMAX,' NODE',NODE
      WRITE(66,69)'ITER#',NI,' VMAX',VMAX,'  UIO ',UIO,!!CROSS,removed 4/5/10
     &' #NODES',NXGT(NI),' TMAX',TMAX,' NODE',NODE
69    FORMAT(' ',A5,I3,A5,F8.3,A6,F10.4,A7,I3,1X,A5,F6.2,A5,I4)!!L3 removed 4/5/10
      IF(NXGT(NI) .GE. INNGTT) THEN
        IA = IA + 1 !count the # times above threshold
C       WRITE(*,*)'IA',IA
        IF(UM(NI) .LT. YMIN)THEN !seek least uio above threshold
          YMIN = UM(NI)
          IF(ABS(YMIN/YMAX) .GT. 1..AND. ABS(YMIN/YMAX).LE. 1.016 !exit iteration process
     &    .AND. IB*IA .GT. 0)THEN
            GOTO 200
          ELSEIF(IB*IA .GT. 0)THEN
            UIO = (YMAX + YMIN)/2.
            GOTO 20
          ELSEIF(IB .EQ. 0)THEN
            UIO = UM(NI)*.5
            GOTO 20
          ENDIF	 
        ENDIF
      ENDIF
      IF(NXGT(NI) .LT. INNGTT)THEN ! NXGT(NI) .LT. INNGTT
        IB = IB + 1
C       WRITE(*,*)'IB',IB
        IF(ABS(UM(NI)) .GT. YMAX)THEN !seek greatest uio below threshold
          YMAX = UM(NI)
          IF(ABS(YMIN/YMAX) .GT. 1..AND. ABS(YMIN/YMAX) .LE. 1.016 !exit iteration process
     &    .AND. IA*IB .GT. 0)THEN
            GOTO 200
          ELSEIF(IA*IB .GT. 0)THEN
            UIO = (YMAX + YMIN)/2.
            GOTO 20
          ELSEIF(IA .EQ. 0)THEN ! IA = 0
            UIO = 2.*UM(NI)
            GOTO 20
          ENDIF !IF CRITERIA MET
        ENDIF !IF UIO .GT. YMAX
      ENDIF !IF #NODES EXCEEDING #NODE CRITERIA
      WRITE(*,*)'FELL THROUGH ITER',NI
      WRITE(*,*)'IA=',IA,' IB=',IB,' YMAX=',YMAX,' YMIN=',YMIN
200   CONTINUE

C PRINT ITERATIVE SUMMARY IF ITHR NE 0

      IF(ITHR.EQ.0) GO TO 301
C      WRITE(*,1509) UIO,XPD,DELAY,UIO2,XPD2,NP
      WRITE(66,250)
      write(*,250)
  250 FORMAT(/,4X,'ITER',4X,'VMAX',5X,'NXGT',8X,'UIO',7X,'TMAX',3X,
     1'NODE')
      DO 260 I=1,NI
        write(*,265) I,VM(I),NXGT(I),UM(I),TM(I),NN(I)
        WRITE(66,265)I,VM(I),NXGT(I),UM(I),TM(I),NN(I)
  260 CONTINUE
  265 FORMAT(' ',3X,I3,1X,F10.4,3X,I3,3X,F11.4,2X,F8.4,2X,I3)
  301 CONTINUE !END OF NP TRIALS
      IRUN=IRUN+1
      WRITE(66,305)IRUN
      WRITE(*,305)IRUN
  305 FORMAT(/' ************* END OF RUN ************* ',I2)
c     CLOSE (6)
      WRITE(17,*)5000,5000
      write(30,*)5000,5000
c      write(40,*)5000,5000
      if(ithr .eq. 1)goto 3333
c      write(31,*)5000,5000
c      write(32,*)5000,5000
c      write(33,*)5000,5000
c      write(34,*)5000,5000
c      write(35,*)5000,5000
      if(iwave .eq. 8.or. iwave .eq. 3.OR. IWAVE .EQ. 9)then
        if(tflag(1))then
          write(*,*)'threshold exceeded in node 1'
        else
          write(*,*)'threshold not exceeded in node 1'
        endif
        if(iwave .eq. 8.OR. IWAVE .EQ. 9)then
          WRITE(*,*)'NEXT UIO2 AT FREQ',FREQ,' ?'
          READ(5,*)UIO2
        endif
        if(iwave .eq. 3)then
          write(*,*)'NEXT UIO AT FREQ',FREQ,' DCOFF ? ',DCOFF
          read(*,*)uio
        endif         
        do i=1,6
          tflag(i)=.false.
        enddo
        GOTO 202
      ENDIF
      GOTO 3333
5000      WRITE(*,*)'HIT EOF ON INPUT'
        STOP
      END




