      SUBROUTINE OUTP(X,Y,DERY,IHLF,NDIM,PRMT)
C****************************************************************
C     MAKES PLOT FILE FOR X VS NODE1                            *
C        STORE DATA FOR PLOTTING AND DATA PRINTOUT              *
C****************************************************************
      REAL LBH,LHN
      real VFLAG(1100),TBT(1100),TAT(1100),VBT(1100),VAT(1100),
     &TTIME(1100)
      logical tflag ,NFOUND
      INTEGER*2 NNODES,NLIN1,NLIN2,NODE1,IWAVE,FS,S,ITHR,IPRNT, pltn 
      DIMENSION Y(5100),DERY(5100),PRMT(5)
      COMMON F,R,T,SODO,SODI,POTO,POTI,PNAB,PKB,PPB,GL,VL,ER
      COMMON VMAX,VTH,IPT,IPRNT,NLIN1,NLIN2,NON,XPD,XPD2,UIO,UIO2
      COMMON NODE,NODE1
      COMMON IWAVE,FREQ,FREQ2,AMP2,ANGLE,ANGLE2,DCOFF,TAUS,DELAY
      COMMON PROD,PROD2,XMULT,PIMULT
      COMMON CA(4,3),CB(4,3),CGA,CGM,CCM,AREA,XA,YA,XC,YC,WIREL,EL,RHOE
      COMMON TIM(1100),EPOT(1100),EPT(1100)
      COMMON UINA(1100),UIK(1100),UIP(1100),UIL(1100)
      COMMON TMAX,TEND,ITHR,INNGTT,NNGTT,NODEZ
      COMMON/SWTCH/FS,S,pltn,TT,tflag(6),Tp,NFOUND
      COMMON/CBPARM/AB,LBH,CMB,GMB !AREA,RESISTANCE,CAPACITANCE,CONDUCTANCE CELL BODY
      COMMON/HPARAM/AH,LHN,CMH,GMH !AREA,RESISTANCE,CAPACITANCE,CONDUCTANCE HILLOCK
      COMMON/WPARM/WB,DIAMB,WH,DIAMH,AN,GAP,GAH,GAB,SUMK(1100)
C         UPDATE VMAX EXAMINING VOLTAGE ARRAY FROM ELECTRODE TO RIGHT
      DATA VMAXO/-999999./,UIOLD/-99999./,NODOLD/99999/,TOLD/-9999./

      JT = 2 * NON + 1
      NNODES= 2*NON +1      
      NODES=NON+1

C VMAX=0 FOR EACH NEW RUN   INITIALIZE VEL INTERPOL PARAMETERS

      IF(VMAX.NE.0.) GO TO 70
      DO 700 I=1,1100
        VFLAG(I)=0.             !voltage flag
        TBT(I)=0.
        TAT(I)=0.
        VBT(I)=0.
        VAT(I)=0.
        TTIME(I)=0.
700	  CONTINUE
      NNGTT=0

   70 CONTINUE

C FIND MAXIMUM VOLTAGE
C STARTING FROM NODE 1 , IE THE 1ST NODE, NOT NODE1 THE IST PRINT NODE
C NFOUND IS TRUE AT THE INSTANT WHEN THE IST NODE EXCEEDING THRESHOLD IS FOUND
C THIS NODE IS CONSIDERED THE EXCITATION NODE
C       DO 5 K = 1,JT
      DO 5 K = 1,NLIN2
        IF (Y(K).LE.VMAX) GO TO 5
        VMAX = Y(K) !NOW THE MAXIMUM VOLTAGE OF THE ENTIRE NODAL STRING IS FOUND
        TMAX=X
        NODE = K
    5 CONTINUE
C        FIND THE EXCITATION NODE
      IF(.not. NFOUND)THEN
        IF(VMAX .GT. VTH)THEN !MAX VOLTAGE FOUND EXCEEDS THRESHOLD
          IF(NODE.GE.NLIN1.and.NODE.LE.NLIN2)THEN
c            write(*,*)'passed node test'
            NODEXCIT=NODE
            NFOUND = .true.
c            write(*,*)'nfound',nfound
          ENDIF
        ENDIF
      ENDIF

C TRAP PARAMETERS TO INTERPOLATE FOR TIMES WHEN VOLTAGES
C  AT NODE1 THROUGH NODE6 REACH VTH
C CHANGED TO START AT NODE1-1 I.E., NODEZ

      m7=6
      if(nodez + 6 .gt. (2*non+1))then
        msd= nodez+6 - (2*non +1)
        m7= 6-msd
      endif
C      DO 75 I=2,M7     !search from y(node1) to y(node2)
C      IF(VFLAG(I).NE.0.) GO TO 75
C      VBT(I)=VAT(I)
C      VAT(I)=Y(NODE1-1+I)
C      TBT(I)=TAT(I)
C      TAT(I)=X
C      IF(Y(NODE1-1+I).LT.VTH) GO TO 75 !if < threhold,inc i & cont search
C      VFLAG(I)=1.     ! set flag(i)=1 
C      NNGTT=NNGTT+1   ! increment # nodes gt threshold
C   75 CONTINUE        ! end of node search
      DO 75 I=NODEXCIT,NLIN2 
        IF(VFLAG(I) .NE. 0.)GOTO 75
        VBT(I)=VAT(I)
        VAT(I)=Y(I)
        TBT(I)=TAT(I)
        TAT(I)=X
        IF(Y(I) .LT. VTH)GOTO 75
        VFLAG(I)=1.
        NNGTT=NNGTT+1
75    CONTINUE     
C CHECK FOR TIME TO END RUN

      IF(X.GT.TEND) GO TO 11
      IF(ITHR.EQ.0) GO TO 15

C  CHECK IF NUMBER OF NODES (NNGTT) EXCEEDING VTH IS > INNGTT

      IF (NNGTT.LT.INNGTT)GO TO 15

C CAUSE END OF RUN BY SETTING PRMT(5)=1.

   11 PRMT(5)=1.0
      NODE2=NODE1+9   ! Changed from 8 to 9, 12/23/08, for data.out
      IF(NODE2 .GT. NNODES)NODE2=NNODES 
c      WRITE(6,500) X,Y(NODEZ),(Y(K),K=NODE1,NODE2)
      WRITE(66,500) X,(Y(K),K=NODE1,NODE2)
c      WRITE(17,*)X,Y(NODEZ)

      WRITE(17,*)X,Y(NODE)   ! ADDED 8/30/97
      write(30,*)x,y(NODE1)  ! CHANGED FROM y(1) 8/30/97. Now first print node.

c        write(31,*)x,y(2)         
c        write(32,*)x,y(3)
c        write(33,*)x,y(4)
c        write(34,*)x,y(5)
c        write(35,*)x,y(6)
      if(.not. tflag(1))then
        tflag(1)=(y(1).ge. vth)
      endif
      if(.not. tflag(2))then
        tflag(2)=(y(2).ge. vth)
      endif
      if(.not. tflag(3))then
        tflag(3)=(y(3).ge. vth)
      endif
      if(.not. tflag(4))then
        tflag(4)=(y(4).ge. vth)
      endif
      if(.not. tflag(5))then
        tflag(5)=(y(5) .ge. vth)
      endif
      if(.not. tflag(6))then
        tflag(6)=(y(6) .ge. vth)
      endif
c      write(22,*)x,uina(pltn)
c      write(23,*)x,uik(pltn)
c      write(24,*)x,uip(pltn)
c      write(25,*)x,uil(pltn)
      utot= uina(pltn)+uik(pltn)+uip(pltn)+uil(pltn)
c      write(26,*)x,utot
c      WRITE(66,501) TIM(NODEZ),(TIM(K),K=NODE1,NODE2)
      WRITE(66,507)'   TIM ' ,(TIM(K),K=NODE1,NODE2)
      WRITE(66,507)'   SUMK', (SUMK(K),K=NODE1,NODE2)
      WRITE(66,507)'   DERY', (DERY(K),K=NODE1,NODE2)
C     MWZ ADDED
      WRITE(66,*) Y(NNODES + 4 * NODE)
C END MWZ

      
      IF(X .EQ. 0)THEN
        WRITE(*,*)'FIRST OUTPUTS OUTP'
      ENDIF
  507 FORMAT('  ',A7,3X,10(2X,E10.3))  ! Change from 8 to 10 (12/17/95)
      IF(VMAX .EQ. VMAXO .AND. UIO .EQ. UIOLD
     &    .AND. NODOLD .EQ. NODE .AND. TOLD .EQ. TMAX)THEN
        GOTO 839
      ELSE
        TOLD=TMAX
        NODOLD=NODE
        VMAXO=VMAX
        UIOLD=UIO
        WRITE(66,503) VMAX,TMAX,NODE,UIO
C       WRITE(5,503) VMAX,TMAX,NODE,UIO
      ENDIF
839   CONTINUE
C
C INTERPOLATE FOR TIMES WHEN VOLTAGES AT NODE1 TO NODE6 EQUAL VTH
C
      WRITE(66,84)
   84 FORMAT('0 TIMES WHEN NODES FIRST REACH VTH ')
C      NODEK=NODE1+1
      NODEK= NODEXCIT+1
c      DO 85 I=2,M7
      IF(NODEXCIT .GT. 0)THEN
        do 85 I=NODEXCIT,NLIN2
          IF(VFLAG(I).EQ.0.) GO TO 85 ! if no threshold,inc i &cont search
          TTIME(I)=TBT(I)+(VTH-VBT(I))*(TAT(I)-TBT(I))/(VAT(I)-VBT(I))
          J=I
C         WRITE(5,83)J,TTIME(I),VMAX,NODEK,Y(NODE1+1),UIO
          WRITE(66,83)J,TTIME(I),VMAX,nodexcit,Y(nodexcit),UIO
   83     FORMAT(' THRESHOLD REACHED AT NODE ' ,I3,' AT TIME = ',F10.6,
     1' FOR VMAX = ',F8.3,'   VN',I3,' = ',F8.3,' UIO = ',F10.4)
   85   CONTINUE
        V12=0.
        V23=0.
        V13=0.
c     If y(node1 +1) and y(node1 + 2) are set, v12,v23,v13 are recomputed      
        IF(VFLAG(2).NE.0.) V12=1./(TTIME(2)-TTIME(1))
        IF(VFLAG(3).NE.0.) V23=1./(TTIME(3)-TTIME(2))
        IF(VFLAG(3).NE.0.) V13=2./(TTIME(3)-TTIME(1))
        IF(VFLAG(2).NE.0) WRITE(66,86)V12,V23,V13
C       IF(VFLAG(2).NE.0) WRITE(5,86)V12,V23,V13
   86   FORMAT(' V12 = ',F10.4,' V23 = ',F10.4,' V13 = ',F12.3)
      ENDIF

      GO TO 60

   15 IPT=IPT+1
C      PRINT IF IPT POSITIVE
C      WRITE(17,*)X,Y(NODEZ)
      WRITE(17,*)X,Y(NODE)

      write(30,*)x,y(NODE1)   ! CHANGED FROM y(1) 8/30/97. Now first print node.

      utot= uina(pltn)+uik(pltn)+uip(pltn)+uil(pltn)
c      write(26,*)x,utot
      if(.not. tflag(1))then
        tflag(1)=(y(1) .ge. vth)
      endif
      if(.not. tflag(2))then
        tflag(2)=(y(2) .ge. vth)
      endif
      if(.not. tflag(3))then
        tflag(3)=(y(3) .ge. vth)
      endif
      if(.not. tflag(4))then
        tflag(4)=(y(4) .ge. vth)
      endif
      if(.not. tflag(5))then
        tflag(5)=(y(5) .ge. vth)
      endif
      if(.not. tflag(6))then
        tflag(6)=(y(6) .ge. vth)
      endif
c enable the next 3 statments to plot h,m,p,n components
c      DO I=1,4
c        WRITE(17+I,*)X,Y(2*NON+I+1) !h,m,p,n
c      ENDDO
 16   IF (IPT-1) 25,25,50
   25 NODE2 = NODE1 + 9   ! Changed from 8 to 9, 12/23/08, for data.out
      if(node2 .gt. nnodes)node2=nnodes  !clamp total #prntd to total# nodes
      WRITE(66,500) X,(Y(K),K=NODE1,NODE2)
      WRITE(66,507)'   TIM ' ,(TIM(K),K=NODE1,NODE2)
      WRITE(66,507)'   SUMK', (SUMK(K),K=NODE1,NODE2)
      WRITE(66,507)'   DERY', (DERY(K),K=NODE1,NODE2)                     
      IF(VMAX .EQ. VMAXO .AND. UIO .EQ. UIOLD
     &    .AND. NODOLD .EQ. NODE .AND. TOLD .EQ. TMAX)THEN
        GOTO 50
      ELSE
        TOLD=TMAX
        NODOLD=NODE
        VMAXO=VMAX
        UIOLD=UIO
        WRITE(66,503) VMAX,TMAX,NODE,UIO
      ENDIF
C      ENDIF
   50 IF(IPT-IPRNT) 60,51,51
   51 IPT=0
   60 RETURN

  500 FORMAT(1H /,2X,F9.4,10(2X,F10.3))  ! Change from 8 to 10 (12/17/95). Effect?
  501 FORMAT(1H ,10X,10(2X,F10.3))  ! Change from 8 to 10 (12/17/95). Effect?
  503 FORMAT(' ',4X,'MAX V(mV) =',F10.3,' AT T =',F9.4,' AT NODE ',I5,
     1' FOR UI0 = ',F9.4)
      END
