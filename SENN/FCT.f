      SUBROUTINE FCT(X,Y,DERY,PRMT)

      include "COMMON.f"
C***********************************************************************
C CALCULATION OF DERIVATIVES                                           *
C PROGRAM MODIFIED FROM THE ORIGINAL PROGRAM WRITTEN BY McNEAL (1976). *
C   THE SUBROUTINE HAS BEEN FURTHER MODIFIED AS FOLLOWS:               *
C   LARGE EXP TEST REPLACED BY 0/0 TEST  L'HOPITALS RULE               *
C   CHANGED DEFN OF END POINTS  25 MAR 85                              *
C   FS SWITCH ADDED 3/17/87                                            *
C***********************************************************************
C      REAL LBH,LHN
C      INTEGER*2 NSINES
C      logical tflag,nfound
C      INTEGER*2 NNODES,NLIN1,NLIN2,NODE1,IWAVE,NP,FS,S,ITHR,IPRNT, pltn 
C      DIMENSION Y(5100),DERY(5100),PRMT(5)
      DIMENSION A(4),B(4)
C      DATA ELD/100./, PERX2/0.0002/, PI/3.141593/
      JT = 2 * NON
      NNODES=2*NON+1
C      K  = JT + 3
       K=JT+2
C      CON= 1000000.
       CON=1.
      XMULT = 1.0
      DELNDX = PRMT(3)    ! Param DELT for indexing input array

      GO TO (16,1,2,3,4,5,444,18,19,12,13,14,15),IWAVE
16    CONTINUE
C      WRITE(37,*) X,XMULT         !! Diagnostic
      GO TO 5
1     IF(X .LT. XPD)THEN
        ARG = PROD*X + ANGLE
        XMULT=SIN(ARG)              ! SINE
      ENDIF    ! SINGLE SINE
C      WRITE(37,*) X,XMULT         !! Diagnostic
      GO TO 5
12    IF(X .LT. XPD)THEN
        ARG = PROD*X + ANGLE        ! FOR FIRST SINE
        ARG2 = PROD2*X + ANGLE2     ! FOR SECOND SINE
        XMULT = SIN(ARG) + AMP2*SIN(ARG2)   ! SUM OF SINES
      ENDIF    ! SUM OF TWO SINES
C      WRITE(37,*) X,XMULT         !! Diagnostic
      GOTO 5
13    IF(X .LT. XPD)THEN
        ARG = PROD*X + ANGLE        ! FOR FIRST SINE
        ARG2 = PROD2*X + ANGLE2     ! FOR SECOND SINE
        XMULT = (1.0+AMP2*COS(ARG2))*SIN(ARG) ! AMP MODULATION
      ENDIF    ! AMPLITUDE MODULATION
C      WRITE(37,*) X,XMULT         !! Diagnostic
      GOTO 5
14    IF(X .LT. XPD)THEN
        XMULT = 0.0    ! INITIALIZING FOR SUM OF INPUT SINES
        DO I = 1,NSINES
          XMULT = XMULT+SINEIN2(I,1)*SIN(SINEIN2(I,2)*X+SINEIN2(I,3))
        END DO   ! SUM OF INPUT SINES
      ENDIF
C      WRITE(37,*) X,XMULT       !! Diagnostic
      GO TO 5
15    IF(X .LT. XPD)THEN
        I = ((X+X)/DELNDX)+1.5 ! the 0.5 for roundoff to nearest integer
        XMULT = YINTERP(I)    ! Y values only from interpolated array.
C                      Myelin calculates relative X values from spacing.
      ENDIF   ! Input wave
C      WRITE(38,*) X,XMULT,I        !! Diagnostic
      GO TO 5


C    2 XMULT=SIN(PROD*X+ANGLE)+DCOFF !SINE ON PEDESTAL !MEMBRANE VOLTAGE
C                                    GOES TO 0 AS UIO APPROACHES 0
2      XMULT=SIN(PROD*X +ANGLE)      !DECOUPLES PEDESTAL FROM UIO
      write(40,*)5000,5000
      if(x .le. xpd)then
        write(40,*)x,xmult
      else
        xmult=0.0
        write(40,*)x,xmult
      endif      
      GO TO 5
    3 XMULT=1.0/(EXP(X/TAUS))         !EXPONENTIAL
C      WRITE(37,*) X,XMULT       !! Diagnostic 
      GO TO 5
    4 XMULT=(SIN(PROD*X+ANGLE))/EXP(X/TAUS) !EXPONENTIAL SINUSOID
C      WRITE(37,*) X,XMULT         !! Diagnostic 
C       CHANGE EXTERNAL POTENTIAL BY APPROPRIATE MULTIPLICATIVE FACTOR
      GOTO 5
444   IF(X .LT. XPD)THEN
        IF(X .GT. 78.E-3)THEN
          ARG = -TAUS*(X-78.E-3)
          IF(ARG .LE. -2.)GOTO 10
          XMULT = .0618*EXP(ARG)
        ELSE
          XMULT = EXP(-2.*X)*SIN(PROD*X + ANGLE)
        ENDIF
      ENDIF	
18    IF(X .LE. XPD)THEN
        XMULTP=uio
        xmult=1.0
        WRITE(40,*)X,XMULTP
      ELSEIF(X .LT. (XPD +XPD2) ) THEN
        XMULTP=UIO2*SIN(PROD*X +ANGLE)
        xmult=sin(PROD*x +angle)
        WRITE(40,*)X,XMULTP
      else
        xmult=0.0
      ENDIF
      GOTO 5
19    IF(X .LT. XPD)THEN
        XMULTP=UIO*SIN(PROD*X +ANGLE)
        xmult=sin(PROD*x +angle)
        WRITE(40,*)X,XMULTP
      ELSEIF(X .LE. (XPD +XPD2) ) THEN
        XMULTP=uio2
        xmult=1.0
        WRITE(40,*)X,XMULTP
      else
        xmult=0.0
      ENDIF
    5 DO  6 I=1,K
        IF(IWAVE .NE. 3)THEN	
          EPT(I) = EPOT(I)*XMULT
        ELSE  ! iwave eq 3
          If(x .lt. Tp)then   !x lt time of pedestal duration
            EPT(I)=EPOT(I)*XMULT +(I-1)*DCOFF*RHOE*EL
          elseif(x .ge. Tp .and. x .le. xpd)then
            ept(i)=epot(i)*xmult
          elseif(x .gt. xpd)then
            ept(i)=0
          endif
        ENDIF
C	  write(41,*)x,ept(i)
    6 CONTINUE
C       TESTING FOR PULSE DURATION
      J=0
      IF (X .LT. XPD) GO TO 8
      IF (X .LE. (XPD+DELAY)) GOTO 10	
      IF(IWAVE .EQ. 8 .AND. X .GT. XPD .AND. X .LE. XPD+XPD2)
     & GOTO 77
      IF(IWAVE .EQ. 9 .AND. X .GT. XPD .AND. X .LE. XPD+XPD2)
     & GOTO 77
C      IF(IWAVE .EQ. 6 .AND. X .GT. XPD .AND. X .LE. XPD+XPD2) !Test for correcting IWAVE=6, 11/5/03
C     & GOTO 77
      DO I =1,NP
        IF(X .GE. PL(I). AND. X .LE. PT(I))THEN
          IK=I
          GOTO 77  ! IWAVE= 6 OR MULTIPLE PULSES
        ENDIF	
      END DO		
      GOTO 10
C       USE NEXT STIMULUS CURRENT FOR EXTERNAL POTENTIALS
 
C**********
77     CONTINUE

C      WRITE (*,*) 'EPT CALC'    !! Diagnostic
	 

      IF(TT .EQ. 2)THEN
        IF(FS .EQ. 0)THEN
          X3= (3-NON-1)*EL
          X2= X3 -LHN
          X1= X3-LHN-LBH
        ENDIF 
      ENDIF

      DO  7 I=1,JT+3  !! "+3" ADDED FOR IWAVE = 6   11/9/03
	    IF(FS .EQ. 0 .OR. FS .EQ. 3)THEN  ! point or wire electrode
C           FKT=I-NON-2 ! old index method
            FKT=I-NON-1
            XI=FKT*EL
	      IF(TT .EQ. 2)THEN
	        IF(I .EQ. 1)THEN
	          RC=SQRT((XC-X1)**2+YC**2)
	          RA=SQRT((XA-X1)**2+YA**2)
                  IF(FS .EQ. 0)THEN  ! point electrode
                    EPT(1)=XMULT*RHOE*UIO2*(1./RA-1./RC)/PIMULT
                  ELSE  ! FS=3, wire electrode
                    EPT(1)=XMULT*RHOE*UIO2*(LOG(RC/RA))/(PIMULT*WIREL)
                  ENDIF
	        ELSEIF(I .EQ. 2)THEN
	          RC=SQRT((XC-X2)**2+YC**2)
	          RA=SQRT((XA-X2)**2+YA**2)
                  IF(FS .EQ. 0)THEN  ! point electrode
                    EPT(2)=XMULT*RHOE*UIO2*(1./RA-1./RC)/PIMULT
                  ELSE  ! FS=3, wire electrode
                    EPT(2)=XMULT*RHOE*UIO2*(LOG(RC/RA))/(PIMULT*WIREL)
                  ENDIF
	        ELSEIF( I .EQ. 3)THEN
	          RC=SQRT((XC-X3)**2+YC**2)
	          RA=SQRT((XA-X3)**2+YA**2)
                  IF(FS .EQ. 0)THEN  ! point electrode
                    EPT(3)=XMULT*RHOE*UIO2*(1./RA-1./RC)/PIMULT
                  ELSE  ! FS=3, wire electrode
                    EPT(3)=XMULT*RHOE*UIO2*(LOG(RC/RA))/(PIMULT*WIREL)
                  ENDIF
	        ELSE         ! I GT 3
	          RC =SQRT((XC-XI)**2+YC**2)
	          RA =SQRT((XA-XI)**2+YA**2)
                  IF(FS .EQ. 0)THEN  ! point electrode
                    EPT(I)=XMULT*RHOE*UIO2*(1./RA-1./RC)/PIMULT
                  ELSE  ! FS=3, wire electrode
                    EPT(I)=XMULT*RHOE*UIO2*(LOG(RC/RA))/(PIMULT*WIREL)
                  ENDIF
	        ENDIF
	      ELSEIF(TT .EQ. 1)THEN         
	        RC =SQRT((XC-XI)**2+YC**2)
	        RA =SQRT((XA-XI)**2+YA**2)
                  IF(FS .EQ. 0)THEN  ! point electrode
                    EPOT(I)=XMULT*RHOE*UIO2*(1./RA-1./RC)/PIMULT  ! EPT??  10/20/03
                  ELSE  ! FS=3, wire electrode
                    EPOT(I)=XMULT*RHOE*UIO2*(LOG(RC/RA))/(PIMULT*WIREL)  ! EPT??  10/20/03
                  ENDIF
	      ENDIF ! FS=0 or FS=3, point or wire electrode

	    ELSEIF(FS .EQ. 1)THEN  ! FS=1  uniform field, generated EPOTs
	      IF(TT .EQ. 2)THEN    ! CELL BODY + HILLOCK MODEL
	        IF(I .EQ.1)THEN
	          EPT(1)=VREF*XMULT
	        ELSE IF(I .EQ.2)THEN
	          EPT(2)=EPT(1) + UIO2*(WB+WH)/2.*RHOE *XMULT
	        ELSE IF(I .EQ. 3)THEN
	          EPT(3)=EPT(2) + UIO2*(WH+GAP)/2.*RHOE*XMULT
	        ELSE ! ALL SUCCESSIVE NODES CELL BODY + HILLOCK
	          EPT(I)= EPT(3) +(I-3)*UIO2*RHOE*EL*XMULT
	        ENDIF !TT=2
	      ELSE ! TT ASSUMED =1 AND FS=1
	        EPT(I) = VREF+(I-1)*UIO2*RHOE*EL*XMULT
	      ENDIF

	    ELSE IF(FS .EQ. 2)THEN  ! FS=2  uniform field, imported EPOTs 
	      EPT(I) = EPOTIN(I)*UIO2*XMULT  ! 11/11/95  !NOT VERIFIED! 
	    
	    ENDIF
 7    CONTINUE
C*********
C      WRITE(66,*)'EPT',(EPT(MM),MM=1,NNODES)  !! Diagnostic 

      IF(IWAVE .EQ. 6.OR. IWAVE .EQ. 8 .OR. IWAVE .EQ. 9)THEN
        IF (X.GE.(XPD+XPD2+DELAY)) GO TO 10 !OUTSIDE OF TOTAL XPD
        GOTO 8
      ENDIF
      if(iwave .ne. 8.AND. IWAVE .NE. 9)then
        DO I =1,NP
          IF(X.GE.PL(I).AND. X .LE. PT(I) )THEN
            IK=I
            GOTO 8 !INSIDE XPD
          ENDIF	
        END DO		
        GOTO 10
      endif

8     IF(TT .EQ. 2)GOTO 88

      TIM(1)=CGA*(Y(2)-Y(1)+EPT(2)-EPT(1))
C      TIM(JT+1)=CGA*(Y(JT)-Y(JT+1)+EPT(JT+1)-EPT(JT+2))*CON !old ept indexing
       TIM(JT+1)=CGA*(Y(JT)-Y(JT+1)+EPT(JT+0)-EPT(JT+1))
      DO 9 I = 2,JT
c   9 TIM(I)=CGA*(Y(I-1)-2.*Y(I)+Y(I+1)+EPT(I+0)-2.*EPT(I+1)+EPT(I+2))
c    1 * CON  !old ept indexing
      TIM(I)=CGA*(Y(I-1)-2.*Y(I)+Y(I+1)+EPT(I-1)-2.*EPT(I+0)+EPT(I+1))
    9 CONTINUE 

      GO TO 20

10    IF(TT.EQ. 2)GOTO 1010
      TIM(1)=CGA*(Y(2)-Y(1))
      TIM(JT+1)=CGA*(Y(JT)-Y(JT+1))
      DO 11 I=2,JT
   11 TIM(I)=CGA*(Y(I-1)-2.*Y(I)+Y(I+1)) !(8) McNeal
   20 CONTINUE
      JT=2*NON+1
      DO 25 I=1,JT
      SUMK(I)=CGM*Y(I)
   25 DERY(I)=(TIM(I)-CGM*Y(I))/CCM
      J = 0
2020  CONTINUE
      JT=2*NON+1
      IF (NLIN1) 50,50,30 !go to 30 if at least 1 nonlinear node,else return

C IN THE EQUATIONS BELOW THERE IS THE FOLLOWING CORRESPONDENCE
C BETWEEN I AND THE SUBSCRIPTS USED IN McNEAL'S PAPER
C     I    SUBSCRIPT
C     1        H
C     2        M
C     3        P
C     4        N
C  THE EQUATIONS BELOW ARE OF THE FORM:
C     F(V)=C(V-A)/(1-EXP((A-V)/B))
C  WHICH ARE INDETERMINATE AT V=A  (0/0)
C  BY L'HOPITALS RULE F(A)=CB
C  THE DERIVATIVE DF/DV IS ALSO INDETERMINATE (0/0) AT A=V
C  BUT AGAIN BY L'HOPITALS RULE
C     DF(A)/DV=.5C
C  WE WANT TO PICK A VALUE DELV SO THAT WHEN
C   V=A+DELV THEN DF(V)=F(V)-F(A) < PER*F(A) AND WE SET F(V)=F(A)
C  .5C*DELV=DF=PER*C*B
C  DELV=2*PER*B
C
C THE PARAMETERS ARE USED IN THE ALPHA SUB X AND BETA SUB X EXPRESSIONS
C  (X=H,M,P,N) AS FOLLOWS:
C
C    ALPHA(X)=CA(X,1)(V-Y)(1-EXP((CA(X,2)-V)/CA(X,3)))**-1
C     SAME FOR BETA(X)   REPLACE CA BY BA
C      EQUATIONS FORMED IN FCT WITH A(X)=ALPHA(X)  B(X)=BETA(X)
C
C  IN ABOVE X=1,2,3,4 FOR H,M,P,N RESPECTIVELY
C      CA(1,1)=0.1  !coeff alpha h
C      CB(1,1)=4.5  !coeff beta h
C      CA(2,1)=0.36 !coeff alpha m
C      CB(2,1)=0.4  !coeff beta m
C      CA(3,1)=0.006!coeff alpha p
C      CB(3,1)=0.09 !coeff beta p
C      CA(4,1)=0.02 !coeff alpha n
C      CB(4,1)=0.05 !coeff beta n
C      CA(1,2)=-10. !term alpha h
C      CB(1,2)=45.  !term beta h
C      CA(2,2)=22.  !term alpha m
C      CB(2,2)=13.  !term beta m
C      CA(3,2)=40.  !term alpha p
C      CB(3,2)=-25. !term beta p
C      CA(4,2)=35.  !term alpha n
C      CB(4,2)=10.  !term beta n
C      CA(1,3)=6.   !denom term alpha h
C      CB(1,3)=10.  !denom term beta h
C      CA(2,3)=3.   !denom term alpha m
C      CB(2,3)=20.  !denom term beta m
C      CA(3,3)=10.  !denom term alpha p
C      CB(3,3)=20.  !denom term beta p
C      CA(4,3)=10.  !denom term alpha n
C      CB(4,3)=10.  !denom term beta n
 30   DO 49 K=NLIN1, NLIN2  ! k is the non-linear node index
        L=JT + 4 * J !jt pointing to last node i.e., 2*non+1 to begin h,m,p,n    
        DELV=PERX2*CA(1,3)
        DEL=CA(1,2)-Y(K) !h computions
C        WRITE(50,*) DEL, Y(K), K     !! DIAGNOSTIC 
        IF(ABS(DEL).GT.DELV) GO TO 31
C  VALUE NEAR 0/0
        A(1)=CA(1,1)*CA(1,3) !h computation
        GO TO 32
31      CONTINUE
        IF(-DEL/CA(1,3) .GT. 87.)THEN
          WRITE(*,*)'Clamping h computation, a(1) clamped to 1.e-36'
          WRITE(*,*)-DEL/CA(1,3)
          A(1)=1.E-36
          GOTO 32
        ENDIF
        A(1)=CA(1,1)*DEL/(1.-EXP(-DEL/CA(1,3))) !alpha h
   32   CONTINUE
        DUM =(CB(1,2) - Y(K))/CB(1,3)
        if(dum .lt. 78.)then
          B(1)=CB(1,1)/(1.+EXP(DUM))              !beta h
        else
          b(1)=0.
        endif
        DO 36 I = 2,4 ! iis m,p,n index
          DELV=PERX2*CA(I,3)                      !alpha m,p,n
          DEL=Y(K)-CA(I,2)
          IF(ABS(DEL).GT.DELV) GO TO 33
C  VALUE NEAR (0/0)
          A(I)=CA(I,1)*CA(I,3)                    !alpha m,p,n
          GO TO 34
33        CONTINUE
          if(-del/ca(i,3) .lt. 87. )then
            A(I)=CA(I,1)*DEL/(1.-EXP(-DEL/CA(I,3))) !alpha m,p,n
          else
            a(i)=1.E-36
          endif
   34     CONTINUE
          DELV=PERX2*CB(I,3)
          DEL=CB(I,2)-Y(K)                        !beta m,p,n
          IF(ABS(DEL).GT.DELV) GO TO 35
C  VALUE NEAR (0/0)
          B(I)=CB(I,1)*CB(I,3)                    !beta m,p,n
          GO TO 37
   35     B(I)=CB(I,1)*DEL/(1.-EXP(-DEL/CB(I,3))) !beta m,p,n
   36   CONTINUE
   37   CONTINUE
        DO 40 I=1,4 ! i is h,m,p,n index,Y(M)=h,m,p,n
          M=L+I
          DERY(M)=A(I)*(1.0-Y(M))-B(I)*Y(M) !dh/dt,dm/dt,dp/dt,dn/dt
   40   CONTINUE
        EFRT=(Y(K)+ER)*F/R/T
        EFRT2=EFRT*F
C       IF (EFRT.LT.174.6) GO TO 44
C       EFRT = 174.6
        IF (EFRT.LT.87.) GO TO 44
        EFRT = 87.
C       WRITE(*,*)'EFRT',EFRT,' Y(I)',Y(I)
   44   EX=EXP(EFRT)
        DEN=1.-EX

        VALH(K)=Y(L+1)
        VALM(K)=Y(L+2)
        VALP(K)=Y(L+3)
        VALN(K)=Y(L+4)
        
        DERH(K)=DERY(L+1)
        DERM(K)=DERY(L+2)
        DERP(K)=DERY(L+3)
        DERN(K)=DERY(L+4)
        
        PNA=PNAB*Y(L+1)*Y(L+2)**2
        PP=PPB*Y(L+3)**2
        PK=PKB*Y(L+4)**2
        
        UINA(k)=PNA*EFRT2*(SODO-SODI*EX)/DEN*1000. !Ina
        UIK(k)=PK*EFRT2*(POTO-POTI*EX)/DEN*1000.   !Ik
        UIP(k)=PP*EFRT2*(SODO-SODI*EX)/DEN*1000.   !Ip
        UIL(k)=GL*(Y(K)-VL)     !Il
        
        SUMK(K)=UINA(K)+UIK(K)+UIP(K)+UIL(K)
        J=J +1
C     dVn/dt per equation (9), McNeal, for non-linear nodes
        IF(TT .EQ.1)THEN
          DERY(K)=(TIM(K)-(UINA(k)+UIK(k)+UIP(k)+UIL(k))*AREA)/CCM
        ENDIF
        IF(TT .EQ. 2)THEN
          IF(K .EQ. 2)THEN !K=1 IS EXCLUDED HERE BECAUSE IT IS DEFINED AS PASSIVE
            DERY(2)=(TIM(2)- SUMK(K)*AH)/CMH
          ELSEIF(K .EQ. 3)THEN
            DERY(3)=(TIM(3)-SUMK(K)*AN)/CCM
          ELSEIF(K .GT. 3)THEN
            DERY(K)=(TIM(K)-SUMK(K)*AN)/CCM
          ENDIF
        ENDIF
   49 CONTINUE
   50 RETURN
88    CONTINUE
      TIM(1)=GAB*(Y(2)-Y(1) +EPT(2)-EPT(1))  !LIKE 8 EXCEPT CELL BODY CONDUCTANCE
      DERY(1)=(TIM(1) - GMB*Y(1))/CMB        !CELL BODY PASSIVE NODE (LINEAR) 
      TIM(2)=GAB*(Y(1)-Y(2) +EPT(1)-EPT(2)) 
     &+GAH*(Y(3)-Y(2)+EPT(3)-EPT(2))
C     DERY(2) NOT NEEDED BECAUSE HILLOCK IS ACTIVE

      DERY(2)=(TIM(2)-Y(2)*GAH)/CMH
      TIM(3)=GAH*(Y(2)-Y(3) +EPT(2)-EPT(3)) 
     &      +CGA*(Y(4)-Y(3) +EPT(4)-EPT(3))
C     DERY(3) NOT NEEDED NODE 3 IS AN ACTIVE REGULAR NODE
         DERY(3)=(TIM(3)-Y(3)*CGM)/CCM
      TIM(NNODES)=CGA*(Y(NNODES-1)-Y(NNODES)
     &             +EPT(NNODES-1)-EPT(NNODES)) !LAST NODE
      DERY(NNODES)=(TIM(NNODES)-CGM*Y(NNODES))/CCM ! IF NTH NODE IS LINEAR
      SUMK(NNODES)=CGM*Y(NNODES)
      DO 255 I=4,JT !JT=NNODES-1
      TIM(I)=CGA*(Y(I-1)-2.*Y(I)+Y(I+1)+EPT(I-1)-2.*EPT(I)+EPT(I+1))
255   DERY(I)=(TIM(I) - CGA*Y(I))/CCM   !ANY ACTIVE NODES WILL BE OVERRIDDEN LATER
      GOTO 2020
1010  CONTINUE
      TIM(1)=GAB*(Y(2)-Y(1))
      DERY(1)=(TIM(1)  - GMB*Y(1))/CMB     !CELL BODY PASSIVE NODE (LINEAR) 
      TIM(2)=GAB*(Y(1)-Y(2) ) +GAH*(Y(3)-Y(2))
C     DERY(2) NOT NEEDED BECAUSE HILLOCK IS ACTIVE
      DERY(2)=(TIM(2)-Y(2)*GAH)/CMH
      TIM(3)=GAH*(Y(2)-Y(3)) +CGA*(Y(4)-Y(3))
C     DERY(3) NOT NEEDED NODE 3 IS AN ACTIVE REGULAR NODE
      TIM(NNODES)=CGA*(Y(NNODES-1)-Y(NNODES)) !LAST NODE
      DERY(NNODES)=(TIM(NNODES) -CGM*Y(NNODES))/CCM ! IF NTH NODE IS LINEAR
      DO 256 I=4,JT !JT=NNODES-1
      TIM(I)=CGA*(Y(I-1)-2.*Y(I)+Y(I+1))
256   DERY(I)=(TIM(I) - CGA*Y(I))/CCM   !ANY ACTIVE NODES WILL BE OVERRIDDEN LATER
      GOTO 2020
      END

