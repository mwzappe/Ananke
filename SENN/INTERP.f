      SUBROUTINE INTERP(XIN,YIN,DELTIN,LENIN,NTRP,XCAL,YCAL,
     &YINTERP,DELTOT,LENOT)
C***********************************************************************
C  TO LINEARLY INTERPOLATE INPUT ARRAY OF (X,Y) POINTS, FOR MODE       *
C  IWAVE = 13. X IS TIME, Y IS CURRENT. XIN AND YIN HAVE THE INPUT     *
C  VALUES, XCAL AND YCAL HAVE THE CALCULATED INTERPOLATED VALUES ONLY. *
C  ARRAY YINTERP HAS Y OUTPUT ONLY. IN THIS IMPLEMENTATION, SENN       *
C  MEASURES X SPACING ONCE, AND CALCULATES RELATIVE X VALUES IN        *
C  SUBROUTINE RKGS. (FOR REFERENCE, OUTPUT FILE XYINTERP HAS THE FULL  *
C  TWO-DIMENSIONAL INTERPOLATED ARRAY, WITH RELATIVE X VALUES, AS      *
C  FORMED BY THIS SUBROUTINE.)                                         *
C***********************************************************************
C      include "common.f"
      REAL*4 XIN(LENIN),YIN(LENIN),XCAL(LENOT),YCAL(LENOT),
     &YINTERP(LENOT),B,M,DELTX
      INTEGER*4 I,J,NTRP,LENIN,LENOT
      REAl*8 DELTIN,DELTOT,VALUE
      IF (LENIN .LT. 2)THEN
        WRITE (*,*) ' INPUT FILE HAS LESS THAN TWO POINTS'
        WRITE (66,*) ' INPUT FILE HAS LESS THAN TWO POINTS'
C        RETURN  ! 6/19/10
        STOP      ! 6/19/10
      END IF
      LENOT = (LENIN -1)*(NTRP + 1) + 1
      DELTOT =  DELTIN/(NTRP +1)
      XCAL(1) = XIN(1)
      YCAL(1) = YIN(1)
C      YINTERP(1) = YIN(1)   ! Don't need?
      VALUE = XIN(1)
      I = 1
      J = 0
	  K = 1   ! Index for interpolated Y array YINTERP
C      WRITE (*,*) XCAL(1),YCAL(1)   !! Diagnostic
      WRITE(2,*) XCAL(1),YCAL(1)
      YINTERP(K) = YCAL(1)   ! or YIN(1)?
      DO WHILE( I .LT. LENIN)
        J = J + 1
        K = K + 1
        VALUE =VALUE + DELTOT
        XCAL(J) = VALUE
        IF (MOD(J,(NTRP +1)).EQ.0) THEN
C          WRITE (*,*) XIN(I+1 ),YIN(I+1 )   !! Diagnostic
          WRITE(2,*) XIN(I+1 ),YIN(I+1 )
          YINTERP(K) = YIN(I+1)
          I =I +1
        ELSE
          DELTX =XIN(I+1)-XIN(I)
          IF (DELTX .LE. 0)THEN
            WRITE (*,*) 'WARNING: INPUT TIME STREAM NOT',
     &' MONOTONIC INCREASING'
            WRITE (66,*) 'WARNING: INPUT TIME STREAM NOT',
     &' MONOTONIC INCREASING'     
C            RETURN  !6/19/10
          END IF
          M = (YIN(I+1) - YIN(I))/DELTIN   !! Orig. DELTX
          B = YIN(I) -M*XIN(I)
          YCAL(J) = M*XCAL(J) + B
C          WRITE(*,*) XCAL(J),YCAL(J)   !! Diagnostic
          WRITE(2,*) XCAL(J),YCAL(J)
          YINTERP(K) = YCAL(J)
        END IF
      END DO
      WRITE (*,*) 'End of new data.    ',
     &LENOT,' (x,y) pairs created'
      WRITE (*,*) 'Output file named XYINTERP'
C      WRITE (*,*) (YCAL(K),K=1,LENOT)   !! Diagnostic
      RETURN
      END
