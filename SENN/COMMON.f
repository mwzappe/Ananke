      REAL LBH,LHN
      REAL*8 DELTIN,DELTOT
      INTEGER MES2(20),LENIN,LENOT
      INTEGER*2 NSINES
      LOGICAL CROSS,tflag,nfound
      INTEGER*2 NI,NNODES,NLIN1,NLIN2,NODE,NODE1,IWAVE,NP,FS,S,ITHR,
     & NTHNODE,IPRNT, pltn  !! 6/5/07 
      DIMENSION Y(5100),DERY(5100),AUX(8,5100),PRMT(5)
      COMMON/PPARAM/ VALH(5100),VALM(5100),VALP(5100),VALN(5100)
      COMMON/PPARAM/ DERH(5100),DERM(5100),DERP(5100),DERN(5100)
      DIMENSION VM(20),TM(20),NN(20),UM(20),NXGT(20),IN(10)
      COMMON F,R,T,SODO,SODI,POTO,POTI,PNAB,PKB,PPB,GL,VL,ER
      COMMON VMAX,VTH,IPT,IPRNT,NLIN1,NLIN2,NON,XPD,XPD2,UIO,UIO2
      COMMON NODE,NODE1
      COMMON IWAVE,FREQ,FREQ2,AMP2,ANGLE,ANGLE2,DCOFF,TAUS,DELAY
      COMMON PROD,PROD2,XMULT,PIMULT
      COMMON CA(4,3),CB(4,3),CGA,CGM,CCM,AREA,XA,YA,XC,YC,WIREL,EL,RHOE
      COMMON TIM(1100),EPOT(1100),EPT(1100)
      COMMON UINA(1100),UIK(1100),UIP(1100),UIL(1100)
      COMMON TMAX,TEND,ITHR,INNGTT,NNGTT,NODEZ
      COMMON/SWTCH/FS,S,pltn,TT,tflag(6),Tp,nfound
      COMMON/FLD/VREF,DIAM,NP,PL(128),PT(128)
      COMMON/CBPARM/AB,LBH,CMB,GMB !AREA,RESISTANCE,CAPACITANCE,CONDUCTANCE CELL BODY
      COMMON/HPARAM/AH,LHN,CMH,GMH !AREA,RESISTANCE,CAPACITANCE,CONDUCTANCE HILLOCK
      COMMON/WPARM/WB,DIAMB,WH,DIAMH,AN,GAP,GAH,GAB,SUMK(1100)
      COMMON/INARRAYS/EPOTIN(1100),NSINES,SINEIN(1100,3),
     &SINEIN2(1100,3),XIN(8001),YIN(8001),XCAL(2392000),
     &YCAL(2392000),YINTERP(2400001)
