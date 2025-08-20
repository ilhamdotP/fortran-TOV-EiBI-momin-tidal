      PROGRAM TOVSOLVER                   
C     *********************************************************
C     
C     RMF BSP parameter set
c      for p-anisotropic slow rotating isotropic NS, 24 Nov 2015
C     R. L. Bowers and E. P. T. Liang, Astrophys. J 88, 657 (1974)
C     LBL=-2 and 2             
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
C      INTEGER I, J, IM, NS, LI, N, IL, IN    
C      DIMENSION YA(10), EK(5,10), Y(10)
       INTEGER IL
       DOUBLE PRECISION LBL,LAM
C--------------------------------------------------------
C result save in file X.dat 
C--------------------

      K=5.D6  !38.D6 !7.D6  ! m^2
      LAM=-2.08D-52   ! m^(-2)

      LBL=0.D0 ! Keep zero for isotropic pressure
      
      !OPEN (unit=5,STATUS='unknown',FILE='cekdata1.dat')
      !OPEN (unit=6,STATUS='unknown',FILE='cekdata2.dat')
      !OPEN (unit=7,STATUS='unknown',FILE='radmassMomIn.dat')
      OPEN (unit=8,STATUS='unknown',FILE='data1_momin.dat')
      !OPEN (unit=9,STATUS='unknown',FILE='data2_momin.dat')

c Set constants: 
c HC: hbar c, 
c PI, 
c GS: Newton constant in m/(m^3 MeV fm^-3)
c MSS: solar mass in m^3 MeV fm^-3
      HC  = 197.327D0
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
  
C---------------------------------------------------------------------
C  PCC is pressure in center
C  TUCT is trial initial condition for metric component nu
c Start looping for different central pressure

       DO 10 IL=600,1,-1
!       IL=600
        
       PCC=1.0D0*IL
       TUCT=0.1D-8

c       PCC=330.D0
C------------------------------------------------------------------------
c Initial steps to calculate the correct initial condition for 
C metric component nu
C   ROS is radius in km
C   GMOS is mass in solar mass at r=R=ROS
C   MNURT is metric component nu at r=R=ROS
C   Rot is the correct initial condition for metric component nu
C
      CALL TOVI(PCC,TUCT,ROS,GMOS,MNURT,K,LAM,LBL,GME,PRS)
      
      IF (PRS .GE. 1.0D9) THEN
            PRINT *,"PRS .GE. 1.0D9"
            GOTO 10
      ELSE IF (ROS .GE. 20.D0) THEN
            PRINT *,"ROS .GE. 20.D0"
            GOTO 10
      ELSE IF (GMOS .GE. 4.D0) THEN
            PRINT *,"GMOS .GE. 4.D0"
            GOTO 10
      ELSE IF (GMOS .NE. GMOS) THEN
            PRINT *,"GMOS .GE. NaN"
            GOTO 10
      ENDIF


C----------------------------------------------------------------------
C In subroutine TOVI the calculation for pressure and mass had been
C done, but initial condition for nu is NOT determined yet 
C Since Exp(nu) -> Exp(nu+constant) does not change the EoM, then
C we can obtain nu(0) -> nu(0)=nu(0)-(nu(R)-Log(1-2GM/R))

      RNS=ROS*1.D3
      DBS=1.D0-2.D0*GS*GMOS*MSS/RNS-LAM*RNS*RNS/(3.D0*L)

      KC=DLOG(DBS)-MNURT
      ROT= TUCT + KC
      
      
C----------------------------------------
C  Check metric nu
C----------------------------------------- 
C   MAMA is metric component nu in R calculate using initial correct nu
c      CALL TOVI(PCC,ROT,ROS,GMOS,MAMA)
c      CMNUR=DLOG(1.D0-2.D0*GS*GMOS*MSS/RNS)
c      WRITE(*,*)MAMA,MNURT,CMNUR,PCC,ROT,ROS,GMOS,MAMA
C--------------------------------------------------
c Moment of Inertia related properties
c----------------------------------------------------
C   ROS2 is radius in km calculate using initial correct nu
C   GMOS2 is mass in solar mass calculate using initial correct nu
C   OMEGA is rotation frequency and KAPPA is nedded to calculate moment
C   of inertia
C------------------------------------------------------
C In subroutine TOVMI, boundary conditions for both 
C omega and kappa are already satisfied using 
C multiplication by a constant zeta to both of them,
C which is similar to what we have done to nu
C------------------------------------------------------
C   MOMIN is I/MR^2 dimensionless of NS
C   MI is moment of inertia in 10^45 g cm^2 of NS
       PMIN=1.0D-9
       CALL TOVMI(PCC,ROT,PMIN,ROS2,GMOS2,MNURT,OMEGA,KAPPA,
     &            K,LAM,LBL,GME)
       RNS2=ROS2*1.D3
       MOMIN=KAPPA/(GS*GMOS*MSS*RNS2*RNS2)
       MI=MOMIN*1.98892D33*1.D10*GMOS*ROS2*ROS2/1.0D45
       C=GS*MSS*GMOS/RNS2
      
c-----------------------------------------------------------
c     CRUST Properties
c----------------------------------------------------------
C     PT pressure at core-crust transition
C     RNST is core radius
c     MOMINC I/MR^2 dimesionless of the core of NS
c     MIC is moment of Inertia    in 10^45 g cm^2 of NS core
c     MICR,MGCR,RPMIC,RPMGC are crust related properties
c--------------------------------------------------------------       
c       PT=2.863D-1 ! for eos BST
       PT=0.5658D0  ! for eos G3
       CALL TOVMI(PCC,ROT,PT,RT,GMC,MNUC,OMEGAC,KAPPAC,
     &            K,LAM,LBL,GMEC)
       RNST=RT*1.D3
       MOMINC=KAPPAC/(GS*GMC*MSS*RNST*RNST)
       MIC=MOMINC*1.98892D33*1.D10*GMC*RT*RT/1.0D45
C------------------------------------------------------------
       MICR=MI-MIC
       MGCR=GMOS2-GMC
       RPMIC=MICR/MI*1D2
       RPMGC= MGCR/GMOS2*1D2
    
      IF (GMOS .NE. GMOS) THEN
            PRINT *,"GMOS .GE. NaN"
            GOTO 10
      ELSE 
            !WRITE(7,*)IL,PCC,FED(PCC),ROS2,GME,GMOS2
            WRITE(8,*)IL,PCC,ROS2,RT,GMOS2,MOMIN,MI,C
            !WRITE(9,*)IL,GMOS2,MGCR,MICR,RPMIC,RPMGC
      ENDIF
       
 

 10   CONTINUE

       
      STOP
      END

C----------------------------------------------------------
C TOV Inertia Moment 
C--------------------------------------------
           
      SUBROUTINE TOVMI(PCC,TUCT,PMIN,ROS,GMOS,MNURT,OMEGA,
     &             KAPPA,K,LAM,LBL,GME)

      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
      INTEGER I, J, IM, NS, LI, N, IN, IL
      DIMENSION YA(10), EK(5,10), Y(10)
      DOUBLE PRECISION LBL,LAM
      HC  = 197.327D0
      PI  = 3.14159265358979D0  
      GS=1.325D-12
      MSS=1.1155D15


C     IM = NUMBER OF EQUATIONS (6)
C     IN =  NUMBER OF FIRST ORDER DIFFERENTIAL EQUATIONS (5)


      IM=6
      IN=IM-1

C---------------------------------------------------
C     Y(1)=Pressure, 
C     Y(2)=Star Mass,
C     Y(3)=Metric Nu
C     Y(4)=Omega, 
C     Y(5)=Kappa,
C----------------------------------------------------
C     Y(6)=Energy Density, 

         
      Y(1)=PCC
      Y(3)=TUCT

      Y(2)=0.1D-11    
      Y(4)=1.D-9
      Y(5)=1.D-9

     
c--------------------------------------------------------------------------
c Compute initial Energy density
c---------------------------------------------------------------------      

      
      Y(6)=FED(PCC)


     
C     PU=INTERVAL OF radius FOR PRINTING, NS=NUMBER OF STEP IN ONE PRINT
C     INTERVAL OF T,XL=MAXIMUM radius TO STOP CALCULATION

      PU=1.0D1  ! default =1.0D-1
      NS=75  ! default =32
      XL=20.0D3 ! default =30.0D3

      H=PU/NS
      XP=1.0D-3
      HH=H/(2.0D0)

C     LINE NUMBER INITIALIZATION

      LI=0
 

 28   LI=LI+1

C     XB=OLD radius, XP=NEW radius, XM=MIDPOINT radius

      DO N=1,NS
         XB=XP
         XP=XP+H
         XM=XB+HH


C    COMPUTE K1, L1, M1, N1, O1

         J=1
         DO I=1,IM
            YA(I)=Y(I)
         END DO
         XA=XB

   
         CALL FUNCX(EK,J,YA,XA,H,K,LAM,LBL)

C    COMPUTE K2, L2, M2, N2, O2

         J=2
         
         DO I=1,IN
            YA(I)=Y(I)+EK(1,I)/(2.D0)
         END DO
            P0=YA(1)
c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        YA(6)=ED
c---------------------------------------------------------------------
                     
         XA=XM

         CALL FUNCX(EK,J,YA,XA,H,K,LAM,LBL)

C    COMPUTE K3, L3, M3, N3, O3

         J=3
         DO I=1,IN
            YA(I)=Y(I)+EK(2,I)/(2.D0)
         END DO
     

          P0=YA(1)

c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        YA(6)=ED    
c---------------------------------------------------------------------- 
                               
         XA=XM


         CALL FUNCX(EK,J,YA,XA,H,K,LAM,LBL)

C    COMPUTE K4, L4, M4, N4, O4

         J=4
         DO I=1,IM
            YA(I)=Y(I)+EK(3,I)
         END DO

          P0=YA(1)

c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        YA(6)=ED

c---------------------------------------------------------------------
                 

         XA=XP

         CALL FUNCX(EK,J,YA,XA,H,K,LAM,LBL)

C    4-TH ORDER RUNGE KUTTA SCHEME

         DO I=1,IN
            Y(I)=Y(I)+(EK(1,I)+2.D0*EK(2,I)+2.D0*EK(3,I)+EK(4,I))/6.D0
         END DO
          P0=Y(1)

c--------------------------------------------------------------------------
c Call Presure vs energy and baryon density parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        Y(6)=ED
c---------------------------------------------------------------------  
                   
          
       END DO




       PS=Y(1)
       
C        IF (ABS(PMIN) .EQ. 0.5658D0) THEN
C             WRITE(6,*)Y(1),Y(4),Y(5)                             !Lagi dicek
C        ELSEIF (ABS(PMIN) .LT. 0.5658D0) THEN
C             WRITE(5,*)Y(1),Y(4),Y(5)
C        ENDIF
         
c       PT=2.863D-1 
c       PMIN=1.0D-9 
      
C        IF ((PS-0.5658D0) .LT. 2.D-4) Y(1)=Y(1)-2.D-4  ! kicking the pressure out of pt=0.5658D0
C        IF ((PS-50.D0) .LT. 2.D-4) Y(1)=Y(1)-2.D-4  ! kicking the pressure out of p=50

       
       

       IF (PS .GT. PMIN .AND. XP .LT. XL) GOTO 28

C --------------------------------------------------
C The initial condition from omet and kapat are not
C satisfied but since: 
C   omet'(r) is a linear function of kapat and 
C   kapat'(r) is a linear function of omet,
C then if zeta is a constant, changing both
C omet -> zeta*omet and kapat -> zeta*kapat
C will not change each EoM. Suppose at r=R we
C obtain omet(R)=(1-2GI/R^3)/zeta and
C kapat(R)=GI/zeta, then the boundary conditions
C will be satisfied by changing both 
C omet(R) -> zeta*omet(R) and kapat -> zeta*kapat(R)
C with zeta=1/(omet(R)+2 kapat(R)/R^3)

C----------------------------------------------------
C In subroutine FUNCMI, I had implemented zeta
C
C      ROS=(XP/1.D3)
C      GMOS=Y(2)
C      MNURT=Y(3)
      OMET=Y(4)
      KAPAT=Y(5)
      ZETA=1.D0/(OMET+2.D0*KAPAT/(XP*XP*XP))
C      OMEGA=ZETA*OMET
C      KAPPA=ZETA*KAPAT


       CALL FUNCMI(XP,Y,MPHY,OMTPHY,KPTPHY,K,LAM,LBL)   

 
      
      ROS=(XP/1.D3)
      GMOS=MPHY
      MNURT=Y(3)
      OMEGA=OMTPHY
      KAPPA=KPTPHY
      GME=Y(2)
      
      IL=PCC
      MOMINEFF=ZETA*KAPAT/(GS*Y(2)*MSS*XP*XP)
      MOMINPHY=KPTPHY/(GS*GMOS*MSS*XP*XP)
      WRITE(*,*) Y(2),GMOS,MOMINEFF,MOMINPHY

   
      RETURN
      END


C--------------------------------------------------------------
C  TOV initial condition for nu
C----------------------------------------------------------

      SUBROUTINE TOVI(PCC,TUCT,ROS,GMOS,MNURT,K,LAM,LBL,GME,PRS)
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
      INTEGER I, J, IM, NS, LI, N, IN, IL
      DIMENSION YA(10), EK(5,10), Y(10)
      DOUBLE PRECISION LBL,LAM

      HC  = 197.327D0
      PI  = 3.14159265358979D0  
      GS=1.325D-12
      MSS=1.1155D15

C     IM = NUMBER OF EQUATIONS
C------------------------------------------
C     IN =  NUMBER OF FIRST ORDER DIFFERENTIAL EQUATIONS
C     Y(1)=Pressure, 
C     Y(2)=Star Mass,
C     Y(3)=Matric Nu,  
C----------------------------------------------
C     Y(4)=Energy Density, 


      IM=4
      IN=IM-1  
      
      XP=1.0D-3

      Y(1)=PCC
      Y(2)=0.1D-11
      Y(3)=TUCT

      
c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      
    
      Y(4)=FED(PCC)
     
C     PU=INTERVAL OF radius FOR PRINTING, NS=NUMBER OF STEP IN ONE PRINT
C     INTERVAL OF T,XL=MAXIMUM radius TO STOP CALCULATION

      PU=1.0D1  ! default =1.0D-1
      NS=75  ! default =16
      XL=20.0D3 ! default =30.0D3

      H=PU/NS
      XP=1.0D-3
      HH=H/(2.0D0)

C     LINE NUMBER INITIALIZATION

      LI=0
 

 28   LI=LI+1

C     XB=OLD radius, XP=NEW radius, XM=MIDPOINT radius

      DO N=1,NS
         XB=XP
         XP=XP+H
         XM=XB+HH

C    COMPUTE K1, L1, M1

         J=1
         DO I=1,IM
            YA(I)=Y(I)
         END DO

         XA=XB
   
         CALL FUNCT(EK,J,YA,XA,H,K,LAM,LBL)
         
C    COMPUTE K2, L2, M2

         J=2
         
         DO I=1,IN
            YA(I)=Y(I)+EK(1,I)/(2.D0)
         END DO
         
        P0=YA(1)
c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        YA(4)=ED
c---------------------------------------------------------------------

        
         XA=XM

         CALL FUNCT(EK,J,YA,XA,H,K,LAM,LBL)

C    COMPUTE K3, L3, M3

         J=3
         DO I=1,IN
            YA(I)=Y(I)+EK(2,I)/(2.D0)
         END DO
     

          P0=YA(1)

c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        YA(4)=ED 
                 
c---------------------------------------------------------------------- 
           
         XA=XM

         CALL FUNCT(EK,J,YA,XA,H,K,LAM,LBL)

C    COMPUTE K4, L4, M4

         J=4
         DO I=1,IM
            YA(I)=Y(I)+EK(3,I)
         END DO

          P0=YA(1)

c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        YA(4)=ED

c---------------------------------------------------------------------
                

         XA=XP

         CALL FUNCT(EK,J,YA,XA,H,K,LAM,LBL)

C    4-TH ORDER RUNGE KUTTA SCHEME

         DO I=1,IN
            Y(I)=Y(I)+(EK(1,I)+2.D0*EK(2,I)+2.D0*EK(3,I)+EK(4,I))/6.D0
         END DO
          P0=Y(1)

c--------------------------------------------------------------------------
c Call Presure vs energy and baryon density parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        Y(4)=ED
c---------------------------------------------------------------------  
                  
          
       END DO

       !PRINT*,XP/1.D3,Y(1),Y(2)   !lagi dicek

       PS=Y(1)
       PMIN=1.0D-9
c       PMIN=2.0D-5
     
C        IF ((PS-0.5658D0) .LT. 2.D-4) Y(1)=Y(1)-2.D-4  ! kicking the pressure out of pt=0.5658D0
C        IF ((PS-50.D0) .LT. 2.D-4) Y(1)=Y(1)-2.D-4  ! kicking the pressure out of p=50
       IF (PS .GT. PMIN .AND. XP .LT. XL) GOTO 28
      PRS=Y(1)
      
      CALL FUNCM(XP,Y,MPHY,K,LAM,LBL)

      ROS=(XP/1.D3)
      GMOS=MPHY
      MNURT=Y(3)
      GME=Y(2)
      
      IL=PCC
      WRITE(*,*) IL
      WRITE(*,*) XP/1.D3
      WRITE(*,*) Y(2),GMOS
           
      
      RETURN
      END


       SUBROUTINE FUNCX(EK,J,YA,XA,H,K,LAM,LBL)
C     *********************************************************
C     DEFINES TOV EQUATIONS AIP BL 
C     P-unisotropic used to calculate Moment of inertia
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z) 
      INTEGER J   
      DIMENSION EK(5,10), YA(10)
      DOUBLE PRECISION LBL,LAM,OMGT,KAPT
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
C ----------------------------------
c    model parameter !! Mmax=2.08 Ms
  
      !LBL=1.82D0   
      !LBL=0.D0 
      !K=7.D6
      L=K*LAM+1.D0
C--------------------------------------       
      PRESS=YA(1)
      MASST=YA(2)
      MNU=YA(3)
      OMGT=YA(4)
      KAPT=YA(5)
      EDEN=YA(6) 
      EXPMNU=EXP(MNU/2.D0)
      EXPMMNU=EXP(-MNU/2.D0)
     
      A=DSQRT(L+8.D0*PI*GS*K*EDEN)
      B=DSQRT(L-8.D0*PI*GS*K*PRESS)
      
      EDENEFF=(A*A-3.D0*B*B+2.D0*A*B*B*B)
     &        /(16.D0*PI*GS*K*A*B*B*B)
      PRESSEFF=(A*A+1.D0*B*B-2.D0*A*B*B*B)
     &        /(16.D0*PI*GS*K*A*B*B*B)
          
      CS2=1.D0/DEDP(PRESS)
          
C Persamaan pressure
      XXX=(4.D0/(A*A-B*B)+3.D0/(B*B)+1.D0/(A*A*CS2))
     &    *(1.D0-2.D0*GS*MSS*MASST/XA-LAM*XA*XA/(3.D0*L))
      
      EK(J,1)=-H*(1.D0/(4.D0*PI*GS*K)*(XA/(2.D0*K)
     &        *(1.D0/(A*B)+A/(B*B*B)-2.D0)
     &        + 2.D0*GS*MSS*MASST/(XA*XA)
     &        + LAM*XA/(3.D0*L))
     &        /XXX)
     
C Persamaan massa
      EK(J,2)=(H/MSS)*XA*XA/(4.D0*GS*K)
     &        *(2.D0/L-3.D0/(A*B)+A/(B*B*B))
      
C Persamaan nu     
      EK(J,3)=H*(-2.D0*EK(J,1)/H*(2.D0*PI*GS*K
     &      *(3.D0/(B*B)+1.D0/(A*A*CS2))
     &      +1.D0/(EDEN+PRESS)))

C Persamaan omega tilde = omet
      EK(J,4)=H*6.D0*EXPMNU/SQRT(1.D0-2.D0*GS*MASST*MSS/XA
     &        -LAM*XA*XA/(3.D0*L))
     &        *KAPT/(XA*XA*XA*XA)
     
     
C Persamaan kappa tilde = kapat
      EK(J,5)=H*(8.D0*PI*GS/3.D0)*(XA*XA*XA*XA)*EXPMMNU
     &        *EDENEFF*(1.D0+PRESSEFF/EDENEFF)
     &        *OMGT/SQRT(1.D0-2.D0*GS*MASST*MSS/XA
     &        -LAM*XA*XA/(3.D0*L))
     
     
      !PRINT*, "tes"
      !PRINT*,OMGT,KAPT
      !STOP

 
      
      RETURN
      END


      SUBROUTINE FUNCT(EK,J,YA,XA,H,K,LAM,LBL)
C     *********************************************************
C     DEFINES TOV EQUATIONS AIP BL 
C     P-anisotropic used to calculate initial condition for nu
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z) 
      INTEGER J   
      DIMENSION EK(5,10), YA(10)
      HC  = 197.327D0
      PI  = 3.14159265358979D0     
      GS=1.325D-12
      MSS=1.1155D15
C ----------------------------------
c    model parameter !! Mmax=2.08 Ms
  
      !LBL=1.82D0   
      !LBL=0.D0 
      !K=7.D6
C-------------------------------------- 
      
      L=LAM*K+1.D0      
      PRESS=YA(1)
      MASST=YA(2)
      EDEN=YA(4)
     
      A=SQRT(L+8.D0*PI*GS*K*EDEN)
      B=SQRT(L-8.D0*PI*GS*K*PRESS)
          
      CS2=1.D0/DEDP(PRESS)
          
C Persamaan pressure
      XXX=(4.D0/(A*A-B*B)+3.D0/(B*B)+1.D0/(A*A*CS2))
     &    *(1.D0-2.D0*GS*MSS*MASST/XA-LAM*XA*XA/(3.D0*L))
      
      EK(J,1)=-H*(1.D0/(4.D0*PI*GS*K)*(XA/(2.D0*K)
     &        *(1.D0/(A*B)+A/(B*B*B)-2.D0)
     &        + 2.D0*GS*MSS*MASST/(XA*XA)
     &        + LAM*XA/(3.D0*L))
     &        /XXX)
     
C Persamaan massa
      EK(J,2)=(H/MSS)*XA*XA/(4.D0*GS*K)
     &        *(2.D0/L-3.D0/(A*B)+A/(B*B*B))
      
C Persamaan nu     
      EK(J,3)=H*(-2.D0*EK(J,1)/H*(2.D0*PI*GS*K
     &      *(3.D0/(B*B)+1.D0/(A*A*CS2))
     &      +1.D0/(EDEN+PRESS)))
      
C Semua persamaan di EK harus dikali H

      !PRINT*,XA/1.D3,PRESS,MASST   !lagi dicek
      !IF (XA .GT. 1.D0) STOP

   
      RETURN
      END

 
C---------------------------------------------------------------------
C  Energy density as a function of pressure
C------------------------------------------------------------------


CC ==============================================================================================

CC        EoS WH 

C       FUNCTION FED(P0)
C       IMPLICIT DOUBLE PRECISION (A-H,K-Z)

C         IF ( P0 .GT. 50.D0 ) THEN


C         FED=        207.22479761675402 + 
C      -  6.279906372059321*P0 - 
C      -  0.021350633175731673*P0**2 + 
C      -  0.00003253273492751005*P0**3 + 
C      -  1.4247301478353118e-7*P0**4 - 
C      -  5.796734043546596e-10*P0**5 + 
C      -  5.80259303998327e-13*P0**6




  

C         ELSE IF ( P0 .GT. 0.5658D0  .AND. P0 .LE. 50.D0 ) THEN


C         FED=        60.35474983120075 + 
C      -  43.358339553599556*P0 - 
C      -  4.694837086435233*P0**2 + 
C      -  0.2844683654092531*P0**3 - 
C      -  0.008749589874350797*P0**4 + 
C      -  0.00013183647143746855*P0**5 - 
C      -  7.738363930359054e-7*P0**6

     
     
C        ELSE IF (P0 .GT. 5.0D-4 .AND. P0 .LE. 0.5658D0 ) THEN


C        FED=        0.46844150803264295 + 
C      -  690.9248196360139*P0 - 
C      -  4130.326006654172*P0**2 + 
C      -  14321.977445145929*P0**3 - 
C      -  26369.24811082957*P0**4 + 
C      -  24451.355662621238*P0**5 - 
C      -  8987.957232232802*P0**6


C         ELSE
        
C        FED=        0.00020663104786834376 + 
C      -  985.755096204799*P0 - 
C      -  6.452649548410587e6*P0**2 + 
C      -  4.045493650683301e10*P0**3 - 
C      -  1.2422017897553961e14*P0**4 + 
C      -  1.7654930953548538e17*P0**5 - 
C      -  9.152941701213566e19*P0**6

      
       
C         END IF

   
C       RETURN
C       END


CC ==============================================================================================

CC        EoS WHSS 

      FUNCTION FED(P0)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)

        IF ( P0 .GT. 50.D0 ) THEN


        FED=        56.336121316129855 + 
     -  11.816705939891484*P0 - 
     -  0.08998007668063691*P0**2 + 
     -  0.0004060689424004221*P0**3 - 
     -  9.838173200421972e-7*P0**4 + 
     -  1.219639435236914e-9*P0**5 - 
     -  6.069131445225517e-13*P0**6




  

        ELSE IF ( P0 .GT. 0.5658D0  .AND. P0 .LE. 50.D0 ) THEN


        FED=        39.0297275662326 + 56.62607981986861*P0 - 
     -  6.9282946501671265*P0**2 + 
     -  0.44222014371897306*P0**3 - 
     -  0.014109491546204781*P0**4 + 
     -  0.00021868926715498178*P0**5 - 
     -  1.3126683542897572e-6*P0**6

     
     
       ELSE IF (P0 .GT. 5.0D-4 .AND. P0 .LE. 0.5658D0 ) THEN


       FED=        0.293210281408273 + 
     -  742.9038848085958*P0 - 
     -  5635.80467355953*P0**2 + 
     -  27146.16793869607*P0**3 - 
     -  71721.53203480988*P0**4 + 
     -  95705.23802777313*P0**5 - 
     -  50358.11224200026*P0**6


        ELSE
        
       FED=        0.00020663104786834376 + 
     -  985.755096204799*P0 - 
     -  6.452649548410587e6*P0**2 + 
     -  4.045493650683301e10*P0**3 - 
     -  1.2422017897553961e14*P0**4 + 
     -  1.7654930953548538e17*P0**5 - 
     -  9.152941701213566e19*P0**6

      
       
        END IF

   
      RETURN
      END


CC ==============================================================================================
      
CC        EoS WoutH 

C       FUNCTION FED(P0)
C       IMPLICIT DOUBLE PRECISION (A-H,K-Z)

C         IF ( P0 .GT. 50.D0 ) THEN


C         FED=        234.1942349585189 + 
C      -  4.244398207462238*P0 - 
C      -  0.013395621016951648*P0**2 + 
C      -  0.000045526410778733156*P0**3 - 
C      -  8.910615387698381e-8*P0**4 + 
C      -  9.208285152708568e-11*P0**5 - 
C      -  3.8734411956242993e-14*P0**6




  

C         ELSE IF ( P0 .GE. 0.5658D0  .AND. P0 .LE. 50.D0 ) THEN


C         FED=        61.15142108633361 + 
C      -  42.41019902820945*P0 - 
C      -  4.486893551483503*P0**2 + 
C      -  0.2701676379200499*P0**3 - 
C      -  0.008407677119287372*P0**4 + 
C      -  0.0001284767514904723*P0**5 - 
C      -  7.626291720876074e-7*P0**6

     
     
C        ELSE IF (P0 .GT. 5.0D-4 .AND. P0 .LT. 0.5658D0 ) THEN


C        FED=        0.4685189238250536 + 
C      -  690.9060074153616*P0 - 
C      -  4129.960666603701*P0**2 + 
C      -  14319.73083001859*P0**3 - 
C      -  26363.06652679627*P0**4 + 
C      -  24443.207464873816*P0**5 - 
C      -  8983.86940730949*P0**6


C         ELSE
        
C        FED=        0.00020663104786834376 + 
C      -  985.755096204799*P0 - 
C      -  6.452649548410587e6*P0**2 + 
C      -  4.045493650683301e10*P0**3 - 
C      -  1.2422017897553961e14*P0**4 + 
C      -  1.7654930953548538e17*P0**5 - 
C      -  9.152941701213566e19*P0**6

      
       
C         END IF

   
C       RETURN
C       END


CC==============================================================================================
      
CC       EoS WoutHSS 

C       FUNCTION FED(P0)
C       IMPLICIT DOUBLE PRECISION (A-H,K-Z)

C         IF ( P0 .GT. 50.D0 ) THEN


C         FED=        293.11268673102506 + 
C      -  3.10114370497162*P0 - 
C      -  0.010901925280260074*P0**2 + 
C      -  0.00004154068779880557*P0**3 - 
C      -  8.163579780831742e-8*P0**4 + 
C      -  7.974755635410319e-11*P0**5 - 
C      -  3.063161763646646e-14*P0**6




  

C         ELSE IF ( P0 .GT. 0.5658D0  .AND. P0 .LE. 50.D0 ) THEN


C         FED=        61.15142108633361 + 
C      -  42.41019902820945*P0 - 
C      -  4.486893551483503*P0**2 + 
C      -  0.2701676379200499*P0**3 - 
C      -  0.008407677119287372*P0**4 + 
C      -  0.0001284767514904723*P0**5 - 
C      -  7.626291720876074e-7*P0**6

     
     
C        ELSE IF (P0 .GT. 5.0D-4 .AND. P0 .LE. 0.5658D0 ) THEN


C        FED=        0.4685189238250536 + 
C      -  690.9060074153616*P0 - 
C      -  4129.960666603701*P0**2 + 
C      -  14319.73083001859*P0**3 - 
C      -  26363.06652679627*P0**4 + 
C      -  24443.207464873816*P0**5 - 
C      -  8983.86940730949*P0**6


C         ELSE
        
C        FED=        0.00020663104786834376 + 
C      -  985.755096204799*P0 - 
C      -  6.452649548410587e6*P0**2 + 
C      -  4.045493650683301e10*P0**3 - 
C      -  1.2422017897553961e14*P0**4 + 
C      -  1.7654930953548538e17*P0**5 - 
C      -  9.152941701213566e19*P0**6

      
       
C         END IF

   
C       RETURN
C       END

CC==============================================================================================
      
CC       EoS naive fitting

C       FUNCTION FED(P0)
C       IMPLICIT DOUBLE PRECISION (A-H,K-Z)
            ! G3 WH
C         FED=        35.38875313226208 + 17.24847480749493*P0 - 
C      -  0.2448168236745378*P0**2 + 
C      -  0.0020467888116389543*P0**3 - 
C      -  8.782852352403798e-6*P0**4 + 
C      -  1.8483442745836703e-8*P0**5 - 
C      -  1.5101080861937848e-11*P0**6
            ! G3 WHSS
C         FED=        38.50467799910013 + 15.685128505945936*P0 - 
C      -  0.18084960995412983*P0**2 + 
C      -  0.001196104238390078*P0**3 - 
C      -  4.1500340315907694e-6*P0**4 + 
C      -  7.155178985377673e-9*P0**5 - 
C      -  4.823487894465202e-12*P0**6
            ! G3 WoutH
C         FED=        47.90062932721776 + 12.225768076115184*P0 - 
C      -  0.11457226593449464*P0**2 + 
C      -  0.000600013940954589*P0**3 - 
C      -  1.569481469885387e-6*P0**4 + 
C      -  1.9891120510167144e-9*P0**5 - 
C      -  9.729037482565278e-13*P0**6
            ! G3 WoutHSS
C         FED=        51.03925729566848 + 11.483991820732193*P0 - 
C      -  0.09784975654592297*P0**2 + 
C      -  0.00043469461945059353*P0**3 - 
C      -  9.522102958267907e-7*P0**4 + 
C      -  1.0081617266058726e-9*P0**5 - 
C      -  4.119793693767378e-13*P0**6
C       RETURN
C       END

CC==============================================================================================
      
CC       EoS linear

C       FUNCTION FED(P0)
C       IMPLICIT DOUBLE PRECISION (A-H,K-Z)
C         FED=        240.D0 + 3.D0*P0
C       RETURN
C       END

C ==============================================================================================

CC     EoS BST

C       FUNCTION FED(P0)
C       IMPLICIT DOUBLE PRECISION (A-H,K-Z)

C         IF ( P0 .GT. 50.D0 ) THEN


C         FED=752.1779439721782D0 - 40.108620203769D0*P0 + 
C      -  1.4752000864096837D0*P0**2 - 
C      -  0.027036299898738337D0*P0**3 + 
C      -  0.0003125837415657684D0*P0**4 - 
C      -  2.441653734444954D-6*P0**5 + 
C      -  1.3308776414303D-8*P0**6 - 
C      -  5.1246996017212506D-11*P0**7 + 
C      -  1.3886907993633848D-13*P0**8 - 
C      -  2.590695651603915D-16*P0**9 + 
C      -  3.166994006343349D-19*P0**10 - 
C      -  2.283247911494129D-22*P0**11 + 
C      -  7.357464623695704D-26*P0**12



  

C         ELSE IF ( P0 .GT.  2.863D-1 .AND. P0 .LE. 50.D0 ) THEN


C         FED=65.94376092512762D0 + 57.952630512664356D0*P0 - 
C      -  15.437797576854498D0*P0**2 + 
C      -  2.7511762122977346D0*P0**3 - 
C      -  0.299643408613215D0*P0**4 + 
C      -  0.02089083577091044D0*P0**5 - 
C      -  0.0009679766327308126D0*P0**6 + 
C      -  0.000030437894581131215D0*P0**7 - 
C      -  6.521271499180282D-7*P0**8 + 
C      -  9.372308503275227D-9*P0**9 - 
C      -  8.643173476154485D-11*P0**10 + 
C      -  4.621195967220976D-13*P0**11 - 
C      -  1.0889239514378402D-15*P0**12
     
C         ELSE IF (P0 .GT. 4.99313436D-4 .AND. P0 .LE. 2.863D-1) THEN


C         FED=0.05015663787134234D0 + 836.2363942486941D0*P0 - 
C      -  9315.146969977652D0*P0**2 + 79689.1930322726D0*P0**3 - 
C      -  412197.6475732246D0*P0**4 + 1.116366190255507D6*P0**5 - 
C      -  1.1988188657021397D6*P0**6


C         ELSE
        
C         FED=0.00020663104786863406D0 + 985.7550962048012D0*P0 - 
C      -  6.452649548410687D6*P0**2 + 4.045493650683396D10*P0**3 - 
C      -  1.2422017897554384D14*P0**4 + 1.765493095354931D17*P0**5 - 
C      -  9.15294170121406D19*P0**6

        
C         END IF

   
C       RETURN
C       END

C ==============================================================================================

      SUBROUTINE FUNCM(XA,YA,MPHY,K,LAM,LBL)
C     *********************************************************
C     Definition for physical mass MPHY
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)   
      DIMENSION YA(10)
      DOUBLE PRECISION LBL,LAM
      HC  = 197.327D0
      PI  = 3.14159265358979D0     
      GS=1.325D-12
      MSS=1.1155D15
C ----------------------------------
      !LBL=1.82D0   
      !LBL=0.D0 
      !K=7.D6
      L=K*LAM+1.D0
      
      PRESS=YA(1)
      EDEN=YA(4)
      MASST=YA(2)
      
      SGM=-(LBL*GS*EDEN*EDEN*XA*XA/3.D0)
     &        *(1.D0+3.D0*PRESS/EDEN)*(1.D0+PRESS/EDEN)
     &        /(1.D0-2.D0*GS*MASST*MSS/XA-LAM*XA*XA/(3.D0*L))

      A=DSQRT(L+8.D0*PI*GS*K*EDEN)
      B=DSQRT(L-8.D0*PI*GS*K*PRESS)
      C=DSQRT(L-8.D0*PI*GS*K*(PRESS-SGM))
      
C       MPHY=XA/(2.D0*GS*MSS*DSQRT(A*B))
C      &     *(1.D0+A*C*C/(B*B)*(2.D0*GS*MASST*MSS/XA-1.D0)
C      &     +(LAM*XA*XA/(3.D0*L))*(A*C*C/B-1.D0/(A*B)))
      MPHY=MASST

   
      RETURN
      END


      SUBROUTINE FUNCMI(XA,YA,MPHY,OMTPHY,KPTPHY,K,LAM,LBL)
C     *********************************************************
C     Definition for physical mass MPHY and moment of inertia MIPHY
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)   
      DIMENSION YA(10)
      DOUBLE PRECISION LBL,LAM
      EXTERNAL FED,DEDP
      
      HC  = 197.327D0
      PI  = 3.14159265358979D0     
      GS=1.325D-12
      MSS=1.1155D15
C ----------------------------------
      !LBL=1.82D0   
      !LBL=0.D0 
      !K=7.D6
      L=K*LAM+1.D0

      PRESS=YA(1)
      MASST=YA(2)
      MNU=YA(3)
      OMET=YA(4)
      KAPAT=YA(5)
      EDEN=YA(6)
      
      

      A=DSQRT(L+8.D0*PI*GS*K*EDEN)
      B=DSQRT(L-8.D0*PI*GS*K*PRESS)
      C=B
      
C      MPHY=XA/(2.D0*GS*MSS*DSQRT(A*B))
C     &     *(1.D0+A*C*C/(B*B)*(2.D0*GS*MASST*MSS/XA-1.D0)
C     &     +(LAM*XA*XA/(3.D0*L))*(A*C*C/B-1.D0/(A*B)))
      MPHY=MASST

C ----------------------------------
C To convert omega & kappa to physical ones,
C I need to fix them by zeta first
      ZETA=1.D0/(OMET+2.D0*KAPAT/(XA*XA*XA))
      
      
      MIA=KAPAT/GS*ZETA
      MIPHY=MIA

C       SOS=1.D0/DEDP(PRESS)
C       XXX=(4.D0/(A*A-B*B)+3.D0/(B*B)+1.D0/(A*A*CS2))
C      &    *(1.D0-2.D0*GS*MSS*MASST/XA-LAM*XA*XA/(3.D0*L))
C       PERSP=-(1.D0/(4.D0*PI*GS*K)*(XA/(2.D0*K)
C      &        *(1.D0/(A*B)+A/(B*B*B)-2.D0)
C      &        + 2.D0*GS*MSS*MASST/(XA*XA)
C      &        + LAM*XA/(3.D0*L))
C      &        /XXX)      
C       MIPHY=MIA/(B*B*B*B)
C      &     +4.D0*PI*K*XA*XA*XA*XA/3.D0
C      &     *OMET*ZETA
C      &     *(1.D0/(SOS*A*A*B*B*B*B)*PERSP
C      &     +1.D0/(B*B*B*B*B*B)*PERSP)
           
      OMTPHY=A*A/(B*B)*ZETA*OMET 
      KPTPHY=GS*MIPHY
      
      
C       !lagi dicek 
C       !MINPHY bisa negatif di core-crust transition saat K terlalu besar 
C
C       PRINT*,"Cek1",MIPHY,MIA/(B*B*B*B),
C      &     (4.D0*PI*K*XA*XA*XA*XA/3.D0
C      &     *OMET*ZETA
C      &     *(1.D0/(SOS*A*A*B*B*B*B)*PERSP
C      &     +1.D0/(B*B*B*B*B*B)*PERSP))
C       PRINT*,"Cek2",4.D0*PI*K*XA*XA*XA*XA/3.D0
C      &     *OMET*ZETA,
C      &     1.D0/(SOS*A*A*B*B*B*B)*PERSP,
C      &     1.D0/(B*B*B*B*B*B)*PERSP
C       PRINT*,"Cek3",1.D0/SOS,PERSP
      
      
C      !Cek hasilnya

   
      RETURN
      END
      
C-----------------------------------------------------------------------
C     DE/DP AS A FUNCTION OF (P0) for NS (ANTO'S VERSION)
C-----------------------------------------------------------------------
      FUNCTION DEDP(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL FED
      h  = 1.D-5
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
   
   
      DEDP = (FED(x4)-8.D0*FED(x2)+8.D0*FED(x1)-FED(x3))/(12.D0*h)
c      WRITE(*,*)xa,DEDP
      RETURN
      END
 
c---------------------------------------------------------------------





C end of the code
