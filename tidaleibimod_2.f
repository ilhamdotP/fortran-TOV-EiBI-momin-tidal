       PROGRAM TOVSOLVER
C     *********************************************************
C     
C      BS 7 Maret 2018, Tidal deformation
C                  
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-M,O-Z) 
      INTEGER I, J, IM, NS, LI, N, IL, IN    
      DIMENSION YA(10), EK(4,10), Y(10)
      double complex HPL1
      

      K=5.D6 !7.D6  ! m^2
      LAM=-2.08D-52    ! m^(-2)
      

      HC  = 197.327D0
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
      
      

      OPEN (unit=5,STATUS='unknown',FILE='love2.dat')
      !OPEN (unit=1,STATUS='unknown',FILE='RadmassNS.dat')
      !OPEN (unit=2,STATUS='unknown',FILE='TIDAL_Y2_CS.dat')
      !OPEN (unit=3,STATUS='unknown',FILE='TIDAL_Y3_CS.dat')
      !OPEN (unit=4,STATUS='unknown',FILE='TIDAL_Y4_CS.dat')

C     IM = NUMBER OF EQUATIONS Y(1)=Pressure, Y(2)=NS Mass and Y(4)=E density
C     Y(3)=tidal function Y=rH'/H

      IM=4
      IN=IM-1

      
      
      DO 10 IL=1,600,1
!       IL=51   ! for testing
      PCC=1.0D0*IL

      DO 20 IY=2,4,1
!      IY=2
      YL=1.D0*IY

      !IL=3
      !PCC=1.D0*IL

      
      Y(1)=PCC
      Y(2)=1.D-12
      Y(3)=YL

      P0=Y(1)
c--------------------------------------------------------------------------
c Call Presure vs energy density parametrization (EOS)
c---------------------------------------------------------------------      


      ED=FED(P0)
c---------------------------------------------------------------------
      EDC=ED
      Y(4)=EDC
    

     
C     PU=INTERVAL OF radius FOR PRINTING, NS=NUMBER OF STEP IN ONE PRINT
C     INTERVAL OF T,XL=MAXIMUM radius TO STOP CALCULATION

      PU=1.0D1  ! default =1.0D-1       ! mine=1.0D1
      NS=75 !75  ! default =32              ! mine=75
      XL=20.0D3 ! default =30.0D3       ! mine=30.D3

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
   
         CALL FUNCT(EK,J,YA,XA,H,YL,K,LAM)

C    COMPUTE K2, L2,M2

         J=2
         
         DO I=1,IN
            YA(I)=Y(I)+EK(1,I)/(2.D0)
         END DO
            P0=YA(1)
c--------------------------------------------------------------------------
c Call Presure vs energy density parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
c---------------------------------------------------------------------

         YA(4)=ED

         XA=XM

         CALL FUNCT(EK,J,YA,XA,H,YL,K,LAM)

C    COMPUTE K3, L3,M3

         J=3
         DO I=1,IN
            YA(I)=Y(I)+EK(2,I)/(2.D0)
         END DO
     

          P0=YA(1)

c--------------------------------------------------------------------------
c Call Presure vs energy density parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
                 
c---------------------------------------------------------------------- 
            YA(4)=ED          

         XA=XM


         CALL FUNCT(EK,J,YA,XA,H,YL,K,LAM)

C    COMPUTE K4, L4,M4

         J=4
         DO I=1,IM
            YA(I)=Y(I)+EK(3,I)
         END DO

          P0=YA(1)

c--------------------------------------------------------------------------
c Call Presure vs energy density parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
c---------------------------------------------------------------------
            YA(4)=ED
     

         XA=XP

         CALL FUNCT(EK,J,YA,XA,H,YL,K,LAM)

C    4-TH ORDER RUNGE KUTTA SCHEME

         DO I=1,IN
            Y(I)=Y(I)+(EK(1,I)+2.D0*EK(2,I)+2.D0*EK(3,I)+EK(4,I))/6.D0
         END DO
          P0=Y(1)

c--------------------------------------------------------------------------
c Call Presure vs energy density parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)

c---------------------------------------------------------------------  
          Y(4)=ED
          
       END DO

C      WRITE(*,*)IL,LI
C      WRITE(5,*)IL,IY,(XP/1.D3),Y(1),Y(4),Y(2),Y(3) 
C      WRITE(*,*)(XP/1.D3),Y(1)  ! for testing
          

       PS=Y(1)
       PMIN=1.0D-8
c       PMIN=2.0D-5
     
      
c      IF ((PS-0.5658D0) .LT. 0.9D-4) Y(1)=Y(1)+0.9D-4  ! kicking the pressure out of pt=0.5658D0
c      IF ((PS-50.D0) .LT. 2.D-4) Y(1)=Y(1)+2.D-4  ! kicking the pressure out of p=50
      IF (PS .GT. PMIN .AND. XP .LT. XL) GOTO 28
      
      IF ( YL .EQ. 2.D0 ) THEN 
        WRITE(*,*)IL
        WRITE(*,*)(XP/1.D3),Y(2),GS*Y(2)*MSS/XP
      ENDIF
      
      IF (Y(2).LE.4.D0 .AND. XP .LT. XL) THEN

      CALL FUNCM(XP,Y,MPHY,K,LAM)
      
      !IF ( YL .EQ. 2.D0 ) THEN 
        !WRITE(1,*)IL,PCC,EDC,(XP/1.D3),Y(2),MPHY
      !ENDIF

    
      IF ( YL .EQ. 2.D0 ) THEN 

       MASST=MPHY
       YR=Y(3)

       CALL FUNCK2(XP,MASST,YR,K,LAM,C,K2)
       !WRITE(2,*)IL,(XP/1.D3), MASST,YR,C,K2
       !WRITE(*,*)YL,YR,K2
       
       LOVE=2.D0*K2/(3.D0*C*C*C*C*C)
	   WRITE(*,*)YL,YR,K2,LOVE
       WRITE(5,*)IL,(XP/1.D3), MASST,YR,C,LOVE,K2

      ELSE IF ( YL .EQ. 3.D0 ) THEN 

       MASST=MPHY
       YR=Y(3)
       CALL FUNCK3(XP,MASST,YR,K,LAM,C,K3)
       !WRITE(3,*)IL,(XP/1.D3), MASST,YR,C,K3
!      WRITE(3,*)IL,PCC,EDC,(XP/1.D3),Y(2),Y(3)
       WRITE(*,*)YL,YR,K3
       
   
      ELSE IF ( YL .EQ. 4.D0 ) THEN 
       MASST=MPHY
       YR=Y(3)
       CALL FUNCK4(XP,MASST,YR,K,LAM,C,K4)
       !WRITE(4,*)IL,(XP/1.D3), MASST,YR,C,K4
!      WRITE(3,*)IL,PCC,EDC,(XP/1.D3),Y(2),Y(3)
       WRITE(*,*)YL,YR,K4

      END IF
      
      END IF

 20   CONTINUE     
 10   CONTINUE
      
      STOP
      END

      SUBROUTINE FUNCT(EK,J,YA,XA,H,YL,K,LAM)
C     *********************************************************
C     DEFINES TOV EQUATIONS 
C
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z) 
      INTEGER J   
      DIMENSION EK(4,10), YA(10)
      HC  = 197.327D0
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
  

         
C--------------------------------------   
      
      L=LAM*K+1.D0    
      PRESS=YA(1)
      MASST=YA(2)
      YTIDAL=YA(3)
      EDEN=YA(4)
      
C sigma isotropik pressure
     
      A=SQRT(ABS(L+8.D0*PI*GS*K*EDEN))
      B=SQRT(ABS(L-8.D0*PI*GS*K*PRESS))
         
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
     
      !PRINT*,XA/1.D3,PRESS,MASST   !Ilham
      !IF (XA .GT. 1.D0) STOP
      
      
      
     
      C=B
      CS2=1.D0/DEDP(PRESS)
      SGM=0.D0
      DSGMDP=0.D0
      EDENEFF=(A*A-B*B-2.D0*C*C+2.D0*A*B*C*C)
     &       /(16.D0*PI*GS*K*A*B*C*C)
      PRESSEFF=(A*A-B*B+2.D0*C*C-2.D0*A*B*C*C)
     &        /(16.D0*PI*GS*K*A*B*C*C)
      
      XXXX=(1.D0/(B*C*C)+B/(A*A*C*C)-2.D0/(A*A*B))
     &       *1.D0/(4.D0*A*CS2)
     &      +(A/(B*B*C*C)+1.D0/(A*C*C)+2.D0/(A*B*B))
     &       *1.D0/(4.D0*B)
     &      +(A/B-B/A)*1.D0/(2.D0*C*C*C*C)
     &       *(1.D0+DSGMDP)
    
      COBAP=H*(-EDENEFF*(1.D0+PRESSEFF/EDENEFF)/2.D0
     &        *(2.D0*GS*MASST*MSS/(XA*XA)
     &        *(1.D0+4.D0*PI*XA*XA*XA*PRESSEFF/(MSS*MASST)
     &        -LAM*XA*XA/(3.D0*L*GS*MSS*MASST))
     &        /(1.D0-2.D0*GS*MASST*MSS/XA-LAM*XA*XA/(3.D0*L)))
     &        -(2.D0/XA)*SGM/(A*B*C*C))
     &        /XXXX
     
      COBAM=H*4.D0*PI*XA*XA*EDENEFF/MSS
      
C      PRINT *, PRESS, EK(J,1), COBAP, (1.D0-EK(J,1)/COBAP)    ! modif
     
C      PRINT *, PRESS, EK(J,2), COBAM, (1.D0-EK(J,2)/COBAM)    ! modif
      
C Persamaan tidal y
      BETAP=-2.D0*EK(J,1)/H*(2.D0*PI*GS*K
     &      *(3.D0/(B*B)+1.D0/(A*A*CS2))
     &      +1.D0/(EDEN+PRESS))
     
      ALPHA=-LOG(1.D0-2.D0*GS*MSS*MASST/XA-LAM*XA*XA/(3.D0*L))
      FTIDAL=1.D0/XA+EXP(ALPHA)/XA
     &       +EXP(ALPHA)*(XA/K)*(1.D0/(A*B)-1.D0)
      GTIDAL=-(YL*(YL+1.D0)*EXP(ALPHA)/(XA*XA)+2.D0*EXP(ALPHA)/K
     &       +BETAP*BETAP)+2.D0*EXP(ALPHA)*A/(K*B*B*B)
     &       *(2.D0-(4.D0/(A*A-B*B))/(4.D0/(A*A-B*B)
     &       + 3.D0/(B*B) + 1.D0/(A*A*CS2) ))
	 
       IF ( K .EQ. 0.D0 ) THEN 
        PERSP=-H*(GS*MASST*MSS*EDEN/(XA*XA)*(1.D0+PRESS/EDEN)
     &            *(4.D0*PI*XA*XA*XA/(MASST*MSS)
     &              +1.D0-LAM*XA*XA*XA/(3.D0*GS*MASST*MSS))
     &             /(1.D0-2.D0*GS*MSS*MASST/XA))
        PERSM=(H/MSS)*4.D0*PI*XA*XA*EDEN
        BETAP=-2.D0*PERSP/H*(2.D0*PI*GS*K*(3.D0/(B*B)+1.D0/(A*A*CS2))
     &      +1.D0/(EDEN+PRESS))
        FTIDAL=1.D0/XA+EXP(ALPHA)*(1.D0/XA+4.D0*PI*GS*(PRESS-EDEN)*XA)
        GTIDAL=-YL*(YL+1.D0)*EXP(ALPHA)/(XA*XA)+4.D0*PI*GS*EXP(ALPHA)
     &         *(5.D0*EDEN+9.D0*PRESS+(EDEN+PRESS)/CS2)-BETAP*BETAP
        EK(J,1)=PERSP
        EK(J,2)=PERSM
       ENDIF
	  
      EK(J,3)=H*(YTIDAL/XA-YTIDAL*YTIDAL/XA
     &        -FTIDAL*YTIDAL-GTIDAL*XA)

      !WRITE(*,*)XA,YTIDAL
      
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
      


      !====================================================================== Ilham

      SUBROUTINE FUNCK2(XP,MASST,YR,K,LAM,C,K2)
      ! Love number k2 from YL=2
      IMPLICIT DOUBLE PRECISION (A-H,K-Z) 
           
    
      PI  = 3.14159265358979D0      
      GS=1.323790281083862D-12
      MSS=1.1154605324408694D15

      YL=2.D0
      L=LAM*K+1.D0
      EPS=(LAM/L)*(GS*MASST*MSS)**2

      C =GS*MASST*MSS/XP
      C2=C*C
      C3=C2*C

      a12=0.D0
      a22=C2/3.D0
      a32=113.D0*C2/84.D0
      a42=13.D0*C2/9.D0

      b12=15.D0/(8.D0*C3)
      b22=0.D0
      b32=1787.D0/(392.D0*C3)
      b42=-(3305.D0/448.D0+5.D0*PI**2/7.D0+15.D0*DLOG(2.D0)**2/7.D0)/C3

      QB2=YR*Q2(C)+C*DQ2(C)
      PB2=YR*P2(C)+C*DP2(C)
      TB2=YR*T2(C)+C*DT2(C)
      SB2=YR*S2(C)+C*DS2(C)

      KN2=(a12 + eps*a32)* Qb2 + (a22 + eps*a42)* Pb2 
     &   + eps* (a12* Tb2 + a22* Sb2)
      KD2=(b12 + eps*b32)* Qb2 + (b22 + eps*b42)* Pb2
     &   + eps* (b12* Tb2 + b22* Sb2)
     
      K2=-DOUBFAC(2.D0*YL-1.D0)/2*KN2/KD2  

C       K2=-(DOUBFAC(2.D0*YL-1.D0)/2)
C      -    *((a22*Pb2 + a12*Qb2)/(b22*Pb2 + b12*Qb2) + 
C      -    eps*((b22*Pb2 + b12*Qb2)*
C      -    (a42*Pb2 + a32*Qb2 + a22*Sb2 + a12*Tb2) - 
C      -    (a22*Pb2 + a12*Qb2)*(b42*Pb2 + b32*Qb2 + b22*Sb2 + b12*Tb2))/
C      -    (b22*Pb2 + b12*Qb2)**2)

      !WRITE(*,*)'KN2',KN2
      !WRITE(*,*)'KD2',KD2
      !WRITE(*,*)'K2',K2
      
C       PRINT *,C,EPS
C       PRINT *,Q2(C),DQ2(C),P2(C),DP2(C)
C       PRINT *,T2(C),DT2(C),S2(C),DS2(C)
C       PRINT *,QB2,PB2,TB2,SB2
C       PRINT *,KN2,KD2
C       PRINT *,b12,b22,b32,b42
C       PRINT *,(b12+EPS*b32)*QB2,(b22+EPS*b42)*PB2,
C      &         EPS*(b12*TB2+b22*SB2)

      RETURN
      END

      SUBROUTINE FUNCK3(XP,MASST,YR,K,LAM,C,K3)
      ! Love number k3 from YL=3
      IMPLICIT DOUBLE PRECISION (A-H,K-Z) 
           
    
      PI  = 3.14159265358979D0      
      GS=1.323790281083862D-12
      MSS=1.1154605324408694D15

      YL=3.D0
      L=LAM*K+1.D0
      EPS=(LAM/L)*(GS*MASST*MSS)**2

      C =GS*MASST*MSS/XP
      C2=C*C
      C3=C2*C
      C4=C3*C

      a13=0.D0
      a23=C3/45.D0
      a33=127.D0*C3/756.D0
      a43=158.D0*C3/945.D0

      b13=35.D0/(8.D0*C4)
      b23=0.D0
      b33=24805.D0/(1512.D0*C4)
      b43=-(13795.D0/576.D0+25.D0*PI**2/9.D0
     &    +25.D0*DLOG(2.D0)**2/3.D0)/C4

      QB3=YR*Q3(C)+C*DQ3(C)
      PB3=YR*P3(C)+C*DP3(C)
      TB3=YR*T3(C)+C*DT3(C)
      SB3=YR*S3(C)+C*DS3(C)

      KN3=(a13+EPS*a33)*QB3 + (a23+EPS*a43)*PB3
     &    +EPS*(a13*TB3+a23*SB3)
      KD3=(b13+EPS*b33)*QB3 + (b23+EPS*b43)*PB3
     &    +EPS*(b13*TB3+b23*SB3)
     
      K3=-DOUBFAC(2.D0*YL-1.D0)/2*KN3/KD3  

C       K3=-(DOUBFAC(2.D0*YL-1.D0)/2)
C      -    *((a23*Pb3 + a13*Qb3)/(b23*Pb3 + b13*Qb3) + 
C      -  (eps*((b23*Pb3 + b13*Qb3)*
C      -        (a43*Pb3 + a33*Qb3 + a23*Sb3 + a13*Tb3)
C      -        - (a23*Pb3 + a13*Qb3)*
C      -        (b43*Pb3 + b33*Qb3 + b23*Sb3 + b13*Tb3)))
C      -    /(b23*Pb3 + b13*Qb3)**2)

      !WRITE(*,*)'KN3',KN3
      !WRITE(*,*)'KD3',KD3
      !WRITE(*,*)'K3',K3
      
C       PRINT *,C,EPS
C       PRINT *,Q3(C),DQ3(C),P3(C),DP3(C)
C       PRINT *,T3(C),DT3(C),S3(C),DS3(C)
C       PRINT *,QB3,PB3,TB3,SB3
C       PRINT *,KN3,KD3
C       PRINT *,b13,b23,b33,b43
C       PRINT *,(b13+EPS*b33)*QB3,(b23+EPS*b43)*PB3,
C      &         EPS*(b13*TB3+b23*SB3)

      RETURN
      END
    
      SUBROUTINE FUNCK4(XP,MASST,YR,K,LAM,C,K4)
        ! Love number k4 from YL=4
      IMPLICIT DOUBLE PRECISION (A-H,K-Z) 
             
      
      PI  = 3.14159265358979D0      
      GS=1.323790281083862D-12
      MSS=1.1154605324408694D15
  
      YL=4.D0
      L=LAM*K+1.D0
      EPS=(LAM/L)*(GS*MASST*MSS)**2
  
      C =GS*MASST*MSS/XP
      C2=C*C
      C3=C2*C
      C4=C3*C
      C5=C4*C
  
      a14=0.D0
      a24=C4/630.D0
      a34=7613.D0*C4/465696.D0
      a44=200077.D0*C4/10866240.D0
  
      b14=735.D0/(64.D0*C5)
      b24=0.D0
      b34=469685.D0/(5808.D0*C5)
      b44=-(1.D0/C5)
     &    *(14680085.D0/202752.D0
     &      +875.D0*PI*PI/88.D0
     &      +2625.D0*LOG(2.D0)*LOG(2.D0)/88.D0)
  
      QB4=YR*Q4(C)+C*DQ4(C)
      PB4=YR*P4(C)+C*DP4(C)
      TB4=YR*T4(C)+C*DT4(C)
      SB4=YR*S4(C)+C*DS4(C)
  
      KN4=(a14+EPS*a34)*QB4 + (a24+EPS*a44)*PB4
     &    +EPS*(a14*TB4+a24*SB4)
      KD4=(b14+EPS*b34)*QB4 + (b24+EPS*b44)*PB4
     &    +EPS*(b14*TB4+b24*SB4)

      K4=-DOUBFAC(2.D0*YL-1.D0)/2*KN4/KD4  
      
C       K4=-(DOUBFAC(2.D0*YL-1.D0)/2)
C      -    *((a24*Pb4 + a14*Qb4)/(b24*Pb4 + b14*Qb4) + 
C      -  (eps*((b24*Pb4 + b14*Qb4)*
C      -        (a44*Pb4 + a34*Qb4 + a24*Sb4 + a14*Tb4)
C      -        - (a24*Pb4 + a14*Qb4)*
C      -        (b44*Pb4 + b34*Qb4 + b24*Sb4 + b14*Tb4)))
C      -    /(b24*Pb4 + b14*Qb4)**2)

      !WRITE(*,*)'KN4',KN4
      !WRITE(*,*)'KD4',KD4
      !WRITE(*,*)'K4',K4
      
C       PRINT *,C,EPS
C       PRINT *,Q4(C),DQ4(C),P4(C),DP4(C)
C       PRINT *,T4(C),DT4(C),S4(C),DS4(C)
C       PRINT *,QB4,PB4,TB4,SB4
C       PRINT *,KN4,KD4
C       PRINT *,b14,b24,b34,b44
C       PRINT *,(b14+EPS*b34)*QB4,(b24+EPS*b44)*PB4,
C      &         EPS*(b14*TB4+b24*SB4)
  
      RETURN
      END

      FUNCTION PLG(C)
      ! PolyLogarithm PolyLog(2,(1-x)/2)) 
	  ! or Li_2 ((1-x)/2) with x=1/c-1
       IMPLICIT NONE
       DOUBLE PRECISION C, x, PLG
        x = C
        IF ( x .GT. 0.01D0 .AND. x .LE. 0.5D0 ) THEN
            PLG = -4.900128181312306 + 
     -            0.00016932349271285354/x**2 - 
     -            0.06607300001104759/x + 
     -            67.81366827244638*x - 
     -            661.3708294454829*x**2 + 
     -            4714.740202946316*x**3 - 
     -            23715.678000133867*x**4 + 
     -            84334.3662190509*x**5 - 
     -            213246.7648131581*x**6 + 
     -            383601.9350809346*x**7 - 
     -            486336.6093429251*x**8 + 
     -            424030.9280758699*x**9 - 
     -            241718.87649564294*x**10 + 
     -            81048.35624350267*x**11 - 
     -            12113.292172158068*x**12
        ELSE IF ( x .GT. 0.D0 .AND. x .LE. 0.01D0 ) THEN
            PLG = 5894.173801180173 - 
     -            5875.9257455803145/
     -             x**0.001
        ELSE
		    PLG = 0.D0
        END IF
		
       !WRITE(*,*)'PLG',PLG
       RETURN
      END

      FUNCTION DOUBFAC(Z)
      ! Double factorial: n!!
      IMPLICIT NONE
      DOUBLE PRECISION Z, Y, DOUBFAC
      DOUBFAC = Z
      Y=Z-2.D0
      DO WHILE (Y.GT.0.D0)
        DOUBFAC=DOUBFAC*Y
        Y=Y-2.D0
      ENDDO
      !WRITE(*,*)'Doubfac',DOUBFAC
      RETURN
      END

      FUNCTION Q2(C)
        ! LegendreQ(2,2,3,XX)
      IMPLICIT NONE
      DOUBLE PRECISION Q2,C,XX
      XX=-1.D0+1.D0/C
      Q2= 6.D0*(-1.D0+XX)*(1.D0+XX)*(-(-1.D0+XX)**2/(48.D0*(1.D0+XX)**2)  
     &    +(-1.D0+XX)/(6.D0*(1.D0+XX))-(1.D0+XX)/(6.D0*(-1.D0+XX))  
     &    +(1.D0+XX)**2/(48.D0*(-1.D0+XX)**2)+(-DLOG(-1.D0+XX) 
     &    +DLOG(1.D0 + XX))/4.D0 )  
      !WRITE(*,*)'Q2',Q2
      RETURN
      END

      FUNCTION Q3(C)
        ! LegendreQ(3,2,3,XX)
      IMPLICIT NONE
      DOUBLE PRECISION Q3,C,XX
      XX=-1.D0+1.D0/C
      Q3= 45.D0*(-1.D0 + XX)*(1.D0 + XX)**2
     &   *(-(-1.D0 + XX)**3/(720.D0*(1.D0 + XX)**3) 
     &   +(-1.D0 + XX)**2/(48.D0*(1.D0 + XX)**2) 
     &   -(1.D0 + XX)/(48.D0*(-1.D0 + XX)) 
     &   +(1.D0 + XX)**2/(720.D0*(-1.D0 + XX)**2) 
     &   +(-4.D0/3.D0 - DLog(-1.D0 + XX) 
     &   +   DLog(1 + XX))/12.D0 
     &   +((-1.D0 + XX)
     &   *   (4.D0/3.D0 - DLog(-1.D0 + XX) 
     &   +     DLog(1.D0 + XX)))/(12.D0*(1.D0 + XX)))
      !WRITE(*,*)'Q3',Q3
      RETURN
      END

      FUNCTION Q4(C)
        ! LegendreQ(4,2,3,XX)
      IMPLICIT NONE
      DOUBLE PRECISION Q4,C,XX
      XX=-1.D0+1.D0/C
      Q4=540.D0*(-1.D0 + XX)*(1.D0 + XX)**3
     &   *(-(-1.D0 + XX)**4/(17280.D0*(1.D0 + XX)**4)
     &   + (-1.D0 + XX)**3/(720.D0*(1.D0 + XX)**3) 
     &   - (1.D0 + XX)/(720.D0*(-1.D0 + XX))  
     &   + (1.D0 + XX)**2/(17280.D0*(-1.D0 + XX)**2)  
     &   +  (-25.D0/12.D0 - Log(-1.D0 + XX)  
     &   +    Log(1.D0 + XX))/96.D0  
     &   + ((-1.D0 + XX)*(-Log(-1.D0 + XX) + Log(1.D0 + XX)))
     &   /  (36.D0*(1.D0 + XX))  
     &   + ((-1.D0 + XX)**2
     &   *    (25.D0/12.D0 - Log(-1.D0 + XX) 
     &   +      Log(1.D0 + XX)))/(96.D0*(1.D0 + XX)**2))
      !WRITE(*,*)'Q4',Q4
      RETURN
      END

      FUNCTION DQ2(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL Q2
      h  = 1.D-5
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
      DQ2 = (Q2(x4)-8.D0*Q2(x2)+8.D0*Q2(x1)-Q2(x3))/(12.D0*h)
      !WRITE(*,*)'Q2',Q2(xa)
      !WRITE(*,*)'DQ2',DQ2
      RETURN
      END

      FUNCTION DQ3(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL Q3
      h  = 1.D-5
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
      DQ3 = (Q3(x4)-8.D0*Q3(x2)+8.D0*Q3(x1)-Q3(x3))/(12.D0*h)
      !WRITE(*,*)'Q3',Q3(xa)
      !WRITE(*,*)'DQ3',DQ3
      RETURN
      END
      
      FUNCTION DQ4(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL Q4
      h  = 1.D-5
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
      DQ4 = (Q4(x4)-8.D0*Q4(x2)+8.D0*Q4(x1)-Q4(x3))/(12.D0*h)
      !WRITE(*,*)'Q4',Q4(xa)
      !WRITE(*,*)'DQ4',DQ4
      RETURN
      END

      FUNCTION P2(C)
        ! LegendreP(2,2,3,XX)
      IMPLICIT NONE
      DOUBLE PRECISION P2,C,XX
      XX=-1.D0+1.D0/C
      P2=(3.D0*(1.D0 - XX)**2*(1.D0 + XX))/(-1.D0 + XX)
      !WRITE(*,*)'P2',P2
      RETURN
      END

      FUNCTION P3(C)
        ! LegendreP(3,2,3,XX)
      IMPLICIT NONE
      DOUBLE PRECISION P3,C,XX
      XX=-1.D0+1.D0/C
      P3=(15.D0*(1.D0 - XX)**2*XX*(1.D0 + XX))/(-1.D0 + XX)
      !WRITE(*,*)'P3',P3
      RETURN
      END

      FUNCTION P4(C)
        ! LegendreP(4,2,3,XX)
      IMPLICIT NONE
      DOUBLE PRECISION P4,C,XX
      XX=-1.D0+1.D0/C
      P4=(15.D0*(1.D0 - XX)**2*(1.D0 + XX)*(-1.D0 + 7.D0*XX**2))
     &   /(2.D0*(-1.D0 + XX))
      !WRITE(*,*)'P4',P4
      RETURN
      END

      FUNCTION DP2(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL P2
      h  = 1.D-5
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
      DP2 = (P2(x4)-8.D0*P2(x2)+8.D0*P2(x1)-P2(x3))/(12.D0*h)
      !WRITE(*,*)'P2',P2(xa)
      !WRITE(*,*)'DP2',DP2
      RETURN
      END

      FUNCTION DP3(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL P3
      h  = 1.D-5
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
      DP3 = (P3(x4)-8.D0*P3(x2)+8.D0*P3(x1)-P3(x3))/(12.D0*h)
      !WRITE(*,*)'P3',P3(xa)
      !WRITE(*,*)'DP3',DP3
      RETURN
      END
        
      FUNCTION DP4(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL P4
      h  = 1.D-5
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
      DP4 = (P4(x4)-8.D0*P4(x2)+8.D0*P4(x1)-P4(x3))/(12.D0*h)
      !WRITE(*,*)'P4',P4(xa)
      !WRITE(*,*)'DP4',DP4
      RETURN
      END


      FUNCTION T2(C)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL PLG
      XX=-1.D0+1.D0/C
      F32=235.D0/56.D0 + (1357.D0*XX)/56.D0 
     &    - (1469.D0*XX**2)/84.D0 - (1987.D0*XX**3)/84.D0 
     &    + (2057.D0*XX**4)/168.D0 + (389.D0*XX**5)/56.D0 
     &    - (8.D0*XX**6)/7.D0
      F42=-(25.D0/14.D0) + (153.D0*XX)/14.D0 + (171.D0*XX**2)/14.D0
     &    - (17.D0*XX**3)/14.D0 - (25.D0*XX**4)/7.D0
     &    - (4.D0*XX**5)/7.D0
      F52=-(25.D0/14.D0) + (34.D0*XX)/7.D0 - (93.D0*XX**2)/14.D0 
     &    + (13.D0*XX**3)/7.D0 + (4.D0*XX**4)/7.D0
      F62=-(24.D0/7.D0)*(DLOG(XX-1.D0)*DLOG((XX+1.D0)/4.D0)
     &    +2.D0*PLG(C))
      T2=F32/((XX+1.D0)*(XX-1.D0)**2)+F42/(XX+1.D0)*DLOG(XX-1.D0)
     &   +F52*(XX+1.D0)/(XX-1.D0)*DLOG(XX+1.D0)+(XX**2-1.D0)*F62
      !WRITE(*,*)'T2',T2
      RETURN
      END

      FUNCTION T3(C)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL PLG
      XX=-1.D0+1.D0/C
      F33=- 32.D0 + (8941.D0*XX)/168.D0 + (7305.D0*XX**2)/56.D0 
     &   -(9977.D0*XX**3)/84.D0 - (35845.D0*XX**4)/252.D0 
     &   +(5225.D0*XX**5)/72.D0 + (8195.D0*XX**6)/168.D0 
     &   -(20.D0*XX**7)/3.D0
      F43=- 635.D0/21.D0 - (845.D0*XX)/42.D0 
     &    +(195.D0*XX**2)/2.D0 + (3425.D0*XX**3)/42.D0 
     &    -(925.D0*XX**4)/42.D0 - (70.D0*XX**5)/3.D0 
     &    -(10.D0*XX**6)/3.D0
      F53=5.D0/21.D0 + (135.D0*XX)/14.D0 + (20.D0*XX**2)/7.D0 
     &    - (1315.D0*XX**3)/42.D0 + (40.D0*XX**4)/3.D0 
     &    + (10.D0*XX**5)/3.D0
      F63=-(200.D0/7.D0)*XX*(DLog(XX- 1.D0)*DLog((1.D0 + XX)/4.D0) 
     &    + 2.D0*PLG(C))
      T3=F33/((XX+1.D0)*(XX-1.D0)**2)+F43*DLOG(XX-1.D0)/(XX+1.D0)
     &   +F53*DLOG(XX+1.D0)*(XX+1.D0)/(XX-1.D0)+(XX**2-1.D0)*F63
      !WRITE(*,*)'T3',T3
      RETURN
      END

      FUNCTION T4(C)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL PLG
      XX=-1.D0+1.D0/C
      F34=-(2160637.D0/44352.D0) - (7430627.D0*XX)/44352.D0 
     &    + (169207.D0*XX**2)/448.D0 
     &    + (7841513.D0*XX**3)/14784.D0 
     &    - (3123931.D0*XX**4)/4928.D0 
     &    - (8024815.D0*XX**5)/14784.D0 
     &    + (14823115.D0*XX**6)/44352.D0 
     &    + (1159715.D0*XX**7)/6336.D0 
     &    - (595.D0*XX**8)/22.D0
      F44=-(3735.D0/308.D0) - (6715.D0*XX)/28.D0 
     &    - (42865.D0*XX**2)/308.D0 + (154025.D0*XX**3)/308.D0 
     &    + (113875.D0*XX**4)/308.D0 - (38235.D0*XX**5)/308.D0 
     &    - (4445.D0*XX**6)/44.D0 - (595.D0*XX**7)/44.D0
      F54=-(3735.D0/308.D0) + (535.D0*XX)/77.D0 
     &    + (30455.D0*XX**2)/308.D0 - (4540.D0*XX**3)/77.D0 
     &    - (34285.D0*XX**4)/308.D0 + (665.D0*XX**5)/11.D0 
     &    + (595.D0*XX**6)/44.D0
      F64=-(1500.D0/77.D0)*(-1.D0 + 7.D0*XX**2)
     &    *(DLog(-1.D0 + XX)*DLog((1.D0 + XX)/4) 
     &    + 2.D0*PLG(C))
      T4=F34/((XX+1.D0)*(XX-1.D0)**2)+F44*DLOG(XX-1.D0)/(XX+1.D0)
     &   +F54*DLOG(XX+1.D0)*(XX+1.D0)/(XX-1.D0)+(XX**2-1.D0)*F64
      !WRITE(*,*)'T4',T4
      RETURN
      END
      
      
      FUNCTION DT2(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL T2
      h  = 1.D-5
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
      DT2 = (T2(x4)-8.D0*T2(x2)+8.D0*T2(x1)-T2(x3))/(12.D0*h)
      !WRITE(*,*)'T2',T2(xa)
      !WRITE(*,*)'DT2',DT2
      RETURN
      END

      FUNCTION DT3(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL T3
      h  = 1.D-5
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
      DT3 = (T3(x4)-8.D0*T3(x2)+8.D0*T3(x1)-T3(x3))/(12.D0*h)
      !WRITE(*,*)'T3',T3(xa)
      !WRITE(*,*)'DT3',DT3
      RETURN
      END

      FUNCTION DT4(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL T4
      h  = 1.D-5
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
      DT4 = (T4(x4)-8.D0*T4(x2)+8.D0*T4(x1)-T4(x3))/(12.D0*h)
      !WRITE(*,*)'T4',T4(xa)
      !WRITE(*,*)'DT4',DT4
      RETURN
      END


      FUNCTION S2(C)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      XX=-1.D0+1.D0/C
      F12=1.D0+21.D0*XX/4.D0-46.D0*XX**2/7.D0-59.D0*XX**3/4.D0
     &    -XX**4/7.D0+6.D0*XX**5+8.D0*XX**6/7.D0
      F22=(3.D0/56.D0)*(113.D0*DLOG(XX-1.D0)+15.D0*DLOG(XX+1.D0))
      S2=F12/(XX**2-1.D0)+(XX**2-1.D0)*F22 
      !WRITE(*,*)'S2',S2
      RETURN
      END

      FUNCTION S3(C)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      XX=-1.D0+1.D0/C
      F13=(20.D0*XX)+(185.D0*XX**2)/4.D0-(220.D0*XX**3)/7.D0
     &    -(375.D0*XX**4)/4.D0-(20.D0*XX**5)/7.D0
     &    +40.D0*XX**6+(20.D0*XX**7)/3.D0
      F23=(25.D0/56.D0)*XX*(127.D0*LOG(XX-1.D0)+LOG(XX+1.D0))
      S3=F13/(XX**2-1.D0)+(XX**2-1.D0)*F23 
      !WRITE(*,*)'S3',S3
      RETURN
      END

      FUNCTION S4(C)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      XX=-1.D0+1.D0/C
      F14=-615.D0/64.D0-(1225.D0/32.D0)*XX+(342515.D0/4928.D0)*XX**2
     &    +(14405.D0/48.D0)*XX**3-(202325.D0/4928.D0)*XX**4
     &    -(14305.D0/32.D0)*XX**5-(277315.D0/4928.D0)*XX**6
     &    +175.D0*XX**7+(595.D0/22.D0)*XX**8
      F24=(25.D0/4928.D0)*(7.D0*XX**2-1.D0)
     &    *(7613.D0*LOG(XX-1.D0)+67.D0*LOG(XX+1.D0))
      S4=F14/(XX**2-1.D0)+(XX**2-1.D0)*F24 
      !WRITE(*,*)'S4',S4
      RETURN
      END

      FUNCTION DS2(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL S2
      h  = 1.D-5
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
      DS2 = (S2(x4)-8.D0*S2(x2)+8.D0*S2(x1)-S2(x3))/(12.D0*h)
      !WRITE(*,*)'S2',S2(xa)
      !WRITE(*,*)'DS2',DS2
      RETURN
      END

      FUNCTION DS3(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL S3
      h  = 1.D-5
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
      DS3 = (S3(x4)-8.D0*S3(x2)+8.D0*S3(x1)-S3(x3))/(12.D0*h)
      !WRITE(*,*)'S3',S3(xa)
      !WRITE(*,*)'DS3',DS3
      RETURN
      END

      FUNCTION DS4(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL S4
      h  = 1.D-5
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
      DS4 = (S4(x4)-8.D0*S4(x2)+8.D0*S4(x1)-S4(x3))/(12.D0*h)
      !WRITE(*,*)'S4',S4(xa)
      !WRITE(*,*)'DS4',DS4
      RETURN
      END

C ==============================================================================================

      SUBROUTINE FUNCM(XA,YA,MPHY,K,LAM)
C     *********************************************************
C     Definition for physical mass MPHY
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)   
      DIMENSION YA(10)
      DOUBLE PRECISION LAM
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
      
      SGM=0.D0

      A=DSQRT(ABS(L+8.D0*PI*GS*K*EDEN))
      B=DSQRT(ABS(L-8.D0*PI*GS*K*PRESS))
      C=DSQRT(ABS(L-8.D0*PI*GS*K*(PRESS-SGM)))
      
C      MPHY=XA/(2.D0*GS*MSS*DSQRT(A*B))
C     &     *(1.D0+A*C*C/(B*B)*(2.D0*GS*MASST*MSS/XA-1.D0)
C     &     +(LAM*XA*XA/(3.D0*L))*(A*C*C/B-1.D0/(A*B)))
      MPHY=MASST

   
      RETURN
      END
      
C===============================================================================================
