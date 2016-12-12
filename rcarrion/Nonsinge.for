************************************************************************
*                                                                      *
*  SUBROTINA QUE CALCULA OS COEFICIENTES DE [H] E [G]. O PONTO DE		 *
*																	 *
*    COLOCA��O EST� FORA DO CONTORNO, ASSIM N�O H� SINGULARIDADE.		 *
*                                                                      *
*         A INTEGRA��O � FEITA USANDO QUADRATURA GAUSSIANA             *
*                                                                      *
************************************************************************
*
      SUBROUTINE NONSINGE(HELEM,GELEM,CO,CXP,CYP,CZP,ETA,
     $N,NP,NPG,GE,RNU,RMU,C1,C2,C3,C4,DELTA)
*
	IMPLICIT REAL*8 (A-H,O-Y)
	IMPLICIT COMPLEX*16 (Z)
	COMMON INP,INQ,IPR,IPS,IPT
	DIMENSION DELTA(3,3)
	DIMENSION HELEM(3,3),GELEM(3,3),U(3,3),T(3,3)
	DIMENSION GI(NPG),OME(NPG)
	DIMENSION CO(4,3),ETA(3),P(2,4),XJ(2,3),F(4)
*
* ZERA AS MATRIZES ELEMENTARES HELEM E GELEM
*
	DO 20 JN=1,3
	DO 20 IN=1,3
	HELEM(IN,JN)=0.D0
 	GELEM(IN,JN)=0.D0
  20	CONTINUE
*
* ACIONA ROTINA QUE CALCULA OS PONTOS E PESOS DE GAUSS
*
	CALL GAULEG(-1.D0,1.D0,GI,OME,NPG)
*
	DO 400 JG=1,NPG
*
      G2=GI(JG)
      P2=OME(JG)
      SP=1.D0+G2
      SM=1.D0-G2
      P(1,1)=-0.25D0*SM
      P(1,2)= 0.25D0*SM
      P(1,3)= 0.25D0*SP
      P(1,4)=-0.25D0*SP
*
      DO 300 IG=1,NPG
*
      G1=GI(IG)
      P1=OME(IG)
      RP=1.D0+G1
      RM=1.D0-G1
      F(1)=0.25D0*RM*SM
      F(2)=0.25D0*RP*SM
      F(3)=0.25D0*RP*SP
      F(4)=0.25D0*RM*SP
      P(2,1)=-0.25D0*RM
      P(2,2)=-0.25D0*RP
      P(2,3)= 0.25D0*RP
      P(2,4)= 0.25D0*RM
*
* CALCULA A RELA��O ENTRE AS COORDENADAS CARTESIANAS E HOMOG�NEAS
*
      DO 90 I=1,2
      DO 80 J=1,3
      TEMP=0.D0
      DO 70 K=1,4
      TEMP=TEMP+P(I,K)*CO(K,J)
  70  CONTINUE
      XJ(I,J)=TEMP
  80  CONTINUE
  90  CONTINUE
*
* CALCULA O JACOBIANO
*
      DET=DSQRT((XJ(1,2)*XJ(2,3)-XJ(2,2)*XJ(1,3))**2 +
     $          (XJ(2,1)*XJ(1,3)-XJ(1,1)*XJ(2,3))**2 +
     $          (XJ(1,1)*XJ(2,2)-XJ(2,1)*XJ(1,2))**2  )
      IF(DET-1.0D-5) 100,110,110
 100  WRITE(IPR,1000) DET
      WRITE(IPR,1100) ((CO(I,J),J=1,3),I=1,4)
 1000 FORMAT(///' NONSING : ERRO, JACOBIANO NULO OU NEGATIVO =',D14.5)
 1100 FORMAT(/1X,'COORDINADAS DOS PONTOS EXTREMOS '/(3D14.5))
      STOP ' NONSING : ERRO, JACOBIANO NULO OU NEGATIVO'
 110  CONTINUE
*
* CALCULA AS COORDENADAS DO PONTO DE INTEGRA��O
*
      CXG=0.D0
      CYG=0.D0
      CZG=0.D0
      DO 130 I=1,4
      CXG=CXG+CO(I,1)*F(I)
      CYG=CYG+CO(I,2)*F(I)
      CZG=CZG+CO(I,3)*F(I)
 130  CONTINUE
*
* ACIONA ROTINA QUE CALCULA A SOLU��O FUNDAMENTAL EST�TICA 3D
*
      CALL SOLFUNE (U,T,CXP,CYP,CZP,CXG,CYG,CZG,ETA,
     $N,NP,NPG,GE,RNU,RMU,C1,C2,C3,C4,DELTA)
*
      P12=P1*P2*DET
      DO 200 JN=1,3
      DO 200 IN=1,3
	HELEM(IN,JN)=HELEM(IN,JN)+T(IN,JN)*P12
      GELEM(IN,JN)=GELEM(IN,JN)+U(IN,JN)*P12
 200  CONTINUE
*
 300  CONTINUE
 400  CONTINUE
*
      RETURN
      END