************************************************************************
*                                                                      *
*         ESTA SUBROTINA CALCULA AS MATRIZES [H] E [G];                *
*                                                                      *
*        AINDA MONTA O SISTEMA DE EQUAÇÕES [A] {x} = {f}               *
*                                                                      *
************************************************************************
*
      SUBROUTINE GHMATECE (CX,CY,CZ,CXM,CYM,CZM,HEST,GEST,NE,NX,NCOX,
     $CONE,N,NBE,NP,NPG,GE,RNU,RMU,DELTA,PI,C1,C2,C3,C4)
*
	IMPLICIT REAL*8 (A-H,O-Y)
	IMPLICIT COMPLEX*16 (Z)
	COMMON INP,INQ,IPR,IPS,IPT
      DIMENSION CX(NCOX),CY(NCOX),CZ(NCOX)
      DIMENSION CXM(NE),CYM(NE),CZM(NE)
	DIMENSION HEST(NX,NX),GEST(NX,NX)
	DIMENSION HELEM(3,3),GELEM(3,3)
	DIMENSION DELTA(3,3)
	DIMENSION CO(4,3),ETA(3)						
	INTEGER CONE(NE,4)
*
      PI=4.D0*DATAN(1.D0)
*
* CONSTANTES USADAS NAS SOLUÇÕES FUNDAMENTAIS ESTÁTICAS
*
      C1=1.D0/(16.D0*PI*RMU*(1.D0-RNU))
      C2=3.D0-(4.D0*RNU)
      C3=-1.D0/(8.D0*PI*(1.D0-RNU))
      C4=1.D0-(2.D0*RNU)
*
      DO 25 I=1,3
      DO 25 J=1,3
      IF(I-J) 10,20,10
  10  DELTA(I,J)=0.D0
      GOTO 25
  20  DELTA(I,J)=1.D0
  25  CONTINUE
*
* ZERANDO AS MATRIZES H E G
*
      DO 50 J=1,N
      DO 50 I=1,N
      II=3*(I-1)
      DO 50 JN=1,3
      DO 50 IN=1,3
	GEST(II+IN,II+JN)=0.D0
  50  HEST(II+IN,II+JN)=0.D0
*
* CÁLCULO DOS COEFICIENTES DAS MATRIZES H E G
*
      DO 220 J=1,N
*
* CÁLCULO DAS COMPONENTES DO VETOR NORMAL
* USANDO O PRODUTO VETORIAL DOS LADOS 1-2 E 1-3
*
      N1=CONE(J,1)
      N2=CONE(J,2)
      N3=CONE(J,3)
      N4=CONE(J,4)
      A=(CY(N2)-CY(N1))*(CZ(N3)-CZ(N1))-(CZ(N2)-CZ(N1))*(CY(N3)-CY(N1))
      B=(CZ(N2)-CZ(N1))*(CX(N3)-CX(N1))-(CX(N2)-CX(N1))*(CZ(N3)-CZ(N1))
      C=(CX(N2)-CX(N1))*(CY(N3)-CY(N1))-(CY(N2)-CY(N1))*(CX(N3)-CX(N1))
      R=DSQRT(A*A+B*B+C*C)
      ETA(1)=A/R
      ETA(2)=B/R
      ETA(3)=C/R
*
* ARMAZENA AS COORDENADAS DOS PONTOS EXTREMOS DO ELEMENTO NO VETOR CO
*
      CO(1,1)=CX(N1)
      CO(1,2)=CY(N1)
      CO(1,3)=CZ(N1)
      CO(2,1)=CX(N2)
      CO(2,2)=CY(N2)
      CO(2,3)=CZ(N2)
      CO(3,1)=CX(N3)
      CO(3,2)=CY(N3)
      CO(3,3)=CZ(N3)
      CO(4,1)=CX(N4)
      CO(4,2)=CY(N4)
      CO(4,3)=CZ(N4)
	JJ=3*(J-1)
*
	DO 200 I=1,NBE
*	 
	IF(I-J)120,140,120
*
* ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE H E G NÃO SINGULAR
*
 120	CALL NONSINGE(HELEM,GELEM,CO,CXM(I),CYM(I),CZM(I),ETA,
     $N,NP,NPG,GE,RNU,RMU,C1,C2,C3,C4,DELTA)
	GO TO 150
*
* ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE G SINGULAR
*
 140	CALL SINGGE(HELEM,GELEM,CXM(I),CYM(I),CZM(I),ETA,CX,CY,CZ,
     $N1,N2,N3,N4,NCOX,N,NP,NPG,GE,RNU,RMU,C1,C2,C3,C4,DELTA)
 150	CONTINUE 
*
	II=3*(I-1)
      DO 190 JN=1,3
      DO 190 IN=1,3
	GEST(II+IN,JJ+JN)=GELEM(IN,JN)
	HEST(II+IN,JJ+JN)=HELEM(IN,JN)
 190  CONTINUE
*
 200  CONTINUE
 220  CONTINUE
*
* ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE H SINGULAR ATRAVÉS
* DA CONSIDERAÇÃO DO MOVIMENTO DO CORPO RÍGIDO
*
	DO 224 MA=1,NBE
	HEST(3*MA-2,3*MA-2)=0.D0
	HEST(3*MA-2,3*MA-1)=0.D0
	HEST(3*MA-2,3*MA)  =0.D0
	HEST(3*MA-1,3*MA-2)=0.D0
	HEST(3*MA-1,3*MA-1)=0.D0
	HEST(3*MA-1,3*MA)  =0.D0
	HEST(3*MA,3*MA-2)  =0.D0
	HEST(3*MA,3*MA-1)  =0.D0
	HEST(3*MA,3*MA)    =0.D0
	DO 227 MB=1,N
	IF (MA-MB)225,227,225
 225	HEST(3*MA-2,3*MA-2)=HEST(3*MA-2,3*MA-2)-HEST(3*MA-2,3*MB-2)
	HEST(3*MA-2,3*MA-1)=HEST(3*MA-2,3*MA-1)-HEST(3*MA-2,3*MB-1)
	HEST(3*MA-2,3*MA)  =HEST(3*MA-2,3*MA)  -HEST(3*MA-2,3*MB)
	HEST(3*MA-1,3*MA-2)=HEST(3*MA-1,3*MA-2)-HEST(3*MA-1,3*MB-2)
	HEST(3*MA-1,3*MA-1)=HEST(3*MA-1,3*MA-1)-HEST(3*MA-1,3*MB-1)
	HEST(3*MA-1,3*MA)  =HEST(3*MA-1,3*MA)  -HEST(3*MA-1,3*MB)
	HEST(3*MA,3*MA-2)  =HEST(3*MA,3*MA-2)  -HEST(3*MA,3*MB-2)
	HEST(3*MA,3*MA-1)  =HEST(3*MA,3*MA-1)  -HEST(3*MA,3*MB-1)
	HEST(3*MA,3*MA)    =HEST(3*MA,3*MA)    -HEST(3*MA,3*MB)
 227  CONTINUE
 224	CONTINUE
*
	RETURN
	END