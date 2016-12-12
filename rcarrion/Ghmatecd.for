************************************************************************
*                                                                      *
*         ESTA SUBROTINA CALCULA AS MATRIZES [H] E [G];                *
*                                                                      *
*        AINDA MONTA O SISTEMA DE EQUA��ES [A] {x} = {f}               *
*                                                                      *
************************************************************************
*
      SUBROUTINE GHMATECD (CX,CY,CZ,CXM,CYM,CZM,HEST,GEST,ZH,ZG,ZFI,DFI,
     $ZDFI,KODE,NE,NX,NCOX,CONE,DELTA,PI,N,NBE,NP,NPG,GE,RNU,RMU,
     $L,FR,DAM,RHO,ZGE,ZCS,ZCP,C1,C2,C3,C4)
	IMPLICIT REAL*8 (A-H,O-Y)
	IMPLICIT COMPLEX*16 (Z)
	COMMON INP,INQ,IPR,IPS,IPT
      DIMENSION CX(NCOX),CY(NCOX),CZ(NCOX)
      DIMENSION CXM(NE),CYM(NE),CZM(NE)
	DIMENSION HEST(NX,NX),GEST(NX,NX)
	DIMENSION ZH(NX,NX),ZG(NX,NX)
	DIMENSION ZHELEM(3,3),ZGELEM(3,3)
C
C quando coloco NX o programa n�o compila
	DIMENSION ZHP(2880,2880),ZGP(2880,2880)
	DIMENSION ZHEST(2880,2880),ZGEST(2880,2880)
C
	DIMENSION CO(4,3),ETA(3)						
	DIMENSION ZFI(NX),ZDFI(NX)
	DIMENSION DFI(NX)
	DIMENSION DELTA(3,3)
	INTEGER CONE(NE,4),KODE(NX)
*
* TRANSFORMA��O DAS CONDI��ES DE CONTORNO EM N�MEROS COMPLEXOS
*
      DO 2 I=1,NBE
      ZDFI(3*I-2)=DCMPLX(DFI(3*I-2),0.D0)
      ZDFI(3*I-1)=DCMPLX(DFI(3*I-1),0.D0)
	ZDFI(3*I)=DCMPLX(DFI(3*I),0.D0)
 2    CONTINUE
*
* TRANSFORMA��O DAS MATRIZES [HEST] E [GEST] EM N�MEROS COMPLEXOS
*
	DO 5 I=1,N
	DO 5 J=1,N
	ZHEST((3*I-2),(3*J-2))=DCMPLX(HEST((3*I-2),(3*J-2)),0.D0)
	ZHEST((3*I-2),(3*J-1))=DCMPLX(HEST((3*I-2),(3*J-1)),0.D0)
	ZHEST((3*I-2),(3*J))  =DCMPLX(HEST((3*I-2),(3*J)),0.D0)
	ZHEST((3*I-1),(3*J-2))=DCMPLX(HEST((3*I-1),(3*J-2)),0.D0)
	ZHEST((3*I-1),(3*J-1))=DCMPLX(HEST((3*I-1),(3*J-1)),0.D0)
	ZHEST((3*I-1),(3*J))  =DCMPLX(HEST((3*I-1),(3*J)),0.D0)
	ZHEST((3*I),(3*J-2))  =DCMPLX(HEST((3*I),(3*J-2)),0.D0)
	ZHEST((3*I),(3*J-1))  =DCMPLX(HEST((3*I),(3*J-1)),0.D0)
	ZHEST((3*I),(3*J))    =DCMPLX(HEST((3*I),(3*J)),0.D0)
*
	ZGEST((3*I-2),(3*J-2))=DCMPLX(GEST((3*I-2),(3*J-2)),0.D0)
	ZGEST((3*I-2),(3*J-1))=DCMPLX(GEST((3*I-2),(3*J-1)),0.D0)
	ZGEST((3*I-2),(3*J))  =DCMPLX(GEST((3*I-2),(3*J)),0.D0)
	ZGEST((3*I-1),(3*J-2))=DCMPLX(GEST((3*I-1),(3*J-2)),0.D0)
	ZGEST((3*I-1),(3*J-1))=DCMPLX(GEST((3*I-1),(3*J-1)),0.D0)
	ZGEST((3*I-1),(3*J))  =DCMPLX(GEST((3*I-1),(3*J)),0.D0)
	ZGEST((3*I),(3*J-2))  =DCMPLX(GEST((3*I),(3*J-2)),0.D0)
	ZGEST((3*I),(3*J-1))  =DCMPLX(GEST((3*I),(3*J-1)),0.D0)
	ZGEST((3*I),(3*J))    =DCMPLX(GEST((3*I),(3*J)),0.D0)
 5	CONTINUE
*
* ZERANDO AS MATRIZES H E G
*
      DO 50 J=1,N
      DO 50 I=1,N
      II=3*(I-1)
      DO 50 JN=1,3
      DO 50 IN=1,3
	ZGP(II+IN,II+JN)=(0.D0,0.D0)
	ZHP(II+IN,II+JN)=(0.D0,0.D0)
	ZG(II+IN,II+JN)=(0.D0,0.D0)
 50	ZH(II+IN,II+JN)=(0.D0,0.D0)
*
* C�LCULO DOS COEFICIENTES DAS MATRIZES H E G
*
      DO 220 J=1,N
*
* C�LCULO DAS COMPONENTES DO VETOR NORMAL
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
* ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE H E G N�O SINGULAR
*
 120	CALL NONSINGD(ZHELEM,ZGELEM,CO,CXM(I),CYM(I),CZM(I),ETA,ZGE,ZCS,
     $ZCP,DELTA,PI,FR,NPG)
*
 	II=3*(I-1)
      DO 190 JN=1,3
      DO 190 IN=1,3
	ZGP(II+IN,JJ+JN)=ZGELEM(IN,JN)
	ZHP(II+IN,JJ+JN)=ZHELEM(IN,JN)
 190  CONTINUE
*
	GO TO 150
*
* ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE H E G SINGULAR
* ATRAV�S DA DIFEREN�A DIN�MICO - EST�TICO
*
 140	CALL SING_DE (ZHELEM,ZGELEM,CO,CXM(I),CYM(I),CZM(I),ETA,ZGE,ZCS,
     $ZCP,C1,C2,C3,C4,DELTA,PI,FR,NPG)
*
	II=3*(I-1)
      DO 191 JN=1,3
      DO 191 IN=1,3
	ZGP(II+IN,JJ+JN)=ZGELEM(IN,JN)+ZGEST(II+IN,JJ+JN)
	ZHP(II+IN,JJ+JN)=ZHELEM(IN,JN)+ZHEST(II+IN,JJ+JN)
 191  CONTINUE
*
 150	CONTINUE 
*
  200 CONTINUE
  220 CONTINUE
*
* TRANSFORMA��O DAS MATRIZES ZHP E ZGP PROVIS�RIAS
* NAS MATRIZES ZH E ZG FINAIS
*
	NN=3*NBE
	DO 400 IMA=1,NN
	DO 400 IMB=1,NN
	ZH(IMA,IMB)=ZHP(IMA,IMB)
	ZG(IMA,IMB)=ZGP(IMA,IMB)
 400	CONTINUE
*
* REORDENA AS COLUNAS DO SISTEMA DE EQUA��ES DE ACORDO COM AS
*
* CONDI��ES DE CONTORNO E FORMA A MATRIZ A QUE � ARMAZENADA EM G
*
      NN=3*NBE
      DO 250 J=1,NN
      IF(KODE(J)) 250,230,250
 230  DO 240 I=1,NN
      ZCH=ZG(I,J)*ZGE
      ZG(I,J)=-ZH(I,J)
 240  ZH(I,J)=-ZCH
 250  CONTINUE
*
* FORMA O LADO DIREITO DO SISTEMA {VETOR f} QUE � ARMAZENADO EM ZFI
*
      DO 300 I=1,NN
      ZFI(I)=(0.D0,0.D0)
      DO 300 J=1,NN
      ZFI(I)=ZFI(I)+ZG(I,J)*ZDFI(J)
  300 CONTINUE
*
	RETURN
	END