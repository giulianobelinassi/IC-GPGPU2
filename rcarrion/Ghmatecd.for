************************************************************************
*                                                                      *
*         ESTA SUBROTINA CALCULA AS MATRIZES [H] E [G];                *
*                                                                      *
*        AINDA MONTA O SISTEMA DE EQUAÇÕES [A] {x} = {f}               *
*                                                                      *
************************************************************************
*
      SUBROUTINE GHMATECD (CX,CY,CZ,CXM,CYM,CZM,HEST,GEST,ZH,ZG,ZFI,DFI,
     $    ZDFI,KODE,NE,NX,NCOX,CONE,DELTA,PI,N,NBE,NP,NPG,GE,RNU,RMU,
     $    L,FR,DAM,RHO,ZGE,ZCS,ZCP,C1,C2,C3,C4)
        IMPLICIT REAL*8 (A-H,O-Y)
        IMPLICIT COMPLEX*16 (Z)
        COMMON INP,INQ,IPR,IPS,IPT
        DIMENSION CX(NCOX),CY(NCOX),CZ(NCOX)
        DIMENSION CXM(NE),CYM(NE),CZM(NE)
        DIMENSION HEST(NX,NX),GEST(NX,NX)
        DIMENSION ZH(NX,NX),ZG(NX,NX)
        DIMENSION ZHELEM(3,3),ZGELEM(3,3)
C
C quando coloco NX o programa não compila
C No gfortran compila :-)
        DIMENSION ZHP(NX,NX),ZGP(NX,NX)
        DIMENSION ZHEST(NX,NX),ZGEST(NX,NX)
C
        DIMENSION CO(4,3),ETA(3)
        DIMENSION ZFI(NX),ZDFI(NX)
        DIMENSION DFI(NX)
        DIMENSION DELTA(3,3)
        INTEGER CONE(NE,4),KODE(NX)
C       REAL t0, t1

*
* TRANSFORMAÇÃO DAS CONDIÇÕES DE CONTORNO EM NÚMEROS COMPLEXOS
*
        DO 2 I=1,NBE
            ZDFI(3*I-2)=DCMPLX(DFI(3*I-2),0.D0)
            ZDFI(3*I-1)=DCMPLX(DFI(3*I-1),0.D0)
            ZDFI(3*I)=DCMPLX(DFI(3*I),0.D0)
 2      CONTINUE
*
* TRANSFORMAÇÃO DAS MATRIZES [HEST] E [GEST] EM NÚMEROS COMPLEXOS
*
        DO 5 J=1,N
            DO 5 I=1,N
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
 5      CONTINUE
*
* ZERANDO AS MATRIZES H E G
*
        DO 50 I=1,N
            II=3*(I-1) + 1 
C           Use a notação de particionamento para o compilador decidir como percorrer a matriz.
            ZGP(II:II+2, II:II+2) = (0.D0,0.D0)
            ZHP(II:II+2, II:II+2) = (0.D0,0.D0)
            ZG (II:II+2, II:II+2) = (0.D0,0.D0)
            ZH (II:II+2, II:II+2) = (0.D0,0.D0)
 50     CONTINUE
*
* CÁLCULO DOS COEFICIENTES DAS MATRIZES H E G
*
*
* CÁLCULO DAS COMPONENTES DO VETOR NORMAL
* USANDO O PRODUTO VETORIAL DOS LADOS 1-2 E 1-3
*

        DO 220 J=1,N
            N1=CONE(J,1)
            N2=CONE(J,2)
            N3=CONE(J,3)
            N4=CONE(J,4)
            A=(CY(N2)-CY(N1))*(CZ(N3)-CZ(N1)) - 
     $           (CZ(N2)-CZ(N1))*(CY(N3)-CY(N1))
            B=(CZ(N2)-CZ(N1))*(CX(N3)-CX(N1)) - 
     $           (CX(N2)-CX(N1))*(CZ(N3)-CZ(N1))
            C=(CX(N2)-CX(N1))*(CY(N3)-CY(N1)) - 
     $           (CY(N2)-CY(N1))*(CX(N3)-CX(N1))
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

            JJ=3*(J-1) + 1
            DO 200 I=1,NBE
                II=3*(I-1) + 1

                IF (I == J) THEN
*                   ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE H E G SINGULAR
*                   ATRAVÉS DA DIFERENÇA DINÂMICO - ESTÁTICO
*
                    CALL SING_DE (ZHELEM,ZGELEM,CO,CXM(I),CYM(I),CZM(I),
     $                  ETA,ZGE,ZCS,ZCP,C1,C2,C3,C4,DELTA,PI,FR,NPG)

                    ZGP(II:II+2, JJ:JJ+2)=ZGELEM+ZGEST(II:II+2, JJ:JJ+2)
                    ZHP(II:II+2, JJ:JJ+2)=ZHELEM+ZHEST(II:II+2, JJ:JJ+2)
                ELSE
*                   ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE H E G NÃO SINGULAR
*
                    CALL NONSINGD(ZHELEM,ZGELEM,CO,CXM(I),CYM(I),CZM(I),
     $                  ETA,ZGE,ZCS,ZCP,DELTA,PI,FR,NPG)

                    ZGP(II:II+2, JJ:JJ+2) = ZGELEM
                    ZHP(II:II+2, JJ:JJ+2) = ZHELEM
                ENDIF
 200        CONTINUE
 220    CONTINUE
*
* TRANSFORMAÇÃO DAS MATRIZES ZHP E ZGP PROVISÓRIAS
* NAS MATRIZES ZH E ZG FINAIS
*
        NN=3*NBE
        ZH(1:NN, 1:NN) = ZHP(1:NN, 1:NN)
        ZG(1:NN, 1:NN) = ZGP(1:NN, 1:NN)
*
* REORDENA AS COLUNAS DO SISTEMA DE EQUAÇÕES DE ACORDO COM AS
*
* CONDIÇÕES DE CONTORNO E FORMA A MATRIZ A QUE É ARMAZENADA EM G
*
        DO 250 J=1,NN
            IF (KODE(J) == 0) THEN
                DO 240 I=1,NN
                    ZCH = ZG(I,J)*ZGE
                    ZG(I,J) = -ZH(I,J)
 240                ZH(I,J) = -ZCH

            ENDIF
 250    CONTINUE

*
* FORMA O LADO DIREITO DO SISTEMA {VETOR f} QUE É ARMAZENADO EM ZFI
*
        
        ZFI(1:NN) = MATMUL(ZG(1:NN, 1:NN), ZDFI(1:NN))
        RETURN
      END
