************************************************************************
*                                                                      *
*         ESTA SUBROTINA CALCULA AS MATRIZES [H] E [G];                *
*                                                                      *
*        AINDA MONTA O SISTEMA DE EQUAÇÕES [A] {x} = {f}               *
*                                                                      *
************************************************************************

      SUBROUTINE GHMATECD (CX,CY,CZ,CXM,CYM,CZM,HEST,GEST,ZH,ZG,ZFI,DFI,
     $    ZDFI,KODE,NE,NX,NCOX,CONE,DELTA,PI,N,NBE,NP,NPG,GE,RNU,RMU,
     $    L,FR,DAM,RHO,ZGE,ZCS,ZCP,C1,C2,C3,C4,ETAS)
        
        USE omp_lib

        IMPLICIT NONE
        INTEGER INP,INQ,IPR,IPS,IPT
        COMMON  INP,INQ,IPR,IPS,IPT
        

        DOUBLE PRECISION, INTENT(IN) :: CX(NCOX),CY(NCOX),CZ(NCOX)
        DOUBLE PRECISION, INTENT(IN) :: CXM(NE),CYM(NE),CZM(NE)
        DOUBLE PRECISION, INTENT(IN) :: HEST(NX,NX), GEST(NX,NX) 
        DOUBLE COMPLEX,   INTENT(OUT):: ZH(NX,NX), ZG(NX,NX)
        DOUBLE COMPLEX,   INTENT(OUT):: ZFI(NX)
        DOUBLE PRECISION, INTENT(IN) :: DFI(NX)
        DOUBLE COMPLEX,   INTENT(OUT):: ZDFI(NX)
        INTEGER,          INTENT(IN) :: KODE(NX),NE,NX,NCOX,CONE(NE,4)
        DOUBLE PRECISION, INTENT(IN) :: DELTA(3,3),PI
        INTEGER,          INTENT(IN) :: N,NBE,NP,NPG
        DOUBLE PRECISION, INTENT(IN) :: GE,RNU,RMU,L,FR,DAM,RHO
        DOUBLE COMPLEX,   INTENT(IN) :: ZGE,ZCS,ZCP
        DOUBLE PRECISION, INTENT(IN) :: C1,C2,C3,C4
        DOUBLE PRECISION, INTENT(IN) :: ETAS(3,NX)

        DOUBLE COMPLEX ZCH!,ZHEST(NX,NX),ZGEST(NX,NX)
        DOUBLE PRECISION CO(4,3)!, ETA(3), A,B,C,R
        INTEGER N1,N2,N3,N4,NN,I,J,II,JJ

#ifdef TEST_GHMATECD_CUDA
        DOUBLE COMPLEX ZHP(NX,NX), ZGP(NX,NX)
#endif

C        DOUBLE PRECISION t0, t1
C        t0 = OMP_GET_WTIME()
*
* TRANSFORMAÇÃO DAS CONDIÇÕES DE CONTORNO EM NÚMEROS COMPLEXOS
*
        DO I=1,NBE
            ZDFI(3*I-2)=DCMPLX(DFI(3*I-2),0.D0)
            ZDFI(3*I-1)=DCMPLX(DFI(3*I-1),0.D0)
            ZDFI(3*I)=DCMPLX(DFI(3*I),0.D0)
        ENDDO
*
* TRANSFORMAÇÃO DAS MATRIZES [HEST] E [GEST] EM NÚMEROS COMPLEXOS
*

!        DO J=1,N
!            DO I=1,N
!               ZHEST((3*I-2),(3*J-2))=DCMPLX(HEST((3*I-2),(3*J-2)),0.D0)
!               ZHEST((3*I-2),(3*J-1))=DCMPLX(HEST((3*I-2),(3*J-1)),0.D0)
!               ZHEST((3*I-2),(3*J))  =DCMPLX(HEST((3*I-2),(3*J)),0.D0)
!               ZHEST((3*I-1),(3*J-2))=DCMPLX(HEST((3*I-1),(3*J-2)),0.D0)
!               ZHEST((3*I-1),(3*J-1))=DCMPLX(HEST((3*I-1),(3*J-1)),0.D0)
!               ZHEST((3*I-1),(3*J))  =DCMPLX(HEST((3*I-1),(3*J)),0.D0)
!               ZHEST((3*I),(3*J-2))  =DCMPLX(HEST((3*I),(3*J-2)),0.D0)
!               ZHEST((3*I),(3*J-1))  =DCMPLX(HEST((3*I),(3*J-1)),0.D0)
!               ZHEST((3*I),(3*J))    =DCMPLX(HEST((3*I),(3*J)),0.D0)
!
!               ZGEST((3*I-2),(3*J-2))=DCMPLX(GEST((3*I-2),(3*J-2)),0.D0)
!               ZGEST((3*I-2),(3*J-1))=DCMPLX(GEST((3*I-2),(3*J-1)),0.D0)
!               ZGEST((3*I-2),(3*J))  =DCMPLX(GEST((3*I-2),(3*J)),0.D0)
!               ZGEST((3*I-1),(3*J-2))=DCMPLX(GEST((3*I-1),(3*J-2)),0.D0)
!               ZGEST((3*I-1),(3*J-1))=DCMPLX(GEST((3*I-1),(3*J-1)),0.D0)
!               ZGEST((3*I-1),(3*J))  =DCMPLX(GEST((3*I-1),(3*J)),0.D0)
!               ZGEST((3*I),(3*J-2))  =DCMPLX(GEST((3*I),(3*J-2)),0.D0)
!               ZGEST((3*I),(3*J-1))  =DCMPLX(GEST((3*I),(3*J-1)),0.D0)
!               ZGEST((3*I),(3*J))    =DCMPLX(GEST((3*I),(3*J)),0.D0)
!            ENDDO
!        ENDDO
*
* ZERANDO AS MATRIZES H E G
*

        ZG(1:3*NBE, 1:3*NBE) = (0.D0, 0.D0)
        ZH(1:3*NBE, 1:3*NBE) = (0.D0, 0.D0)
*
* CÁLCULO DOS COEFICIENTES DAS MATRIZES H E G
*
*
* CÁLCULO DAS COMPONENTES DO VETOR NORMAL
* USANDO O PRODUTO VETORIAL DOS LADOS 1-2 E 1-3
*

!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(N1,N2,N3,N4,J,I,CO,II,JJ)
        DO J=1,N
            N1=CONE(J,1)
            N2=CONE(J,2)
            N3=CONE(J,3)
            N4=CONE(J,4)
C            A=(CY(N2)-CY(N1))*(CZ(N3)-CZ(N1)) - 
C     $           (CZ(N2)-CZ(N1))*(CY(N3)-CY(N1))
C            B=(CZ(N2)-CZ(N1))*(CX(N3)-CX(N1)) - 
C     $           (CX(N2)-CX(N1))*(CZ(N3)-CZ(N1))
C            C=(CX(N2)-CX(N1))*(CY(N3)-CY(N1)) - 
C     $           (CY(N2)-CY(N1))*(CX(N3)-CX(N1))
C            R=DSQRT(A*A+B*B+C*C)
C            ETA(1)=A/R
C            ETA(2)=B/R
C            ETA(3)=C/R
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
            DO I=1,NBE
                II=3*(I-1) + 1

                IF (I == J) THEN
*                   ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE H E G SINGULAR
*                   ATRAVÉS DA DIFERENÇA DINÂMICO - ESTÁTICO
*
                    CALL SING_DE (ZH(II:II+2, JJ:JJ+2),
     $                  ZG(II:II+2, JJ:JJ+2),
     $                  CO,CXM(I),CYM(I),CZM(I),
     $                  ETAS(1:3,J),
     $                  ZGE,ZCS,ZCP,C1,C2,C3,C4,DELTA,PI,FR,NPG)


                    ZG(II:II+2, JJ:JJ+2) = ZG(II:II+2, JJ:JJ+2) +
     $                  GEST(II:II+2, JJ:JJ+2)
                    ZH(II:II+2, JJ:JJ+2) = ZH(II:II+2, JJ:JJ+2) + 
     $                  HEST(II:II+2, JJ:JJ+2)

                ELSE
*                   ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE H E G NÃO SINGULAR
*
                    CALL NONSINGD(ZH(II:II+2, JJ:JJ+2),
     $                  ZG(II:II+2, JJ:JJ+2),
     $                  CO,CXM(I),CYM(I),CZM(I),
     $                  ETAS(1:3,J),ZGE,ZCS,ZCP,DELTA,PI,FR,NPG)
                ENDIF
            ENDDO
        ENDDO
!$OMP END PARALLEL DO
*
* TRANSFORMAÇÃO DAS MATRIZES ZHP E ZGP PROVISÓRIAS
* NAS MATRIZES ZH E ZG FINAIS
*
        NN=3*NBE

! MATRIZES PROVISÓRIAS ELIMINADAS. 
!        ZH(1:NN, 1:NN) = ZHP(1:NN, 1:NN)
!        ZG(1:NN, 1:NN) = ZGP(1:NN, 1:NN)

#ifdef TEST_GHMATECD_CUDA
        ZHP = (0.D0,0.D0)
        ZGP = (0.D0,0.D0)

        CALL cuda_ghmatecd(NE,
     $                      NBE,
     $                      NX,
     $                      NPG,
     $                      NCOX,
     $                      N,
     $                      CONE,
     $                      CX,
     $                      CY,
     $                      CZ,
     $                      CXM,
     $                      CYM,
     $                      CZM,
     $                      ETAS,
     $                      ZGE,
     $                      ZCS,
     $                      ZCP,
     $                      C1,
     $                      C2,
     $                      C3,
     $                      C4,
     $                      DELTA,
     $                      PI,
     $                      FR,
     $                      HEST,
     $                      GEST,
     $                      ZGP,
     $                      ZHP
     $                      )
        CALL ASSERT_GHMATECD_ZH_ZG(ZH, ZHP, ZG, ZGP, NX, NN)
#endif
*
* REORDENA AS COLUNAS DO SISTEMA DE EQUAÇÕES DE ACORDO COM AS
*
* CONDIÇÕES DE CONTORNO E FORMA A MATRIZ A QUE É ARMAZENADA EM G
*
        DO  J=1,NN
            IF (KODE(J) == 0) THEN
                DO I=1,NN
                    ZCH = ZG(I,J)*ZGE
                    ZG(I,J) = -ZH(I,J)
                    ZH(I,J) = -ZCH
                ENDDO
            ENDIF
        ENDDO

*
* FORMA O LADO DIREITO DO SISTEMA {VETOR f} QUE É ARMAZENADO EM ZFI
*
        
        ZFI(1:NN) = MATMUL(ZG(1:NN, 1:NN), ZDFI(1:NN))
C        t1 = OMP_GET_WTIME()
C        PRINT *, "Tempo gasto em Ghmatecd: ", (t1-t0)
        RETURN
      END

      SUBROUTINE ASSERT_GHMATECD_ZH_ZG(ZH, ZHP, ZG, ZGP, NX, NN)
        DOUBLE COMPLEX, INTENT(IN), DIMENSION(NX,NX) :: ZH,ZHP,ZG,ZGP
        INTEGER, INTENT(IN) :: NX, NN

        INTEGER :: i, j
        LOGICAL :: ghmatecd_asserted = .TRUE.

        DO j = 1, NN
            DO i = 1, NN
                IF (ZHP(i,j) /= ZH(i,j)) THEN
                    ghmatecd_asserted = .FALSE.
                ENDIF
            ENDDO
        ENDDO
       WRITE(0,"(A)") "As matrizes ZH e ZG em Ghmatecd_cu sao iguais as"
        WRITE(0,"(A)") "calculadas em Ghmatecd?"
        IF (ghmatecd_asserted .EQV. .TRUE.) THEN
            WRITE(0,"(A)")"[OK]"
        ELSE
            WRITE(0,"(A)"),"[FALHOU]"
        ENDIF
        WRITE(0,"(A)"), ""

      END SUBROUTINE ASSERT_GHMATECD_ZH_ZG
