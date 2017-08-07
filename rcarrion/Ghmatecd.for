************************************************************************
*                                                                      *
*         ESTA SUBROTINA CALCULA AS MATRIZES [H] E [G];                *
*                                                                      *
*        AINDA MONTA O SISTEMA DE EQUAÇÕES [A] {x} = {f}               *
*                                                                      *
************************************************************************



      SUBROUTINE GHMATECD (CX,CY,CZ,CXM,CYM,CZM,HEST,GEST,ZH,ZG,ZFI,DFI,
     $    ZDFI,KODE,NE,NX,NCOX,CONE,DELTA,PI,N,NBE,NP,NPG,GE,RNU,RMU,
     $    L,FR,DAM,RHO,ZGE,ZCS,ZCP,C1,C2,C3,C4,ETAS,GI,OME)
        
        USE omp_lib


        IMPLICIT NONE
       
#ifdef USE_GPU
        INCLUDE 'kernels/Ghmatecd_cu.fd'
#endif
        
        INTEGER INP,INQ,IPR,IPS,IPT
        COMMON  INP,INQ,IPR,IPS,IPT
        

        REAL, DIMENSION(NP), INTENT(IN) :: CX, CY, CZ
        REAL, DIMENSION(N), INTENT(IN) :: CXM, CYM, CZM
        REAL, DIMENSION(3*NBE, 3*N), INTENT(IN) :: HEST, GEST

        COMPLEX, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ZH, ZG
        COMPLEX, DIMENSION(:), INTENT(OUT), ALLOCATABLE:: ZFI
        REAL, INTENT(IN) :: DFI(3*NBE)
        COMPLEX, INTENT(OUT), ALLOCATABLE :: ZDFI(:)
        INTEGER, INTENT(IN) :: KODE(3*NBE),NE,NX,NCOX,CONE(N,4)
        REAL, INTENT(IN) :: DELTA(3,3),PI
        INTEGER, INTENT(IN) :: L,N,NBE,NP,NPG
        REAL, INTENT(IN) :: GE,RNU,RMU,FR,DAM,RHO
        COMPLEX,   INTENT(IN) :: ZGE,ZCS,ZCP
        REAL, INTENT(IN) :: C1,C2,C3,C4
        REAL, INTENT(IN) :: ETAS(3,N)

        COMPLEX ZCH
        REAL, INTENT(IN) :: GI(NPG), OME(NPG)
        DOUBLE PRECISION :: t1, t2
        INTEGER NN,I,J, stats1, stats2

#define USE_CPU

#ifdef USE_GPU
#undef USE_CPU
#endif

#ifdef TEST_CUDA
#undef  USE_CPU
#undef  USE_GPU
#define USE_CPU
#define USE_GPU
        COMPLEX, ALLOCATABLE :: ZHP(:,:), ZGP(:,:)
#endif

#ifdef  USE_CPU
        INTEGER N1,N2,N3,N4,II,JJ
        REAL :: CO(4,3)
#endif
#ifdef USE_GPU
        INTEGER RET
#endif

*
* TRANSFORMAÇÃO DAS CONDIÇÕES DE CONTORNO EM NÚMEROS COMPLEXOS
*
        NN = 3*NBE

        ALLOCATE(ZDFI(NN), STAT = stats1)
        ALLOCATE(ZFI(NN) , STAT = stats2)
        IF (stats1 /= 0 .OR. stats2 /= 0) THEN
            PRINT*, "MEMÓRIA INSUFICIENTE"
        ENDIF

        DO I=1,NBE
            ZDFI(3*I-2) = CMPLX(DFI(3*I-2),0.0)
            ZDFI(3*I-1) = CMPLX(DFI(3*I-1),0.0)
            ZDFI(3*I  ) = CMPLX(DFI(3*I),0.0)
        ENDDO
*
* CÁLCULO DOS COEFICIENTES DAS MATRIZES H E G
*
*
        ALLOCATE(ZH(3*NBE, 3*N), STAT = stats1)
        ALLOCATE(ZG(3*NBE, 3*N), STAT = stats2)

        IF (stats1 /= 0 .or. stats2 /= 0) THEN
            PRINT*, "MEMÓRIA INSUFICIENTE"
            STOP
        ENDIF

        t1 = OMP_GET_WTIME()

#ifdef USE_CPU
* ZERANDO AS MATRIZES H E G
*

!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(N1,N2,N3,N4,J,I,CO,II,JJ)
        DO J=1,N
            N1=CONE(J,1)
            N2=CONE(J,2)
            N3=CONE(J,3)
            N4=CONE(J,4)
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
C                   ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE H E G SINGULAR
C                   ATRAVÉS DA DIFERENÇA DINÂMICO - ESTÁTICO

                    CALL SING_DE (ZH(II:II+2, JJ:JJ+2),
     $                  ZG(II:II+2, JJ:JJ+2),
     $                  CO,CXM(I),CYM(I),CZM(I),
     $                  ETAS(1:3,J),
     $                  ZGE,ZCS,ZCP,C1,C2,C3,C4,DELTA,PI,FR,GI,OME,NPG)


                    ZG(II:II+2, JJ:JJ+2) = ZG(II:II+2, JJ:JJ+2) +
     $                  GEST(II:II+2, JJ:JJ+2)
                    ZH(II:II+2, JJ:JJ+2) = ZH(II:II+2, JJ:JJ+2) + 
     $                  HEST(II:II+2, JJ:JJ+2)

                ELSE
*                   ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE H E G NÃO SINGULAR
*
                    CALL NONSINGD (ZH(II:II+2, JJ:JJ+2),
     $                  ZG(II:II+2, JJ:JJ+2),
     $                  CO,CXM(I),CYM(I),CZM(I),
     $                  ETAS(1:3,J),
     $                  ZGE,ZCS,ZCP,DELTA,PI,FR,GI,OME,NPG)
                ENDIF
            ENDDO
        ENDDO
!$OMP END PARALLEL DO
        t2 = OMP_GET_WTIME()
        PRINT *, "GHMATECD: Tempo na CPU: ", (t2-t1)
#endif

#ifdef TEST_CUDA
!FAÇA UMA CÓPIA DAS MATRIZES ZG E ZH PARA COMPARAÇÃO COM O RESULTADO DA GPU.
        ALLOCATE(ZHP(3*NBE, 3*N))
        ALLOCATE(ZGP(3*NBE, 3*N))
        ZGP = ZG
        ZHP = ZH
#endif

#ifdef USE_GPU
        t1 = OMP_GET_WTIME()

        CALL cuda_ghmatecd(
     $                      NBE,
     $                      NPG,
     $                      N,
     $                      NP,
     $                      ZGE,
     $                      ZCS,
     $                      ZCP,
     $                      C1,
     $                      C2,
     $                      C3,
     $                      C4,
     $                      FR,
     $                      HEST,
     $                      GEST,
     $                      ZG,
     $                      ZH,
     $                      RET
     $                      )

        
!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(N1,N2,N3,N4,J,I,CO,II,JJ)
        DO J=1,NBE
            N1=CONE(J,1)
            N2=CONE(J,2)
            N3=CONE(J,3)
            N4=CONE(J,4)
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
            I = J
            II=3*(I-1) + 1


            CALL SING_DE (ZH(II:II+2, JJ:JJ+2),
     $          ZG(II:II+2, JJ:JJ+2),
     $          CO,CXM(I),CYM(I),CZM(I),
     $          ETAS(1:3,J),
     $          ZGE,ZCS,ZCP,C1,C2,C3,C4,DELTA,PI,FR,GI,OME,NPG)


            ZG(II:II+2, JJ:JJ+2) = ZG(II:II+2, JJ:JJ+2) +
     $          GEST(II:II+2, JJ:JJ+2)
            ZH(II:II+2, JJ:JJ+2) = ZH(II:II+2, JJ:JJ+2) + 
     $          HEST(II:II+2, JJ:JJ+2)

        ENDDO
!$OMP END PARALLEL DO

        IF (RET /= 0) THEN
            PRINT*, "GHMATECD: Erro: Matriz Singular."
        ENDIF

        t2 = OMP_GET_WTIME()
    
        PRINT *, "GHMATECD: Tempo na GPU: ", (t2 - t1)
#endif
        
#ifdef TEST_CUDA
        CALL ASSERT_GHMATECD_ZH_ZG(ZH, ZHP, ZG, ZGP, NBE, N)
        DEALLOCATE(ZHP)
        DEALLOCATE(ZGP)
#endif
*
* REORDENA AS COLUNAS DO SISTEMA DE EQUAÇÕES DE ACORDO COM AS
*
* CONDIÇÕES DE CONTORNO E FORMA A MATRIZ A QUE É ARMAZENADA EM G
*
        t1 = OMP_GET_WTIME()
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
        
!        ZFI = MATMUL(ZG(1:NN, 1:NN), ZDFI)
!        PRINT*, NN, 3*N
        CALL CGEMV('N', NN, NN, (1.0,0), ZG, NN, ZDFI, 1, (0,0), ZFI, 1)
        t2 = OMP_GET_WTIME()
        PRINT*, "GHMATECD: Tempo gasto no restante: ", (t2-t1)


C        t1 = OMP_GET_WTIME()
C        PRINT *, "Tempo gasto em Ghmatecd: ", (t1-t0)
        RETURN
      END


      SUBROUTINE ASSERT_GHMATECD_ZH_ZG(ZH, ZHP, ZG, ZGP, NBE, N)
        IMPLICIT NONE
        COMPLEX, INTENT(IN), DIMENSION(3*NBE,3*N) :: ZH,ZHP,ZG,ZGP
        INTEGER, INTENT(IN) :: NBE, N  

        INTEGER :: i, j
        INTEGER :: N3, NBE3
        LOGICAL :: ghmatecd_asserted = .TRUE.
        REAL :: sum_norms = 0, eps
        REAL :: local_sum = 0, max_local_sum = 0

        eps = 1.5E-6*nbe
        N3   = N*3
        NBE3 = NBE*3 

        sum_norms = 0.0
        

        DO j = 1, N3
            local_sum = 0
            DO i = 1, NBE3
                local_sum = local_sum + ABS(ZHP(i,j) - ZH(i,j))
            ENDDO
            sum_norms = sum_norms + local_sum
            max_local_sum = MAX(local_sum, max_local_sum)
        ENDDO

!        DO j = 1, N3
!            DO i = 1, NBE3
!                local_sum = ABS(ZHP(i,j) - ZH(i,j))
!                if (local_sum > max_local_sum) then
!                    max_local_sum = local_sum
!                endif
!            ENDDO
!        ENDDO
        PRINT*, max_local_sum
       

        IF (max_local_sum > eps) THEN
            ghmatecd_asserted = .FALSE.
            PRINT*, "||ZH||_inf = ", max_local_sum
        ENDIF

        sum_norms = 0.0
        max_local_sum = 0
        
        DO j = 1, N3
            local_sum = 0
            DO i = 1, NBE3
                local_sum = local_sum + ABS(ZGP(i,j) - ZG(i,j))
            ENDDO
            sum_norms = sum_norms + local_sum
            max_local_sum = MAX(local_sum, max_local_sum)
        ENDDO

        IF (max_local_sum > eps) THEN
            ghmatecd_asserted = .FALSE.
            PRINT*, "||ZG||_inf = ", max_local_sum
        ENDIF

 200    FORMAT (A,ES7.1)     
       WRITE(0,"(A)") "As matrizes ZH e ZG em Ghmatecd_cu sao iguais as"
        WRITE(0,200) "calculadas em Ghmatecd com um erro de ", eps
        IF (ghmatecd_asserted .EQV. .TRUE.) THEN
            WRITE(0,"(A)")"[OK]"
        ELSE
            WRITE(0,"(A)")"[FALHOU]"
        ENDIF
        WRITE(0,"(A)") ""

      END SUBROUTINE ASSERT_GHMATECD_ZH_ZG
