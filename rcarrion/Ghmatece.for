************************************************************************
*                                                                      *
*         ESTA SUBROTINA CALCULA AS MATRIZES [H] E [G];                *
*                                                                      *
*        AINDA MONTA O SISTEMA DE EQUAÇÕES [A] {x} = {f}               *
*                                                                      *
************************************************************************
*
      SUBROUTINE GHMATECE (CX,CY,CZ,CXM,CYM,CZM,HEST,GEST,NE,NX,NCOX,
     $    CONE,N,NBE,NP,NPG,GE,RNU,RMU,DELTA,PI,C1,C2,C3,C4,ETAS)
*
        USE omp_lib        

        IMPLICIT REAL (A-H,O-Y)
        IMPLICIT COMPLEX (Z)
        COMMON INP,INQ,IPR,IPS,IPT
        REAL, DIMENSION(NP), INTENT(IN) :: CX, CY, CZ
        REAL, DIMENSION(N) , INTENT(IN) :: CXM, CYM, CZM
        INTEGER, DIMENSION(N, 4), INTENT(IN) :: CONE
        REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: HEST, GEST
        REAL, DIMENSION(NPG) :: GI, OME 

        DIMENSION HELEM(3,3),GELEM(3,3)
        DIMENSION DELTA(3,3)
        DIMENSION CO(4,3)
        REAL, INTENT(IN) :: ETAS(3,n)
        INTEGER stats1, stats2

        DOUBLE PRECISION :: t1, t2

#ifdef TEST_GHMATECE_CUDA
#undef  GHMATECE_USE_CPU
#undef  GHMATECE_USE_GPU
#define GHMATECE_USE_CPU
#define GHMATECE_USE_GPU
        REAL, ALLOCATABLE, DIMENSION(:,:) :: HP, GP
#endif

        PI=4.0*ATAN(1.0)
* ACIONA ROTINA QUE CALCULA OS PONTOS E PESOS DE GAUSS
*
        CALL GAULEG(-1.0,1.0,GI,OME,NPG)
*
* CONSTANTES USADAS NAS SOLUÇÕES FUNDAMENTAIS ESTÁTICAS
*
        C1=1.0/(16.0*PI*RMU*(1.0-RNU))
        C2=3.0-(4.0*RNU)
        C3=-1.0/(8.0*PI*(1.0-RNU))
        C4=1.0-(2.0*RNU)
*
        DO 25 J=1,3
            DO 25 I=1,3
                IF (I == J) THEN
                    DELTA(I, J) = 1.0
                ELSE
                    DELTA(I, J) = 0.0
                ENDIF
  25    CONTINUE
*
* ZERANDO AS MATRIZES H E G
*
        ALLOCATE(HEST(3*NBE, 3*N), STAT = stats1)        
        ALLOCATE(GEST(3*NBE, 3*N), STAT = stats2)
        
        IF (stats1 /= 0 .or. stats2 /= 0) THEN
            PRINT*, "MEMÓRIA INSUFICIENTE!"
            STOP
        ENDIF

*
* CÁLCULO DOS COEFICIENTES DAS MATRIZES H E G
*

#ifdef GHMATECE_USE_CPU
        t1 = OMP_GET_WTIME() 
        GEST = 0
        HEST = 0
!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(N1,N2,N3,N4,J,I,CO,II,JJ,HELEM,GELEM)
        DO J=1,N
*
* CÁLCULO DAS COMPONENTES DO VETOR NORMAL
*
            N1=CONE(J,1)
            N2=CONE(J,2)
            N3=CONE(J,3)
            N4=CONE(J,4)
C            A=(CY(N2)-CY(N1))*(CZ(N3)-CZ(N1)) - 
C     $          (CZ(N2)-CZ(N1))*(CY(N3)-CY(N1))
C            B=(CZ(N2)-CZ(N1))*(CX(N3)-CX(N1)) -
C     $          (CX(N2)-CX(N1))*(CZ(N3)-CZ(N1))
C            C=(CX(N2)-CX(N1))*(CY(N3)-CY(N1)) - 
C     $          (CY(N2)-CY(N1))*(CX(N3)-CX(N1))
C            R=DSQRT(A*A+B*B+C*C)
C            ETAS(1)=A/R
C            ETAS(2)=B/R
C            ETAS(3)=C/R
            
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
*	 
                IF (I == J) THEN
!                   ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE G SINGULAR
                    CALL SINGGE(HELEM,GELEM,CXM(I),CYM(I),CZM(I),
     $                      ETAS(1:3,J),CX,
     $                  CY,CZ,N1,N2,N3,N4,NCOX,N,NP,NPG,GE,RNU,RMU,C1,
     $                  C2,C3,C4,DELTA,GI,OME)
                ELSE
*                   ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE H E G NÃO SINGULAR
                    
                    CALL NONSINGE(HELEM,GELEM,CO,CXM(I),CYM(I),CZM(I),
     $                  ETAS(1:3,J),N,NP,NPG,GE,RNU,RMU,C1,C2,C3,C4,
     $                  DELTA,GI,OME)

                ENDIF

                II=3*(I-1) + 1
                GEST(II:II+2, JJ:JJ+2) = GELEM
                HEST(II:II+2, JJ:JJ+2) = HELEM
*
            ENDDO
        ENDDO
!$OMP END PARALLEL DO
        t2 = OMP_GET_WTIME()
        PRINT *, "GHMATECE: Tempo na CPU: ", (t2-t1)
#endif

#ifdef TEST_GHMATECE_CUDA
        ALLOCATE(HP(3*NBE, 3*N))
        ALLOCATE(GP(3*NBE, 3*N))
        GP = GEST
        HP = HEST
#endif

#ifdef GHMATECE_USE_GPU
        t1 = OMP_GET_WTIME()
       
        CALL cuda_ghmatece(
     $          NBE,
     $          NPG,
     $          N,
     $          NP,
     $          CONE,
     $          CX,
     $          CY,
     $          CZ,
     $          CXM,
     $          CYM,
     $          CZM,
     $          ETAS,
     $          C1,
     $          C2,
     $          C3,
     $          C4,
     $          FR,
     $          HEST,
     $          GEST,
     $          OME,
     $          GI,
     $          STATS1
     $          )

        CALL GHMATECE_SINGULAR(
     $          HEST,
     $          GEST,
     $          CXM,
     $          CYM, 
     $          CZM,
     $          ETAS,
     $          CX,
     $          CY,
     $          CZ,
     $          NCOX,
     $          N,
     $          NP,
     $          NPG,
     $          GE,
     $          RNU,
     $          RMU,
     $          C1,
     $          C2,
     $          C3,
     $          C4,
     $          DELTA,
     $          GI,
     $          OME,
     $          CONE,
     $          NBE
     $          )
        t2 = OMP_GET_WTIME()
        PRINT *, "GHMATECE: Tempo na GPU: ", (t2-t1)
#endif


#ifdef TEST_GHMATECE_CUDA
        CALL ASSERT_GHMATECE_H_G(HEST, HP, GEST, GP, NBE, N)
        DEALLOCATE(HP)
        DEALLOCATE(GP)
#endif

*
* ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE H SINGULAR ATRAVÉS
* DA CONSIDERAÇÃO DO MOVIMENTO DO CORPO RÍGIDO
*
        DO 224 MA=1,NBE
            HEST(3*MA-2,3*MA-2)=0.0
            HEST(3*MA-2,3*MA-1)=0.0
            HEST(3*MA-2,3*MA)  =0.0
            HEST(3*MA-1,3*MA-2)=0.0
            HEST(3*MA-1,3*MA-1)=0.0
            HEST(3*MA-1,3*MA)  =0.0
            HEST(3*MA,3*MA-2)  =0.0
            HEST(3*MA,3*MA-1)  =0.0
            HEST(3*MA,3*MA)    =0.0
            DO 227 MB=1, N
                IF (MA /= MB) THEN
                    HEST(3*MA-2,3*MA-2)=HEST(3*MA-2,3*MA-2) - 
     $                  HEST(3*MA-2,3*MB-2)
                    HEST(3*MA-2,3*MA-1)=HEST(3*MA-2,3*MA-1) - 
     $                  HEST(3*MA-2,3*MB-1)
                    HEST(3*MA-2,3*MA)  =HEST(3*MA-2,3*MA) - 
     $                  HEST(3*MA-2,3*MB)
                    HEST(3*MA-1,3*MA-2)=HEST(3*MA-1,3*MA-2) - 
     $                  HEST(3*MA-1,3*MB-2)
                    HEST(3*MA-1,3*MA-1)=HEST(3*MA-1,3*MA-1) - 
     $                  HEST(3*MA-1,3*MB-1)
                    HEST(3*MA-1,3*MA)  =HEST(3*MA-1,3*MA) - 
     $                  HEST(3*MA-1,3*MB)
                    HEST(3*MA,3*MA-2)  =HEST(3*MA,3*MA-2) - 
     $                  HEST(3*MA,3*MB-2)
                    HEST(3*MA,3*MA-1)  =HEST(3*MA,3*MA-1) - 
     $                  HEST(3*MA,3*MB-1)
                    HEST(3*MA,3*MA)    =HEST(3*MA,3*MA) - 
     $                  HEST(3*MA,3*MB)
                ENDIF
            
 227        CONTINUE
 224    CONTINUE
*
        RETURN
      END

!     Note que há cálculos muito diferentes dos demais no problema
!     singular, e portanto decidi tratá-lo de maneira diferenciada.
      SUBROUTINE GHMATECE_SINGULAR(HEST, GEST, CXM, CYM, CZM, ETAS,
     $      CX, CY, CZ, NCOX, N, NP, NPG, GE, RNU, RMU, 
     $      C1, C2, C3, C4, DELTA, GI, OME, CONE, NBE)
       
        IMPLICIT NONE
        REAL, DIMENSION(NP), INTENT(IN) :: CX, CY, CZ
        REAL, DIMENSION(N) , INTENT(IN) :: CXM, CYM, CZM
        INTEGER, DIMENSION(N, 4), INTENT(IN) :: CONE
        REAL, DIMENSION(3*NBE, 3*N), INTENT(OUT) :: HEST, GEST
        REAL, DIMENSION(NPG) :: GI, OME
        REAL, INTENT(IN) :: ETAS(3,n)
        INTEGER, INTENT(IN) :: NBE, N, NP, NPG, NCOX
        REAL, INTENT(IN) :: GE, RNU, RMU, C1,C2,C3,C4
        REAL, INTENT(IN) :: DELTA(3,3)
        INTEGER :: J, JJ, N1, N2, N3, N4

        DO J=1, NBE  
            JJ = 3*(J-1) + 1
            HEST(JJ:JJ+2, JJ:JJ+2) = 0
            GEST(JJ:JJ+2, JJ:JJ+2) = 0
        ENDDO

!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(N1,N2,N3,N4,J,JJ)
        DO J=1,NBE
            N1=CONE(J,1)
            N2=CONE(J,2)
            N3=CONE(J,3)
            N4=CONE(J,4)
            
            JJ=3*(J-1) + 1
            CALL SINGGE(HEST(JJ:JJ+2, JJ:JJ+2),GEST(JJ:JJ+2, JJ:JJ+2),
     $                  CXM(J),CYM(J),CZM(J),
     $                  ETAS(1:3,J),CX,
     $                  CY,CZ,N1,N2,N3,N4,NCOX,N,NP,NPG,GE,RNU,RMU,C1,
     $                  C2,C3,C4,DELTA,GI,OME)
        ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE GHMATECE_SINGULAR

      SUBROUTINE ASSERT_GHMATECE_H_G(H, HP, G, GP, NBE, N)
        IMPLICIT NONE
        REAL, INTENT(IN), DIMENSION(3*NBE,3*N) :: H,HP,G,GP
        INTEGER, INTENT(IN) :: NBE, N  

        INTEGER :: i, j
        INTEGER :: N3, NBE3
        LOGICAL :: ghmatecd_asserted = .TRUE.
        REAL :: sum_norms = 0, eps
        REAL :: local_sum = 0, max_local_sum = 0

        eps = 0

        N3   = N*3
        NBE3 = NBE*3 

        sum_norms = 0.0

        DO j = 1, N3
            local_sum = 0
            DO i = 1, NBE3
                local_sum = local_sum + ABS(HP(i,j) - H(i,j))
            ENDDO
            sum_norms = sum_norms + local_sum
            max_local_sum = MAX(local_sum, max_local_sum)
        ENDDO

        IF (max_local_sum > eps) THEN
            ghmatecd_asserted = .FALSE.
            PRINT*, "||H||_inf = ", max_local_sum
        ENDIF

        sum_norms = 0.0
        max_local_sum = 0
        
        DO j = 1, N3
            local_sum = 0
            DO i = 1, NBE3
                local_sum=local_sum + ABS(GP(i,j) - G(i,j))
            ENDDO
            sum_norms = sum_norms + local_sum
            max_local_sum = MAX(local_sum, max_local_sum)
        ENDDO

        IF (max_local_sum > eps) THEN
            ghmatecd_asserted = .FALSE.
            PRINT*, "||G||_inf = ", max_local_sum
        ENDIF

 200    FORMAT (A,ES7.1)     
       WRITE(0,"(A)") "As matrizes H e G em Ghmatecd_cu sao iguais as"
        WRITE(0,200) "calculadas em Ghmatece com um erro de ", eps
        IF (ghmatecd_asserted .EQV. .TRUE.) THEN
            WRITE(0,"(A)")"[OK]"
        ELSE
            WRITE(0,"(A)")"[FALHOU]"
        ENDIF
        WRITE(0,"(A)") ""

      END SUBROUTINE ASSERT_GHMATECE_H_G
