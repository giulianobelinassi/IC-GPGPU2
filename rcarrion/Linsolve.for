      SUBROUTINE LINSOLVE(NN, N, ZH, ZG, ZFI, DFI, ZDFI)
        USE omp_lib
    
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NN, N
        COMPLEX, INTENT(INOUT) :: ZH(NN, NN), ZG(NN, NN)
        REAL, INTENT(IN) :: DFI(NN)
        COMPLEX, INTENT(INOUT), ALLOCATABLE :: ZFI(:)
        COMPLEX, INTENT(INOUT), ALLOCATABLE :: ZDFI(:)
		INTEGER :: stats1, stats2
        DOUBLE PRECISION :: t1, t2
#define USE_CPU

#ifdef USE_GPU
#undef USE_CPU
#endif

#ifdef  TEST_CUDA
#undef  USE_CPU
#undef  USE_GPU
#define USE_CPU
#define USE_GPU
        COMPLEX :: ZFI_ORIG(NN), ZFIP(NN)
        COMPLEX, ALLOCATABLE :: ZH_ORIG(:,:), ZHP(:,:)
#endif

#ifdef USE_CPU
        INTEGER :: stats 
        INTEGER, ALLOCATABLE :: PIV(:)
#endif

        ALLOCATE(ZDFI(NN), STAT = stats1)
        ALLOCATE(ZFI(NN) , STAT = stats2)
        IF (stats1 /= 0 .OR. stats2 /= 0) THEN
            PRINT*, "MEMÓRIA INSUFICIENTE"
        ENDIF

!
! TRANSFORMAÇÃO DAS CONDIÇÕES DE CONTORNO EM NÚMEROS COMPLEXOS
!
        ZDFI = DFI


! FORMA O LADO DIREITO DO SISTEMA {VETOR f} QUE É ARMAZENADO EM ZFI

#ifdef TEST_CUDA
        ALLOCATE(ZH_ORIG(NN, NN))
        ALLOCATE(ZHP(NN, NN))
        ZH_ORIG = ZH
        ZFI_ORIG = ZFI
#endif

#ifdef USE_CPU

        t1 = OMP_GET_WTIME()
        
        IF (SIZEOF(1.0) == 8) THEN
            CALL ZGEMV('N',NN,NN,(1.0,0),ZG,NN,ZDFI,1,(0,0),ZFI,1)
        ELSEIF (SIZEOF(1.0) == 4) THEN
            CALL CGEMV('N',NN,NN,(1.0,0),ZG,NN,ZDFI,1,(0,0),ZFI,1)
        ELSE
            PRINT*, "ERRO FATAL: Precisão desconhecida"
            STOP
        ENDIF
        
        ALLOCATE(PIV(NN), STAT = stats)
        IF (stats /= 0) THEN
            PRINT*, "MEMORIA INSUFICIENTE"
            STOP
        ENDIF

        IF (SIZEOF(1.0) == 8) THEN
            CALL ZGESV(NN,1,ZH,NN,PIV,ZFI,NN,stats)
        ELSEIF (SIZEOF(1.0) == 4) THEN
            CALL CGESV(NN,1,ZH,NN,PIV,ZFI,NN,stats)
        ELSE
            PRINT*, "ERRO FATAL: Precisão de float desconhecida"
        ENDIF

        IF (stats < 0) THEN
            PRINT *, "Erro em ZGESV :-("
        ELSE IF (stats > 0) THEN
            PRINT *, "Matriz Singular :-|"
        ENDIF
        DEALLOCATE(PIV)
        t2 = OMP_GET_WTIME()
        PRINT*, "LINSOLVE: Tempo na CPU: ", (t2-t1)
#endif

#ifdef  TEST_CUDA
        
        ZFIP = ZFI
        ZH = ZHP
        ZFI = ZFI_ORIG
        ZH = ZH_ORIG
#endif

#ifdef USE_GPU
        t1 = OMP_GET_WTIME()  
        CALL cuda_linsolve(NN, N, ZFI, ZDFI)
        t2 = OMP_GET_WTIME()
        PRINT*, "LINSOLVE: Tempo na GPU: ", (t2-t1)
#endif

#ifdef  TEST_CUDA
        CALL ASSERT_ZFI(ZFI, ZFIP, NN)
        DEALLOCATE(ZHP)
        DEALLOCATE(ZH_ORIG)
#endif

      END SUBROUTINE


      SUBROUTINE ASSERT_ZFI(ZFI, ZFIP, NN)
        IMPLICIT NONE
        COMPLEX, DIMENSION(NN), INTENT(IN) :: ZFI, ZFIP
        INTEGER :: NN, i
        LOGICAL :: asserted = .TRUE. 
        REAL, PARAMETER :: eps = 1.2E-6
        REAL :: maxentry = 0

        DO i = 1, NN
            maxentry = MAX(maxentry, ABS(ZFI(i) - ZFIP(i)))
!            PRINT*, ZFI(i), ZFIP(i)
        ENDDO

        IF (maxentry > eps) THEN
            asserted = .FALSE.
        ENDIF

        PRINT*, "||ZFI||_inf = ", maxentry


 200    FORMAT (A,ES7.1)     
        WRITE(0,"(A)") "O vetor ZFI calculado em Interec1_cu e igual ao"
        WRITE(0,200) "calculado em Interec.for com um erro de ", eps
        IF (asserted .EQV. .TRUE.) THEN
            WRITE(0,"(A)")"[OK]"
        ELSE
            WRITE(0,"(A)")"[FALHOU]"
        ENDIF
        WRITE(0,"(A)") ""
      END SUBROUTINE ASSERT_ZFI
