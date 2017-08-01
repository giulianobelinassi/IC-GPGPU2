      SUBROUTINE LINSOLVE(NN, N, ZH, ZFI)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NN, N
        COMPLEX, INTENT(IN) :: ZH(NN, 3*N)
        COMPLEX, INTENT(INOUT) :: ZFI(NN)
        
#define USE_CPU

#ifdef USE_GPU
#undef USE_CPU
#endif

#ifdef  TEST_CUDA
#undef  USE_CPU
#undef  USE_GPU
#define USE_CPU
#define USE_GPU
        COMPLEX :: ZFIP(NN)
#endif

#ifdef USE_CPU
        INTEGER :: stats
        INTEGER, ALLOCATABLE :: PIV(:)

        ALLOCATE(PIV(NN), STAT = stats)
        IF (stats /= 0) THEN
            PRINT*, "MEMORIA INSUFICIENTE"
            STOP
        ENDIF

        CALL CGESV(NN,1,ZH,NN,PIV,ZFI,NN,stats)
        IF (stats < 0) THEN
            PRINT *, "Erro em ZGESV :-("
        ELSE IF (stats > 0) THEN
            PRINT *, "Matriz Singular :-|"
        ENDIF
        DEALLOCATE(PIV)
#endif

#ifdef  TEST_CUDA
        ZFIP = ZFI
#endif

#ifdef USE_GPU
        CALL cuda_linsolve(NN, N, ZH, ZFI) 
#endif

#ifdef  TEST_CUDA
        CALL ASSERT_ZFI(ZFI, ZFIP, NN) 
#endif

      END SUBROUTINE


      SUBROUTINE ASSERT_ZFI(ZFI, ZFIP, NN)
        IMPLICIT NONE
        COMPLEX, DIMENSION(NN), INTENT(IN) :: ZFI, ZFIP
        INTEGER :: NN, i
        LOGICAL :: asserted = .TRUE. 
        REAL, PARAMETER :: eps = 1.0E-6
        REAL :: maxentry = 0

        DO i = 1, NN
            maxentry = MAX(maxentry, ABS(ZFI(i) - ZFIP(i)))
!            PRINT*, ZFI(i), ZFIP(i)
        ENDDO

        IF (maxentry > eps) THEN
            asserted = .FALSE.
            PRINT*, "||ZFI||_inf = ", maxentry
        ENDIF

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
