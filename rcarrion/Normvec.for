! Subrotina que calcula o produto vetorial usado originalmente em Ghmatecd.for,
! Ghmatece.for e Interec.for, evitando cálculos redundantes.

      SUBROUTINE NORMVEC(cone, cx, cy, cz, nx, ne, ncox, n, etas)
        IMPLICIT NONE
        INTEGER, DIMENSION(NE,4), INTENT(IN) :: CONE
        DOUBLE PRECISION, DIMENSION(NCOX), INTENT(IN) :: CX, CY, CZ
        INTEGER, INTENT(IN) :: NX, NE, NCOX, N
        DOUBLE PRECISION, DIMENSION(3,nx), INTENT(INOUT) :: ETAS

        INTEGER J, N1,N2,N3,N4
        DOUBLE PRECISION A, B, C, R

!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(N1,N2,N3,N4,J,A,B,C,R)
        DO J=1, N
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
            ETAS(1,j)=A/R
            ETAS(2,j)=B/R
            ETAS(3,j)=C/R
        ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE NORMVEC
