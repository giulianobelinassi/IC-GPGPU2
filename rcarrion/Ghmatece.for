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
        IMPLICIT REAL (A-H,O-Y)
        IMPLICIT COMPLEX (Z)
        COMMON INP,INQ,IPR,IPS,IPT
        DIMENSION CX(NCOX),CY(NCOX),CZ(NCOX)
        DIMENSION CXM(NE),CYM(NE),CZM(NE)
        DIMENSION HEST(NX,NX),GEST(NX,NX)
        DIMENSION HELEM(3,3),GELEM(3,3)
        DIMENSION DELTA(3,3)
        DIMENSION CO(4,3)!, ETA(3)
        INTEGER CONE(NE,4)
        REAL, INTENT(IN) :: ETAS(3,NX)
*
        PI=4.0*ATAN(1.0)
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
        DO 50 I=1,N
                II=3*(I-1) + 1
                GEST(II:II+2, II:II+2) = 0.0
                HEST(II:II+2, II:II+2) = 0.0
  50    CONTINUE
*
* CÁLCULO DOS COEFICIENTES DAS MATRIZES H E G
*
   
!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(N1,N2,N3,N4,J,I,CO,II,JJ,HELEM,GELEM)
        DO J=1,N
*
* CÁLCULO DAS COMPONENTES DO VETOR NORMAL
* USANDO O PRODUTO VETORIAL DOS LADOS 1-2 E 1-3
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
*                   ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE G SINGULAR
                    CALL SINGGE(HELEM,GELEM,CXM(I),CYM(I),CZM(I),
     $                      ETAS(1:3,J),CX,
     $                  CY,CZ,N1,N2,N3,N4,NCOX,N,NP,NPG,GE,RNU,RMU,C1,
     $                  C2,C3,C4,DELTA)
                ELSE

*                   ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE H E G NÃO SINGULAR
                    CALL NONSINGE(HELEM,GELEM,CO,CXM(I),CYM(I),CZM(I),
     $                  ETAS(1:3,J),N,NP,NPG,GE,RNU,RMU,C1,C2,C3,C4,
     $                  DELTA)
                ENDIF

                II=3*(I-1) + 1
                GEST(II:II+2, JJ:JJ+2) = GELEM
                HEST(II:II+2, JJ:JJ+2) = HELEM
*
            ENDDO
        ENDDO
!$OMP END PARALLEL DO
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
