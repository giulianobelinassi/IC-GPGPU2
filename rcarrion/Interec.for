************************************************************************
*                                                                      *
*       ESTA SUBROTINA CALCULA OS VALORES DOS DESLOCAMENTOS            *
*                                                                      *
*                        NOS PONTOS INTERNOS                           *
*                                                                      *
************************************************************************
*
      SUBROUTINE INTEREC(ZFI,ZDFI,KODE,CX,CY,CZ,CXI,CYI,CZI,ZDSOL,ZSSOL,
     $  NE,NX,NCOX,NPIX,CONE,ZGE,ZCS,ZCP,DELTA,PI,FR,NPG,L,N,NBE,RHO)
*
        IMPLICIT REAL*8 (A-H,O-Y)
        IMPLICIT COMPLEX*16 (Z)
        COMMON INP,INQ,IPR,IPS,IPT
        DIMENSION CX(NCOX),CY(NCOX),CZ(NCOX)
        DIMENSION CXI(NPIX),CYI(NPIX),CZI(NPIX)
        DIMENSION ZHELEM(3,3),ZGELEM(3,3)
        DIMENSION CO(4,3),ETA(3)
        DIMENSION DELTA(3,3)
        DIMENSION ZFI(NX),ZDFI(NX),ZDSOL(3*NPIX),ZSSOL(9*NPIX)
        INTEGER CONE(NE,4),KODE(NX)
        DIMENSION ZD(3,3,3),ZS(3,3,3)
*
* REARRANJA OS VETORES ZFI AND ZDFI PARA ARMAZENAR TODOS OS VALORES DOS
*
* DESLOCAMENTOS EM ZFI E TODOS OS VALORES DE FORÇAS DE SUPERFÍCIE EM ZDFI
*
        NN=3*NBE
        DO 20 I=1,NN

            IF (KODE(I) == 0) THEN
                ZCH = ZFI(I)*ZGE
                ZFI(I)=ZDFI(I)
                ZDFI(I)=ZCH
            ENDIF
  20    CONTINUE
*
* CALCULA OS VALORES DOS DESLOCAMENTOS EM PONTOS INTERNOS
*
      IF(L <= 0) RETURN
*
        DO 30 K=1,L*3
            ZDSOL(K)=(0.D0,0.D0)
  30    CONTINUE



!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(N1,N2,N3,N4,J,JJ,K,KK,CO,ETA,A,B,C,R,ZHELEM,ZGELEM)
!$OMP& REDUCTION(+:ZDSOL)
        DO J=1,N
*
            N1=CONE(J,1)
            N2=CONE(J,2)
            N3=CONE(J,3)
            N4=CONE(J,4)
            A=(CY(N2)-CY(N1))*(CZ(N3)-CZ(N1)) - 
     $          (CZ(N2)-CZ(N1))*(CY(N3)-CY(N1))
            B=(CZ(N2)-CZ(N1))*(CX(N3)-CX(N1)) - 
     $          (CX(N2)-CX(N1))*(CZ(N3)-CZ(N1))
            C=(CX(N2)-CX(N1))*(CY(N3)-CY(N1)) - 
     $          (CY(N2)-CY(N1))*(CX(N3)-CX(N1))
            R=DSQRT(A*A+B*B+C*C)
      
            ETA(1)=A/R
            ETA(2)=B/R
            ETA(3)=C/R
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
            DO K=1,L
*
                CALL NONSINGD(ZHELEM,ZGELEM,CO,CXI(K),CYI(K),CZI(K),ETA,
     $              ZGE,ZCS,ZCP,DELTA,PI,FR,NPG)
*
                KK=3*(K-1)
                ZDSOL(KK+1:KK+3)=ZDSOL(KK+1:KK+3) +
     $              MATMUL(ZGELEM,ZDFI(JJ+1:JJ+3))-
     $              MATMUL(ZHELEM,ZFI(JJ+1:JJ+3))

C                DO IN=1,3
C                    DO JN=1,3
C                        ZDSOL(KK+IN)=ZDSOL(KK+IN)+
C     $                      ZDFI(JJ+JN)*ZGELEM(IN,JN)-
C     $                      ZFI(JJ+JN)*ZHELEM(IN,JN)
C                    ENDDO
C                ENDDO
            ENDDO
        ENDDO
!$OMP END PARALLEL DO
*
* ACRESCENTADO POSTERIORMANTE (APÓS O PROGRAMA ESTAR RODANDO ATÉ
* O CÁLCULO PARA OS DESLOCAMENTOS EM PONTOS INTERNOS).
*
* CÁLCULO DAS TENSÕES EM PONTOS INTERNOS
*
       
        ZSSOL(1:L*9) = 0.D0
       
!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(N1,N2,N3,N4,J,K,CO,ETA,A,B,C,R,ZD,ZS)
!$OMP& REDUCTION(+:ZDSOL)
        DO K=1,L
            DO J=1,N
*
                N1=CONE(J,1)
                N2=CONE(J,2)
                N3=CONE(J,3)
                N4=CONE(J,4)
                A=(CY(N2)-CY(N1))*(CZ(N3)-CZ(N1)) - 
     $              (CZ(N2)-CZ(N1))*(CY(N3)-CY(N1))
                B=(CZ(N2)-CZ(N1))*(CX(N3)-CX(N1)) - 
     $              (CX(N2)-CX(N1))*(CZ(N3)-CZ(N1))
                C=(CX(N2)-CX(N1))*(CY(N3)-CY(N1)) - 
     $              (CY(N2)-CY(N1))*(CX(N3)-CX(N1))
                
                R=DSQRT(A*A+B*B+C*C)
                ETA(1)=A/R
                ETA(2)=B/R
                ETA(3)=C/R
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
*
                CALL SIGMAEC(CO,CXI(K),CYI(K),CZI(K),ETA,DELTA,PI,FR,
     $              ZGE,RHO,ZCS,ZCP,NPG,ZD,ZS)
*
                ZSSOL(9*K-8)=ZSSOL(9*K-8)+ZDFI(3*J-2)*ZD(1,1,1)+
     $                                    ZDFI(3*J-1)*ZD(2,1,1)+
     $                                    ZDFI(3*J)  *ZD(3,1,1)-
     $                                    ZFI (3*J-2)*ZS(1,1,1)-
     $                                    ZFI (3*J-1)*ZS(2,1,1)-
     $                                    ZFI (3*J)  *ZS(3,1,1)
                ZSSOL(9*K-7)=ZSSOL(9*K-7)+ZDFI(3*J-2)*ZD(1,1,2)+
     $                                    ZDFI(3*J-1)*ZD(2,1,2)+
     $                                    ZDFI(3*J)  *ZD(3,1,2)-
     $                                    ZFI (3*J-2)*ZS(1,1,2)-
     $                                    ZFI (3*J-1)*ZS(2,1,2)-
     $                                    ZFI (3*J)  *ZS(3,1,2)
                ZSSOL(9*K-6)=ZSSOL(9*K-6)+ZDFI(3*J-2)*ZD(1,1,3)+
     $                                    ZDFI(3*J-1)*ZD(2,1,3)+
     $                                    ZDFI(3*J)  *ZD(3,1,3)-
     $                                    ZFI (3*J-2)*ZS(1,1,3)-
     $                                    ZFI (3*J-1)*ZS(2,1,3)-
     $                                    ZFI (3*J)  *ZS(3,1,3)
                ZSSOL(9*K-5)=ZSSOL(9*K-5)+ZDFI(3*J-2)*ZD(1,2,1)+
     $                                    ZDFI(3*J-1)*ZD(2,2,1)+
     $                                    ZDFI(3*J)  *ZD(3,2,1)-
     $                                    ZFI (3*J-2)*ZS(1,2,1)-
     $                                    ZFI (3*J-1)*ZS(2,2,1)-
     $                                    ZFI (3*J)  *ZS(3,2,1)
                ZSSOL(9*K-4)=ZSSOL(9*K-4)+ZDFI(3*J-2)*ZD(1,2,2)+
     $                                    ZDFI(3*J-1)*ZD(2,2,2)+
     $                                    ZDFI(3*J)  *ZD(3,2,2)-
     $                                    ZFI (3*J-2)*ZS(1,2,2)-
     $                                    ZFI (3*J-1)*ZS(2,2,2)-
     $                                    ZFI (3*J)  *ZS(3,2,2)
                ZSSOL(9*K-3)=ZSSOL(9*K-3)+ZDFI(3*J-2)*ZD(1,2,3)+
     $                                    ZDFI(3*J-1)*ZD(2,2,3)+
     $                                    ZDFI(3*J)  *ZD(3,2,3)-
     $                                    ZFI (3*J-2)*ZS(1,2,3)-
     $                                    ZFI (3*J-1)*ZS(2,2,3)-
     $                                    ZFI (3*J)  *ZS(3,2,3)
                ZSSOL(9*K-2)=ZSSOL(9*K-2)+ZDFI(3*J-2)*ZD(1,3,1)+
     $                                    ZDFI(3*J-1)*ZD(2,3,1)+
     $                                    ZDFI(3*J)  *ZD(3,3,1)-
     $                                    ZFI (3*J-2)*ZS(1,3,1)-
     $                                    ZFI (3*J-1)*ZS(2,3,1)-
     $                                    ZFI (3*J)  *ZS(3,3,1)
                ZSSOL(9*K-1)=ZSSOL(9*K-1)+ZDFI(3*J-2)*ZD(1,3,2)+
     $                                    ZDFI(3*J-1)*ZD(2,3,2)+
     $                                    ZDFI(3*J)  *ZD(3,3,2)-
     $                                    ZFI (3*J-2)*ZS(1,3,2)-
     $                                    ZFI (3*J-1)*ZS(2,3,2)-
     $                                    ZFI (3*J)  *ZS(3,3,2)
                ZSSOL(9*K)=  ZSSOL(9*K)  +ZDFI(3*J-2)*ZD(1,3,3)+
     $                                    ZDFI(3*J-1)*ZD(2,3,3)+
     $                                    ZDFI(3*J)  *ZD(3,3,3)-
     $                                    ZFI (3*J-2)*ZS(1,3,3)-
     $                                    ZFI (3*J-1)*ZS(2,3,3)-
     $                                    ZFI (3*J)  *ZS(3,3,3)
            ENDDO
        ENDDO        
!$OMP END PARALLEL DO
*
        RETURN
      END
