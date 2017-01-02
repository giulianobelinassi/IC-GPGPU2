************************************************************************
*                                                                      *
*  SUBROTINA QUE CALCULA OS COEFICIENTES DE [H] E [G]. O PONTO DE      *
*																	   *
*    COLOCA��O EST� FORA DO ELEMENTO, ASSIM N�O H� SINGULARIDADE.      *
*                                                                      *
*         A INTEGRA��O � FEITA USANDO QUADRATURA GAUSSIANA             *
*                                                                      *
************************************************************************
*
      SUBROUTINE NONSINGD(ZHELEM,ZGELEM,CO,CXP,CYP,CZP,ETA,ZGE,ZCS,ZCP,
     $DELTA,PI,FR,NPG)
*
        IMPLICIT NONE
        COMMON INP,INQ,IPR,IPS,IPT
        DOUBLE COMPLEX, INTENT(OUT) :: ZHELEM(3,3), ZGELEM(3,3)
        DOUBLE PRECISION, INTENT(IN):: CO(4,3), CXP, CYP, CZP, ETA(3)
        DOUBLE COMPLEX, INTENT(IN ) :: ZGE, ZCS, ZCP
        DOUBLE PRECISION, INTENT(IN):: DELTA, PI, FR
        INTEGER, INTENT(IN)         :: NPG
        INTEGER INP, INQ, IPR, IPS, IPT

        DOUBLE COMPLEX ZU(3,3), ZT(3,3)
        DOUBLE PRECISION GI(NPG), OME(NPG), P(2,4), XJ(2,3), F(4)
        DOUBLE PRECISION G1, G2, P1, P2, P12, SP, SM, RP, RM, TEMP, DET
        DOUBLE PRECISION CXG, CYG, CZG
        INTEGER i, j, k, ig, jg

*

* ZERA AS MATRIZES ELEMENTARES HELEM E GELEM
*       
        ZHELEM = 0.D0
        ZGELEM = 0.D0
*
* ACIONA ROTINA QUE CALCULA OS PONTOS E PESOS DE GAUSS
*
        CALL GAULEG(-1.D0,1.D0,GI,OME,NPG)
*
C!$OMP  PARALLEL DO DEFAULT(PRIVATE)
C!$OMP& SHARED(INP,INQ,IPR,IPS,IPT, GI, OME, CO, CXP, CYP, CZP, NPG, ETA,
C!$OMP& DELTA, ZGE, PI, FR, ZCS, ZCP)
C!$OMP& REDUCTION(+:ZHELEM, ZGELEM)
        DO JG=1,NPG
*
            G2=GI(JG)
            P2=OME(JG)
            SP=1.D0+G2
            SM=1.D0-G2
            P(1,1)=-0.25D0*SM
            P(1,2)= 0.25D0*SM
            P(1,3)= 0.25D0*SP
            P(1,4)=-0.25D0*SP
*
            DO IG=1,NPG
*
                IF (IG > NPG) THEN
                    PRINT*, "Vai crashar :P"
                ENDIF

                G1=GI(IG)
                P1=OME(IG)
                RP=1.D0+G1
                RM=1.D0-G1
                F(1)=0.25D0*RM*SM
                F(2)=0.25D0*RP*SM
                F(3)=0.25D0*RP*SP
                F(4)=0.25D0*RM*SP
                P(2,1)=-0.25D0*RM
                P(2,2)=-0.25D0*RP
                P(2,3)= 0.25D0*RP
                P(2,4)= 0.25D0*RM
*
* CALCULA A RELA��O ENTRE AS COORDENADAS CARTESIANAS E HOMOG�NEAS
*
                DO  I=1,2
                    DO J=1,3
                        TEMP=0.D0
                        DO K=1,4
                            TEMP=TEMP+P(I,K)*CO(K,J)
                        ENDDO
                        XJ(I,J)=TEMP
                    ENDDO
                ENDDO
*
* CALCULA O JACOBIANO
*
                DET=DSQRT((XJ(1,2)*XJ(2,3)-XJ(2,2)*XJ(1,3))**2 +
     $              (XJ(2,1)*XJ(1,3)-XJ(1,1)*XJ(2,3))**2 +
     $              (XJ(1,1)*XJ(2,2)-XJ(2,1)*XJ(1,2))**2  )

                IF (DET < 1.0D-5) THEN
                    WRITE(IPR,1000) DET
                    WRITE(IPR,1100) ((CO(I,J),J=1,3),I=1,4)
 1000   FORMAT(///' NONSING : ERRO, JACOBIANO NULO OU NEGATIVO =',D14.5)
 1100            FORMAT(/1X,'COORDINADAS DOS PONTOS EXTREMOS '/(3D14.5))
                    STOP ' NONSING : ERRO, JACOBIANO NULO OU NEGATIVO'
                ENDIF
*
* CALCULA AS COORDENADAS DO PONTO DE INTEGRA��O
*
                CXG=0.D0
                CYG=0.D0
                CZG=0.D0
                DO I=1,4
                    CXG=CXG+CO(I,1)*F(I)
                    CYG=CYG+CO(I,2)*F(I)
                    CZG=CZG+CO(I,3)*F(I)
                ENDDO
*
* ACIONA ROTINA QUE CALCULA A SOLU��O FUNDAMENTAL EST�TICA 3D
*
                CALL SOLFUND(ZU,ZT,CXP,CYP,CZP,CXG,CYG,CZG,ETA,ZGE,ZCS,
     $ZCP,DELTA,PI,FR)
*
                P12=P1*P2*DET
                ZHELEM = ZHELEM + ZT*P12
                ZGELEM = ZGELEM + ZU*P12
*
            ENDDO
        ENDDO
C!$OMP END PARALLEL DO

        RETURN
      END
