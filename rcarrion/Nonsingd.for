************************************************************************
*                                                                      *
*  SUBROTINA QUE CALCULA OS COEFICIENTES DE [H] E [G]. O PONTO DE      *
*																	   *
*    COLOCAÇÃO ESTÁ FORA DO ELEMENTO, ASSIM NÃO HÁ SINGULARIDADE.      *
*                                                                      *
*         A INTEGRAÇÃO É FEITA USANDO QUADRATURA GAUSSIANA             *
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
C        DIMENSION CO_cp(4,3),ETA_cp(3),DELTA_cp(3,3)    

C Guarda as entradas
C        CO_cp  = CO
C        CXP_cp = CXP
C        CYP_cp = CYP
C        CZP_cp = CZP
C        ETA_cp = ETA
C        ZGE_cp = ZGE
C        ZSC_cp = ZSC
C        ZCP_cp = ZCP
C        DELTA_cp = DELTA
C        PI_cp = PI
C        FR_co = FR
C        NPG_cp = NPG

* ZERA AS MATRIZES ELEMENTARES HELEM E GELEM
*       
        ZHELEM = 0.D0
        ZGELEM = 0.D0
*
* ACIONA ROTINA QUE CALCULA OS PONTOS E PESOS DE GAUSS
*
        CALL GAULEG(-1.D0,1.D0,GI,OME,NPG)
*
!$OMP  PARALLEL DEFAULT(PRIVATE) 
!$OMP& PRIVATE(JG, IG)
!$OMP& SHARED(INP,INQ,IPR,IPS,IPT,GI,OME)
!$OMP& REDUCTION(+:ZHELEM, ZGELEM)

!$OMP DO
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
* CALCULA A RELAÇÃO ENTRE AS COORDENADAS CARTESIANAS E HOMOGÊNEAS
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
* CALCULA AS COORDENADAS DO PONTO DE INTEGRAÇÃO
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
* ACIONA ROTINA QUE CALCULA A SOLUÇÃO FUNDAMENTAL ESTÁTICA 3D
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
!$OMP END DO
!$OMP END PARALLEL

*
C        DO 321 j=1,3
C            DO 321 i=1,4
C                IF (CO_cp(i, j) /= CO(i, j)) THEN           
C                    PRINT*, "Variável Modificada: CO"
C                ENDIF
C 321    CONTINUE
C        DO 123 j=1,3
C            DO 123, i=1,3
C                IF (DELTA_cp(i, j) /= DELTA(i, j)) THEN 
C                    PRINT*, "Variável modificada: DELTA"
C                ENDIF
C 123    CONTINUE
C
C        DO 213 i=1,3
C            IF (ETA_cp(i) /= ETA_cp(i)) THEN
C                PRINT*, "Variável Modificada: ETA"
C            ENDIF
C 213    CONTINUE
C
C        IF (CXP_cp /= CXP) THEN
C            PRINT*, "Variável Modificada: CXP"
C        ENDIF
C        IF (CYP_cp /= CYP) THEN
C            PRINT*, "Variável Modificada: CYP"
C        ENDIF
C        IF (CZP_cp /= CZP) THEN
C            PRINT*, "Variável Modificada: CZP"
C        ENDIF
C        IF (ZGE_cp /= ZGE) THEN
C            PRINT*, "Variável Modificada: ZGE"
C        ENDIF
C        IF (ZSC_cp /= ZSC) THEN
C            PRINT*, "Variável Modificada: ZSC"
C        ENDIF
C        IF (ZCP_cp /= ZCP) THEN
C            PRINT*, "Variável Modificada: ZCP"
C        ENDIF
C        IF (PI_cp /= PI) THEN
C            PRINT*, "Variável Modificada: PI"
C        ENDIF
C        IF (FR_co /= FR) THEN
C            PRINT*, "Variável Modificada: FR"
C        ENDIF
C        IF (NPG_cp /= NPG) THEN
C            PRINT*, "Variável Modificada: NPG"
C        ENDIF
        
        
        
        
        
        
        
        
        
        
        
        
        RETURN
      END
