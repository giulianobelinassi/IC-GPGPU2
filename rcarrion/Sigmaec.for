************************************************************************
*                                                                      *
*       ESTA SUBROTINA CALCULA OS VALORES DAS MATRIZES "S" E "D"       *
*                                                                      *
*			USADAS NO CÁLCULO DAS TENSÕES EM PONTOS INTERNOS           *
*                                                                      *
************************************************************************
*
       SUBROUTINE SIGMAEC(CO,CXP,CYP,CZP,RN,DELTA,PI,FR,ZGE,RHO,ZCS,
     $ZCP,NPG,ZD,ZS)
*
        IMPLICIT REAL*8 (A-H,O-Y)
        IMPLICIT COMPLEX*16 (Z)
        COMMON INP,INQ,IPR,IPS,IPT
        DIMENSION DELTA(3,3)
        DIMENSION GI(NPG),OME(NPG)
        DIMENSION CO(4,3),P(2,4),XJ(2,3),F(4)
        DIMENSION RD(3),RN(3)
        DIMENSION ZD(3,3,3),ZS(3,3,3),ZDP(3,3,3),ZSP(3,3,3)
C       DIMENSION ETA(3)

*
* ZERA AS MATRIZES ZD E ZS
*
        ZD = (0.D0, 0.D0)        
        ZS = (0.D0, 0.D0)
*
* ACIONA ROTINA QUE CALCULA OS PONTOS E PESOS DE GAUSS
*
        CALL GAULEG(-1.D0,1.D0,GI,OME,NPG)
*
        DO 400 JG=1,NPG
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
            DO 300 IG=1,NPG
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
                DO 90 I=1,2
                    DO 80 J=1,3
                        TEMP=0.D0
                        DO 70 K=1,4
                            TEMP=TEMP+P(I,K)*CO(K,J)
  70                    CONTINUE
                        XJ(I,J)=TEMP
  80                CONTINUE
  90            CONTINUE
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
 1100   FORMAT(/1X,'COORDINADAS DOS PONTOS EXTREMOS '/(3D14.5))
                     STOP ' NONSING : ERRO, JACOBIANO NULO OU NEGATIVO'
                ENDIF

* CALCULA AS COORDENADAS DO PONTO DE INTEGRAÇÃO
*
                CXG=0.D0
                CYG=0.D0
                CZG=0.D0
                DO 130 I=1,4
                    CXG=CXG+CO(I,1)*F(I)
                    CYG=CYG+CO(I,2)*F(I)
                    CZG=CZG+CO(I,3)*F(I)
 130            CONTINUE
*
* COMPRIMENTO DO RAIO
                R1=CXG-CXP
                R2=CYG-CYP
                R3=CZG-CZP
                R=DSQRT(R1**2+R2**2+R3**2)
*
* DERIVADA DO RAIO EM RELAÇÃO À NORMAL 
                DRN=(R1*RN(1)+R2*RN(2)+R3*RN(3))/R
*
* DERIVADAS DIRECIONAIS DO RAIO
                RD(1)=R1/R
                RD(2)=R2/R
                RD(3)=R3/R
*
* CONSTANTES USADAS NAS SOLUÇÕES FUNDAMENTAIS DINÂMICAS
*
                ZWI=DCMPLX(0.D0,FR)
                ZC0=(1.D0,0.D0)/(DCMPLX(4.D0*PI,0.D0)*ZGE)
                ZC1=(ZCP/ZCS)**2
                ZC2=(ZCS/ZCP)**2
                ZKP=-ZWI/ZCP
                ZKS=-ZWI/ZCS
                ZZP=ZKP*R
                ZZS=ZKS*R
                ZEZP=CDEXP(ZZP)
                ZEZS=CDEXP(ZZS)
                ZP2=ZZP*ZZP
                ZS2=ZZS*ZZS
*
                ZFHI=(1.D0+1.D0/ZS2-1.D0/ZZS)*ZEZS/R-
     @              ZC2*(1.D0/ZP2-1.D0/ZZP)*ZEZP/R
                ZCAPPA=(1.D0+3.D0/ZS2-3.D0/ZZS)*ZEZS/R-
     @              ZC2*(1.D0+3.D0/ZP2-3.D0/ZZP)*ZEZP/R
                ZFHIDR=(-2.D0+ZZS+3.D0/ZZS-3.D0/ZS2)*ZEZS/R**2-
     @              ZC2*(-1.D0+3.D0/ZZP-3.D0/ZP2)*ZEZP/R**2
                ZCAPPADR=(ZZS-4.D0+9.D0/ZZS-9.D0/ZS2)*ZEZS/R**2-
     @              ZC2*(ZZP-4.D0+9.D0/ZZP-9.D0/ZP2)*ZEZP/R**2
*
* ESTA PARTE FOI TIRADA DE UM TRABALHO DO PIOTR FEDELINSK
* QUE USA AS COSTANTES DE LAMÈ (LAMBDA E MI) E
* NÚMERO DE ONDAS (K1 E K2)
*
* VOU FAZER COMO ESTÁ NESTE TRABALHO POIS AS EQUAÇÕES SÃO MUITO GRANDES,
* EMBORA ALGUMAS VARIAVÉIS CRIADAS DAQUI PRA FRENTE JÁ EXISTAM
*
                ZK1=ZWI/ZCP
                ZK2=ZWI/ZCS
                ZEK1=CDEXP(-ZK1*R)
                ZEK2=CDEXP(-ZK2*R)
                ZMU=ZGE
                ZLAMB=(ZCP*ZCP*RHO)-(2.D0*ZMU)
*
                ZFHIDRR=-(ZCS/ZCP)**2*
     $              ZEK1*(12.D0+12.D0*ZK1*R+5.D0*(ZK1*R)**2+(ZK1*R)**3)*
     $              (1.D0/((ZK1*R)**2*R**3))+
     $          ZEK2*(12.D0+12.D0*ZK2*R+7.D0*(ZK2*R)**2+3.D0*(ZK2*R)**3+
     $              (ZK2*R)**4)*
     $              (1.D0/((ZK2*R)**2*R**3))
*
                ZCAPPADRR=-(ZCS/ZCP)**2*
     $         ZEK1*(36.D0+36.D0*ZK1*R+17.D0*(ZK1*R)**2+5.D0*(ZK1*R)**3+
     $              (ZK1*R)**4)*
     $              (1.D0/((ZK1*R)**2*R**3))+
     $         ZEK2*(36.D0+36.D0*ZK2*R+17.D0*(ZK2*R)**2+5.D0*(ZK2*R)**3+
     $              (ZK2*R)**4)*
     $              (1.D0/((ZK2*R)**2*R**3))
*
                ZCTE1=2.D0*(ZCAPPA/R)
                ZCTE2=(ZLAMB/ZMU)*(ZFHIDR-ZCAPPADR-ZCTE1)
                ZCTE3=(ZFHIDR-ZCAPPA/R)
                ZCTE4=2.D0*(ZCAPPADR-ZCTE1)

             ZCTE5=4.D0*(ZCAPPADRR-5.D0*(ZCAPPADR/R)+8.D0*(ZCAPPA/R**2))
           ZCTE6=ZFHIDRR-(ZFHIDR/R)-3.D0*(ZCAPPADR/R)+6.D0*(ZCAPPA/R**2)
                ZCTE7=2.D0*(ZCAPPADR/R)-4.D0*(ZCAPPA/R**2)
           ZCTE8=(ZLAMB/ZMU)*(ZCAPPADRR+(ZCAPPADR/R)-4.D0*(ZCAPPA/R**2)-
     $ZFHIDRR+(ZFHIDR/R))
                ZCTE9=2.D0*(ZCTE7+ZCTE8)
                ZCTE10=4.D0*(ZCAPPA/R**2)+4.D0*(ZLAMB/ZMU)*
     $((ZCAPPADR/R)+2.D0*(ZCAPPA/R**2)-(ZFHIDR/R))
                ZCTE11=((ZLAMB/ZMU)**2)*(ZCAPPADRR+4.D0*(ZCAPPADR/R)+
     $2.D0*(ZCAPPA/R**2)-ZFHIDRR-2.D0*(ZFHIDR/R))
                ZCTE12=2.D0*((ZFHIDR/R)-(ZCAPPA/R**2))

                DO 40 K=1,3
                    DO 40 I=1,3
                        DO 40 J=1,3
                            ZDP(K,I,J)=(1.D0/(4.D0*PI))*
     $                      (
     $                        (ZCTE1-ZCTE2)*DELTA(I,J)*RD(K)-
     $                        ZCTE3*(DELTA(I,K)*RD(J)+DELTA(J,K)*RD(I))+
     $                        ZCTE4*RD(I)*RD(J)*RD(K)
     $                      )
*
                ZSP(K,I,J)=(ZMU/(4.D0*PI))*
     $          (
     $              (ZCTE5*RD(I)*RD(J)*RD(K)-
     $                  ZCTE6*(DELTA(I,K)*RD(J)+DELTA(J,K)*RD(I))+
     $                  ZCTE9*DELTA(I,J)*RD(K))*DRN+
     $                  ZCTE9*RD(I)*RD(J)*RN(K)-
     $                  ZCTE6*(RD(J)*RN(I)+RD(I)*RN(J))*RD(K)+
     $                  (ZCTE10+ZCTE11)*DELTA(I,J)*RN(K)-
     $                  ZCTE12*(DELTA(J,K)*RN(I)+DELTA(I,K)*RN(J))
     $          )
 40             CONTINUE
*
                P12=P1*P2*DET
                DO 200 KN=1,3
                    DO 200 IN=1,3
                        DO 200 JN=1,3
                            ZD(KN,IN,JN)=ZD(KN,IN,JN)+ZDP(KN,IN,JN)*P12
                            ZS(KN,IN,JN)=ZS(KN,IN,JN)+ZSP(KN,IN,JN)*P12
 200            CONTINUE
*
 300        CONTINUE
 400    CONTINUE
*
        RETURN
      END
