************************************************************************
*                                                                      *
*  SUBROTINA QUE CALCULA AS SOLUÇÕES FUNDAMENTAIS DINÂMICAS 3D		 *
*																	 *
************************************************************************
*
      SUBROUTINE SOLFUND(ZU,ZT,CXP,CYP,CZP,CXG,CYG,CZG,RN,ZGE,ZCS,ZCP,
     $DELTA,PI,FR)
*
	IMPLICIT REAL*8 (A-H,O-Y)
	IMPLICIT COMPLEX*16 (Z)
	COMMON INP,INQ,IPR,IPS,IPT
	DIMENSION ZU(3,3),ZT(3,3)
	DIMENSION RD(3)
	DIMENSION RN(3)
	DIMENSION DELTA(3,3)
*
      R1=CXG-CXP
      R2=CYG-CYP
      R3=CZG-CZP
      R=DSQRT(R1**2+R2**2+R3**2)
      DRN=(R1*RN(1)+R2*RN(2)+R3*RN(3))/R
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
     @    ZC2*(1.D0/ZP2-1.D0/ZZP)*ZEZP/R
      ZCAPPA=(1.D0+3.D0/ZS2-3.D0/ZZS)*ZEZS/R-
     @      ZC2*(1.D0+3.D0/ZP2-3.D0/ZZP)*ZEZP/R
      ZFHIDR=(-2.D0+ZZS+3.D0/ZZS-3.D0/ZS2)*ZEZS/R**2-
     @       ZC2*(-1.D0+3.D0/ZZP-3.D0/ZP2)*ZEZP/R**2
      ZCAPPADR=(ZZS-4.D0+9.D0/ZZS-9.D0/ZS2)*ZEZS/R**2-
     @        ZC2*(ZZP-4.D0+9.D0/ZZP-9.D0/ZP2)*ZEZP/R**2
*
      ZAA=ZFHIDR-ZCAPPA/R
      ZBB=4.D0*ZCAPPA/R-2.D0*ZCAPPADR
      ZCC=(ZC1-2.D0)*(ZAA+5.D-1*ZBB-3.D0*ZCAPPA/R)-2.D0*ZCAPPA/R      
*
      DO 40 I=1,3
      DO 40 J=1,3
      ZU(I,J)=ZC0*(ZFHI*DELTA(I,J)-ZCAPPA*RD(J)*RD(I))
      ZT(I,J)=(1.D0/(4.D0*PI))*((ZAA*(DRN*DELTA(I,J)+RD(J)*RN(I)))+
     $RD(I)*RD(J)*DRN*ZBB+RD(I)*RN(J)*ZCC)
 40   CONTINUE
*
	RETURN
	END