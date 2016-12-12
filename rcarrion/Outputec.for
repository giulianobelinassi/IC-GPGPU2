************************************************************************
*                                                                      *
*        ESTA SUBROTINA IMPRIME OS VALORES DOS DESLOCAMENTOS E         *
*                                                                      *
*              FOR«AS DE SUPERFÕCIE NOS N”S DO CONTORNO                *
*                                                                      *
************************************************************************
*
      SUBROUTINE OUTPUTEC(ZFI,ZDFI,ZDSOL,ZSSOL,NPIX,NX,N,NBE,L,FR)
*
	IMPLICIT REAL*8 (A-H,O-Y)
	IMPLICIT COMPLEX*16 (Z)
	COMMON INP,INQ,IPR,IPS,IPT
	DIMENSION ZFI(NX),ZDFI(NX),ZDSOL(3*NPIX),ZSSOL(9*NPIX)
	DIMENSION DELTA(3,3)
*
      WRITE(IPS,50)
  50  FORMAT(/' ',79('*')/)
*
      WRITE(IPS,700),FR
 700  FORMAT(30X,'FREQ‹ NCIA = ',D14.7/)
*
      WRITE(IPS,100)
 100  FORMAT(' ',79('*')//1X,'RESULTADOS'//2X,'N”S DO CONTORNO'//' ELEM',
     $6X,'DESLOCAMENTO X',11X,'DESLOCAMENTOT Y',11X,'DESLOCAMENTO Z'/)
      DO 10 I=1,NBE
      WRITE(IPS,200) I,(ZFI(3*(I-1)+J),J=1,3)
  10  CONTINUE
*
      WRITE(IPS,150)
 150  FORMAT(/' ELEM',8X,'TRACTION X',15X,'TRACTION Y',15X
     $,'TRACTION Z'/)
      DO 15 I=1,NBE
      WRITE(IPS,200) I,(ZDFI(3*(I-1)+J),J=1,3)
  15  CONTINUE   
 200  FORMAT(I5,3(2X,D11.4,',',D11.4))
*
      IF(L.LE.0) GO TO 60
      WRITE(IPS,300)
 300  FORMAT(//2X,'PONTOS INTERNOS'//' PONTO',5X,
     $ 'DESLOCAMENTO X',11X,'DESLOCAMENTO Y',11X,'DESLOCAMENTO Z')
      DO 20 I=1,L
      WRITE(IPS,200) I,(ZDSOL(3*(I-1)+J),J=1,3)
  20  CONTINUE
*
      IF(L.LE.0) GO TO 60
      WRITE(IPS,400)
 400  FORMAT(//2X,'PONTOS INTERNOS'//' PONTO',8X,
     $'SIGMAXX',17X,'SIGMAXY',18X,'SIGMAXZ',18X,
     $'SIGMAYX',18X,'SIGMAYY',18X,'SIGMAYZ',18X,
     $'SIGMAZX',18X,'SIGMAZY',18X,'SIGMAZZ')
      DO 25 I=1,L
      WRITE(IPS,350) I,(ZSSOL(9*(I-1)+J),J=1,9)
 350	FORMAT(I4,9(2X,D11.4,',',D11.4))
  25  CONTINUE

*
	WRITE(IPT,*)FR,CDABS(ZFI(157)),CDABS(ZFI(166))
*
  60  WRITE(IPS,500)
 500  FORMAT(' ',230('*'))
*
      RETURN
      END