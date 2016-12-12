************************************************************************
************************************************************************
**                                                                    **
**     ESTE PROGRAMA RESOLVE PROBLEMAS VISCOELASTICOS, DIN�MICOS, 	**
**                                                                    **
**   TRIDIMENSIONAIS ATRAV�S DO M�TODO DOS ELEMENTOS DE CONTORNO NUM	**
**																	**
**						SEMI ESPA�O INFINITO						**
**																	**
************************************************************************
**                                                                    **
** PROGRAMA "MECECCS" = (M�TODO DOS ELEMENTOS DE CONTORNO USANDO      **
**                                                                    **
**					   ELEMENTOS CONSTANTES COM SINGULARIDADE)		**
**																	**
**								DIN�MICO							**
**																	**
** RONALDO CARRION / EUCLIDES DE MESQUITA NETO                        **
**                                                                    **
** DATA: 18/07/01														**
**																	**
************************************************************************
************************************************************************
*
	PROGRAM MAIN
C		USE MSIMSL
*
		IMPLICIT REAL*8 (A-H,O-Y)
		IMPLICIT COMPLEX*16 (Z)
		PARAMETER (NE=960,NX=3*NE,NCOX=962,NPIX=10,NFRX=7)
		COMMON INP,INQ,IPR,IPS,IPT
		DIMENSION CX(NCOX),CY(NCOX),CZ(NCOX)
		DIMENSION CXM(NE),CYM(NE),CZM(NE)
		DIMENSION CXI(NPIX),CYI(NPIX),CZI(NPIX)
		DIMENSION HEST(NX,NX),GEST(NX,NX)
		DIMENSION ZH(NX,NX),ZG(NX,NX)
		DIMENSION ZDFI(NX),ZFI(NX),ZVETSOL(NX),ZDSOL(3*NPIX),ZSSOL(9*NPIX)
		DIMENSION DELTA(3,3)
		DIMENSION BC(NX),DFI(NX),AFR(NFRX)
		INTEGER CONE(NE,4),KODE(NX)
		REAL INSTANTE_INICIAL,INSTANTE_FINAL
C		Declara a matriz de pivoteamento para o LAPACK.
		INTEGER P(1000)
		COMPLEX*16 WORK(1000)
		COMPLEX   SWORK(1000000)
		DOUBLE PRECISION RWORK(1000)
		INTEGER ITER

*
* NE = N�MERO M�XIMO DE ELEMENTOS DA MALHA (CONTORNO + ENCLOSING)
* NX = DIMENS�O M�XIMA DO SISTEMA DE EQUA��ES		
* NCOX = N�MERO M�XIMO DE PONTOS EXTREMOS DOS ELEMENTOS (NE+2)
* NPIX = N�MERO M�XIMO DE PONTOS INTERNOS
*
* PAR�METROS ESPECIFICADORES DOS N�MEROS DE UNIDADE (ARQUIVO) DE E/S
*
		INP=5
		INQ=6
		IPR=7
		IPS=8
		IPT=9
      
		OPEN(UNIT=INP,FILE='ESOLO240E_-5+5_e10.DAT',STATUS='OLD')
      	OPEN(UNIT=INQ,FILE='ESOLO240D_-5+5_e10.DAT',STATUS='OLD')
      	OPEN(UNIT=IPR,FILE='SSOLO240E_-5+5_e10.DAT',STATUS='UNKNOWN')
      	OPEN(UNIT=IPS,FILE='SSOLO240D_-5+5_e10.DAT',STATUS='UNKNOWN')
C     	OPEN(UNIT=IPT,FILE='GSOLO240_a5b5.DAT',STATUS='UNKNOWN')
*
		CALL CPU_TIME(INSTANTE_INICIAL)
*
* ACIONA ROTINA QUE L� OS DADOS DE ENTRADA PARA O PROBLEMA EST�TICO
*
      
      	CALL INPUTECE (CX,CY,CZ,NE,NCOX,CONE,CXM,CYM,CZM,
     $N,NBE,NP,NPG,GE,RNU,RMU)
*
* ACIONA ROTINA QUE CALCULA AS MATRIZES [H] E [G] DO PROBLEMA EST�TICO
*
      	CALL GHMATECE (CX,CY,CZ,CXM,CYM,CZM,HEST,GEST,NE,NX,NCOX,CONE,
     $N,NBE,NP,NPG,GE,RNU,RMU,DELTA,PI,C1,C2,C3,C4)
*
* ACIONA ROTINA QUE L� OS DADOS DE ENTRADA PARA O PROBLEMA DIN�MICO
*
      	CALL INPUTECD (CX,CY,CZ,CXI,CYI,CZI,KODE,BC,NFR,AFR,NE,NX,
     $NCOX,NPIX,NFRX,CONE,CXM,CYM,CZM,N,NBE,NP,NPG,GE,RNU,RMU,
     $L,FR,DAM,RHO,ZGE,ZCS,ZCP)
	
		NN=3*NBE
		DO 40 J=1,NN
			DFI(J)=BC(J)
 40		CONTINUE
		DO 12 I=1,NFR
			FR=AFR(I)
		WRITE(*,*)'rodando para frequencia ',i
*
* ACIONA ROTINA QUE CALCULA AS MATRIZES [H] E [G] DO PROBLEMA EST�TICO
*
      	CALL GHMATECD (CX,CY,CZ,CXM,CYM,CZM,HEST,GEST,ZH,ZG,ZFI,DFI,ZDFI,
     $KODE,NE,NX,NCOX,CONE,DELTA,PI,N,NBE,NP,NPG,GE,RNU,RMU,
     $L,FR,DAM,RHO,ZGE,ZCS,ZCP,C1,C2,C3,C4)
*
* ACIONA ROTINA QUE RESOLVE O SISTEMA DE EQUA��ES
*
* (ROTINA PRONTA DA BIBLIOTECA DO FORTRAN PS 4.0)	
*
C      	CALL DLSACG(NN,ZH,NX,ZFI,1,ZVETSOL)
		CALL ZCGESV(NN,1,ZH,NN,P,ZFI,NN,VETSOL,NN,WORK,SWORK,RWORK,ITER)
		IF (ITER /= 0) THEN
			PRINT *, "Erro em ZCGESV :-("
		ENDIF
		ZFI=ZVETSOL
*
* ACIONA ROTINA QUE CALCULA OS DESLOCAMENTOS EM PONTOS INTERNOS
*
		CALL INTEREC(ZFI,ZDFI,KODE,CX,CY,CZ,CXI,CYI,CZI,ZDSOL,ZSSOL,NE,NX,
     $NCOX,NPIX,CONE,ZGE,ZCS,ZCP,DELTA,PI,FR,NPG,L,N,NBE,RHO)
*
* ACIONA ROTINA QUE IMPRIME AS VARI�VEIS NO CONTORNO E PONTOS INTERNOS
*
     	 CALL OUTPUTEC(ZFI,ZDFI,ZDSOL,ZSSOL,NPIX,NX,N,NBE,L,FR)
 12	CONTINUE
*
      	CLOSE (INP)          
      	CLOSE (INQ)
      	CLOSE (IPR)          
      	CLOSE (IPS)
C       CLOSE (IPT)
		CALL CPU_TIME(INSTANTE_FINAL)
		WRITE (*,*)'TEMPO DE PROCESSAMENTO: ',
     $INSTANTE_FINAL-INSTANTE_INICIAL,' SEGUNDOS'
		PAUSE
    
		STOP
      	END
*
*#######################################################################
*#######################################################################
*                   FIM DO PROGRAMA PRINCIPAL                          #
*#######################################################################
*#######################################################################
