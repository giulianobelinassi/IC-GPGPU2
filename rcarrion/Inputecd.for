************************************************************************
*                                                                      *
*      ESTA SUBROTINA L� E IMPRIME OS DADOS DO ARQUIVO DE ENTRADA      *
*                                                                      *
************************************************************************
*
      SUBROUTINE INPUTECD (CX,CY,CZ,CXI,CYI,CZI,KODE,BC,NFR,AFR,NE,NX,
     $      NCOX,NPIX,NFRX,CONE,CXM,CYM,CZM,N,NBE,NP,NPG,GE,RNU,RMU,
     $      L,FR,DAM,RHO,ZGE,ZCS,ZCP)
*
        IMPLICIT REAL*8 (A-H,O-Y)
        IMPLICIT COMPLEX*16 (Z)
        COMMON INP,INQ,IPR,IPS,IPT
        CHARACTER*75 TITULO
        DIMENSION CX(NCOX),CY(NCOX),CZ(NCOX)
        DIMENSION CXM(NE),CYM(NE),CZM(NE)
        DIMENSION CXI(NPIX),CYI(NPIX),CZI(NPIX)
        DIMENSION BC(NX)
        DIMENSION AFR(NFRX)
C       DIMENSION DELTA(3,3)
        INTEGER CONE(NE,4),KODE(NX)
*
        WRITE(IPS,100)
 100    FORMAT(/' ',79('*')/)
*
* LEITURA DO T�TULO DO TRABALHO
*
        READ(INQ,'(A)')TITULO
        WRITE(IPS,'(A)')TITULO
*
* LEITURA DO N�MERO DE FREQ��NCIAS E DAS PROPRIEDADES DO MATERIAL
*
        READ(INQ,*)L,NFR,DAM,RHO 
        WRITE(IPS,200)N,NBE,NP,L,NPG,NFR,GE,RNU,DAM,RHO     
 200    FORMAT(//'DADOS'//2X,'N�MERO DE ELEMENTOS DA MALHA= ',I4/
     $      2X,'N�MERO DE ELEMENTOS DE CONTORNO= 'I4/
     $      2X,'N�MERO DE PONTOS EXTREMOS DOS ELEMENTOS= 'I4/
     $      2X,'N�MERO DE PONTOS INTERNOS= ',I3/
     $      2X,'N�MERO DE PONTOS DE GAUSS= ',I3/    
     $      2X,'N�MERO DE FREQ��NCIAS= ',I3/    
     $      2X,'M�DULO DE CISALHAMENTO (REAL)= 'D14.7/
     $      2X,'COEFICIENTE DE POISSON= ',D14.7/
     $      2X,'COEFICIENTE DE AMORTECIMENTO= ',D14.7/
     $      2X,'DENSIDADE DE MASSA= ',D14.7)
*
        ZGE=DCMPLX(GE,GE*2.D0*DAM)
        ZCS=CDSQRT(ZGE/RHO)
        ZCP=ZCS*DSQRT((2.D0-2.D0*RNU)/(1.D0-2.D0*RNU))
        WRITE(IPS,230)ZGE,ZCS,ZCP     
 230        FORMAT(//'DADOS DERIVADOS'//
     $          2X,'M�DULO DE CISALHAMENTO (COMPLEXO)= '2(D14.7)/
     $          2X,'VELOCIDADE DA ONDA-S (COMPLEXO)= ',2(D14.7)/
     $          2X,'VELOCIDADE DA ONDA-P (COMPLEXO)= ',2(D14.7))

*
* LEITURA DAS FREQ��NCIAS
*
        READ(INQ,*) (AFR(I),I=1,NFR)
        WRITE(IPS,250) (AFR(I),I=1,NFR)
 250    FORMAT(/2X,'FREQUENCIAS'//5(2X,D14.7))
*
* IMPRIME AS COORDENADAS DOS PONTOS EXTREMOS DOS ELEMENTOS DE CONTORNO
* NAS DIRE��ES 'X','Y' E 'Z'
*
        WRITE(IPS,300)
 300      FORMAT(//2X,'COORDENADAS DOS PONTOS EXTREMOS DOS ELEMENTOS DA 
     $MALHA'//4X,'PONTO',10X,'X',18X,'Y',18X,'Z')
C      READ(INP,*)(CX(I),CY(I),CZ(I),I=1,NP)
        DO 10 I=1,NP
 10         WRITE(IPS,400)I,CX(I),CY(I),CZ(I)
 400            FORMAT(4X,I4,3(5X,D14.7))
*
* IMPRIME A CONECTIVIDADE ENTRE OS ELEMENTOS
*
        WRITE(IPS,500)
 500        FORMAT(//2X,'CONECTIVIDADE DOS ELEMENTOS DA MALHA'//
     $      1X,'ELEMENTO',3X,'EXTR I',3X,'EXTR J',3X,'EXTR K'
     $      ,3X,'EXTR L',5X,'XMED',8X,'YMED',8X,'ZMED')
        DO 18 I=1,N
C       READ(INP,*) (CONE(I,J),J=1,4)
*
* CALCULA AS COORDENADAS NODAIS E ARMAZENA NOS VETORES XM, YM E ZM
*
C           N1=CONE(I,1)
C           N2=CONE(I,2)
C           N3=CONE(I,3)
C           N4=CONE(I,4)
*
* ELEMENTOS QUADRILATERAIS DE QUATRO N�S
*
C           CXM(I)=(CX(N1)+CX(N2)+CX(N3)+CX(N4))*0.25D0
C           CYM(I)=(CY(N1)+CY(N2)+CY(N3)+CY(N4))*0.25D0
C           CZM(I)=(CZ(N1)+CZ(N2)+CZ(N3)+CZ(N4))*0.25D0
            WRITE(IPS,550) I,(CONE(I,J),J=1,4),CXM(I),CYM(I),CZM(I)
 550        FORMAT(I6,4I9,2X,3D12.4)
  18    CONTINUE

*
* LEITURA DAS CONDI��ES DE CONTORNO ARMAZENADAS NO VETOR BC(I)
* SE KODE(I)=0; BC(I) CONT�M O VALOR DO DESLOCAMENTO
* SE KODE(I)=1; BC(I) CONT�M O VALOR DA "TRACTION"
*
        WRITE(IPS,600)
 600  FORMAT(//2X,'CONDI��ES DE CONTORNO'//8X,3('VALOR PRESCRITO',10X)
     $/1X,'ELEM',6X,'DIRE��O X',4X,'C�DIGO',6X,'DIRE��O Y',4X,'C�DIGO',
     $6X,'DIRE��O Z',4X,'C�DIGO')

        DO 20 I=1,NBE
            READ(INQ,*) (KODE(3*(I-1)+J),BC(3*(I-1)+J),J=1,3)
  20        WRITE(IPS,650) I,(BC(3*(I-1)+J),KODE(3*(I-1)+J),J=1,3)
 650        FORMAT(I4,3(4X,D10.3,6X,I3,2X))
*
* LEITURA DAS COORDENADAS DOS PONTOS INTERNOS
*
        IF(L == 0) THEN
            RETURN
        ENDIF

        WRITE(IPS,700)
 700    FORMAT(//2X,'COORDENADAS DOS PONTOS INTERNOS ',
     $  //4X,'PONTO',10X,'X',18X,'Y',17X,' Z')
        READ(INQ,*)(CXI(I),CYI(I),CZI(I),I=1,L) 
        WRITE(IPS,750) (I,CXI(I),CYI(I),CZI(I),I=1,L)
 750    FORMAT(4X,I3,3D19.7)
*
        RETURN
      END
