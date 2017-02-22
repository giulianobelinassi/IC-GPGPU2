************************************************************************
*                                                                      *
*      ESTA SUBROTINA L� E IMPRIME OS DADOS DO ARQUIVO DE ENTRADA      *
*                                                                      *
************************************************************************
*
      SUBROUTINE INPUTECE (CX,CY,CZ,NE,NCOX,CONE,CXM,CYM,CZM,
     $    N,NBE,NP,NPG,GE,RNU,RMU)
*
        IMPLICIT REAL (A-H,O-Y)
        IMPLICIT COMPLEX (Z)
        COMMON INP,INQ,IPR,IPS,IPT
        CHARACTER*75 TITULO
        DIMENSION CX(NCOX),CY(NCOX),CZ(NCOX)
        DIMENSION CXM(NE),CYM(NE),CZM(NE)
        INTEGER CONE(NE,4)
        REAL cone_in(4)
*
        WRITE(IPR,100)
 100    FORMAT(/' ',79('*')/)
*
* LEITURA DO T�TULO DO TRABALHO
*
        READ(INP,'(A)')TITULO
        WRITE(IPR,'(A)')TITULO
*
* LEITURA DO N�MERO DE N�S E DAS PROPRIEDADES DO MATERIAL
*
        READ(INP,*)N,NBE,NP,NPG,GE,RNU 
        WRITE(IPR,200)N,NBE,NP,NPG,GE,RNU     
 200    FORMAT(//'DADOS'//2X,'N�MERO DE ELEMENTOS DA MALHA= ',I4/
     $      2X,'N�MERO DE ELEMENTOS DE CONTORNO= 'I4/
     $      2X,'N�MERO DE PONTOS EXTREMOS DOS ELEMENTOS= 'I4/
     $      2X,'N�MERO DE PONTOS DE GAUSS= ',I3/    
     $      2X,'M�DULO DE CISALHAMENTO= 'D14.7/
     $      2X,'COEFICIENTE DE POISSON= ',D14.7)
*
        RMU=GE 
*
* LEITURA DAS COORDENADAS DOS PONTOS EXTREMOS DOS ELEMENTOS DE CONTORNO
* NAS DIRE��ES 'X','Y' E 'Z'
*
        WRITE(IPR,300)
 300    FORMAT(//2X,'COORDENADAS DOS PONTOS EXTREMOS DOS ELEMENTOS DE CO
     $NTORNO'//4X,'PONTO',10X,'X',18X,'Y',18X,'Z')
        READ(INP,*)(CX(I),CY(I),CZ(I),I=1,NP)
        DO 10 I=1,NP
 10         WRITE(IPR,400)I,CX(I),CY(I),CZ(I)
 400    FORMAT(4X,I4,3(5X,D14.7))
*
* LEITURA DA CONECTIVIDADE ENTRE OS ELEMENTOS
*
        WRITE(IPR,500)
 500    FORMAT(//2X,'CONECTIVIDADE DOS ELEMENTOS DE CONTORNO'//
     $  1X,'ELEMENTO',3X,'EXTR I',3X,'EXTR J',3X,'EXTR K',
     $  3X,'EXTR L',5X,'XMED',8X,'YMED',8X,'ZMED')
        DO 18 I=1,N
            READ(INP,*) (cone_in(J),J=1,4)
            CONE(I, 1:4) = INT(cone_in(1:4))
*
* CALCULA AS COORDENADAS NODAIS E ARMAZENA NOS VETORES XM, YM E ZM
*
            N1=CONE(I,1)
            N2=CONE(I,2)
            N3=CONE(I,3)
            N4=CONE(I,4)
*
* ELEMENTOS QUADRILATERAIS DE QUATRO N�S
*
            CXM(I)=(CX(N1)+CX(N2)+CX(N3)+CX(N4))*0.250
            CYM(I)=(CY(N1)+CY(N2)+CY(N3)+CY(N4))*0.250
            CZM(I)=(CZ(N1)+CZ(N2)+CZ(N3)+CZ(N4))*0.250
            WRITE(IPR,550) I,(CONE(I,J),J=1,4),CXM(I),CYM(I),CZM(I)
 550        FORMAT(I6,4I9,2X,3D12.4)
  18    CONTINUE
*
        RETURN
      END
