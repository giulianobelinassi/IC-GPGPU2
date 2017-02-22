************************************************************************
*                                                                      *
*    SUBROTINA QUE CALCULA OS COEFICIENTES DA DIAGONAL DE [G] POR      *
*																	   *	
*  INTEGRAÇÃO NUMÉRICA. O PONTO DE COLOCAÇÃO ESTÁ NO ELEMENTO A SER	   *
*                                                                      *
*                       INTEGRADO                                      *
*                                                                      *
************************************************************************
*
      SUBROUTINE SINGGE(HELEM,GELEM,CXM,CYM,CZM,ETA,CX,CY,CZ,
     $  N1,N2,N3,N4,NCOX,N,NP,NPG,GE,RNU,RMU,C1,C2,C3,C4,DELTA)
*
        IMPLICIT REAL (A-H,O-Y)
        IMPLICIT COMPLEX (Z)
        COMMON INP,INQ,IPR,IPS,IPT
        DIMENSION DELTA(3,3)
        DIMENSION HELEM(3,3),GELEM(3,3),HP(3,3),GP(3,3)
        DIMENSION CX(NCOX),CY(NCOX),CZ(NCOX)
        DIMENSION CO(4,3),ETA(3)
*
* INTEGRAÇÃO NUMÉRICA SOBRE TRIÂNGULOS INTERNOS COLAPSANDO OS CANTOS 3-4
* DE UM ELEMENTO QUADRILATERAL
*
        CO(3,1)=CXM
        CO(3,2)=CYM
        CO(3,3)=CZM
        CO(4,1)=CXM
        CO(4,2)=CYM
        CO(4,3)=CZM
*


        HELEM = 0.D0
        GELEM = 0.D0
        
        DO 500 IT=1,4

*           SWITCH IT
            GO TO (20,30,40,50),IT
  20            CO(1,1)=CX(N1)
                CO(1,2)=CY(N1)
                CO(1,3)=CZ(N1)
                CO(2,1)=CX(N2)
                CO(2,2)=CY(N2)
                CO(2,3)=CZ(N2)
            GO TO 60
  30            CO(1,1)=CX(N2)
                CO(1,2)=CY(N2)
                CO(1,3)=CZ(N2)
                CO(2,1)=CX(N3)
                CO(2,2)=CY(N3)
                CO(2,3)=CZ(N3)
            GO TO 60
  40            CO(1,1)=CX(N3)
                CO(1,2)=CY(N3)
                CO(1,3)=CZ(N3)
                CO(2,1)=CX(N4)
                CO(2,2)=CY(N4)
                CO(2,3)=CZ(N4)
            GO TO 60
  50            CO(1,1)=CX(N4)
                CO(1,2)=CY(N4)
                CO(1,3)=CZ(N4)
                CO(2,1)=CX(N1)
                CO(2,2)=CY(N1)
                CO(2,3)=CZ(N1)
   60       CONTINUE
*
* INTEGRAÇÃO SOBRE O ELEMENTO QUADRILATERAL COLAPSADO
*
            CALL NONSINGE(HP,GP,CO,CXM,CYM,CZM,ETA,N,NP,NPG,GE,RNU,RMU,
     $          C1,C2,C3,C4,DELTA)
*
            
            GELEM = GELEM + GP
            HELEM = HELEM + HP
*
  500   CONTINUE
*
        RETURN
      END
