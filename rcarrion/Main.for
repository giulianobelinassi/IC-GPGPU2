************************************************************************
************************************************************************
**                                                                    **
**     ESTE PROGRAMA RESOLVE PROBLEMAS VISCOELASTICOS, DIN�MICOS,     **
**                                                                    **
**   TRIDIMENSIONAIS ATRAV�S DO M�TODO DOS ELEMENTOS DE CONTORNO NUM  **
**                                                                    **
**                      SEMI ESPA�O INFINITO                          **
**                                                                    **
************************************************************************
**                                                                    **
** PROGRAMA "MECECCS" = (M�TODO DOS ELEMENTOS DE CONTORNO USANDO      **
**                                                                    **
**                     ELEMENTOS CONSTANTES COM SINGULARIDADE)        **
**                                                                    **
**                              DIN�MICO                              **
**                                                                    **
** RONALDO CARRION / EUCLIDES DE MESQUITA NETO                        **
**                                                                    **
** DATA: 18/07/01                                                     **
**                                                                    **
** MODIFICADO POR: GIULIANO AUGUSTO FAULIN BELINASSI                  ** 
** DATA: 12/2016                                                      **
************************************************************************
************************************************************************
*
      PROGRAM MAIN
*
        IMPLICIT REAL (A-H,O-Y)
        IMPLICIT COMPLEX (Z)
        INTERFACE 
          SUBROUTINE INPUTECE(CX,CY,CZ,NE,NCOX,CONE,CXM,CYM,CZM,
     $        N,NBE,NP,NPG,GE,RNU,RMU)
            IMPLICIT REAL(A-H, O-Y)
            IMPLICIT COMPLEX(Z)
            REAL, DIMENSION(:), INTENT(OUT), ALLOCATABLE :: CX, CY, CZ
            REAL, DIMENSION(:), INTENT(OUT), ALLOCATABLE :: CXM,CYM, CZM
            INTEGER, DIMENSION(:,:), INTENT(OUT), ALLOCATABLE :: CONE
          END
        END INTERFACE 

        INTERFACE
          FUNCTION NORMVEC(cone, cx, cy, cz, n, np) RESULT(ETAS)
            IMPLICIT NONE
            INTEGER, DIMENSION(n,4), INTENT(IN) :: CONE
            REAL, DIMENSION(np), INTENT(IN) :: CX, CY, CZ
            INTEGER, INTENT(IN) :: n, np
            REAL, DIMENSION(:,:), ALLOCATABLE:: ETAS
          END
        END INTERFACE 

        INTERFACE
          SUBROUTINE GHMATECE(CX,CY,CZ,CXM,CYM,CZM,HEST,GEST,NE,NX,NCOX,
     $        CONE,N,NBE,NP,NPG,GE,RNU,RMU,DELTA,PI,C1,C2,C3,C4,ETAS)
            IMPLICIT REAL (A-H,O-Y)
            IMPLICIT COMPLEX (Z)
            COMMON INP,INQ,IPR,IPS,IPT
            REAL, DIMENSION(NP), INTENT(IN) :: CX, CY, CZ
            REAL, DIMENSION(N) , INTENT(IN) :: CXM, CYM, CZM
            INTEGER, DIMENSION(N, 4), INTENT(IN) :: CONE
            REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: HEST, GEST
            REAL :: DELTA(3,3), ETAS(3, n)
          END
        END INTERFACE

        INTERFACE
          SUBROUTINE INPUTECD (CX,CY,CZ,CXI,CYI,CZI,KODE,BC,NFR,AFR,NE,
     $      NX,NCOX,NPIX,NFRX,CONE,CXM,CYM,CZM,N,NBE,NP,NPG,GE,RNU,RMU,
     $      L,FR,DAM,RHO,ZGE,ZCS,ZCP)
              
            IMPLICIT REAL (A-H,O-Y)
            IMPLICIT COMPLEX (Z)
            REAL, DIMENSION(NP), INTENT(IN) :: CX, CY, CZ
            REAL, DIMENSION(N) , INTENT(IN) :: CXM, CYM, CZM
            REAL, DIMENSION(:) , INTENT(OUT), ALLOCATABLE:: CXI, CYI
            REAL, DIMENSION(:) , INTENT(OUT), ALLOCATABLE:: CZI, BC
            REAL, DIMENSION(:) , INTENT(OUT), ALLOCATABLE:: AFR
            INTEGER, DIMENSION(:), INTENT(OUT), ALLOCATABLE:: KODE
            INTEGER, DIMENSION(N, 4), INTENT(IN) :: CONE
          END
        END INTERFACE

        INTERFACE
          SUBROUTINE GHMATECD (CX,CY,CZ,CXM,CYM,CZM,HEST,GEST,ZH,ZG,ZFI,
     $      DFI,ZDFI,KODE,NE,NX,NCOX,CONE,DELTA,PI,N,NBE,NP,NPG,GE,RNU,
     $      RMU,L,FR,DAM,RHO,ZGE,ZCS,ZCP,C1,C2,C3,C4,ETAS)
          
            IMPLICIT NONE
            REAL, DIMENSION(NP), INTENT(IN) :: CX, CY, CZ
            REAL, DIMENSION(N), INTENT(IN) :: CXM, CYM, CZM
            REAL, DIMENSION(3*NBE, 3*N), INTENT(IN) :: HEST, GEST

            COMPLEX, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ZH
            COMPLEX, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ZG
            COMPLEX, DIMENSION(:), INTENT(OUT), ALLOCATABLE:: ZFI
            REAL, INTENT(IN) :: DFI(3*NBE)
            COMPLEX, INTENT(OUT), ALLOCATABLE :: ZDFI(:)
            INTEGER, INTENT(IN) :: KODE(3*NBE),NE,NX,NCOX,CONE(N,4)
            REAL, INTENT(IN) :: DELTA(3,3),PI
            INTEGER, INTENT(IN) :: N,NBE,NP,NPG, L
            REAL, INTENT(IN) :: GE,RNU,RMU,FR,DAM,RHO
            COMPLEX,   INTENT(IN) :: ZGE,ZCS,ZCP
            REAL, INTENT(IN) :: C1,C2,C3,C4
            REAL, INTENT(IN) :: ETAS(3,N)
          END
        END INTERFACE

        INTERFACE
          SUBROUTINE INTEREC(ZFI,ZDFI,KODE,CX,CY,CZ,CXI,CYI,CZI,ZDSOL,
     $      ZSSOL,NE,NX,NCOX,NPIX,CONE,ZGE,ZCS,ZCP,DELTA,PI,FR,NPG,L,N,
     $      NBE,RHO,ETAS)
            
            IMPLICIT REAL (A-H,O-Y)
            IMPLICIT COMPLEX (Z)
            REAL, DIMENSION(:), INTENT(IN) :: CX, CY, CZ
            REAL, DIMENSION(L),  INTENT(IN) :: CXI, CYI, CZI
            COMPLEX, DIMENSION(3,3) :: ZHELEM, ZGELEM 
            REAL, INTENT(IN) :: DELTA(3,3)
            COMPLEX, DIMENSION(3*NBE), INTENT(INOUT) :: ZDFI, ZFI
            COMPLEX, INTENT(OUT), ALLOCATABLE :: ZDSOL(:), ZSSOL(:)
            INTEGER, INTENT(IN) ::  CONE(N,4),KODE(3*NBE)
            DIMENSION ZD(3,3,3),ZS(3,3,3)
            REAL, INTENT(IN) :: ETAS(3,NX)
          END 
        END INTERFACE
       
        INTERFACE
          SUBROUTINE OUTPUTEC(ZFI,ZDFI,ZDSOL,ZSSOL,NPIX,NX,N,NBE,L,FR)
            IMPLICIT REAL (A-H,O-Y)
            IMPLICIT COMPLEX (Z)
            COMMON INP,INQ,IPR,IPS,IPT
            DIMENSION ZFI(3*NBE),ZDFI(3*NBE),ZDSOL(3*L),ZSSOL(9*L)
          END
        END INTERFACE

C       PARAMETER (NE=960,NX=3*NE,NCOX=962,NPIX=10,NFRX=7)
        PARAMETER (NE=2200,NX=3*NE,NCOX=2200,NPIX=10,NFRX=7)
        COMMON INP,INQ,IPR,IPS,IPT
        REAL, DIMENSION(:), ALLOCATABLE :: CX, CY, CZ, CXM, CYM, CZM
        INTEGER, ALLOCATABLE :: CONE(:,:), KODE(:)
        REAL, DIMENSION(:), ALLOCATABLE :: CXI,CYI, CZI
        REAL, DIMENSION(:), ALLOCATABLE :: BC, AFR, DFI
        REAL, DIMENSION(:,:), ALLOCATABLE :: HEST, GEST
        COMPLEX, DIMENSION(:,:), ALLOCATABLE :: ZH, ZG
        COMPLEX, DIMENSION(:), ALLOCATABLE :: ZDFI, ZFI, ZDSOL, ZSSOL
     

        DIMENSION DELTA(3,3)
        INTEGER, DIMENSION(:), ALLOCATABLE :: PIV
        REAL, DIMENSION(:,:), ALLOCATABLE :: ETAS
        INTEGER stats1
        !REAL, DIMENSION(3,NX) :: ETAS
*
* NE = N�MERO M�XIMO DE ELEMENTOS DA MALHA (CONTORNO + ENCLOSING)
* NX = DIMENS�O M�XIMA DO SISTEMA DE EQUA��ES       
* NCOX = N�MERO M�XIMO DE PONTOS EXTREMOS DOS ELEMENTOS (NE+2)
* NPIX = N�MERO M�XIMO DE PONTOS INTERNOS
*
* PAR�METROS ESPECIFICADORES DOS N�MEROS DE UNIDADE (ARQUIVO) DE E/S
*
        INP=100
        INQ=101
        IPR=102
        IPS=103
        IPT=104
*
  
C       OPEN(UNIT=INP,FILE='ESOLO240E_-5+5.DAT',STATUS='OLD')
C       OPEN(UNIT=INQ,FILE='ESOLO240D_-5+5.DAT',STATUS='OLD')
C       OPEN(UNIT=IPR,FILE='SSOLO240E_-5+5.DAT',STATUS='UNKNOWN')
C       OPEN(UNIT=IPS,FILE='SSOLO240D_-5+5.DAT',STATUS='UNKNOWN')

        OPEN(UNIT=INP,FILE='ESOLO2160E_-5+5.DAT',STATUS='OLD')
        OPEN(UNIT=INQ,FILE='ESOLO2160D_-5+5.DAT',STATUS='OLD')
        OPEN(UNIT=IPR,FILE='SSOLO2160E_-5+5.DAT',STATUS='UNKNOWN')
        OPEN(UNIT=IPS,FILE='SSOLO2160D_-5+5.DAT',STATUS='UNKNOWN')

C       OPEN(UNIT=IPT,FILE='GSOLO240_a5b5.DAT',STATUS='UNKNOWN')

C ACIONA ROTINA QUE REMOVE AS RESTRI��ES IMPOSTAS PELO SISTEMA
C OPERACIONAL SOBRE A STACK
!        CALL request_unlimited_stack(ITER)
!        IF (ITER /= 0) THEN
!            PRINT*, "Falha ao pedir uma stack de tamanho ilimitado! :-("
!            STOP
!        ENDIF

*
*
* ACIONA ROTINA QUE L� OS DADOS DE ENTRADA PARA O PROBLEMA EST�TICO
*
      
        CALL INPUTECE (CX,CY,CZ,NE,NCOX,CONE,CXM,CYM,CZM,
     $      N,NBE,NP,NPG,GE,RNU,RMU)


! Aciona a rotina que calcula as normas originalmente usadas em
! Ghmatece.for Ghmatecd.for e Interec.for. Isto evita calculos
! reduntantes.
        ETAS = NORMVEC(CONE, CX, CY, CZ, N, NP)
*
* ACIONA ROTINA QUE CALCULA AS MATRIZES [H] E [G] DO PROBLEMA EST�TICO
*
        CALL GHMATECE (CX,CY,CZ,CXM,CYM,CZM,HEST,GEST,NE,NX,NCOX,CONE,
     $      N,NBE,NP,NPG,GE,RNU,RMU,DELTA,PI,C1,C2,C3,C4,ETAS)
*
* ACIONA ROTINA QUE L� OS DADOS DE ENTRADA PARA O PROBLEMA DIN�MICO
*
        CALL INPUTECD (CX,CY,CZ,CXI,CYI,CZI,KODE,BC,NFR,AFR,NE,NX,
     $     NCOX,NPIX,NFRX,CONE,CXM,CYM,CZM,N,NBE,NP,NPG,GE,RNU,RMU,
     $     L,FR,DAM,RHO,ZGE,ZCS,ZCP)
       
        
        NN=3*NBE
        ALLOCATE(DFI(NN), STAT = ITER)
        IF (ITER /= 0) THEN
            PRINT*, "Mem�ria insuficiente"
        ENDIF

        DFI = BC
        
        DO 12 I=1,NFR
            FR=AFR(I)
            PRINT *, 'rodando para frequencia ', i

*
* ACIONA ROTINA QUE CALCULA AS MATRIZES [H] E [G] DO PROBLEMA EST�TICO
*
            CALL GHMATECD (CX,CY,CZ,CXM,CYM,CZM,HEST,GEST,ZH,ZG,ZFI,DFI,
     $          ZDFI,KODE,NE,NX,NCOX,CONE,DELTA,PI,N,NBE,NP,NPG,GE,RNU,
     $          RMU,L,FR,DAM,RHO,ZGE,ZCS,ZCP,C1,C2,C3,C4,ETAS
     $      )
*
* ACIONA ROTINA QUE RESOLVE O SISTEMA DE EQUA��ES (LAPACK)
*
            ALLOCATE(PIV(NN), STAT = stats1)
            IF (stats1 /= 0) THEN
                PRINT*, "MEMORIA INSUFICIENTE"
                STOP
            ENDIF

            CALL CGESV(NN,1,ZH,NN,PIV,ZFI,NN,ITER)
            IF (ITER < 0) THEN
                PRINT *, "Erro em ZGESV :-("
            ELSE IF (ITER > 0) THEN
                PRINT *, "Matriz Singular :-|"
            ENDIF

*
* ACIONA ROTINA QUE CALCULA OS DESLOCAMENTOS EM PONTOS INTERNOS
*

            CALL INTEREC(ZFI,ZDFI,KODE,CX,CY,CZ,CXI,CYI,CZI,ZDSOL,ZSSOL,
     $          NE,NX,NCOX,NPIX,CONE,ZGE,ZCS,ZCP,DELTA,PI,FR,NPG,L,N,
     $          NBE,RHO,ETAS
     $      )
     
*
* ACIONA ROTINA QUE IMPRIME AS VARI�VEIS NO CONTORNO E PONTOS INTERNOS
*
            CALL OUTPUTEC(ZFI,ZDFI,ZDSOL,ZSSOL,NPIX,NX,N,NBE,L,FR)
 12     CONTINUE
*
*
            DEALLOCATE(CX)
            DEALLOCATE(CY)
            DEALLOCATE(CZ)

            DEALLOCATE(CXM)
            DEALLOCATE(CYM)
            DEALLOCATE(CZM)
           
            DEALLOCATE(CXI)
            DEALLOCATE(CYI)
            DEALLOCATE(CZI)

            DEALLOCATE(BC)
            DEALLOCATE(KODE)
            DEALLOCATE(AFR)

            DEALLOCATE(CONE) 
            
            DEALLOCATE(ETAS)

            DEALLOCATE(HEST)
            DEALLOCATE(GEST)
    
            DEALLOCATE(ZH)
            DEALLOCATE(ZG)

            DEALLOCATE(ZDFI)
            DEALLOCATE(ZFI)

            DEALLOCATE(DFI)
            DEALLOCATE(ZDSOL)
            DEALLOCATE(ZSSOL)
            DEALLOCATE(PIV)


            CLOSE (INP)          
            CLOSE (INQ)
            CLOSE (IPR)          
            CLOSE (IPS)
C           CLOSE (IPT)
            STOP
      END
*
*#######################################################################
*#######################################################################
*                   FIM DO PROGRAMA PRINCIPAL                          #
*#######################################################################
*#######################################################################
