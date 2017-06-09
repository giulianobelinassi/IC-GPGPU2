************************************************************************
************************************************************************
**                                                                    **
**     ESTE PROGRAMA RESOLVE PROBLEMAS VISCOELASTICOS, DINÂMICOS,     **
**                                                                    **
**   TRIDIMENSIONAIS ATRAVÉS DO MÉTODO DOS ELEMENTOS DE CONTORNO NUM  **
**                                                                    **
**                      SEMI ESPAÇO INFINITO                          **
**                                                                    **
************************************************************************
**                                                                    **
** PROGRAMA "MECECCS" = (MÉTODO DOS ELEMENTOS DE CONTORNO USANDO      **
**                                                                    **
**                     ELEMENTOS CONSTANTES COM SINGULARIDADE)        **
**                                                                    **
**                              DINÂMICO                              **
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
C       PARAMETER (NE=960,NX=3*NE,NCOX=962,NPIX=10,NFRX=7)
        PARAMETER (NE=2200,NX=3*NE,NCOX=2200,NPIX=10,NFRX=7)
        COMMON INP,INQ,IPR,IPS,IPT
        REAL, DIMENSION(:), ALLOCATABLE :: CX, CY, CZ, CXM, CYM, CZM
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: CONE, KODE
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
* NE = NÚMERO MÁXIMO DE ELEMENTOS DA MALHA (CONTORNO + ENCLOSING)
* NX = DIMENSÃO MÁXIMA DO SISTEMA DE EQUAÇÕES       
* NCOX = NÚMERO MÁXIMO DE PONTOS EXTREMOS DOS ELEMENTOS (NE+2)
* NPIX = NÚMERO MÁXIMO DE PONTOS INTERNOS
*
* PARÂMETROS ESPECIFICADORES DOS NÚMEROS DE UNIDADE (ARQUIVO) DE E/S
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

C ACIONA ROTINA QUE REMOVE AS RESTRIÇÕES IMPOSTAS PELO SISTEMA
C OPERACIONAL SOBRE A STACK
        CALL request_unlimited_stack(ITER)
        IF (ITER /= 0) THEN
            PRINT*, "Falha ao pedir uma stack de tamanho ilimitado! :-("
            STOP
        ENDIF

*
*
* ACIONA ROTINA QUE LÊ OS DADOS DE ENTRADA PARA O PROBLEMA ESTÁTICO
*
      
        CALL INPUTECE (CX,CY,CZ,NE,NCOX,CONE,CXM,CYM,CZM,
     $      N,NBE,NP,NPG,GE,RNU,RMU)


! Aciona a rotina que calcula as normas originalmente usadas em
! Ghmatece.for Ghmatecd.for e Interec.for. Isto evita calculos
! reduntantes.
        ETAS = NORMVEC(CONE, CX, CY, CZ, NX, NE, NCOX, N)
*
* ACIONA ROTINA QUE CALCULA AS MATRIZES [H] E [G] DO PROBLEMA ESTÁTICO
*
        CALL GHMATECE (CX,CY,CZ,CXM,CYM,CZM,HEST,GEST,NE,NX,NCOX,CONE,
     $      N,NBE,NP,NPG,GE,RNU,RMU,DELTA,PI,C1,C2,C3,C4,ETAS)
*
* ACIONA ROTINA QUE LÊ OS DADOS DE ENTRADA PARA O PROBLEMA DINÂMICO
*
        CALL INPUTECD (CX,CY,CZ,CXI,CYI,CZI,KODE,BC,NFR,AFR,NE,NX,
     $     NCOX,NPIX,NFRX,CONE,CXM,CYM,CZM,N,NBE,NP,NPG,GE,RNU,RMU,
     $     L,FR,DAM,RHO,ZGE,ZCS,ZCP)
       
        
        NN=3*NBE
        ALLOCATE(DFI(NN), STAT = ITER)
        IF (ITER == 0) THEN
            PRINT*, "Memória insuficiente"
        ENDIF

        DFI = BC
        
        DO 12 I=1,NFR
            FR=AFR(I)
            PRINT *, 'rodando para frequencia ', i

*
* ACIONA ROTINA QUE CALCULA AS MATRIZES [H] E [G] DO PROBLEMA ESTÁTICO
*
            CALL GHMATECD (CX,CY,CZ,CXM,CYM,CZM,HEST,GEST,ZH,ZG,ZFI,DFI,
     $          ZDFI,KODE,NE,NX,NCOX,CONE,DELTA,PI,N,NBE,NP,NPG,GE,RNU,
     $          RMU,L,FR,DAM,RHO,ZGE,ZCS,ZCP,C1,C2,C3,C4,ETAS
     $      )
*
* ACIONA ROTINA QUE RESOLVE O SISTEMA DE EQUAÇÕES (LAPACK)
*
            ALLOCATE(PIV(NN), STAT = stats1)
            IF (stats1 == 0) THEN
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
* ACIONA ROTINA QUE IMPRIME AS VARIÁVEIS NO CONTORNO E PONTOS INTERNOS
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
