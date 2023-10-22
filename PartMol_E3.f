* Program PartMolE3 - Self-Energy Decomposition
* Authors: Augusto C. Moreira
*          Julio C. S. Da Silva
*          Gabriela Balbino

*      Open Overlap Matrix        *
       PROGRAM QUEBRAMOL3
       PARAMETER (NMAX=500000)
       PARAMETER (NLI=500)
       PARAMETER (NLK=7)
       PARAMETER (NLJ=5)
       CHARACTER ARQ3*35,V1*8,V*23,ARQOV*35
       INTEGER FAT,I,J,K,L1,L2,L3(NLI),M,N,O,P,BIT(NLK),W,NU
       INTEGER SUB,OM,IA(NLK),IB(NLK)
       INTEGER D,E1,E2,E3,E4,E5,CAM,E,F,C1,C2
3      FORMAT (1A23,5F10.5)
4      FORMAT (1A8,5D14.6)
6      FORMAT (1I3,5F11.5,1F13.6)

       REAL X(NLK,NLI,NLJ)
       REAL A(NMAX,5),B(NLI,NLI),CO(NLI,NLI,NLJ)
       REAL COL(NLK,NLI,NLI,NLJ)
       REAL COLC(NLK,NLI,NLI,NLJ)
       REAL TSM(NLI,NLJ),MAT(NLI,NLI)
       REAL AV(NLI,NLJ)
       REAL LAP(NLK,NLK,NLI,NLJ)
*       REAL SL(NLK,NLI,NLJ)
       REAL TS(NLI,NLJ)
       

******************************************************************************
******************************************************************************
******************************************************************************
         L1=359291

       ARQOV='Smatriz_anth_cis.DAT'
       ARQ3='SAI_C60_anth_cis.DAT'

         CAM=1
         N=1892
         OM=40

* NUMERO DE PEDAÄOS *****************************************************

         NU=3
         DO W=1,CAM
         BIT(W)=385
         ENDDO

          IA(1)=  1
          IB(1)=  212
          IA(2)=  213
          IB(2)=  1052
          IA(3)=  1053
          IB(3)=  N

******************************************************************************
******************************************************************************
******************************************************************************
         L2=L1+3

*         ARQ(1)='COF00.DAT'
*         ARQ(2)='ANT_Z+18.log'
*         ARQ(3)='ANT_Z+19.log'
*         ARQ(4)='ANT_Z+20.log'

*       NU3,A3,B3 = ELETRODO 1
*       NU2,A2,B2 = ELETRODO 2
*       NU1,A1,B1 = MOLECULA

*         ARQOV='BZ_Pb10_ETIL3_A.out'
         PRINT*,'LENDO OVERLAP'
         OPEN (1,FILE=ARQOV,STATUS='OLD')

         DO J=1,L1
         READ(1,4) V1,A(J,1),A(J,2),A(J,3),A(J,4),A(J,5)
         ENDDO
         CLOSE(1)


                                                           DO W=1,CAM
       PRINT*,'LENDO COEFICIENTES'
       OPEN (2,FILE=ARQOV,STATUS='OLD')

         DO J=1,L2
         READ(2,*)
         ENDDO

          F=(OM/5)
*             DO L=L3(W)+1,(L3(W)+1)+(N+3)*F
         DO D=1,F
         E=5*D-5
       E1=E+1
       E2=E+2
       E3=E+3
       E4=E+4
       E5=E+5
       READ(2,3)V,AV(E1,W),AV(E2,W),AV(E3,W),AV(E4,W),AV(E5,W)
       DO J=1,N
       READ(2,3)V,CO(J,E1,W),CO(J,E2,W),CO(J,E3,W),CO(J,E4,W),CO(J,E5,W)
       ENDDO
       READ(2,*)
       READ(2,*)
         ENDDO
*             ENDDO

       CLOSE(2)

                                                            ENDDO



        !PRINT*,'DIGITE O NOME DO ARQUIVO DE SA÷DA                  '
*        READ(*,*)ARQ3

                                                            
                                                            DO W=1,CAM
       J=0
10     SUB=N+1
       J=J+1
       K=0
       O=J
20            DO M=1,5
              IF ((K+M).GT.J) THEN
                 DO FAT=(K+M),N
                 B(J,FAT)=0
                 ENDDO
                 GOTO 10
              ELSE
                B(J,K+M)=A(O,M)
              ENDIF
              ENDDO
              SUB=SUB-5
              O=O+SUB
              K=K+5
              IF (J.LE.N) THEN
                 GOTO 20
              ENDIF

              DO I=1,N
              DO J=1,I
                MAT(J,I)=B(I,J)
                MAT(I,J)=MAT(J,I)
              ENDDO
              ENDDO
              

*      ZERANDO OS COEFICIENTES DOS SUB ESPAÄOS                           *
*      I - êSIMA FUNÄ«O DE BASE
*      J - êSIMO INTERVALO A SER "ZERADO"
*      K - êSIMO ORBITAL MOLECULAR
           DO K=1,OM
           DO J=1,NU
           DO I=1,N
               COLC(J,I,K,W)=CO(I,K,W)
           ENDDO
           ENDDO
           ENDDO
*       ELETRODO 1                                                       *

               DO K=1,OM
               DO J=1,NU
               C1=IA(J)
               C2=IB(J)
*               PRINT*,C1,C2,J
               DO I=C1,C2
               COLC(J,I,K,W)=0
               ENDDO
               ENDDO
               ENDDO

               DO K=1,OM
               DO J=1,NU
               DO I=1,N
               COL(J,I,K,W)=COLC(J,I,K,W)-CO(I,K,W)
               ENDDO
               ENDDO
               ENDDO

       PRINT*,'                                                 '

       PRINT*,'CSC           CSC           CSC             CSC    '
       
        DO I=1,OM
        DO J=1,NU
        DO L=1,NU
        TSM(I,W)=0
        X(J,I,W)=0
        LAP(J,L,I,W)=0
        TS(I,W)=0
        ENDDO
        ENDDO
        ENDDO

       PRINT*,'                                                   '

              DO I=1,OM
              DO K=1,N
              DO P=1,N
              TSM(I,W)=TSM(I,W)+CO(K,I,W)*MAT(K,P)*CO(P,I,W)
              ENDDO
              ENDDO
              ENDDO

        DO I=1,OM
          DO J=1,NU
          DO L=1,NU
            DO K=1,N
            DO P=1,N
            LAP(J,L,I,W)=LAP(J,L,I,W)+COL(J,K,I,W)*MAT(K,P)*COL(L,P,I,W)
            ENDDO
            ENDDO
*           PRINT*,J,L,I,W,LAP(J,L,I,W)
          ENDDO
          ENDDO
       ENDDO


        PRINT*,'                                                  '
        DO I=1,OM
        DO J=1,NU
        DO L=1,NU
          X(J,I,W)=X(J,I,W)+LAP(J,L,I,W)
        ENDDO
*        PRINT*,I,J,W,' ',X(J,I,W)
        TS(I,W)=TS(I,W)+abs(X(J,I,W))
*        TS(I,W)=TS(I,W)+X(J,I,W)
        ENDDO
           DO J=1,NU
           X(J,I,W)=abs(X(J,I,W))/TS(I,W)
           ENDDO
           TS(I,W)=0
           DO J=1,NU
           TS(I,W)=TS(I,W)+abs(X(J,I,W))
           ENDDO
        ENDDO

        DO I=1,OM
        PRINT*,I,W,' ',TS(I,W)
        ENDDO

                                                               ENDDO

        PRINT*,'                                                   '
        PRINT*,' FIM ----- FIM ----- FIM ----- FIM ------ FIM '
        PRINT*,'                                                   '

       OPEN (4,FILE=ARQ3,STATUS='UNKNOWN')
       WRITE(4,*) 'PROBABILIDADES COM OVERLAP INCLUSAS                 '
       WRITE(4,*) '                                                    '
       WRITE(4,*) 'ORB   CONT1  CONT2  PONTE   MOLLD1  MOLLD2   AU.VAL.'
       WRITE(4,*) '                                                    '

       DO W=1,CAM
       DO I=1,OM
       BIT(W)=BIT(W)+1
       WRITE(4,6)BIT(W),X(1,I,W),X(2,I,W),X(3,I,W),TS(I,W),AV(I,W)
       ENDDO
       WRITE(4,*) '                                                    '
       ENDDO

       DO W=1,CAM
       BIT(W)=245
       DO I=1,OM
       BIT(W)=BIT(W)+1
       WRITE(4,6)BIT(W),X(1,I,W)*X(3,I,W),TS(I,W),AV(I,W)
       ENDDO
       WRITE(4,*) '                                                    '
       ENDDO

       DO W=1,CAM
       BIT(W)=245
       DO I=1,OM
       BIT(W)=BIT(W)+1
       WRITE(4,6)BIT(W),LAP(1,2,I,W),LAP(2,3,I,W),LAP(3,2,I,W),AV(I,W)
       ENDDO
       WRITE(4,*) '                                                    '
       ENDDO

       DO W=1,CAM
       BIT(W)=245
       DO I=1,OM
       BIT(W)=BIT(W)+1
       WRITE(4,6)BIT(W),LAP(1,2,I,W)*AV(I,W),LAP(2,3,I,W)*AV(I,W)
       ENDDO
       WRITE(4,*) '                                                    '
       ENDDO

       DO I=1,N
       WRITE(4,*) MAT(I,I)
       ENDDO
       READ(*,*)
       END
