        subroutine pdf_tr_size_setup(s1,s2,PDFs1,PDFs2,ds1,ds2,Ns1,Ns2) 
c
        save
c
        INCLUDE 'parameters.include'
c
        integer I
        integer Ns1,Ns2

        real s1(maxtrs1), s2(maxtrs2),ds1(maxtrs1), ds2(maxtrs2)
        real PDFs1(maxtrs1),PDFs2(maxtrs2),TRUNK_HGHT_FUNC

        INTEGER I_PDF_LEAF_SIZE, I_PDF_BR_1_SIZE, I_PDF_BR_2_SIZE,
     &          I_PDF_BR_3_SIZE, I_PDF_BR_4_SIZE, I_PDF_BR_5_SIZE, 
     &          I_PDF_BR_6_SIZE, I_PDF_NDL_SIZE, I_PDF_TRUNK_SIZE

        COMMON /I_PDF_SIZE/ I_PDF_LEAF_SIZE,I_PDF_BR_1_SIZE,
     &               I_PDF_BR_2_SIZE, I_PDF_BR_3_SIZE, I_PDF_BR_4_SIZE, 
     &               I_PDF_BR_5_SIZE, I_PDF_BR_6_SIZE, I_PDF_NDL_SIZE,
     &               I_PDF_TRUNK_SIZE


c
        REAL MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
        REAL DENSITY, CROWN_HGHT, TRUNK_HGHT
        REAL TRUNK_Dsig, TRUNK_Hsig 

        COMMON /R_TRUNK/ MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
        COMMON /R_CANOPY/ DENSITY, CROWN_HGHT, TRUNK_HGHT
        COMMON /R_TRUNK_SIG/ TRUNK_Dsig, TRUNK_Hsig 
C
        REAL TRUNK_DIAM_HIST_X(MAXTRDHIST),TRUNK_DIAM_HIST_Y(MAXTRDHIST)
        REAL TR_DIAM_HIST_DX(MAXTRDHIST)
        INTEGER I_TR_DIAM_HIST
C
        COMMON /TRUNK_HIST/ TRUNK_DIAM_HIST_X, TRUNK_DIAM_HIST_Y,
     &                      TR_DIAM_HIST_DX, I_TR_DIAM_HIST
c
        IF(I_PDF_TRUNK_SIZE.EQ.0)THEN
            Ns1 = 1
            Ns2 = 1
            s1(1) = TRUNK_HGHT
            s2(1) = TRUNK_DIAM
            ds1(1) = 1.0
            ds2(1) = 1.0
            PDFs1(1) = 1.0
            PDFs2(1) = 1.0
        ELSE IF(I_PDF_TRUNK_SIZE.EQ.1)THEN
            Ns1 = 1
            Ns2 = I_TR_DIAM_HIST
            DO 10 I=1,Ns2
C
                s2(I) = TRUNK_DIAM_HIST_X(I)
                PDFs2(I) = TRUNK_DIAM_HIST_Y(I)
                ds2(I) = TR_DIAM_HIST_DX(I)
C
                PDFs1(I) = 1.0
                ds1(I) = 1.0
                s1(I) = TRUNK_HGHT_FUNC(s2(I))
c
10          CONTINUE
C
C        ELSE IF(I_PDF_TRUNK_SIZE.EQ.2)THEN
CC
C            Ns1 = INT(TRUNK_Hsig/TRUNK_HGHT)
C            Ns2 = INT(TRUNK_Dsig/TRUNK_DIAM)
C            do 30 i=1,Ns1
C                x = 
C                do 20 j=1,Ns2
C                    y = 
C                    f(i,j) = bivnorm(x,y,TRUNK_HGHT,TRUNK_DIAM,
C     &                           TRUNK_Hsig,TRUNK_Dsig,RHOtrHD)
C20              continue
C30          continue
C
C
C        ELSE IF(I_PDF_TRUNK_SIZE.EQ.3)THEN
C            GAUSS_DIST(X,TRUNK_HGHT,TRUNK_Hsig)
C        ELSE IF(I_PDF_TRUNK_SIZE.EQ.4)THEN
C            GAUSS_DIST(X,TRUNK_DIAM,TRUNK_Dsig)
        ELSE
            CALL WRITE_ERROR(115)
            STOP
        ENDIF
c
        return
        end
