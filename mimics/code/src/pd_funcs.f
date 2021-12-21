        REAL FUNCTION PDF_NDL(THETAN)
        save
c
        REAL THETAN

        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &          I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &          I_PDF_NDL,I_PDF_TRUNK

        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &                 I_PDF_NDL,I_PDF_TRUNK



c
        include 'constants.include'
c
        IF(I_PDF_NDL.EQ.1)THEN
            IF(thetan.GT.PI/2.)THEN
                PDF_NDL = 0.0
            ELSE IF(thetan.EQ.PI/2.)THEN
                PDF_NDL = 0.5
            ELSE
                PDF_NDL = SIN(thetan)
            ENDIF
c
        ELSE IF(I_PDF_NDL.EQ.2)THEN
            IF((thetan.GT.0.0).AND.(thetan.LT.(PI/2.0)))THEN
                  DUM = SIN(2*thetan)
                  DUM = DUM*DUM
                  PDF_NDL = (16.0/(3.0*PI))*DUM*DUM
            ELSE
                  PDF_NDL = 0.0
            ENDIF
c
        ELSE IF(I_PDF_NDL.EQ.3)THEN
            IF((thetan.GT.0).AND.(thetan.LT.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (4.0/PI)*DUM*DUM
            ELSE IF((thetan.EQ.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (2.0/PI)*DUM*DUM
            ELSE
                 PDF_NDL = 0.0
            ENDIF
c
        ELSE IF(I_PDF_NDL.EQ.4)THEN
            IF((thetan.GT.0).AND.(thetan.LT.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (1.5)*DUM*DUM*DUM
            ELSE IF((thetan.EQ.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (0.75)*DUM*DUM*DUM
            ELSE
                 PDF_NDL = 0.0
            ENDIF
c
        ELSE IF(I_PDF_NDL.EQ.5)THEN
            IF((thetan.GT.0).AND.(thetan.LT.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (1.697653)*(DUM**4)
            ELSE IF((thetan.EQ.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (0.848826)*(DUM**4)
            ELSE
                 PDF_NDL = 0.0
            ENDIF
c
        ELSE IF(I_PDF_NDL.EQ.6)THEN
            IF((thetan.GT.0).AND.(thetan.LT.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (15./8.)*(DUM**5)
            ELSE IF((thetan.EQ.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (15./16.)*(DUM**5)
            ELSE
                 PDF_NDL = 0.0
            ENDIF
c
        ELSE IF(I_PDF_NDL.EQ.7)THEN
            IF((thetan.GT.0).AND.(thetan.LT.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (6.4/pi)*(DUM**6)
            ELSE IF((thetan.EQ.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (3.2/pi)*(DUM**6)
            ELSE
                 PDF_NDL = 0.0
            ENDIF
c
        ELSE IF(I_PDF_NDL.EQ.8)THEN
            IF((thetan.GT.0).AND.(thetan.LT.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (2.187499)*(DUM**7)
            ELSE IF((thetan.EQ.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (1.093750)*(DUM**7)
            ELSE
                 PDF_NDL = 0.0
            ENDIF
C
        ELSE IF(I_PDF_NDL.EQ.9)THEN
            IF((thetan.GT.0).AND.(thetan.LT.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (2.328210)*(DUM**8)
            ELSE IF((thetan.EQ.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (1.164105)*(DUM**8)
            ELSE
                 PDF_NDL = 0.0
            ENDIF
c
        ELSE IF(I_PDF_NDL.EQ.10)THEN
            IF((thetan.GT.0).AND.(thetan.LT.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (2.460938)*(DUM**9)
            ELSE IF((thetan.EQ.PI/2.))THEN
                 DUM = SIN(thetan)
                 PDF_NDL = (1.230469)*(DUM**9)
            ELSE
                 PDF_NDL = 0.0
            ENDIF
c
        ELSE IF(I_PDF_NDL.EQ.11)THEN
            IF((thetan.GT.0).AND.(thetan.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(thetan-OFFSET))**9 
     &               + ABS(SIN(PI-thetan-OFFSET))**9 
                 PDF_NDL = (1.23046)*(DUM)
            ELSE
                 PDF_NDL = 0.0
            ENDIF
c
        ELSE IF(I_PDF_NDL.EQ.12)THEN
            IF((thetan.GT.0).AND.(thetan.LE.PI/2.))THEN
                 OFFSET = PI/3.
                 DUM = ABS(SIN(thetan+OFFSET))**9 
     &               + ABS(SIN(PI-thetan+OFFSET))**9 
                 PDF_NDL = (1.23)*(DUM)
            ELSE
                 PDF_NDL = 0.0
            ENDIF
c
        ELSE IF(I_PDF_NDL.EQ.13)THEN
            IF((thetan.GT.0).AND.(thetan.LE.PI/2.))THEN
                 PDF_NDL = (3.61311)*(sin(thetan))**20
            ELSE IF((thetan.EQ.PI/2.))THEN
                 PDF_NDL = (1.80656)*(sin(thetan))**20
            ELSE
                 PDF_NDL = 0.0
            ENDIF
c
        ELSE IF(I_PDF_NDL.EQ.14)THEN
            IF((thetan.GT.0).AND.(thetan.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(thetan+OFFSET))**20
     &               + ABS(SIN(PI-thetan+OFFSET))**20
                 PDF_NDL = (1.80656)*(DUM)
            ELSE
                 PDF_NDL = 0.0
            ENDIF
c
        ELSE IF(I_PDF_NDL.EQ.17)THEN
c
         IF(thetan.EQ.0.0)THEN
            dum = abs(cos(thetan))
            dum = dum**16.0
            PDF_NDL = 1.62088*dum
         ELSE IF((thetan.GT.0.0).AND.(thetan.LT.PI/2.))THEN
            dum = abs(cos(thetan))
            dum = dum**16.0
            PDF_NDL = 3.24176*dum
         ELSE
            PDF_NDL = 0.0
         ENDIF
c

        ELSE
            print*,' *** Bad value for needle pdf function ***'
            STOP
        ENDIF

C
        RETURN
        END
c
        REAL FUNCTION PDF_LF(THETAd)
        save
c
        REAL THETAd, MU, NU ,TOL
        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &          I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &          I_PDF_NDL,I_PDF_TRUNK

        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &                 I_PDF_NDL,I_PDF_TRUNK


c
        include 'constants.include'
c
        TOL = 8.727E-7
c
        IF(I_PDF_LEAF.EQ.1)THEN
            IF(THETAd.GT.PI/2.)THEN
                PDF_LF = 0.0
            ELSE IF(THETAd.EQ.PI/2.)THEN
                PDF_LF = 0.5
            ELSE
                PDF_LF = SIN(THETAd)
            ENDIF
c
        ELSE IF(I_PDF_LEAF.EQ.2)THEN
            IF(THETAd.GE.(PI/2.-TOL))THEN
                PDF_LF = 0.0
            else if(THETAd.LE.TOL)then
                PDF_LF = 0.0
            else
                nu = 2.770
                mu = 1.172
                call leaf_beta_dist(mu,nu,THETAd,PDF_LF)
            endif
c
        ELSE IF(I_PDF_LEAF.EQ.3)THEN
            IF(THETAd.GE.(PI/2.-TOL))THEN
                PDF_LF = 0.0
            else if(THETAd.LE.TOL)then
                PDF_LF = 0.0
            else
                nu = 1.1720
                mu = 2.7700
                call leaf_beta_dist(mu,nu,THETAd,PDF_LF)
            endif
c
        ELSE IF(I_PDF_LEAF.EQ.4)THEN
            IF(THETAd.GE.(PI/2.-TOL))THEN
                PDF_LF = 0.0
            else if(THETAd.LE.TOL)then
                PDF_LF = 0.0
            else
                nu = 3.326
                mu = 3.326
                call leaf_beta_dist(mu,nu,THETAd,PDF_LF)
            endif
c
        ELSE IF(I_PDF_LEAF.EQ.5)THEN
            IF(THETAd.GE.(PI/2.-TOL))THEN
                PDF_LF = 0.0
            else if(THETAd.LE.TOL)then
                PDF_LF = 0.0
            else
                nu = 0.433
                mu = 0.433
                call leaf_beta_dist(mu,nu,THETAd,PDF_LF)
            endif
c
        ELSE IF(I_PDF_LEAF.EQ.6)THEN
            IF(THETAd.GE.(PI/2.-TOL))THEN
                PDF_LF = 0.0
            else if(THETAd.LE.TOL)then
                PDF_LF = 0.0
            else
                nu = 1.000
                mu = 1.000
                call leaf_beta_dist(mu,nu,THETAd,PDF_LF)
            endif
c
        ELSE IF(I_PDF_LEAF.EQ.7)THEN
            IF(THETAd.GE.(PI/2.-TOL))THEN
                PDF_LF = 0.0
            else if(THETAd.LE.TOL)then
                PDF_LF = 0.0
            else
                nu = 1.101
                mu = 1.930
                call leaf_beta_dist(mu,nu,THETAd,PDF_LF)
            endif
c
        ELSE
            CALL WRITE_ERROR(111)
            STOP
        ENDIF
C
        RETURN
        END
c
        REAL FUNCTION PDF_BR(THETAc)
        save
c
        REAL THETAc, DUM, OFFSET
        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &          I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &          I_PDF_NDL,I_PDF_TRUNK

        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &                 I_PDF_NDL,I_PDF_TRUNK


c
        include 'constants.include'
c
        IF(I_PDF_BR_1.EQ.1)THEN
            IF(THETAc.GT.PI/2.)THEN
                PDF_BR = 0.0
            ELSE IF(THETAc.EQ.PI/2.)THEN
                PDF_BR = 0.5
            ELSE
                PDF_BR = SIN(THETAc)
            ENDIF
c
        ELSE IF(I_PDF_BR_1.EQ.2)THEN
            IF((THETAC.GT.0.0).AND.(THETAC.LT.(PI/2.0)))THEN
                  DUM = SIN(2*THETAC)
                  DUM = DUM*DUM
                  PDF_BR = (16.0/(3.0*PI))*DUM*DUM
            ELSE
                  PDF_BR = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_1.EQ.3)THEN
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (4.0/PI)*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (2.0/PI)*DUM*DUM
            ELSE
                 PDF_BR = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_1.EQ.4)THEN
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (1.5)*DUM*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (0.75)*DUM*DUM*DUM
            ELSE
                 PDF_BR = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_1.EQ.5)THEN
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (1.697653)*(DUM**4)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (0.848826)*(DUM**4)
            ELSE
                 PDF_BR = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_1.EQ.6)THEN
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (15./8.)*(DUM**5)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (15./16.)*(DUM**5)
            ELSE
                 PDF_BR = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_1.EQ.7)THEN
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (6.4/pi)*(DUM**6)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (3.2/pi)*(DUM**6)
            ELSE
                 PDF_BR = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_1.EQ.8)THEN
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (2.187499)*(DUM**7)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (1.093750)*(DUM**7)
            ELSE
                 PDF_BR = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_1.EQ.9)THEN
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (2.328210)*(DUM**8)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (1.164105)*(DUM**8)
            ELSE
                 PDF_BR = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_1.EQ.10)THEN
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (2.460938)*(DUM**9)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR = (1.230469)*(DUM**9)
            ELSE
                 PDF_BR = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_1.EQ.11)THEN
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc-OFFSET))**9 
     &               + ABS(SIN(PI-THETAc-OFFSET))**9 
                 PDF_BR = (1.23046)*(DUM)
            ELSE
                 PDF_BR = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_1.EQ.12)THEN
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/3.
                 DUM = ABS(SIN(THETAc+OFFSET))**9 
     &               + ABS(SIN(PI-THETAc+OFFSET))**9 
                 PDF_BR = (1.23)*(DUM)
            ELSE
                 PDF_BR = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_1.EQ.13)THEN
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 PDF_BR = (3.61311)*(sin(thetac))**20
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 PDF_BR = (1.80656)*(sin(thetac))**20
            ELSE
                 PDF_BR = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_1.EQ.14)THEN
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc+OFFSET))**20
     &               + ABS(SIN(PI-THETAc+OFFSET))**20
                 PDF_BR = (1.80656)*(DUM)
            ELSE
                 PDF_BR = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_1.EQ.17)THEN
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 64.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**16
     &               + ABS(SIN(PI-THETAc+OFFSET))**16
                 PDF_BR = (1.62088)*(DUM)
            ELSE
                 PDF_BR = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_1.EQ.18)THEN
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 41.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**8
     &               + ABS(SIN(PI-THETAc+OFFSET))**8
                 PDF_BR = (1.1641)*(DUM)
            ELSE
                 PDF_BR = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_1.EQ.23)THEN
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = COS(THETAc)
                 PDF_BR = (4.0/PI)*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = COS(THETAc)
                 PDF_BR = (2.0/PI)*DUM*DUM
            ELSE
                 PDF_BR = 0.0
            ENDIF
c

        ELSE
            CALL WRITE_ERROR(110)
            STOP
        ENDIF
C  
        RETURN
        END
c
        REAL FUNCTION PDF_BR2(THETAc)
        save
c
        REAL THETAc, DUM, OFFSET
        INTEGER I_PDF_LEAF, I_PDF_BR_1, I_PDF_BR_2,
     &          I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &          I_PDF_NDL,I_PDF_TRUNK

        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &                 I_PDF_NDL,I_PDF_TRUNK

c
        include 'constants.include'
c
        IF(I_PDF_BR_2.EQ.1)THEN
            IF(THETAc.GT.PI/2.)THEN
                PDF_BR2 = 0.0
            ELSE IF(THETAc.EQ.PI/2.)THEN
                PDF_BR2 = 0.5
            ELSE
                PDF_BR2 = SIN(THETAc)
            ENDIF
c
        ELSE IF(I_PDF_BR_2.EQ.2)THEN
            IF((THETAC.GT.0.0).AND.(THETAC.LT.(PI/2.0)))THEN
                  DUM = SIN(2*THETAC)
                  DUM = DUM*DUM
                  PDF_BR2 = (16.0/(3.0*PI))*DUM*DUM
            ELSE
                  PDF_BR2 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_2.EQ.3)THEN
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (4.0/PI)*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (2.0/PI)*DUM*DUM
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_2.EQ.4)THEN
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (1.5)*DUM*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (0.75)*DUM*DUM*DUM
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_2.EQ.5)THEN
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (1.697653)*(DUM**4)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (0.848826)*(DUM**4)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_2.EQ.6)THEN
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (15./8.)*(DUM**5)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (15./16.)*(DUM**5)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_2.EQ.7)THEN
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (6.4/pi)*(DUM**6)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (3.2/pi)*(DUM**6)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_2.EQ.8)THEN
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (2.187499)*(DUM**7)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (1.093750)*(DUM**7)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_2.EQ.9)THEN
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (2.328210)*(DUM**8)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (1.164105)*(DUM**8)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_2.EQ.10)THEN
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (2.460938)*(DUM**9)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR2 = (1.230469)*(DUM**9)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_2.EQ.11)THEN
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc-OFFSET))**9 
     &               + ABS(SIN(PI-THETAc-OFFSET))**9 
                 PDF_BR2 = (1.23046)*(DUM)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_2.EQ.12)THEN
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/3.
                 DUM = ABS(SIN(THETAc+OFFSET))**9 
     &               + ABS(SIN(PI-THETAc+OFFSET))**9 
                 PDF_BR2 = (1.23)*(DUM)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_2.EQ.13)THEN
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 PDF_BR2 = (3.61311)*(sin(thetac))**20
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 PDF_BR2 = (1.80656)*(sin(thetac))**20
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_2.EQ.14)THEN
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc+OFFSET))**20
     &               + ABS(SIN(PI-THETAc+OFFSET))**20
                 PDF_BR2 = (1.80656)*(DUM)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
        ELSE IF(I_PDF_BR_2.EQ.17)THEN
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 64.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**16
     &               + ABS(SIN(PI-THETAc+OFFSET))**16
                 PDF_BR2 = (1.62088)*(DUM)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF

        ELSE IF(I_PDF_BR_2.EQ.18)THEN
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 41.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**8
     &               + ABS(SIN(PI-THETAc+OFFSET))**8
                 PDF_BR2 = (1.1641)*(DUM)
            ELSE
                 PDF_BR2 = 0.0
            ENDIF
c

        ELSE
            CALL WRITE_ERROR(110)
            STOP
        ENDIF
C  
        RETURN
        END
c
        REAL FUNCTION PDF_TR(THETAC)
        save
c
        REAL THETAC, dum
        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &          I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &          I_PDF_NDL,I_PDF_TRUNK

        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &                 I_PDF_NDL,I_PDF_TRUNK

c
        include 'constants.include'
c
        IF(I_PDF_TRUNK.EQ.1)THEN
            IF(THETAc.eq.0.0)THEN
                PDF_TR = 11.4592
            ELSE 
                PDF_TR = 0.0
            ENDIF
C
        ELSE IF(I_PDF_TRUNK.EQ.2)THEN
c
         IF(THETAC.EQ.0.0)THEN
            dum = abs(cos(thetac))
            dum = dum**8.0
            PDF_TR = 1.1641*dum
         ELSE IF((THETAC.GT.0.0).AND.(THETAC.LT.PI/2.))THEN
            dum = abs(cos(thetac))
            dum = dum**8.0
            PDF_TR = 2.32821*dum
         ELSE
            PDF_TR = 0.0
         ENDIF
C
        ELSE IF(I_PDF_TRUNK.EQ.3)THEN
c
         IF(THETAC.EQ.0.0)THEN
            dum = abs(cos(thetac))
            dum = dum**6.0
            PDF_TR = 1.01859*dum
         ELSE IF((THETAC.GE.0.0).AND.(THETAC.LT.PI/2.))THEN
            dum = abs(cos(thetac))
            dum = dum**6.0
            PDF_TR = 2.03718*dum
         ELSE
            PDF_TR = 0.0
         ENDIF
C
        ELSE IF(I_PDF_TRUNK.EQ.4)THEN
c
         IF((THETAC.GE.0.0).AND.(THETAC.LT.PI/2.))THEN
            PDF_TR = 8.0*(EXP(-8.0*THETAC))
         ELSE
            PDF_TR = 0.0
         ENDIF
c
        ELSE
            PRINT*,'IMPROPER TRUNK LAYER PDF'
            STOP
        ENDIF
C
        RETURN
        END
c
        subroutine leaf_beta_dist(mu,nu,theta,pdf)
c
        real mu,nu,pi,pdf,theta,gammln
c
        pi = 3.141592654
c
        pdf = (2./pi)*exp(gammln(mu+nu) - gammln(mu) - gammln(nu))
     &        *((1.0-theta/(pi/2.))**(nu-1.))*
     &          ((theta/(pi/2.))**(mu-1.))
c
        return
        end
c
      FUNCTION GAMMLN(XX)
      SAVE
c
      INTEGER J
      REAL*4 GAMMLN,XX
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END

c***********************************************************************
