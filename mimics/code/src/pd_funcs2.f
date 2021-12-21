        REAL FUNCTION PDF_BR3(THETAc)
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
        IF(I_PDF_BR_3.EQ.1)THEN
c
            IF(THETAc.GT.PI/2.)THEN
                PDF_BR3 = 0.0
            ELSE IF(THETAc.EQ.PI/2.)THEN
                PDF_BR3 = 0.5
            ELSE
                PDF_BR3 = SIN(THETAc)
            ENDIF
c
        ELSE IF(I_PDF_BR_3.EQ.2)THEN
c
            IF((THETAC.GT.0.0).AND.(THETAC.LT.(PI/2.0)))THEN
                  DUM = SIN(2*THETAC)
                  DUM = DUM*DUM
                  PDF_BR3 = (16.0/(3.0*PI))*DUM*DUM
            ELSE
                  PDF_BR3 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_3.EQ.3)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (4.0/PI)*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (2.0/PI)*DUM*DUM
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_3.EQ.4)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (1.5)*DUM*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (0.75)*DUM*DUM*DUM
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_3.EQ.5)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (1.697653)*(DUM**4)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (0.848826)*(DUM**4)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_3.EQ.6)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (15./8.)*(DUM**5)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (15./16.)*(DUM**5)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_3.EQ.7)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (6.4/pi)*(DUM**6)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (3.2/pi)*(DUM**6)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_3.EQ.8)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (2.187499)*(DUM**7)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (1.093750)*(DUM**7)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_3.EQ.9)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (2.328210)*(DUM**8)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (1.164105)*(DUM**8)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_3.EQ.10)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (2.460938)*(DUM**9)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR3 = (1.230469)*(DUM**9)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_3.EQ.11)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc-OFFSET))**9 
     &               + ABS(SIN(PI-THETAc-OFFSET))**9 
                 PDF_BR3 = (1.23046)*(DUM)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_3.EQ.12)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/3.
                 DUM = ABS(SIN(THETAc+OFFSET))**9 
     &               + ABS(SIN(PI-THETAc+OFFSET))**9 
                 PDF_BR3 = (1.23)*(DUM)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_3.EQ.13)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 PDF_BR3 = (3.61311)*(sin(thetac))**20
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 PDF_BR3 = (1.80656)*(sin(thetac))**20
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_3.EQ.14)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc+OFFSET))**20
     &               + ABS(SIN(PI-THETAc+OFFSET))**20
                 PDF_BR3 = (1.80656)*(DUM)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_3.EQ.17)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 64.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**16
     &               + ABS(SIN(PI-THETAc+OFFSET))**16
                 PDF_BR3 = (1.62088)*(DUM)
            ELSE
                 PDF_BR3 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_3.EQ.18)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 41.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**8
     &               + ABS(SIN(PI-THETAc+OFFSET))**8
                 PDF_BR3 = (1.1641)*(DUM)
            ELSE
                 PDF_BR3 = 0.0
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
        REAL FUNCTION PDF_BR4(THETAc)
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
        IF(I_PDF_BR_4.EQ.1)THEN
c
            IF(THETAc.GT.PI/2.)THEN
                PDF_BR4 = 0.0
            ELSE IF(THETAc.EQ.PI/2.)THEN
                PDF_BR4 = 0.5
            ELSE
                PDF_BR4 = SIN(THETAc)
            ENDIF
c
        ELSE IF(I_PDF_BR_4.EQ.2)THEN
c
            IF((THETAC.GT.0.0).AND.(THETAC.LT.(PI/2.0)))THEN
                  DUM = SIN(2*THETAC)
                  DUM = DUM*DUM
                  PDF_BR4 = (16.0/(3.0*PI))*DUM*DUM
            ELSE
                  PDF_BR4 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_4.EQ.3)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (4.0/PI)*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (2.0/PI)*DUM*DUM
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_4.EQ.4)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (1.5)*DUM*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (0.75)*DUM*DUM*DUM
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_4.EQ.5)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (1.697653)*(DUM**4)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (0.848826)*(DUM**4)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_4.EQ.6)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (15./8.)*(DUM**5)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (15./16.)*(DUM**5)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_4.EQ.7)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (6.4/pi)*(DUM**6)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (3.2/pi)*(DUM**6)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_4.EQ.8)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (2.187499)*(DUM**7)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (1.093750)*(DUM**7)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_4.EQ.9)THEN
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (2.328210)*(DUM**8)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (1.164105)*(DUM**8)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_4.EQ.10)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (2.460938)*(DUM**9)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR4 = (1.230469)*(DUM**9)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_4.EQ.11)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc-OFFSET))**9 
     &               + ABS(SIN(PI-THETAc-OFFSET))**9 
                 PDF_BR4 = (1.23046)*(DUM)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_4.EQ.12)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/3.
                 DUM = ABS(SIN(THETAc+OFFSET))**9 
     &               + ABS(SIN(PI-THETAc+OFFSET))**9 
                 PDF_BR4 = (1.23)*(DUM)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_4.EQ.13)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 PDF_BR4 = (3.61311)*(sin(thetac))**20
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 PDF_BR4 = (1.80656)*(sin(thetac))**20
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_4.EQ.14)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc+OFFSET))**20
     &               + ABS(SIN(PI-THETAc+OFFSET))**20
                 PDF_BR4 = (1.80656)*(DUM)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_4.EQ.17)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 64.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**16
     &               + ABS(SIN(PI-THETAc+OFFSET))**16
                 PDF_BR4 = (1.62088)*(DUM)
            ELSE
                 PDF_BR4 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_4.EQ.18)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 41.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**8
     &               + ABS(SIN(PI-THETAc+OFFSET))**8
                 PDF_BR4 = (1.1641)*(DUM)
            ELSE
                 PDF_BR4 = 0.0
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
        REAL FUNCTION PDF_BR5(THETAc)
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
        IF(I_PDF_BR_5.EQ.1)THEN
c
            IF(THETAc.GT.PI/2.)THEN
                PDF_BR5 = 0.0
            ELSE IF(THETAc.EQ.PI/2.)THEN
                PDF_BR5 = 0.5
            ELSE
                PDF_BR5 = SIN(THETAc)
            ENDIF
c
        ELSE IF(I_PDF_BR_5.EQ.2)THEN
c
            IF((THETAC.GT.0.0).AND.(THETAC.LT.(PI/2.0)))THEN
                  DUM = SIN(2*THETAC)
                  DUM = DUM*DUM
                  PDF_BR5 = (16.0/(3.0*PI))*DUM*DUM
            ELSE
                  PDF_BR5 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_5.EQ.3)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (4.0/PI)*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (2.0/PI)*DUM*DUM
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_5.EQ.4)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (1.5)*DUM*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (0.75)*DUM*DUM*DUM
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_5.EQ.5)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (1.697653)*(DUM**4)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (0.848826)*(DUM**4)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_5.EQ.6)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (15./8.)*(DUM**5)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (15./16.)*(DUM**5)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_5.EQ.7)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (6.4/pi)*(DUM**6)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (3.2/pi)*(DUM**6)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_5.EQ.8)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (2.187499)*(DUM**7)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (1.093750)*(DUM**7)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_5.EQ.9)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (2.328210)*(DUM**8)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (1.164105)*(DUM**8)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_5.EQ.10)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (2.460938)*(DUM**9)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR5 = (1.230469)*(DUM**9)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_5.EQ.11)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc-OFFSET))**9 
     &               + ABS(SIN(PI-THETAc-OFFSET))**9 
                 PDF_BR5 = (1.23046)*(DUM)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_5.EQ.12)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/3.
                 DUM = ABS(SIN(THETAc+OFFSET))**9 
     &               + ABS(SIN(PI-THETAc+OFFSET))**9 
                 PDF_BR5 = (1.23)*(DUM)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_5.EQ.13)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 PDF_BR5 = (3.61311)*(sin(thetac))**20
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 PDF_BR5 = (1.80656)*(sin(thetac))**20
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_5.EQ.14)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc+OFFSET))**20
     &               + ABS(SIN(PI-THETAc+OFFSET))**20
                 PDF_BR5 = (1.80656)*(DUM)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_5.EQ.17)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 64.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**16
     &               + ABS(SIN(PI-THETAc+OFFSET))**16
                 PDF_BR5 = (1.62088)*(DUM)
            ELSE
                 PDF_BR5 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_5.EQ.18)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 41.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**8
     &               + ABS(SIN(PI-THETAc+OFFSET))**8
                 PDF_BR5 = (1.1641)*(DUM)
            ELSE
                 PDF_BR5 = 0.0
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
        REAL FUNCTION PDF_BR6(THETAc)
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
        IF(I_PDF_BR_6.EQ.1)THEN
c
            IF(THETAc.GT.PI/2.)THEN
                PDF_BR6 = 0.0
            ELSE IF(THETAc.EQ.PI/2.)THEN
                PDF_BR6 = 0.5
            ELSE
                PDF_BR6 = SIN(THETAc)
            ENDIF
c
        ELSE IF(I_PDF_BR_6.EQ.2)THEN
c
            IF((THETAC.GT.0.0).AND.(THETAC.LT.(PI/2.0)))THEN
                  DUM = SIN(2*THETAC)
                  DUM = DUM*DUM
                  PDF_BR6 = (16.0/(3.0*PI))*DUM*DUM
            ELSE
                  PDF_BR6 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_6.EQ.3)THEN
C
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (4.0/PI)*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (2.0/PI)*DUM*DUM
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_6.EQ.4)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (1.5)*DUM*DUM*DUM
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (0.75)*DUM*DUM*DUM
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_6.EQ.5)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (1.697653)*(DUM**4)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (0.848826)*(DUM**4)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_6.EQ.6)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (15./8.)*(DUM**5)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (15./16.)*(DUM**5)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_6.EQ.7)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (6.4/pi)*(DUM**6)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (3.2/pi)*(DUM**6)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_6.EQ.8)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (2.187499)*(DUM**7)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (1.093750)*(DUM**7)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_6.EQ.9)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (2.328210)*(DUM**8)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (1.164105)*(DUM**8)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_6.EQ.10)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LT.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (2.460938)*(DUM**9)
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 DUM = SIN(THETAc)
                 PDF_BR6 = (1.230469)*(DUM**9)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_6.EQ.11)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc-OFFSET))**9 
     &               + ABS(SIN(PI-THETAc-OFFSET))**9 
                 PDF_BR6 = (1.23046)*(DUM)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_6.EQ.12)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/3.
                 DUM = ABS(SIN(THETAc+OFFSET))**9 
     &               + ABS(SIN(PI-THETAc+OFFSET))**9 
                 PDF_BR6 = (1.23)*(DUM)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_6.EQ.13)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 PDF_BR6 = (3.61311)*(sin(thetac))**20
            ELSE IF((THETAC.EQ.PI/2.))THEN
                 PDF_BR6 = (1.80656)*(sin(thetac))**20
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_6.EQ.14)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = PI/6.
                 DUM = ABS(SIN(THETAc+OFFSET))**20
     &               + ABS(SIN(PI-THETAc+OFFSET))**20
                 PDF_BR6 = (1.80656)*(DUM)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_6.EQ.17)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 64.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**16
     &               + ABS(SIN(PI-THETAc+OFFSET))**16
                 PDF_BR6 = (1.62088)*(DUM)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
c
        ELSE IF(I_PDF_BR_6.EQ.18)THEN
c
            IF((THETAC.GT.0).AND.(THETAC.LE.PI/2.))THEN
                 OFFSET = 41.*(pi/180.)
                 DUM = ABS(SIN(THETAc+OFFSET))**8
     &               + ABS(SIN(PI-THETAc+OFFSET))**8
                 PDF_BR6 = (1.1641)*(DUM)
            ELSE
                 PDF_BR6 = 0.0
            ENDIF
c
        ELSE
            CALL WRITE_ERROR(110)
            STOP
        ENDIF
C  
        RETURN
        END
