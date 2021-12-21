        REAL FUNCTION TRUNK_HGHT_FUNC(DBH)
        save

        REAL DBH,DBHCU,MEAND,CUTOFF
      
        MEAND=10.0 
        CUTOFF=8.0
	DBHCU=0.00221*DBH*DBH*DBH

	IF (DBH.LT.CUTOFF) THEN
       		TRUNK_HGHT_FUNC = 0.1437+0.683*DBH+0.01517*DBH*DBH-DBHCU
        ELSE IF (MEAND.LT.CUTOFF) THEN
	        TRUNK_HGHT_FUNC = 7.3132 
	ELSE
                TRUNK_HGHT_FUNC = 10**(0.8688*ALOG10(MEAND+0.0795))        
        ENDIF 

        RETURN
        END
C***********************************************************************
