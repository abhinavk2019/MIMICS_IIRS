        SUBROUTINE WRITE_ERROR(IFLAG)
        save
c
        INTEGER IFLAG
c
         OPEN(UNIT=1,FILE='../output_dir/error_messages.output')
        if(iflag.eq.2) then
            WRITE(1,2)
        else if(iflag.eq.3) then
            WRITE(1,3)
        else if(iflag.eq.10) then
            WRITE(1,10)
        else if(iflag.eq.15) then
            WRITE(1,15)
        else if(iflag.eq.20) then
            WRITE(1,20)
        else if(iflag.eq.25) then
            WRITE(1,25)
        else if(iflag.eq.30) then
            WRITE(1,30)
        else if(iflag.eq.35) then
            WRITE(1,35)
        else if(iflag.eq.40) then
            WRITE(1,40)
        else if(iflag.eq.50) then
            WRITE(1,50)
        else if(iflag.eq.100) then
            WRITE(1,100)
        else if(iflag.eq.110) then
            WRITE(1,110)
        else if(iflag.eq.111) then
            WRITE(1,111)
        else if(iflag.eq.115) then
            WRITE(1,115)
        else if(iflag.eq.120) then
            WRITE(1,120)
        else if(iflag.eq.125) then
            WRITE(1,125)
        endif
        CLOSE(1)
c
2       FORMAT('*** ERROR -- Secondary branches are not yet ',
     &         'incorporated in the model *** (routine read_configure)')
3       FORMAT('*** ERROR -- VARIATION IN PHI OR SIZE MUST BE DEFAULT ',
     &         'P.D.F. *** (routine read_configure)')
10      FORMAT('***ERROR,I not equal to 0,1 or 2 in SET_TWO_FLAGS**')
15      FORMAT('*** ERROR -- I not equal to 0 or 1 in SET_FLAG ***')
20      FORMAT(' *** ERROR *** DELTA value is equal to zero ')
25      FORMAT('*** ERROR  -- SURFACE TYPE IS OUT OF BOUNDS ***')
30      FORMAT('*** ERROR  VALUE OF J OUT OF BOUNDS IN ROUTINE STEP')
35      FORMAT('*** ERROR  VALUE OF J OUT OF BOUNDS IN ROUTINE RESET')
40      FORMAT(' *** ERROR ***  soil gravometric moiture not',
     &            ' yet incorporated for frozen soil.')
50      FORMAT('*** ERROR ** ONLY SOIL SURFACE ALLOWED FOR BACKSCATTER')
100     FORMAT('*** ERROR ** Value of CHAR not set in WRITE_DATA ***')
110     FORMAT('*** ERROR ** PDF flag for primary branches ',
     &          'is out of range.')
111     FORMAT('*** ERROR ** PDF flag for leaves ',
     &          'is out of range.')
115     FORMAT('*** ERROR ** PDF flag for trunk size ',
     &          'is out of range.')
120     FORMAT('*** ERROR ** DIELECTRIC INPUT FLAG IS OUT OF RANGE',
     &          'in READ_EPS_TABLE')
125     FORMAT('*** ERROR ** HISTOGRAM INPUT FLAG IS OUT OF RANGE',
     &          'in READ_HISTOGRAM')
C
        RETURN
        END
