      REAL V(4)
      REAL dbh, h
      integer i
      CHARACTER*1 c1


 5    CONTINUE
      print *,' ' 
      PRINT *,' This program allows you to specify the coefficients'
      PRINT *,' of the polynomial that relates trunk diameter in cm'
      PRINT *,' to trunk height in meters. '
      print *,' '
      print *,' Enter the constant term:'
      READ  *,V(1)
      print *,' Enter the coeff for DIAM:'
      READ  *,V(2)
      print *,' Enter the coeff for DIAM^2:'
      READ  *,V(3)
      print *,' Enter the coeff for DIAM^3:'
      READ  *,V(4)
      print *,' Enter the coeff for DIAM^4:'
      READ  *,V(5)
      print *,' '
      print *,' The function you have defined will look like:'
      print *,' '
      print *,V(1),' + ',V(2),'*(DBH)',' + ',V(3),'*(DBH^2)',
     &             ' + ',V(4),'*(DBH^3)',' + ',V(5),'*(DBH^4)'
      print *,' '
      print *,' resulting in the following:'
      print *,' '
      print *,' diameter    calculated height'
      DO i=1,5
       dbh = 20. + (i-1)*50.
       h   = V(1) + V(2)*dbh + V(3)*dbh**2 + V(4)*dbh**3 + V(5)*dbh**4
       WRITE(*,10)dbh,h
 10    FORMAT(1X,F4.0,' cm',6X,F6.2,' meters')
      ENDDO
      print *,' '
      print *,' Is this satisfactory? '
      print *,'(Y will write the file, anything else will start over)'
      READ(*,'(a1)')c1
      IF(c1.NE.'Y')GOTO 5
       

               open(unit=25,file='trunk_hght_f.f',STATUS='UNKNOWN')
               write(25,301)
 301           format('        REAL FUNCTION TRUNK_HGHT_FUNC(DBH)')
               write(25,302)
 302           format('C')
               write(25,303)
 303           format('C   Function to compute trunk heights ',
     &                     'given the diameters')
               write(25,304)
 304           format('C')
               write(25,305)
 305           format('        REAL DBH')
               write(25,306)
 306           format('C******************************************')
               write(25,307)
 307           format(' ')
               write(25,308)v(1)
 308           format('        TRUNK_HGHT_FUNC=',e20.10,'+')
               write(25,309)v(2)
 309           format('     &                  ',e20.10,'*DBH+')
               write(25,310)v(3)
 310           format('     &                  ',e20.10,'*DBH*DBH +')
               write(25,311)v(4)
 311           format('     &                  ',e20.10,'*(DBH**3) +')
               write(25,312)v(5)
 312           format('     &                  ',e20.10,'*(DBH**4)')
               write(25,313)
 313           format(' ')
               write(25,314)
 314           format('           ')
               write(25,315)
 315           format('        RETURN')
               write(25,316)
 316           format('        END')
               write(25,317)
 317           format('C*******************************************')

               close(unit=25)

         print *,'wrote trunk_hght_f.f'
         stop
         end
