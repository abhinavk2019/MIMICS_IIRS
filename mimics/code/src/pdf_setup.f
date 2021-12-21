        SUBROUTINE PDF_NDL_SETUP(n_theta_1,t_rad_start_1,t_rad_stop_1,
     & delta_t_rad_1,n_theta_2,t_rad_start_2,t_rad_stop_2,delta_t_rad_2,
     &               n_theta_3,t_rad_start_3,t_rad_stop_3,delta_t_rad_3,
     &               n_phi_1,p_rad_start_1,p_rad_stop_1,delta_p_rad_1,
     &               n_phi_2,p_rad_start_2,p_rad_stop_2,delta_p_rad_2,
     &               n_phi_3,p_rad_start_3,p_rad_stop_3,delta_p_rad_3)
c
        save
C
        INTEGER n_theta_1,n_theta_2,n_theta_3, n_phi_1,n_phi_2,n_phi_3
        REAL t_rad_start_1, t_rad_stop_1, delta_t_rad_1
        REAL t_rad_start_2, t_rad_stop_2, delta_t_rad_2
        REAL t_rad_start_3, t_rad_stop_3, delta_t_rad_3
        REAL p_rad_start_1, p_rad_stop_1, delta_p_rad_1
        REAL p_rad_start_2, p_rad_stop_2, delta_p_rad_2
        REAL p_rad_start_3, p_rad_stop_3, delta_p_rad_3
C
        REAL t_deg_start_1, t_deg_stop_1, delta_t_deg_1
        REAL t_deg_start_2, t_deg_stop_2, delta_t_deg_2
        REAL t_deg_start_3, t_deg_stop_3, delta_t_deg_3
        REAL p_deg_start_1, p_deg_stop_1, delta_p_deg_1
        REAL p_deg_start_2, p_deg_stop_2, delta_p_deg_2
        REAL p_deg_start_3, p_deg_stop_3, delta_p_deg_3
C
        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &          I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &          I_PDF_NDL,I_PDF_TRUNK

        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &                 I_PDF_NDL,I_PDF_TRUNK


C
        REAL WORK
c
        INCLUDE 'constants.include'
c
        WORK = PI/180.
c
        if((i_pdf_ndl.eq.1).or.(i_pdf_ndl.eq.2))then
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 27.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 6
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 30.0
            t_deg_stop_2  = 150.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 121
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 152.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 6
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if(i_pdf_ndl.eq.3)then
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 32.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 7
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 35.0
            t_deg_stop_2  = 145.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 111
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 147.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 7
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if((i_pdf_ndl.eq.4).or.(i_pdf_ndl.eq.5))then
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 42.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 9
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 45.0
            t_deg_stop_2  = 135.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 91
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 147.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 7
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if(i_pdf_ndl.eq.6)then
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 47.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 10
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 50.0
            t_deg_stop_2  = 130.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 81
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 132.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 10
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if((i_pdf_ndl.eq.7).or.(i_pdf_ndl.eq.8).or.
     &                                (i_pdf_ndl.eq.9))then
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 52.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 11
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 55.0
            t_deg_stop_2  = 125.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 71
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 127.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 11
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if(i_pdf_ndl.eq.10)then
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 57.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 12
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 60.0
            t_deg_stop_2  = 120.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 61
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 122.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 12
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if(i_pdf_ndl.eq.11)then
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 27.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 6
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 30.0
            t_deg_stop_2  = 90.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 61
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 92.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 18
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if(i_pdf_ndl.eq.12)then
            t_deg_start_1 = 1.0
            t_deg_stop_1  = 60.0
            delta_t_deg_1 = 1.0
            n_theta_1 = 60
c
            p_deg_start_1 = 0.5
            p_deg_stop_1  = 179.5
            delta_p_deg_1 = 1.0
            n_phi_1 = 180
c
            t_deg_start_2 = 62.5
            t_deg_stop_2  = 177.5
            delta_t_deg_2 = 5.0
            n_theta_2 = 24
c
            p_deg_start_2 = 2.5
            p_deg_stop_2  = 177.5
            delta_p_deg_2 = 5.0
            n_phi_2 = 36
c
            t_deg_start_3 = 0.0
            t_deg_stop_3  = 0.0
            delta_t_deg_3 = 0.0
            n_theta_3 = 0
c
            p_deg_start_3 = 0.0
            p_deg_stop_3  = 0.0
            delta_p_deg_3 = 0.0
            n_phi_3 = 0
c
        else if(i_pdf_ndl.eq.13)then
c
            t_deg_start_1 = 7.5
            t_deg_stop_1  = 57.5
            delta_t_deg_1 = 10.0
            n_theta_1 = 0
c
            p_deg_start_1 = 0.0
            p_deg_stop_1  = 360.0
            delta_p_deg_1 = 10.0
            n_phi_1 = 0
c
            t_deg_start_2 = 60.0
            t_deg_stop_2  = 120.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 61
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 122.5
            t_deg_stop_3  = 172.5
            delta_t_deg_3 = 10.0
            n_theta_3 = 0
c
            p_deg_start_3 = 0.0
            p_deg_stop_3  = 360.0
            delta_p_deg_3 = 10.0
            n_phi_3 = 0
c
        else if(i_pdf_ndl.eq.14)then
c
            t_deg_start_1 = 0.0
            t_deg_stop_1  = 30.0
            delta_t_deg_1 = 10.0
            n_theta_1 = 0
c
            p_deg_start_1 = 0.0
            p_deg_stop_1  = 360.0
            delta_p_deg_1 = 10.0
            n_phi_1 = 0
c
            t_deg_start_2 = 30.0
            t_deg_stop_2  = 90.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 61
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 90.
            t_deg_stop_3  = 180.0
            delta_t_deg_3 = 10.0
            n_theta_3 = 0
c
            p_deg_start_3 = 0.0
            p_deg_stop_3  = 360.0
            delta_p_deg_3 = 10.0
            n_phi_3 = 0
c
        else if(i_pdf_ndl.eq.17)then
c
            t_deg_start_1 =  1.0
            t_deg_stop_1  = 30.0
            delta_t_deg_1 =  1.0
            n_theta_1 = 30
c
            p_deg_start_1 =   0.5
            p_deg_stop_1  = 179.5
            delta_p_deg_1 =   1.0
            n_phi_1 = 180
c
            t_deg_start_2 = 35.0
            t_deg_stop_2  = 85.0
            delta_t_deg_2 =  5.0
            n_theta_2 = 11
c
            p_deg_start_2 =   5.0
            p_deg_stop_2  = 175.0
            delta_p_deg_2 =  10.0
            n_phi_2 = 18
c
            t_deg_start_3 = 90.
            t_deg_stop_3  = 180.0
            delta_t_deg_3 = 10.0
            n_theta_3 = 0
c
            p_deg_start_3 = 0.0
            p_deg_stop_3  = 360.0
            delta_p_deg_3 = 10.0
            n_phi_3 = 0
c
        else
c
            t_deg_start_1 = 2.5 
            t_deg_stop_1  = 177.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 35
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 0.0
            t_deg_stop_2  = 0.0
            delta_t_deg_2 = 0.0
            n_theta_2 = 0
c
            p_deg_start_2 = 0.0
            p_deg_stop_2  = 0.0
            delta_p_deg_2 = 0.0
            n_phi_2 = 0
c
            t_deg_start_3 = 0.0
            t_deg_stop_3  = 0.0
            delta_t_deg_3 = 0.0
            n_theta_3 = 0
c
            p_deg_start_3 = 0.0
            p_deg_stop_3  = 0.0
            delta_p_deg_3 = 0.0
            n_phi_3 = 0
c
        endif
c
c-----------------------------------------------------------------------
c
        t_rad_start_1 =  WORK*t_deg_start_1
        t_rad_stop_1  =  WORK*t_deg_stop_1 
        delta_t_rad_1 =  WORK*delta_t_deg_1
                                     
        p_rad_start_1 =  WORK*p_deg_start_1
        p_rad_stop_1  =  WORK*p_deg_stop_1 
        delta_p_rad_1 =  WORK*delta_p_deg_1
                                     
        t_rad_start_2 =  WORK*t_deg_start_2
        t_rad_stop_2  =  WORK*t_deg_stop_2 
        delta_t_rad_2 =  WORK*delta_t_deg_2
                                     
        p_rad_start_2 =  WORK*p_deg_start_2
        p_rad_stop_2  =  WORK*p_deg_stop_2 
        delta_p_rad_2 =  WORK*delta_p_deg_2
                                     
        t_rad_start_3 =  WORK*t_deg_start_3
        t_rad_stop_3  =  WORK*t_deg_stop_3 
        delta_t_rad_3 =  WORK*delta_t_deg_3
                                     
        p_rad_start_3 =  WORK*p_deg_start_3
        p_rad_stop_3  =  WORK*p_deg_stop_3 
        delta_p_rad_3 =  WORK*delta_p_deg_3
c
        RETURN
        END
c
        SUBROUTINE PDF_LF_SETUP(n_theta_1,t_rad_start_1,t_rad_stop_1,
     & delta_t_rad_1,n_theta_2,t_rad_start_2,t_rad_stop_2,delta_t_rad_2,
     &               n_theta_3,t_rad_start_3,t_rad_stop_3,delta_t_rad_3,
     &               n_phi_1,p_rad_start_1,p_rad_stop_1,delta_p_rad_1,
     &               n_phi_2,p_rad_start_2,p_rad_stop_2,delta_p_rad_2,
     &               n_phi_3,p_rad_start_3,p_rad_stop_3,delta_p_rad_3)
c
        save
C
        INTEGER n_theta_1,n_theta_2,n_theta_3, n_phi_1,n_phi_2,n_phi_3
        REAL t_rad_start_1, t_rad_stop_1, delta_t_rad_1
        REAL t_rad_start_2, t_rad_stop_2, delta_t_rad_2
        REAL t_rad_start_3, t_rad_stop_3, delta_t_rad_3
        REAL p_rad_start_1, p_rad_stop_1, delta_p_rad_1
        REAL p_rad_start_2, p_rad_stop_2, delta_p_rad_2
        REAL p_rad_start_3, p_rad_stop_3, delta_p_rad_3
C
        REAL t_deg_start_1, t_deg_stop_1, delta_t_deg_1
        REAL t_deg_start_2, t_deg_stop_2, delta_t_deg_2
        REAL t_deg_start_3, t_deg_stop_3, delta_t_deg_3
        REAL p_deg_start_1, p_deg_stop_1, delta_p_deg_1
        REAL p_deg_start_2, p_deg_stop_2, delta_p_deg_2
        REAL p_deg_start_3, p_deg_stop_3, delta_p_deg_3
C
        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &          I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &          I_PDF_NDL,I_PDF_TRUNK

        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &                 I_PDF_NDL,I_PDF_TRUNK


C
        REAL WORK
c
        INCLUDE 'constants.include'
c
        WORK = PI/180.
C
        t_deg_start_1 = 2.5
        t_deg_stop_1  = 177.5
        delta_t_deg_1 = 5.0
        n_theta_1 = 35
c
        p_deg_start_1 = 2.5
        p_deg_stop_1  = 177.5
        delta_p_deg_1 = 5.0
        n_phi_1 = 36
c
c
        t_deg_start_2 = 0.0
        t_deg_stop_2  = 0.0
        delta_t_deg_2 = 0.0
        n_theta_2 = 0
c
        p_deg_start_2 = 0.0
        p_deg_stop_2  = 0.0
        delta_p_deg_2 = 0.0
        n_phi_2 = 0
c
c
        t_deg_start_3 = 0.0
        t_deg_stop_3  = 0.0
        delta_t_deg_3 = 0.0
        n_theta_3 = 0
c
        p_deg_start_3 = 0.0
        p_deg_stop_3  = 0.0
        delta_p_deg_3 = 0.0
        n_phi_3 = 0
C
C
        t_rad_start_1 =  WORK*t_deg_start_1
        t_rad_stop_1  =  WORK*t_deg_stop_1 
        delta_t_rad_1 =  WORK*delta_t_deg_1
                                     
        p_rad_start_1 =  WORK*p_deg_start_1
        p_rad_stop_1  =  WORK*p_deg_stop_1 
        delta_p_rad_1 =  WORK*delta_p_deg_1
                                     
        t_rad_start_2 =  WORK*t_deg_start_2
        t_rad_stop_2  =  WORK*t_deg_stop_2 
        delta_t_rad_2 =  WORK*delta_t_deg_2
                                     
        p_rad_start_2 =  WORK*p_deg_start_2
        p_rad_stop_2  =  WORK*p_deg_stop_2 
        delta_p_rad_2 =  WORK*delta_p_deg_2
                                     
        t_rad_start_3 =  WORK*t_deg_start_3
        t_rad_stop_3  =  WORK*t_deg_stop_3 
        delta_t_rad_3 =  WORK*delta_t_deg_3
                                     
        p_rad_start_3 =  WORK*p_deg_start_3
        p_rad_stop_3  =  WORK*p_deg_stop_3 
        delta_p_rad_3 =  WORK*delta_p_deg_3
C
        RETURN
        END
c
        SUBROUTINE PDF_BR_SETUP(n_theta_1,t_rad_start_1,t_rad_stop_1,
     & delta_t_rad_1,n_theta_2,t_rad_start_2,t_rad_stop_2,delta_t_rad_2,
     &               n_theta_3,t_rad_start_3,t_rad_stop_3,delta_t_rad_3,
     &               n_phi_1,p_rad_start_1,p_rad_stop_1,delta_p_rad_1,
     &               n_phi_2,p_rad_start_2,p_rad_stop_2,delta_p_rad_2,
     &               n_phi_3,p_rad_start_3,p_rad_stop_3,delta_p_rad_3)
c
        save
C
        INTEGER n_theta_1,n_theta_2,n_theta_3, n_phi_1,n_phi_2,n_phi_3
        REAL t_rad_start_1, t_rad_stop_1, delta_t_rad_1
        REAL t_rad_start_2, t_rad_stop_2, delta_t_rad_2
        REAL t_rad_start_3, t_rad_stop_3, delta_t_rad_3
        REAL p_rad_start_1, p_rad_stop_1, delta_p_rad_1
        REAL p_rad_start_2, p_rad_stop_2, delta_p_rad_2
        REAL p_rad_start_3, p_rad_stop_3, delta_p_rad_3
C
        REAL t_deg_start_1, t_deg_stop_1, delta_t_deg_1
        REAL t_deg_start_2, t_deg_stop_2, delta_t_deg_2
        REAL t_deg_start_3, t_deg_stop_3, delta_t_deg_3
        REAL p_deg_start_1, p_deg_stop_1, delta_p_deg_1
        REAL p_deg_start_2, p_deg_stop_2, delta_p_deg_2
        REAL p_deg_start_3, p_deg_stop_3, delta_p_deg_3
C
        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &          I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &          I_PDF_NDL,I_PDF_TRUNK

        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &                 I_PDF_NDL,I_PDF_TRUNK
C
        REAL WORK
C
c
        INCLUDE 'constants.include'
c
        WORK = PI/180.
c
        if((i_pdf_br_1.eq.1).or.(i_pdf_br_1.eq.2))then
c
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 27.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 6
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 30.0
            t_deg_stop_2  = 150.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 121
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 152.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 6
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if(i_pdf_br_1.eq.3)then
c
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 32.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 7
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 35.0
            t_deg_stop_2  = 145.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 111
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 147.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 7
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if((i_pdf_br_1.eq.4).or.(i_pdf_br_1.eq.5))then
c
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 42.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 9
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 45.0
            t_deg_stop_2  = 135.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 91
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 147.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 7
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if(i_pdf_br_1.eq.6)then
c
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 47.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 10
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 50.0
            t_deg_stop_2  = 130.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 81
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 132.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 10
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if((i_pdf_br_1.eq.7).or.(i_pdf_br_1.eq.8).or.
     &                                (i_pdf_br_1.eq.9))then
c
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 52.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 11
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 55.0
            t_deg_stop_2  = 125.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 71
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 127.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 11
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if(i_pdf_br_1.eq.10)then
c
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 57.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 12
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 60.0
            t_deg_stop_2  = 120.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 61
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 122.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 12
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if(i_pdf_br_1.eq.11)then
c
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 27.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 6
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 30.0
            t_deg_stop_2  = 90.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 61
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 92.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 18
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if(i_pdf_br_1.eq.12)then
c
            t_deg_start_1 = 1.0
            t_deg_stop_1  = 60.0
            delta_t_deg_1 = 1.0
            n_theta_1 = 60
c
            p_deg_start_1 = 0.5
            p_deg_stop_1  = 179.5
            delta_p_deg_1 = 1.0
            n_phi_1 = 180
c
            t_deg_start_2 = 62.5
            t_deg_stop_2  = 177.5
            delta_t_deg_2 = 5.0
            n_theta_2 = 24
c
            p_deg_start_2 = 2.5
            p_deg_stop_2  = 177.5
            delta_p_deg_2 = 5.0
            n_phi_2 = 36
c
            t_deg_start_3 = 0.0
            t_deg_stop_3  = 0.0
            delta_t_deg_3 = 0.0
            n_theta_3 = 0
c
            p_deg_start_3 = 0.0
            p_deg_stop_3  = 0.0
            delta_p_deg_3 = 0.0
            n_phi_3 = 0
c
        else if(i_pdf_br_1.eq.13)then
c
            t_deg_start_1 = 7.5
            t_deg_stop_1  = 57.5
            delta_t_deg_1 = 10.0
            n_theta_1 = 0
c
            p_deg_start_1 = 0.0
            p_deg_stop_1  = 360.0
            delta_p_deg_1 = 10.0
            n_phi_1 = 0
c
            t_deg_start_2 = 60.0
            t_deg_stop_2  = 120.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 61
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 122.5
            t_deg_stop_3  = 172.5
            delta_t_deg_3 = 10.0
            n_theta_3 = 0
c
            p_deg_start_3 = 0.0
            p_deg_stop_3  = 360.0
            delta_p_deg_3 = 10.0
            n_phi_3 = 0
c
        else if(i_pdf_br_1.eq.14)then
c
            t_deg_start_1 = 0.0
            t_deg_stop_1  = 30.0
            delta_t_deg_1 = 10.0
            n_theta_1 = 0
c
            p_deg_start_1 = 0.0
            p_deg_stop_1  = 360.0
            delta_p_deg_1 = 10.0
            n_phi_1 = 0
c
            t_deg_start_2 = 30.0
            t_deg_stop_2  = 90.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 61
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 90.
            t_deg_stop_3  = 180.0
            delta_t_deg_3 = 10.0
            n_theta_3 = 0
c
            p_deg_start_3 = 0.0
            p_deg_stop_3  = 360.0
            delta_p_deg_3 = 10.0
            n_phi_3 = 0
c
        else if(i_pdf_br_1.eq.17)then
c
            t_deg_start_1 =  1.0
            t_deg_stop_1  = 50.0
            delta_t_deg_1 =  1.0
            n_theta_1 = 50
c
            p_deg_start_1 =   0.5
            p_deg_stop_1  = 179.5
            delta_p_deg_1 =   1.0
            n_phi_1 = 180
c
            t_deg_start_2 = 55.0
            t_deg_stop_2  = 85.0
            delta_t_deg_2 =  5.0
            n_theta_2 = 7
c
            p_deg_start_2 =   5.0
            p_deg_stop_2  = 175.0
            delta_p_deg_2 =   5.0
            n_phi_2 = 35
c
            t_deg_start_3 = 90.
            t_deg_stop_3  = 180.0
            delta_t_deg_3 = 10.0
            n_theta_3 = 0
c
            p_deg_start_3 = 0.0
            p_deg_stop_3  = 360.0
            delta_p_deg_3 = 10.0
            n_phi_3 = 0

c
        else if(i_pdf_br_1.eq.18)then
c
            t_deg_start_1 = 5.0
            t_deg_stop_1  = 15.0
            delta_t_deg_1 = 10.0
            n_theta_1 = 2
c
            p_deg_start_1 =   5.0
            p_deg_stop_1  = 175.0
            delta_p_deg_1 =  10.0
            n_phi_1 = 18
c
            t_deg_start_2 = 20.0
            t_deg_stop_2  = 90.0
            delta_t_deg_2 =  1.0
            n_theta_2 = 71
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 90.
            t_deg_stop_3  = 180.0
            delta_t_deg_3 = 10.0
            n_theta_3 = 0
c
            p_deg_start_3 = 0.0
            p_deg_stop_3  = 360.0
            delta_p_deg_3 = 10.0
            n_phi_3 = 0
c
        else
c
            t_deg_start_1 = 2.5 
            t_deg_stop_1  = 177.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 35
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 0.0
            t_deg_stop_2  = 0.0
            delta_t_deg_2 = 0.0
            n_theta_2 = 0
c
            p_deg_start_2 = 0.0
            p_deg_stop_2  = 0.0
            delta_p_deg_2 = 0.0
            n_phi_2 = 0
c
            t_deg_start_3 = 0.0
            t_deg_stop_3  = 0.0
            delta_t_deg_3 = 0.0
            n_theta_3 = 0
c
            p_deg_start_3 = 0.0
            p_deg_stop_3  = 0.0
            delta_p_deg_3 = 0.0
            n_phi_3 = 0
c
        endif
c
        t_rad_start_1 =  WORK*t_deg_start_1
        t_rad_stop_1  =  WORK*t_deg_stop_1 
        delta_t_rad_1 =  WORK*delta_t_deg_1
                                     
        p_rad_start_1 =  WORK*p_deg_start_1
        p_rad_stop_1  =  WORK*p_deg_stop_1 
        delta_p_rad_1 =  WORK*delta_p_deg_1
                                     
        t_rad_start_2 =  WORK*t_deg_start_2
        t_rad_stop_2  =  WORK*t_deg_stop_2 
        delta_t_rad_2 =  WORK*delta_t_deg_2
                                     
        p_rad_start_2 =  WORK*p_deg_start_2
        p_rad_stop_2  =  WORK*p_deg_stop_2 
        delta_p_rad_2 =  WORK*delta_p_deg_2
                                     
        t_rad_start_3 =  WORK*t_deg_start_3
        t_rad_stop_3  =  WORK*t_deg_stop_3 
        delta_t_rad_3 =  WORK*delta_t_deg_3
                                     
        p_rad_start_3 =  WORK*p_deg_start_3
        p_rad_stop_3  =  WORK*p_deg_stop_3 
        delta_p_rad_3 =  WORK*delta_p_deg_3
c
        RETURN
        END
c
        SUBROUTINE PDF_BR2_SETUP(n_theta_1,t_rad_start_1,t_rad_stop_1,
     & delta_t_rad_1,n_theta_2,t_rad_start_2,t_rad_stop_2,delta_t_rad_2,
     &               n_theta_3,t_rad_start_3,t_rad_stop_3,delta_t_rad_3,
     &               n_phi_1,p_rad_start_1,p_rad_stop_1,delta_p_rad_1,
     &               n_phi_2,p_rad_start_2,p_rad_stop_2,delta_p_rad_2,
     &               n_phi_3,p_rad_start_3,p_rad_stop_3,delta_p_rad_3)
c
        save
C
        INTEGER n_theta_1,n_theta_2,n_theta_3, n_phi_1,n_phi_2,n_phi_3
        REAL t_rad_start_1, t_rad_stop_1, delta_t_rad_1
        REAL t_rad_start_2, t_rad_stop_2, delta_t_rad_2
        REAL t_rad_start_3, t_rad_stop_3, delta_t_rad_3
        REAL p_rad_start_1, p_rad_stop_1, delta_p_rad_1
        REAL p_rad_start_2, p_rad_stop_2, delta_p_rad_2
        REAL p_rad_start_3, p_rad_stop_3, delta_p_rad_3
C
        REAL t_deg_start_1, t_deg_stop_1, delta_t_deg_1
        REAL t_deg_start_2, t_deg_stop_2, delta_t_deg_2
        REAL t_deg_start_3, t_deg_stop_3, delta_t_deg_3
        REAL p_deg_start_1, p_deg_stop_1, delta_p_deg_1
        REAL p_deg_start_2, p_deg_stop_2, delta_p_deg_2
        REAL p_deg_start_3, p_deg_stop_3, delta_p_deg_3
C
        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &          I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &          I_PDF_NDL,I_PDF_TRUNK

        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &                 I_PDF_NDL,I_PDF_TRUNK

C
        REAL WORK
c
        INCLUDE 'constants.include'
c
        WORK = PI/180.
C
        if((i_pdf_br_2.eq.1).or.(i_pdf_br_2.eq.2))then
c
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 27.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 6
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 30.0
            t_deg_stop_2  = 150.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 121
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 152.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 6
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if(i_pdf_br_2.eq.3)then
c
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 32.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 7
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 35.0
            t_deg_stop_2  = 145.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 111
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 147.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 7
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if((i_pdf_br_2.eq.4).or.(i_pdf_br_2.eq.5))then
c
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 42.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 9
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 45.0
            t_deg_stop_2  = 135.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 91
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 147.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 7
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if(i_pdf_br_2.eq.6)then
c
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 47.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 10
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 50.0
            t_deg_stop_2  = 130.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 81
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 132.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 10
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if((i_pdf_br_2.eq.7).or.(i_pdf_br_2.eq.8).or.
     &                                (i_pdf_br_2.eq.9))then
c
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 52.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 11
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 55.0
            t_deg_stop_2  = 125.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 71
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 127.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 11
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if(i_pdf_br_2.eq.10)then
c
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 57.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 12
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 60.0
            t_deg_stop_2  = 120.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 61
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 122.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 12
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if(i_pdf_br_2.eq.11)then
c
            t_deg_start_1 = 2.5
            t_deg_stop_1  = 27.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 6
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 30.0
            t_deg_stop_2  = 90.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 61
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 92.5
            t_deg_stop_3  = 177.5
            delta_t_deg_3 = 5.0
            n_theta_3 = 18
c
            p_deg_start_3 = 2.5
            p_deg_stop_3  = 177.5
            delta_p_deg_3 = 5.0
            n_phi_3 = 36
c
        else if(i_pdf_br_2.eq.12)then
c
            t_deg_start_1 = 1.0
            t_deg_stop_1  = 60.0
            delta_t_deg_1 = 1.0
            n_theta_1 = 60
c
            p_deg_start_1 = 0.5
            p_deg_stop_1  = 179.5
            delta_p_deg_1 = 1.0
            n_phi_1 = 180
c
            t_deg_start_2 = 62.5
            t_deg_stop_2  = 177.5
            delta_t_deg_2 = 5.0
            n_theta_2 = 24
c
            p_deg_start_2 = 2.5
            p_deg_stop_2  = 177.5
            delta_p_deg_2 = 5.0
            n_phi_2 = 36
c
            t_deg_start_3 = 0.0
            t_deg_stop_3  = 0.0
            delta_t_deg_3 = 0.0
            n_theta_3 = 0
c
            p_deg_start_3 = 0.0
            p_deg_stop_3  = 0.0
            delta_p_deg_3 = 0.0
            n_phi_3 = 0
c
        else if(i_pdf_br_2.eq.13)then
c
            t_deg_start_1 = 7.5
            t_deg_stop_1  = 57.5
            delta_t_deg_1 = 10.0
            n_theta_1 = 0
c
            p_deg_start_1 = 0.0
            p_deg_stop_1  = 360.0
            delta_p_deg_1 = 10.0
            n_phi_1 = 0
c
            t_deg_start_2 = 60.0
            t_deg_stop_2  = 120.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 61
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 122.5
            t_deg_stop_3  = 172.5
            delta_t_deg_3 = 10.0
            n_theta_3 = 0
c
            p_deg_start_3 = 0.0
            p_deg_stop_3  = 360.0
            delta_p_deg_3 = 10.0
            n_phi_3 = 0
c
        else if(i_pdf_br_2.eq.14)then
c
            t_deg_start_1 = 0.0
            t_deg_stop_1  = 30.0
            delta_t_deg_1 = 10.0
            n_theta_1 = 0
c
            p_deg_start_1 = 0.0
            p_deg_stop_1  = 360.0
            delta_p_deg_1 = 10.0
            n_phi_1 = 0
c
            t_deg_start_2 = 30.0
            t_deg_stop_2  = 90.0
            delta_t_deg_2 = 1.0
            n_theta_2 = 61
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 90.
            t_deg_stop_3  = 180.0
            delta_t_deg_3 = 10.0
            n_theta_3 = 0
c
            p_deg_start_3 = 0.0
            p_deg_stop_3  = 360.0
            delta_p_deg_3 = 10.0
            n_phi_3 = 0
c
        else if(i_pdf_br_2.eq.17)then
c
            t_deg_start_1 =  1.0
            t_deg_stop_1  = 50.0
            delta_t_deg_1 =  1.0
            n_theta_1 = 50
c
            p_deg_start_1 =   0.5
            p_deg_stop_1  = 179.5
            delta_p_deg_1 =   1.0
            n_phi_1 = 180
c
            t_deg_start_2 = 55.0
            t_deg_stop_2  = 85.0
            delta_t_deg_2 =  5.0
            n_theta_2 = 7
c
            p_deg_start_2 =   5.0
            p_deg_stop_2  = 175.0
            delta_p_deg_2 =   5.0
            n_phi_2 = 35
c
            t_deg_start_3 = 90.
            t_deg_stop_3  = 180.0
            delta_t_deg_3 = 10.0
            n_theta_3 = 0
c
            p_deg_start_3 = 0.0
            p_deg_stop_3  = 360.0
            delta_p_deg_3 = 10.0
            n_phi_3 = 0
c
        else if(i_pdf_br_2.eq.18)then
c
            t_deg_start_1 = 5.0
            t_deg_stop_1  = 15.0
            delta_t_deg_1 = 10.0
            n_theta_1 = 2
c
            p_deg_start_1 =   5.0
            p_deg_stop_1  = 175.0
            delta_p_deg_1 =  10.0
            n_phi_1 = 18
c
            t_deg_start_2 = 20.0
            t_deg_stop_2  = 90.0
            delta_t_deg_2 =  1.0
            n_theta_2 = 71
c
            p_deg_start_2 = 0.5
            p_deg_stop_2  = 179.5
            delta_p_deg_2 = 1.0
            n_phi_2 = 180
c
            t_deg_start_3 = 90.
            t_deg_stop_3  = 180.0
            delta_t_deg_3 = 10.0
            n_theta_3 = 0
c
            p_deg_start_3 = 0.0
            p_deg_stop_3  = 360.0
            delta_p_deg_3 = 10.0
            n_phi_3 = 0
c
        else
c
            t_deg_start_1 = 2.5 
            t_deg_stop_1  = 177.5
            delta_t_deg_1 = 5.0
            n_theta_1 = 35
c
            p_deg_start_1 = 2.5
            p_deg_stop_1  = 177.5
            delta_p_deg_1 = 5.0
            n_phi_1 = 36
c
            t_deg_start_2 = 0.0
            t_deg_stop_2  = 0.0
            delta_t_deg_2 = 0.0
            n_theta_2 = 0
c
            p_deg_start_2 = 0.0
            p_deg_stop_2  = 0.0
            delta_p_deg_2 = 0.0
            n_phi_2 = 0
c
            t_deg_start_3 = 0.0
            t_deg_stop_3  = 0.0
            delta_t_deg_3 = 0.0
            n_theta_3 = 0
c
            p_deg_start_3 = 0.0
            p_deg_stop_3  = 0.0
            delta_p_deg_3 = 0.0
            n_phi_3 = 0
c
        endif
c
        t_rad_start_1 =  WORK*t_deg_start_1
        t_rad_stop_1  =  WORK*t_deg_stop_1 
        delta_t_rad_1 =  WORK*delta_t_deg_1
                                     
        p_rad_start_1 =  WORK*p_deg_start_1
        p_rad_stop_1  =  WORK*p_deg_stop_1 
        delta_p_rad_1 =  WORK*delta_p_deg_1
                                     
        t_rad_start_2 =  WORK*t_deg_start_2
        t_rad_stop_2  =  WORK*t_deg_stop_2 
        delta_t_rad_2 =  WORK*delta_t_deg_2
                                     
        p_rad_start_2 =  WORK*p_deg_start_2
        p_rad_stop_2  =  WORK*p_deg_stop_2 
        delta_p_rad_2 =  WORK*delta_p_deg_2
                                     
        t_rad_start_3 =  WORK*t_deg_start_3
        t_rad_stop_3  =  WORK*t_deg_stop_3 
        delta_t_rad_3 =  WORK*delta_t_deg_3
                                     
        p_rad_start_3 =  WORK*p_deg_start_3
        p_rad_stop_3  =  WORK*p_deg_stop_3 
        delta_p_rad_3 =  WORK*delta_p_deg_3
c
        RETURN
        END
c
        SUBROUTINE PDF_TR_SETUP(n_theta_1,t_rad_start_1,t_rad_stop_1,
     & delta_t_rad_1,n_theta_2,t_rad_start_2,t_rad_stop_2,delta_t_rad_2,
     &               n_theta_3,t_rad_start_3,t_rad_stop_3,delta_t_rad_3,
     &               n_phi_1,p_rad_start_1,p_rad_stop_1,delta_p_rad_1,
     &               n_phi_2,p_rad_start_2,p_rad_stop_2,delta_p_rad_2,
     &               n_phi_3,p_rad_start_3,p_rad_stop_3,delta_p_rad_3)
c
        save
C
        INTEGER n_theta_1,n_theta_2,n_theta_3, n_phi_1,n_phi_2,n_phi_3
        REAL t_rad_start_1, t_rad_stop_1, delta_t_rad_1
        REAL t_rad_start_2, t_rad_stop_2, delta_t_rad_2
        REAL t_rad_start_3, t_rad_stop_3, delta_t_rad_3
        REAL p_rad_start_1, p_rad_stop_1, delta_p_rad_1
        REAL p_rad_start_2, p_rad_stop_2, delta_p_rad_2
        REAL p_rad_start_3, p_rad_stop_3, delta_p_rad_3
C
        REAL t_deg_start_1, t_deg_stop_1, delta_t_deg_1
        REAL t_deg_start_2, t_deg_stop_2, delta_t_deg_2
        REAL t_deg_start_3, t_deg_stop_3, delta_t_deg_3
        REAL p_deg_start_1, p_deg_stop_1, delta_p_deg_1
        REAL p_deg_start_2, p_deg_stop_2, delta_p_deg_2
        REAL p_deg_start_3, p_deg_stop_3, delta_p_deg_3
C
        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &          I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &          I_PDF_NDL,I_PDF_TRUNK

        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &                 I_PDF_NDL,I_PDF_TRUNK
C
        REAL WORK
c
        INCLUDE 'constants.include'
c
        WORK = PI/180.
C
        t_deg_start_1 = 0.0 
        t_deg_stop_1  = 5.0
        delta_t_deg_1 = 5.0
        n_theta_1 = 1
c
        p_deg_start_1 = 0.0
        p_deg_stop_1  = 360.0
        delta_p_deg_1 = 360.0
        n_phi_1 = 1
c
c
        t_deg_start_2 = t_deg_stop_1
        t_deg_stop_2  = 180.0
        delta_t_deg_2 = 5.0
        n_theta_2 = 35
c
        p_deg_start_2 = 2.5
        p_deg_stop_2  = 177.5
        delta_p_deg_2 = 5.0
        n_phi_2 = 36
c
c
        t_deg_start_3 = 0.0
        t_deg_stop_3  = 0.0
        delta_t_deg_3 = 0.0
        n_theta_3 = 0
c
        p_deg_start_3 = 0.0
        p_deg_stop_3  = 0.0
        delta_p_deg_3 = 0.0
        n_phi_3 = 0
C
C
        t_rad_start_1 =  WORK*t_deg_start_1
        t_rad_stop_1  =  WORK*t_deg_stop_1 
        delta_t_rad_1 =  WORK*delta_t_deg_1
                                     
        p_rad_start_1 =  WORK*p_deg_start_1
        p_rad_stop_1  =  WORK*p_deg_stop_1 
        delta_p_rad_1 =  WORK*delta_p_deg_1
                                     
        t_rad_start_2 =  WORK*t_deg_start_2
        t_rad_stop_2  =  WORK*t_deg_stop_2 
        delta_t_rad_2 =  WORK*delta_t_deg_2
                                     
        p_rad_start_2 =  WORK*p_deg_start_2
        p_rad_stop_2  =  WORK*p_deg_stop_2 
        delta_p_rad_2 =  WORK*delta_p_deg_2
                                     
        t_rad_start_3 =  WORK*t_deg_start_3
        t_rad_stop_3  =  WORK*t_deg_stop_3 
        delta_t_rad_3 =  WORK*delta_t_deg_3
                                     
        p_rad_start_3 =  WORK*p_deg_start_3
        p_rad_stop_3  =  WORK*p_deg_stop_3 
        delta_p_rad_3 =  WORK*delta_p_deg_3
c
        RETURN
        END
