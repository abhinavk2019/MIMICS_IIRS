        SUBROUTINE NEEDLE_PHASE(P_NEEDLE)
        save
c
        include 'parameters.include'
c
        INTEGER II,I,L,K,JJ
                               
        REAL LENGTHCM, diamcm, K0A
        REAL WORK, NORM ,PDF, PDF_NDL
        REAL THETAI,PHII,THETAS,PHIS,THETAN,PHIN
        REAL P_NEEDLE(4,4,4),STOKES_MAT(4,4)
        REAL P_NEEDLE_tmp_1(4,4,4), P_NEEDLE_tmp_2(4,4,4)
        REAL P_NEEDLE_tmp_3(4,4,4)
C
        REAL MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
        REAL FREQ_HERTZ, WAVELENGTH, k0
        REAL THETA, CTHETA, STHETA
        REAL DENSITY, CROWN_HGHT, TRUNK_HGHT
        REAL EXP_KAPPA_NDL_P(4,4),EXP_KAPPA_NDL_M(4,4)
        real lambdacm, lengthresh
        logical rayleigh
C
        COMPLEX EPSILONR(N_EPS),EPSILONRC(N_EPS)
C
        COMPLEX J, CWORK, EPSR
        COMPLEX SMAT_AVG(2,2,2),M_NEEDLE_P(2,2),M_NEEDLE_M(2,2)
        COMPLEX SMAT_AVG_tmp_1(2,2,2), SMAT_AVG_tmp_2(2,2,2)
        COMPLEX SMAT_AVG_tmp_3(2,2,2)
        COMPLEX Smat(2,2)
c
        integer n_theta_1, n_theta_2, n_theta_3
        integer n_phi_1, n_phi_2, n_phi_3
        real t_rad_start_1,t_rad_stop_1,delta_t_rad_1
        real t_rad_start_2,t_rad_stop_2,delta_t_rad_2
        real t_rad_start_3,t_rad_stop_3,delta_t_rad_3
        real p_rad_start_1,p_rad_stop_1,delta_p_rad_1
        real p_rad_start_2,p_rad_stop_2,delta_p_rad_2
        real p_rad_start_3,p_rad_stop_3,delta_p_rad_3
        real delta_p1_t1, delta_p2_t2, delta_p3_t3
C
        COMMON /R_NDL/ MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
        COMMON /RADAR_FREQ/ FREQ_HERTZ, WAVELENGTH, k0
        COMMON /RADAR_ANGLE/ THETA, CTHETA, STHETA
        COMMON /R_CANOPY/ DENSITY, CROWN_HGHT, TRUNK_HGHT
        COMMON /C_DIELECTRIC/ EPSILONR,EPSILONRC
C
        common /ndldatar/ diamcm,lengthcm,k0A
        common /ndldatac/ epsr
        common /ndlscatmat/ Smat
c
        COMMON /NDL_EXT/ EXP_KAPPA_NDL_p, EXP_KAPPA_NDL_m
        COMMON /NDL_M_MAT/ M_NEEDLE_P, M_NEEDLE_M
c
        include 'constants.include'
c
        J = CMPLX(0.0,1.0)
        EPSR = EPSILONRC(3)
        diamcm = NDL_DIAM
        lengthcm = NDL_LNG
        k0a = k0
        lambdacm = 100.*(2.*pi/k0)
        lengthresh = 1.5*lambdacm
        if(lengthcm.lt.lengthresh)then
           rayleigh = .true.
        else
           rayleigh = .false.
        endif
        print*,' needle rayleigh = ',rayleigh
c
        DO i=1,4
         DO k=1,4
          DO l=1,4
             P_NEEDLE(l,k,i) = 0.0
             P_NEEDLE_tmp_1(l,k,i) = 0.0
             P_NEEDLE_tmp_2(l,k,i) = 0.0
             P_NEEDLE_tmp_3(l,k,i) = 0.0
	  END DO
	 END DO
	END DO
C
        DO i=1,2
         DO k=1,2
          DO l=1,2
             Smat_avg(l,k,i) = CMPLX(0.0,0.0)
             Smat_avg_tmp_1(l,k,i) = CMPLX(0.0,0.0)
             Smat_avg_tmp_2(l,k,i) = CMPLX(0.0,0.0)
             Smat_avg_tmp_3(l,k,i) = CMPLX(0.0,0.0)
	  END DO
	 END DO
	END DO
c
        NORM = 2.0*PI
        PHII = 0.0
c
        call pdf_ndl_setup(n_theta_1,t_rad_start_1,t_rad_stop_1,
     &       delta_t_rad_1,n_theta_2,t_rad_start_2,t_rad_stop_2,
     &       delta_t_rad_2,n_theta_3,t_rad_start_3,t_rad_stop_3,
     &       delta_t_rad_3,n_phi_1,p_rad_start_1,p_rad_stop_1,
     &  delta_p_rad_1,n_phi_2,p_rad_start_2,p_rad_stop_2,delta_p_rad_2,
     &       n_phi_3,p_rad_start_3,p_rad_stop_3,delta_p_rad_3)
c
        delta_p1_t1 = delta_p_rad_1*delta_t_rad_1
        delta_p2_t2 = delta_p_rad_2*delta_t_rad_2
        delta_p3_t3 = delta_p_rad_3*delta_t_rad_3
c
      IF((n_theta_1.ne.0).and.(n_phi_1.ne.0))THEN
c
      do 150 jj=1,n_theta_1
       thetan = t_rad_start_1 + delta_t_rad_1*float(jj-1)
       PDF = PDF_NDL(THETAN)
       IF(PDF.GT.0.0)THEN
        print*,' thetan, pdf = ',thetan,pdf
c
        DO 145 II = 1,n_phi_1
         PHIN = p_rad_start_1 + delta_p_rad_1*float(ii-1)
c
        THETAI = THETA
        THETAS = PI - THETA
        PHIS   = PI + PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_1(L,K,1) =
     &               P_NEEDLE_tmp_1(L,K,1) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)

        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_1(L,K,1) =
     &               P_NEEDLE_tmp_1(L,K,1) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

c
        THETAI = PI - THETA
        THETAS = PI - THETA
        PHIS   = PI + PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_1(L,K,2) =
     &       P_NEEDLE_tmp_1(L,K,2) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_1(L,K,2) =
     &       P_NEEDLE_tmp_1(L,K,2) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO
c
        THETAI = THETA
        THETAS = THETA
        PHIS = PI + PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_1(L,K,3) =
     &             P_NEEDLE_tmp_1(L,K,3) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_1(L,K,3) =
     &             P_NEEDLE_tmp_1(L,K,3) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO
c
        THETAI = PI - THETA
        THETAS = THETA
        PHIS   = PI + PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_1(L,K,4) =
     &            P_NEEDLE_tmp_1(L,K,4) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_1(L,K,4) =
     &            P_NEEDLE_tmp_1(L,K,4) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

c
        THETAI = THETA
        THETAS = THETA
        PHIS   = PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        DO K=1,2
         DO L=1,2
          Smat_avg_tmp_1(L,K,1) =
     &        Smat_avg_tmp_1(L,K,1) + PDF*Smat(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        DO K=1,2
         DO L=1,2
          Smat_avg_tmp_1(L,K,1) =
     &        Smat_avg_tmp_1(L,K,1) + PDF*Smat(L,K)
         ENDDO
        ENDDO

c
        THETAI = PI - THETA
        THETAS = PI - THETA
        PHIS   = PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        DO K=1,2
         DO L=1,2
          Smat_avg_tmp_1(L,K,2) =
     &          Smat_avg_tmp_1(L,K,2) + PDF*Smat(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        DO K=1,2
         DO L=1,2
          Smat_avg_tmp_1(L,K,2) =
     &          Smat_avg_tmp_1(L,K,2) + PDF*Smat(L,K)
         ENDDO
        ENDDO
C
145   CONTINUE
C
      ENDIF
150   CONTINUE
      endif
c
      IF((n_theta_2.ne.0).and.(n_phi_2.ne.0))THEN
c
      do 170 jj=1,n_theta_2
       thetan = t_rad_start_2 + delta_t_rad_2*float(jj-1)
       PDF = PDF_NDL(THETAN)
       IF(PDF.GT.0.0)THEN
        print*,' thetan, pdf = ',thetan,pdf
c
        DO 165 II = 1,n_phi_2
         phin = p_rad_start_2 + delta_p_rad_2*float(ii-1)
c
        THETAI = THETA
        THETAS = PI - THETA
        PHIS   = PI + PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_2(L,K,1) =
     &               P_NEEDLE_tmp_2(L,K,1) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_2(L,K,1) =
     &               P_NEEDLE_tmp_2(L,K,1) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO
c
        THETAI = PI - THETA
        THETAS = PI - THETA
        PHIS   = PI + PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_2(L,K,2) =
     &       P_NEEDLE_tmp_2(L,K,2) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_2(L,K,2) =
     &       P_NEEDLE_tmp_2(L,K,2) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO
c
        THETAI = THETA
        THETAS = THETA
        PHIS = PI + PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_2(L,K,3) =
     &             P_NEEDLE_tmp_2(L,K,3) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_2(L,K,3) =
     &             P_NEEDLE_tmp_2(L,K,3) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

c
        THETAI = PI - THETA
        THETAS = THETA
        PHIS   = PI + PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_2(L,K,4) =
     &            P_NEEDLE_tmp_2(L,K,4) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_2(L,K,4) =
     &            P_NEEDLE_tmp_2(L,K,4) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO
c
        THETAI = THETA
        THETAS = THETA
        PHIS   = PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        DO K=1,2
         DO L=1,2
          Smat_avg_tmp_2(L,K,1) =
     &        Smat_avg_tmp_2(L,K,1) + PDF*Smat(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        DO K=1,2
         DO L=1,2
          Smat_avg_tmp_2(L,K,1) =
     &        Smat_avg_tmp_2(L,K,1) + PDF*Smat(L,K)
         ENDDO
        ENDDO

c
        THETAI = PI - THETA
        THETAS = PI - THETA
        PHIS   = PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        DO K=1,2
         DO L=1,2
          Smat_avg_tmp_2(L,K,2) =
     &          Smat_avg_tmp_2(L,K,2) + PDF*Smat(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        DO K=1,2
         DO L=1,2
          Smat_avg_tmp_2(L,K,2) =
     &          Smat_avg_tmp_2(L,K,2) + PDF*Smat(L,K)
         ENDDO
        ENDDO
C
165   CONTINUE
C
      ENDIF
170   CONTINUE
      endif
c
      IF((n_theta_2.ne.0).and.(n_phi_2.ne.0))THEN
      do 200 jj=1,n_theta_3
       thetan = t_rad_start_3 + delta_t_rad_3*float(jj-1)
       PDF = PDF_NDL(THETAN)
       IF(PDF.GT.0.0)THEN
        print*,' thetan, pdf = ',thetan,pdf
c
        DO 195 II = 1,n_phi_3
         phin = p_rad_start_3 + delta_p_rad_3*float(ii-1)
c
        THETAI = THETA
        THETAS = PI - THETA
        PHIS   = PI + PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_3(L,K,1) =
     &               P_NEEDLE_tmp_3(L,K,1) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_3(L,K,1) =
     &               P_NEEDLE_tmp_3(L,K,1) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO
c
        THETAI = PI - THETA
        THETAS = PI - THETA
        PHIS   = PI + PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_3(L,K,2) =
     &       P_NEEDLE_tmp_3(L,K,2) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_3(L,K,2) =
     &       P_NEEDLE_tmp_3(L,K,2) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO
c
        THETAI = THETA
        THETAS = THETA
        PHIS = PI + PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_3(L,K,3) =
     &             P_NEEDLE_tmp_3(L,K,3) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_3(L,K,3) =
     &             P_NEEDLE_tmp_3(L,K,3) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO
c
        THETAI = PI - THETA
        THETAS = THETA
        PHIS   = PI + PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        CALL STOKES_SUB(SMAT,STOKES_MAT)
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_3(L,K,4) =
     &            P_NEEDLE_tmp_3(L,K,4) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        CALL STOKES_SUB(SMAT,STOKES_MAT)
        DO K=1,4
         DO L=1,4
          P_NEEDLE_tmp_3(L,K,4) =
     &            P_NEEDLE_tmp_3(L,K,4) + PDF*STOKES_MAT(L,K)
         ENDDO
        ENDDO
c
        THETAI = THETA
        THETAS = THETA
        PHIS   = PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        DO K=1,2
         DO L=1,2
          Smat_avg_tmp_3(L,K,1) =
     &        Smat_avg_tmp_3(L,K,1) + PDF*Smat(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        DO K=1,2
         DO L=1,2
          Smat_avg_tmp_3(L,K,1) =
     &        Smat_avg_tmp_3(L,K,1) + PDF*Smat(L,K)
         ENDDO
        ENDDO
c
        THETAI = PI - THETA
        THETAS = PI - THETA
        PHIS   = PHII
C
        if(rayleigh)then
          CALL NEEDLE_RLGH_ELPSD_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetan,phin)
        else
          call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetan,phin,diamcm,lengthcm,k0,epsr,Smat)
        endif

        DO K=1,2
         DO L=1,2
          Smat_avg_tmp_3(L,K,2) =
     &          Smat_avg_tmp_3(L,K,2) + PDF*Smat(L,K)
         ENDDO
        ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

        DO K=1,2
         DO L=1,2
          Smat_avg_tmp_3(L,K,2) =
     &          Smat_avg_tmp_3(L,K,2) + PDF*Smat(L,K)
         ENDDO
        ENDDO
C
195   CONTINUE
C
      ENDIF
200   CONTINUE
      endif
c
        WORK = 1.0/NORM
C
        DO i=1,4
         DO k=1,4
          DO l=1,4
           P_NEEDLE(l,k,i) = WORK*(delta_p1_t1*P_NEEDLE_tmp_1(l,k,i) +
     &                             delta_p2_t2*P_NEEDLE_tmp_2(l,k,i) +
     &                             delta_p3_t3*P_NEEDLE_tmp_3(l,k,i))
	  END DO
	 END DO
	END DO
C
        CWORK = J*2.0*PI/(NORM*k0)
C
        DO k=1,2
         DO l=1,2
           M_NEEDLE_P(l,k) = CWORK*(delta_p1_t1*Smat_avg_tmp_1(l,k,1) +
     &                              delta_p2_t2*Smat_avg_tmp_2(l,k,1) +
     &                              delta_p3_t3*Smat_avg_tmp_3(l,k,1))
	 END DO
	END DO
        DO k=1,2
         DO l=1,2
           M_NEEDLE_M(l,k) = CWORK*(delta_p1_t1*Smat_avg_tmp_1(l,k,2) +
     &                              delta_p2_t2*Smat_avg_tmp_2(l,k,2) +
     &                              delta_p3_t3*Smat_avg_tmp_3(l,k,2))
	 END DO
	END DO
C
        RETURN
        END
