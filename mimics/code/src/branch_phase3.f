        SUBROUTINE BRANCH_PHASE3(P_BRANCH_3)
        save
c
        include 'parameters.include'
c
        REAL FREQ_HERTZ, WAVELENGTH, k0
        REAL THETA, CTHETA, STHETA
        REAL MG_BR3, RHO_BR3, BR3_DENS, BR3_DIAM, BR3_LNG
        REAL DENSITY, CROWN_HGHT, TRUNK_HGHT
        COMPLEX EPSILONR(N_EPS),EPSILONRC(N_EPS)

        COMMON /R_CONSTANTS/ PI,MU_0,EPS_0,LIGHT
        COMMON /RADAR_FREQ/ FREQ_HERTZ, WAVELENGTH, k0
        COMMON /RADAR_ANGLE/ THETA, CTHETA, STHETA

        COMMON /R_BR3/ MG_BR3, RHO_BR3, BR3_DENS, BR3_DIAM, BR3_LNG
        COMMON /R_CANOPY/ DENSITY, CROWN_HGHT, TRUNK_HGHT
        COMMON /C_DIELECTRIC/ EPSILONR,EPSILONRC

        INTEGER II,I,L,K,JJ
C
        REAL diamcm, lengthm, k0a
        REAL THETAI,PHII,THETAS,PHIS,THETAC,PHIC
        REAL P_BRANCH_3(4,4,4),STOKES_MAT(4,4)
        REAL P_BRANCH_tmp_1(4,4,4), P_BRANCH_tmp_2(4,4,4)
        REAL P_BRANCH_tmp_3(4,4,4)
        REAL EXP_KAPPA_BR_P3(4,4),EXP_KAPPA_BR_M3(4,4)
        REAL PDF, PDF_BR3, NORM, WORK
C
        COMPLEX J, CWORK
        COMPLEX SMAT_AVG(2,2,2),M_BRANCH_P3(2,2),M_BRANCH_M3(2,2)
        COMPLEX SMAT_AVG_tmp_1(2,2,2), SMAT_AVG_tmp_2(2,2,2)
        COMPLEX SMAT_AVG_tmp_3(2,2,2), smat_tmp1(2,2), smat_tmp2(2,2)
        COMPLEX epsr
        COMPLEX Smat(2,2)

        real lengthcm, lengthresh, diamthresh, lambdacm
        logical resonant
C   
        LOGICAL LFLAG
C
        integer n_theta_1, n_theta_2, n_theta_3
        integer n_phi_1, n_phi_2, n_phi_3
c
        real t_rad_start_1, t_rad_stop_1, delta_t_rad_1
        real t_rad_start_2, t_rad_stop_2, delta_t_rad_2
        real t_rad_start_3, t_rad_stop_3, delta_t_rad_3
        real p_rad_start_1, p_rad_stop_1, delta_p_rad_1
        real p_rad_start_2, p_rad_stop_2, delta_p_rad_2
        real p_rad_start_3, p_rad_stop_3, delta_p_rad_3
        real delta_p1_t1, delta_p2_t2, delta_p3_t3
c
        common /rcyldata/ diamcm,lengthm,k0a
        common /ccyldata/ epsr
        common /cylscatmat/ Smat
c
        COMMON /BR3_EXT/ EXP_KAPPA_BR_P3, EXP_KAPPA_BR_M3
        COMMON /BR3_M_MAT/ M_BRANCH_P3, M_BRANCH_M3
c
        include 'constants.include'
c
        J = CMPLX(0.0,1.0)
        EPSR = EPSILONRC(9)
        diamcm = BR3_DIAM
        lengthm = BR3_LNG
        k0a = k0

        lengthcm = 100.*lengthm
        lambdacm = 100.*(2.*pi/k0)
        lengthresh = 2.0*lambdacm
        diamthresh = 0.1*lambdacm/(2.*pi*cabs(csqrt(epsr)))

        if((diamcm.le.diamthresh).and.(lengthcm.ge.lengthresh))then
            resonant = .false.
        else
            resonant = .true.
        endif

        print*,' 3rd branch, resonant = ',resonant
c
        DO i=1,4
         DO k=1,4
          DO l=1,4
             P_branch_3(l,k,i) = 0.0
             P_branch_tmp_1(l,k,i) = 0.0
             P_branch_tmp_2(l,k,i) = 0.0
             P_branch_tmp_3(l,k,i) = 0.0
	  END DO
	 END DO
	END DO
c
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
C
        call pdf_br3_setup(n_theta_1,t_rad_start_1,t_rad_stop_1,
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
       thetac = t_rad_start_1 + delta_t_rad_1*float(jj-1)
       pdf = pdf_br3(thetac)
       IF(PDF.GT.0.0)THEN
        print*,' branch 3a thetac,pdf = ',thetac,pdf
c
        DO 145 II = 1,n_phi_1
         PHIC = p_rad_start_1 + delta_p_rad_1*float(ii-1)
c
         THETAI = THETA
         THETAS = PI - THETA
         PHIS   = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)

         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_1(L,K,1) =
     &          P_BRANCH_tmp_1(L,K,1) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)

         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_1(L,K,1) =
     &          P_BRANCH_tmp_1(L,K,1) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
c
         THETAI = PI - THETA
         THETAS = PI - THETA
         PHIS   = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
            smat_tmp1(1,1) = smat(1,1)
            smat_tmp1(1,2) = -smat(2,1)
            smat_tmp1(2,1) = -smat(1,2)
            smat_tmp1(2,2) = smat(2,2)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_1(L,K,2) =
     &          P_BRANCH_tmp_1(L,K,2) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_1(L,K,2) =
     &          P_BRANCH_tmp_1(L,K,2) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
c
         THETAI = THETA
         THETAS = THETA
         PHIS = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
            smat_tmp2(1,1) =  smat(1,1)
            smat_tmp2(1,2) = -smat(2,1)
            smat_tmp2(2,1) = -smat(1,2)
            smat_tmp2(2,2) =  smat(2,2)
          else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_1(L,K,3) =
     &         P_BRANCH_tmp_1(L,K,3) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_1(L,K,3) =
     &         P_BRANCH_tmp_1(L,K,3) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
c
         if(resonant)then
           CALL STOKES_SUB(smat_tmp1,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_1(L,K,3) =
     &         P_BRANCH_tmp_1(L,K,3) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO
           smat_tmp1(1,2) = -smat_tmp1(1,2)
           smat_tmp1(2,1) = -smat_tmp1(2,1)
           CALL STOKES_SUB(smat_tmp1,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_1(L,K,3) =
     &         P_BRANCH_tmp_1(L,K,3) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO

           CALL STOKES_SUB(smat_tmp2,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_1(L,K,2) =
     &         P_BRANCH_tmp_1(L,K,2) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO
           smat_tmp2(1,2) = -smat_tmp2(1,2)
           smat_tmp2(2,1) = -smat_tmp2(2,1)
           CALL STOKES_SUB(smat_tmp2,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_1(L,K,2) =
     &         P_BRANCH_tmp_1(L,K,2) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO
         endif
c
         THETAI = PI - THETA
         THETAS = THETA
         PHIS   = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_1(L,K,4) =
     &         P_BRANCH_tmp_1(L,K,4) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_1(L,K,4) =
     &         P_BRANCH_tmp_1(L,K,4) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
c
         THETAI = THETA
         THETAS = THETA
         PHIS   = PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         DO K=1,2
          DO L=1,2
           Smat_avg_tmp_1(L,K,1) =
     &           Smat_avg_tmp_1(L,K,1) + PDF*Smat(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         DO K=1,2
          DO L=1,2
           Smat_avg_tmp_1(L,K,1) =
     &           Smat_avg_tmp_1(L,K,1) + PDF*Smat(L,K)
          ENDDO
         ENDDO
c
         THETAI = PI - THETA
         THETAS = PI - THETA
         PHIS   = PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
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
c
      endif
c
      IF((n_theta_2.ne.0).and.(n_phi_2.ne.0))THEN
c
      do 170 jj=1,n_theta_2
       thetac = t_rad_start_2 + delta_t_rad_2*float(jj-1)
       pdf = pdf_br3(thetac)
       IF(PDF.GT.0.0)THEN
        print*,' branch 3b thetac,pdf = ',thetac,pdf
c
        DO 165 II = 1,n_phi_2
         PHIC = p_rad_start_2 + delta_p_rad_2*float(ii-1)
c
         THETAI = THETA
         THETAS = PI - THETA
         PHIS   = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)

         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_2(L,K,1) =
     &          P_BRANCH_tmp_2(L,K,1) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)

         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_2(L,K,1) =
     &          P_BRANCH_tmp_2(L,K,1) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

c
         THETAI = PI - THETA
         THETAS = PI - THETA
         PHIS   = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
            smat_tmp1(1,1) =  smat(1,1)
            smat_tmp1(1,2) = -smat(2,1)
            smat_tmp1(2,1) = -smat(1,2)
            smat_tmp1(2,2) =  smat(2,2)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_2(L,K,2) =
     &          P_BRANCH_tmp_2(L,K,2) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_2(L,K,2) =
     &          P_BRANCH_tmp_2(L,K,2) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
c
         THETAI = THETA
         THETAS = THETA
         PHIS = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
            smat_tmp2(1,1) =  smat(1,1)
            smat_tmp2(1,2) = -smat(2,1)
            smat_tmp2(2,1) = -smat(1,2)
            smat_tmp2(2,2) =  smat(2,2)
          else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_2(L,K,3) =
     &         P_BRANCH_tmp_2(L,K,3) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_2(L,K,3) =
     &         P_BRANCH_tmp_2(L,K,3) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
c
         if(resonant)then
           CALL STOKES_SUB(smat_tmp1,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_2(L,K,3) =
     &         P_BRANCH_tmp_2(L,K,3) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO
           smat_tmp1(1,2) = -smat_tmp1(1,2)
           smat_tmp1(2,1) = -smat_tmp1(2,1)
           CALL STOKES_SUB(smat_tmp1,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_2(L,K,3) =
     &         P_BRANCH_tmp_2(L,K,3) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO

           CALL STOKES_SUB(smat_tmp2,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_2(L,K,2) =
     &         P_BRANCH_tmp_2(L,K,2) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO
           smat_tmp2(1,2) = -smat_tmp2(1,2)
           smat_tmp2(2,1) = -smat_tmp2(2,1)
           CALL STOKES_SUB(smat_tmp2,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_2(L,K,2) =
     &         P_BRANCH_tmp_2(L,K,2) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO
         endif
c
         THETAI = PI - THETA
         THETAS = THETA
         PHIS   = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_2(L,K,4) =
     &         P_BRANCH_tmp_2(L,K,4) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_2(L,K,4) =
     &         P_BRANCH_tmp_2(L,K,4) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
c
         THETAI = THETA
         THETAS = THETA
         PHIS   = PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         DO K=1,2
          DO L=1,2
           Smat_avg_tmp_2(L,K,1) =
     &           Smat_avg_tmp_2(L,K,1) + PDF*Smat(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         DO K=1,2
          DO L=1,2
           Smat_avg_tmp_2(L,K,1) =
     &           Smat_avg_tmp_2(L,K,1) + PDF*Smat(L,K)
          ENDDO
         ENDDO
c
         THETAI = PI - THETA
         THETAS = PI - THETA
         PHIS   = PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
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
c
      endif
c
      IF((n_theta_3.ne.0).and.(n_phi_3.ne.0))THEN
c
      do 200 jj=1,n_theta_3
       thetac = t_rad_start_3 + delta_t_rad_3*float(jj-1)
       pdf = pdf_br3(thetac)
       IF(PDF.GT.0.0)THEN
        print*,' branch 3c thetac,pdf = ',thetac,pdf
c
        DO 195 II = 1,n_phi_3
         PHIC = p_rad_start_3 + delta_p_rad_3*float(ii-1)
c
         THETAI = THETA
         THETAS = PI - THETA
         PHIS   = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)

         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_3(L,K,1) =
     &          P_BRANCH_tmp_3(L,K,1) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)

         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_3(L,K,1) =
     &          P_BRANCH_tmp_3(L,K,1) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
c
         THETAI = PI - THETA
         THETAS = PI - THETA
         PHIS   = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
            smat_tmp1(1,1) =  smat(1,1)
            smat_tmp1(1,2) = -smat(2,1)
            smat_tmp1(2,1) = -smat(1,2)
            smat_tmp1(2,2) =  smat(2,2)
          else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_3(L,K,2) =
     &          P_BRANCH_tmp_3(L,K,2) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_3(L,K,2) =
     &          P_BRANCH_tmp_3(L,K,2) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
c
         THETAI = THETA
         THETAS = THETA
         PHIS = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
            smat_tmp2(1,1) =  smat(1,1)
            smat_tmp2(1,2) = -smat(2,1)
            smat_tmp2(2,1) = -smat(1,2)
            smat_tmp2(2,2) =  smat(2,2)
          else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_3(L,K,3) =
     &         P_BRANCH_tmp_3(L,K,3) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)
C         
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_3(L,K,3) =
     &         P_BRANCH_tmp_3(L,K,3) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
c
         if(resonant)then
           CALL STOKES_SUB(smat_tmp1,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_3(L,K,3) =
     &         P_BRANCH_tmp_3(L,K,3) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO
           smat_tmp1(1,2) = -smat_tmp1(1,2)
           smat_tmp1(2,1) = -smat_tmp1(2,1)
           CALL STOKES_SUB(smat_tmp1,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_3(L,K,3) =
     &         P_BRANCH_tmp_3(L,K,3) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO

           CALL STOKES_SUB(smat_tmp2,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_3(L,K,2) =
     &         P_BRANCH_tmp_3(L,K,2) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO
           smat_tmp2(1,2) = -smat_tmp2(1,2)
           smat_tmp2(2,1) = -smat_tmp2(2,1)
           CALL STOKES_SUB(smat_tmp2,STOKES_MAT)
           DO K=1,4
            DO L=1,4
             P_BRANCH_tmp_3(L,K,2) =
     &         P_BRANCH_tmp_3(L,K,2) + PDF*STOKES_MAT(L,K)
            ENDDO
           ENDDO
         endif
c
         THETAI = PI - THETA
         THETAS = THETA
         PHIS   = PI + PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         CALL STOKES_SUB(SMAT,STOKES_MAT)
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_3(L,K,4) =
     &         P_BRANCH_tmp_3(L,K,4) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         CALL STOKES_SUB(SMAT,STOKES_MAT)
         DO K=1,4
          DO L=1,4
           P_BRANCH_tmp_3(L,K,4) =
     &         P_BRANCH_tmp_3(L,K,4) + PDF*STOKES_MAT(L,K)
          ENDDO
         ENDDO
c
         THETAI = THETA
         THETAS = THETA
         PHIS   = PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
         endif

         DO K=1,2
          DO L=1,2
           Smat_avg_tmp_3(L,K,1) =
     &           Smat_avg_tmp_3(L,K,1) + PDF*Smat(L,K)
          ENDDO
         ENDDO

        smat(1,2) = -smat(1,2)
        smat(2,1) = -smat(2,1)

         DO K=1,2
          DO L=1,2
           Smat_avg_tmp_3(L,K,1) =
     &           Smat_avg_tmp_3(L,K,1) + PDF*Smat(L,K)
          ENDDO
         ENDDO
c
         THETAI = PI - THETA
         THETAS = PI - THETA
         PHIS   = PHII
C
         if(resonant)then
             CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
         else
             call long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,smat)
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
c
      endif
c
        WORK = 1.0/NORM
C
        DO i=1,4
         DO k=1,4
          DO l=1,4
           P_branch_3(l,k,i) = WORK*(delta_p1_t1*P_branch_tmp_1(l,k,i) +
     &                               delta_p2_t2*P_branch_tmp_2(l,k,i) +
     &                               delta_p3_t3*P_branch_tmp_3(l,k,i))
	  END DO
	 END DO
	END DO
C
        CWORK = J*2.0*PI/(NORM*k0)
C
        DO k=1,2
         DO l=1,2
          M_BRANCH_P3(l,k) = CWORK*(delta_p1_t1*Smat_avg_tmp_1(l,k,1) +
     &                              delta_p2_t2*Smat_avg_tmp_2(l,k,1) +
     &                              delta_p3_t3*Smat_avg_tmp_3(l,k,1))
	 END DO
	END DO
        DO k=1,2
         DO l=1,2
          M_BRANCH_M3(l,k) = CWORK*(delta_p1_t1*Smat_avg_tmp_1(l,k,2) +
     &                              delta_p2_t2*Smat_avg_tmp_2(l,k,2) +
     &                              delta_p3_t3*Smat_avg_tmp_3(l,k,2))

	 END DO
	END DO
c
        if(resonant)then
         DO k=1,4
          DO l=1,4
           P_branch_3(l,k,2) = 0.5*P_branch_3(l,k,2)
           P_branch_3(l,k,3) = 0.5*P_branch_3(l,k,3)
          enddo
         enddo
        endif
c
        RETURN
        END
c
