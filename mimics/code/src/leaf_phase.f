        SUBROUTINE LEAF_PHASE(P_LEAF)
        save
c
        include 'parameters.include'
c
        INTEGER II,I,L,K,JJ

        REAL  WORK, NORM ,PDF, PDF_LF
        REAL THETAI,PHII,THETAS,PHIS,THETAD,PHID
        REAL P_LEAF(4,4,4), STOKES_MAT(4,4)
        REAL P_LEAF_TMP_1(4,4,4),P_LEAF_TMP_2(4,4,4),P_LEAF_TMP_3(4,4,4)
        REAL EXP_KAPPA_LF_P(4,4),EXP_KAPPA_LF_M(4,4)
C
        REAL MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
        REAL DENSITY, CROWN_HGHT, TRUNK_HGHT
        REAL FREQ_HERTZ, WAVELENGTH, k0
        REAL THETA, CTHETA, STHETA
C
        COMPLEX EPSILONR(N_EPS),EPSILONRC(N_EPS)
C
        COMPLEX J, CWORK
        COMPLEX SMAT_AVG(2,2,2),M_LEAF_P(2,2),M_LEAF_M(2,2)
        COMPLEX SMAT_AVG_TMP_1(2,2,2), SMAT_AVG_TMP_2(2,2,2)
        COMPLEX SMAT_AVG_TMP_3(2,2,2)
        COMPLEX Smat(2,2)
c
        integer n_theta_1, n_theta_2, n_theta_3
        integer n_phi_1, n_phi_2, n_phi_3
        real t_rad_start_1, t_rad_stop_1, delta_t_rad_1
        real t_rad_start_2, t_rad_stop_2, delta_t_rad_2
        real t_rad_start_3, t_rad_stop_3, delta_t_rad_3
        real p_rad_start_1, p_rad_stop_1, delta_p_rad_1
        real p_rad_start_2, p_rad_stop_2, delta_p_rad_2
        real p_rad_start_3, p_rad_stop_3, delta_p_rad_3
        real delta_p1_t1, delta_p2_t2, delta_p3_t3
c
        logical L_PO
C
        COMMON /R_LEAF/ MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
        COMMON /RADAR_FREQ/ FREQ_HERTZ, WAVELENGTH, k0
        COMMON /RADAR_ANGLE/ THETA, CTHETA, STHETA
        COMMON /R_CANOPY/ DENSITY, CROWN_HGHT, TRUNK_HGHT
        COMMON /C_DIELECTRIC/ EPSILONR,EPSILONRC
C
        COMMON /LEAF_SCAT_MAT/ Smat
c
        COMMON /LEAF_EXT/ EXP_KAPPA_LF_p, EXP_KAPPA_LF_m
        COMMON /LEAF_M_MAT/ M_LEAF_P, M_LEAF_M
c
        include 'constants.include'
c
        J = CMPLX(0.0,1.0)
c
        if((WAVELENGTH/(LEAF_DIAM/100.)).lt.1.5)then
            L_PO = .true.
        else
            L_PO = .false.
        endif
        print*,' Leaf Phys Optics = ',L_PO
c
        DO i=1,4
         DO k=1,4
          DO l=1,4
             P_LEAF(l,k,i) = 0.0
             P_LEAF_tmp_1(l,k,i) = 0.0
             P_LEAF_tmp_2(l,k,i) = 0.0
             P_LEAF_tmp_3(l,k,i) = 0.0
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
c
        call pdf_lf_setup(n_theta_1,t_rad_start_1,t_rad_stop_1,
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
        thetad = t_rad_start_1 + delta_t_rad_1*float(jj-1)
         PDF = PDF_LF(THETAD)
         IF(PDF.GT.0.0)THEN
           print*,' thetad, pdf = ',thetad, pdf
c
          do 145 ii=1,n_phi_1
            PHId = p_rad_start_1 + delta_p_rad_1*float(ii-1)
c
            THETAI = THETA
            THETAS = PI-THETA
            PHII   = PI
            PHIS   = 0.0
c
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_1(L,K,1) =
     &               P_LEAF_tmp_1(L,K,1) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_1(L,K,1) =
     &               P_LEAF_tmp_1(L,K,1) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
c
            THETAI = PI - THETA
            THETAS = PI - THETA
            PHII = 0.0
            PHIS   = PI
c
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_1(L,K,2) =
     &             P_LEAF_tmp_1(L,K,2)+PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_1(L,K,2) =
     &             P_LEAF_tmp_1(L,K,2)+PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
c
            THETAI = THETA
            THETAS = THETA
            PHII = PI
            PHIS = 0.0
c
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_1(L,K,3) =
     &             P_LEAF_tmp_1(L,K,3) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_1(L,K,3) =
     &             P_LEAF_tmp_1(L,K,3) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
c
            THETAI = PI - THETA
            THETAS = THETA
            PHII   = 0.0
            PHIS   = PI
c               { PI + PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF

            CALL STOKES_SUB(SMAT,STOKES_MAT)

            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_1(L,K,4) =
     &             P_LEAF_tmp_1(L,K,4) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)

            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_1(L,K,4) =
     &             P_LEAF_tmp_1(L,K,4) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
c
            THETAI = THETA
            THETAS = THETAI
            PHII   = PI
            PHIS   = PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_1(L,K,1) =
     &             Smat_avg_tmp_1(L,K,1) + PDF*Smat(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_1(L,K,1) =
     &             Smat_avg_tmp_1(L,K,1) + PDF*Smat(L,K)
             ENDDO
            ENDDO

c
            THETAI = PI - THETA
            THETAS = PI - THETA
            PHII   = 0.0
            PHIS   = PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_1(L,K,2) =
     &            Smat_avg_tmp_1(L,K,2) + PDF*Smat(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_1(L,K,2) =
     &            Smat_avg_tmp_1(L,K,2) + PDF*Smat(L,K)
             ENDDO
            ENDDO

C
145       CONTINUE
C
         ENDIF
150   CONTINUE
c
      endif
c
      IF((n_theta_2.ne.0).and.(n_phi_2.ne.0))THEN
c
      do 170 jj=1,n_theta_2
        thetad = t_rad_start_2 + delta_t_rad_2*float(jj-1)
         PDF = PDF_LF(THETAD)
         IF(PDF.GT.0.0)THEN
           print*,' thetad, pdf = ',thetad,pdf
c
          do 165 ii=1,n_phi_2
            PHId = p_rad_start_2 + delta_p_rad_2*float(ii-1)
c
            THETAI = THETA
            THETAS = PI-THETA
            PHII   = PI
            PHIS   = 0.0
c
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_2(L,K,1) =
     &               P_LEAF_tmp_2(L,K,1) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_2(L,K,1) =
     &               P_LEAF_tmp_2(L,K,1) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

c
            THETAI = PI - THETA
            THETAS = PI - THETA
            PHII = 0.0
            PHIS   = PI
c
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_2(L,K,2) =
     &             P_LEAF_tmp_2(L,K,2)+PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_2(L,K,2) =
     &             P_LEAF_tmp_2(L,K,2)+PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

c
            THETAI = THETA
            THETAS = THETA
            PHII = PI
            PHIS = 0.0
c
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_2(L,K,3) =
     &             P_LEAF_tmp_2(L,K,3) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_2(L,K,3) =
     &             P_LEAF_tmp_2(L,K,3) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
c
            THETAI = PI - THETA
            THETAS = THETA
            PHII   = 0.0
            PHIS   = PI
c
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_2(L,K,4) =
     &             P_LEAF_tmp_2(L,K,4) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_2(L,K,4) =
     &             P_LEAF_tmp_2(L,K,4) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

c
            THETAI = THETA
            THETAS = THETAI
            PHII   = PI
            PHIS   = PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_2(L,K,1) =
     &             Smat_avg_tmp_2(L,K,1) + PDF*Smat(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_2(L,K,1) =
     &             Smat_avg_tmp_2(L,K,1) + PDF*Smat(L,K)
             ENDDO
            ENDDO

c
            THETAI = PI - THETA
            THETAS = PI - THETA
            PHII   = 0.0
            PHIS   = PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_2(L,K,2) =
     &            Smat_avg_tmp_2(L,K,2) + PDF*Smat(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_2(L,K,2) =
     &            Smat_avg_tmp_2(L,K,2) + PDF*Smat(L,K)
             ENDDO
            ENDDO
C
165       CONTINUE
C
         ENDIF
170   CONTINUE
c
      endif
c
      IF((n_theta_3.ne.0).and.(n_phi_3.ne.0))THEN
c
      do 200 jj=1,n_theta_3
        thetad = t_rad_start_3 + delta_t_rad_3*float(jj-1)
         PDF = PDF_LF(THETAD)
         IF(PDF.GT.0.0)THEN
           print*,' thetad, pdf = ',thetad,pdf
c
          do 195 ii=1,n_phi_3
            PHId = p_rad_start_3 + delta_p_rad_3*float(ii-1)
c
            THETAI = THETA
            THETAS = PI-THETA
            PHII   = PI
            PHIS   = 0.0
c
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_3(L,K,1) =
     &               P_LEAF_tmp_3(L,K,1) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_3(L,K,1) =
     &               P_LEAF_tmp_3(L,K,1) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
c
            THETAI = PI - THETA
            THETAS = PI - THETA
            PHII = 0.0
            PHIS   = PI
c                 {PI + PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_3(L,K,2) =
     &             P_LEAF_tmp_3(L,K,2)+PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_3(L,K,2) =
     &             P_LEAF_tmp_3(L,K,2)+PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
c
            THETAI = THETA
            THETAS = THETA
            PHII = PI
            PHIS = 0.0
c                     {PI + PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_3(L,K,3) =
     &             P_LEAF_tmp_3(L,K,3) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)
C
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_3(L,K,3) =
     &             P_LEAF_tmp_3(L,K,3) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
c
            THETAI = PI - THETA
            THETAS = THETA
            PHII   = 0.0
            PHIS   = PI
c
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            CALL STOKES_SUB(SMAT,STOKES_MAT)
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_3(L,K,4) =
     &             P_LEAF_tmp_3(L,K,4) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            CALL STOKES_SUB(SMAT,STOKES_MAT)
            DO K=1,4
             DO L=1,4
              P_LEAF_tmp_3(L,K,4) =
     &             P_LEAF_tmp_3(L,K,4) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
c
            THETAI = THETA
            THETAS = THETAI
            PHII   = PI
            PHIS   = PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_3(L,K,1) =
     &             Smat_avg_tmp_3(L,K,1) + PDF*Smat(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_3(L,K,1) =
     &             Smat_avg_tmp_3(L,K,1) + PDF*Smat(L,K)
             ENDDO
            ENDDO
c
            THETAI = PI - THETA
            THETAS = PI - THETA
            PHII   = 0.0
            PHIS   = PHII
C
           IF(L_PO) then
            CALL LEAF_PHYS_OPTICS_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetad,phid)
           ELSE
            CALL leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &             thetad,phid)
           ENDIF
C
            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_3(L,K,2) =
     &            Smat_avg_tmp_3(L,K,2) + PDF*Smat(L,K)
             ENDDO
            ENDDO

            smat(1,2) = -smat(1,2)
            smat(2,1) = -smat(2,1)

            DO K=1,2
             DO L=1,2
              Smat_avg_tmp_3(L,K,2) =
     &            Smat_avg_tmp_3(L,K,2) + PDF*Smat(L,K)
             ENDDO
            ENDDO

C
195       CONTINUE
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
             P_LEAF(l,k,i) = WORK*(delta_p1_t1*P_LEAF_tmp_1(l,k,i) +
     &                             delta_p2_t2*P_LEAF_tmp_2(l,k,i) +
     &                             delta_p3_t3*P_LEAF_tmp_3(l,k,i))
	  END DO
	 END DO
	END DO
C
        CWORK = J*2.0*PI/(NORM*k0)
C
        DO k=1,2
         DO l=1,2
             M_LEAF_P(l,k) = CWORK*(delta_p1_t1*Smat_avg_tmp_1(l,k,1) +
     &                              delta_p2_t2*Smat_avg_tmp_2(l,k,1) +
     &                              delta_p3_t3*Smat_avg_tmp_3(l,k,1))
	 END DO
	END DO
        DO k=1,2
         DO l=1,2
             M_LEAF_M(l,k) = CWORK*(delta_p1_t1*Smat_avg_tmp_1(l,k,2) +
     &                              delta_p2_t2*Smat_avg_tmp_2(l,k,2) +
     &                              delta_p3_t3*Smat_avg_tmp_3(l,k,2))
	 END DO
	END DO
C
        RETURN
        END
