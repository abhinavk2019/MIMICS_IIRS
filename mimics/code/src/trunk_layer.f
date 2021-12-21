        SUBROUTINE TRUNK_LAYER(LOG_LAYER)
        save
c
        include 'parameters.include'
c
        INTEGER I,J,K
C
        REAL Z, dum(4,4)
        REAL DENSITY, CROWN_HGHT, TRUNK_HGHT
        REAL P_TRUNK(4,4,2), P_TRUNK_LAYER(4,4,2)
        REAL EXP_KAPPA_T_p(4,4), EXP_KAPPA_T_m(4,4)
        real kappa_t_p(4,4), kappa_t_m(4,4)
        real trunk_prop_diff_p(3), trunk_prop_diff_m(3)
        real trunk_spec_p_diff(3,2)
        real phi_hv_vv,phi_vh_vv,phi_hh_vv
C
        COMPLEX M_trunk_p(2,2), M_trunk_m(2,2)
        COMPLEX M_TRUNK_LAYER_p(2,2), M_TRUNK_LAYER_m(2,2)
        COMPLEX LAMBDA_T_p(4),Q_T_p(4,4),Q_T_p_INV(4,4)
        COMPLEX LAMBDA_T_m(4),Q_T_m(4,4),Q_T_m_INV(4,4)
C
        COMMON /R_CANOPY/ DENSITY, CROWN_HGHT, TRUNK_HGHT

        COMMON /TRUNK_PHASE/ P_TRUNK, P_TRUNK_LAYER
        COMMON /TRUNK_M_MAT/ M_trunk_p, M_trunk_m
        COMMON /TRUNK_EXT/ EXP_KAPPA_T_p, EXP_KAPPA_T_m
        COMMON /TRUNK_KAPPA/ kappa_t_p, kappa_t_m
        COMMON /TRUNK_QS/ Q_T_p, Q_T_p_INV, Q_T_m, Q_T_m_INV
        COMMON /TRUNK_EIGEN/ LAMBDA_T_p, LAMBDA_T_m
c
        common /trunk_prop_p/ trunk_prop_diff_p, trunk_prop_diff_m
        common /trunk_spec_p/ trunk_spec_p_diff
C
        REAL THETA, CTHETA, STHETA
        COMMON /RADAR_ANGLE/ THETA, CTHETA, STHETA
C
        LOGICAL OPEN, STEP_VARIABLE(N_VARIABLES)
        LOGICAL STEP_THIS_TIME(N_VARIABLES)
        LOGICAL  CALL_SUB(N_CALLS,N_SUB_CALLS)
        COMMON /L_FLAGS/ STEP_VARIABLE, STEP_THIS_TIME, CALL_SUB, OPEN
C
        LOGICAL LOG_LAYER

        common /mmats_tr/ M_TRUNK_LAYER_p,M_TRUNK_LAYER_m
c
        IF(CALL_SUB(3,2)) THEN
            CALL TRUNK_PHASE_SUB(P_TRUNK)
C
        DO J=1,2
         DO I=1,4
          DO K=1,4
            P_TRUNK_LAYER(K,I,J) = DENSITY*P_TRUNK(K,I,J)
	  END DO
	 END DO
	END DO

c
        do i=1,4
          do j=1,4
            dum(j,i) = P_TRUNK_LAYER(j,i,1)
          enddo
        enddo
        call phase_diff(dum,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        trunk_spec_p_diff(1,1) = phi_hh_vv
        trunk_spec_p_diff(2,1) = phi_vh_vv
        trunk_spec_p_diff(3,1) = phi_hv_vv
c
        do i=1,4
          do j=1,4
            dum(j,i) = P_TRUNK_LAYER(j,i,2)
          enddo
        enddo
        call phase_diff(dum,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        trunk_spec_p_diff(1,2) = phi_hh_vv
        trunk_spec_p_diff(2,2) = phi_vh_vv
        trunk_spec_p_diff(3,2) = phi_hv_vv
C
        DO J=1,2
         DO I=1,2
            M_TRUNK_LAYER_p(i,j) = DENSITY*M_trunk_p(i,j)
            M_TRUNK_LAYER_m(i,j) = DENSITY*M_trunk_m(i,j)
	 END DO
	END DO
C
        ELSE IF(.NOT.LOG_LAYER)THEN
         do j=1,2
          do i=1,4
           do k=1,4
            P_TRUNK_LAYER(K,I,J) = 0.0
            P_TRUNK(K,I,J) = 0.0
           enddo
          enddo
         enddo
         do j=1,2
          do i=1,2
            M_TRUNK_LAYER_p(i,j) = CMPLX(0.0,0.0)
            M_TRUNK_LAYER_m(i,j) = CMPLX(0.0,0.0)
          enddo
         enddo
         do j=1,2
          do i=1,3
            trunk_spec_p_diff(i,j) = 0.0
          enddo
         enddo
        ENDIF
c
        CALL MAT_EIGEN_SUB(M_TRUNK_LAYER_p,LAMBDA_T_p,Q_T_p,Q_T_p_INV)
        CALL MAT_EIGEN_SUB(M_TRUNK_LAYER_m,LAMBDA_T_m,Q_T_m,Q_T_m_INV)
c
        CALL KAPPA_MAT_SUB(M_trunk_layer_p, kappa_t_p)
        CALL KAPPA_MAT_SUB(M_trunk_layer_m, kappa_t_m)
c
        z = -TRUNK_HGHT/CTHETA
       CALL EXP_MAT(LAMBDA_T_p,Q_T_p,Q_T_p_INV,z,EXP_KAPPA_T_p)
       CALL EXP_MAT(LAMBDA_T_m,Q_T_m,Q_T_m_INV,z,EXP_KAPPA_T_m)
c
        call phase_diff(EXP_KAPPA_T_p,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        trunk_prop_diff_p(1) = phi_hh_vv
        trunk_prop_diff_p(2) = phi_vh_vv
        trunk_prop_diff_p(3) = phi_hv_vv
c
        call phase_diff(EXP_KAPPA_T_m,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        trunk_prop_diff_m(1) = phi_hh_vv
        trunk_prop_diff_m(2) = phi_vh_vv
        trunk_prop_diff_m(3) = phi_hv_vv
C
        RETURN
        END
c
        SUBROUTINE TRUNK_PHASE_SUB(P_TRUNK)
        save
c
        include 'parameters.include'
c
        INTEGER I, K, L, II, JJ
C
        REAL WORK, PDF, PDF_TR, NORM
        REAL DIAMCM,LENGTHM, K0A
        real thetai,phii,thetas,phis,thetac,phic

        REAL FREQ_HERTZ, WAVELENGTH, k0
        REAL THETA, CTHETA, STHETA
        REAL MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
        REAL DENSITY, CROWN_HGHT, TRUNK_HGHT

c
        COMPLEX EPSILONR(N_EPS),EPSILONRC(N_EPS),EPSR
        COMPLEX J,CWORK
        COMPLEX SMAT(2,2),Smat_avg(2,2,2)
        COMPLEX Smat_avg_tmp_1(2,2,2), Smat_avg_tmp_2(2,2,2)
        COMPLEX Smat_avg_tmp_3(2,2,2), smat_tmp1(2,2), smat_tmp2(2,2)
        COMPLEX M_trunk_p(2,2), M_trunk_m(2,2)
        REAL P_TRUNK(4,4,2), STOKES_MAT(4,4), STOKES_MAT2(4,4)
        REAL P_TRUNK_tmp_1(4,4,2), P_TRUNK_tmp_2(4,4,2)
        REAL P_TRUNK_tmp_3(4,4,2)
c
        integer n_theta_1, n_theta_2, n_theta_3
        integer n_phi_1, n_phi_2, n_phi_3
        real t_rad_start_1, t_rad_stop_1, t_rad_start_2, t_rad_stop_2
        real t_rad_start_3, t_rad_stop_3
        real p_rad_start_1, p_rad_stop_1, p_rad_start_2, p_rad_stop_2
        real p_rad_start_3, p_rad_stop_3
        real delta_p_rad_1, delta_p_rad_2, delta_p_rad_3
        real delta_t_rad_1, delta_t_rad_2, delta_t_rad_3
        real delta_p1_t1, delta_p2_t2, delta_p3_t3
c
        LOGICAL LFLAG
c
        COMMON /RADAR_FREQ/ FREQ_HERTZ, WAVELENGTH, k0
        COMMON /RADAR_ANGLE/ THETA, CTHETA, STHETA
        COMMON /C_DIELECTRIC/ EPSILONR,EPSILONRC
        COMMON /R_TRUNK/ MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
        COMMON /R_CANOPY/ DENSITY, CROWN_HGHT, TRUNK_HGHT
c
        COMMON /TRUNK_M_MAT/ M_trunk_p, M_trunk_m
c
        common /rcyldata/ diamcm,lengthm,k0a
        common /ccyldata/ epsr
        common /cylscatmat/ Smat
c
        LOGICAL LOG_CONSTITUENT(N_CONSTITUENTS)
        LOGICAL LOG_EPS_TABLE(N_CONSTITUENTS)
        LOGICAL LOG_DRY_DENSITY(N_CONSTITUENTS) 
        LOGICAL LOG_PDF_TYPE(N_CONSTITUENTS,N_CONST_VARY)
        LOGICAL LOG_HIST(N_CONSTITUENTS)
C
        COMMON /L_CONFIGURE/ LOG_CONSTITUENT, LOG_EPS_TABLE,
     &                       LOG_DRY_DENSITY, LOG_PDF_TYPE, LOG_HIST
C
        REAL STOKES_MAT_INT(4,4)
        COMPLEX SMAT_INT(2,2)
c
        include 'constants.include'
c
        DO i=1,2
         DO k=1,4
          DO l=1,4
             P_TRUNK(l,k,i) = 0.0
             P_TRUNK_tmp_1(l,k,i) = 0.0
             P_TRUNK_tmp_2(l,k,i) = 0.0
             P_TRUNK_tmp_3(l,k,i) = 0.0
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
        k0a = k0
C
        NORM = 2.0*PI
c
        diamcm = TRUNK_DIAM
        lengthm = TRUNK_HGHT
        epsr = EPSILONRC(5)
c
        phii = 0.0
c
        call pdf_tr_setup(n_theta_1,t_rad_start_1,t_rad_stop_1,
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
        pdf = pdf_tr(thetac)
        IF(PDF.GT.0.0)THEN
          print*,' thetac, pdf = ',thetac,pdf
c
          do 140 ii=1,n_phi_1
            PHIc = p_rad_start_1 + delta_p_rad_1*float(ii-1)
c
            THETAI = THETA
            THETAS = THETA
            PHIS = PI + PHII
C
            IF(LOG_PDF_TYPE(1,1).OR.LOG_HIST(1))THEN
                CALL CYL_SCAT_MAT_SIZE_INT(thetai,phii,thetas,phis,
     &             thetac,phic,STOKES_MAT_INT,SMAT_INT)
                DO K=1,4
                 DO L=1,4
                  STOKES_MAT(L,K) = STOKES_MAT_INT(L,K)
                 ENDDO
                ENDDO
                DO K=1,2
                 DO L=1,2
                  SMAT(L,K) = SMAT_INT(L,K)
                 ENDDO
                ENDDO
            ELSE
                CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
                CALL STOKES_SUB(SMAT,STOKES_MAT)

                smat_tmp1(1,1) =  smat(1,1)
                smat_tmp1(1,2) = -smat(2,1)
                smat_tmp1(2,1) = -smat(1,2)
                smat_tmp1(2,2) =  smat(2,2)

                if(n_phi_1.gt.1)then
                 smat(1,2) = -smat(1,2)
                 smat(2,1) = -smat(2,1)

                 CALL STOKES_SUB(SMAT,STOKES_MAT2)

                 DO K=1,4
                  DO L=1,4
                   STOKES_MAT(L,K)=(STOKES_MAT(L,K)+STOKES_MAT2(L,K))
                  ENDDO
                 ENDDO
                endif

            ENDIF
C
            DO K=1,4
             DO L=1,4
              P_TRUNK_tmp_1(L,K,1) =
     &             P_TRUNK_tmp_1(L,K,1) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO

c
            THETAI = PI - THETA
            THETAS = PI - THETA
            PHIS   = PI + PHII
C
            IF(LOG_PDF_TYPE(1,1).OR.LOG_HIST(1))THEN
                CALL CYL_SCAT_MAT_SIZE_INT(thetai,phii,thetas,phis,
     &             thetac,phic,STOKES_MAT_INT,SMAT_INT)
                DO K=1,4
                 DO L=1,4
                  STOKES_MAT(L,K) = STOKES_MAT_INT(L,K)
                 ENDDO
                ENDDO
                DO K=1,2
                 DO L=1,2
                  SMAT(L,K) = SMAT_INT(L,K)
                 ENDDO
                ENDDO
            ELSE
                CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
                CALL STOKES_SUB(SMAT,STOKES_MAT)

                smat_tmp2(1,1) =  smat(1,1)
                smat_tmp2(1,2) = -smat(2,1)
                smat_tmp2(2,1) = -smat(1,2)
                smat_tmp2(2,2) =  smat(2,2)

                if(n_phi_1.gt.1)then
                 smat(1,2) = -smat(1,2)
                 smat(2,1) = -smat(2,1)

                 CALL STOKES_SUB(SMAT,STOKES_MAT2)

                 DO K=1,4
                  DO L=1,4
                  STOKES_MAT(L,K)=(STOKES_MAT(L,K)+STOKES_MAT2(L,K))
                  ENDDO
                 ENDDO
                endif

            ENDIF
C
            DO K=1,4
             DO L=1,4
              P_TRUNK_tmp_1(L,K,2) =
     &             P_TRUNK_tmp_1(L,K,2) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
c
            IF(LOG_PDF_TYPE(1,1).OR.LOG_HIST(1)) then
             continue
            else 
             CALL STOKES_SUB(smat_tmp1,STOKES_MAT)
             if(n_phi_1.gt.1)then
              smat_tmp1(1,2) = -smat_tmp1(1,2)
              smat_tmp1(2,1) = -smat_tmp1(2,1)
              CALL STOKES_SUB(smat_tmp1,STOKES_MAT2)
              DO K=1,4
               DO L=1,4
                STOKES_MAT(L,K) = STOKES_MAT(L,K) + STOKES_MAT2(L,K) 
               ENDDO
              ENDDO
             endif

             DO K=1,4
              DO L=1,4
               P_TRUNK_tmp_1(L,K,2) =
     &             P_TRUNK_tmp_1(L,K,2) + PDF*STOKES_MAT(L,K)
              ENDDO
             ENDDO

             CALL STOKES_SUB(smat_tmp2,STOKES_MAT)
              if(n_phi_1.gt.1)then
               smat_tmp2(1,2) = -smat_tmp2(1,2)
               smat_tmp2(2,1) = -smat_tmp2(2,1)
               CALL STOKES_SUB(smat_tmp2,STOKES_MAT2)
               DO K=1,4
                DO L=1,4
                 STOKES_MAT(L,K) = STOKES_MAT(L,K) + STOKES_MAT2(L,K) 
                ENDDO
               ENDDO
              endif

              DO K=1,4
               DO L=1,4
                P_TRUNK_tmp_1(L,K,1) =
     &             P_TRUNK_tmp_1(L,K,1) + PDF*STOKES_MAT(L,K)
               ENDDO
              ENDDO
            endif
c
            THETAI = THETA
            THETAS = THETA
            PHIS   = PHII
C
            IF(LOG_PDF_TYPE(1,1).OR.LOG_HIST(1))THEN
                CALL CYL_SCAT_MAT_SIZE_INT(thetai,phii,thetas,phis,
     &             thetac,phic,STOKES_MAT_INT,SMAT_INT)
                DO K=1,4
                 DO L=1,4
                  STOKES_MAT(L,K) = STOKES_MAT_INT(L,K)
                 ENDDO
                ENDDO
                DO K=1,2
                 DO L=1,2
                  SMAT(L,K) = SMAT_INT(L,K)
                 ENDDO
                ENDDO

                DO K=1,2
                 DO L=1,2
                  Smat_avg_tmp_1(L,K,1) =
     &            Smat_avg_tmp_1(L,K,1) + PDF*Smat(L,K)
                 ENDDO
                ENDDO

            ELSE
                CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)

                DO K=1,2
                 DO L=1,2
                  Smat_avg_tmp_1(L,K,1) =
     &            Smat_avg_tmp_1(L,K,1) + PDF*Smat(L,K)
                 ENDDO
                ENDDO

                if(n_phi_1.gt.1)then
                  smat(1,2) = -smat(1,2)
                  smat(2,1) = -smat(2,1)

                  DO K=1,2
                   DO L=1,2
                    Smat_avg_tmp_1(L,K,1) =
     &                Smat_avg_tmp_1(L,K,1) + PDF*Smat(L,K)
                   ENDDO
                  ENDDO
                endif

            ENDIF
c
            THETAI = PI - THETA
            THETAS = PI - THETA
            PHIS   = PHII
C
            IF(LOG_PDF_TYPE(1,1).OR.LOG_HIST(1))THEN
                CALL CYL_SCAT_MAT_SIZE_INT(thetai,phii,thetas,phis,
     &             thetac,phic,STOKES_MAT_INT,SMAT_INT)
                DO K=1,4
                 DO L=1,4
                  STOKES_MAT(L,K) = STOKES_MAT_INT(L,K)
                 ENDDO
                ENDDO
                DO K=1,2
                 DO L=1,2
                  SMAT(L,K) = SMAT_INT(L,K)
                 ENDDO
                ENDDO

                DO K=1,2
                 DO L=1,2
                  Smat_avg_tmp_1(L,K,2) =
     &             Smat_avg_tmp_1(L,K,2) + PDF*Smat(L,K)
                 ENDDO
                ENDDO

            ELSE
                CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)

                DO K=1,2
                 DO L=1,2
                  Smat_avg_tmp_1(L,K,2) =
     &             Smat_avg_tmp_1(L,K,2) + PDF*Smat(L,K)
                 ENDDO
                ENDDO

                if(n_phi_1.gt.1)then
                  smat(1,2) = -smat(1,2)
                  smat(2,1) = -smat(2,1)

                  DO K=1,2
                   DO L=1,2
                    Smat_avg_tmp_1(L,K,2) =
     &             Smat_avg_tmp_1(L,K,2) + PDF*Smat(L,K)
                   ENDDO
                  ENDDO
                endif

            ENDIF
C
C
140        CONTINUE
c
        endif
150    continue
c
       ENDIF
c
      IF((n_theta_2.ne.0).and.(n_phi_2.ne.0))THEN
c
      do 170 jj=1,n_theta_2
        thetac = t_rad_start_2 + delta_t_rad_2*float(jj-1)
        pdf = pdf_tr(thetac)
        IF(PDF.GT.0.0)THEN
          print*,' thetac, pdf = ',thetac,pdf
c
          do 165 ii=1,n_phi_2
            PHIc = p_rad_start_2 + delta_p_rad_2*float(ii-1)
c
            THETAI = THETA
            THETAS = THETA
            PHIS = PI + PHII
C
            IF(LOG_PDF_TYPE(1,1).OR.LOG_HIST(1))THEN
                CALL CYL_SCAT_MAT_SIZE_INT(thetai,phii,thetas,phis,
     &             thetac,phic,STOKES_MAT_INT,SMAT_INT)
                DO K=1,4
                 DO L=1,4
                  STOKES_MAT(L,K) = STOKES_MAT_INT(L,K)
                 ENDDO
                ENDDO
                DO K=1,2
                 DO L=1,2
                  SMAT(L,K) = SMAT_INT(L,K)
                 ENDDO
                ENDDO
            ELSE
                CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
                CALL STOKES_SUB(SMAT,STOKES_MAT)

                smat_tmp1(1,1) =  smat(1,1)
                smat_tmp1(1,2) = -smat(2,1)
                smat_tmp1(2,1) = -smat(1,2)
                smat_tmp1(2,2) =  smat(2,2)

                if(n_phi_2.gt.1)then
                 smat(1,2) = -smat(1,2)
                 smat(2,1) = -smat(2,1)

                 CALL STOKES_SUB(SMAT,STOKES_MAT2)

                 DO K=1,4
                  DO L=1,4
                   STOKES_MAT(L,K)=(STOKES_MAT(L,K)+STOKES_MAT2(L,K))
                  ENDDO
                 ENDDO
                endif

            ENDIF
C
            DO K=1,4
             DO L=1,4
              P_TRUNK_tmp_2(L,K,1) =
     &             P_TRUNK_tmp_2(L,K,1) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
c
            THETAI = PI - THETA
            THETAS = PI - THETA
            PHIS   = PI + PHII
C
            IF(LOG_PDF_TYPE(1,1).OR.LOG_HIST(1))THEN
                CALL CYL_SCAT_MAT_SIZE_INT(thetai,phii,thetas,phis,
     &             thetac,phic,STOKES_MAT_INT,SMAT_INT)
                DO K=1,4
                 DO L=1,4
                  STOKES_MAT(L,K) = STOKES_MAT_INT(L,K)
                 ENDDO
                ENDDO
                DO K=1,2
                 DO L=1,2
                  SMAT(L,K) = SMAT_INT(L,K)
                 ENDDO
                ENDDO
            ELSE
                CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
                CALL STOKES_SUB(SMAT,STOKES_MAT)

                smat_tmp2(1,1) =  smat(1,1)
                smat_tmp2(1,2) = -smat(2,1)
                smat_tmp2(2,1) = -smat(1,2)
                smat_tmp2(2,2) =  smat(2,2)

                if(n_phi_2.gt.1)then
                 smat(1,2) = -smat(1,2)
                 smat(2,1) = -smat(2,1)

                 CALL STOKES_SUB(SMAT,STOKES_MAT2)

                 DO K=1,4
                  DO L=1,4
                   STOKES_MAT(L,K)=(STOKES_MAT(L,K)+STOKES_MAT2(L,K))
                  ENDDO
                 ENDDO
                endif

            ENDIF
C
            DO K=1,4
             DO L=1,4
              P_TRUNK_tmp_2(L,K,2) =
     &             P_TRUNK_tmp_2(L,K,2) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
c
            IF(LOG_PDF_TYPE(1,1).OR.LOG_HIST(1)) then
             continue
            else 
             CALL STOKES_SUB(smat_tmp1,STOKES_MAT)
             if(n_phi_2.gt.1)then
              smat_tmp1(1,2) = -smat_tmp1(1,2)
              smat_tmp1(2,1) = -smat_tmp1(2,1)
              CALL STOKES_SUB(smat_tmp1,STOKES_MAT2)
              DO K=1,4
               DO L=1,4
                STOKES_MAT(L,K) = STOKES_MAT(L,K) + STOKES_MAT2(L,K) 
               ENDDO
              ENDDO
             endif

             DO K=1,4
              DO L=1,4
               P_TRUNK_tmp_2(L,K,2) =
     &             P_TRUNK_tmp_2(L,K,2) + PDF*STOKES_MAT(L,K)
              ENDDO
             ENDDO

             CALL STOKES_SUB(smat_tmp2,STOKES_MAT)
             if(n_phi_2.gt.1)then
              smat_tmp2(1,2) = -smat_tmp2(1,2)
              smat_tmp2(2,1) = -smat_tmp2(2,1)
              CALL STOKES_SUB(smat_tmp2,STOKES_MAT2)
              DO K=1,4
               DO L=1,4
                STOKES_MAT(L,K) = STOKES_MAT(L,K) + STOKES_MAT2(L,K) 
               ENDDO
              ENDDO
             endif

             DO K=1,4
              DO L=1,4
               P_TRUNK_tmp_2(L,K,1) =
     &             P_TRUNK_tmp_2(L,K,1) + PDF*STOKES_MAT(L,K)
              ENDDO
             ENDDO
            endif

c
            THETAI = THETA
            THETAS = THETA
            PHIS   = PHII
C
            IF(LOG_PDF_TYPE(1,1).OR.LOG_HIST(1))THEN
                CALL CYL_SCAT_MAT_SIZE_INT(thetai,phii,thetas,phis,
     &             thetac,phic,STOKES_MAT_INT,SMAT_INT)
                DO K=1,4
                 DO L=1,4
                  STOKES_MAT(L,K) = STOKES_MAT_INT(L,K)
                 ENDDO
                ENDDO
                DO K=1,2
                 DO L=1,2
                  SMAT(L,K) = SMAT_INT(L,K)
                 ENDDO
                ENDDO

                DO K=1,2
                 DO L=1,2
                  Smat_avg_tmp_2(L,K,1) =
     &              Smat_avg_tmp_2(L,K,1) + PDF*Smat(L,K)
                 ENDDO
                ENDDO

            ELSE
                CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)

                DO K=1,2
                 DO L=1,2
                  Smat_avg_tmp_2(L,K,1) =
     &              Smat_avg_tmp_2(L,K,1) + PDF*Smat(L,K)
                 ENDDO
                ENDDO 

                if(n_phi_2.gt.1)then
                  smat(1,2) = -smat(1,2)
                  smat(2,1) = -smat(2,1)

                  DO K=1,2
                   DO L=1,2
                    Smat_avg_tmp_2(L,K,1) =
     &             Smat_avg_tmp_2(L,K,1) + PDF*Smat(L,K)
                   ENDDO
                  ENDDO
                endif


            ENDIF
c
            THETAI = PI - THETA
            THETAS = PI - THETA
            PHIS   = PHII
C
            IF(LOG_PDF_TYPE(1,1).OR.LOG_HIST(1))THEN
                CALL CYL_SCAT_MAT_SIZE_INT(thetai,phii,thetas,phis,
     &             thetac,phic,STOKES_MAT_INT,SMAT_INT)
                DO K=1,4
                 DO L=1,4
                  STOKES_MAT(L,K) = STOKES_MAT_INT(L,K)
                 ENDDO
                ENDDO
                DO K=1,2
                 DO L=1,2
                  SMAT(L,K) = SMAT_INT(L,K)
                 ENDDO
                ENDDO 

                DO K=1,2
                 DO L=1,2
                  Smat_avg_tmp_2(L,K,2) =
     &             Smat_avg_tmp_2(L,K,2) + PDF*Smat(L,K)
                 ENDDO
                ENDDO
            ELSE
                CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)

                DO K=1,2
                 DO L=1,2
                  Smat_avg_tmp_2(L,K,2) =
     &             Smat_avg_tmp_2(L,K,2) + PDF*Smat(L,K)
                 ENDDO
                ENDDO

                if(n_phi_2.gt.1)then
                  smat(1,2) = -smat(1,2)
                  smat(2,1) = -smat(2,1)

                  DO K=1,2
                   DO L=1,2
                    Smat_avg_tmp_2(L,K,2) =
     &             Smat_avg_tmp_2(L,K,2) + PDF*Smat(L,K)
                   ENDDO
                  ENDDO
                endif
            ENDIF
C
165        CONTINUE
c
        endif
170    continue
c
       ENDIF
c
      IF((n_theta_3.ne.0).and.(n_phi_3.ne.0))THEN
c
      do 200 jj=1,n_theta_3
        thetac = t_rad_start_3 + delta_t_rad_3*float(jj-1)
        pdf = pdf_tr(thetac)
        IF(PDF.GT.0.0)THEN
          print*,' thetac, pdf = ',thetac,pdf
c
          do 195 ii=1,n_phi_3
            PHIc = p_rad_start_3 + delta_p_rad_3*float(ii-1)
c
            THETAI = THETA
            THETAS = THETA
            PHIS = PI + PHII
C
            IF(LOG_PDF_TYPE(1,1).OR.LOG_HIST(1))THEN
                CALL CYL_SCAT_MAT_SIZE_INT(thetai,phii,thetas,phis,
     &             thetac,phic,STOKES_MAT_INT,SMAT_INT)
                DO K=1,4
                 DO L=1,4
                  STOKES_MAT(L,K) = STOKES_MAT_INT(L,K)
                 ENDDO
                ENDDO
                DO K=1,2
                 DO L=1,2
                  SMAT(L,K) = SMAT_INT(L,K)
                 ENDDO
                ENDDO
            ELSE
                CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
                CALL STOKES_SUB(SMAT,STOKES_MAT)

                smat_tmp1(1,1) =  smat(1,1)
                smat_tmp1(1,2) = -smat(2,1)
                smat_tmp1(2,1) = -smat(1,2)
                smat_tmp1(2,2) =  smat(2,2)

                if(n_phi_3.gt.1)then
                 smat(1,2) = -smat(1,2)
                 smat(2,1) = -smat(2,1)

                 CALL STOKES_SUB(SMAT,STOKES_MAT2)

                 DO K=1,4
                  DO L=1,4
                   STOKES_MAT(L,K)=(STOKES_MAT(L,K)+STOKES_MAT2(L,K))
                  ENDDO
                 ENDDO
                endif

            ENDIF
C
            DO K=1,4
             DO L=1,4
              P_TRUNK_tmp_2(L,K,1) =
     &             P_TRUNK_tmp_3(L,K,1) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
c
            THETAI = PI - THETA
            THETAS = PI - THETA
            PHIS   = PI + PHII
C
            IF(LOG_PDF_TYPE(1,1).OR.LOG_HIST(1))THEN
                CALL CYL_SCAT_MAT_SIZE_INT(thetai,phii,thetas,phis,
     &             thetac,phic,STOKES_MAT_INT,SMAT_INT)
                DO K=1,4
                 DO L=1,4
                  STOKES_MAT(L,K) = STOKES_MAT_INT(L,K)
                 ENDDO
                ENDDO
                DO K=1,2
                 DO L=1,2
                  SMAT(L,K) = SMAT_INT(L,K)
                 ENDDO
                ENDDO
            ELSE
                CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
                CALL STOKES_SUB(SMAT,STOKES_MAT)

                smat_tmp2(1,1) =  smat(1,1)
                smat_tmp2(1,2) = -smat(2,1)
                smat_tmp2(2,1) = -smat(1,2)
                smat_tmp2(2,2) =  smat(2,2)

                if(n_phi_3.gt.1)then
                 smat(1,2) = -smat(1,2)
                 smat(2,1) = -smat(2,1)

                 CALL STOKES_SUB(SMAT,STOKES_MAT2)

                 DO K=1,4
                  DO L=1,4
                   STOKES_MAT(L,K)=(STOKES_MAT(L,K)+STOKES_MAT2(L,K))
                  ENDDO
                 ENDDO
                endif

            ENDIF
C
            DO K=1,4
             DO L=1,4
              P_TRUNK_tmp_3(L,K,2) =
     &             P_TRUNK_tmp_3(L,K,2) + PDF*STOKES_MAT(L,K)
             ENDDO
            ENDDO
c
            IF(LOG_PDF_TYPE(1,1).OR.LOG_HIST(1)) then
             continue
            else 
             CALL STOKES_SUB(smat_tmp1,STOKES_MAT) 
             if(n_phi_3.gt.1)then
              smat_tmp1(1,2) = -smat_tmp1(1,2)
              smat_tmp1(2,1) = -smat_tmp1(2,1)
              CALL STOKES_SUB(smat_tmp1,STOKES_MAT2)
              DO K=1,4
               DO L=1,4
                STOKES_MAT(L,K) = STOKES_MAT(L,K) + STOKES_MAT2(L,K) 
               ENDDO
              ENDDO
             endif

             DO K=1,4
              DO L=1,4
               P_TRUNK_tmp_3(L,K,2) =
     &             P_TRUNK_tmp_3(L,K,2) + PDF*STOKES_MAT(L,K)
              ENDDO
             ENDDO

             CALL STOKES_SUB(smat_tmp2,STOKES_MAT)
             if(n_phi_3.gt.1)then
              smat_tmp2(1,2) = -smat_tmp2(1,2)
              smat_tmp2(2,1) = -smat_tmp2(2,1)
              CALL STOKES_SUB(smat_tmp2,STOKES_MAT2)
              DO K=1,4
               DO L=1,4
                STOKES_MAT(L,K) = STOKES_MAT(L,K) + STOKES_MAT2(L,K) 
               ENDDO
              ENDDO
             endif

             DO K=1,4
              DO L=1,4
               P_TRUNK_tmp_3(L,K,1) =
     &             P_TRUNK_tmp_3(L,K,1) + PDF*STOKES_MAT(L,K)
              ENDDO
             ENDDO
            endif

c
            THETAI = THETA
            THETAS = THETA
            PHIS   = PHII
C
            IF(LOG_PDF_TYPE(1,1).OR.LOG_HIST(1))THEN
                CALL CYL_SCAT_MAT_SIZE_INT(thetai,phii,thetas,phis,
     &             thetac,phic,STOKES_MAT_INT,SMAT_INT)
                DO K=1,4
                 DO L=1,4
                  STOKES_MAT(L,K) = STOKES_MAT_INT(L,K)
                 ENDDO
                ENDDO
                DO K=1,2
                 DO L=1,2
                  SMAT(L,K) = SMAT_INT(L,K)
                 ENDDO
                ENDDO

                DO K=1,2
                 DO L=1,2
                  Smat_avg_tmp_3(L,K,1) =
     &              Smat_avg_tmp_3(L,K,1) + PDF*Smat(L,K)
                 ENDDO
                ENDDO

            ELSE
                CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)

                DO K=1,2
                 DO L=1,2
                  Smat_avg_tmp_3(L,K,1) =
     &              Smat_avg_tmp_3(L,K,1) + PDF*Smat(L,K)
                 ENDDO
                ENDDO

                if(n_phi_3.gt.1)then
                  smat(1,2) = -smat(1,2)
                  smat(2,1) = -smat(2,1)

                  DO K=1,2
                   DO L=1,2
                    Smat_avg_tmp_3(L,K,2) =
     &             Smat_avg_tmp_3(L,K,2) + PDF*Smat(L,K)
                   ENDDO
                  ENDDO
                endif

            ENDIF
c
            THETAI = PI - THETA
            THETAS = PI - THETA
            PHIS   = PHII
C
            IF(LOG_PDF_TYPE(1,1).OR.LOG_HIST(1))THEN
                CALL CYL_SCAT_MAT_SIZE_INT(thetai,phii,thetas,phis,
     &             thetac,phic,STOKES_MAT_INT,SMAT_INT)
                DO K=1,4
                 DO L=1,4
                  STOKES_MAT(L,K) = STOKES_MAT_INT(L,K)
                 ENDDO
                ENDDO
                DO K=1,2
                 DO L=1,2
                  SMAT(L,K) = SMAT_INT(L,K)
                 ENDDO
                ENDDO

                DO K=1,2
                 DO L=1,2
                  Smat_avg_tmp_3(L,K,2) =
     &              Smat_avg_tmp_3(L,K,2) + PDF*Smat(L,K)
                 ENDDO
                ENDDO

            ELSE
                CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)

                DO K=1,2
                 DO L=1,2
                  Smat_avg_tmp_3(L,K,2) =
     &              Smat_avg_tmp_3(L,K,2) + PDF*Smat(L,K)
                 ENDDO
                ENDDO

                if(n_phi_3.gt.1)then
                  smat(1,2) = -smat(1,2)
                  smat(2,1) = -smat(2,1)

                  DO K=1,2
                   DO L=1,2
                    Smat_avg_tmp_3(L,K,2) =
     &             Smat_avg_tmp_3(L,K,2) + PDF*Smat(L,K)
                   ENDDO
                  ENDDO
                endif

            ENDIF
C
195        CONTINUE
c
        endif
200    continue
c
       ENDIF
c
        WORK = 1.0/(NORM*TRUNK_HGHT)
C
        DO i=1,2
         DO k=1,4
          DO l=1,4
           IF(LOG_PDF_TYPE(1,1).OR.LOG_HIST(1)) then
            P_TRUNK(l,k,i) = WORK*(delta_p1_t1*P_TRUNK_tmp_1(l,k,i) +
     &            delta_p2_t2*P_TRUNK_tmp_2(l,k,i) +
     &            delta_p3_t3*P_TRUNK_tmp_3(l,k,i))
           else 
            P_TRUNK(l,k,i) = 0.5*WORK*(delta_p1_t1*P_TRUNK_tmp_1(l,k,i)+
     &            delta_p2_t2*P_TRUNK_tmp_2(l,k,i) +
     &            delta_p3_t3*P_TRUNK_tmp_3(l,k,i))
           endif
	  END DO
	 END DO
	END DO
C
        J = CMPLX(0.0,1.0)
        CWORK = J*2.0*PI/(NORM*TRUNK_HGHT*k0)
C
        DO k=1,2
         DO l=1,2
           M_TRUNK_P(l,k) = CWORK*(delta_p1_t1*Smat_avg_tmp_1(l,k,1) +
     &            delta_p2_t2*Smat_avg_tmp_2(l,k,1) +
     &            delta_p3_t3*Smat_avg_tmp_3(l,k,1))
	 END DO
	END DO
C
        DO k=1,2
         DO l=1,2
           M_TRUNK_M(l,k) = CWORK*(delta_p1_t1*Smat_avg_tmp_1(l,k,2) +
     &            delta_p2_t2*Smat_avg_tmp_2(l,k,2) +
     &            delta_p3_t3*Smat_avg_tmp_3(l,k,2))
	 END DO
	END DO
c
        RETURN
        END
