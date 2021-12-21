        SUBROUTINE CROWN_LAYER(LOG_LAYER)
        save
c
        INCLUDE '../src/parameters.include'
c
        INTEGER I, J, K
C
        REAL Z
        REAL P_crown(4,4,4), P_leaf(4,4,4)
        REAL P_needle(4,4,4),P_branch(4,4,4),P_branch_2(4,4,4)
        REAL P_branch_3(4,4,4),P_branch_4(4,4,4),
     &       P_branch_5(4,4,4),P_branch_6(4,4,4)
        REAL EXP_KAPPA_C_p(4,4), EXP_KAPPA_C_m(4,4)
        real kappa_c_p(4,4), kappa_c_m(4,4)
        real crown_prop_diff_p(3), crown_prop_diff_m(3)
        real crown_p_diff(3,4),dum(4,4)
        real phi_hv_vv,phi_vh_vv,phi_hh_vv
C
        COMPLEX M_crown_P(2,2),M_leaf_P(2,2)
        COMPLEX M_needle_P(2,2),M_branch_P(2,2),M_branch_P2(2,2)
        COMPLEX M_branch_P3(2,2),M_branch_P4(2,2),
     &          M_branch_P5(2,2),M_branch_P6(2,2)
        COMPLEX M_crown_M(2,2),M_leaf_M(2,2)
        COMPLEX M_needle_M(2,2),M_branch_M(2,2),M_branch_M2(2,2)
        COMPLEX M_branch_M3(2,2),M_branch_M4(2,2),
     &          M_branch_M5(2,2),M_branch_M6(2,2)
        COMPLEX LAMBDA_C_p(4),Q_C_p(4,4),Q_C_p_INV(4,4)
        COMPLEX LAMBDA_C_m(4),Q_C_m(4,4),Q_C_m_INV(4,4)
C
        REAL THETA, CTHETA, STHETA
        COMMON /RADAR_ANGLE/ THETA, CTHETA, STHETA
C
        REAL DENSITY, CROWN_HGHT, TRUNK_HGHT
        COMMON /R_CANOPY/ DENSITY, CROWN_HGHT, TRUNK_HGHT
        common /mmats/ M_crown_P, M_crown_M
C
        COMMON /CROWN_PHASE/ P_CROWN,P_NEEDLE,P_LEAF,P_BRANCH,P_BRANCH_2
     &                      ,P_BRANCH_3,P_BRANCH_4,P_BRANCH_5,P_BRANCH_6
        COMMON /CROWN_KAPPA/ kappa_c_p, kappa_c_m
        COMMON /CROWN_EXT/ EXP_KAPPA_C_p, EXP_KAPPA_C_m
        COMMON /NDL_M_MAT/ M_NEEDLE_P, M_NEEDLE_M
        COMMON /LEAF_M_MAT/ M_LEAF_P, M_LEAF_M
        COMMON /BR1_M_MAT/ M_BRANCH_P, M_BRANCH_M
        COMMON /BR2_M_MAT/ M_BRANCH_P2, M_BRANCH_M2

        COMMON /BR3_M_MAT/ M_BRANCH_P3, M_BRANCH_M3
        COMMON /BR4_M_MAT/ M_BRANCH_P4, M_BRANCH_M4
        COMMON /BR5_M_MAT/ M_BRANCH_P5, M_BRANCH_M5
        COMMON /BR6_M_MAT/ M_BRANCH_P6, M_BRANCH_M6
C
        COMMON /CROWN_QS/ Q_C_p,Q_C_p_INV,Q_C_m,Q_C_m_INV
        COMMON /CROWN_EIGEN/ LAMBDA_C_p, LAMBDA_C_m
c
        common /crown_prop_p/ crown_prop_diff_p, crown_prop_diff_m
        common /crown_spec_p/ crown_p_diff
C
        LOGICAL OPEN, STEP_VARIABLE(N_VARIABLES)
        LOGICAL STEP_THIS_TIME(N_VARIABLES)
        LOGICAL  CALL_SUB(N_CALLS,N_SUB_CALLS)
        COMMON /L_FLAGS/ STEP_VARIABLE, STEP_THIS_TIME, CALL_SUB, OPEN
C
        LOGICAL LOG_LAYER
C
        REAL MG_BR1, RHO_BR1, BR1_DENS, BR1_DIAM, BR1_LNG
        REAL MG_BR2, RHO_BR2, BR2_DENS, BR2_DIAM, BR2_LNG
        REAL MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
        REAL MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
        COMMON /R_BR1/ MG_BR1, RHO_BR1, BR1_DENS, BR1_DIAM, BR1_LNG
        COMMON /R_BR2/ MG_BR2, RHO_BR2, BR2_DENS, BR2_DIAM, BR2_LNG
        COMMON /R_LEAF/ MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
        COMMON /R_NDL/ MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG

        REAL MG_BR3, RHO_BR3, BR3_DENS, BR3_DIAM, BR3_LNG
        COMMON /R_BR3/ MG_BR3, RHO_BR3, BR3_DENS, BR3_DIAM, BR3_LNG
        REAL MG_BR4, RHO_BR4, BR4_DENS, BR4_DIAM, BR4_LNG
        COMMON /R_BR4/ MG_BR4, RHO_BR4, BR4_DENS, BR4_DIAM, BR4_LNG
        REAL MG_BR5, RHO_BR5, BR5_DENS, BR5_DIAM, BR5_LNG
        COMMON /R_BR5/ MG_BR5, RHO_BR5, BR5_DENS, BR5_DIAM, BR5_LNG
        REAL MG_BR6, RHO_BR6, BR6_DENS, BR6_DIAM, BR6_LNG
        COMMON /R_BR6/ MG_BR6, RHO_BR6, BR6_DENS, BR6_DIAM, BR6_LNG
c
      IF(LOG_LAYER)THEN

          IF(CALL_SUB(4,2)) CALL BRANCH_PHASE(P_BRANCH)
          IF(CALL_SUB(4,3)) CALL LEAF_PHASE(P_LEAF)
          IF(CALL_SUB(4,4)) CALL NEEDLE_PHASE(P_NEEDLE)
          IF(CALL_SUB(4,5)) CALL BRANCH_PHASE2(P_BRANCH_2)
          IF(CALL_SUB(4,6)) CALL BRANCH_PHASE3(P_BRANCH_3)
          IF(CALL_SUB(4,7)) CALL BRANCH_PHASE4(P_BRANCH_4)
          IF(CALL_SUB(4,8)) CALL BRANCH_PHASE5(P_BRANCH_5)
          IF(CALL_SUB(4,9)) CALL BRANCH_PHASE6(P_BRANCH_6)
c
        DO K=1,4
         DO J=1,4
          DO I=1,4
           P_crown(i,j,k) = LEAF_DENS*P_leaf(i,j,k) 
     &          + NDL_DENS*P_needle(i,j,k) + BR1_DENS*P_branch(i,j,k)
     &          + BR2_DENS*P_branch_2(i,j,k)
     &          + BR3_DENS*P_branch_3(i,j,k)
     &          + BR4_DENS*P_branch_4(i,j,k)
     &          + BR5_DENS*P_branch_5(i,j,k)
     &          + BR6_DENS*P_branch_6(i,j,k)
	  END DO
	 END DO
	END DO
C
        DO J=1,2
         DO I=1,2
          M_crown_P(i,j) = LEAF_DENS*M_leaf_P(i,j) +
     &           NDL_DENS*M_needle_P(i,j)+BR1_DENS*M_branch_P(i,j)
     &           + BR2_DENS*M_branch_P2(i,j)
     &           + BR3_DENS*M_branch_P3(i,j)
     &           + BR4_DENS*M_branch_P4(i,j)
     &           + BR5_DENS*M_branch_P5(i,j)
     &           + BR6_DENS*M_branch_P6(i,j)
          M_crown_M(i,j) = LEAF_DENS*M_leaf_M(i,j) +
     &           NDL_DENS*M_needle_M(i,j)+BR1_DENS*M_branch_M(i,j)
     &           + BR2_DENS*M_branch_M2(i,j)
     &           + BR3_DENS*M_branch_M3(i,j)
     &           + BR4_DENS*M_branch_M4(i,j)
     &           + BR5_DENS*M_branch_M5(i,j)
     &           + BR6_DENS*M_branch_M6(i,j)
	 END DO
	END DO
c
        do i=1,4
          do j=1,4
            dum(j,i) = P_CROWN(j,i,1)
          enddo
        enddo
        call flip_stokes(dum,dum)
        call phase_diff(dum,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        crown_p_diff(1,1) = phi_hh_vv
        crown_p_diff(2,1) = phi_vh_vv
        crown_p_diff(3,1) = phi_hv_vv
c
        do i=1,4
          do j=1,4
            dum(j,i) = P_CROWN(j,i,2)
          enddo
        enddo
        call phase_diff(dum,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        crown_p_diff(1,2) = phi_hh_vv
        crown_p_diff(2,2) = phi_vh_vv
        crown_p_diff(3,2) = phi_hv_vv
c
        do i=1,4
          do j=1,4
            dum(j,i) = P_CROWN(j,i,3)
          enddo
        enddo
        call phase_diff(dum,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        crown_p_diff(1,3) = phi_hh_vv
        crown_p_diff(2,3) = phi_vh_vv
        crown_p_diff(3,3) = phi_hv_vv
c
        do i=1,4
          do j=1,4
            dum(j,i) = P_CROWN(j,i,4)
          enddo
        enddo
        call flip_stokes(dum,dum)
        call phase_diff(dum,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        crown_p_diff(1,4) = phi_hh_vv
        crown_p_diff(2,4) = phi_vh_vv
        crown_p_diff(3,4) = phi_hv_vv
C
      ELSE
C
         do j=1,2
          do i=1,4
           do k=1,4
            P_CROWN(K,I,J) = 0.0
           enddo
          enddo
         enddo
         do j=1,2
          do i=1,2
            M_CROWN_p(i,j) = CMPLX(0.0,0.0)
            M_CROWN_m(i,j) = CMPLX(0.0,0.0)
          enddo
         enddo
         do i=1,4
           do j=1,3
             crown_p_diff(j,i) = 0.0
           enddo
         enddo
      ENDIF
c
        CALL MAT_EIGEN_SUB(M_crown_p,LAMBDA_C_p,Q_C_p,Q_C_p_INV)
        CALL MAT_EIGEN_SUB(M_crown_m,LAMBDA_C_m,Q_C_m,Q_C_m_INV)
c
        CALL KAPPA_MAT_SUB(M_crown_p, kappa_c_p)
        CALL KAPPA_MAT_SUB(M_crown_m, kappa_c_m)
c
        z = -CROWN_HGHT/CTHETA
        CALL EXP_MAT(LAMBDA_C_p,Q_C_p,Q_C_p_INV,z,EXP_KAPPA_C_p)
        CALL EXP_MAT(LAMBDA_C_m,Q_C_m,Q_C_m_INV,z,EXP_KAPPA_C_m)
c
        call phase_diff(EXP_KAPPA_C_p,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        crown_prop_diff_p(1) = phi_hh_vv
        crown_prop_diff_p(2) = phi_vh_vv
        crown_prop_diff_p(3) = phi_hv_vv
c
        call phase_diff(EXP_KAPPA_C_m,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        crown_prop_diff_m(1) = phi_hh_vv
        crown_prop_diff_m(2) = phi_vh_vv
        crown_prop_diff_m(3) = phi_hv_vv
C
        RETURN
        END
c
