        SUBROUTINE SOLVE_CANOPY(LOG_TRUNK,LOG_CROWN)
        save
c
        INCLUDE 'parameters.include'
c
        INTEGER I, J, K
C
        REAL THETA, CTHETA, STHETA, DENSITY, CROWN_HGHT, TRUNK_HGHT
        REAL PHASE_PASS(4,4)
        REAL EXP_KAPPA_T_p(4,4), EXP_KAPPA_T_m(4,4)
        REAL EXP_KAPPA_C_p(4,4), EXP_KAPPA_C_m(4,4)
        REAL TRUNK_PHASE_MAT_m(4,4),TRUNK_PHASE_MAT_p(4,4)
        REAL P_crown(4,4,4),P_leaf(4,4,4),P_needle(4,4,4)
        REAL P_branch(4,4,4)
        REAL P_branch_2(4,4,4), P_branch_3(4,4,4), P_branch_4(4,4,4)
        REAL P_branch_5(4,4,4), P_branch_6(4,4,4)
        REAL GRND_REFLECT_MAT(4,4), GRND_BACK_MAT(4,4)
        REAL Rprime0(4,4), RprimePI(4,4), RDUMMAT1(4,4), RDUMMAT2(4,4)
        REAL RDUMMAT(4,4)
        REAL BACKTERMS(4,4,N_SCAT_TERMS), T(4,4)
        REAL EXP_CANOPY_p(4,4), EXP_CANOPY_m(4,4)
        real phi_hv_vv,phi_vh_vv,phi_hh_vv
        real phase_diff_terms(3,n_scat_terms), phase_diff_total(3)
        real a_phase(3,n_scat_terms)
C
        COMPLEX Q_T_p(4,4),Q_T_p_INV(4,4),Q_T_m(4,4),Q_T_m_INV(4,4)
        COMPLEX Q_C_p(4,4),Q_C_p_INV(4,4),Q_C_m(4,4),Q_C_m_INV(4,4)
        COMPLEX LAMBDA_T_p(4), LAMBDA_T_m(4)
        COMPLEX LAMBDA_C_p(4), LAMBDA_C_m(4)
        COMPLEX CDUMMAT1(4,4), WORK1, WORK2, WORK3
        COMPLEX A_1(4,4),A_2(4,4),A_3(4,4),A_4(4,4),A_5(4,4),A_6(4,4)
C
        LOGICAL LOG_TRUNK,LOG_CROWN
C
        COMMON /RADAR_ANGLE/ THETA, CTHETA, STHETA
        COMMON /R_CANOPY/ DENSITY, CROWN_HGHT, TRUNK_HGHT
        COMMON /R_GROUND_MATS/ GRND_REFLECT_MAT, GRND_BACK_MAT
        COMMON /TRUNK_EXT/ EXP_KAPPA_T_p, EXP_KAPPA_T_m

        REAL P_TRUNK(4,4,2), P_TRUNK_LAYER(4,4,2)
        COMMON /TRUNK_PHASE/ P_TRUNK, P_TRUNK_LAYER

        COMMON /TRUNK_QS/ Q_T_p, Q_T_p_INV, Q_T_m, Q_T_m_INV
        COMMON /TRUNK_EIGEN/ LAMBDA_T_p, LAMBDA_T_m
        COMMON /CROWN_EXT/ EXP_KAPPA_C_p, EXP_KAPPA_C_m
        COMMON /CROWN_PHASE/ P_CROWN,P_NEEDLE,P_LEAF,P_BRANCH,P_BRANCH_2
     &                      ,P_BRANCH_3,P_BRANCH_4,P_BRANCH_5,P_BRANCH_6


        COMMON /CROWN_QS/ Q_C_p, Q_C_p_INV, Q_C_m, Q_C_m_INV
        COMMON /CROWN_EIGEN/ LAMBDA_C_p, LAMBDA_C_m
        COMMON /STOKES_MATS/ T, BACKTERMS
        COMMON /A_MATS/ A_1, A_2, A_3, A_4, A_5, A_6
        COMMON /CANOPY_EXT/ EXP_CANOPY_p, EXP_CANOPY_m
        COMMON /phase_diff_mats/ phase_diff_total, phase_diff_terms
        common /a_phase_diffs/ a_phase
c
        CALL MATMULT2(EXP_KAPPA_C_p,EXP_KAPPA_T_p,EXP_CANOPY_p)
        CALL MATMULT2(EXP_KAPPA_T_m,EXP_KAPPA_C_m,EXP_CANOPY_m)
c
      IF(LOG_CROWN)THEN
C
        DO I=1,4
         DO J=1,4
          PHASE_PASS(I,J) = P_CROWN(I,J,1)
         ENDDO
        ENDDO
C
        CALL CrMATMULT3(Q_C_m_INV,PHASE_PASS,Q_C_p,CDUMMAT1)
C
        DO J=1,4
         DO I=1,4
            WORK1 = (LAMBDA_C_m(I) + LAMBDA_C_p(J))/CTHETA
            WORK2 = (1.0-EXP(-WORK1*CROWN_HGHT))/WORK1
            A_1(I,J) = WORK2*CDUMMAT1(I,J)
	 END DO
	END DO
c
        DO I=1,4
         DO J=1,4
          PHASE_PASS(I,J) = P_CROWN(I,J,2)
         ENDDO
        ENDDO
C
        CALL CrMATMULT3(Q_C_m_INV,PHASE_PASS,Q_C_m,CDUMMAT1)
C
        DO J=1,4
         DO I=1,4
            IF((I.EQ.J).OR.
     &         (CABS((LAMBDA_C_m(I)-LAMBDA_C_m(J))/LAMBDA_C_m(I))
     &            .LT.1.0E-5))THEN
                A_2(I,J) = CEXP(-LAMBDA_C_m(I)*CROWN_HGHT/CTHETA)
     &                     *CDUMMAT1(I,J)*CROWN_HGHT
            ELSE
                WORK1 = (-LAMBDA_C_m(I) + LAMBDA_C_m(J))/CTHETA
                WORK2 = CEXP(-LAMBDA_C_m(I)*CROWN_HGHT/CTHETA) - 
     &                  CEXP(-LAMBDA_C_m(J)*CROWN_HGHT/CTHETA) 
                WORK3 = WORK2/WORK1
                A_2(I,J) = WORK3*CDUMMAT1(I,J)
            ENDIF
	 END DO
	END DO
c
        DO I=1,4
         DO J=1,4
          PHASE_PASS(I,J) = P_CROWN(I,J,3)
         ENDDO
        ENDDO
C
        CALL CrMATMULT3(Q_C_p_INV,PHASE_PASS,Q_C_p,CDUMMAT1)
C
        DO J=1,4
         DO I=1,4
            IF((I.EQ.J).OR.
     &         (CABS((LAMBDA_C_p(I)-LAMBDA_C_p(J))/LAMBDA_C_p(I))
     &            .LT.1.0E-5))THEN
                A_3(I,J) = CEXP(-LAMBDA_C_p(I)*CROWN_HGHT/CTHETA)
     &                     *CDUMMAT1(I,J)*CROWN_HGHT
            ELSE
                WORK1 = (LAMBDA_C_p(I) - LAMBDA_C_p(J))/CTHETA
                WORK2 = CEXP(-LAMBDA_C_p(J)*CROWN_HGHT/CTHETA)- 
     &                  CEXP(-LAMBDA_C_p(I)*CROWN_HGHT/CTHETA) 
                WORK3 = WORK2/WORK1
                A_3(I,J) = WORK3*CDUMMAT1(I,J)
            ENDIF
	 END DO
	END DO
c
        DO I=1,4
         DO J=1,4
          PHASE_PASS(I,J) = P_CROWN(I,J,4)
         ENDDO
        ENDDO
C
        CALL CrMATMULT3(Q_C_p_INV,PHASE_PASS,Q_C_m,CDUMMAT1)
C
        DO J=1,4
         DO I=1,4
            WORK1 = (LAMBDA_C_p(I) + LAMBDA_C_m(J))/CTHETA
            WORK2 = (1.0-CEXP(-WORK1*CROWN_HGHT))/WORK1
            A_4(I,J) = WORK2*CDUMMAT1(I,J)
	 END DO
	END DO
C
      ELSE
C
        DO 110 J=1,4
         DO 109 I=1,4
            A_1(I,J) = CMPLX(0.0,0.0)
            A_2(I,J) = CMPLX(0.0,0.0)
            A_3(I,J) = CMPLX(0.0,0.0)
            A_4(I,J) = CMPLX(0.0,0.0)
109      CONTINUE
110     CONTINUE
C
      ENDIF
C
      IF(LOG_TRUNK)THEN
C
        DO J=1,4
         DO I=1,4
            TRUNK_PHASE_MAT_p(I,J) = P_TRUNK_LAYER(I,J,1)
            TRUNK_PHASE_MAT_m(I,J) = P_TRUNK_LAYER(I,J,2)
         ENDDO
        ENDDO
c
        CALL CrMATMULT3(Q_T_m_INV,TRUNK_PHASE_MAT_m,Q_T_m,CDUMMAT1)
C
        DO J=1,4
         DO I=1,4
            IF((I.EQ.J).OR.
     &         (CABS((LAMBDA_T_m(I)-LAMBDA_T_m(J))/LAMBDA_T_m(I))
     &            .LT.1.0E-5))THEN
                A_5(I,J) = CEXP(-LAMBDA_T_m(I)*TRUNK_HGHT/CTHETA)
     &                     *CDUMMAT1(I,J)*TRUNK_HGHT
            ELSE
                WORK1 = (-LAMBDA_T_m(I) + LAMBDA_T_m(J))/CTHETA
                WORK2 = CEXP(-LAMBDA_T_m(I)*TRUNK_HGHT/CTHETA) - 
     &                  CEXP(-LAMBDA_T_m(J)*TRUNK_HGHT/CTHETA) 
                WORK3 = WORK2/WORK1
                A_5(I,J) = WORK3*CDUMMAT1(I,J)
            ENDIF
	 END DO
	END DO
c
        CALL CrMATMULT3(Q_T_p_INV,TRUNK_PHASE_MAT_p,Q_T_p,CDUMMAT1)
C
        DO J=1,4
         DO I=1,4
            IF((I.EQ.J).OR.
     &         (CABS((LAMBDA_T_p(I)-LAMBDA_T_p(J))/LAMBDA_T_p(I))
     &            .LT.1.0E-5))THEN
                A_6(I,J) = CEXP(-LAMBDA_t_p(I)*TRUNK_HGHT/CTHETA)
     &                     *CDUMMAT1(I,J)*TRUNK_HGHT
            ELSE
                WORK1 = (LAMBDA_t_p(I) - LAMBDA_t_p(J))/CTHETA
                WORK2 = CEXP(-LAMBDA_t_p(J)*TRUNK_HGHT/CTHETA) - 
     &                  CEXP(-LAMBDA_t_p(I)*TRUNK_HGHT/CTHETA) 
                WORK3 = WORK2/WORK1
                A_6(I,J) = WORK3*CDUMMAT1(I,J)
            ENDIF
	 END DO
	END DO
C
      ELSE
C
        DO J=1,4
         DO I=1,4
            A_5(I,J) = CMPLX(0.0,0.0)
            A_6(I,J) = CMPLX(0.0,0.0)
	 END DO
	END DO
      ENDIF
C
c
          CALL MATMULT3(EXP_KAPPA_T_p, GRND_REFLECT_MAT, EXP_KAPPA_T_m,
     &               Rprime0)
          CALL MATMULT3(EXP_KAPPA_T_m, GRND_REFLECT_MAT, EXP_KAPPA_T_p,
     &               RprimePI)
c
        CALL MATMULT3Cr(Q_C_m,A_1,Q_C_p_INV,RDUMMAT1)
c
        call flip_stokes(RDUMMAT1,RDUMMAT)
c
        call phase_diff(RDUMMAT,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        a_phase(1,1) = phi_hh_vv
        a_phase(2,1) = phi_vh_vv
        a_phase(3,1) = phi_hv_vv
c
        CALL MATMULT3(RprimePI,RDUMMAT1,Rprime0,RDUMMAT2)
        CALL MATMULT3(EXP_KAPPA_C_p,RDUMMAT2,EXP_KAPPA_C_m,RDUMMAT1)
c
        call flip_stokes(RDUMMAT1,RDUMMAT)
c
        CALL SET_TERM(RDUMMAT,BACKTERMS,CTHETA,1)
c
        call phase_diff(RDUMMAT,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        phase_diff_terms(1,1) = phi_hh_vv
        phase_diff_terms(2,1) = phi_vh_vv
        phase_diff_terms(3,1) = phi_hv_vv
c
        CALL MATMULT3Cr(Q_C_p,A_3,Q_C_p_INV,RDUMMAT1)
c
        call phase_diff(RDUMMAT1,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        a_phase(1,2) = phi_hh_vv
        a_phase(2,2) = phi_vh_vv
        a_phase(3,2) = phi_hv_vv
c
        CALL MATMULT3(RDUMMAT1,RprimePI,EXP_KAPPA_C_m,RDUMMAT2)
        CALL SET_TERM(RDUMMAT2,BACKTERMS,CTHETA,2)
c
        call phase_diff(RDUMMAT2,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        phase_diff_terms(1,2) = phi_hh_vv
        phase_diff_terms(2,2) = phi_vh_vv
        phase_diff_terms(3,2) = phi_hv_vv
c
        CALL MATMULT3Cr(Q_C_m,A_2,Q_C_m_INV,RDUMMAT1)
c
        call phase_diff(RDUMMAT1,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        a_phase(1,3) = phi_hh_vv
        a_phase(2,3) = phi_vh_vv
        a_phase(3,3) = phi_hv_vv
c
        CALL MATMULT3(EXP_KAPPA_C_p,Rprime0,RDUMMAT1,RDUMMAT2)
        CALL SET_TERM(RDUMMAT2,BACKTERMS,CTHETA,3)
c
        call phase_diff(RDUMMAT2,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        phase_diff_terms(1,3) = phi_hh_vv
        phase_diff_terms(2,3) = phi_vh_vv
        phase_diff_terms(3,3) = phi_hv_vv
c
        CALL MATMULT3Cr(Q_C_p,A_4,Q_C_m_INV,RDUMMAT1)
c
        call flip_stokes(RDUMMAT1,RDUMMAT)
c
        CALL SET_TERM(RDUMMAT,BACKTERMS,CTHETA,4)
c
        call phase_diff(RDUMMAT,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        phase_diff_terms(1,4) = phi_hh_vv
        phase_diff_terms(2,4) = phi_vh_vv
        phase_diff_terms(3,4) = phi_hv_vv
        a_phase(1,4) = phi_hh_vv
        a_phase(2,4) = phi_vh_vv
        a_phase(3,4) = phi_hv_vv
c
        CALL MATMULT3Cr(Q_T_m,A_5,Q_T_m_INV,RDUMMAT1)
c
        call phase_diff(RDUMMAT1,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        a_phase(1,5) = phi_hh_vv
        a_phase(2,5) = phi_vh_vv
        a_phase(3,5) = phi_hv_vv
c
        CALL MATMULT3(GRND_REFLECT_MAT,RDUMMAT1,EXP_KAPPA_C_m,RDUMMAT2)
        CALL MATMULT3(EXP_KAPPA_C_p,EXP_KAPPA_T_p,RDUMMAT2,RDUMMAT1)
        CALL SET_TERM(RDUMMAT1,BACKTERMS,CTHETA,5)
c
        call phase_diff(RDUMMAT1,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        phase_diff_terms(1,5) = phi_hh_vv
        phase_diff_terms(2,5) = phi_vh_vv
        phase_diff_terms(3,5) = phi_hv_vv
c
        CALL MATMULT3Cr(Q_T_p,A_6,Q_T_p_INV,RDUMMAT1)
c
        call phase_diff(RDUMMAT1,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        a_phase(1,6) = phi_hh_vv
        a_phase(2,6) = phi_vh_vv
        a_phase(3,6) = phi_hv_vv
c
        CALL MATMULT3(EXP_KAPPA_C_p,RDUMMAT1,GRND_REFLECT_MAT,RDUMMAT2)
        CALL MATMULT3(RDUMMAT2,EXP_KAPPA_T_m,EXP_KAPPA_C_m,RDUMMAT1)
        CALL SET_TERM(RDUMMAT1,BACKTERMS,CTHETA,6)
c
        call phase_diff(RDUMMAT1,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        phase_diff_terms(1,6) = phi_hh_vv
        phase_diff_terms(2,6) = phi_vh_vv
        phase_diff_terms(3,6) = phi_hv_vv
c
        call phase_diff(GRND_BACK_MAT,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        a_phase(1,7) = phi_hh_vv
        a_phase(2,7) = phi_vh_vv
        a_phase(3,7) = phi_hv_vv
c
       CALL MATMULT3(EXP_KAPPA_T_p,GRND_BACK_MAT,EXP_KAPPA_T_m,RDUMMAT1)
       CALL MATMULT3(EXP_KAPPA_C_p,RDUMMAT1,EXP_KAPPA_C_m,RDUMMAT2)
       CALL SET_TERM(RDUMMAT2,BACKTERMS,ctheta,7)
c
        call phase_diff(RDUMMAT2,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        phase_diff_terms(1,7) = phi_hh_vv
        phase_diff_terms(2,7) = phi_vh_vv
        phase_diff_terms(3,7) = phi_hv_vv
c
        DO I=1,4
         DO J=1,4
            T(J,I) = 0.0
	 END DO
	END DO
C
        DO K=1,7
         DO I=1,4
          DO J=1,4
            T(J,I) = T(J,I) + BACKTERMS(J,I,K)
	  END DO
	 END DO
	END DO
c
        call phase_diff(T,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        phase_diff_total(1) = phi_hh_vv
        phase_diff_total(2) = phi_vh_vv
        phase_diff_total(3) = phi_hv_vv
c
        RETURN
        END
c
        SUBROUTINE SET_TERM(INPUT,BACKTERMS,MU0,K)
        save
c
        INCLUDE 'parameters.include'
c
        INTEGER I,J,K
        REAL MU0, BACKTERMS(4,4,N_SCAT_TERMS),INPUT(4,4)
c
        DO I=1,4
          DO J=1,4
            BACKTERMS(J,I,K) = INPUT(J,I)/MU0
	 END DO
	END DO
C
        RETURN
        END
c
        subroutine phase_diff(stokes,phi_hv_vv,phi_vh_vv,phi_hh_vv)
        save
c
        real stokes(4,4), phi_hv_vv, phi_vh_vv, phi_hh_vv
        real x, y, pi
c
        pi = 3.141592654
c
        if((stokes(1,4).ne.0.0).or.(stokes(1,3).ne.0.0))then
            phi_vh_vv = atan2(stokes(1,4),stokes(1,3))
        else
            phi_vh_vv = 0.0
        endif  
c
        if((stokes(4,1).ne.0.0).or.(stokes(3,1).ne.0.0))then
            phi_hv_vv = atan2(-stokes(4,1),stokes(3,1))
        else
            phi_hv_vv = 0.0
        endif 
c
        x = -2.0*stokes(1,1)*stokes(4,3) + stokes(3,1)*stokes(1,4) +
     &        stokes(1,3)*stokes(4,1)
        y = 2.0*stokes(1,1)*stokes(3,3) - stokes(3,1)*stokes(1,3) +
     &        stokes(1,4)*stokes(4,1)
c
        phi_hh_vv = atan2(x,y)
c
        phi_hh_vv = phi_hh_vv*(180./pi)
        phi_vh_vv = phi_vh_vv*(180./pi)
        phi_hv_vv = phi_hv_vv*(180./pi)
c
        return
        end
c
        subroutine flip_stokes(stokes_in,stokes_out)
        save
c
        real stokes_in(4,4), stokes_out(4,4)

        stokes_out(1,1) =  stokes_in(1,1)
        stokes_out(1,2) =  stokes_in(1,2)
        stokes_out(1,3) =  stokes_in(1,3)
        stokes_out(1,4) =  stokes_in(1,4)
        stokes_out(2,1) =  stokes_in(2,1)
        stokes_out(2,2) =  stokes_in(2,2)
        stokes_out(2,3) =  stokes_in(2,3)
        stokes_out(2,4) =  stokes_in(2,4)
        stokes_out(3,1) = -stokes_in(3,1)
        stokes_out(3,2) = -stokes_in(3,2)
        stokes_out(3,3) = -stokes_in(3,3)
        stokes_out(3,4) = -stokes_in(3,4)
        stokes_out(4,1) = -stokes_in(4,1)
        stokes_out(4,2) = -stokes_in(4,2)
        stokes_out(4,3) = -stokes_in(4,3)
        stokes_out(4,4) = -stokes_in(4,4)
c
        return
        end
