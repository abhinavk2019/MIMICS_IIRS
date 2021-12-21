        SUBROUTINE WRITE_DATA(WRITE_HEADER)
        save
c
        INCLUDE 'parameters.include'
c
        integer numfamily, k ,I, j, INT_VAR, ipoint, icount, i_surf
        LOGICAL WRITE_HEADER
        CHARACTER*12 CHAR_FORM
        CHARACTER*11 CHAR_OUT
        character*9 char, char2
        CHARACTER*5 CHAR_VAR
C
        INTEGER PARAM_NUM(N_VARIABLES), INEST(N_VARIABLES)
        INTEGER LOOP_NUM(N_VARIABLES), LOOP_COUNT(N_VARIABLES)
C
        COMMON /I_NEST/ PARAM_NUM, INEST
        COMMON /I_COUNT/ LOOP_NUM, LOOP_COUNT
C
        INTEGER PARAM_OUT(N_VARIABLES), NEST_OUT(N_VARIABLES)
        INTEGER LOOP_OUT(N_VARIABLES)
        logical nest_set(N_VARIABLES),STEP_SET(N_VARIABLES)
C
        LOGICAL LOG_CONSTITUENT(N_CONSTITUENTS)
        LOGICAL LOG_EPS_TABLE(N_CONSTITUENTS)
        LOGICAL LOG_DRY_DENSITY(N_CONSTITUENTS) 
        LOGICAL LOG_PDF_TYPE(N_CONSTITUENTS,N_CONST_VARY)
        LOGICAL LOG_HIST(N_CONSTITUENTS)
        LOGICAL SOIL_SURFACE,WATER_SURFACE,SNOW_SURFACE,ICE_SURFACE
C
        COMMON /L_CONFIGURE/ LOG_CONSTITUENT, LOG_EPS_TABLE,
     &                       LOG_DRY_DENSITY, LOG_PDF_TYPE, LOG_HIST
        COMMON /L_SURFACE_TYPE/ SOIL_SURFACE,WATER_SURFACE,
     &                          SNOW_SURFACE,ICE_SURFACE
        common /surface_mod/ i_surf

        LOGICAL L_POL, L_PHASE, L_KAPPA, L_TRANS
        COMMON /L_FILES/ L_POL, L_PHASE, L_KAPPA, L_TRANS


        complex M_crown_P(2,2), M_crown_M(2,2)
        complex M_TRUNK_LAYER_p(2,2),M_TRUNK_LAYER_m(2,2)
        common /mmats_tr/ M_TRUNK_LAYER_p,M_TRUNK_LAYER_m
        common /mmats/ M_crown_P, M_crown_M

        REAL P_TRUNK(4,4,2), P_TRUNK_LAYER(4,4,2)
        COMMON /TRUNK_PHASE/ P_TRUNK, P_TRUNK_LAYER

        real kappa_t_p(4,4), kappa_t_m(4,4)
        COMMON /TRUNK_KAPPA/ kappa_t_p, kappa_t_m

c
        REAL THETA_DEGREES, FREQ_GHZ
        REAL MV_SOIL, RMS_SOIL, LS_SOIL, SAND, CLAY, SALT
        REAL MG_TRUNK,DENSITY,RHO_TRUNK,CROWN_HGHT,TRUNK_DIAM,TRUNK_HGHT
        REAL MG_LEAF, RHO_LEAF, LEAF_DENS, LEAF_DIAM, LEAF_TAU
        REAL MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
        REAL MG_BR1, RHO_BR1, BR1_DENS, BR1_DIAM, BR1_LNG
        REAL MG_BR2, RHO_BR2, BR2_DENS, BR2_DIAM, BR2_LNG
        REAL MG_BR3, RHO_BR3, BR3_DENS, BR3_DIAM, BR3_LNG
        REAL MG_BR4, RHO_BR4, BR4_DENS, BR4_DIAM, BR4_LNG
        REAL MG_BR5, RHO_BR5, BR5_DENS, BR5_DIAM, BR5_LNG
        REAL MG_BR6, RHO_BR6, BR6_DENS, BR6_DIAM, BR6_LNG
C
        REAL T_SOIL, T_VEG, T_WATER
        REAL T_SNOW
C
        COMMON /R_SENSOR/ THETA_DEGREES, FREQ_GHZ
        COMMON /R_SOIL/ MV_SOIL, SAND, CLAY
        COMMON /R_WATER/ SALT
        COMMON /R_SURFACE/ RMS_SOIL, LS_SOIL
        COMMON /R_TRUNK/ MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
        COMMON /R_CANOPY/ DENSITY, CROWN_HGHT, TRUNK_HGHT
        COMMON /R_LEAF/ MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
        COMMON /R_NDL/ MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
        COMMON /R_BR1/ MG_BR1, RHO_BR1, BR1_DENS, BR1_DIAM, BR1_LNG
        COMMON /R_BR2/ MG_BR2, RHO_BR2, BR2_DENS, BR2_DIAM, BR2_LNG
        COMMON /R_BR3/ MG_BR3, RHO_BR3, BR3_DENS, BR3_DIAM, BR3_LNG
        COMMON /R_BR4/ MG_BR4, RHO_BR4, BR4_DENS, BR4_DIAM, BR4_LNG
        COMMON /R_BR5/ MG_BR5, RHO_BR5, BR5_DENS, BR5_DIAM, BR5_LNG
        COMMON /R_BR6/ MG_BR6, RHO_BR6, BR6_DENS, BR6_DIAM, BR6_LNG
C
        COMMON /R_ENVIRON/ T_SOIL, T_WATER, T_VEG
        COMMON /R_SNOW/ T_SNOW
C
        REAL VARIABLE(N_VARIABLES)
C
        COMPLEX EPSILONR(N_EPS),EPSILONRC(N_EPS)
        COMMON /C_DIELECTRIC/ EPSILONR,EPSILONRC
c
        REAL GRND_REFLECT_MAT(4,4), GRND_BACK_MAT(4,4)
        REAL EXP_KAPPA_T_p(4,4), EXP_KAPPA_T_m(4,4)
        REAL EXP_KAPPA_C_p(4,4), EXP_KAPPA_C_m(4,4)
        REAL EXP_CANOPY_p(4,4), EXP_CANOPY_m(4,4)
        REAL T(4,4), BACKTERMS(4,4,N_SCAT_TERMS)
        REAL T_SIG0(4,4), SIG0(4,4,N_SCAT_TERMS)
        REAL T_SIG0_dB(4,4), SIG0_dB(4,4,N_SCAT_TERMS)
        real phase_diff_terms(3,n_scat_terms), phase_diff_total(3)
        real gnd_spec_p_diff(3), gnd_back_p_diff(3)
        real trunk_prop_diff_p(3), trunk_prop_diff_m(3)
        real trunk_spec_p_diff(3,2)
        real crown_prop_diff_p(3), crown_prop_diff_m(3)
        real crown_p_diff(3,4)
        real a_phase(3,n_scat_terms)
        REAL P_CROWN(4,4,4),P_NEEDLE(4,4,4),P_LEAF(4,4,4)
        REAL P_BRANCH(4,4,4)
        REAL P_BRANCH_2(4,4,4),P_BRANCH_3(4,4,4),P_BRANCH_4(4,4,4),
     &          P_BRANCH_5(4,4,4),P_BRANCH_6(4,4,4)

c
        COMPLEX Q_T_p(4,4), Q_T_p_INV(4,4), Q_T_m(4,4), Q_T_m_INV(4,4)
        COMPLEX Q_C_p(4,4), Q_C_p_INV(4,4), Q_C_m(4,4), Q_C_m_INV(4,4)
        COMPLEX A_1(4,4),A_2(4,4),A_3(4,4),A_4(4,4),A_5(4,4),A_6(4,4)
c
        COMMON /R_GROUND_MATS/ GRND_REFLECT_MAT, GRND_BACK_MAT
        COMMON /TRUNK_EXT/ EXP_KAPPA_T_p, EXP_KAPPA_T_m



        COMMON /TRUNK_QS/ Q_T_p, Q_T_p_INV, Q_T_m, Q_T_m_INV
        real kappa_c_p(4,4), kappa_c_m(4,4)
        COMMON /CROWN_KAPPA/ kappa_c_p, kappa_c_m

        COMMON /CROWN_EXT/ EXP_KAPPA_C_p, EXP_KAPPA_C_m
        COMMON /CROWN_PHASE/ P_CROWN,P_NEEDLE,P_LEAF,P_BRANCH,P_BRANCH_2
     &                      ,P_BRANCH_3,P_BRANCH_4,P_BRANCH_5,P_BRANCH_6
        COMMON /CROWN_QS/ Q_C_p, Q_C_p_INV, Q_C_m, Q_C_m_INV
        COMMON /STOKES_MATS/ T, BACKTERMS
        COMMON /A_MATS/ A_1, A_2, A_3, A_4, A_5, A_6
        COMMON /CANOPY_EXT/ EXP_CANOPY_p, EXP_CANOPY_m
        COMMON /SIGMA_0/ T_SIG0,T_SIG0_dB, SIG0, SIG0_dB
c
        COMMON /phase_diff_mats/ phase_diff_total, phase_diff_terms
        common /grnd_p_diff/ gnd_spec_p_diff, gnd_back_p_diff
        common /trunk_prop_p/ trunk_prop_diff_p, trunk_prop_diff_m
        common /trunk_spec_p/ trunk_spec_p_diff
        common /crown_prop_p/ crown_prop_diff_p, crown_prop_diff_m
        common /crown_spec_p/ crown_p_diff
        common /a_phase_diffs/ a_phase
c
        LOGICAL OPEN, STEP_VARIABLE(N_VARIABLES)
        LOGICAL STEP_THIS_TIME(N_VARIABLES)
        LOGICAL  CALL_SUB(N_CALLS,N_SUB_CALLS)
        COMMON /L_FLAGS/ STEP_VARIABLE, STEP_THIS_TIME, CALL_SUB, OPEN
C
        LOGICAL CHANGE_VAR(N_VARIABLES),STEP_LAST_TIME(N_VARIABLES)
        COMMON /L_FLAG_SWAP/ CHANGE_VAR,STEP_LAST_TIME
C
        INTEGER I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &          I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &          I_PDF_NDL,I_PDF_TRUNK

        COMMON /I_PDF/ I_PDF_LEAF,I_PDF_BR_1,I_PDF_BR_2,
     &                 I_PDF_BR_3, I_PDF_BR_4, I_PDF_BR_5, I_PDF_BR_6,
     &                 I_PDF_NDL,I_PDF_TRUNK


        INTEGER I_PDF_LEAF_SIZE, I_PDF_BR_1_SIZE, I_PDF_BR_2_SIZE,
     &          I_PDF_BR_3_SIZE, I_PDF_BR_4_SIZE, I_PDF_BR_5_SIZE, 
     &          I_PDF_BR_6_SIZE, I_PDF_NDL_SIZE, I_PDF_TRUNK_SIZE

        COMMON /I_PDF_SIZE/ I_PDF_LEAF_SIZE,I_PDF_BR_1_SIZE,
     &               I_PDF_BR_2_SIZE, I_PDF_BR_3_SIZE, I_PDF_BR_4_SIZE, 
     &               I_PDF_BR_5_SIZE, I_PDF_BR_6_SIZE, I_PDF_NDL_SIZE,
     &               I_PDF_TRUNK_SIZE


C
        CHARACTER*30 C_PDF_NDL,C_PDF_BR1,C_PDF_BR2,C_PDF_BR3,C_PDF_BR4,
     &               C_PDF_BR5,C_PDF_BR6,C_PDF_LF,C_PDF_TR
        CHARACTER*30 C_PDF_NDL_SZ,C_PDF_BR1_SZ,C_PDF_BR2_SZ,
     &               C_PDF_BR3_SZ,C_PDF_BR4_SZ,C_PDF_BR5_SZ,
     &               C_PDF_BR6_SZ,C_PDF_LF_SZ,C_PDF_TR_SZ
        character*30 c_surf
C
        INTEGER ILOW
        LOGICAL CHECK_NEST(N_VARIABLES)

        INTEGER TRUNK, BR1, BR2, BR3,BR4,BR5,BR6, LEAVES, NEEDLE, GROUND
c
        INCLUDE 'constants.include'
c
        TRUNK = 1
        BR1   = 2
        BR2   = 3
        BR3   = 4
        BR4   = 5
        BR5   = 6
        BR6   = 7
        LEAVES= 8
        NEEDLE= 9
        GROUND=10
c
        if (i_surf.eq.1) then
            c_surf = 'Geometrical Optics'
        else if (i_surf.eq.2) then
            c_surf = 'Physical Optics'
        else if (i_surf.eq.3) then
            c_surf = 'Small Perturbation'
        else if (i_surf.eq.0) then
            c_surf = 'Specular'
        else if (i_surf.eq.4) then
            c_surf = 'UMich Empirical'
        else
            print*,' ERROR -- Invalid surface model (in write_data.f)'
            stop
        endif
C
        IF(.NOT.OPEN)THEN
c
         open(unit=8,file='../results/forest_sigma_like.out')
         open(unit=7,file='../results/forest_sigma_cross.out')
         open(unit=11,file='../results/forest_trans.out')
         open(unit=12,file='../results/forest_phase_terms.out')
         open(unit=13,file='../results/forest_phase_crown.out')
         open(unit=14,file='../results/forest_phase_trunk.out')

         if(l_pol)open(unit=15,file='../results/polarimetric.out')
         if(l_phase)
     &       open(unit=16,file='../results/forest_phase_mats.out')
         if(l_kappa)
     &       open(unit=17,file='../results/forest_kappa_mats.out')
         if(l_trans)
     &       open(unit=18,file='../results/forest_trans_mats.out')
         open(unit=19,file='../results/mimics2.data')


c
            numfamily = 0
c
c
            CHAR_VAR = ' 1. 5a'
            write(8,5) CHAR_VAR
            write(7,5) CHAR_VAR
            write(11,5) CHAR_VAR
            write(12,5) CHAR_VAR
            write(13,5) CHAR_VAR
            write(14,5) CHAR_VAR
            if(l_pol) write(15,5) CHAR_VAR
            if(l_phase)write(16,5) CHAR_VAR
            if(l_kappa)write(17,5) CHAR_VAR
            if(l_trans)write(18,5) CHAR_VAR
            write(19,5) CHAR_VAR


            write(8,10)
            write(7,15)  
            write(11,60)  
            write(12,65)  
            write(13,66)  
            write(14,67)  
            if(l_pol)write(15,68) 
            if(l_phase)write(16,69)  
            if(l_kappa)write(17,70)
            if(l_trans)write(18,71)
            write(19,72)

 
c
            write(8,100) (LOG_CONSTITUENT(i),i=1,N_CONSTITUENTS-1)
            write(7,100) (LOG_CONSTITUENT(i),i=1,N_CONSTITUENTS-1)
            write(11,100) (LOG_CONSTITUENT(i),i=1,N_CONSTITUENTS-1)
            write(12,100) (LOG_CONSTITUENT(i),i=1,N_CONSTITUENTS-1)
            write(13,100) (LOG_CONSTITUENT(i),i=1,N_CONSTITUENTS-1)
            write(14,100) (LOG_CONSTITUENT(i),i=1,N_CONSTITUENTS-1)
        if(l_pol)write(15,100) (LOG_CONSTITUENT(i),i=1,N_CONSTITUENTS-1)
        if(l_phase)write(15,100) 
     &             (LOG_CONSTITUENT(i),i=1,N_CONSTITUENTS-1)
        if(l_kappa)write(15,100) 
     &             (LOG_CONSTITUENT(i),i=1,N_CONSTITUENTS-1)
        if(l_trans)write(15,100) 
     &             (LOG_CONSTITUENT(i),i=1,N_CONSTITUENTS-1)
            write(19,100) (LOG_CONSTITUENT(i),i=1,N_CONSTITUENTS-1)
c
            write(8,105) SOIL_SURFACE,WATER_SURFACE,
     &                   SNOW_SURFACE,ICE_SURFACE
            write(7,105) SOIL_SURFACE,WATER_SURFACE,
     &                   SNOW_SURFACE,ICE_SURFACE
            write(11,105) SOIL_SURFACE,WATER_SURFACE,
     &                    SNOW_SURFACE,ICE_SURFACE
            write(12,105) SOIL_SURFACE,WATER_SURFACE,
     &                    SNOW_SURFACE,ICE_SURFACE
            write(13,105) SOIL_SURFACE,WATER_SURFACE,
     &                    SNOW_SURFACE,ICE_SURFACE
            write(14,105) SOIL_SURFACE,WATER_SURFACE,
     &                    SNOW_SURFACE,ICE_SURFACE
            if(l_pol) write(15,105) SOIL_SURFACE,WATER_SURFACE,
     &                    SNOW_SURFACE,ICE_SURFACE

            if(l_phase) write(15,105) SOIL_SURFACE,WATER_SURFACE,
     &                    SNOW_SURFACE,ICE_SURFACE

            if(l_kappa) write(15,105) SOIL_SURFACE,WATER_SURFACE,
     &                    SNOW_SURFACE,ICE_SURFACE

            if(l_trans) write(15,105) SOIL_SURFACE,WATER_SURFACE,
     &                    SNOW_SURFACE,ICE_SURFACE

            write(19,105) SOIL_SURFACE,WATER_SURFACE,
     &                    SNOW_SURFACE,ICE_SURFACE
c
            write(7,107) (LOG_EPS_TABLE(j),j=1,N_CONSTITUENTS)
            write(8,107) (LOG_EPS_TABLE(j),j=1,N_CONSTITUENTS)
            write(11,107) (LOG_EPS_TABLE(j),j=1,N_CONSTITUENTS)
            write(12,107) (LOG_EPS_TABLE(j),j=1,N_CONSTITUENTS)
            write(13,107) (LOG_EPS_TABLE(j),j=1,N_CONSTITUENTS)
            write(14,107) (LOG_EPS_TABLE(j),j=1,N_CONSTITUENTS)
            if(l_pol)write(15,107) (LOG_EPS_TABLE(j),j=1,N_CONSTITUENTS)
            if(l_phase)write(15,107) 
     &                 (LOG_EPS_TABLE(j),j=1,N_CONSTITUENTS)
            if(l_kappa)write(15,107) 
     &                 (LOG_EPS_TABLE(j),j=1,N_CONSTITUENTS)
            if(l_trans)write(15,107) 
     &                 (LOG_EPS_TABLE(j),j=1,N_CONSTITUENTS)
            write(19,107) (LOG_EPS_TABLE(j),j=1,N_CONSTITUENTS)
 
C
            IF(I_PDF_LEAF.EQ.1)THEN
                C_PDF_LF = '0.5*SIN(Theta_d) - Uniform'
            ELSE IF(I_PDF_LEAF.EQ.2)THEN
                C_PDF_LF = 'PLANOPHILE   -- BETA DIST.'
            ELSE IF(I_PDF_LEAF.EQ.3)THEN
                C_PDF_LF = 'ERECTOPHILE  -- BETA DIST.'
            ELSE IF(I_PDF_LEAF.EQ.4)THEN
                C_PDF_LF = 'PLAGIOPHILE  -- BETA DIST.'
            ELSE IF(I_PDF_LEAF.EQ.5)THEN
                C_PDF_LF = 'EXTREMOPHILE -- BETA DIST.'
            ELSE IF(I_PDF_LEAF.EQ.6)THEN
                C_PDF_LF = 'UNIFORM      -- BETA DIST.'
            ELSE IF(I_PDF_LEAF.EQ.7)THEN
                C_PDF_LF = 'SPHERICAL    -- BETA DIST.'
            ELSE
                C_PDF_LF = 'UNKNOWN'
            ENDIF

            IF(I_PDF_BR_1.EQ.1)THEN
                C_PDF_BR1 = '0.5*SIN(Theta_c) - Uniform'
            ELSE IF(I_PDF_BR_1.EQ.2)THEN
                C_PDF_BR1 = '(16/(3*pi))*SIN(2*THETA_c)**4'
            ELSE IF(I_PDF_BR_1.EQ.3)THEN
                C_PDF_BR1 = '(2/pi)*SIN(THETA_c)**2'
            ELSE IF(I_PDF_BR_1.EQ.4)THEN
                C_PDF_BR1 = '(3/4)*SIN(THETA_c)**3'
            ELSE IF(I_PDF_BR_1.EQ.5)THEN
                C_PDF_BR1 = '(8/(3*pi))*SIN(THETA_c)**4'
            ELSE IF(I_PDF_BR_1.EQ.6)THEN
                C_PDF_BR1 = '(15/16)**SIN(THETA_c)**5'
            ELSE IF(I_PDF_BR_1.EQ.7)THEN
                C_PDF_BR1 = '(3.2/pi)*SIN(THETA_c)**6'
            ELSE IF(I_PDF_BR_1.EQ.8)THEN
                C_PDF_BR1 = '(1.093750)*SIN(THETA_c)**7'
            ELSE IF(I_PDF_BR_1.EQ.9)THEN
                C_PDF_BR1 = '(1.164105)*SIN(THETA_c)**8'
            ELSE IF(I_PDF_BR_1.EQ.10)THEN
                C_PDF_BR1 = '(1.230469)*SIN(THETA_c)**9'
            ELSE IF(I_PDF_BR_1.EQ.11)THEN
                C_PDF_BR1 = '(1.230469)*SIN(THETA_c-30)**9'
            ELSE IF(I_PDF_BR_1.EQ.12)THEN
                C_PDF_BR1 = '(1.23)*SIN(THETA_c+60)**9'
            ELSE IF(I_PDF_BR_1.EQ.13)THEN
                C_PDF_BR1 = '(1.81)*SIN(THETA_c)**20'
            ELSE IF(I_PDF_BR_1.EQ.14)THEN
                C_PDF_BR1 = '(1.23)*SIN(THETA_c+30)**20'
            ELSE IF(I_PDF_BR_1.EQ.17)THEN
                C_PDF_BR1 = '(1.62)*SIN(THETA_c+64)**16'
            ELSE IF(I_PDF_BR_1.EQ.18)THEN
                C_PDF_BR1 = '(1.16)*SIN(THETA_c+41)**8'
            ELSE
                C_PDF_BR1 = 'UNKNOWN'
            ENDIF

            IF(I_PDF_BR_2.EQ.1)THEN
                C_PDF_BR2 = '0.5*SIN(Theta_c) - Uniform'
            ELSE IF(I_PDF_BR_2.EQ.2)THEN
                C_PDF_BR2 = '(16/(3*pi))*SIN(2*THETA_c)**4'
            ELSE IF(I_PDF_BR_2.EQ.3)THEN
                C_PDF_BR2 = '(2/pi)*SIN(THETA_c)**2'
            ELSE IF(I_PDF_BR_2.EQ.4)THEN
                C_PDF_BR2 = '(3/4)*SIN(THETA_c)**3'
            ELSE IF(I_PDF_BR_2.EQ.5)THEN
                C_PDF_BR2 = '(8/(3*pi))*SIN(THETA_c)**4'
            ELSE IF(I_PDF_BR_2.EQ.6)THEN
                C_PDF_BR2 = '(15/16)**SIN(THETA_c)**5'
            ELSE IF(I_PDF_BR_2.EQ.7)THEN
                C_PDF_BR2 = '(3.2/pi)*SIN(THETA_c)**6'
            ELSE IF(I_PDF_BR_2.EQ.8)THEN
                C_PDF_BR2 = '(1.093750)*SIN(THETA_c)**7'
            ELSE IF(I_PDF_BR_2.EQ.9)THEN
                C_PDF_BR2 = '(1.164105)*SIN(THETA_c)**8'
            ELSE IF(I_PDF_BR_2.EQ.10)THEN
                C_PDF_BR2 = '(1.230469)*SIN(THETA_c)**9'
            ELSE IF(I_PDF_BR_2.EQ.11)THEN
                C_PDF_BR2 = '(1.230469)*SIN(THETA_c-30)**9'
            ELSE IF(I_PDF_BR_2.EQ.12)THEN
                C_PDF_BR2 = '(1.23)*SIN(THETA_c+60)**9'
            ELSE IF(I_PDF_BR_2.EQ.13)THEN
                C_PDF_BR2 = '(1.81)*SIN(THETA_c)**20'
            ELSE IF(I_PDF_BR_2.EQ.14)THEN
                C_PDF_BR2 = '(1.23)*SIN(THETA_c+30)**20'
            ELSE IF(I_PDF_BR_2.EQ.17)THEN
                C_PDF_BR2 = '(1.62)*SIN(THETA_c+64)**16'
            ELSE IF(I_PDF_BR_2.EQ.18)THEN
                C_PDF_BR2 = '(1.16)*SIN(THETA_c+41)**8'

            ELSE
                C_PDF_BR2 = 'UNKNOWN'
            ENDIF

            IF(I_PDF_BR_3.EQ.1)THEN
                C_PDF_BR3 = '0.5*SIN(Theta_c) - Uniform'
            ELSE IF(I_PDF_BR_3.EQ.2)THEN
                C_PDF_BR3 = '(16/(3*pi))*SIN(2*THETA_c)**4'
            ELSE IF(I_PDF_BR_3.EQ.3)THEN
                C_PDF_BR3 = '(2/pi)*SIN(THETA_c)**2'
            ELSE IF(I_PDF_BR_3.EQ.4)THEN
                C_PDF_BR3 = '(3/4)*SIN(THETA_c)**3'
            ELSE IF(I_PDF_BR_3.EQ.5)THEN
                C_PDF_BR3 = '(8/(3*pi))*SIN(THETA_c)**4'
            ELSE IF(I_PDF_BR_3.EQ.6)THEN
                C_PDF_BR3 = '(15/16)**SIN(THETA_c)**5'
            ELSE IF(I_PDF_BR_3.EQ.7)THEN
                C_PDF_BR3 = '(3.2/pi)*SIN(THETA_c)**6'
            ELSE IF(I_PDF_BR_3.EQ.8)THEN
                C_PDF_BR3 = '(1.093750)*SIN(THETA_c)**7'
            ELSE IF(I_PDF_BR_3.EQ.9)THEN
                C_PDF_BR3 = '(1.164105)*SIN(THETA_c)**8'
            ELSE IF(I_PDF_BR_3.EQ.10)THEN
                C_PDF_BR3 = '(1.230469)*SIN(THETA_c)**9'
            ELSE IF(I_PDF_BR_3.EQ.11)THEN
                C_PDF_BR3 = '(1.230469)*SIN(THETA_c-30)**9'
            ELSE IF(I_PDF_BR_3.EQ.12)THEN
                C_PDF_BR3 = '(1.23)*SIN(THETA_c+60)**9'
            ELSE IF(I_PDF_BR_3.EQ.13)THEN
                C_PDF_BR3 = '(1.81)*SIN(THETA_c)**20'
            ELSE IF(I_PDF_BR_3.EQ.14)THEN
                C_PDF_BR3 = '(1.23)*SIN(THETA_c+30)**20'
            ELSE IF(I_PDF_BR_3.EQ.17)THEN
                C_PDF_BR3 = '(1.62)*SIN(THETA_c+64)**16'
            ELSE IF(I_PDF_BR_3.EQ.18)THEN
                C_PDF_BR3 = '(1.16)*SIN(THETA_c+41)**8'

            ELSE
                C_PDF_BR3 = 'UNKNOWN'
            ENDIF

            IF(I_PDF_BR_4.EQ.1)THEN
                C_PDF_BR4 = '0.5*SIN(Theta_c) - Uniform'
            ELSE IF(I_PDF_BR_4.EQ.2)THEN
                C_PDF_BR4 = '(16/(3*pi))*SIN(2*THETA_c)**4'
            ELSE IF(I_PDF_BR_4.EQ.3)THEN
                C_PDF_BR4 = '(2/pi)*SIN(THETA_c)**2'
            ELSE IF(I_PDF_BR_4.EQ.4)THEN
                C_PDF_BR4 = '(3/4)*SIN(THETA_c)**3'
            ELSE IF(I_PDF_BR_4.EQ.5)THEN
                C_PDF_BR4 = '(8/(3*pi))*SIN(THETA_c)**4'
            ELSE IF(I_PDF_BR_4.EQ.6)THEN
                C_PDF_BR4 = '(15/16)**SIN(THETA_c)**5'
            ELSE IF(I_PDF_BR_4.EQ.7)THEN
                C_PDF_BR4 = '(3.2/pi)*SIN(THETA_c)**6'
            ELSE IF(I_PDF_BR_4.EQ.8)THEN
                C_PDF_BR4 = '(1.093750)*SIN(THETA_c)**7'
            ELSE IF(I_PDF_BR_4.EQ.9)THEN
                C_PDF_BR4 = '(1.164105)*SIN(THETA_c)**8'
            ELSE IF(I_PDF_BR_4.EQ.10)THEN
                C_PDF_BR4 = '(1.230469)*SIN(THETA_c)**9'
            ELSE IF(I_PDF_BR_4.EQ.11)THEN
                C_PDF_BR4 = '(1.230469)*SIN(THETA_c-30)**9'
            ELSE IF(I_PDF_BR_4.EQ.12)THEN
                C_PDF_BR4 = '(1.23)*SIN(THETA_c+60)**9'
            ELSE IF(I_PDF_BR_4.EQ.13)THEN
                C_PDF_BR4 = '(1.81)*SIN(THETA_c)**20'
            ELSE IF(I_PDF_BR_4.EQ.14)THEN
                C_PDF_BR4 = '(1.23)*SIN(THETA_c+30)**20'
            ELSE IF(I_PDF_BR_4.EQ.17)THEN
                C_PDF_BR4 = '(1.62)*SIN(THETA_c+64)**16'
            ELSE IF(I_PDF_BR_4.EQ.18)THEN
                C_PDF_BR4 = '(1.16)*SIN(THETA_c+41)**8'

            ELSE
                C_PDF_BR4 = 'UNKNOWN'
            ENDIF

            IF(I_PDF_BR_5.EQ.1)THEN
                C_PDF_BR5 = '0.5*SIN(Theta_c) - Uniform'
            ELSE IF(I_PDF_BR_5.EQ.2)THEN
                C_PDF_BR5 = '(16/(3*pi))*SIN(2*THETA_c)**4'
            ELSE IF(I_PDF_BR_5.EQ.3)THEN
                C_PDF_BR5 = '(2/pi)*SIN(THETA_c)**2'
            ELSE IF(I_PDF_BR_5.EQ.4)THEN
                C_PDF_BR5 = '(3/4)*SIN(THETA_c)**3'
            ELSE IF(I_PDF_BR_5.EQ.5)THEN
                C_PDF_BR5 = '(8/(3*pi))*SIN(THETA_c)**4'
            ELSE IF(I_PDF_BR_5.EQ.6)THEN
                C_PDF_BR5 = '(15/16)**SIN(THETA_c)**5'
            ELSE IF(I_PDF_BR_5.EQ.7)THEN
                C_PDF_BR5 = '(3.2/pi)*SIN(THETA_c)**6'
            ELSE IF(I_PDF_BR_5.EQ.8)THEN
                C_PDF_BR5 = '(1.093750)*SIN(THETA_c)**7'
            ELSE IF(I_PDF_BR_5.EQ.9)THEN
                C_PDF_BR5 = '(1.164105)*SIN(THETA_c)**8'
            ELSE IF(I_PDF_BR_5.EQ.10)THEN
                C_PDF_BR5 = '(1.230469)*SIN(THETA_c)**9'
            ELSE IF(I_PDF_BR_5.EQ.11)THEN
                C_PDF_BR5 = '(1.230469)*SIN(THETA_c-30)**9'
            ELSE IF(I_PDF_BR_5.EQ.12)THEN
                C_PDF_BR5 = '(1.23)*SIN(THETA_c+60)**9'
            ELSE IF(I_PDF_BR_5.EQ.13)THEN
                C_PDF_BR5 = '(1.81)*SIN(THETA_c)**20'
            ELSE IF(I_PDF_BR_5.EQ.14)THEN
                C_PDF_BR5 = '(1.23)*SIN(THETA_c+30)**20'
            ELSE IF(I_PDF_BR_5.EQ.17)THEN
                C_PDF_BR5 = '(1.62)*SIN(THETA_c+64)**16'
            ELSE IF(I_PDF_BR_5.EQ.18)THEN
                C_PDF_BR5 = '(1.16)*SIN(THETA_c+41)**8'

            ELSE
                C_PDF_BR5 = 'UNKNOWN'
            ENDIF

            IF(I_PDF_BR_6.EQ.1)THEN
                C_PDF_BR6 = '0.5*SIN(Theta_c) - Uniform'
            ELSE IF(I_PDF_BR_6.EQ.2)THEN
                C_PDF_BR6 = '(16/(3*pi))*SIN(2*THETA_c)**4'
            ELSE IF(I_PDF_BR_6.EQ.3)THEN
                C_PDF_BR6 = '(2/pi)*SIN(THETA_c)**2'
            ELSE IF(I_PDF_BR_6.EQ.4)THEN
                C_PDF_BR6 = '(3/4)*SIN(THETA_c)**3'
            ELSE IF(I_PDF_BR_6.EQ.5)THEN
                C_PDF_BR6 = '(8/(3*pi))*SIN(THETA_c)**4'
            ELSE IF(I_PDF_BR_6.EQ.6)THEN
                C_PDF_BR6 = '(15/16)**SIN(THETA_c)**5'
            ELSE IF(I_PDF_BR_6.EQ.7)THEN
                C_PDF_BR6 = '(3.2/pi)*SIN(THETA_c)**6'
            ELSE IF(I_PDF_BR_6.EQ.8)THEN
                C_PDF_BR6 = '(1.093750)*SIN(THETA_c)**7'
            ELSE IF(I_PDF_BR_6.EQ.9)THEN
                C_PDF_BR6 = '(1.164105)*SIN(THETA_c)**8'
            ELSE IF(I_PDF_BR_6.EQ.10)THEN
                C_PDF_BR6 = '(1.230469)*SIN(THETA_c)**9'
            ELSE IF(I_PDF_BR_6.EQ.11)THEN
                C_PDF_BR6 = '(1.230469)*SIN(THETA_c-30)**9'
            ELSE IF(I_PDF_BR_6.EQ.12)THEN
                C_PDF_BR6 = '(1.23)*SIN(THETA_c+60)**9'
            ELSE IF(I_PDF_BR_6.EQ.13)THEN
                C_PDF_BR6 = '(1.81)*SIN(THETA_c)**20'
            ELSE IF(I_PDF_BR_6.EQ.14)THEN
                C_PDF_BR6 = '(1.23)*SIN(THETA_c+30)**20'
            ELSE IF(I_PDF_BR_6.EQ.17)THEN
                C_PDF_BR6 = '(1.62)*SIN(THETA_c+64)**16'
            ELSE IF(I_PDF_BR_6.EQ.18)THEN
                C_PDF_BR6 = '(1.16)*SIN(THETA_c+41)**8'

            ELSE
                C_PDF_BR6 = 'UNKNOWN'
            ENDIF

            IF(I_PDF_NDL.EQ.1)THEN
                C_PDF_NDL = '0.5*SIN(Theta_n) - Uniform'
            ELSE IF(I_PDF_NDL.EQ.2)THEN
                C_PDF_NDL = '(16/(3*pi))*SIN(2*THETA_c)**4'
            ELSE IF(I_PDF_NDL.EQ.3)THEN
                C_PDF_NDL = '(2/pi)*SIN(THETA_c)**2'
            ELSE IF(I_PDF_NDL.EQ.4)THEN
                C_PDF_NDL = '(3/4)*SIN(THETA_c)**3'
            ELSE IF(I_PDF_NDL.EQ.5)THEN
                C_PDF_NDL = '(8/(3*pi))*SIN(THETA_c)**4'
            ELSE IF(I_PDF_NDL.EQ.6)THEN
                C_PDF_NDL = '(15/16)**SIN(THETA_c)**5'
            ELSE IF(I_PDF_NDL.EQ.7)THEN
                C_PDF_NDL = '(3.2/pi)*SIN(THETA_c)**6'
            ELSE IF(I_PDF_NDL.EQ.8)THEN
                C_PDF_NDL = '(1.093750)*SIN(THETA_c)**7'
            ELSE IF(I_PDF_NDL.EQ.9)THEN
                C_PDF_NDL = '(1.164105)*SIN(THETA_c)**8'
            ELSE IF(I_PDF_NDL.EQ.10)THEN
                C_PDF_NDL = '(1.230469)*SIN(THETA_c)**9'
            ELSE IF(I_PDF_NDL.EQ.11)THEN
                C_PDF_NDL = '(1.230469)*SIN(THETA_c-30)**9'
            ELSE IF(I_PDF_NDL.EQ.12)THEN
                C_PDF_NDL = '(1.23)*SIN(THETA_c+60)**9'
            ELSE IF(I_PDF_NDL.EQ.13)THEN
                C_PDF_NDL = '(1.81)*SIN(THETA_c)**20'
            ELSE IF(I_PDF_NDL.EQ.14)THEN
                C_PDF_NDL = '(1.23)*SIN(THETA_c+30)**20'
            ELSE IF(I_PDF_NDL.EQ.17)THEN
                C_PDF_NDL = '(1.62)*COS(THETA_c)**16'
            ELSE
                C_PDF_NDL = 'UNKNOWN'
            ENDIF

            IF(I_PDF_TRUNK.EQ.1)THEN
                C_PDF_TR = 'VERTICAL'
            ELSE IF(I_PDF_TRUNK.EQ.2)THEN 
                C_PDF_TR = '(1/0.85903)*COS(THETA_c)**8'
            ELSE IF(I_PDF_TRUNK.EQ.3)THEN 
                C_PDF_TR = '(1/0.98175)*COS(THETA_c)**6'
            ELSE IF(I_PDF_TRUNK.EQ.4)THEN 
                C_PDF_TR = '8*EXP(-8*THETA_c)'
            ELSE
                C_PDF_TR = 'UNKNOWN'
            ENDIF
C
            IF(I_PDF_LEAF_SIZE.EQ.0)THEN
                C_PDF_LF_SZ = 'DEFAULT VALUES'
            ELSE
                C_PDF_LF_SZ = 'UNKNOWN'
            ENDIF

            IF(I_PDF_BR_1_SIZE.EQ.0)THEN
                C_PDF_BR1_SZ = 'DEFAULT VALUES'
            ELSE
                C_PDF_BR1_SZ = 'UNKNOWN'
            ENDIF

            IF(I_PDF_BR_2_SIZE.EQ.0)THEN
                C_PDF_BR2_SZ = 'DEFAULT VALUES'
            ELSE
                C_PDF_BR2_SZ = 'UNKNOWN'
            ENDIF

            IF(I_PDF_BR_3_SIZE.EQ.0)THEN
                C_PDF_BR3_SZ = 'DEFAULT VALUES'
            ELSE
                C_PDF_BR3_SZ = 'UNKNOWN'
            ENDIF

            IF(I_PDF_BR_4_SIZE.EQ.0)THEN
                C_PDF_BR4_SZ = 'DEFAULT VALUES'
            ELSE
                C_PDF_BR4_SZ = 'UNKNOWN'
            ENDIF

            IF(I_PDF_BR_5_SIZE.EQ.0)THEN
                C_PDF_BR5_SZ = 'DEFAULT VALUES'
            ELSE
                C_PDF_BR5_SZ = 'UNKNOWN'
            ENDIF

            IF(I_PDF_BR_6_SIZE.EQ.0)THEN
                C_PDF_BR6_SZ = 'DEFAULT VALUES'
            ELSE
                C_PDF_BR6_SZ = 'UNKNOWN'
            ENDIF

            IF(I_PDF_NDL_SIZE.EQ.0)THEN
                C_PDF_NDL_SZ = 'DEFAULT VALUES'
            ELSE
                C_PDF_NDL_SZ = 'UNKNOWN'
            ENDIF

            IF(I_PDF_TRUNK_SIZE.EQ.0)THEN
                C_PDF_TR_SZ = 'DEFAULT VALUES'
            ELSE IF(I_PDF_TRUNK_SIZE.EQ.1)THEN
                C_PDF_TR_SZ = 'HISTOGRAM DATA'
            ELSE
                C_PDF_TR_SZ = 'UNKNOWN'
            ENDIF

            WRITE(7,199)
            WRITE(8,199)
            WRITE(11,199)
            WRITE(12,199)
            WRITE(13,199)
            WRITE(14,199)
            if(l_pol)WRITE(15,199)
            if(l_phase)WRITE(16,199)
            if(l_kappa)WRITE(17,199)
            if(l_trans)WRITE(18,199)
            WRITE(19,199)

            if(LOG_CONSTITUENT(TRUNK)) THEN
              WRITE(7,200) I_PDF_TRUNK,C_PDF_TR, I_PDF_TRUNK_SIZE,
     &                     C_PDF_TR_SZ
              WRITE(8,200) I_PDF_TRUNK,C_PDF_TR, I_PDF_TRUNK_SIZE,
     &                     C_PDF_TR_SZ
              WRITE(11,200) I_PDF_TRUNK,C_PDF_TR, I_PDF_TRUNK_SIZE,
     &                     C_PDF_TR_SZ
              WRITE(12,200) I_PDF_TRUNK,C_PDF_TR, I_PDF_TRUNK_SIZE,
     &                     C_PDF_TR_SZ
              WRITE(13,200) I_PDF_TRUNK,C_PDF_TR, I_PDF_TRUNK_SIZE,
     &                     C_PDF_TR_SZ
              WRITE(14,200) I_PDF_TRUNK,C_PDF_TR, I_PDF_TRUNK_SIZE,
     &                     C_PDF_TR_SZ
              if(l_pol)WRITE(15,200) I_PDF_TRUNK, C_PDF_TR, 
     &                     I_PDF_TRUNK_SIZE, C_PDF_TR_SZ
              if(l_phase)WRITE(16,200) I_PDF_TRUNK, C_PDF_TR, 
     &                     I_PDF_TRUNK_SIZE, C_PDF_TR_SZ
              if(l_kappa)WRITE(17,200) I_PDF_TRUNK, C_PDF_TR, 
     &                     I_PDF_TRUNK_SIZE, C_PDF_TR_SZ
              if(l_trans)WRITE(18,200) I_PDF_TRUNK, C_PDF_TR, 
     &                     I_PDF_TRUNK_SIZE, C_PDF_TR_SZ
              WRITE(19,200) I_PDF_TRUNK,C_PDF_TR, I_PDF_TRUNK_SIZE,
     &                     C_PDF_TR_SZ
            ENDIF
            if(LOG_CONSTITUENT(LEAVES)) THEN
              WRITE(7,201) I_PDF_LEAF,C_PDF_LF, I_PDF_LEAF_SIZE,
     &                     C_PDF_LF_SZ
              WRITE(8,201) I_PDF_LEAF,C_PDF_LF, I_PDF_LEAF_SIZE,
     &                     C_PDF_LF_SZ
              WRITE(11,201) I_PDF_LEAF,C_PDF_LF, I_PDF_LEAF_SIZE,
     &                     C_PDF_LF_SZ
              WRITE(12,201) I_PDF_LEAF,C_PDF_LF, I_PDF_LEAF_SIZE,
     &                     C_PDF_LF_SZ
              WRITE(13,201) I_PDF_LEAF,C_PDF_LF, I_PDF_LEAF_SIZE,
     &                     C_PDF_LF_SZ
              WRITE(14,201) I_PDF_LEAF,C_PDF_LF, I_PDF_LEAF_SIZE,
     &                     C_PDF_LF_SZ
              if(l_pol) WRITE(15,201) I_PDF_LEAF,C_PDF_LF, 
     &                     I_PDF_LEAF_SIZE, C_PDF_LF_SZ
              if(l_phase)WRITE(16,201) I_PDF_LEAF,C_PDF_LF, 
     &                     I_PDF_LEAF_SIZE, C_PDF_LF_SZ
              if(l_kappa) WRITE(17,201) I_PDF_LEAF,C_PDF_LF, 
     &                     I_PDF_LEAF_SIZE, C_PDF_LF_SZ
              if(l_trans)WRITE(18,201) I_PDF_LEAF,C_PDF_LF, 
     &                     I_PDF_LEAF_SIZE, C_PDF_LF_SZ
               WRITE(19,201) I_PDF_LEAF,C_PDF_LF, 
     &                     I_PDF_LEAF_SIZE, C_PDF_LF_SZ
            ENDIF
            if(LOG_CONSTITUENT(NEEDLE)) THEN
              WRITE(7,202) I_PDF_NDL,C_PDF_NDL, I_PDF_NDL_SIZE,
     &                     C_PDF_NDL_SZ
              WRITE(8,202) I_PDF_NDL,C_PDF_NDL, I_PDF_NDL_SIZE,
     &                     C_PDF_NDL_SZ
              WRITE(11,202) I_PDF_NDL,C_PDF_NDL, I_PDF_NDL_SIZE,
     &                     C_PDF_NDL_SZ
              WRITE(12,202) I_PDF_NDL,C_PDF_NDL, I_PDF_NDL_SIZE,
     &                     C_PDF_NDL_SZ
              WRITE(13,202) I_PDF_NDL,C_PDF_NDL, I_PDF_NDL_SIZE,
     &                     C_PDF_NDL_SZ
              WRITE(14,202) I_PDF_NDL,C_PDF_NDL, I_PDF_NDL_SIZE,
     &                     C_PDF_NDL_SZ
              if(l_pol) WRITE(15,202) I_PDF_NDL,C_PDF_NDL,
     &                     I_PDF_NDL_SIZE, C_PDF_NDL_SZ
              if(l_phase)WRITE(16,202) I_PDF_NDL,C_PDF_NDL,
     &                     I_PDF_NDL_SIZE, C_PDF_NDL_SZ
              if(l_kappa) WRITE(17,202) I_PDF_NDL,C_PDF_NDL,
     &                     I_PDF_NDL_SIZE, C_PDF_NDL_SZ
              if(l_trans)WRITE(18,202) I_PDF_NDL,C_PDF_NDL,
     &                     I_PDF_NDL_SIZE, C_PDF_NDL_SZ
              WRITE(19,202) I_PDF_NDL,C_PDF_NDL,
     &                     I_PDF_NDL_SIZE, C_PDF_NDL_SZ
            ENDIF
            if(LOG_CONSTITUENT(BR1)) THEN
              WRITE(7,203) I_PDF_BR_1, C_PDF_BR1, I_PDF_BR_1_SIZE,
     &                     C_PDF_BR1_SZ
              WRITE(8,203) I_PDF_BR_1, C_PDF_BR1, I_PDF_BR_1_SIZE,
     &                     C_PDF_BR1_SZ
              WRITE(11,203) I_PDF_BR_1, C_PDF_BR1, I_PDF_BR_1_SIZE,
     &                     C_PDF_BR1_SZ
              WRITE(12,203) I_PDF_BR_1, C_PDF_BR1, I_PDF_BR_1_SIZE,
     &                     C_PDF_BR1_SZ
              WRITE(13,203) I_PDF_BR_1, C_PDF_BR1, I_PDF_BR_1_SIZE,
     &                     C_PDF_BR1_SZ
              WRITE(14,203) I_PDF_BR_1, C_PDF_BR1, I_PDF_BR_1_SIZE,
     &                     C_PDF_BR1_SZ
              if(l_pol) WRITE(15,203) I_PDF_BR_1, C_PDF_BR1,
     &                      I_PDF_BR_1_SIZE, C_PDF_BR1_SZ
              if(l_phase) WRITE(16,203) I_PDF_BR_1, C_PDF_BR1,
     &                      I_PDF_BR_1_SIZE, C_PDF_BR1_SZ
              if(l_kappa) WRITE(17,203) I_PDF_BR_1, C_PDF_BR1,
     &                      I_PDF_BR_1_SIZE, C_PDF_BR1_SZ
              if(l_trans) WRITE(18,203) I_PDF_BR_1, C_PDF_BR1,
     &                      I_PDF_BR_1_SIZE, C_PDF_BR1_SZ
              WRITE(19,203) I_PDF_BR_1, C_PDF_BR1,
     &                      I_PDF_BR_1_SIZE, C_PDF_BR1_SZ
            ENDIF
            if(LOG_CONSTITUENT(BR2)) THEN
              WRITE(7,204) I_PDF_BR_2, C_PDF_BR2, I_PDF_BR_2_SIZE,
     &                     C_PDF_BR2_SZ
              WRITE(8,204) I_PDF_BR_2, C_PDF_BR2, I_PDF_BR_2_SIZE,
     &                     C_PDF_BR2_SZ
              WRITE(11,204) I_PDF_BR_2, C_PDF_BR2, I_PDF_BR_2_SIZE,
     &                     C_PDF_BR2_SZ
              WRITE(12,204) I_PDF_BR_2, C_PDF_BR2, I_PDF_BR_2_SIZE,
     &                     C_PDF_BR2_SZ
              WRITE(13,204) I_PDF_BR_2, C_PDF_BR2, I_PDF_BR_2_SIZE,
     &                     C_PDF_BR2_SZ
              WRITE(14,204) I_PDF_BR_2, C_PDF_BR2, I_PDF_BR_2_SIZE,
     &                     C_PDF_BR2_SZ
              if(l_pol) WRITE(15,204) I_PDF_BR_2, C_PDF_BR2,
     &                      I_PDF_BR_2_SIZE, C_PDF_BR2_SZ
              if(l_phase)WRITE(16,204) I_PDF_BR_2, C_PDF_BR2,
     &                      I_PDF_BR_2_SIZE, C_PDF_BR2_SZ
              if(l_kappa) WRITE(17,204) I_PDF_BR_2, C_PDF_BR2,
     &                      I_PDF_BR_2_SIZE, C_PDF_BR2_SZ
              if(l_trans)WRITE(18,204) I_PDF_BR_2, C_PDF_BR2,
     &                      I_PDF_BR_2_SIZE, C_PDF_BR2_SZ
              WRITE(19,204) I_PDF_BR_2, C_PDF_BR2,
     &                      I_PDF_BR_2_SIZE, C_PDF_BR2_SZ
            ENDIF
            if(LOG_CONSTITUENT(BR3)) THEN
              WRITE(7,205) I_PDF_BR_3, C_PDF_BR3, I_PDF_BR_3_SIZE,
     &                     C_PDF_BR3_SZ
              WRITE(8,205) I_PDF_BR_3, C_PDF_BR3, I_PDF_BR_3_SIZE,
     &                     C_PDF_BR3_SZ
              WRITE(11,205) I_PDF_BR_3, C_PDF_BR3, I_PDF_BR_3_SIZE,
     &                     C_PDF_BR3_SZ
              WRITE(12,205) I_PDF_BR_3, C_PDF_BR3, I_PDF_BR_3_SIZE,
     &                     C_PDF_BR3_SZ
              WRITE(13,205) I_PDF_BR_3, C_PDF_BR3, I_PDF_BR_3_SIZE,
     &                     C_PDF_BR3_SZ
              WRITE(14,205) I_PDF_BR_3, C_PDF_BR3, I_PDF_BR_3_SIZE,
     &                     C_PDF_BR3_SZ
              if(l_pol) WRITE(15,205) I_PDF_BR_3, C_PDF_BR3,
     &                      I_PDF_BR_3_SIZE, C_PDF_BR3_SZ
              if(l_phase) WRITE(16,205) I_PDF_BR_3, C_PDF_BR3,
     &                      I_PDF_BR_3_SIZE, C_PDF_BR3_SZ
               if(l_kappa) WRITE(17,205) I_PDF_BR_3, C_PDF_BR3,
     &                      I_PDF_BR_3_SIZE, C_PDF_BR3_SZ
              if(l_trans) WRITE(18,205) I_PDF_BR_3, C_PDF_BR3,
     &                      I_PDF_BR_3_SIZE, C_PDF_BR3_SZ
              WRITE(19,205) I_PDF_BR_3, C_PDF_BR3,
     &                      I_PDF_BR_3_SIZE, C_PDF_BR3_SZ
           ENDIF
            if(LOG_CONSTITUENT(BR4)) THEN
              WRITE(7,206) I_PDF_BR_4, C_PDF_BR4, I_PDF_BR_4_SIZE,
     &                     C_PDF_BR4_SZ
              WRITE(8,206) I_PDF_BR_4, C_PDF_BR4, I_PDF_BR_4_SIZE,
     &                     C_PDF_BR4_SZ
              WRITE(11,206) I_PDF_BR_4, C_PDF_BR4, I_PDF_BR_4_SIZE,
     &                     C_PDF_BR4_SZ
              WRITE(12,206) I_PDF_BR_4, C_PDF_BR4, I_PDF_BR_4_SIZE,
     &                     C_PDF_BR4_SZ
              WRITE(13,206) I_PDF_BR_4, C_PDF_BR4, I_PDF_BR_4_SIZE,
     &                     C_PDF_BR4_SZ
              WRITE(14,206) I_PDF_BR_4, C_PDF_BR4, I_PDF_BR_4_SIZE,
     &                     C_PDF_BR4_SZ
              if(l_pol) WRITE(15,206) I_PDF_BR_4, C_PDF_BR4,
     &                      I_PDF_BR_4_SIZE, C_PDF_BR4_SZ
              if(l_phase) WRITE(16,205) I_PDF_BR_4, C_PDF_BR4,
     &                      I_PDF_BR_4_SIZE, C_PDF_BR4_SZ
               if(l_kappa) WRITE(17,205) I_PDF_BR_4, C_PDF_BR4,
     &                      I_PDF_BR_4_SIZE, C_PDF_BR4_SZ
              if(l_trans) WRITE(18,205) I_PDF_BR_4, C_PDF_BR4,
     &                      I_PDF_BR_4_SIZE, C_PDF_BR4_SZ
              WRITE(19,205) I_PDF_BR_4, C_PDF_BR4,
     &                      I_PDF_BR_4_SIZE, C_PDF_BR4_SZ
            ENDIF
            if(LOG_CONSTITUENT(BR5)) THEN
              WRITE(7,207) I_PDF_BR_5, C_PDF_BR5, I_PDF_BR_5_SIZE,
     &                     C_PDF_BR5_SZ
              WRITE(8,207) I_PDF_BR_5, C_PDF_BR5, I_PDF_BR_5_SIZE,
     &                     C_PDF_BR5_SZ
              WRITE(11,207) I_PDF_BR_5, C_PDF_BR5, I_PDF_BR_5_SIZE,
     &                     C_PDF_BR5_SZ
              WRITE(12,207) I_PDF_BR_5, C_PDF_BR5, I_PDF_BR_5_SIZE,
     &                     C_PDF_BR5_SZ
              WRITE(13,207) I_PDF_BR_5, C_PDF_BR5, I_PDF_BR_5_SIZE,
     &                     C_PDF_BR5_SZ
              WRITE(14,207) I_PDF_BR_5, C_PDF_BR5, I_PDF_BR_5_SIZE,
     &                     C_PDF_BR5_SZ
              if(l_pol) WRITE(15,207) I_PDF_BR_5, C_PDF_BR5,
     &                      I_PDF_BR_5_SIZE, C_PDF_BR5_SZ
              if(l_phase) WRITE(16,207) I_PDF_BR_5, C_PDF_BR5,
     &                      I_PDF_BR_5_SIZE, C_PDF_BR5_SZ
              if(l_kappa) WRITE(17,207) I_PDF_BR_5, C_PDF_BR5,
     &                      I_PDF_BR_5_SIZE, C_PDF_BR5_SZ
              if(l_trans) WRITE(18,207) I_PDF_BR_5, C_PDF_BR5,
     &                      I_PDF_BR_5_SIZE, C_PDF_BR5_SZ
              WRITE(19,207) I_PDF_BR_5, C_PDF_BR5,
     &                      I_PDF_BR_5_SIZE, C_PDF_BR5_SZ
            ENDIF
            if(LOG_CONSTITUENT(BR6)) THEN
              WRITE(7,208) I_PDF_BR_6, C_PDF_BR6, I_PDF_BR_6_SIZE,
     &                     C_PDF_BR6_SZ
              WRITE(8,208) I_PDF_BR_6, C_PDF_BR6, I_PDF_BR_6_SIZE,
     &                     C_PDF_BR6_SZ
              WRITE(11,208) I_PDF_BR_6, C_PDF_BR6, I_PDF_BR_6_SIZE,
     &                     C_PDF_BR6_SZ
              WRITE(12,208) I_PDF_BR_6, C_PDF_BR6, I_PDF_BR_6_SIZE,
     &                     C_PDF_BR6_SZ
              WRITE(13,208) I_PDF_BR_6, C_PDF_BR6, I_PDF_BR_6_SIZE,
     &                     C_PDF_BR6_SZ
              WRITE(14,208) I_PDF_BR_6, C_PDF_BR6, I_PDF_BR_6_SIZE,
     &                     C_PDF_BR6_SZ
              if(l_pol) WRITE(15,208) I_PDF_BR_6, C_PDF_BR6,
     &                      I_PDF_BR_6_SIZE, C_PDF_BR6_SZ
              if(l_phase) WRITE(16,208) I_PDF_BR_6, C_PDF_BR6,
     &                      I_PDF_BR_6_SIZE, C_PDF_BR6_SZ
              if(l_kappa) WRITE(17,208) I_PDF_BR_6, C_PDF_BR6,
     &                      I_PDF_BR_6_SIZE, C_PDF_BR6_SZ
              if(l_trans) WRITE(18,208) I_PDF_BR_6, C_PDF_BR6,
     &                      I_PDF_BR_6_SIZE, C_PDF_BR6_SZ
              WRITE(19,208) I_PDF_BR_6, C_PDF_BR6,
     &                      I_PDF_BR_6_SIZE, C_PDF_BR6_SZ
            ENDIF
C
            K = 0
            DO I=1,N_VARIABLES
                IF(STEP_VARIABLE(I))THEN
                    K = K+1
                    PARAM_OUT(k) = I
                    NEST_OUT(k) = INEST(I)
                    LOOP_OUT(k) = LOOP_NUM(I)
                ENDIF
            ENDDO
C
            DO I=1,K
             DO J=1,K
              IF(I.NE.J)THEN
               IF(NEST_OUT(I).EQ.NEST_OUT(J))THEN
                STEP_SET(I) = .TRUE.
               ENDIF
              ENDIF
             ENDDO
            ENDDO
C
C            IF(NEST_OUT(1).NE.NEST_OUT(2))THEN
C
                do i=1,k
                  nest_set(i) = .false.
                enddo
c
                icount = 0
900             do i=1,k
                    if(.not.nest_set(i))then
                        ipoint = i
                        do j=1,k
                            if(nest_out(j).lt.nest_out(ipoint))then
                                if(.not.nest_set(j))then
                                    ipoint = j
                                endif
                            endif
                        enddo
                        icount = icount+1
c                        IF(.NOT.STEP_SET(IPOINT))THEN
c                         nest_out(ipoint) = icount
c                        ENDIF
                        nest_set(ipoint) = .true.
                        if(icount.ge.k)then
                            goto 901
                        else
                            goto 900
                        endif
                    endif
                enddo
901             continue
c
            DO J=1,K
                CHECK_NEST(J) = .FALSE.
            ENDDO
            DO J=1,K
              ILOW = N_VARIABLES + 1
              DO I=1,K
                IF((NEST_OUT(I).LT.ILOW).AND.(.NOT.CHECK_NEST(I)))THEN
                    ILOW = NEST_OUT(I)
                ENDIF
              ENDDO
C
              DO I=1,K
                IF(NEST_OUT(I).EQ.ILOW)THEN
                    NEST_OUT(I) = J
                    CHECK_NEST(I) = .TRUE.
                ENDIF
              ENDDO
            ENDDO
C
            write(8,110) (PARAM_OUT(I), I=1,K)
            write(8,111)  (NEST_OUT(I), I=1,K)
            write(8,112) (LOOP_OUT(I), I=1,K)

            write(7,110) (PARAM_OUT(I), I=1,K)
            write(7,111)  (NEST_OUT(I), I=1,K)
            write(7,112) (LOOP_OUT(I), I=1,K)

            write(11,110) (PARAM_OUT(I), I=1,K)
            write(11,111)  (NEST_OUT(I), I=1,K)
            write(11,112) (LOOP_OUT(I), I=1,K)

            write(12,110) (PARAM_OUT(I), I=1,K)
            write(12,111)  (NEST_OUT(I), I=1,K)
            write(12,112) (LOOP_OUT(I), I=1,K)

            write(13,110) (PARAM_OUT(I), I=1,K)
            write(13,111)  (NEST_OUT(I), I=1,K)
            write(13,112) (LOOP_OUT(I), I=1,K)

            write(14,110) (PARAM_OUT(I), I=1,K)
            write(14,111)  (NEST_OUT(I), I=1,K)
            write(14,112) (LOOP_OUT(I), I=1,K)

            if(l_pol) then
                write(15,110) (PARAM_OUT(I), I=1,K)
                write(15,111)  (NEST_OUT(I), I=1,K)
                write(15,112) (LOOP_OUT(I), I=1,K)
            endif

            if(l_phase) then
                write(16,110) (PARAM_OUT(I), I=1,K)
                write(16,111)  (NEST_OUT(I), I=1,K)
                write(16,112) (LOOP_OUT(I), I=1,K)
            endif

            if(l_kappa) then
                write(17,110) (PARAM_OUT(I), I=1,K)
                write(17,111)  (NEST_OUT(I), I=1,K)
                write(17,112) (LOOP_OUT(I), I=1,K)
            endif

            if(l_trans) then
                write(18,110) (PARAM_OUT(I), I=1,K)
                write(18,111)  (NEST_OUT(I), I=1,K)
                write(18,112) (LOOP_OUT(I), I=1,K)
            endif

            write(19,110) (PARAM_OUT(I), I=1,K)
            write(19,111)  (NEST_OUT(I), I=1,K)
            write(19,112) (LOOP_OUT(I), I=1,K)

c
            open = .true.
            write_header = .true.
c              
        endif
c
        VARIABLE(1)  = THETA_DEGREES
        VARIABLE(2)  = FREQ_GHZ
        if(LOG_EPS_TABLE(10))then
            VARIABLE(3)  = real(epsilonrc(2))
        else
            VARIABLE(3)  = MV_SOIL
        endif

        VARIABLE(4)  = RMS_SOIL
        VARIABLE(5)  = LS_SOIL
        VARIABLE(6)  = SAND
        VARIABLE(7)  = CLAY
        VARIABLE(8)  = SALT
        if(LOG_EPS_TABLE(1))then
            VARIABLE(9)  = real(epsilonrc(5))
        else
            VARIABLE(9)  = MG_TRUNK
        endif

        VARIABLE(10) = DENSITY
        VARIABLE(11) = RHO_TRUNK
        VARIABLE(12) = CROWN_HGHT
        VARIABLE(13) = TRUNK_DIAM
        VARIABLE(14) = TRUNK_HGHT
        if(LOG_EPS_TABLE(8))then
            VARIABLE(15)  = real(epsilonrc(4))
        else
            VARIABLE(15) = MG_LEAF
        endif

        VARIABLE(16) = RHO_LEAF
        VARIABLE(17) = LEAF_DENS
        VARIABLE(18) = LEAF_DIAM
        VARIABLE(19) = LEAF_TAU

        if(LOG_EPS_TABLE(9))then
            VARIABLE(20)  = real(epsilonrc(3))
        else
            VARIABLE(20) = MG_NDL
        endif
        VARIABLE(21) = RHO_NDL
        VARIABLE(22) = NDL_DENS
        VARIABLE(23) = NDL_DIAM
        VARIABLE(24) = NDL_LNG
        if(LOG_EPS_TABLE(2))then
            VARIABLE(25)  = real(epsilonrc(6))
        else
            VARIABLE(25) = MG_BR1
        endif

        VARIABLE(26) = RHO_BR1
        VARIABLE(27) = BR1_DENS
        VARIABLE(28) = BR1_DIAM
        VARIABLE(29) = BR1_LNG
        VARIABLE(30) = T_SOIL
        VARIABLE(31) = T_VEG
        VARIABLE(32) = T_WATER

        if(LOG_EPS_TABLE(3))then
            VARIABLE(33)  = real(epsilonrc(7))
        else
            VARIABLE(33) = MG_BR2
        endif
        VARIABLE(34) = RHO_BR2
        VARIABLE(35) = BR2_DENS
        VARIABLE(36) = BR2_DIAM
        VARIABLE(37) = BR2_LNG
C
        if(LOG_EPS_TABLE(4))then
            VARIABLE(38)  = real(epsilonrc(9))
        else
            VARIABLE(38) = MG_BR3
        endif
        VARIABLE(39) = RHO_BR3
        VARIABLE(40) = BR3_DENS
        VARIABLE(41) = BR3_DIAM
        VARIABLE(42) = BR3_LNG
C
        if(LOG_EPS_TABLE(5))then
            VARIABLE(43)  = real(epsilonrc(10))
        else
            VARIABLE(43) = MG_BR4
        endif
        VARIABLE(44) = RHO_BR4
        VARIABLE(45) = BR4_DENS
        VARIABLE(46) = BR4_DIAM
        VARIABLE(47) = BR4_LNG
C
        if(LOG_EPS_TABLE(6))then
            VARIABLE(48)  = real(epsilonrc(11))
        else
            VARIABLE(48) = MG_BR5
        endif
        VARIABLE(48) = MG_BR5
        VARIABLE(49) = RHO_BR5
        VARIABLE(50) = BR5_DENS
        VARIABLE(51) = BR5_DIAM
        VARIABLE(52) = BR5_LNG
C
        if(LOG_EPS_TABLE(7))then
            VARIABLE(53)  = real(epsilonrc(12))
        else
            VARIABLE(53) = MG_BR6
        endif
        VARIABLE(53) = MG_BR6
        VARIABLE(54) = RHO_BR6
        VARIABLE(55) = BR6_DENS
        VARIABLE(56) = BR6_DIAM
        VARIABLE(57) = BR6_LNG
C
        VARIABLE(60) = T_SNOW
        VARIABLE(61) = real(epsilonrc(8))
C
        if (WRITE_HEADER) then
c
        if(STEP_LAST_TIME(1))then
            char =  '  Theta  '
            char2 = '  (deg)  '
            INT_VAR = 1
            CHAR_FORM = '(1x,f6.3,4x)'
        else if(STEP_LAST_TIME(2))then
            char  = 'Frequency'
            char2 = '  (GHz)  '
            INT_VAR = 2
            CHAR_FORM = '(1x,f6.3,4x)'
        else if(STEP_LAST_TIME(3))then
          if(LOG_EPS_TABLE(10))then
            char  = '  Soil   '
            char2 = ' RE(epsr)'
            INT_VAR = 3
            CHAR_FORM = '(1x,f6.3,4x)'
          else
            char  = ' Soil Mv '
            char2 = '         '
            INT_VAR = 3
            CHAR_FORM = '(1x,f6.3,4x)'
          endif
        else if(STEP_LAST_TIME(4))then
            char  = 'RMS Hght '
            char2 = '  (cm)   '
            INT_VAR = 4
            CHAR_FORM = '(1x,f6.3,4x)'
        else if(STEP_LAST_TIME(5))then
            char =  'Corr Lng '
            char2 = '  (cm)   '
            INT_VAR = 5
            CHAR_FORM = '(1x,f6.2,4x)'
        else if(STEP_LAST_TIME(6))then
            char =  'Sand Frac'
            char2 = '         '
            INT_VAR = 6
            CHAR_FORM = '(1x,f6.2,4x)'
        else if(STEP_LAST_TIME(7))then
            char =  'Clay Frac'
            char2 = '         '
            INT_VAR = 7
            CHAR_FORM = '(1x,f6.2,4x)'
        else if(STEP_LAST_TIME(8))then
            char =  '  Slt   '
            char2 = '  (ppt)  '
            INT_VAR = 8
            CHAR_FORM = '(1x,f6.3,4x)'
        else if(STEP_LAST_TIME(9))then
          if(LOG_EPS_TABLE(1))then
            char =  '  Trunk  '
            char2 = ' RE(epsr)'
            INT_VAR = 9
            CHAR_FORM = '(1x,f6.3,4x)'
          else
            char =  'Trunk Mg '
            char2 = '         '
            INT_VAR = 9
            CHAR_FORM = '(1x,f6.3,4x)'
          endif
        else if(STEP_LAST_TIME(10))then
            char =  ' Density '
            char2 = 'trees/m*m'
            INT_VAR = 10
            CHAR_FORM = '(1x,f6.4,4x)'
        else if(STEP_LAST_TIME(11))then
            char =  'Trunk Rho'
            char2 = '         '
            INT_VAR = 11
            CHAR_FORM = '(1x,f6.4,4x)'
        else if(STEP_LAST_TIME(12))then
            char =  'Crown Hgt'
            char2 = ' (meters)'
            INT_VAR = 12
            CHAR_FORM = '(1x,f6.3,4x)'
        else if(STEP_LAST_TIME(13))then
            char =  'Trnk Diam'
            char2 = '   (cm)  '
            INT_VAR = 13
            CHAR_FORM = '(1x,f6.3,4x)'
        else if(STEP_LAST_TIME(14))then
            char =  'Trunk Hgt'
            char2 = ' (meters)'
            INT_VAR = 14
            CHAR_FORM = '(1x,f6.3,4x)'
        else if(STEP_LAST_TIME(15))then
          if(LOG_EPS_TABLE(8))then
            char =  '   Leaf  '
            char2 = ' RE(epsr)'
            INT_VAR = 15
            CHAR_FORM = '(1x,f6.3,4x)'
          else
            char =  ' Leaf Mg '
            char2 = '         '
            INT_VAR = 15
            CHAR_FORM = '(1x,f6.3,4x)'
          endif
        else if(STEP_LAST_TIME(16))then
            char =  ' Leaf Rho'
            char2 = '         '
            INT_VAR = 16
            CHAR_FORM = '(1x,f6.4,4x)'
        else if(STEP_LAST_TIME(17))then
            char =  'Leaf Dens'
            char2 = 'leaf/m**3'
            INT_VAR = 17
            CHAR_FORM = '(1x,f6.1,4x)'
        else if(STEP_LAST_TIME(18))then
            char =  'Leaf Diam'
            char2 = '  (cm)   '
            INT_VAR = 18
            CHAR_FORM = '(1x,f6.3,4x)'
        else if(STEP_LAST_TIME(19))then
            char =  'Leaf Thck'
            char2 = '  (cm)   '
            INT_VAR = 19
            CHAR_FORM = '(1x,f6.4,4x)'
        else if(STEP_LAST_TIME(20))then
          if(LOG_EPS_TABLE(9))then
            char =  ' Needle  '
            char2 = ' RE(epsr)'
            INT_VAR = 20
            CHAR_FORM = '(1x,f6.3,4x)'
          else
            char =  'Needle Mg'
            char2 = '         '
            INT_VAR = 20
            CHAR_FORM = '(1x,f6.3,4x)'
          endif
        else if(STEP_LAST_TIME(21))then
            char =  ' Ndle Rho'
            char2 = '         '
            INT_VAR = 21
            CHAR_FORM = '(1x,f6.4,4x)'
        else if(STEP_LAST_TIME(22))then
            char =  'Ndle Dens'
            char2 = 'ndle/m**3'
            INT_VAR = 22
            CHAR_FORM = '(E11.4)'
        else if(STEP_LAST_TIME(23))then
            char =  'Ndle Diam'
            char2 = '  (cm)   '
            INT_VAR = 23
            CHAR_FORM = '(1x,f6.3,4x)'
        else if(STEP_LAST_TIME(24))then
            char =  'Ndle Lng '
            char2 = '  (cm)   '
            INT_VAR = 24
            CHAR_FORM = '(1x,f6.3,4x)'
        else if(STEP_LAST_TIME(25))then
          if(LOG_EPS_TABLE(2))then
            char =  ' Branch  '
            char2 = ' RE(epsr)'
            INT_VAR = 25
            CHAR_FORM = '(1x,f6.3,4x)'
          else
            char =  'Branch Mg'
            char2 = '         '
            INT_VAR = 25
            CHAR_FORM = '(1x,f6.3,4x)'
          endif
        else if(STEP_LAST_TIME(26))then
            char =  'Brnch Rho'
            char2 = '         '
            INT_VAR = 26
            CHAR_FORM = '(1x,f6.4,4x)'
        else if(STEP_LAST_TIME(27))then
            char =  'Brnch Den'
            char2 = ' brn/m**3'
            INT_VAR = 27
            CHAR_FORM = '(1x,f6.2,4x)'
        else if(STEP_LAST_TIME(28))then
            char =  'Brnch Dia'
            char2 = '  (cm)   '
            INT_VAR = 28
            CHAR_FORM = '(1x,f6.3,4x)'
        else if(STEP_LAST_TIME(29))then
            char =  'Brnch Lng'
            char2 = ' (meter) '
            INT_VAR = 29
            CHAR_FORM = '(1x,f6.3,4x)'
        else if(STEP_LAST_TIME(30))then
            char =  'Soil Temp'
            char2 = '  Deg C. '
            INT_VAR = 30
            CHAR_FORM = '(1x,f6.2,4x)'
        else if(STEP_LAST_TIME(31))then
            char =  'Veg. Temp'
            char2 = '  Deg C. '
            INT_VAR = 31
            CHAR_FORM = '(1x,f6.2,4x)'
        else if(STEP_LAST_TIME(32))then
            char =  'Water Tmp'
            char2 = '  Deg C. '
            INT_VAR = 32
            CHAR_FORM = '(1x,f6.2,4x)'
C
        else if(STEP_LAST_TIME(33))then
          if(LOG_EPS_TABLE(3))then
            char =  ' Brnch2  '
            char2 = ' RE(epsr)'
            INT_VAR = 33
            CHAR_FORM = '(1x,f6.3,4x)'
          else
            char =  'Brnch2 Mg'
            char2 = '         '
            INT_VAR = 33
            CHAR_FORM = '(1x,f6.3,4x)'
          endif
        else if(STEP_LAST_TIME(34))then
            char =  'Brch2 Rho'
            char2 = '         '
            INT_VAR = 34
            CHAR_FORM = '(1x,f6.4,4x)'
        else if(STEP_LAST_TIME(35))then
            char =  'Brch2 Den'
            char2 = ' brn/m**3'
            INT_VAR = 35
            CHAR_FORM = '(1x,f6.2,4x)'
        else if(STEP_LAST_TIME(36))then
            char =  'Brch2 Dia'
            char2 = '  (cm)   '
            INT_VAR = 36
            CHAR_FORM = '(1x,f6.3,4x)'
        else if(STEP_LAST_TIME(37))then
            char =  'Brch2 Lng'
            char2 = ' (meter) '
            INT_VAR = 37
            CHAR_FORM = '(1x,f6.3,4x)'
C
        else if(STEP_LAST_TIME(38))then
          if(LOG_EPS_TABLE(4))then
            char =  ' Brnch3  '
            char2 = ' RE(epsr)'
            INT_VAR = 38
            CHAR_FORM = '(1x,f6.3,4x)'
          else
            char =  'Brnch3 Mg'
            char2 = '         '
            INT_VAR = 38
            CHAR_FORM = '(1x,f6.3,4x)'
          endif
        else if(STEP_LAST_TIME(39))then
            char =  'Brch3 Rho'
            char2 = '         '
            INT_VAR = 39
            CHAR_FORM = '(1x,f6.4,4x)'
        else if(STEP_LAST_TIME(40))then
            char =  'Brch3 Den'
            char2 = ' brn/m**3'
            INT_VAR = 40
            CHAR_FORM = '(1x,f6.2,4x)'
        else if(STEP_LAST_TIME(41))then
            char =  'Brch3 Dia'
            char2 = '  (cm)   '
            INT_VAR = 41
            CHAR_FORM = '(1x,f6.3,4x)'
        else if(STEP_LAST_TIME(42))then
            char =  'Brch3 Lng'
            char2 = ' (meter) '
            INT_VAR = 42
            CHAR_FORM = '(1x,f6.3,4x)'
C
        else if(STEP_LAST_TIME(43))then
          if(LOG_EPS_TABLE(5))then
            char =  ' Brnch4  '
            char2 = ' RE(epsr)'
            INT_VAR = 43
            CHAR_FORM = '(1x,f6.3,4x)'
          else
            char =  'Brnch4 Mg'
            char2 = '         '
            INT_VAR = 43
            CHAR_FORM = '(1x,f6.3,4x)'
          endif
        else if(STEP_LAST_TIME(44))then
            char =  'Brch4 Rho'
            char2 = '         '
            INT_VAR = 44
            CHAR_FORM = '(1x,f6.4,4x)'
        else if(STEP_LAST_TIME(45))then
            char =  'Brch4 Den'
            char2 = ' brn/m**3'
            INT_VAR = 45
            CHAR_FORM = '(1x,f6.2,4x)'
        else if(STEP_LAST_TIME(46))then
            char =  'Brch4 Dia'
            char2 = '  (cm)   '
            INT_VAR = 46
            CHAR_FORM = '(1x,f6.3,4x)'
        else if(STEP_LAST_TIME(47))then
            char =  'Brch4 Lng'
            char2 = ' (meter) '
            INT_VAR = 47
            CHAR_FORM = '(1x,f6.3,4x)'
C
        else if(STEP_LAST_TIME(48))then
          if(LOG_EPS_TABLE(6))then
            char =  ' Brnch5  '
            char2 = ' RE(epsr)'
            INT_VAR = 48
            CHAR_FORM = '(1x,f6.3,4x)'
          else
            char =  'Brnch5 Mg'
            char2 = '         '
            INT_VAR = 48
            CHAR_FORM = '(1x,f6.3,4x)'
          endif
        else if(STEP_LAST_TIME(49))then
            char =  'Brch5 Rho'
            char2 = '         '
            INT_VAR = 49
            CHAR_FORM = '(1x,f6.4,4x)'
        else if(STEP_LAST_TIME(50))then
            char =  'Brch5 Den'
            char2 = ' brn/m**3'
            INT_VAR = 50
            CHAR_FORM = '(1x,f6.2,4x)'
        else if(STEP_LAST_TIME(51))then
            char =  'Brch5 Dia'
            char2 = '  (cm)   '
            INT_VAR = 51
            CHAR_FORM = '(1x,f6.3,4x)'
        else if(STEP_LAST_TIME(52))then
            char =  'Brch5 Lng'
            char2 = ' (meter) '
            INT_VAR = 52
            CHAR_FORM = '(1x,f6.3,4x)'
C
        else if(STEP_LAST_TIME(53))then
          if(LOG_EPS_TABLE(7))then
            char =  ' Brnch6  '
            char2 = ' RE(epsr)'
            INT_VAR = 53
            CHAR_FORM = '(1x,f6.3,4x)'
          else
            char =  'Brnch6 Mg'
            char2 = '         '
            INT_VAR = 53
            CHAR_FORM = '(1x,f6.3,4x)'
          endif
        else if(STEP_LAST_TIME(54))then
            char =  'Brch6 Rho'
            char2 = '         '
            INT_VAR = 54
            CHAR_FORM = '(1x,f6.4,4x)'
        else if(STEP_LAST_TIME(55))then
            char =  'Brch6 Den'
            char2 = ' brn/m**3'
            INT_VAR = 55
            CHAR_FORM = '(1x,f6.2,4x)'
        else if(STEP_LAST_TIME(56))then
            char =  'Brch6 Dia'
            char2 = '  (cm)   '
            INT_VAR = 56
            CHAR_FORM = '(1x,f6.3,4x)'
        else if(STEP_LAST_TIME(57))then
            char =  'Brch6 Lng'
            char2 = ' (meter) '
            INT_VAR = 57
            CHAR_FORM = '(1x,f6.3,4x)'

C
        else if(STEP_LAST_TIME(60))then
            char =  'Snow Dpth'
            char2 = ' (meter) '
            INT_VAR = 60
            CHAR_FORM = '(1x,f6.3,4x)'
        else
            char  = ' Theta '
            char2 = ' (deg) '
            INT_VAR = 1
            CHAR_FORM = '(1x,f6.3,4x)'
        endif
c
            numfamily = numfamily + 1
C
            WRITE(8,140) NUMFAMILY
            WRITE(7,140) NUMFAMILY
            WRITE(11,140) NUMFAMILY
            WRITE(12,140) NUMFAMILY
            WRITE(13,140) NUMFAMILY
            WRITE(14,140) NUMFAMILY
            if(l_pol) WRITE(15,140) NUMFAMILY
            if(l_phase)WRITE(16,140) NUMFAMILY
            if(l_kappa) WRITE(17,140) NUMFAMILY
            if(l_trans)WRITE(18,140) NUMFAMILY
            WRITE(19,140) NUMFAMILY


            write(8,150) THETA_DEGREES, FREQ_GHZ
            write(8,151) DENSITY, CROWN_HGHT, TRUNK_HGHT, T_VEG
            write(7,150) THETA_DEGREES, FREQ_GHZ
            write(7,151) DENSITY, CROWN_HGHT, TRUNK_HGHT, T_VEG
            write(11,150) THETA_DEGREES, FREQ_GHZ
            write(11,151) DENSITY, CROWN_HGHT, TRUNK_HGHT, T_VEG
            write(12,150) THETA_DEGREES, FREQ_GHZ
            write(12,151) DENSITY, CROWN_HGHT, TRUNK_HGHT, T_VEG
            write(13,150) THETA_DEGREES, FREQ_GHZ
            write(13,151) DENSITY, CROWN_HGHT, TRUNK_HGHT, T_VEG
            write(14,150) THETA_DEGREES, FREQ_GHZ
            write(14,151) DENSITY, CROWN_HGHT, TRUNK_HGHT, T_VEG
            if(l_pol) then
                write(15,150) THETA_DEGREES, FREQ_GHZ
                write(15,151) DENSITY, CROWN_HGHT, TRUNK_HGHT, T_VEG
            endif
            if(l_phase) then
                write(16,150) THETA_DEGREES, FREQ_GHZ
                write(16,151) DENSITY, CROWN_HGHT, TRUNK_HGHT, T_VEG
            endif
            if(l_kappa) then
                write(17,150) THETA_DEGREES, FREQ_GHZ
                write(17,151) DENSITY, CROWN_HGHT, TRUNK_HGHT, T_VEG
            endif
            if(l_trans) then
                write(18,150) THETA_DEGREES, FREQ_GHZ
                write(18,151) DENSITY, CROWN_HGHT, TRUNK_HGHT, T_VEG
            endif
            write(19,150) THETA_DEGREES, FREQ_GHZ
            write(19,151) DENSITY, CROWN_HGHT, TRUNK_HGHT, T_VEG




c
        if(LOG_EPS_TABLE(10).and.SOIL_SURFACE)then
          write(7,160) real(epsilonrc(2)),aimag(epsilonrc(2)),RMS_SOIL,
     &                   LS_SOIL,c_surf
          write(8,160) real(epsilonrc(2)),aimag(epsilonrc(2)),RMS_SOIL,
     &                   LS_SOIL,c_surf
          write(11,160) real(epsilonrc(2)),aimag(epsilonrc(2)),RMS_SOIL,
     &                   LS_SOIL,c_surf
          write(12,160) real(epsilonrc(2)),aimag(epsilonrc(2)),RMS_SOIL,
     &                   LS_SOIL,c_surf
          write(13,160) real(epsilonrc(2)),aimag(epsilonrc(2)),RMS_SOIL,
     &                   LS_SOIL,c_surf
          write(14,160) real(epsilonrc(2)),aimag(epsilonrc(2)),RMS_SOIL,
     &                   LS_SOIL,c_surf
          if(l_pol)write(15,160) real(epsilonrc(2)),aimag(epsilonrc(2)),
     &                   RMS_SOIL, LS_SOIL,c_surf
          if(l_phase) write(16,160) real(epsilonrc(2)),
     &                   aimag(epsilonrc(2)),RMS_SOIL,LS_SOIL,c_surf
          if(l_kappa) write(17,160) real(epsilonrc(2)),
     &                   aimag(epsilonrc(2)),RMS_SOIL,LS_SOIL,c_surf
          if(l_trans) write(18,160) real(epsilonrc(2)),
     &                   aimag(epsilonrc(2)),RMS_SOIL,LS_SOIL,c_surf
          write(19,160) real(epsilonrc(2)),aimag(epsilonrc(2)),RMS_SOIL,
     &                   LS_SOIL,c_surf
         else if(LOG_EPS_TABLE(10).and.WATER_SURFACE)then
            write(7,161) real(epsilonr(2)),aimag(epsilonr(2))
            write(8,161) real(epsilonr(2)),aimag(epsilonr(2))
            write(11,161) real(epsilonr(2)),aimag(epsilonr(2))
            write(12,161) real(epsilonr(2)),aimag(epsilonr(2))
            write(13,161) real(epsilonr(2)),aimag(epsilonr(2))
            write(14,161) real(epsilonr(2)),aimag(epsilonr(2))
            if(l_pol) write(15,161) real(epsilonr(2)),aimag(epsilonr(2))
            if(l_phase)write(16,161)real(epsilonr(2)),aimag(epsilonr(2))
            if(l_kappa)write(17,161)real(epsilonr(2)),aimag(epsilonr(2))
            if(l_trans)write(18,161)real(epsilonr(2)),aimag(epsilonr(2))
            write(19,161)real(epsilonr(2)),aimag(epsilonr(2))
         else if(SOIL_SURFACE)then
         write(8,152) MV_SOIL,RMS_SOIL,LS_SOIL,SAND,CLAY,T_SOIL,c_surf
         write(7,152) MV_SOIL,RMS_SOIL,LS_SOIL,SAND,CLAY,T_SOIL,c_surf
         write(11,152) MV_SOIL,RMS_SOIL,LS_SOIL,SAND,CLAY,T_SOIL,c_surf
         write(12,152) MV_SOIL,RMS_SOIL,LS_SOIL,SAND,CLAY,T_SOIL,c_surf
         write(13,152) MV_SOIL,RMS_SOIL,LS_SOIL,SAND,CLAY,T_SOIL,c_surf
         write(14,152) MV_SOIL,RMS_SOIL,LS_SOIL,SAND,CLAY,T_SOIL,c_surf
         if(l_pol) write(15,152) MV_SOIL,RMS_SOIL,LS_SOIL,SAND,CLAY,
     &                           T_SOIL,c_surf
         if(l_phase) write(16,152) MV_SOIL,RMS_SOIL,LS_SOIL,SAND,CLAY,
     &                           T_SOIL,c_surf
         if(l_kappa) write(17,152) MV_SOIL,RMS_SOIL,LS_SOIL,SAND,CLAY,
     &                           T_SOIL,c_surf
         if(l_trans) write(18,152) MV_SOIL,RMS_SOIL,LS_SOIL,SAND,CLAY,
     &                           T_SOIL,c_surf
         write(19,152) MV_SOIL,RMS_SOIL,LS_SOIL,SAND,CLAY,T_SOIL,c_surf
         else if(WATER_SURFACE) then
            write(8,153) SALT, T_WATER
            write(7,153) SALT, T_WATER
            write(11,153) SALT, T_WATER
            write(12,153) SALT, T_WATER
            write(13,153) SALT, T_WATER
            write(14,153) SALT, T_WATER
            if(l_pol) write(15,153) SALT, T_WATER
            if(l_phase) write(16,153) SALT, T_WATER
            if(l_kappa) write(17,153) SALT, T_WATER
            if(l_trans) write(18,153) SALT, T_WATER
            write(19,153) SALT, T_WATER
        endif
c
        if(snow_surface)then
            write(8,159) REAL(EPSILONRC(8)),AIMAG(EPSILONRC(8)),T_SNOW
            write(7,159) REAL(EPSILONRC(8)),AIMAG(EPSILONRC(8)),T_SNOW
            write(11,159) REAL(EPSILONRC(8)),AIMAG(EPSILONRC(8)),T_SNOW
            write(12,159) REAL(EPSILONRC(8)),AIMAG(EPSILONRC(8)),T_SNOW
            write(13,159) REAL(EPSILONRC(8)),AIMAG(EPSILONRC(8)),T_SNOW
            write(14,159) REAL(EPSILONRC(8)),AIMAG(EPSILONRC(8)),T_SNOW
            if(l_pol) write(15,159) REAL(EPSILONRC(8)),
     &                  AIMAG(EPSILONRC(8)),T_SNOW
            if(l_phase) write(16,159) REAL(EPSILONRC(8)),
     &                  AIMAG(EPSILONRC(8)),T_SNOW
            if(l_kappa) write(17,159) REAL(EPSILONRC(8)),
     &                  AIMAG(EPSILONRC(8)),T_SNOW
            if(l_trans) write(18,159) REAL(EPSILONRC(8)),
     &                  AIMAG(EPSILONRC(8)),T_SNOW
            write(19,159) REAL(EPSILONRC(8)),AIMAG(EPSILONRC(8)),T_SNOW
        endif
c
        if(LOG_CONSTITUENT(TRUNK))then
          if(LOG_EPS_TABLE(TRUNK))then
            write(8,164) REAL(EPSILONRC(5)),AIMAG(EPSILONRC(5)),
     &                    TRUNK_DIAM
            write(7,164) REAL(EPSILONRC(5)),AIMAG(EPSILONRC(5)),
     &                    TRUNK_DIAM
            write(11,164) REAL(EPSILONRC(5)),AIMAG(EPSILONRC(5)),
     &                    TRUNK_DIAM
            write(12,164) REAL(EPSILONRC(5)),AIMAG(EPSILONRC(5)),
     &                    TRUNK_DIAM
            write(13,164) REAL(EPSILONRC(5)),AIMAG(EPSILONRC(5)),
     &                    TRUNK_DIAM
            write(14,164) REAL(EPSILONRC(5)),AIMAG(EPSILONRC(5)),
     &                    TRUNK_DIAM
            if(l_pol) write(15,164) REAL(EPSILONRC(5)),
     &                    AIMAG(EPSILONRC(5)), TRUNK_DIAM
            if(l_phase) write(16,164) REAL(EPSILONRC(5)),
     &                    AIMAG(EPSILONRC(5)), TRUNK_DIAM
            if(l_kappa) write(17,164) REAL(EPSILONRC(5)),
     &                    AIMAG(EPSILONRC(5)), TRUNK_DIAM
            if(l_trans) write(18,164) REAL(EPSILONRC(5)),
     &                    AIMAG(EPSILONRC(5)), TRUNK_DIAM
            write(19,164) REAL(EPSILONRC(5)),AIMAG(EPSILONRC(5)),
     &                    TRUNK_DIAM
          else
            write(8,154) MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
            write(7,154) MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
            write(11,154) MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
            write(12,154) MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
            write(13,154) MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
            write(14,154) MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
            if(l_pol) write(15,154) MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
            if(l_phase) write(16,154) MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
            if(l_kappa) write(17,154) MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
            if(l_trans) write(18,154) MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
            write(19,154) MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
           endif
        endif
c
        if(LOG_CONSTITUENT(LEAVES)) then
          if(LOG_EPS_TABLE(LEAVES))then
            write(8,165) REAL(EPSILONRC(4)),AIMAG(EPSILONRC(4)),
     &                 LEAF_DENS,LEAF_DIAM,LEAF_TAU
            write(7,165) REAL(EPSILONRC(4)),AIMAG(EPSILONRC(4)),
     &                 LEAF_DENS,LEAF_DIAM,LEAF_TAU
            write(11,165)REAL(EPSILONRC(4)),AIMAG(EPSILONRC(4)),
     &                 LEAF_DENS,LEAF_DIAM,LEAF_TAU
            write(12,165)REAL(EPSILONRC(4)),AIMAG(EPSILONRC(4)),
     &                 LEAF_DENS,LEAF_DIAM,LEAF_TAU
            write(13,165)REAL(EPSILONRC(4)),AIMAG(EPSILONRC(4)),
     &                 LEAF_DENS,LEAF_DIAM,LEAF_TAU
            write(14,165)REAL(EPSILONRC(4)),AIMAG(EPSILONRC(4)),
     &                 LEAF_DENS,LEAF_DIAM,LEAF_TAU
            if(l_pol) write(15,165)REAL(EPSILONRC(4)),
     &                 AIMAG(EPSILONRC(4)),LEAF_DENS,LEAF_DIAM,LEAF_TAU
            if(l_phase) write(16,165)REAL(EPSILONRC(4)),
     &                 AIMAG(EPSILONRC(4)),LEAF_DENS,LEAF_DIAM,LEAF_TAU
            if(l_kappa) write(17,165)REAL(EPSILONRC(4)),
     &                 AIMAG(EPSILONRC(4)),LEAF_DENS,LEAF_DIAM,LEAF_TAU
            if(l_trans) write(18,165)REAL(EPSILONRC(4)),
     &                 AIMAG(EPSILONRC(4)),LEAF_DENS,LEAF_DIAM,LEAF_TAU
            write(19,165)REAL(EPSILONRC(4)),AIMAG(EPSILONRC(4)),
     &                 LEAF_DENS,LEAF_DIAM,LEAF_TAU
          else
            write(8,155) MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
            write(7,155) MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
            write(11,155) MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
            write(12,155) MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
            write(13,155) MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
            write(14,155) MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
            if(l_pol) write(15,155) MG_LEAF,RHO_LEAF,
     &                        LEAF_DENS,LEAF_DIAM,LEAF_TAU
             if(l_phase) write(16,155) MG_LEAF,RHO_LEAF,LEAF_DENS,
     &                       LEAF_DIAM,LEAF_TAU
            if(l_kappa) write(17,155) MG_LEAF,RHO_LEAF,
     &                        LEAF_DENS,LEAF_DIAM,LEAF_TAU
            if(l_trans) write(18,155) MG_LEAF,RHO_LEAF,LEAF_DENS,
     &                       LEAF_DIAM,LEAF_TAU
            write(19,155) MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
          endif
        endif
c
        if(LOG_CONSTITUENT(NEEDLE)) then
          if(LOG_EPS_TABLE(NEEDLE))then
            write(8,166) REAL(EPSILONRC(3)),AIMAG(EPSILONRC(3)),
     &                  NDL_DENS, NDL_DIAM, NDL_LNG
            write(7,166) REAL(EPSILONRC(3)),AIMAG(EPSILONRC(3)),
     &                  NDL_DENS, NDL_DIAM, NDL_LNG
            write(11,166)REAL(EPSILONRC(3)),AIMAG(EPSILONRC(3)),
     &                  NDL_DENS, NDL_DIAM, NDL_LNG
            write(12,166)REAL(EPSILONRC(3)),AIMAG(EPSILONRC(3)),
     &                  NDL_DENS, NDL_DIAM, NDL_LNG
            write(13,166)REAL(EPSILONRC(3)),AIMAG(EPSILONRC(3)),
     &                  NDL_DENS, NDL_DIAM, NDL_LNG
            write(14,166)REAL(EPSILONRC(3)),AIMAG(EPSILONRC(3)),
     &                  NDL_DENS, NDL_DIAM, NDL_LNG
            if(l_pol) write(15,166)REAL(EPSILONRC(3)),
     &                AIMAG(EPSILONRC(3)), NDL_DENS, NDL_DIAM, NDL_LNG
            if(l_phase) write(16,166)REAL(EPSILONRC(3)),
     &                AIMAG(EPSILONRC(3)),NDL_DENS, NDL_DIAM, NDL_LNG
            if(l_kappa) write(17,166)REAL(EPSILONRC(3)),
     &                AIMAG(EPSILONRC(3)), NDL_DENS, NDL_DIAM, NDL_LNG
            if(l_trans) write(18,166)REAL(EPSILONRC(3)),
     &                AIMAG(EPSILONRC(3)),NDL_DENS, NDL_DIAM, NDL_LNG
            write(19,166)REAL(EPSILONRC(3)),AIMAG(EPSILONRC(3)),
     &                  NDL_DENS, NDL_DIAM, NDL_LNG
          else
            write(8,156) MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
            write(7,156) MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
            write(11,156) MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
            write(12,156) MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
            write(13,156) MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
            write(14,156) MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
            if(l_pol) write(15,156) MG_NDL, RHO_NDL, NDL_DENS, 
     &                             NDL_DIAM, NDL_LNG
            if(l_phase) write(16,156) MG_NDL, RHO_NDL, NDL_DENS,
     &                             NDL_DIAM, NDL_LNG
            if(l_kappa) write(17,156) MG_NDL, RHO_NDL, NDL_DENS, 
     &                             NDL_DIAM, NDL_LNG
            if(l_trans) write(18,156) MG_NDL, RHO_NDL, NDL_DENS,
     &                             NDL_DIAM, NDL_LNG
            write(19,156) MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
          endif
        endif
c
        if(LOG_CONSTITUENT(BR1)) then
          if(LOG_EPS_TABLE(BR1))then
            write(8,167) REAL(EPSILONRC(6)),AIMAG(EPSILONRC(6)),
     &                  BR1_DENS, BR1_DIAM, BR1_LNG
            write(7,167) REAL(EPSILONRC(6)),AIMAG(EPSILONRC(6)),
     &                  BR1_DENS, BR1_DIAM, BR1_LNG
            write(11,167)REAL(EPSILONRC(6)),AIMAG(EPSILONRC(6)),
     &                  BR1_DENS, BR1_DIAM, BR1_LNG
            write(12,167)REAL(EPSILONRC(6)),AIMAG(EPSILONRC(6)),
     &                  BR1_DENS, BR1_DIAM, BR1_LNG
            write(13,167)REAL(EPSILONRC(6)),AIMAG(EPSILONRC(6)),
     &                  BR1_DENS, BR1_DIAM, BR1_LNG
            write(14,167)REAL(EPSILONRC(6)),AIMAG(EPSILONRC(6)),
     &                  BR1_DENS, BR1_DIAM, BR1_LNG
            if(l_pol) write(15,167)REAL(EPSILONRC(6)),
     &              AIMAG(EPSILONRC(6)),BR1_DENS, BR1_DIAM, BR1_LNG
            if(l_phase) write(16,167)REAL(EPSILONRC(6)),
     &              AIMAG(EPSILONRC(6)),BR1_DENS, BR1_DIAM, BR1_LNG
            if(l_kappa) write(17,167)REAL(EPSILONRC(6)),
     &              AIMAG(EPSILONRC(6)),BR1_DENS, BR1_DIAM, BR1_LNG
            if(l_trans) write(18,167)REAL(EPSILONRC(6)),
     &              AIMAG(EPSILONRC(6)),BR1_DENS, BR1_DIAM, BR1_LNG
            write(19,167)REAL(EPSILONRC(6)),AIMAG(EPSILONRC(6)),
     &                  BR1_DENS, BR1_DIAM, BR1_LNG
          else
            write(8,157) MG_BR1, RHO_BR1, BR1_DENS, BR1_DIAM, BR1_LNG
            write(7,157) MG_BR1, RHO_BR1, BR1_DENS, BR1_DIAM, BR1_LNG
            write(11,157) MG_BR1, RHO_BR1, BR1_DENS, BR1_DIAM, BR1_LNG
            write(12,157) MG_BR1, RHO_BR1, BR1_DENS, BR1_DIAM, BR1_LNG
            write(13,157) MG_BR1, RHO_BR1, BR1_DENS, BR1_DIAM, BR1_LNG
            write(14,157) MG_BR1, RHO_BR1, BR1_DENS, BR1_DIAM, BR1_LNG
            if(l_pol) write(15,157) MG_BR1, RHO_BR1, BR1_DENS,
     &                             BR1_DIAM, BR1_LNG
            if(l_phase) write(16,157) MG_BR1, RHO_BR1, BR1_DENS,
     &                             BR1_DIAM, BR1_LNG
            if(l_kappa) write(17,157) MG_BR1, RHO_BR1, BR1_DENS,
     &                             BR1_DIAM, BR1_LNG
            if(l_trans) write(18,157) MG_BR1, RHO_BR1, BR1_DENS,
     &                             BR1_DIAM, BR1_LNG
            write(19,157) MG_BR1, RHO_BR1, BR1_DENS, BR1_DIAM, BR1_LNG
          endif
        endif
c
        if(LOG_CONSTITUENT(BR2)) then
          if(LOG_EPS_TABLE(BR2))then
            write(8,168) REAL(EPSILONRC(7)),AIMAG(EPSILONRC(7)),
     &                  BR2_DENS, BR2_DIAM, BR2_LNG
            write(7,168) REAL(EPSILONRC(7)),AIMAG(EPSILONRC(7)),
     &                  BR2_DENS, BR2_DIAM, BR2_LNG
            write(11,168)REAL(EPSILONRC(7)),AIMAG(EPSILONRC(7)),
     &                  BR2_DENS, BR2_DIAM, BR2_LNG
            write(12,168)REAL(EPSILONRC(7)),AIMAG(EPSILONRC(7)),
     &                  BR2_DENS, BR2_DIAM, BR2_LNG
            write(13,168)REAL(EPSILONRC(7)),AIMAG(EPSILONRC(7)),
     &                  BR2_DENS, BR2_DIAM, BR2_LNG
            write(14,168)REAL(EPSILONRC(7)),AIMAG(EPSILONRC(7)),
     &                  BR2_DENS, BR2_DIAM, BR2_LNG
            if(l_pol) write(15,168)REAL(EPSILONRC(7)),
     &                AIMAG(EPSILONRC(7)), BR2_DENS, BR2_DIAM, BR2_LNG
            if(l_phase) write(16,168)REAL(EPSILONRC(7)),
     &                AIMAG(EPSILONRC(7)), BR2_DENS, BR2_DIAM, BR2_LNG
            if(l_kappa) write(17,168)REAL(EPSILONRC(7)),
     &                AIMAG(EPSILONRC(7)), BR2_DENS, BR2_DIAM, BR2_LNG
            if(l_trans) write(18,168)REAL(EPSILONRC(7)),
     &                AIMAG(EPSILONRC(7)), BR2_DENS, BR2_DIAM, BR2_LNG
            write(19,168)REAL(EPSILONRC(7)),AIMAG(EPSILONRC(7)),
     &                  BR2_DENS, BR2_DIAM, BR2_LNG

          else
            write(8,158) MG_BR2, RHO_BR2, BR2_DENS, BR2_DIAM, BR2_LNG
            write(7,158) MG_BR2, RHO_BR2, BR2_DENS, BR2_DIAM, BR2_LNG
            write(11,158) MG_BR2, RHO_BR2, BR2_DENS, BR2_DIAM, BR2_LNG
            write(12,158) MG_BR2, RHO_BR2, BR2_DENS, BR2_DIAM, BR2_LNG
            write(13,158) MG_BR2, RHO_BR2, BR2_DENS, BR2_DIAM, BR2_LNG
            write(14,158) MG_BR2, RHO_BR2, BR2_DENS, BR2_DIAM, BR2_LNG
            if(l_pol) write(15,158) MG_BR2, RHO_BR2, BR2_DENS, 
     &                             BR2_DIAM, BR2_LNG
            if(l_phase) write(16,158) MG_BR2, RHO_BR2, BR2_DENS,
     &                             BR2_DIAM, BR2_LNG
            if(l_kappa) write(17,158) MG_BR2, RHO_BR2, BR2_DENS, 
     &                             BR2_DIAM, BR2_LNG
            if(l_trans) write(18,158) MG_BR2, RHO_BR2, BR2_DENS,
     &                             BR2_DIAM, BR2_LNG
            write(19,158) MG_BR2, RHO_BR2, BR2_DENS, BR2_DIAM, BR2_LNG
           endif
        endif
c
        if(LOG_CONSTITUENT(BR3)) then
          if(LOG_EPS_TABLE(BR3))then
            write(8,312) REAL(EPSILONRC(9)),AIMAG(EPSILONRC(9)),
     &                  BR3_DENS, BR3_DIAM, BR3_LNG
            write(7,312) REAL(EPSILONRC(9)),AIMAG(EPSILONRC(9)),
     &                  BR3_DENS, BR3_DIAM, BR3_LNG
            write(11,312)REAL(EPSILONRC(9)),AIMAG(EPSILONRC(9)),
     &                  BR3_DENS, BR3_DIAM, BR3_LNG
            write(12,312)REAL(EPSILONRC(9)),AIMAG(EPSILONRC(9)),
     &                  BR3_DENS, BR3_DIAM, BR3_LNG
            write(13,312)REAL(EPSILONRC(9)),AIMAG(EPSILONRC(9)),
     &                  BR3_DENS, BR3_DIAM, BR3_LNG
            write(14,312)REAL(EPSILONRC(9)),AIMAG(EPSILONRC(9)),
     &                  BR3_DENS, BR3_DIAM, BR3_LNG
            if(l_pol) write(15,312)REAL(EPSILONRC(9)),
     &                AIMAG(EPSILONRC(9)), BR3_DENS, BR3_DIAM, BR3_LNG
            if(l_phase) write(16,312)REAL(EPSILONRC(9)),
     &                AIMAG(EPSILONRC(9)), BR3_DENS, BR3_DIAM, BR3_LNG
            if(l_kappa) write(17,312)REAL(EPSILONRC(9)),
     &                AIMAG(EPSILONRC(9)), BR3_DENS, BR3_DIAM, BR3_LNG
            if(l_trans) write(18,312)REAL(EPSILONRC(9)),
     &                AIMAG(EPSILONRC(9)), BR3_DENS, BR3_DIAM, BR3_LNG
            write(19,312)REAL(EPSILONRC(9)),
     &                AIMAG(EPSILONRC(9)), BR3_DENS, BR3_DIAM, BR3_LNG
          else
            write(8,302) MG_BR3, RHO_BR3, BR3_DENS, BR3_DIAM, BR3_LNG
            write(7,302) MG_BR3, RHO_BR3, BR3_DENS, BR3_DIAM, BR3_LNG
            write(11,302) MG_BR3, RHO_BR3, BR3_DENS, BR3_DIAM, BR3_LNG
            write(12,302) MG_BR3, RHO_BR3, BR3_DENS, BR3_DIAM, BR3_LNG
            write(13,302) MG_BR3, RHO_BR3, BR3_DENS, BR3_DIAM, BR3_LNG
            write(14,302) MG_BR3, RHO_BR3, BR3_DENS, BR3_DIAM, BR3_LNG
            if(l_pol) write(15,302) MG_BR3, RHO_BR3, BR3_DENS, 
     &                             BR3_DIAM, BR3_LNG
            if(l_phase) write(16,302) MG_BR3, RHO_BR3, BR3_DENS, 
     &                             BR3_DIAM, BR3_LNG
            if(l_kappa) write(17,302) MG_BR3, RHO_BR3, BR3_DENS, 
     &                             BR3_DIAM, BR3_LNG
            if(l_trans) write(18,302) MG_BR3, RHO_BR3, BR3_DENS, 
     &                             BR3_DIAM, BR3_LNG
            write(19,302) MG_BR3, RHO_BR3, BR3_DENS, 
     &                             BR3_DIAM, BR3_LNG
          endif
        endif
c
        if(LOG_CONSTITUENT(BR4)) then
          if(LOG_EPS_TABLE(BR4))then
            write(8,314) REAL(EPSILONRC(10)),AIMAG(EPSILONRC(10)),
     &                  BR4_DENS, BR4_DIAM, BR4_LNG
            write(7,314) REAL(EPSILONRC(10)),AIMAG(EPSILONRC(10)),
     &                  BR4_DENS, BR4_DIAM, BR4_LNG
            write(11,314)REAL(EPSILONRC(10)),AIMAG(EPSILONRC(10)),
     &                  BR4_DENS, BR4_DIAM, BR4_LNG
            write(12,314)REAL(EPSILONRC(10)),AIMAG(EPSILONRC(10)),
     &                  BR4_DENS, BR4_DIAM, BR4_LNG
            write(13,314)REAL(EPSILONRC(10)),AIMAG(EPSILONRC(10)),
     &                  BR4_DENS, BR4_DIAM, BR4_LNG
            write(14,314)REAL(EPSILONRC(10)),AIMAG(EPSILONRC(10)),
     &                  BR4_DENS, BR4_DIAM, BR4_LNG
            if(l_pol) write(15,314)REAL(EPSILONRC(10)),
     &                AIMAG(EPSILONRC(10)), BR4_DENS, BR4_DIAM, BR4_LNG
            if(l_phase) write(16,314)REAL(EPSILONRC(10)),
     &                AIMAG(EPSILONRC(10)), BR4_DENS, BR4_DIAM, BR4_LNG
            if(l_kappa) write(17,314)REAL(EPSILONRC(10)),
     &                AIMAG(EPSILONRC(10)), BR4_DENS, BR4_DIAM, BR4_LNG
            if(l_trans) write(18,314)REAL(EPSILONRC(10)),
     &                AIMAG(EPSILONRC(10)), BR4_DENS, BR4_DIAM, BR4_LNG
            write(19,314)REAL(EPSILONRC(10)),
     &                AIMAG(EPSILONRC(10)), BR4_DENS, BR4_DIAM, BR4_LNG
          else
            write(8,304) MG_BR4, RHO_BR4, BR4_DENS, BR4_DIAM, BR4_LNG
            write(7,304) MG_BR4, RHO_BR4, BR4_DENS, BR4_DIAM, BR4_LNG
            write(11,304) MG_BR4, RHO_BR4, BR4_DENS, BR4_DIAM, BR4_LNG
            write(12,304) MG_BR4, RHO_BR4, BR4_DENS, BR4_DIAM, BR4_LNG
            write(13,304) MG_BR4, RHO_BR4, BR4_DENS, BR4_DIAM, BR4_LNG
            write(14,304) MG_BR4, RHO_BR4, BR4_DENS, BR4_DIAM, BR4_LNG
            if(l_pol) write(15,304) MG_BR4, RHO_BR4, BR4_DENS, 
     &                             BR4_DIAM, BR4_LNG
            if(l_phase) write(16,304) MG_BR4, RHO_BR4, BR4_DENS, 
     &                             BR4_DIAM, BR4_LNG
            if(l_kappa) write(17,304) MG_BR4, RHO_BR4, BR4_DENS, 
     &                             BR4_DIAM, BR4_LNG
            if(l_trans) write(18,304) MG_BR4, RHO_BR4, BR4_DENS, 
     &                             BR4_DIAM, BR4_LNG
            write(19,304) MG_BR4, RHO_BR4, BR4_DENS, 
     &                             BR4_DIAM, BR4_LNG
          endif
        endif
c
        if(LOG_CONSTITUENT(BR5)) then
          if(LOG_EPS_TABLE(BR5))then
            write(8,316) REAL(EPSILONRC(11)),AIMAG(EPSILONRC(11)),
     &                  BR5_DENS, BR5_DIAM, BR5_LNG
            write(7,316) REAL(EPSILONRC(11)),AIMAG(EPSILONRC(11)),
     &                  BR5_DENS, BR5_DIAM, BR5_LNG
            write(11,316)REAL(EPSILONRC(11)),AIMAG(EPSILONRC(11)),
     &                  BR5_DENS, BR5_DIAM, BR5_LNG
            write(12,316)REAL(EPSILONRC(11)),AIMAG(EPSILONRC(11)),
     &                  BR5_DENS, BR5_DIAM, BR5_LNG
            write(13,316)REAL(EPSILONRC(11)),AIMAG(EPSILONRC(11)),
     &                  BR5_DENS, BR5_DIAM, BR5_LNG
            write(14,316)REAL(EPSILONRC(11)),AIMAG(EPSILONRC(11)),
     &                  BR5_DENS, BR5_DIAM, BR5_LNG
            if(l_pol) write(15,316)REAL(EPSILONRC(11)),
     &                AIMAG(EPSILONRC(11)), BR5_DENS, BR5_DIAM, BR5_LNG
            if(l_phase) write(16,316)REAL(EPSILONRC(11)),
     &                AIMAG(EPSILONRC(11)), BR5_DENS, BR5_DIAM, BR5_LNG
            if(l_kappa) write(17,316)REAL(EPSILONRC(11)),
     &                AIMAG(EPSILONRC(11)), BR5_DENS, BR5_DIAM, BR5_LNG
            if(l_trans) write(18,316)REAL(EPSILONRC(11)),
     &                AIMAG(EPSILONRC(11)), BR5_DENS, BR5_DIAM, BR5_LNG
            write(19,316)REAL(EPSILONRC(11)),
     &                AIMAG(EPSILONRC(11)), BR5_DENS, BR5_DIAM, BR5_LNG
          else
            write(8,306) MG_BR5, RHO_BR5, BR5_DENS, BR5_DIAM, BR5_LNG
            write(7,306) MG_BR5, RHO_BR5, BR5_DENS, BR5_DIAM, BR5_LNG
            write(11,306) MG_BR5, RHO_BR5, BR5_DENS, BR5_DIAM, BR5_LNG
            write(12,306) MG_BR5, RHO_BR5, BR5_DENS, BR5_DIAM, BR5_LNG
            write(13,306) MG_BR5, RHO_BR5, BR5_DENS, BR5_DIAM, BR5_LNG
            write(14,306) MG_BR5, RHO_BR5, BR5_DENS, BR5_DIAM, BR5_LNG
            if(l_pol) write(15,306) MG_BR5, RHO_BR5, BR5_DENS, 
     &                             BR5_DIAM, BR5_LNG
            if(l_phase) write(16,306) MG_BR5, RHO_BR5, BR5_DENS, 
     &                             BR5_DIAM, BR5_LNG
            if(l_kappa) write(17,306) MG_BR5, RHO_BR5, BR5_DENS, 
     &                             BR5_DIAM, BR5_LNG
            if(l_trans) write(18,306) MG_BR5, RHO_BR5, BR5_DENS, 
     &                             BR5_DIAM, BR5_LNG
            write(19,306) MG_BR5, RHO_BR5, BR5_DENS, 
     &                             BR5_DIAM, BR5_LNG
          endif
        endif
c
        if(LOG_CONSTITUENT(BR6)) then
          if(LOG_EPS_TABLE(BR6))then
            write(8,316) REAL(EPSILONRC(12)),AIMAG(EPSILONRC(12)),
     &                  BR6_DENS, BR6_DIAM, BR6_LNG
            write(7,316) REAL(EPSILONRC(12)),AIMAG(EPSILONRC(12)),
     &                  BR6_DENS, BR6_DIAM, BR6_LNG
            write(11,316)REAL(EPSILONRC(12)),AIMAG(EPSILONRC(12)),
     &                  BR6_DENS, BR6_DIAM, BR6_LNG
            write(12,316)REAL(EPSILONRC(12)),AIMAG(EPSILONRC(12)),
     &                  BR6_DENS, BR6_DIAM, BR6_LNG
            write(13,316)REAL(EPSILONRC(12)),AIMAG(EPSILONRC(12)),
     &                  BR6_DENS, BR6_DIAM, BR6_LNG
            write(14,316)REAL(EPSILONRC(12)),AIMAG(EPSILONRC(12)),
     &                  BR6_DENS, BR6_DIAM, BR6_LNG
            if(l_pol) write(15,316)REAL(EPSILONRC(12)),
     &                AIMAG(EPSILONRC(12)), BR6_DENS, BR6_DIAM, BR6_LNG
            if(l_phase) write(16,316)REAL(EPSILONRC(12)),
     &                AIMAG(EPSILONRC(12)), BR6_DENS, BR6_DIAM, BR6_LNG
            if(l_kappa) write(17,316)REAL(EPSILONRC(12)),
     &                AIMAG(EPSILONRC(12)), BR6_DENS, BR6_DIAM, BR6_LNG
            if(l_trans) write(18,316)REAL(EPSILONRC(12)),
     &                AIMAG(EPSILONRC(12)), BR6_DENS, BR6_DIAM, BR6_LNG
            write(19,316)REAL(EPSILONRC(12)),
     &                AIMAG(EPSILONRC(12)), BR6_DENS, BR6_DIAM, BR6_LNG
          else
            write(8,308) MG_BR6, RHO_BR6, BR6_DENS, BR6_DIAM, BR6_LNG
            write(7,308) MG_BR6, RHO_BR6, BR6_DENS, BR6_DIAM, BR6_LNG
            write(11,308) MG_BR6, RHO_BR6, BR6_DENS, BR6_DIAM, BR6_LNG
            write(12,308) MG_BR6, RHO_BR6, BR6_DENS, BR6_DIAM, BR6_LNG
            write(13,308) MG_BR6, RHO_BR6, BR6_DENS, BR6_DIAM, BR6_LNG
            write(14,308) MG_BR6, RHO_BR6, BR6_DENS, BR6_DIAM, BR6_LNG
            if(l_pol) write(15,308) MG_BR6, RHO_BR6, BR6_DENS, 
     &                             BR6_DIAM, BR6_LNG
            if(l_phase) write(16,308) MG_BR6, RHO_BR6, BR6_DENS, 
     &                             BR6_DIAM, BR6_LNG
            if(l_kappa) write(17,308) MG_BR6, RHO_BR6, BR6_DENS, 
     &                             BR6_DIAM, BR6_LNG
            if(l_trans) write(18,308) MG_BR6, RHO_BR6, BR6_DENS, 
     &                             BR6_DIAM, BR6_LNG
            write(19,308) MG_BR6, RHO_BR6, BR6_DENS, 
     &                             BR6_DIAM, BR6_LNG
          endif
        endif
C
            write(8,500) char,char2
            write(7,501) char,char2
            write(11,504) char,char2
            write(12,505) char,char2
            write(13,506) char,char2
            write(14,507) char,char2
            if(l_pol) write(15,508) char,char2
            if(l_phase) write(16,509) char,char2
            if(l_kappa) write(17,510) char,char2
            if(l_trans) write(18,511) char,char2
            write(19,512) char,char2
        ENDIF
c
        OPEN(UNIT=1,STATUS='SCRATCH')
        WRITE(1,CHAR_FORM) VARIABLE(INT_VAR)
        REWIND(1)
        READ(1,1) CHAR_OUT
        CLOSE(1)

c
        WRITE(8,300) CHAR_OUT,
     &                T_SIG0_dB(1,1),T_SIG0_dB(2,2),
     &               (SIG0_dB(1,1,K), SIG0_dB(2,2,K),K=1,N_SCAT_TERMS)
C
        WRITE(7,300) CHAR_OUT,
     &                T_SIG0_dB(1,2),T_SIG0_dB(2,1),
     &               (SIG0_dB(1,2,K), SIG0_dB(2,1,K),K=1,N_SCAT_TERMS)
C
        WRITE(11,450) CHAR_OUT, 
     &        EXP_CANOPY_p(1,1), EXP_CANOPY_p(2,2),
     &        EXP_CANOPY_m(1,1), EXP_CANOPY_m(2,2),
     &        EXP_KAPPA_C_p(1,1), EXP_KAPPA_C_p(2,2),
     &        EXP_KAPPA_C_m(1,1), EXP_KAPPA_C_m(2,2),
     &        EXP_KAPPA_T_p(1,1), EXP_KAPPA_T_p(2,2),
     &        EXP_KAPPA_T_m(1,1), EXP_KAPPA_T_m(2,2)
C
        WRITE(12,460) CHAR_OUT,
     &        (phase_diff_total(i),i=1,3),
     &        ((phase_diff_terms(i,j),i=1,3),j=1,n_scat_terms)
C
c        WRITE(13,470) CHAR_OUT, 
c     &       (crown_prop_diff_p(i),i=1,3),(crown_prop_diff_m(i),i=1,3),
c     &       ((crown_p_diff(i,j),i=1,3),j=1,4)
        WRITE(13,470) CHAR_OUT, 
     &        (crown_prop_diff_p(i),i=1,3),(crown_prop_diff_m(i),i=1,3),
     &        ((a_phase(i,j),i=1,3),j=1,4)
C
        WRITE(14,470) CHAR_OUT, 
     &        (gnd_back_p_diff(i),i=1,3),(gnd_spec_p_diff(i),i=1,3),
     &        (trunk_prop_diff_p(i),i=1,3),(trunk_prop_diff_m(i),i=1,3),
     &        ((trunk_spec_p_diff(i,j),i=1,3),j=1,2)

        if(l_pol)then
         WRITE(15,480) CHAR_OUT,
     &       (T(1,i),i=1,4),((BACKTERMS(1,i,j),i=1,4),j=1,3),
     &       (T(2,i),i=1,4),((BACKTERMS(2,i,j),i=1,4),j=1,3),
     &       (T(3,i),i=1,4),((BACKTERMS(3,i,j),i=1,4),j=1,3),
     &       (T(4,i),i=1,4),((BACKTERMS(4,i,j),i=1,4),j=1,3)
         WRITE(15,485) (((BACKTERMS(k,i,j),i=1,4),j=4,7),k=1,4)
        endif
        if(l_phase)then
          write(16,490) CHAR_OUT,
     &          (p_crown(1,j,1),j=1,4),(p_crown(1,j,2),j=1,4),
     &          (p_crown(1,j,3),j=1,4),(p_crown(1,j,4),j=1,4),
     &          (0.5*(p_trunk_layer(1,j,1)+p_trunk_layer(1,j,2)),j=1,4),
     &          (p_crown(2,j,1),j=1,4),(p_crown(2,j,2),j=1,4),
     &          (p_crown(2,j,3),j=1,4),(p_crown(2,j,4),j=1,4),
     &          (0.5*(p_trunk_layer(2,j,1)+p_trunk_layer(2,j,2)),j=1,4),
     &           (p_crown(3,j,1),j=1,4),(p_crown(3,j,2),j=1,4),
     &           (p_crown(3,j,3),j=1,4),(p_crown(3,j,4),j=1,4),
     &          (0.5*(p_trunk_layer(3,j,1)+p_trunk_layer(3,j,2)),j=1,4),
     &           (p_crown(4,j,1),j=1,4),(p_crown(4,j,2),j=1,4),
     &           (p_crown(4,j,3),j=1,4),(p_crown(4,j,4),j=1,4),
     &          (0.5*(p_trunk_layer(4,j,1)+p_trunk_layer(4,j,2)),j=1,4)
        endif

        if(l_kappa)then
          write(17,492) CHAR_OUT,
     &        (kappa_t_p(1,j),j=1,4), (kappa_t_m(1,j),j=1,4),
     &         (kappa_c_p(1,j),j=1,4), (kappa_c_m(1,j),j=1,4),
     &        (kappa_t_p(2,j),j=1,4), (kappa_t_m(2,j),j=1,4),
     &         (kappa_c_p(2,j),j=1,4), (kappa_c_m(2,j),j=1,4),
     &        (kappa_t_p(3,j),j=1,4), (kappa_t_m(3,j),j=1,4),
     &         (kappa_c_p(3,j),j=1,4), (kappa_c_m(3,j),j=1,4),
     &        (kappa_t_p(4,j),j=1,4), (kappa_t_m(4,j),j=1,4),
     &         (kappa_c_p(4,j),j=1,4), (kappa_c_m(4,j),j=1,4)
        endif

        if(l_trans)then
          write(18,494) CHAR_OUT,
     &            (EXP_CANOPY_p(1,j),j=1,4),(EXP_KAPPA_T_p(1,j),j=1,4),
     &                                      (EXP_KAPPA_C_p(1,j),j=1,4),
     &            (EXP_CANOPY_p(2,j),j=1,4),(EXP_KAPPA_T_p(2,j),j=1,4),
     &                                      (EXP_KAPPA_C_p(2,j),j=1,4),
     &            (EXP_CANOPY_p(3,j),j=1,4),(EXP_KAPPA_T_p(3,j),j=1,4),
     &                                      (EXP_KAPPA_C_p(3,j),j=1,4),
     &            (EXP_CANOPY_p(4,j),j=1,4),(EXP_KAPPA_T_p(4,j),j=1,4),
     &                                      (EXP_KAPPA_C_p(4,j),j=1,4)
        endif

        write(19,495) CHAR_OUT,
     &    (GRND_REFLECT_MAT(1,j),j=1,4),M_crown_p(1,1),M_crown_p(1,2),
     &    (GRND_REFLECT_MAT(2,j),j=1,4),M_crown_p(2,1),M_crown_p(2,2),
     &    (GRND_REFLECT_MAT(3,j),j=1,4),
     &                     M_TRUNK_LAYER_p(1,1),M_TRUNK_LAYER_p(1,2),
     &    (GRND_REFLECT_MAT(4,j),j=1,4),
     &                     M_TRUNK_LAYER_p(2,1),M_TRUNK_LAYER_p(2,2)

c
1       FORMAT(A11)
5       FORMAT(' MIMICS Version ',a5)
c
10          format(/,' ****  Forest Canopy Backscatter Model Output File 
     & -- Like-polarized Sigma0 Values ****')
15          format(/,' ****  Forest Canopy Backscatter Model Output File 
     & -- Cross-polarized Sigma0 Values (receive-transmit) ****')
20          format(/,' ****  Forest Canopy Backscatter Model Diagnostic 
     &Output File ****')
30          format(/,' ****  Forest Canopy Backscatter Model Output File 
     & -- Phase Difference Values ****')
40          format(/,' ****  Forest Canopy Backscatter Model Output File 
     & -- Stokes Matrix Values ****')
50          format(/,' ****  Forest Canopy Backscatter Model Output File 
     & -- Extinction Matrix Values ****')
60          format(/,' ****  Forest Canopy Backscatter Model Output File 
     & -- Transmissivity Values ****')
65          format(/,' ****  Forest Canopy Backscatter Model Output File 
     & -- Phase Difference Values -- Total Contributions from Individual
     & Terms ****')
66          format(/,' ****  Forest Canopy Backscatter Model Output File 
     & -- Phase Difference Values -- Individual Constituent Contribution
     &s of the Crown ****')
c
67          format(/,' ****  Forest Canopy Backscatter Model Output File 
     & -- Phase Difference Values -- Individual Constituent Contribution
     &s of the Trunk and Ground ****')
68          format(/,' ****  Forest Canopy Backscatter Model Output File 
     & -- Canopy Transformation Matrices -- Total and Indiviual Terms **
     &**')
69          format(/,' ****  Forest Canopy Backscatter Model Output File 
     & -- Phase Matrices of Indiviual Layers ****')
70          format(/,' ****  Forest Canopy Backscatter Model Output File 
     & -- Kappa Matrices of Indiviual Layers ****')
71          format(/,' ****  Forest Canopy Backscatter Model Output File 
     & -- Transmissivity Matrices of Indiviual Layers ****')
72          format(/,' ****  Forest Canopy Backscatter Model Output File 
     & -- Misc. Values for MIMICS II input ****')
c

c
100     format(' Trunks = ',L1,
     &         ', Primary Branches = ',L1,
     &         ', Secondary Branches = ',L1,
     &         ', 3rd Branches = ',L1,
     &         ', 4th Branches = ',L1,
     &         ', 5th Branches = ',L1,
     &         ', 6th Branches = ',L1,
     &         ', Leaves = ',L1,
     &         ', Needles = ',L1)
c
105     format(' SOIL_SURFACE = ',L1,', WATER_SURFACE = ',L1,
     &         ', SNOW_SURFACE = ',L1,', ICE_SURFACE = ',L1)
c
107     format(' Dielectric Lookup Tables: Trunk = ',L1,
     &         ', Primary Branch = ',L1,
     &         ', Secondary Branch = ',L1,
     &         ', 3rd Branch = ',L1,
     &         ', 4th Branch = ',L1,
     &         ', 5th Branch = ',L1,
     &         ', 6th Branch = ',L1,
     &         ', Leaf = ',L1,
     &         ', Needle = ',L1,
     &         ', Ground = ',L1)
C
110     FORMAT(/,/,' Nesting Data:',/,
     &         '  Variable Number = ',30(1x,i3))
111     FORMAT('    Nesting Order = ',30(1x,i3))
112     FORMAT('  Number of Loops = ',30(1x,i3))
c
140     format(/,' Family Number ',i5,/)
c
150     FORMAT(' RADAR Look Angle = ',f4.1,' Degrees,',
     &         ' Radar Frequency = ',f6.2,' GHz')
151     FORMAT(' Canopy Density = ',f5.2,' Trees per square meter,',
     &         ' Crown Height = ',f5.2,' meters,',
     &         ' Trunk Height = ',f5.2,' meters,',
     &         ' Vegetation Temperature = ',f6.2,' degrees C.')
152     FORMAT(' Soil Volumetric Moisture = ',f5.3,
     &    ', Soil RMS Roughness = ',f5.2,' cm',
     &    ', Soil Correlation Length = ',f6.2,' cm',
     &    ', Soil % Sand = ',f5.2,' cm, Soil % Clay = ',f5.2,' cm',
     &    ', Soil Temperature = ',f6.2,' degrees C., Model Type = ',a30)
153     FORMAT(' Surface Water Salt Content = ',f6.3,' ppt',
     &         ', Water Temperature = ',f6.2,' degrees C.')
154     FORMAT(' Trunk Gravimetric Moisture = ',f5.3,
     &         ', Trunk Dry Density = ',f5.3,' cm,',
     &         ' Trunk Diameter = ',f6.3,' cm')
155     FORMAT(' Leaf Gravimetric Moisture = ',f5.3,
     &         ', Leaf Dry Density = ',f5.3,' cm',
     &         ', Leaf Density = ',f7.2,' leaves per cubic meter',
     &         ', Leaf Diameter = ',f6.3,' cm',
     &         ', Leaf Thickness = ',f5.3,' cm')
156     FORMAT(' Needle Gravimetric Moisture = ',f5.3,
     &         ', Needle Dry Density = ',f5.3,' cm',
     &         ', Needle Density = ',e11.4,' needles per cubic meter',
     &         ', Needle Diameter = ',f5.3,' cm',
     &         ', Needle Length = ',f6.3,' cm')
157     FORMAT(' Branch Gravimetric Moisture = ',f5.3,
     &         ', Branch Dry Density = ',f5.3,' cm',
     &         ', Branch Density = ',f7.3,' branches per cubic meter',
     &         ', Branch Diameter = ',f6.3,' cm',
     &         ', Branch Length = ',f7.3,' meters')
158     FORMAT(' Secondary Branch Grav Moist = ',f5.3,
     &         ', Sec Branch Dry Dens= ',f5.3,' cm',
     &         ', Sec Branch Dens= ',f7.3,' branches per cubic meter',
     &         ', Sec Branch Diam = ',f6.3,' cm',
     &         ', Sec Branch Lng= ',f7.3,' meters')
C
302     FORMAT(' 3rd Branch Grav Moist = ',f5.3,
     &         ', 3rd Branch Dry Dens= ',f5.3,' cm',
     &         ', 3rd Branch Dens= ',f7.3,' branches per cubic meter',
     &         ', 3rd Branch Diam = ',f6.3,' cm',
     &         ', 3rd Branch Lng= ',f7.3,' meters')
C
304     FORMAT(' 4th Branch Grav Moist = ',f5.3,
     &         ', 4th Branch Dry Dens= ',f5.3,' cm',
     &         ', 4th Branch Dens= ',f7.3,' branches per cubic meter',
     &         ', 4th Branch Diam = ',f6.3,' cm',
     &         ', 4th Branch Lng= ',f7.3,' meters')
C
306     FORMAT(' 5th Branch Grav Moist = ',f5.3,
     &         ', 5th Branch Dry Dens= ',f5.3,' cm',
     &         ', 5th Branch Dens= ',f7.3,' branches per cubic meter',
     &         ', 5th Branch Diam = ',f6.3,' cm',
     &         ', 5th Branch Lng= ',f7.3,' meters')
C
308     FORMAT(' 6th Branch Grav Moist = ',f5.3,
     &         ', 6th Branch Dry Dens= ',f5.3,' cm',
     &         ', 6th Branch Dens= ',f7.3,' branches per cubic meter',
     &         ', 6th Branch Diam = ',f6.3,' cm',
     &         ', 6th Branch Lng= ',f7.3,' meters')
C
159     FORMAT(' Snow Dielectric = ',f6.3,' -j',f6.3,
     &         ' Snow Layer Thickness = ',f6.3,' m')
160     FORMAT(' Surface Dielectric = ',f6.3,' -j',f6.3,
     &   ', Soil RMS Roughness = ',f5.2,' cm',
     &   ', Soil Correlation Length = ',f6.2,' cm, Model Type = ',a30)
161     FORMAT(' Surface Dielectric = ',f6.3,' -j',f6.3)
164     FORMAT(' Trunk Dielectric = ',f6.3,' -j',f6.3,
     &         ' Trunk Diameter = ',f6.3,' cm')
165     FORMAT(' Leaf  Dielectric = ',f6.3,' -j',f6.3,
     &         ', Leaf Density = ',f7.2,' leaves per cubic meter',
     &         ', Leaf Diameter = ',f6.3,' cm',
     &         ', Leaf Thickness = ',f5.3,' cm')
166     FORMAT(' Needle Dielectric = ',f6.3,' -j',f6.3,
     &         ', Needle Density = ',e11.4,' needles per cubic meter',
     &         ', Needle Diameter = ',f5.3,' cm',
     &         ', Needle Length = ',f6.3,' cm')
167     FORMAT(' Branch Dielectric = ',f6.3,' -j',f6.3,
     &         ', Branch Density = ',f7.3,' branches per cubic meter',
     &         ', Branch Diameter = ',f6.3,' cm',
     &         ', Branch Length = ',f7.3,' meters')
168     FORMAT(' Secondary Branch Dielectric = ',f6.3,' -j',f6.3,
     &         ', Sec Branch Dens= ',f7.3,' branches per cubic meter',
     &         ', Sec Branch Diam = ',f6.3,' cm',
     &         ', Sec Branch Lng= ',f7.3,' meters')
C
312     FORMAT(' 3rd Branch Dielectric = ',f6.3,' -j',f6.3,
     &         ', 3rd Branch Dens= ',f7.3,' branches per cubic meter',
     &         ', 3rd Branch Diam = ',f6.3,' cm',
     &         ', 3rd Branch Lng= ',f7.3,' meters')
C
314     FORMAT(' 4th Branch Dielectric = ',f6.3,' -j',f6.3,
     &         ', 4th Branch Dens= ',f7.3,' branches per cubic meter',
     &         ', 4th Branch Diam = ',f6.3,' cm',
     &         ', 4th Branch Lng= ',f7.3,' meters')
C
316     FORMAT(' 5th Branch Dielectric = ',f6.3,' -j',f6.3,
     &         ', 5th Branch Dens= ',f7.3,' branches per cubic meter',
     &         ', 5th Branch Diam = ',f6.3,' cm',
     &         ', 5th Branch Lng= ',f7.3,' meters')
C
318     FORMAT(' 6th Branch Dielectric = ',f6.3,' -j',f6.3,
     &         ', 6th Branch Dens= ',f7.3,' branches per cubic meter',
     &         ', 6th Branch Diam = ',f6.3,' cm',
     &         ', 6th Branch Lng= ',f7.3,' meters')
C
199     FORMAT(/,' CONSTITUENT ORIENTATION AND SIZE DISTRIBUTIONS:')
200     FORMAT(1X,'  Trunk Orient PDF = ',i2,', Form = ',a30,
     &        '  Size PDF = ',i2,', Form = ',a30)
201     FORMAT(1X,'   Leaf Orient PDF = ',i2,', Form = ',a30,
     &        '  Size PDF = ',i2,', Form = ',a30)
202     FORMAT(1X,' Needle Orient PDF = ',i2,', Form = ',a30,
     &        '  Size PDF = ',i2,', Form = ',a30)
203     FORMAT(1X,'Prim Br Orient PDF = ',i2,', Form = ',a30,
     &        '  Size PDF = ',i2,', Form = ',a30)
204     FORMAT(1X,' Sec Br Orient PDF = ',i2,', Form = ',a30,
     &        '  Size PDF = ',i2,', Form = ',a30)
C
205     FORMAT(1X,' 3rd Br Orient PDF = ',i2,', Form = ',a30,
     &        '  Size PDF = ',i2,', Form = ',a30)
C
206     FORMAT(1X,' 4th Br Orient PDF = ',i2,', Form = ',a30,
     &        '  Size PDF = ',i2,', Form = ',a30)
C
207     FORMAT(1X,' 5th Br Orient PDF = ',i2,', Form = ',a30,
     &        '  Size PDF = ',i2,', Form = ',a30)
C
208     FORMAT(1X,' 6th Br Orient PDF = ',i2,', Form = ',a30,
     &        '  Size PDF = ',i2,', Form = ',a30)
C
500     format(' **** ','Backscattering Cross Section - Sigma0 (dB)',
     &      ' ****',/,
     &      1x,a9,5x,
     &      'Total Backscatter', 3x,'Gnd-Crown-Gnd',7x,
     &      'Crown-Ground',      7x,'Ground-Crown', 7x,
     &      'Direct Crown',      7x,'Trunk-Ground', 7x,
     &      'Ground-Trunk',      7x,'Direct Ground',/,
     &      1x,a9,2x,8(6x,'VV',6x,'HH',3x))
c
501     format(' **** ','Backscattering Cross Section - Sigma0 (dB)',
     &      ' ****',/,
     &      1x,a9,5x,
     &      'Total Backscatter', 3x,'Gnd-Crown-Gnd',7x,
     &      'Crown-Ground',      7x,'Ground-Crown', 7x,
     &      'Direct Crown',      7x,'Trunk-Ground', 7x,
     &      'Ground-Trunk',      7x,'Direct Ground',/,
     &      1x,a9,2x,8(6x,'VH',6x,'HV',3x))
C
502     format(' **** ','Backscattering Stokes Matrices',' ****',/,
     &      1x,a9,3x,
     &      'Total Backscatter', 7x,'Ground-Crown-Ground',10x,
     &      'Crown-Ground',      15x,'Ground-Crown',13x,
     &      'Direct Crown',      15x,'Trunk-Ground',14x,
     &      'Ground-Trunk',      13x,'Direct Ground',/,
     &      1x,a9)
c
503     format (' **** ','Extinction Matrices',' ****',/,
     &       1x,a9,16x,
     &      ' Total Canopy',40x,'Crown Layer',41x,'Trunk Layer',/,
     &      1x,a9,5x,
     &      ' + going ',17x,' - going ',17x,
     &      ' + going ',17x,' - going ',17x,
     &      ' + going ',17x,' - going ')
c
504     format(' **** ','Transmissivity Values (Tau)',' ****',/,
     &       1x,a9,20x,
     &      ' Total Canopy', 40x,'Crown Layer',41x,'Trunk Layer',/,
     &      1x,a9,
     &      3(7x,'Tau_v+',6x,'Tau_h+',8x,'Tau_v-',6x,'Tau_h-',1x))
c
505     format(' **** ',
     &      'Backscattering Phase Difference (tr-VV) (Degrees)',
     &      ' ****',/,1x,a9,5x,
     &      'Total Backscatter', 7x,'Ground-Crown-Ground',9x,
     &      'Crown-Ground',11x,'Ground-Crown',13x,
     &      'Direct Crown',11x,'Trunk-Ground',13x,
     &      'Ground-Trunk',10x,'Direct Ground',/,
     &      1x,a9,8(6x,'HH',5x,'VH',5x,'HV',2x))
c
506    format(' **** Backscattering Phase Difference',
     &      '  -- Crown Contributions (tr-VV) (Degrees) ****',/,
     &      1x,a9,5x,
     &      'Crown Propagation +',5x,'Crown Propagation -',4x,
     &      'Ground-Crown-Ground',8x,'Specular Crown -',8x,
     &      'Specular Crown +',10x,'Direct Crown',
     &       /,1x,a9,6(7x,'HH',5x,'VH',5x,'HV',1x))
c
507    format(' **** Backscattering Phase Difference',
     &      '  -- Trunk and Ground Contributions (tr-VV) (Degrees) ****'
     &      ,/,1x,a9,10x,'Ground Back',10x,'Specular Ground',7x,
     &      'Trunk Propagation +',4x,'Trunk Propagation -',7x,
     &      'Specular Trunk +',8x,'Specular Trunk -',
     &       /,1x,a9,6(7x,'HH',5x,'VH',5x,'HV',1x))
508    format(' **** Backscattering Transformation Matrices ****'
     &      ,/,1x,a9,/,1x,a9)
509    format(' **** Phase Matrices ****'
     &   ,/,1x,a9,15x,'Ground-Crown-Ground',16x,19x,'Crown-Ground',19x,
     &            19x,'Ground-Crown',19x,19x,'Direct Crown',19x,
     &            18x,'Specular Trunk'
     &      ,/,1x,a9)
510    format(' **** Extinction Matrices ****'
     &   ,/,1x,a9,44x,'Trunk Layer',45x,44x,'Crown Layer',
     &    /,1x,a9,21x,'positive-going',19x,18x,'negative-going',18x,
     &            18x,'positive-going',18x,18x,'negative-going')
511    format(' **** Transmissivity Matrices  (Positive-Going)****'
     &   ,/,1x,a9,19x,'Total Canopy',19x,19x,'Trunk Layer',20x,
     &            21x,'Crown Layer',/,1x,a9)
512    format(' **** Ground Reflection Matrix *** Crown M matrices ****'
     &   ,/,1x,a9,19x,'Ground Reflection',19x,19x,
     &   ' M-matrices (crown & trunk)',/,1x,a9)
c

c
300     format(A11,8(3X,F7.2,2X,F7.2))
400     format(1x,f6.3)
401     format(6x,8(3x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4))
411     format(6x,6(3x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4))
450     format(A11,6(3x,E11.4,1x,E11.4))
460     format(A11,8(3x,3(f6.1,1x)))
470     format(A11,6(3x,3(f6.1,1x)))
480     format(11x,21x,'Total',38x,'Ground-Crown-Ground',35x,
     &     'Ground-Crown',38x,'Crown-Ground',/,
     &        A11,4(4(E11.4,1x),2x),/,
     &      3(11x,4(4(E11.4,1x),2x),/))
485     format(11x,17x,'Direct Crown',39x,'Trunk-Ground',37x,
     &       'Ground-Trunk',38x,'Direct Ground',/,
     &      4(11x,4(4(E11.4,1x),2x),/))
c
c
490     format(A11,5(2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4),/,
     &         11x,5(2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4),/,
     &         11x,5(2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4),/,
     &         11x,5(2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4),/)
c
492     format(A11,4(3x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4),/,
     &         11x,4(3x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4),/,
     &         11x,4(3x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4),/,
     &         11x,4(3x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4),/)
c
494     format(A11,3(3x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4),/,
     &         11x,3(3x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4),/,
     &         11x,3(3x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4),/,
     &         11x,3(3x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4),/)
c
495     format(A11,2(3x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4),/,
     &         11x,2(3x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4),/,
     &         11x,2(3x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4),/,
     &         11x,2(3x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4),/)

        return
        end
c
        SUBROUTINE FINISH
        save

c
        LOGICAL L_POL, L_PHASE, L_KAPPA, L_TRANS
        COMMON /L_FILES/ L_POL, L_PHASE, L_KAPPA, L_TRANS

c
        close(unit=7)
        close(unit=8)
        close(unit=11)
        close(unit=12)
        close(unit=13)
        close(unit=14)
        if(l_pol) close(unit=15)
        if(l_phase) close(unit=16)
        if(l_kappa) close(unit=17)
        if(l_trans) close(unit=18)
        close(unit=19)

C
        RETURN
        END
