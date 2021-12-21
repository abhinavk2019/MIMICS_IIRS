        SUBROUTINE DIELECTRIC(EPSILONR,EPSILONRC)
        save
c
        include 'parameters.include'
c
        REAL MV_SOIL, SAND, CLAY
        REAL SALT
        REAL MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
        REAL MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
        REAL MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
        REAL MG_BR1, RHO_BR1, BR1_DENS, BR1_DIAM, BR1_LNG
        REAL MG_BR2, RHO_BR2, BR2_DENS, BR2_DIAM, BR2_LNG
        REAL T_SOIL, T_WATER, T_VEG
C
        COMPLEX EPSILONR(N_EPS),EPSILONRC(N_EPS)
C
        LOGICAL LOG_CONSTITUENT(N_CONSTITUENTS)
        LOGICAL LOG_EPS_TABLE(N_CONSTITUENTS)
        LOGICAL LOG_DRY_DENSITY(N_CONSTITUENTS) 
        LOGICAL LOG_PDF_TYPE(N_CONSTITUENTS,N_CONST_VARY)
        LOGICAL LOG_HIST(N_CONSTITUENTS)
C
        COMMON /L_CONFIGURE/ LOG_CONSTITUENT, LOG_EPS_TABLE,
     &                       LOG_DRY_DENSITY, LOG_PDF_TYPE, LOG_HIST

C
        LOGICAL OPEN, STEP_VARIABLE(N_VARIABLES)
        LOGICAL STEP_THIS_TIME(N_VARIABLES)
        LOGICAL  CALL_SUB(N_CALLS,N_SUB_CALLS)
        COMMON /L_FLAGS/ STEP_VARIABLE, STEP_THIS_TIME, CALL_SUB, OPEN
C
        COMMON /R_ENVIRON/ T_SOIL, T_WATER, T_VEG
C
        COMMON /R_SOIL/ MV_SOIL, SAND, CLAY
        COMMON /R_WATER/ SALT
C
        COMMON /R_TRUNK/ MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
        COMMON /R_LEAF/ MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
        COMMON /R_NDL/ MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
        COMMON /R_BR1/ MG_BR1, RHO_BR1, BR1_DENS, BR1_DIAM, BR1_LNG
        COMMON /R_BR2/ MG_BR2, RHO_BR2, BR2_DENS, BR2_DIAM, BR2_LNG

        REAL MG_BR3, RHO_BR3, BR3_DENS, BR3_DIAM, BR3_LNG
        COMMON /R_BR3/ MG_BR3, RHO_BR3, BR3_DENS, BR3_DIAM, BR3_LNG

        REAL MG_BR4, RHO_BR4, BR4_DENS, BR4_DIAM, BR4_LNG
        COMMON /R_BR4/ MG_BR4, RHO_BR4, BR4_DENS, BR4_DIAM, BR4_LNG

        REAL MG_BR5, RHO_BR5, BR5_DENS, BR5_DIAM, BR5_LNG
        COMMON /R_BR5/ MG_BR5, RHO_BR5, BR5_DENS, BR5_DIAM, BR5_LNG

        REAL MG_BR6, RHO_BR6, BR6_DENS, BR6_DIAM, BR6_LNG
        COMMON /R_BR6/ MG_BR6, RHO_BR6, BR6_DENS, BR6_DIAM, BR6_LNG

c
        IF(CALL_SUB(1,2)) THEN
            CALL EPS_WATER(EPSILONR(1),T_WATER)
            IF(AIMAG(EPSILONR(1)).GT.0.) EPSILONR(1)=CONJG(EPSILONR(1))
            EPSILONRC(1) = CONJG(EPSILONR(1))
        ENDIF

c
        if(CALL_SUB(1,3)) THEN
            CALL EPS_SOIL(EPSILONR(2), T_SOIL)
            IF(AIMAG(EPSILONR(2)).GT.0.) EPSILONR(2)=CONJG(EPSILONR(2))
            EPSILONRC(2) = CONJG(EPSILONR(2))
        ENDIF

c
        if(CALL_SUB(1,4)) THEN
            CALL EPS_NEEDLE(EPSILONR(3))
            IF(AIMAG(EPSILONR(3)).GT.0.) EPSILONR(3)=CONJG(EPSILONR(3))
            EPSILONRC(3) = CONJG(EPSILONR(3))
        ENDIF

c
        if(CALL_SUB(1,5)) THEN
            CALL EPS_LEAF(EPSILONR(4))
            IF(AIMAG(EPSILONR(4)).GT.0.) EPSILONR(4)=CONJG(EPSILONR(4))
            EPSILONRC(4) = CONJG(EPSILONR(4))
        ENDIF

c
        if(CALL_SUB(1,6)) THEN
            CALL EPS_TRUNK(EPSILONR(5))
            IF(AIMAG(EPSILONR(5)).GT.0.) EPSILONR(5)=CONJG(EPSILONR(5))
            EPSILONRC(5) = CONJG(EPSILONR(5))
        ENDIF

c
        if(CALL_SUB(1,7)) THEN
            CALL EPS_BRANCH(EPSILONR(6))
            IF(AIMAG(EPSILONR(6)).GT.0.) EPSILONR(6)=CONJG(EPSILONR(6))
            EPSILONRC(6) = CONJG(EPSILONR(6))
        ENDIF
c
        if(CALL_SUB(1,8)) THEN
            CALL EPS_BRANCH2(EPSILONR(7))
            IF(AIMAG(EPSILONR(7)).GT.0.) EPSILONR(7)=CONJG(EPSILONR(7))
            EPSILONRC(7) = CONJG(EPSILONR(7))
        ENDIF
c
        if(CALL_SUB(1,9)) THEN
            IF(AIMAG(EPSILONR(8)).GT.0.) EPSILONR(8)=CONJG(EPSILONR(8))
            EPSILONRC(8) = CONJG(EPSILONR(8))
        ENDIF
c
        if(CALL_SUB(1,10)) THEN
            CALL EPS_BRANCH3(EPSILONR(9))
            IF(AIMAG(EPSILONR(9)).GT.0.) EPSILONR(9)=CONJG(EPSILONR(9))
            EPSILONRC(9) = CONJG(EPSILONR(9))
        ENDIF
c
        if(CALL_SUB(1,11)) THEN
          CALL EPS_BRANCH4(EPSILONR(10))
          IF(AIMAG(EPSILONR(10)).GT.0.) EPSILONR(10)=CONJG(EPSILONR(10))
          EPSILONRC(10) = CONJG(EPSILONR(10))
        ENDIF
c
        if(CALL_SUB(1,12)) THEN
          CALL EPS_BRANCH5(EPSILONR(11))
          IF(AIMAG(EPSILONR(11)).GT.0.) EPSILONR(11)=CONJG(EPSILONR(11))
          EPSILONRC(11) = CONJG(EPSILONR(11))
        ENDIF
c
        if(CALL_SUB(1,13)) THEN
          CALL EPS_BRANCH6(EPSILONR(12))
          IF(AIMAG(EPSILONR(12)).GT.0.) EPSILONR(12)=CONJG(EPSILONR(12))
          EPSILONRC(12) = CONJG(EPSILONR(12))
        ENDIF
c
        RETURN
        END
c
        SUBROUTINE EPS_LEAF(EPSR)
        save
c
        include 'parameters.include'
c
        REAL MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
        REAL T_SOIL, T_WATER, T_VEG
        COMPLEX EPSR
        LOGICAL LOG_CONSTITUENT(N_CONSTITUENTS)
        LOGICAL LOG_EPS_TABLE(N_CONSTITUENTS)
        LOGICAL LOG_DRY_DENSITY(N_CONSTITUENTS) 
        LOGICAL LOG_PDF_TYPE(N_CONSTITUENTS,N_CONST_VARY)
        LOGICAL LOG_HIST(N_CONSTITUENTS)
C
        COMMON /L_CONFIGURE/ LOG_CONSTITUENT, LOG_EPS_TABLE,
     &                       LOG_DRY_DENSITY, LOG_PDF_TYPE, LOG_HIST
C
        COMMON /R_ENVIRON/ T_SOIL, T_WATER, T_VEG
        COMMON /R_LEAF/ MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
c
        IF(LOG_DRY_DENSITY(8))THEN
            CALL EPSVEG_Mg_2(MG_LEAF,RHO_LEAF,T_VEG,EPSR)
        ELSE
            CALL EPSVEG_Mg_1(MG_LEAF,T_VEG,EPSR)
        ENDIF
C
        RETURN
        END
c
        SUBROUTINE EPS_NEEDLE(EPSR)
        save
c
        include 'parameters.include'
c
        REAL MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
        REAL T_SOIL, T_WATER, T_VEG
        COMPLEX EPSR
        LOGICAL LOG_CONSTITUENT(N_CONSTITUENTS)
        LOGICAL LOG_EPS_TABLE(N_CONSTITUENTS)
        LOGICAL LOG_DRY_DENSITY(N_CONSTITUENTS) 
        LOGICAL LOG_PDF_TYPE(N_CONSTITUENTS,N_CONST_VARY)
        LOGICAL LOG_HIST(N_CONSTITUENTS)
C
        COMMON /L_CONFIGURE/ LOG_CONSTITUENT, LOG_EPS_TABLE,
     &                       LOG_DRY_DENSITY, LOG_PDF_TYPE, LOG_HIST
C
        COMMON /R_ENVIRON/ T_SOIL, T_WATER, T_VEG
        COMMON /R_NDL/ MG_NDL, RHO_NDL, NDL_DENS, NDL_DIAM, NDL_LNG
c
        IF(LOG_DRY_DENSITY(9))then
            CALL EPSVEG_Mg_2(MG_NDL,RHO_NDL,T_VEG,EPSR)
        ELSE
            CALL EPSVEG_Mg_1(MG_NDL,T_VEG,EPSR)
        ENDIF
C
        RETURN
        END
c
        SUBROUTINE EPS_BRANCH(EPSR)
        save
c
        include 'parameters.include'
c
        REAL MG_BR1, RHO_BR1, BR1_DENS, BR1_DIAM, BR1_LNG
        REAL T_SOIL, T_WATER, T_VEG
        COMPLEX EPSR
        LOGICAL LOG_CONSTITUENT(N_CONSTITUENTS)
        LOGICAL LOG_EPS_TABLE(N_CONSTITUENTS)
        LOGICAL LOG_DRY_DENSITY(N_CONSTITUENTS) 
        LOGICAL LOG_PDF_TYPE(N_CONSTITUENTS,N_CONST_VARY)
        LOGICAL LOG_HIST(N_CONSTITUENTS)
C
        COMMON /L_CONFIGURE/ LOG_CONSTITUENT, LOG_EPS_TABLE,
     &                       LOG_DRY_DENSITY, LOG_PDF_TYPE, LOG_HIST
C
        COMMON /R_ENVIRON/ T_SOIL, T_WATER, T_VEG
        COMMON /R_BR1/ MG_BR1, RHO_BR1, BR1_DENS, BR1_DIAM, BR1_LNG
c
        IF(LOG_DRY_DENSITY(2))THEN
            CALL EPSVEG_Mg_2(MG_BR1,RHO_BR1,T_VEG,EPSR)
        ELSE
            CALL EPSVEG_Mg_1(MG_BR1,T_VEG,EPSR)
        ENDIF
C
        RETURN
        END
c
        SUBROUTINE EPS_BRANCH2(EPSR)
        save
c
        include 'parameters.include'
c
        REAL MG_BR2, RHO_BR2, BR2_DENS, BR2_DIAM, BR2_LNG
        REAL T_SOIL, T_WATER, T_VEG
        COMPLEX EPSR
        LOGICAL LOG_CONSTITUENT(N_CONSTITUENTS)
        LOGICAL LOG_EPS_TABLE(N_CONSTITUENTS)
        LOGICAL LOG_DRY_DENSITY(N_CONSTITUENTS) 
        LOGICAL LOG_PDF_TYPE(N_CONSTITUENTS,N_CONST_VARY)
        LOGICAL LOG_HIST(N_CONSTITUENTS)
C
        COMMON /L_CONFIGURE/ LOG_CONSTITUENT, LOG_EPS_TABLE,
     &                       LOG_DRY_DENSITY, LOG_PDF_TYPE, LOG_HIST
C
        COMMON /R_ENVIRON/ T_SOIL, T_WATER, T_VEG
        COMMON /R_BR2/ MG_BR2, RHO_BR2, BR2_DENS, BR2_DIAM, BR2_LNG
c
        IF(LOG_DRY_DENSITY(3))THEN
            CALL EPSVEG_Mg_2(MG_BR2,RHO_BR2,T_VEG,EPSR)
        ELSE
            CALL EPSVEG_Mg_1(MG_BR2,T_VEG,EPSR)
        ENDIF
C
        RETURN
        END





c
        SUBROUTINE EPS_BRANCH3(EPSR)
        save
c
        include 'parameters.include'
c
        REAL MG_BR3, RHO_BR3, BR3_DENS, BR3_DIAM, BR3_LNG
        REAL T_SOIL, T_WATER, T_VEG
        COMPLEX EPSR
        LOGICAL LOG_CONSTITUENT(N_CONSTITUENTS)
        LOGICAL LOG_EPS_TABLE(N_CONSTITUENTS)
        LOGICAL LOG_DRY_DENSITY(N_CONSTITUENTS) 
        LOGICAL LOG_PDF_TYPE(N_CONSTITUENTS,N_CONST_VARY)
        LOGICAL LOG_HIST(N_CONSTITUENTS)
C
        COMMON /L_CONFIGURE/ LOG_CONSTITUENT, LOG_EPS_TABLE,
     &                       LOG_DRY_DENSITY, LOG_PDF_TYPE, LOG_HIST
C
        COMMON /R_ENVIRON/ T_SOIL, T_WATER, T_VEG
        COMMON /R_BR3/ MG_BR3, RHO_BR3, BR3_DENS, BR3_DIAM, BR3_LNG
c
        IF(LOG_DRY_DENSITY(4))THEN
            CALL EPSVEG_Mg_2(MG_BR3,RHO_BR3,T_VEG,EPSR)
        ELSE
            CALL EPSVEG_Mg_1(MG_BR3,T_VEG,EPSR)
        ENDIF
C
        RETURN
        END
c
        SUBROUTINE EPS_BRANCH4(EPSR)
        save
c
        include 'parameters.include'
c
        REAL MG_BR4, RHO_BR4, BR4_DENS, BR4_DIAM, BR4_LNG
        REAL T_SOIL, T_WATER, T_VEG
        COMPLEX EPSR
        LOGICAL LOG_CONSTITUENT(N_CONSTITUENTS)
        LOGICAL LOG_EPS_TABLE(N_CONSTITUENTS)
        LOGICAL LOG_DRY_DENSITY(N_CONSTITUENTS) 
        LOGICAL LOG_PDF_TYPE(N_CONSTITUENTS,N_CONST_VARY)
        LOGICAL LOG_HIST(N_CONSTITUENTS)
C
        COMMON /L_CONFIGURE/ LOG_CONSTITUENT, LOG_EPS_TABLE,
     &                       LOG_DRY_DENSITY, LOG_PDF_TYPE, LOG_HIST
C
        COMMON /R_ENVIRON/ T_SOIL, T_WATER, T_VEG
        COMMON /R_BR4/ MG_BR4, RHO_BR4, BR4_DENS, BR4_DIAM, BR4_LNG
c
        IF(LOG_DRY_DENSITY(5))THEN
            CALL EPSVEG_Mg_2(MG_BR4,RHO_BR4,T_VEG,EPSR)
        ELSE
            CALL EPSVEG_Mg_1(MG_BR4,T_VEG,EPSR)
        ENDIF
C
        RETURN
        END
c
        SUBROUTINE EPS_BRANCH5(EPSR)
        save
c
        include 'parameters.include'
c
        REAL MG_BR5, RHO_BR5, BR5_DENS, BR5_DIAM, BR5_LNG
        REAL T_SOIL, T_WATER, T_VEG
        COMPLEX EPSR
        LOGICAL LOG_CONSTITUENT(N_CONSTITUENTS)
        LOGICAL LOG_EPS_TABLE(N_CONSTITUENTS)
        LOGICAL LOG_DRY_DENSITY(N_CONSTITUENTS) 
        LOGICAL LOG_PDF_TYPE(N_CONSTITUENTS,N_CONST_VARY)
        LOGICAL LOG_HIST(N_CONSTITUENTS)
C
        COMMON /L_CONFIGURE/ LOG_CONSTITUENT, LOG_EPS_TABLE,
     &                       LOG_DRY_DENSITY, LOG_PDF_TYPE, LOG_HIST
C
        COMMON /R_ENVIRON/ T_SOIL, T_WATER, T_VEG
        COMMON /R_BR5/ MG_BR5, RHO_BR5, BR5_DENS, BR5_DIAM, BR5_LNG
c
        IF(LOG_DRY_DENSITY(6))THEN
            CALL EPSVEG_Mg_2(MG_BR5,RHO_BR5,T_VEG,EPSR)
        ELSE
            CALL EPSVEG_Mg_1(MG_BR5,T_VEG,EPSR)
        ENDIF
C
        RETURN
        END
c
        SUBROUTINE EPS_BRANCH6(EPSR)
        save
c
        include 'parameters.include'
c
        REAL MG_BR6, RHO_BR6, BR6_DENS, BR6_DIAM, BR6_LNG
        REAL T_SOIL, T_WATER, T_VEG
        COMPLEX EPSR
        LOGICAL LOG_CONSTITUENT(N_CONSTITUENTS)
        LOGICAL LOG_EPS_TABLE(N_CONSTITUENTS)
        LOGICAL LOG_DRY_DENSITY(N_CONSTITUENTS) 
        LOGICAL LOG_PDF_TYPE(N_CONSTITUENTS,N_CONST_VARY)
        LOGICAL LOG_HIST(N_CONSTITUENTS)
C
        COMMON /L_CONFIGURE/ LOG_CONSTITUENT, LOG_EPS_TABLE,
     &                       LOG_DRY_DENSITY, LOG_PDF_TYPE, LOG_HIST
C
        COMMON /R_ENVIRON/ T_SOIL, T_WATER, T_VEG
        COMMON /R_BR6/ MG_BR6, RHO_BR6, BR6_DENS, BR6_DIAM, BR6_LNG
c
        IF(LOG_DRY_DENSITY(7))THEN
            CALL EPSVEG_Mg_2(MG_BR6,RHO_BR6,T_VEG,EPSR)
        ELSE
            CALL EPSVEG_Mg_1(MG_BR6,T_VEG,EPSR)
        ENDIF
C
        RETURN
        END








c
        SUBROUTINE EPS_TRUNK(EPSR)
        save
c
        include 'parameters.include'
c
        REAL MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
        REAL T_SOIL, T_WATER, T_VEG
        COMPLEX EPSR
        LOGICAL LOG_CONSTITUENT(N_CONSTITUENTS)
        LOGICAL LOG_EPS_TABLE(N_CONSTITUENTS)
        LOGICAL LOG_DRY_DENSITY(N_CONSTITUENTS) 
        LOGICAL LOG_PDF_TYPE(N_CONSTITUENTS,N_CONST_VARY)
        LOGICAL LOG_HIST(N_CONSTITUENTS)
C
        COMMON /L_CONFIGURE/ LOG_CONSTITUENT, LOG_EPS_TABLE,
     &                       LOG_DRY_DENSITY, LOG_PDF_TYPE, LOG_HIST
C
        COMMON /R_ENVIRON/ T_SOIL, T_WATER, T_VEG
        COMMON /R_TRUNK/ MG_TRUNK, RHO_TRUNK, TRUNK_DIAM
c
        IF(LOG_DRY_DENSITY(1))THEN
            CALL EPSVEG_Mg_2(MG_TRUNK,RHO_TRUNK,T_VEG,EPSR)
        ELSE
            CALL EPSVEG_Mg_1(MG_TRUNK,T_VEG,EPSR)
        ENDIF
C
        RETURN
        END
c
        SUBROUTINE EPSVEG_Mg_1(MG,TdegC,EPSR)
        save
c
        REAL MG, TdegC
        REAL A,B,C,D
        REAL TOPITAU,FNAUGHT,ENAUGHT,EINFIN
        REAL THETA_DEGREES, FREQ_GHZ
        COMPLEX EPSR, E, F, G, CWORK
C
        COMMON /R_SENSOR/ THETA_DEGREES, FREQ_GHZ
c
        a = 1.7 - 0.74*mg + 6.16*mg*mg
        b = mg*(0.55*mg - 0.076)
        c = 4.64*mg*mg/(7.36*mg*mg +1.0)
        d = 22.74
c
        topitau = 1.1109e-10 + TdegC*(-3.824e-12 + TdegC*
     &                (6.938e-14 - TdegC*5.096e-16))
        fnaught = 1.0/(topitau*1.0e9)
        enaught = 88.045 + TdegC*(-0.4147 + TdegC*(6.295e-4 + 
     &                 TdegC*1.075e-5))    
        einfin = 4.9
c
        e = cmplx(1.0,(freq_ghz/fnaught))
        f = cmplx(0.0,(d/freq_ghz))
c        cwork = csqrt(cmplx(0.0,freq_ghz/0.18))
        cwork = (sqrt(2.)/2.)*cmplx(1.0,1.0)*sqrt(freq_ghz/0.18)
        g = 1.0 + cwork
c
        epsr = a + b*(4.9 + (enaught-einfin)/e - f) + c*(2.9 + 55.0/g)
c
c
        return
        end
c
        SUBROUTINE EPSVEG_Mg_2(Mg,RHO,TdegC,EPSR)
        save
c
        real Mg,RHO,Mv
        real a,b,c,d
        REAL THETA_DEGREES, FREQ_GHZ
        real TdegC  
        real topitau,fnaught,enaught,einfin
        complex epsr
        complex e,f,g
        complex cwork
c
        COMMON /R_SENSOR/ THETA_DEGREES, FREQ_GHZ
c
        Mv = Mg*RHO/(1.0-Mg*(1.0-RHO))
C
        a = 1.7 + 3.20*MV + 6.5*MV*MV
        b = MV*(0.82*MV + 0.166)
        c = 31.4*MV*MV/(59.5*MV*MV +1.0)
        d = 22.74
c
        topitau = 1.1109e-10 + TdegC*(-3.824e-12 + TdegC*
     &                (6.938e-14 - TdegC*5.096e-16))
        fnaught = 1.0/(topitau*1.0e9)
        enaught = 88.045 + TdegC*(-0.4147 + TdegC*(6.295e-4 + 
     &                 TdegC*1.075e-5))    
        einfin = 4.9
c
        e = cmplx(1.0,(freq_ghz/fnaught))
        f = cmplx(0.0,(d/freq_ghz))
c        cwork = csqrt(cmplx(0.0,freq_ghz/0.18))
        cwork = (sqrt(2.)/2.)*cmplx(1.0,1.0)*sqrt(freq_ghz/0.18)
        g = 1.0 + cwork
c
        epsr = a + b*(4.9 + (enaught-einfin)/e - f) + c*(2.9 + 55.0/g)
c
c
        return
        end
c
        SUBROUTINE EPS_SOIL(EPSR,T_degC_s)
        save
c
        real T_degC_s
c        real epsp,epspp
        REAL THETA_DEGREES, FREQ_GHZ
        REAL MV_SOIL, SAND, CLAY
C
        complex epsr
        complex a0(9),a1(9),a2(9),b0(9),b1(9),b2(9)
        complex c0(9),c1(9),c2(9)
        complex aa0,aa1,aa2,bb0,bb1,bb2,cc0,cc1,cc2
C
        COMMON /R_SOIL/ MV_SOIL, SAND, CLAY
        COMMON /R_SENSOR/ THETA_DEGREES, FREQ_GHZ
C
        if (T_degC_s.gt.0.0) then
c
            a0(1) = cmplx(2.862,0.356)
            a1(1) = cmplx(-0.012,-0.003)
            a2(1) = cmplx(0.001,-0.008)
c
            b0(1) = cmplx(3.803,5.507)
            b1(1) = cmplx(0.462,0.044)
            b2(1) = cmplx(-0.341,-0.002)
c
            c0(1) = cmplx(119.006,17.753)
            c1(1) = cmplx(-0.500,-0.313)
            c2(1) = cmplx(0.633,0.206)
c
            a0(2) = cmplx(2.927,0.004)
            a1(2) = cmplx(-0.012,0.001)
            a2(2) = cmplx(-0.001,0.002)
c         
            b0(2) = cmplx(5.505,0.951)
            b1(2) = cmplx(0.371,0.005)
            b2(2) = cmplx(0.062,-0.010)
c         
            c0(2) = cmplx(114.826,16.759)
            c1(2) = cmplx(-0.389,0.192)
            c2(2) = cmplx(-0.547,0.290)
c
            a0(3) = cmplx(1.993,-0.123)
            a1(3) = cmplx(0.002,0.002)
            a2(3) = cmplx(0.015,0.003)
c         
            b0(3) = cmplx(38.086,7.502)
            b1(3) = cmplx(-0.176,-0.058)
            b2(3) = cmplx(-0.633,-0.116)
c         
            c0(3) = cmplx(10.72,2.942)
            c1(3) = cmplx(1.256,0.452)
            c2(3) = cmplx(1.522,0.543)
c
            a0(4) = cmplx(1.997,-0.201)
            a1(4) = cmplx(0.002,0.003)
            a2(4) = cmplx(0.018,0.003)
c       
            b0(4) = cmplx(25.579,11.266)
            b1(4) = cmplx(-0.17,-0.085)
            b2(4) = cmplx(-0.412,-0.155)
c         
            c0(4) = cmplx(39.793,0.194)
            c1(4) = cmplx(0.723,0.584)
            c2(4) = cmplx(0.941,0.581)
c
            a0(5) = cmplx(2.502,-0.070)
            a1(5) = cmplx(-0.003,0.000)
            a2(5) = cmplx(-0.003,0.001)
c                     
            b0(5) = cmplx(10.101,6.620)
            b1(5) = cmplx(0.221,0.015)
            b2(5) = cmplx(-0.004,-0.081)
c         
            c0(5) = cmplx(77.482,21.578)
            c1(5) = cmplx(-0.061,0.293)
            c2(5) = cmplx(-0.135,0.332)
c
            a0(6) = cmplx(2.200,-0.142)
            a1(6) = cmplx(-0.001,0.001)
            a2(6) = cmplx(0.012,0.003)
c                     
            b0(6) = cmplx(26.473,11.868)
            b1(6) = cmplx(0.013,-0.059)
            b2(6) = cmplx(-0.523,-0.225)
c                     
            c0(6) = cmplx(34.333,7.817)
            c1(6) = cmplx(0.284,0.570)
            c2(6) = cmplx(1.062,0.801)
c                     
            a0(7) = cmplx(2.301,-0.096)
            a1(7) = cmplx(0.001,0.001)
            a2(7) = cmplx(0.009,0.002)
c                     
            b0(7) = cmplx(17.918,8.583)
            b1(7) = cmplx(0.084,-0.005)
            b2(7) = cmplx(-0.282,-0.153)
c         
            c0(7) = cmplx(50.149,28.707)
            c1(7) = cmplx(0.012,0.297)
            c2(7) = cmplx(0.387,0.357)
c
            a0(8) = cmplx(2.237,-0.027)
            a1(8) = cmplx(0.002,-0.001)
            a2(8) = cmplx(0.009,0.003)
c                     
            b0(8) = cmplx(15.505,6.179)
            b1(8) = cmplx(0.076,0.074)
            b2(8) = cmplx(-0.217,-0.086)
c                     
            c0(8) = cmplx(48.260,34.126)
            c1(8) = cmplx(0.168,0.143)
            c2(8) = cmplx(0.289,0.206)
c                     
            a0(9) = cmplx(1.912,-0.071)
            a1(9) = cmplx(0.007,0.000)
            a2(9) = cmplx(0.021,0.003)
c                     
            b0(9) = cmplx(29.123,6.938)
            b1(9) = cmplx(-0.190,0.029)
            b2(9) = cmplx(-0.545,-0.128)
c         
            c0(9) = cmplx(6.960,29.945)
            c1(9) = cmplx(0.822,0.275)
            c2(9) = cmplx(1.195,0.377)
c  
            if (freq_ghz.le.1.4)then                 
                aa0 = a0(1)
                aa1 = a1(1)
                aa2 = a2(1)
c
                bb0 = b0(1)
                bb1 = b1(1)
                bb2 = b2(1)
c               
                cc0 = c0(1)
                cc1 = c1(1)
                cc2 = c2(1)
            endif
c            
            if ((freq_ghz.gt.1.4).and.(freq_ghz.le.4.))then                 
                aa0 = a0(1)+((a0(1)-a0(2))/(1.4-4.0))*(freq_ghz-1.4)
                aa1 = a1(1)+((a1(1)-a1(2))/(1.4-4.0))*(freq_ghz-1.4)
                aa2 = a2(1)+((a2(1)-a2(2))/(1.4-4.0))*(freq_ghz-1.4)
c                 
                bb0 = b0(1)+((b0(1)-b0(2))/(1.4-4.0))*(freq_ghz-1.4)
                bb1 = b1(1)+((b1(1)-b1(2))/(1.4-4.0))*(freq_ghz-1.4)
                bb2 = b2(1)+((b2(1)-b2(2))/(1.4-4.0))*(freq_ghz-1.4)
c               
                cc0 = c0(1)+((c0(1)-c0(2))/(1.4-4.0))*(freq_ghz-1.4)
                cc1 = c1(1)+((c1(1)-c1(2))/(1.4-4.0))*(freq_ghz-1.4)
                cc2 = c2(1)+((c2(1)-c2(2))/(1.4-4.0))*(freq_ghz-1.4)
            endif
c  
            if ((freq_ghz.gt.4.0).and.(freq_ghz.le.6.))then                 
                aa0 = a0(2)+((a0(2)-a0(3))/(4.0-6.0))*(freq_ghz-4.0)
                aa1 = a1(2)+((a1(2)-a1(3))/(4.0-6.0))*(freq_ghz-4.0)
                aa2 = a2(2)+((a2(2)-a2(3))/(4.0-6.0))*(freq_ghz-4.0)
c                                                         
                bb0 = b0(2)+((b0(2)-b0(3))/(4.0-6.0))*(freq_ghz-4.0)
                bb1 = b1(2)+((b1(2)-b1(3))/(4.0-6.0))*(freq_ghz-4.0)
                bb2 = b2(2)+((b2(2)-b2(3))/(4.0-6.0))*(freq_ghz-4.0)
c               
                cc0 = c0(2)+((c0(2)-c0(3))/(4.0-6.0))*(freq_ghz-4.0)
                cc1 = c1(2)+((c1(2)-c1(3))/(4.0-6.0))*(freq_ghz-4.0)
                cc2 = c2(2)+((c2(2)-c2(3))/(4.0-6.0))*(freq_ghz-4.0)
            endif
c  
            if ((freq_ghz.gt.6.0).and.(freq_ghz.le.8.))then                 
                aa0 = a0(3)+((a0(3)-a0(4))/(6.-8.))*(freq_ghz-6.0)
                aa1 = a1(3)+((a1(3)-a1(4))/(6.-8.))*(freq_ghz-6.0)
                aa2 = a2(3)+((a2(3)-a2(4))/(6.-8.))*(freq_ghz-6.0)
c                                                       
                bb0 = b0(3)+((b0(3)-b0(4))/(6.-8.))*(freq_ghz-6.0)
                bb1 = b1(3)+((b1(3)-b1(4))/(6.-8.))*(freq_ghz-6.0)
                bb2 = b2(3)+((b2(3)-b2(4))/(6.-8.))*(freq_ghz-6.0)
c                                                       
                cc0 = c0(3)+((c0(3)-c0(4))/(6.-8.))*(freq_ghz-6.0)
                cc1 = c1(3)+((c1(3)-c1(4))/(6.-8.))*(freq_ghz-6.0)
                cc2 = c2(3)+((c2(3)-c2(4))/(6.-8.))*(freq_ghz-6.0)
            endif
c
            if ((freq_ghz.gt.8.0).and.(freq_ghz.le.10.))then                
                aa0 = a0(4)+((a0(4)-a0(5))/(8.-10.))*(freq_ghz-8.0)
                aa1 = a1(4)+((a1(4)-a1(5))/(8.-10.))*(freq_ghz-8.0)
                aa2 = a2(4)+((a2(4)-a2(5))/(8.-10.))*(freq_ghz-8.0)
c                 
                bb0 = b0(4)+((b0(4)-b0(5))/(8.-10.))*(freq_ghz-8.0)
                bb1 = b1(4)+((b1(4)-b1(5))/(8.-10.))*(freq_ghz-8.0)
                bb2 = b2(4)+((b2(4)-b2(5))/(8.-10.))*(freq_ghz-8.0)
c                                                            
                cc0 = c0(4)+((c0(4)-c0(5))/(8.-10.))*(freq_ghz-8.0)
                cc1 = c1(4)+((c1(4)-c1(5))/(8.-10.))*(freq_ghz-8.0)
                cc2 = c2(4)+((c2(4)-c2(5))/(8.-10.))*(freq_ghz-8.0)
            endif   
c
            if ((freq_ghz.gt.10.0).and.(freq_ghz.le.12.))then                
                aa0 = a0(5)+((a0(5)-a0(6))/(10.-12.))*(freq_ghz-10.0)
                aa1 = a1(5)+((a1(5)-a1(6))/(10.-12.))*(freq_ghz-10.0)
                aa2 = a2(5)+((a2(5)-a2(6))/(10.-12.))*(freq_ghz-10.0)
c                 
                bb0 = b0(5)+((b0(5)-b0(6))/(10.-12.))*(freq_ghz-10.0)
                bb1 = b1(5)+((b1(5)-b1(6))/(10.-12.))*(freq_ghz-10.0)
                bb2 = b2(5)+((b2(5)-b2(6))/(10.-12.))*(freq_ghz-10.0)
c                                                              
                cc0 = c0(5)+((c0(5)-c0(6))/(10.-12.))*(freq_ghz-10.0)
                cc1 = c1(5)+((c1(5)-c1(6))/(10.-12.))*(freq_ghz-10.0)
                cc2 = c2(5)+((c2(5)-c2(6))/(10.-12.))*(freq_ghz-10.0)
            endif   
c       
            if ((freq_ghz.gt.12.0).and.(freq_ghz.le.14.))then                
                aa0 = a0(6)+((a0(6)-a0(7))/(12.-14.))*(freq_ghz-12.0)
                aa1 = a1(6)+((a1(6)-a1(7))/(12.-14.))*(freq_ghz-12.0)
                aa2 = a2(6)+((a2(6)-a2(7))/(12.-14.))*(freq_ghz-12.0)
c                     
                bb0 = b0(6)+((b0(6)-b0(7))/(12.-14.))*(freq_ghz-12.0)
                bb1 = b1(6)+((b1(6)-b1(7))/(12.-14.))*(freq_ghz-12.0)
                bb2 = b2(6)+((b2(6)-b2(7))/(12.-14.))*(freq_ghz-12.0)
c                                                                
                cc0 = c0(6)+((c0(6)-c0(7))/(12.-14.))*(freq_ghz-12.0)
                cc1 = c1(6)+((c1(6)-c1(7))/(12.-14.))*(freq_ghz-12.0)
                cc2 = c2(6)+((c2(6)-c2(7))/(12.-14.))*(freq_ghz-12.0)
            endif   
c       
            if ((freq_ghz.gt.14.0).and.(freq_ghz.le.16.))then                
                aa0 = a0(7)+((a0(7)-a0(8))/(14.-16.))*(freq_ghz-14.0)
                aa1 = a1(7)+((a1(7)-a1(8))/(14.-16.))*(freq_ghz-14.0)
                aa2 = a2(7)+((a2(7)-a2(8))/(14.-16.))*(freq_ghz-14.0)
c                     
                bb0 = b0(7)+((b0(7)-b0(8))/(14.-16.))*(freq_ghz-14.0)
                bb1 = b1(7)+((b1(7)-b1(8))/(14.-16.))*(freq_ghz-14.0)
                bb2 = b2(7)+((b2(7)-b2(8))/(14.-16.))*(freq_ghz-14.0)
c                                                               
                cc0 = c0(7)+((c0(7)-c0(8))/(14.-16.))*(freq_ghz-14.0)
                cc1 = c1(7)+((c1(7)-c1(8))/(14.-16.))*(freq_ghz-14.0)
                cc2 = c2(7)+((c2(7)-c2(8))/(14.-16.))*(freq_ghz-14.0)
            endif   
c       
            if ((freq_ghz.gt.16.0).and.(freq_ghz.le.18.))then                
                aa0 = a0(8)+((a0(8)-a0(9))/(16.-18.))*(freq_ghz-16.0)
                aa1 = a1(8)+((a1(8)-a1(9))/(16.-18.))*(freq_ghz-16.0)
                aa2 = a2(8)+((a2(8)-a2(9))/(16.-18.))*(freq_ghz-16.0)
c                                           
                bb0 = b0(8)+((b0(8)-b0(9))/(16.-18.))*(freq_ghz-16.0)
                bb1 = b1(8)+((b1(8)-b1(9))/(16.-18.))*(freq_ghz-16.0)
                bb2 = b2(8)+((b2(8)-b2(9))/(16.-18.))*(freq_ghz-16.0)
c                                                               
                cc0 = c0(8)+((c0(8)-c0(9))/(16.-18.))*(freq_ghz-16.0)
                cc1 = c1(8)+((c1(8)-c1(9))/(16.-18.))*(freq_ghz-16.0)
                cc2 = c2(8)+((c2(8)-c2(9))/(16.-18.))*(freq_ghz-16.0)
            endif   
c       
            if (freq_ghz.gt.18.)then                 
                aa0 = a0(9)
                aa1 = a1(9)
                aa2 = a2(9)
c       
                bb0 = b0(9)
                bb1 = b1(9)
                bb2 = b2(9)
c                   
                cc0 = c0(9)
                cc1 = c1(9)
                cc2 = c2(9)
            endif
c
c
         epsr = (aa0+aa1*Sand+aa2*Clay)+(bb0+bb1*Sand+bb2*Clay)*Mv_SOIL
     &             + (cc0+cc1*Sand+cc2*Clay)*Mv_SOIL*Mv_SOIL
c
        else
c
            CALL WRITE_ERROR(40)
c
C            epsp = 3.0 + Mg_S*(20.0+(2./3.)*T_degC_s)
C            epspp = (1.0 + T_degC_s/50.)*2.*Mg_s/0.15
c
C            epsr = cmplx(epsp,-epspp)
            stop
c
        endif
c
        return
        end
c
        SUBROUTINE EPS_WATER(EPSR, TdegC)
        save
c
        REAL TdegC  
        REAL SALT
        REAL F,T,S
        REAL THETA_DEGREES, FREQ_GHZ
        REAL ESWI,ESW0,ESW0S,A,B,ESWP,PHI,SIG,TAW0,DELT
        REAL ESWD,TAWS,SIG25
c
        COMPLEX EPSR, ESW
c
        COMMON /R_SENSOR/ THETA_DEGREES, FREQ_GHZ
        COMMON /R_WATER/ SALT
c
        include 'constants.include'
c
        F = FREQ_GHZ
        T = TdegC
        S = SALT
C
       ESWI=4.9
       ESW0=87.134-1.949E-1*T-1.276E-2*(T**2)+2.491E-4*(T**3)
       A=1.0+1.613E-5*T*S-3.656E-3*S+3.21E-5*(S**2)-4.232E-7*(S**3)
       ESW0S=ESW0*A
       TAW0=(1./(2.*PI))*(1.1109E-10-3.824E-12*T+6.938E-14*(T**2)-
     & 5.096E-16*(T**3))
       B=1.0+2.282E-5*T*S-7.638E-4*S-7.76E-6*(S**2)+1.105E-8*(S**3)
       TAWS=TAW0*B
       DELT=25.0-T
       PHI=DELT*(2.033E-2+1.266E-4*DELT+2.464E-6*(DELT**2)-
     & S*(1.849E-5-2.551E-7*DELT+2.551E-8*(DELT**2)))
       SIG25=S*(0.18252-1.4619E-3*S+2.093E-5*(S**2)-1.282E-7*(S**3))
       SIG=SIG25*EXP(-PHI)
       ESWP=ESWI+((ESW0S-ESWI)/(1.+((2.*PI*F*1.E9*TAWS)**2)))
       ESWD=2.*PI*F*1.E9*TAWS*(ESW0S-ESWI)/(1.+((2.*PI*F*1.E9*TAWS)**2))
     + +SIG/(2.*PI*EPS_0*F*1.E9)
       ESW=CMPLX(ESWP,-ESWD)
c
        EPSR = ESW
c
       RETURN
       END
c
