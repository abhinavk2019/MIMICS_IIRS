        SUBROUTINE leaf_rlgh_scat_mat(thetai,phii,thetas,phis,
     &                      thetad,phid)
        save
c
        include 'parameters.include'
c
        real dum,dum2,sumcheck
        real k0, t_lf, d_lf
        real thetai,phii,thetad,phid,thetas,phis
        real xhat(3),yhat(3),zhat(3),xhatl(3),yhatl(3),zhatl(3)
        real khati(3),khats(3)
        real hhati(3),vhati(3),hhats(3),vhats(3)
        real vhsdxhl,vhsdyhl,vhsdzhl,hhsdxhl,hhsdyhl,hhsdzhl
        real vhidxhl,vhidyhl,vhidzhl,hhidxhl,hhidyhl,hhidzhl
c
        real Ae,Be,Ce
c
        complex epsr,cdum
        COMPLEX EPSILONR(N_EPS), EPSILONRC(N_EPS)
        complex Smat(2,2)

        real Ac, Ab, Aa, Vo
        complex Vd, cduma, cdumb, cdumc
c
        real FREQ_HERTZ, WAVELENGTH
        real MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
c
        COMMON /RADAR_FREQ/ FREQ_HERTZ, WAVELENGTH, k0
        COMMON /R_LEAF/ MG_LEAF,RHO_LEAF,LEAF_DENS,LEAF_DIAM,LEAF_TAU
        COMMON /C_DIELECTRIC/ EPSILONR, EPSILONRC
        COMMON /LEAF_SCAT_MAT/ Smat
c
        include 'constants.include'
c
         t_lf = LEAF_TAU
         d_lf = LEAF_DIAM
c
         epsr =  EPSILONRC(4)
c
        xhat(1) = 1.0
        xhat(2) = 0.0
        xhat(3) = 0.0
c
        yhat(1) = 0.0
        yhat(2) = 1.0
        yhat(3) = 0.0
c
        zhat(1) = 0.0
        zhat(2) = 0.0
        zhat(3) = 1.0
c
        khati(1) = sin(thetai)*cos(phii)
        khati(2) = sin(thetai)*sin(phii)
        khati(3) = cos(thetai)
c
        call vctsmth(khati,khati)
c
        khats(1) = sin(thetas)*cos(phis)
        khats(2) = sin(thetas)*sin(phis)
        khats(3) = cos(thetas)
c
        call vctsmth(khats,khats)
c
        hhati(1) = -sin(phii)
        hhati(2) = cos(phii)
        hhati(3) = 0.0
c
        call vctsmth(hhati,hhati)
c
        call cross(hhati,khati,vhati)
        call vctsmth(vhati,vhati)
c         
        hhats(1) = -sin(phis)
        hhats(2) = cos(phis)
        hhats(3) = 0.0
c
        call vctsmth(hhats,hhats)
c
        call cross(hhats,khats,vhats)
        call vctsmth(vhats,vhats)
c
        xhatl(1) = cos(thetad)*cos(phid)
        xhatl(2) = cos(thetad)*sin(phid)
        xhatl(3) = -sin(thetad)
c
        call vctsmth(xhatl,xhatl)
c
        yhatl(1) = -sin(phid)
        yhatl(2) = cos(phid)
        yhatl(3) = 0.0
c
        call vctsmth(yhatl,yhatl)
c
        zhatl(1) = sin(thetad)*cos(phid)
        zhatl(2) = sin(thetad)*sin(phid)
        zhatl(3) = cos(thetad)
c
        call vctsmth(zhatl,zhatl)
c
        call dot(vhats,xhatl,vhsdxhl)
        call dot(vhats,yhatl,vhsdyhl)
        call dot(vhats,zhatl,vhsdzhl)

        call dot(hhats,xhatl,hhsdxhl)
        call dot(hhats,yhatl,hhsdyhl)
        call dot(hhats,zhatl,hhsdzhl)

        call dot(vhati,xhatl,vhidxhl)
        call dot(vhati,yhatl,vhidyhl)
        call dot(vhati,zhatl,vhidzhl)

        call dot(hhati,xhatl,hhidxhl)
        call dot(hhati,yhatl,hhidyhl)
        call dot(hhati,zhatl,hhidzhl)

c
            DUM = (3.0/2.0)**(1.0/3.0)
c
            Ce = (t_lf/200.)*dum
            Ae = (d_lf/200.)*dum
            Be = Ae
c
            dum = Ae*Ae-Ce*Ce
            dum2 = (sqrt(dum))/Ce
            Ac = (2.0/(dum**(3./2.)))*(dum2-atan(dum2))
            Ab = (2./(Ae*Be*Ce) - Ac)/2.0
            Aa = Ab
c
            Vo = 4.*pi*Ae*Be*Ce/3.
            Vd = (Ae*Be*Ce/2.)*(epsr-1.)
c
        cdum = (k0*k0/(4.*pi))*Vo*(epsr-1.)
        cduma = 1.+ Vd*Aa
        cdumb = cduma
        cdumc = 1.+ Vd*Ac
c
          Smat(1,1) = cdum*(vhsdxhl*vhidxhl/cduma
     &                 + vhsdyhl*vhidyhl/cdumb + vhsdzhl*vhidzhl/cdumc)

          Smat(2,1) = cdum*(hhsdxhl*vhidxhl/cduma
     &                 + hhsdyhl*vhidyhl/cdumb + hhsdzhl*vhidzhl/cdumc)

          Smat(1,2) = cdum*(vhsdxhl*hhidxhl/cduma
     &                 + vhsdyhl*hhidyhl/cdumb + vhsdzhl*hhidzhl/cdumc)

          Smat(2,2) = cdum*(hhsdxhl*hhidxhl/cduma
     &                 + hhsdyhl*hhidyhl/cdumb + hhsdzhl*hhidzhl/cdumc)
c
        RETURN
        END
