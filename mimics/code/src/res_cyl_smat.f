        SUBROUTINE resonant_cyl_scat_mat(thetai,phii,thetas,phis,
     &                 thetac,phic,lflag)
        save
c
        real pi,tol,tol2,dum,dum2,dum3,sumcheck,spsii,spsis
        real k0,diamcm,lengthm,fhz,thresh
        real thetai,phii,thetac,phic,thetas,phis
        real xhat(3),yhat(3),zhat(3),xhatl(3),yhatl(3),zhatl(3)
        real khati(3),khats(3),rhatc(3),phihatc(3),zhatc(3)
        real hhati(3),vhati(3),hhats(3),vhats(3),hhatc(3),vhatc(3)
        real hhatcs(3),vhatcs(3),hhcinf(3),vhcinf(3),hhinf(3),vhinf(3)
        real psii,psis,phip,khspinf(3),khsinf(3)
        real khatip(3),khatsp(3),khipxy(3),khspxy(3)
        real vhidvhc,hhidvhc,vsdzhc,vsdrhc,vsdphc,hsdzhc,hsdrhc,hsdphc
        real vhidhhc,hhidhhc
        real vhsdvhsc, vhsdhhsc, hhsdvhsc, hhsdhhsc
c
        complex sumpolte,sumpoltm,sumdepol,j,q
        complex Smat(2,2)
        complex epsr

        complex SMATP(2,2), CDUMMAT(2,2)
        real T(2,2),S(2,2)
c
        logical lflag
c
        common /cylsum/ sumpolte,sumpoltm,sumdepol
        common /rcyldata/ diamcm,lengthm,k0
        common /ccyldata/ epsr
        common /cylscatmat/ Smat
c
        pi = 3.141592654
c        pi = 3.141593
        j = cmplx(0.0,1.0)
c        tol = 1.0e-2
        tol = 0.035
        tol2 = 0.01
        thresh = 0.1
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
        xhatl(1) = cos(thetac)*cos(phic)
        xhatl(2) = cos(thetac)*sin(phic)
        xhatl(3) = -sin(thetac)
        call vctsmth(xhatl,xhatl)
c
        yhatl(1) = -sin(phic)
        yhatl(2) = cos(phic)
        yhatl(3) = 0.0
        call vctsmth(yhatl,yhatl)
c
        zhatl(1) = sin(thetac)*cos(phic)
        zhatl(2) = sin(thetac)*sin(phic)
        zhatl(3) = cos(thetac)
        call vctsmth(zhatl,zhatl)
c
        khati(1) = sin(thetai)*cos(phii)
        khati(2) = sin(thetai)*sin(phii)
        khati(3) = cos(thetai)
c
        khats(1) = sin(thetas)*cos(phis)
        khats(2) = sin(thetas)*sin(phis)
        khats(3) = cos(thetas)
c        call vctsmth(khats,khats)
c        print*,'khats = ',khats
c
cc        if(zhatl(3).lt.0.0)then
c        call dot(zhatl,khati,dum)
c        if(dum.gt.0.0)then
c            zhatl(1) = -zhatl(1)
c            zhatl(2) = -zhatl(2)
c            zhatl(3) = -zhatl(3)
c                xhatl(1) = -xhatl(1)
c                xhatl(2) = -xhatl(2)
c                xhatl(3) = -xhatl(3)
c        endif
c
c---------- set polarization vectors ----------------------------------c
c------ incident direction ---------
c
        hhati(1) = -sin(phii)
        hhati(2) = cos(phii)
        hhati(3) = 0.0
c
        call vctsmth(hhati,hhati)
c
        vhati(1) = cos(phii)*cos(thetai)
        vhati(2) = sin(phii)*cos(thetai)
        vhati(3) = -sin(thetai)
c
c        call cross(hhati,khati,vhati)
c        call vctsmth(vhati,vhati)
          
        hhats(1) = -sin(phis)
        hhats(2) = cos(phis)
        hhats(3) = 0.0
c
        call vctsmth(hhats,hhats)
c
        vhats(1) = cos(phis)*cos(thetas)
        vhats(2) = sin(phis)*cos(thetas)
        vhats(3) = -sin(thetas)
c
c        call cross(hhats,khats,vhats)
c        call vctsmth(vhats,vhats)

        zhatc(1) = zhatl(1)
        zhatc(2) = zhatl(2)
        zhatc(3) = zhatl(3)
c
        call dot(zhatc,khati,dum)
c
        rhatc(1) = -(khati(1) - dum*zhatc(1))
        rhatc(2) = -(khati(2) - dum*zhatc(2))
        rhatc(3) = -(khati(3) - dum*zhatc(3))
        call vctsmth(rhatc,rhatc)
c
        call cross(zhatc,rhatc,phihatc)
        call vctsmth(phihatc,phihatc)
c
        call cross(zhatl,khati,hhatc)
        call vctsmth(hhatc,hhatc)
c
        call cross(hhatc,khati,vhatc)
        call vctsmth(vhatc,vhatc)
c
        sumcheck = hhatc(1)*hhatc(1)+hhatc(2)*hhatc(2)+hhatc(3)*hhatc(3)
c
        if (sumcheck.eq.0.0) then
            hhatc(1) = hhati(1)
            hhatc(2) = hhati(2)
            hhatc(3) = hhati(3)
            vhatc(1) = vhati(1)
            vhatc(2) = vhati(2)
            vhatc(3) = vhati(3)
        endif
c
        call cross(zhatl,khats,hhatcs)
        call vctsmth(hhatcs,hhatcs)
c
        call cross(hhatcs,khats,vhatcs)
        call vctsmth(vhatcs,vhatcs)
c
        sumcheck = hhatcs(1)*hhatcs(1) + hhatcs(2)*hhatcs(2)
     &             + hhatcs(3)*hhatcs(3)
c
        if (sumcheck.eq.0.0) then
            hhatcs(1) = hhats(1)
            hhatcs(2) = hhats(2)
            hhatcs(3) = hhats(3)
            vhatcs(1) = vhats(1)
            vhatcs(2) = vhats(2)
            vhatcs(3) = vhats(3)
        endif
c
        call dot(khati,zhatl,dum)
        if(dum.gt.1.0)dum = 1.0
        if(dum.lt.-1.0)dum = -1.0
        dum2 = acos(-dum)
c
           psii = pi/2. - dum2
c
            call dot(khats,zhatl,dum)
            if(dum.gt.1.0)dum = 1.0
            if(dum.lt.-1.0)dum = -1.0
            dum3 = acos(dum)
c
                    psis = (pi/2.0 - dum3)
c
            spsii = sin(psii)
            spsis = sin(psis)
            dum = spsii+spsis
            if(abs((abs(psii)-(pi/2.))).lt.tol)then
                Smat(1,1) = cmplx(0.0,0.0)
                Smat(1,2) = cmplx(0.0,0.0)
                Smat(2,1) = cmplx(0.0,0.0)
                Smat(2,2) = cmplx(0.0,0.0)
                lflag = .false.
c
                RETURN
c
            else if(abs(dum).lt.tol2)then
                dum2 = 1
                lflag = .true.
            else
                dum3 = k0*dum*lengthm/2.0
                dum2 = sin(dum3)/dum3
                if(abs(dum2).lt.thresh)then
                    Smat(1,1) = cmplx(0.0,0.0)
                    Smat(1,2) = cmplx(0.0,0.0)
                    Smat(2,1) = cmplx(0.0,0.0)
                    Smat(2,2) = cmplx(0.0,0.0)
                    lflag = .false.
c
                    RETURN
c
                else
                    lflag = .true.
                endif
            endif
c
            call dot(khati,xhatl,khatip(1))
            call dot(khati,yhatl,khatip(2))
            call dot(khati,zhatl,khatip(3))
            call vctsmth(khatip,khatip)
            khipxy(1) = -khatip(1)
            khipxy(2) = -khatip(2)
            khipxy(3) = 0.0
            call vctsmth(khipxy,khipxy)
c
            call dot(khats,xhatl,khatsp(1))
            call dot(khats,yhatl,khatsp(2))
            call dot(khats,zhatl,khatsp(3))
            call vctsmth(khatsp,khatsp)
            khspxy(1) = khatsp(1)
            khspxy(2) = khatsp(2)
            khspxy(3) = 0.0
            call vctsmth(khspxy,khspxy)
c
            call dot(khspxy,khipxy,dum)
            if(dum.le.-1.0) dum = -1.0
            if(dum.ge.1.0) dum = 1.0
            phip = acos(dum)

            dum = khipxy(1)*khspxy(2) - khipxy(2)*khspxy(1)
            if(dum.lt.0.0) phip = -phip
c
            khspinf(1) = khspxy(1)*cos(psii)
            khspinf(2) = khspxy(2)*cos(psii)
c
            call dot(khati,zhatl,dum)
            if(dum.gt.0.0)then
                khspinf(3) = sin(psii)
            else
                khspinf(3) = -sin(psii)
            endif
c
            call vctsmth(khspinf,khspinf)
c
            khsinf(1) = khspinf(1)*xhatl(1) + khspinf(2)*yhatl(1) +
     &                  khspinf(3)*zhatl(1)
            khsinf(2) = khspinf(1)*xhatl(2) + khspinf(2)*yhatl(2) +
     &                  khspinf(3)*zhatl(2)
            khsinf(3) = khspinf(1)*xhatl(3) + khspinf(2)*yhatl(3) +
     &                  khspinf(3)*zhatl(3)
            call vctsmth(khsinf,khsinf)
c
        zhatc(1) = zhatl(1)
        zhatc(2) = zhatl(2)
        zhatc(3) = zhatl(3)
c
        call dot(zhatc,khsinf,dum)
c
        rhatc(1) = khsinf(1) - dum*zhatc(1)
        rhatc(2) = khsinf(2) - dum*zhatc(2)
        rhatc(3) = khsinf(3) - dum*zhatc(3)
        call vctsmth(rhatc,rhatc)
c
        call cross(zhatc,rhatc,phihatc)
        call vctsmth(phihatc,phihatc)
c
        hhcinf(1) = hhatcs(1) 
        hhcinf(2) = hhatcs(2) 
        hhcinf(3) = hhatcs(3)
c
        call cross(hhcinf,khsinf,vhcinf)
c
        call dot(vhats,vhatcs,dum)
        call dot(vhats,hhatcs,dum2)
        if(dum.ge.1.0) dum=1.0
        if(dum2.ge.1.0) dum2=1.0
        if(dum.le.-1.0) dum=-1.0
        if(dum2.le.-1.0) dum2=-1.0
        dum = acos(dum)
        dum2 = acos(dum2)
c
        vhinf(1) = vhcinf(1)*cos(dum) + hhcinf(1)*cos(dum2)
        vhinf(2) = vhcinf(2)*cos(dum) + hhcinf(2)*cos(dum2)
        vhinf(3) = vhcinf(3)*cos(dum) + hhcinf(3)*cos(dum2)
c
        call vctsmth(vhinf,vhinf)
c
        call dot(hhats,vhatcs,dum)
        call dot(hhats,hhatcs,dum2)
        dum = acos(dum)
        dum2 = acos(dum2)
c
        hhinf(1) = vhcinf(1)*cos(dum) + hhcinf(1)*cos(dum2)
        hhinf(2) = vhcinf(2)*cos(dum) + hhcinf(2)*cos(dum2)
        hhinf(3) = vhcinf(3)*cos(dum) + hhcinf(3)*cos(dum2)
c
        call vctsmth(hhinf,hhinf)
c
        call dot(vhati,vhatc,vhidvhc)
        call dot(hhati,vhatc,hhidvhc)
        call dot(vhati,hhatc,vhidhhc)
        call dot(hhati,hhatc,hhidhhc)
c
        call dot(zhatc,vhinf,vsdzhc)
        call dot(rhatc,vhinf,vsdrhc)
        call dot(phihatc,vhinf,vsdphc)
c
        call dot(phihatc,hhinf,hsdphc)
        call dot(zhatc,hhinf,hsdzhc)
        call dot(rhatc,hhinf,hsdrhc)
c
        call dot(vhinf,vhcinf,vhsdvhsc)
        call dot(vhinf,hhcinf,vhsdhhsc)
        call dot(hhinf,vhcinf,hhsdvhsc)
        call dot(hhinf,hhcinf,hhsdhhsc)
c
         fhz = (k0*2.9979e8)/(2.0*pi)
       call cylinf(fhz,diamcm,epsr,psii,phip,sumpolte,sumpoltm,sumdepol)
c
        call dot(vhati,zhat,dum)
        call dot(khati,zhatc,dum2)
c
            spsii = sin(psii)
            spsis = sin(psis)
            dum = spsii+spsis
            if(abs(dum).lt.tol2)then
                dum2 = 1
            else
                dum3 = k0*dum*lengthm/2.0
                dum2 = sin(dum3)/dum3
            endif
            Q = ((lengthm*cos(psis))/(j*pi*cos(psii)))*dum2
c
            SMATP(1,1) = sumpoltm
c
            SMATP(2,1) = -sumdepol*cmplx(0.0,1.0)
            SMATP(1,2) = sumdepol*cmplx(0.0,1.0)
c
            SMATP(2,2) = sumpolte
c
            T(1,1) = vhidvhc 
            T(1,2) = hhidvhc
            T(2,1) = vhidhhc
            T(2,2) = hhidhhc
c
            S(1,1) = vhsdvhsc 
            S(1,2) = vhsdhhsc
            S(2,1) = hhsdvhsc
            S(2,2) = hhsdhhsc
c
            CDUMMAT(1,1) = SMATP(1,1)*T(1,1) + SMATP(1,2)*T(2,1)
            CDUMMAT(2,1) = SMATP(2,1)*T(1,1) + SMATP(2,2)*T(2,1)
            CDUMMAT(1,2) = SMATP(1,1)*T(1,2) + SMATP(1,2)*T(2,2)
            CDUMMAT(2,2) = SMATP(2,1)*T(1,2) + SMATP(2,2)*T(2,2)
c
            Smat(1,1) = Q*(S(1,1)*CDUMMAT(1,1) + S(1,2)*CDUMMAT(2,1))
            Smat(2,1) = Q*(S(2,1)*CDUMMAT(1,1) + S(2,2)*CDUMMAT(2,1))
            Smat(1,2) = Q*(S(1,1)*CDUMMAT(1,2) + S(1,2)*CDUMMAT(2,2))
            Smat(2,2) = Q*(S(2,1)*CDUMMAT(1,2) + S(2,2)*CDUMMAT(2,2))
c
        RETURN
        END
