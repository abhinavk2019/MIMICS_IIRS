        SUBROUTINE long_thin_cyl_scatmat(thetai,phii,thetas,phis,
     &                      thetac,phic,diamcm,lengthcm,k0,epsr,Smat)
        save
c
        real pi
        real k0,diamcm,radiuscm,lengthcm
        real thetai,phii,thetac,phic,thetas,phis
        real xhat(3),yhat(3),zhat(3),xhatl(3),yhatl(3),zhatl(3)
        real khati(3),khats(3)
        real hhati(3),vhati(3),hhats(3),vhats(3)
        real khatip(3),khatsp(3)
        real vhatip(3),hhatip(3),vhatsp(3),hhatsp(3)
        real xhatp(3),yhatp(3),zhatp(3)
c
        real area, U, sinfact, mfactor, khsdzhl, cosbeta
c
        complex epsr, Pxx, Pyy, Pzz, pdotv(3), pdoth(3)
        complex Smat(2,2)
        complex Svv,Svh,Shv,Shh, Sh(3), Sv(3)
c
        pi = 3.141592654
        radiuscm = diamcm/2.0
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
        xhatl(1) = cos(thetac)*cos(phic)
        xhatl(2) = cos(thetac)*sin(phic)
        xhatl(3) = -sin(thetac)
c
        call vctsmth(xhatl,xhatl)
c
        yhatl(1) = -sin(phic)
        yhatl(2) = cos(phic)
        yhatl(3) = 0.0
c
        call vctsmth(yhatl,yhatl)
c
        zhatl(1) = sin(thetac)*cos(phic)
        zhatl(2) = sin(thetac)*sin(phic)
        zhatl(3) = cos(thetac)
c
        call vctsmth(zhatl,zhatl)
c
        xhatp(1) = xhatl(1)
        xhatp(2) = yhatl(1)
        xhatp(3) = zhatl(1)

        yhatp(1) = xhatl(2)
        yhatp(2) = yhatl(2)
        yhatp(3) = zhatl(2)

        zhatp(1) = xhatl(3)
        zhatp(2) = yhatl(3)
        zhatp(3) = zhatl(3)
c
        call dot(vhati,xhatl,vhatip(1))
        call dot(vhati,yhatl,vhatip(2))
        call dot(vhati,zhatl,vhatip(3))

        call dot(vhats,xhatl,vhatsp(1))
        call dot(vhats,yhatl,vhatsp(2))
        call dot(vhats,zhatl,vhatsp(3))

        call dot(hhati,xhatl,hhatip(1))
        call dot(hhati,yhatl,hhatip(2))
        call dot(hhati,zhatl,hhatip(3))

        call dot(hhats,xhatl,hhatsp(1))
        call dot(hhats,yhatl,hhatsp(2))
        call dot(hhats,zhatl,hhatsp(3))
c
        call cross(hhatip,vhatip,khatip)
        call cross(hhatsp,vhatsp,khatsp)
c
        call dot(khats,zhatl,khsdzhl)
        call dot(khati,zhatl,cosbeta)
c
        area = pi*radiuscm*radiuscm/(100.*100.)

        Pzz = area*(epsr - 1.0)
        Pxx = 2.*Pzz/(epsr + 1.)
        Pyy = Pxx 

        U = 0.5*k0*(lengthcm/100.)*(khsdzhl-cosbeta)
        if(abs(U).gt. 0.01) then
            sinfact = sin(u)/u
        else
            sinfact = 1.0
        endif

        mfactor = -sinfact*0.25*k0*k0*lengthcm/(100.*4.*pi)
c
        pdotv(1) = Pxx*vhatip(1) 
        pdotv(2) = Pyy*vhatip(2) 
        pdotv(3) = Pzz*vhatip(3)

        Sv(1) = pdotv(1)*(-(khatsp(2)*khatsp(2) + khatsp(3)*khatsp(3)))
     &        + pdotv(2)*(khatsp(1)*khatsp(2))
     &        + pdotv(3)*(khatsp(1)*khatsp(3))

        Sv(2) = pdotv(1)*(khatsp(1)*khatsp(2))
     &        + pdotv(2)*(-(khatsp(3)*khatsp(3) + khatsp(1)*khatsp(1)))
     &        + pdotv(3)*(khatsp(3)*khatsp(2))

        Sv(3) = pdotv(1)*(khatsp(1)*khatsp(3))
     &        + pdotv(2)*(khatsp(2)*khatsp(3))
     &        + pdotv(3)*(-(khatsp(1)*khatsp(1) + khatsp(2)*khatsp(2)))
c
        pdoth(1) = Pxx*hhatip(1) 
        pdoth(2) = Pyy*hhatip(2) 
        pdoth(3) = Pzz*hhatip(3)

        Sh(1) = pdoth(1)*(-(khatsp(2)*khatsp(2) + khatsp(3)*khatsp(3)))
     &        + pdoth(2)*(khatsp(1)*khatsp(2))
     &        + pdoth(3)*(khatsp(1)*khatsp(3))

        Sh(2) = pdoth(1)*(khatsp(1)*khatsp(2))
     &        + pdoth(2)*(-(khatsp(3)*khatsp(3) + khatsp(1)*khatsp(1)))
     &        + pdoth(3)*(khatsp(3)*khatsp(2))

        Sh(3) = pdoth(1)*(khatsp(1)*khatsp(3))
     &        + pdoth(2)*(khatsp(2)*khatsp(3))
     &        + pdoth(3)*(-(khatsp(1)*khatsp(1) + khatsp(2)*khatsp(2)))
c
        Svv = Sv(1)*vhatsp(1) + Sv(2)*vhatsp(2) + Sv(3)*vhatsp(3)
        Shv = Sv(1)*hhatsp(1) + Sv(2)*hhatsp(2) + Sv(3)*hhatsp(3) 
        Svh = Sh(1)*vhatsp(1) + Sh(2)*vhatsp(2) + Sh(3)*vhatsp(3)
        Shh = Sh(1)*hhatsp(1) + Sh(2)*hhatsp(2) + Sh(3)*hhatsp(3)


        Smat(1,1) = Svv*mfactor
        Smat(2,1) = Shv*mfactor
        Smat(1,2) = Svh*mfactor
        Smat(2,2) = Shh*mfactor

c
        RETURN
        END
