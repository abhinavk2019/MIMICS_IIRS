        subroutine cylinf(fhz,dcm,epsr,psi,phipr,
     &                    sumpolte,sumpoltm,sumdepol)

        save
c
c
        integer bessord,rmaxord,cmaxord
        integer ordlim,sumlim
        integer i
c
        real fghz,dcm,phipr,psi
        real spi
        real fhz

c
        complex epsr,sumpoltm,sumpolte,sumdepol
c
c            
        real*8 pi,mu0,eps0
        real*8 rj(500),ry(500)
        real*8 diam
        real*8 ccc,sss
        real*8 wvlth,k0,a0,x0,s0
        real*8 error
        real*8 dphipr,dpsi,repsrr,repsri
c
c 
        real*8 cjr(500),cyr(500)
        real*8 cji(500),cyi(500)
        real*8 x1r,qr,s1r,r1r
        real*8 x1i,qi,s1i,r1i
        real*8 cwork1r
        real*8 cwork1i
        real*8 sumpolmr,sumpoler,sumdpolr
        real*8 sumpolmi,sumpolei,sumdpoli
        real*8 workr,worki
        real*8 mydble
c
c
        common /rbessel/ rj,ry
        common /cbessel/ cjr,cji,cyr,cyi
        common /besparr/ x0,x1r,s0,s1r,r1r,qr
        common /bespari/ x1i,s1i,r1i,qi
        common /clnsumsr/ sumpolmr,sumpoler,sumdpolr
        common /clnsumsi/ sumpolmi,sumpolei,sumdpoli
        common /limit/ sumlim
c             
        pi = 3.14159265358979324d0
        spi = 3.1415926535897932
        dphipr = mydble(phipr)
        dpsi = mydble(psi)
c
        mu0 = 4.d0*pi*1.0d-7
        eps0 = 8.854d-12
        error = 1.0d-12
c
        fghz = fhz/1.0e9
        wvlth = 0.29979d0/mydble(fghz)
        k0 = 2.d0*pi/wvlth
        diam = mydble(dcm)/100.d0 
        a0 = diam/2.d0
c
        ccc = dcos(dpsi)
        sss = dsin(dpsi)
c
        repsrr = mydble(real(epsr))
        if (aimag(epsr).lt.0.) then
            repsri = -mydble(aimag(epsr))
        else
            repsri = mydble(aimag(epsr))
        endif
c  
        bessord = 30
c
        x0 = k0*a0*ccc
        cwork1r = repsrr - sss*sss
        cwork1i = repsri
        call cdiv(1.d0,0.d0,cwork1r,cwork1i,workr,worki)
        qr = (sss/(k0*a0))*(workr - 1.d0/(ccc*ccc))
        qi = (sss/(k0*a0))*(worki)
        call csqrta(cwork1r,cwork1i,workr,worki)
        cwork1r = workr
        cwork1i = worki
        x1r = k0*a0*cwork1r
        x1i = k0*a0*cwork1i
        s0 = 1.d0/ccc
        call cdiv(repsrr,repsri,cwork1r,cwork1i,s1r,s1i)
        call cdiv(1.d0,0.d0,cwork1r,cwork1i,r1r,r1i)
c
c
20      do 30 i=1,500
            rj(i) = 0.d0
            ry(i) = 0.d0
            cjr(i) = 0.d0
            cji(i) = 0.d0
            cyr(i) = 0.d0
            cyi(i) = 0.d0

30      continue
c
        call rbesrec(x0,bessord,rj,ry,rmaxord)
c
        do 50 i=1,500
            if ((rj(i).eq.0.d0).or.(ry(i).eq.0.d0)) go to 55
50      continue
c
55      if ((i-2).lt.bessord) then
          if ((i-2).ge.0) then
            ordlim = i-2
          else
           print*,' *** Only ',i,' real bessel functions computed ***'
           stop
          endif
        else
            ordlim = bessord
        endif          
c        print*,' ordlim = ',ordlim
c
        call cbesrec(x1r,x1i,bessord,cjr,cji,cyr,cyi,cmaxord)
c
        do 60 i=1,500
            if(((dabs(cjr(i)).eq.0.d0).and.
     &         (dabs(cji(i)).eq.0.d0)).or.
     &         ((dabs(cyr(i)).eq.0.d0).and.
     &         (dabs(cyi(i)).eq.0.d0)))go to 65
60      continue
c
65      if ((i-2).lt.ordlim) then
          if ((i-2).ge.0) then
            ordlim = i-2
          else
           print*,' *** Only ',i,
     &         ' complex bessel functions computed ***'
           stop
          endif
        endif
c
        if (ordlim.gt.bessord) then
            ordlim = bessord
        endif
c
c        print*,' ordlim = ',ordlim
        call sumcons(dphipr,error,ordlim,sumlim)
c
        if (sumlim.gt.ordlim) then
            if (ordlim.lt.bessord) then
                print*,'*** Summation of Cn`s does not converge',
     &                  ' after ',ordlim,' iterations. ***'
            else
                bessord = bessord + 10
                go to 20
            endif
        endif
c
c
        sumpoltm = cmplx(sngl(sumpolmr),sngl(sumpolmi))
        sumpolte = cmplx(sngl(sumpoler),sngl(sumpolei))
        sumdepol = cmplx(sngl(sumdpolr),sngl(sumdpoli))
c
        return
        end
c
        subroutine sumcons(dphipr,error,ordlim,sumlim)
c
        save

c
        integer ordlim,sumlim,n
c
        real*8 dphipr,error,mult
        real*8 temp1,temp2
c
        real*8 sumpolmr,sumpoler,sumdpolr
        real*8 sumpolmi,sumpolei,sumdpoli
        real*8 cntmr,cnter,cndr
        real*8 cntmi,cntei,cndi
        real*8 cmultr, cmulti, workr, worki, work1r
        real*8 dfloat
c
        common /clnsumsr/ sumpolmr,sumpoler,sumdpolr
        common /clnsumsi/ sumpolmi,sumpolei,sumdpoli
c
        sumpolmr = 0.d0
        sumpolmi = 0.d0
        sumpoler = 0.d0
        sumpolei = 0.d0
        sumdpolr = 0.d0
        sumdpoli = 0.d0
        n = 0
c
        call cntmd(n,ordlim,sumpolmr,sumpolmi,
     &            sumpoler,sumpolei,sumdpolr,sumdpoli)
c
        mult = 1.d0
        do 500 n=1,ordlim
            mult = -mult
c 
            call cntmd(n,ordlim,cntmr,cntmi,cnter,cntei,cndr,cndi)
c
            workr = 0.d0
            worki = 1.d0*dfloat(n)*dphipr
            call cexpa(workr,worki,cmultr,cmulti)
c
            call cmult(cmultr,cmulti,cntmr,cntmi,workr,worki)
            sumpolmr = sumpolmr + mult*workr
            sumpolmi = sumpolmi + mult*worki
c
            call cmult(cmultr,cmulti,cnter,cntei,workr,worki)
            sumpoler = sumpoler + mult*workr
            sumpolei = sumpolei + mult*worki
c
            call cmult(cmultr,cmulti,cndr,cndi,workr,worki)
            sumdpolr = sumdpolr + mult*workr
            sumdpoli = sumdpoli + mult*worki
c
c
            workr = 0.d0
            worki = -1.d0*dfloat(n)*dphipr
            call cexpa(workr,worki,cmultr,cmulti)
c
            call cmult(cmultr,cmulti,cntmr,cntmi,workr,worki)
            sumpolmr = sumpolmr + mult*workr
            sumpolmi = sumpolmi + mult*worki
c
            call cmult(cmultr,cmulti,cnter,cntei,workr,worki)
            sumpoler = sumpoler + mult*workr
            sumpolei = sumpolei + mult*worki
c
            call cmult(cmultr,cmulti,cndr,cndi,workr,worki)
            sumdpolr = sumdpolr - mult*workr
            sumdpoli = sumdpoli - mult*worki
c
            call cabsa(cntmr,cntmi,workr)
            call cabsa(sumpolmr,sumpolmi,work1r)
            temp1= workr/work1r
            call cabsa(cnter,cntei,workr)
            call cabsa(sumpoler,sumpolei,work1r)
            temp2 = workr/work1r
c            call cabsa(cndr,cndi,workr)
c            call cabsa(sumdpolr,sumdpoli,work1r)
c            temp3 = workr/work1r
c            print*,'temp3 = ',temp3
c
c            WRITE(6,111)N,CNTM,TEMP1,CNTE,TEMP2
c
            if(TEMP1.LT.ERROR)then
                if(TEMP2.LT.ERROR)then
                    sumlim = n
                    return
                endif
            endif
500     continue
c
        sumlim = ordlim + 1
c
        end                     
c
        subroutine cntmd(n,ordlim,cntmr,cntmi,
     &                   cnter,cntei,cndr,cndi)
c
        save

c
        integer n,n1,n2,ordlim
c
        real*8 x0,s0,rj(500),ry(500)
        real*8 rjnp,rynp,pi
        real*8 workr,worki,work1r,work1i,work2r,work2i
c
        real*8 x1r,s1r,r1r,qr,cjr(500),cyr(500)
        real*8 x1i,s1i,r1i,qi,cji(500),cyi(500)
        real*8 cntmr,cnter,cndr
        real*8 cntmi,cntei,cndi
        real*8 cjnpr,chnpr,hnx0r
        real*8 cjnpi,chnpi,hnx0i
        real*8 ct1r,ct2r
        real*8 ct1i,ct2i
c
        real*8 qnr,vnr,pnr,nnr,mnr
        real*8 qni,vni,pni,nni,mni
        real*8 dfloat
c
        common /rbessel/ rj,ry
        common /cbessel/ cjr,cji,cyr,cyi
        common /besparr/ x0,x1r,s0,s1r,r1r,qr
        common /bespari/ x1i,s1i,r1i,qi
c
        pi = 3.1415926535897d0
        cntmr = 0.d0
        cntmi = 0.d0
        cnter = 0.d0
        cntei = 0.d0
c
        n1=n+1
        n2=n+2    
c
        hnx0r = rj(n1)
        hnx0i = ry(n1)
c
        if(n.lt.ordlim)then
            RJNP = dfloat(N)*RJ(N1)/X0-RJ(N2)
            rynp = dfloat(n)*ry(n1)/x0-ry(n2)
            call cdiv(cjr(n1),cji(n1),x1r,x1i,workr,worki)
            cjnpr = dfloat(n)*workr - cjr(n2)
            cjnpi = dfloat(n)*worki - cji(n2)
            chnpr = rjnp           
            chnpi = rynp
        else if(n.eq.ordlim)then
            rjnp = rj(n) - dfloat(n)*rj(n1)/x0
            rynp = ry(n) - dfloat(n)*ry(n1)/x0
            call cdiv(cjr(n1),cji(n1),x1r,x1i,workr,worki)
            cjnpr = cjr(n1) - dfloat(n)*workr
            cjnpi = cji(n1) - dfloat(n)*worki
            chnpr = rjnp
            chnpi = rynp
        else
            print*,'  *** ERROR ***   Requested n is greater than',
     &             ' available order.'
        endif
c
        call cmult(qr,qi,dfloat(n),0.d0,qnr,qni)
c
        call cmult(s1r,s1i,cjnpr,cjnpi,workr,worki)
        vnr = rj(n1)*workr - s0*rjnp*cjr(n1)
        vni = rj(n1)*worki - s0*rjnp*cji(n1)
c
        call cmult(r1r,r1i,hnx0r,hnx0i,workr,worki)
        call cmult(workr,worki,cjnpr,cjnpi,work1r,work1i)
        call cmult(chnpr,chnpi,cjr(n1),cji(n1),workr,worki)
        pnr = work1r - s0*workr
        pni = work1i - s0*worki
c
        call cmult(s1r,s1i,hnx0r,hnx0i,work2r,work2i)
        call cmult(work2r,work2i,cjnpr,cjnpi,work1r,work1i)
        nnr = work1r - s0*workr
        nni = work1i - s0*worki
c
        call cmult(r1r,r1i,cjnpr,cjnpi,workr,worki)
        mnr = rj(n1)*workr - s0*rjnp*cjr(n1)
        mni = rj(n1)*worki - s0*rjnp*cji(n1)
c
        call cmult(qnr,qni,qnr,qni,workr,worki)
        call cmult(workr,worki,hnx0r,hnx0i,work1r,work1i)
        call cmult(work1r,work1i,cjr(n1),cji(n1),workr,worki)
        call cmult(workr,worki,cjr(n1),cji(n1),work1r,work1i)
        ct1r = work1r*rj(n1)
        ct1i = work1i*rj(n1)
c
        call cmult(pnr,pni,nnr,nni,workr,worki)
        call cmult(qnr,qni,hnx0r,hnx0i,work1r,work1i)
        call cmult(work1r,work1i,cjr(n1),cji(n1),work2r,work2i)
        call cmult(work2r,work2i,work2r,work2i,work1r,work1i)
        ct2r = workr - work1r
        ct2i = worki - work1i
c
        call cmult(vnr,vni,pnr,pni,workr,worki)
        work1r = ct1r-workr
        work1i = ct1i-worki
        call cdiv(work1r,work1i,ct2r,ct2i,cntmr,cntmi)
c
        call cmult(mnr,mni,nnr,nni,workr,worki)
        work1r = ct1r-workr
        work1i = ct1i-worki
        call cdiv(work1r,work1i,ct2r,ct2i,cnter,cntei)
C
        if (n.eq.0) then
            cndr = 0.d0
            cndi = 0.d0
        else
            call cmult(cjr(n1),cji(n1),cjr(n1),cji(n1),workr,worki)
            call cdiv(workr,worki,ct2r,ct2i,work1r,work1i)
            call cmult(qnr,qni,work1r,work1i,workr,worki)
            work1r = (2.d0/(pi*x0))*workr*s0
            work1i = (2.d0/(pi*x0))*worki*s0
            call cmult(0.d0,1.d0,work1r,work1i,cndr,cndi)
        endif
c
      RETURN
      END
c
        subroutine rbesrec(x,order,jn,yn,maxord) 
c
        save

c
        integer order,imagx,nsig,l,n
        integer k,kk,maxord
c
c       
        real*8 x
        real*8 dx,magx
        real*8 jn(500),yn(500)
        real*8 work1,work2,work3
        real*8 p(500),test,mtest,magp
        real*8 y(500),mu
        real*8 sum
        real*8 y0
        real*8 dfloat,pi
c
        do 20 k=1,500
            jn(k) = 0.d0
            yn(k) = 0.d0
            p(k) = 0.d0
            y(k) = 0.d0
20      continue
c      
        pi = 3.14159265358979324d0
        nsig = 16
c
        dx = x
c
        magx = dabs(dx)
        imagx = idint(magx)
        l = max0((imagx+1),(order-1))
c
        p(imagx+1) = 0.d0
        p(imagx+2) = 1.d0
c
        do 100 k=imagx+2,500
c
            work3 = 2.d0*dfloat(k-1)/dx
            work2 = p(k)*work3
c
	    IF(k<500) THEN
            p(k+1) = work2 - p(k-1)
            endif
c
            if(k.gt.(l+1)) then
c
                work1 = p(l+1)*p(l+2)
c
                work3 = 2.d0*(10.d0**nsig)
                work2 = work1*work3
                work1 = dsqrt(work2)
c      
                test = work1
c
                mtest = dabs(test)
                magp = dabs(p(k+1))
c
                if(mtest.lt.magp)then
                    n = k
                    go to 105
                endif
            endif
100     continue
c
        print*,' ** No maximum value for N reached **  '
c
105     continue
c
        y(n+1) = 0.d0
        y(n) = 1.d0/p(n+1)
c
        do 200 k=1,n-1
            kk = n - k 
c
            work3 = 2.d0*dfloat(kk)/dx
            work2 = work3*y(kk+1)
            y(kk) = work2 - y(kk+2)
c
200     continue
c
        sum = 0.d0
c
        do 240 k=1,int(n/2)
            sum = sum + y(2*k+1)
240     continue
c
        mu = 2.d0*sum + y(1)
c
        do 300 k=1,n+1
            jn(k) = y(k)/mu
300     continue
c
        call yzeror(dx,y,mu,n,y0)
c
        yn(1) = y0
c
        do 350 k=1,n
c
            work3 = jn(k+1)*yn(k)
            work2 = 2.d0/(pi*dx)
c
            work1 = work3 - work2
c
            yn(k+1) = work1/jn(k)
c
350     continue
c
            maxord = n
c
            return
            end
c
        subroutine yzeror(dx,y,mu,n,y0)
        save
c
        integer n,k,j,i
c
        real*8 pi,euler
        real*8 dx,y(500),y0
        real*8 sum,mult,mu
        real*8 dfloat
        real*8 dlnx
        real*8 work1,work2
c       
        pi = 3.14159265358979324d0
        euler = 0.57721566490153286061d0
c
        sum = 0.d0
        mult = 1.d0
c
        i = 1 + 2*int(n/2)
        do 100 j=1,int(n/2)
            k = i - 2*(j-1)
            mult = -mult
            sum = sum + mult*y(k)/dfloat(int((k-1)/2))
100     continue      
c
        if(mult.gt.0.d0)then
            sum = -sum
        endif
c
        sum = sum*(4.d0/pi)
c
        dlnx = dlog(0.5d0*dx)
c
        work1 = dlnx + euler
        work1 = (2.d0/pi)*work1*y(1)
        work2 = work1 - sum
c
         y0 = work2/mu
c
        return
        end
c
        subroutine cbesrec(dzr,dzi,order,jnr,jni,ynr,yni,maxord) 
        save
c
        integer order,imagz,nsig,l,n
        integer k,kk,maxord
c       
        real*8 dzr,dzi,magz
        real*8 jnr(500),jni(500),ynr(500),yni(500)
        real*8 work1r,work1i,work2r,work2i,work3r,work3i
        real*8 rp(500),ip(500),testr,testi,mtest,magp
        real*8 ry(500),iy(500),mur,mui
        real*8 sumr,sumi,rfact,ifact
        real*8 y0r,y0i
        real*8 dfloat,pi
c
        do 20 k=1,500
            jnr(k) = 0.d0
            jni(k) = 0.d0
            ynr(k) = 0.d0
            yni(k) = 0.d0
            rp(k) = 0.d0
            ip(k) = 0.d0
            ry(k) = 0.d0
            iy(k) = 0.d0 
20      continue
c
c       
        pi = 3.14159265358979324d0
        nsig = 16
c
        call cabsa(dzr,dzi,magz)
c
        imagz = idint(magz)
        l = max0((imagz+1),(order-1))
c
        rp(imagz+1) = 0.d0
        ip(imagz+1) = 0.d0
        rp(imagz+2) = 1.d0
        ip(imagz+2) = 0.d0
c
        do 100 k=imagz+2,500
c
            work1r = 2.d0*dfloat(k-1)
            work1i = 0.d0
            work2r = dzr
            work2i = dzi
c
            call cdiv(work1r,work1i,work2r,work2i,work3r,work3i)
c
            work1r = rp(k)
            work1i = ip(k)
c
            call cmult(work3r,work3i,work1r,work1i,work2r,work2i)
c

	    if(k<500) then
            rp(k+1) = work2r - rp(k-1)
            ip(k+1) = work2i - ip(k-1)
            endif
c
            if(k.gt.(l+1)) then
c
                work2r = rp(l+1)
                work2i = ip(l+1)
                work3r = rp(l+2)
                work3i = ip(l+2)
c
                call cmult(work2r,work2i,work3r,work3i,work1r,work1i)
c
                work3r = 2.d0*(10.d0**nsig)
                work2r = work1r*work3r
                work2i = work1i*work3r
c
                call csqrta(work2r,work2i,work1r,work1i)
c      
                testr = work1r
                testi = work1i
c
                call cabsa(testr,testi,mtest)
                call cabsa(rp(k+1),ip(k+1),magp)
c
                if(mtest.lt.magp)then
                    n = k
                    go to 105
                endif
            endif
100     continue
c
        print*,' ** No maximum value for N reached **  '
c
105     continue
c
        ry(n+1) = 0.d0
        iy(n+1) = 0.d0
c
        work1r = rp(n+1)
        work1i = ip(n+1)
c
        call cdiv(1.d0,0.d0,work1r,work1i,ry(n),iy(n))
c
        do 200 k=1,n-1
            kk = n - k 
c
            work1r = 2.d0*dfloat(kk)
            work1i = 0.d0
            work2r = dzr
            work2i = dzi
c
            call cdiv(work1r,work1i,work2r,work2i,work3r,work3i)
c
            work1r = ry(kk+1)
            work1i = iy(kk+1)
c
            call cmult(work3r,work3i,work1r,work1i,work2r,work2i)
c
            ry(kk) = work2r - ry(kk+2)
            iy(kk) = work2i - iy(kk+2)
c
200     continue
c
        if (dzi.eq.0.d0) then
c
            sumr = 0.d0
            sumi = 0.d0
c
            do 240 k=1,int(n/2)
                sumr = sumr + ry(2*k+1)
                sumi = sumi + iy(2*k+1)
240         continue
c
            mur = 2.d0*sumr + ry(1)
            mui = 2.d0*sumi + iy(1)
c
        endif
c
c
        if (dzi.gt.0.d0) then
c
            sumr = 0.d0
            sumi = 0.d0
            rfact = 1.d0
            ifact = 0.d0
c
            do 260 k=1,n
c
                work1r = rfact
                work1i = ifact
c
                call cmult(work1r,work1i,0.d0,-1.d0,work3r,work3i)
c
                rfact = work3r
                ifact = work3i
c
                work1r = ry(k+1)
                work1i = iy(k+1)
c
                call cmult(work3r,work3i,work1r,work1i,work2r,work2i)
c
                sumr = sumr + work2r
                sumi = sumi + work2i
c
260         continue
c
            sumr = 2.d0*sumr + ry(1)
            sumi = 2.d0*sumi + iy(1)
c
            work1r = dzr
            work1i = dzi
c
            call cmult(0.d0,1.d0,work1r,work1i,work2r,work2i)
c
            call cexpa(work2r,work2i,work1r,work1i)
c
            call cmult(work1r,work1i,sumr,sumi,work2r,work2i)
c
            mur = work2r
            mui = work2i
c
        endif
c
c
        if (dzi.lt.0.d0) then
c
            sumr = 0.d0
            sumi = 0.d0
            rfact = 1.d0
            ifact = 0.d0
c
            do 280 k=1,n
c
                work1r = rfact
                work1i = ifact
c
                call cmult(work1r,work1i,0.d0,1.d0,work3r,work3i)
c
                rfact = work3r
                ifact = work3i
c
                work1r = ry(k+1)
                work1i = iy(k+1)
c
                call cmult(work3r,work3i,work1r,work1i,work2r,work2i)
c
                sumr = sumr + work2r
                sumi = sumi + work2i
c
280         continue
c
            sumr = 2.d0*sumr + ry(1)
            sumi = 2.d0*sumi + iy(1)
c
            work1r = dzr
            work1i = dzi
c
            call cmult(0.d0,-1.d0,work1r,work1i,work2r,work2i)
c
            call cexpa(work2r,work2i,work1r,work1i)
c
            call cmult(work1r,work1i,sumr,sumi,work2r,work2i)
c
            mur = work2r
            mui = work2i
c
        endif
c
        do 300 k=1,n+1
c
            work1r = ry(k)
            work1i = iy(k)
            work2r = mur
            work2i = mui
c
            call cdiv(work1r,work1i,work2r,work2i,work3r,work3i)
c
            jnr(k) = work3r
            jni(k) = work3i
c
300     continue
c
        call yzeroc(dzr,dzi,ry,iy,mur,mui,n,y0r,y0i)
c
        ynr(1) = y0r
        yni(1) = y0i
c
        do 350 k=1,n
            work1r = jnr(k+1)
            work1i = jni(k+1)
            work2r = ynr(k)
            work2i = yni(k)
c
            call cmult(work1r,work1i,work2r,work2i,work3r,work3i)
c
            work1r = 2.d0/pi
            work1i = 0.d0
c
            call cdiv(work1r,work1i,dzr,dzi,work2r,work2i)
c
            work1r = work3r - work2r
            work1i = work3i - work2i
c
            work2r = jnr(k)
            work2i = jni(k)
c
            call cdiv(work1r,work1i,work2r,work2i,work3r,work3i)
c
            ynr(k+1) = work3r
            yni(k+1) = work3i
c
350     continue
c
            maxord = n
c
            return
            end
c
        subroutine yzeroc(dzr,dzi,ry,iy,mur,mui,n,y0r,y0i)
        save
c
        integer n,k,j,i
c
        real*8 pi,euler
        real*8 dzr,dzi,ry(500),iy(500),y0r,y0i
        real*8 sumr,sumi,mult,mur,mui
        real*8 dfloat 
        real*8 dlnzr,dlnzi
        real*8 work1r,work1i,work2r,work2i
c       
        pi = 3.14159265358979324d0
        euler = 0.57721566490153286061d0
c
        sumr = 0.d0
        sumi = 0.d0
        mult = 1.d0
c
        i = 1 + 2*int(n/2)
        do 100 j=1,int(n/2)
            k = i - 2*(j-1)
            mult = -mult
            sumr = sumr + mult*ry(k)/dfloat(int((k-1)/2))
            sumi = sumi + mult*iy(k)/dfloat(int((k-1)/2))
100     continue      
c
        if(mult.gt.0.d0)then
            sumr = -sumr
            sumi = -sumi
        endif
c
        sumr = sumr*(4.d0/pi)
        sumi = sumi*(4.d0/pi)
c
        work1r = 0.5d0*dzr
        work1i = 0.5d0*dzi                  
c
        call clna(work1r,work1i,dlnzr,dlnzi)
c
        work1r = dlnzr + euler
        work1i = dlnzi
c
        call cmult(work1r,work1i,ry(1),iy(1),work2r,work2i)
c
        work1r = (2.d0/pi)*work2r
        work1i = (2.d0/pi)*work2i
c
        work2r = work1r - sumr
        work2i = work1i - sumi
c
        call cdiv(work2r,work2i,mur,mui,work1r,work1i)
c
         y0r = work1r
         y0i = work1i
c
        return
        end
c
