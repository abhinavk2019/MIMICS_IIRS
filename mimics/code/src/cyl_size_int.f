        SUBROUTINE CYL_SCAT_MAT_SIZE_INT(thetai,phii,thetas,phis,
     &             thetac,phic,PHASE_MAT,SMAT_avg)
        save
c
        INCLUDE 'parameters.include'
c
        integer i, j, k, l
        integer Ns1,Ns2
c
        real STOKES_MAT(4,4), phase_mat(4,4)
        real s1(maxtrs1), s2(maxtrs2),ds1(maxtrs1), ds2(maxtrs2)
        real PDFs1(maxtrs1),PDFs2(maxtrs2)
        real PDF, delta
c
        complex SMAT_avg(2,2)
c
        real thetai,phii,thetas,phis,thetac,phic
        real diamcm,lengthm,k0a
        complex epsr, SMAT(2,2)
        logical lflag
c
        common /rcyldata/ diamcm,lengthm,k0a
        common /ccyldata/ epsr
        common /cylscatmat/ Smat
c
        do 10 i=1,2
            do 9 j=1,2
                smat_avg(j,i) = cmplx(0.0,0.0)
9           continue
10      continue
c
        do 20 i=1,4
            do 19 j=1,4
                phase_mat(j,i) = 0.0
19          continue
20      continue
c
        call pdf_tr_size_setup(s1,s2,PDFs1,PDFs2,ds1,ds2,Ns1,Ns2)
c
        do 150 i=1,Ns1
C            lengthm = s1(i)
            do 140 j=1,Ns2
               lengthm = s1(j)
                diamcm = s2(j)
                delta = ds1(i)*ds2(j)
                PDF = PDFs1(i)*PDFs2(j)

c         print*,' i,j,PDFs1(i),PDFs2(j),PDFs1(i)*PDFs2(j),',
c     &              'ds1(i),ds2(j),ds1(i)*ds2(j) = ',
c     &            i,j,PDFs1(i),PDFs2(j),PDFs1(i)*PDFs2(j),
c     &               ds1(i),ds2(j),ds1(i)*ds2(j) 
c
                CALL RESONANT_CYL_SCAT_MAT(thetai,phii,thetas,phis,
     &             thetac,phic,lflag)
                CALL STOKES_SUB(SMAT,STOKES_MAT)
c
                do 110 k=1,4
                 do 100 l=1,4
                   phase_mat(l,k) = phase_mat(l,k) + 
     &                               PDF*delta*STOKES_MAT(l,k)
100              continue
110             continue
c
                do 130 k=1,2
                 do 120 l=1,2
                   smat_avg(l,k) = smat_avg(l,k) + PDF*delta*SMAT(l,k)
120              continue
130             continue
c
140         continue
150     continue
c
        return
        end
c
