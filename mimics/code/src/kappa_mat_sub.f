        subroutine kappa_mat_sub(m, kappa)  
        save
c
        real kappa(4,4)
        complex m(2,2)
c
        kappa(1,1) = -2.0*real(m(1,1))
        kappa(2,1) =  0.0
        kappa(3,1) = -2.0*real(m(2,1))
        kappa(4,1) =  2.0*aimag(m(2,1))

        kappa(1,2) =  0.0
        kappa(2,2) = -2.0*real(m(2,2))
        kappa(3,2) = -2.0*real(m(1,2))
        kappa(4,2) = -2.0*aimag(m(1,2))

        kappa(1,3) = -real(m(1,2))
        kappa(2,3) = -real(m(2,1))
        kappa(3,3) = -real(m(1,1) + m(2,2))
        kappa(4,3) = -aimag(m(1,1) - m(2,2))

        kappa(1,4) = -aimag(m(1,2))
        kappa(2,4) =  aimag(m(2,1))
        kappa(3,4) =  aimag(m(1,1) - m(2,2))
        kappa(4,4) = -real(m(1,1) + m(2,2))
c
        return
        end
