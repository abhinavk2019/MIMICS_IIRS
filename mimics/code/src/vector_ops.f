        subroutine dot(v1,v2,dprod)
        save
c                                         
        real v1(3),v2(3)
        real dprod
c
        dprod = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
c
        return
        end
c
        subroutine cross(v1,v2,cprod)
        save
c
        real v1(3),v2(3),cprod(3)
c
        cprod(1) = v1(2)*v2(3) - v1(3)*v2(2)
        cprod(2) = v1(3)*v2(1) - v1(1)*v2(3)
        cprod(3) = v1(1)*v2(2) - v1(2)*v2(1)
c
        return
        end
c
        subroutine vctsmth(inv,outv)
        save
c
        integer i
        real inv(3), outv(3), mag
c
        mag = sqrt(inv(1)*inv(1)+inv(2)*inv(2)+inv(3)*inv(3))
c
        if(mag.ne.0.0)then
            do 5 i=1,3
                inv(i) = inv(i)/mag
5           continue
            call checkmag(inv,inv)
            mag = sqrt(inv(1)*inv(1)+inv(2)*inv(2)+inv(3)*inv(3))
            do 10 i=1,3
                outv(i) = inv(i)/mag
10          continue
        else
            do 20 i=1,3
                outv(i) = inv(i)
20          continue
        endif
c
        return
        end
c
        subroutine checkmag(inv,outv)
        save
c
        integer i
        real inv(3), outv(3), tol
c                       
        tol = 1.0e-5
c
        do 10 i=1,3
            if(abs(inv(i)).lt.tol) then
                outv(i) = 0.0
            else
                outv(i) = inv(i)
            endif
10      continue
c
        return
        end
c
        subroutine crosscmplx(v1,v2,cprod)
        save
c
        real v1(3)
        complex v2(3),cprod(3)
c
        cprod(1) = v1(2)*v2(3) - v1(3)*v2(2)
        cprod(2) = v1(3)*v2(1) - v1(1)*v2(3)
        cprod(3) = v1(1)*v2(2) - v1(2)*v2(1)
c
        return
        end
