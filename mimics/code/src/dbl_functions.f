        double precision function mydble(a)
        save
        real a
        character*30 chara

        write(chara,10)a
 10     format(e30.17)

        read(chara,10)mydble

        return
        end
c
        double precision function dfloat(n)
        save
        integer n
        real a
        character*30 chara
c
        a = float(n)
        write(chara,10)a
 10     format(e30.17)

        read(chara,10)dfloat

        return
        end
c
