        subroutine  snow_layer(eps_snow,theta_loc,two_way_tau)
        save
c
        real two_way_tau, kappa_s, theta_loc
        complex eps_snow

        REAL FREQ_HERTZ, WAVELENGTH, k0
        REAL T_SNOW
c
        COMMON /RADAR_FREQ/ FREQ_HERTZ, WAVELENGTH, k0
        COMMON /R_SNOW/ T_SNOW

        kappa_s = 2.0*k0*abs(aimag(csqrt(eps_snow)))
        two_way_tau = exp(-2.0*kappa_s*t_snow/cos(theta_loc))
c
        return
        end
