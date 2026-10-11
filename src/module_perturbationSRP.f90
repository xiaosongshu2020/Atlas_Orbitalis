
    module module_perturbationSRP
    use module_planets
    use module_const
    use module_time
    use module_elements
    
    integer, public :: shadowModel = 1 ! 1 - 圆柱形地影；2 - 锥形地影；
    
    double precision :: rho_s = 0.3166d-17 ! Solar radiative constant (normalized unit)
    double precision :: rho_s_SI = 4.5605d-6 ! N/m^2
    
    ! amr(m^2/kg) -> k(normalized acc unit)
    double PRECISION :: amr2k = 1d0/(length_unit*1d3)**2*mass_unit*0.3169d-17
    double PRECISION :: k2amr = 1d0/(1d0/(length_unit*1d3)**2*mass_unit*0.3169d-17)
    
    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This subroutine is written to calculate the perturbation force of Solar radiative pressure
    !   F_epsilon = - kappa*S/m*rho_s*r_sun_unit, If the satellite in not in the shadow of the Earth
    !last modified: 2021/3/23    Shawn Hu
    subroutine solar_pressure_perturbation(UTC, r_sat, reflectivity, amrSRP, &
        &   SRP_model, ECOM, F)

    implicit none
    double precision, intent(in) :: UTC, r_sat(6), reflectivity, amrSRP
    integer, intent(in) :: SRP_Model
    double precision, intent(in) :: ECOM(9)
    double precision, intent(out) :: F(3)
    logical :: inshadow = .false.

    double precision k, r_unit(3)
    double precision r_sun(3), solarPV(6), kappa, TT
    double precision rsu(3), amrSRPu
    ! ECOM model variables
    
    double precision kpl(6), r, E, fo, u, r_sat_u(6)
    double precision d_unit(3), y_unit(3), b_unit(3)
    double precision D0, Dc, Ds, Y0, Yc, Ys, B0, Bc, Bs
    double precision f_d(3), f_y(3), f_b(3)
    !--------------------------------------------------------
    
    rsu = r_sat(1:3)/length_unit
    r_sat_u(1:3) = rsu
    r_sat_u(4:6) = r_sat(4:6)/length_unit*time_unit
    amrSRPu = amrSRP/(length_unit*1d3)**2*mass_unit
    
    call UTC2TT(UTC, TT)
    call planetsPV(TT, nEarth, nSun, SolarPV)
    r_sun = solarPV(1:3)/length_unit
    
    call in_shadow(rsu, r_sun, inshadow)
    
    ! If sat is in the shadow, F_SRP = 0
    if (inshadow) then
        F = 0d0
        return
    endif
    
    ! calculate SRP
    if (SRP_Model == 1) then ! canon ball
        
        kappa = 1d0 + reflectivity
        k = kappa*amrSRPu*rho_s*(1.49591492d8/6378.1363d0/norm2(rsu - r_sun))**2
        r_unit = (rsu - r_sun)/norm2(rsu - r_sun)
        F = k*r_unit
        
    elseif (SRP_Model == 2) then ! ECOM (9 paras)
        
        ! get u
        call cts2kpl(r_sat_u, kpl)
        call kepler2rEf(kpl, r, E, fo)
        u = kpl(5) + fo
        
        ! get d, y, b directions
        d_unit = (rsu - r_sun)/norm2(rsu - r_sun)
        call cross_3d(r_sun, rsu, y_unit)
        y_unit = y_unit/norm2(y_unit)
        call cross_3d(d_unit, y_unit, b_unit)
        
        ! get force
        D0 = ECOM(1)
        Y0 = ECOM(2)
        B0 = ECOM(3)
        Bc = ECOM(4)
        Bs = ECOM(5)
        Dc = ECOM(6)
        Ds = ECOM(7)
        Yc = ECOM(8)
        Ys = ECOM(9)
        
        f_d = (D0 + Dc*dcos(u) + Ds*dsin(u)) * d_unit
        f_y = (Y0 + Yc*dcos(u) + Ys*dsin(u)) * y_unit
        f_b = (B0 + Bc*dcos(u) + Bs*dsin(u)) * b_unit
        F = f_d + f_y + f_b
        
    end if
        
    F = F*acc_unit
    
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine in_shadow(r_sat, r_sun, flag_shadow)
    implicit none
    logical flag_shadow
    double precision r_sat(3), r_sun(3), r_sun_unit(3), r_sat_unit(3)
    double precision sin_phi, cos_phi
    !---------------------------------------------------------------
    r_sun_unit = r_sun/norm2(r_sun)
    r_sat_unit = r_sat/norm2(r_sat)
    call dot(3, r_sun_unit, r_sat_unit, cos_phi)
    sin_phi = dsqrt(1d0 - cos_phi**2)
    if ((cos_phi < 0) .and. (sin_phi < 1d0/norm2(r_sat))) then
        flag_shadow = .true.
    else
        flag_shadow = .false.
    end if
    end subroutine

    
    endmodule