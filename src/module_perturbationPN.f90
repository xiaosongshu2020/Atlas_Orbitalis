
    module module_perturbationPN
    use module_toolKits
    implicit none
    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine post_newton_effect(r_sat, v_sat, F)
    use module_const
    implicit none
    double precision r_sat(3), v_sat(3), F(3)
    double precision c, r_sat_norm, v_sat_norm, r_rd
    double precision rsu(3), vsu(3)
    !--------------------------------------------------
    rsu = r_sat/length_unit
    vsu = v_sat/length_unit*time_unit
    
    c = 299792458d0/(length_unit*1d3)*time_unit
    r_sat_norm = norm2(rsu)
    v_sat_norm = norm2(vsu)
    call dot(3, rsu, vsu, r_rd)

    F = 1d0/c**2/r_sat_norm**3*((4d0/r_sat_norm - v_sat_norm**2)*rsu + &
    & 4d0*r_rd*vsu)
    
    F = F*acc_unit

    end subroutine
    endmodule