
    module module_perturbationThirdBody
    use module_planets
    use module_const    
    use module_time
    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine third_body_perturbation(UTC, rs, Moon, MoonTide, Sun, SunTide, &
                                       Mercury, Venus, Mars, Jupiter, Saturn, Uranus, Pluto, F)
    implicit none

    double precision, intent(in) :: UTC, rs(3)
    integer, intent(in) :: Moon, MoonTide, Sun, SunTide
    integer, intent(in) :: Mercury, Venus, Mars, Jupiter, Saturn, Uranus, Pluto
    double precision, intent(out) :: F(3)
    double precision TT
    double precision LunarPV(6), SolarPV(6), PlanetPV(6)
    double precision HM
    double precision r_lunar(3), r_sun(3), r_planet(3), r_lunar_norm, r_sun_norm, r_planet_norm
    double precision delta_lunar(3), delta_sun(3), delta_planet(3), delta_lunar_norm, delta_sun_norm, delta_planet_norm
    double precision mass_lunar, mass_sun, mass_planet, F_moon(3), F_sun(3), F_planet(3), F_lunar_tide(3), F_solar_tide(3)
    double precision rsu(3)
    !-----------------------------------------------------------------------------
    F = 0d0
    !get TT
    call UTC2TT(UTC, TT)
    
    rsu = rs/length_unit
    
    ! For the lunar gravity
    if (Moon == 1) then
        mass_lunar = 1.2300037226627106D-2
        
        call planetsPV(TT, nEarth, nMoon, LunarPV)
        r_lunar = LunarPV(1:3)/length_unit
        
        r_lunar_norm = norm2(r_lunar)
        delta_lunar=rsu - r_lunar
        delta_lunar_norm = norm2(delta_lunar)
        F_moon=-mass_lunar*(delta_lunar/delta_lunar_norm**3+r_lunar/r_lunar_norm**3)
        F = F + F_moon
        ! For lunar tidal deformation
        if(MoonTide == 1) then
            !call get_F_tide(rsu, r_lunar, mass_lunar, F_lunar_tide)
            call get_F_tide_Liulin(rsu, r_lunar, mass_lunar, F_lunar_tide)
            F = F + F_lunar_tide
        endif
    endif
    
    ! For solar gravity
    if(Sun == 1) then
        mass_sun = 332946.03864607110D0

        call planetsPV(TT, nEarth, nSun, SolarPV)
        r_sun = solarPV(1:3)/length_unit
        
        !print*, solarPV(1:3)/length_unit
        
        r_sun_norm = norm2(r_sun)
        delta_sun=rsu - r_sun 
        delta_sun_norm = norm2(delta_sun)
        F_sun=-mass_sun*(delta_sun/delta_sun_norm**3+r_sun/r_sun_norm**3)
        F = F + F_sun
        ! For solar tide
        if (SunTide == 1) then
            call get_F_tide_Liulin(rsu, r_sun, mass_sun, F_solar_tide)
            F = F + F_solar_tide
        endif
    endif
    
    ! For Mercury gravity
    if (Mercury == 1) then
        mass_planet = mercuryEarthMassRatio
        
        call planetsPV(TT, nEarth, nMercury, PlanetPV)
        r_planet = PlanetPV(1:3)/length_unit
        
        r_planet_norm = norm2(r_planet)
        delta_planet = rsu - r_planet
        delta_planet_norm = norm2(delta_planet)
        F_planet = -mass_planet*(delta_planet/delta_planet_norm**3 + r_planet/r_planet_norm**3)
        F = F + F_planet
    endif
    
    ! For Venus gravity
    if (Venus == 1) then
        mass_planet = venusEarthMassRatio
        
        call planetsPV(TT, nEarth, nVenus, PlanetPV)
        r_planet = PlanetPV(1:3)/length_unit
        
        r_planet_norm = norm2(r_planet)
        delta_planet = rsu - r_planet
        delta_planet_norm = norm2(delta_planet)
        F_planet = -mass_planet*(delta_planet/delta_planet_norm**3 + r_planet/r_planet_norm**3)
        F = F + F_planet
    endif
    
    ! For Mars gravity
    if (Mars == 1) then
        mass_planet = marsEarthMassRatio
        
        call planetsPV(TT, nEarth, nMars, PlanetPV)
        r_planet = PlanetPV(1:3)/length_unit
        
        r_planet_norm = norm2(r_planet)
        delta_planet = rsu - r_planet
        delta_planet_norm = norm2(delta_planet)
        F_planet = -mass_planet*(delta_planet/delta_planet_norm**3 + r_planet/r_planet_norm**3)
        F = F + F_planet
    endif
    
    ! For Jupiter gravity
    if (Jupiter == 1) then
        mass_planet = jupiterEarthMassRatio
        
        call planetsPV(TT, nEarth, nJupiter, PlanetPV)
        r_planet = PlanetPV(1:3)/length_unit
        
        r_planet_norm = norm2(r_planet)
        delta_planet = rsu - r_planet
        delta_planet_norm = norm2(delta_planet)
        F_planet = -mass_planet*(delta_planet/delta_planet_norm**3 + r_planet/r_planet_norm**3)
        F = F + F_planet
    endif
    
    ! For Saturn gravity
    if (Saturn == 1) then
        mass_planet = saturnEarthMassRatio
        
        call planetsPV(TT, nEarth, nSaturn, PlanetPV)
        r_planet = PlanetPV(1:3)/length_unit
        
        r_planet_norm = norm2(r_planet)
        delta_planet = rsu - r_planet
        delta_planet_norm = norm2(delta_planet)
        F_planet = -mass_planet*(delta_planet/delta_planet_norm**3 + r_planet/r_planet_norm**3)
        F = F + F_planet
    endif
    
    ! For Uranus gravity
    if (Uranus == 1) then
        mass_planet = uranusEarthMassRatio
        
        call planetsPV(TT, nEarth, nUranus, PlanetPV)
        r_planet = PlanetPV(1:3)/length_unit
        
        r_planet_norm = norm2(r_planet)
        delta_planet = rsu - r_planet
        delta_planet_norm = norm2(delta_planet)
        F_planet = -mass_planet*(delta_planet/delta_planet_norm**3 + r_planet/r_planet_norm**3)
        F = F + F_planet
    endif
    
    ! For Pluto gravity
    if (Pluto == 1) then
        mass_planet = plutoEarthMassRatio
        
        call planetsPV(TT, nEarth, nPluto, PlanetPV)
        r_planet = PlanetPV(1:3)/length_unit
        
        r_planet_norm = norm2(r_planet)
        delta_planet = rsu - r_planet
        delta_planet_norm = norm2(delta_planet)
        F_planet = -mass_planet*(delta_planet/delta_planet_norm**3 + r_planet/r_planet_norm**3)
        F = F + F_planet
    endif
    
    F = F*acc_unit
    
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_F_tide(r_sat, r_body, m_body, F)
    implicit none
    double precision r_sat(3), r_body(3), F(3)
    double precision r_sat_norm, r_body_norm
    double precision sin_phi, cos_phi, cos_lambda, sin_lambda
    double precision sin_phi_p, cos_phi_p, cos_lambda_p, sin_lambda_p
    double precision Moon, nu, beta, m_body, lambda_p, lambda_star
    double precision P2, P2p, P21, P21p, P22, P22p
    double precision P2b, P2pb, P21b, P21pb, P22b, P22pb
    double precision k(3), G(3), F20(3), F21(3), F22(3)
    !-------------------------------------------------------

    r_sat_norm = norm2(r_sat)
    r_body_norm = norm2(r_body)
    
    sin_phi = r_sat(3)/r_sat_norm
    cos_phi = norm2(r_sat(1:2))/r_sat_norm
    cos_lambda = r_sat(1)/norm2(r_sat(1:2))
    sin_lambda = r_sat(2)/norm2(r_sat(1:2))

    sin_phi_p = r_body(3)/r_body_norm
    cos_phi_p = norm2(r_body(1:2))/r_body_norm
    cos_lambda_p = r_body(1)/norm2(r_body(1:2))
    sin_lambda_p = r_body(2)/norm2(r_body(1:2))

    Moon = 0.299d0
    nu = 0.0436d0 
    beta = m_body/r_body_norm**3

    call get_angle(sin_lambda_p, cos_lambda_p, lambda_p)

    lambda_star = lambda_p + nu

    P2 = 3d0/2d0*sin_phi**2 - 0.5d0
    P2p = 3d0*sin_phi
    P21 = 1.5d0*2d0*sin_phi*cos_phi
    P21p = 3d0*(cos_phi**2 - sin_phi**2)/cos_phi
    P22 = 3d0*cos_phi**2
    P22p = -6d0*sin_phi
    P2b = 3d0/2d0*sin_phi_p**2 - 0.5d0
    P2pb = 3d0*sin_phi_p
    P21b = 1.5d0*2d0*sin_phi_p*cos_phi_p
    P21pb = 3d0*(cos_phi_p**2 - sin_phi_p**2)/cos_phi_p
    P22b = 3d0*cos_phi_p**2
    P22pb = -6d0*sin_phi_p
    
    k = (/0d0, 0d0, 1d0/)
    G = (/-sin_lambda, cos_lambda, 0d0/)

    F20 = -Moon*beta*P2b*1d0/r_sat_norm**5*((3d0*P2 + sin_phi*P2p)*r_sat - r_sat_norm*P2p*k)

    F21 = -Moon*beta*2d0/6d0*P21b*(1d0/r_sat_norm**5*((3d0*P21 + sin_phi*P21p)*r_sat - r_sat_norm*P21p*k)&
    &*(cos_lambda*dcos(lambda_star) + sin_lambda*dsin(lambda_star)) + 1d0/norm2(r_sat(1:2))*1d0/r_sat_norm**3&
    &*(P21*(sin_lambda*dcos(lambda_star) - dsin(lambda_star)*cos_lambda))*G)

    F22 = -Moon*beta*2d0/24d0*P22b*(1d0/r_sat_norm**5*((3d0*P22 + sin_phi*P22p)*r_sat - r_sat_norm*P22p*k)*&
    &((cos_lambda**2 - sin_lambda**2)*(dcos(lambda_star)**2 - dsin(lambda_star)**2) + 2d0*sin_lambda*cos_lambda*&
    &2d0*dsin(lambda_star)*dcos(lambda_star)) + 2d0/norm2(r_sat(1:2))*1d0/r_sat_norm**3*(P22*(2d0*(sin_lambda*dcos(lambda_star) &
    & - dsin(lambda_star)*cos_lambda)*(cos_lambda*dcos(lambda_star) + sin_lambda*dsin(lambda_star))))*G)

    F = F20 + F21 + F22
    end subroutine
    
    subroutine get_F_tide_Liulin(r_sat, r_body, m_body, F)
    implicit none
    double precision, intent(in) :: r_sat(3), r_body(3), m_body
    double precision, intent(out) :: F(3)
    double precision :: rp(3), rp1, rp3, r(3), Gm1, tlag
    rp = r_body
    rp1 = norm2(r_body)
    rp3 = rp1**3
    r = r_sat
    Gm1 = m_body
    tlag = 0.3d0*deg2rad
    
    call F34tX(rp, rp1, rp3, r, Gm1, tlag, F)
    
    endsubroutine get_F_tide_Liulin
    
    SUBROUTINE F34tX(rp, rp1, rp3, r, Gm1, tlag, F34t)
      IMPLICIT NONE
      ! 参数声明
      DOUBLE PRECISION, INTENT(IN)  :: rp(3), rp1, rp3, r(3), Gm1, tlag
      DOUBLE PRECISION, INTENT(OUT) :: F34t(3)

      ! 局部变量
      DOUBLE PRECISION :: Tidek2, r1, r3_local, r5
      DOUBLE PRECISION :: Sfai, Sfai2, Cfai, Cfai2, S2fai, C2fai
      DOUBLE PRECISION :: Sfaip, Sfaip2, Cfaip, Cfaip2, S2faip, C2faip
      DOUBLE PRECISION :: rr, SLamda, CLamda, S2Lamda, C2Lamda
      DOUBLE PRECISION :: rrp, WS, WC, WLamda, WX
      DOUBLE PRECISION :: SLx, CLx, S2Lx, C2Lx
      DOUBLE PRECISION :: P2fai, P21fai, P22fai, P2p, P21p, P22p
      DOUBLE PRECISION :: P2faip, P21faip, P22faip
      DOUBLE PRECISION :: SLL, CLL, S2LL, C2LL
      DOUBLE PRECISION :: FJ2X, FJ2Y, FJ2Z
      DOUBLE PRECISION :: FJ21X, FJ21Y, FJ21Z
      DOUBLE PRECISION :: FJ22X, FJ22Y, FJ22Z

      ! 基本量计算
      Tidek2 = 0.3D0 * Gm1 / rp3

      r1 = SQRT(r(1)**2 + r(2)**2 + r(3)**2)
      r3_local = r1 * r1 * r1
      r5 = r1 * r1 * r3_local

      Sfai  = r(3) / r1
      Sfai2 = Sfai**2
      Cfai  = SQRT(1.0D0 - Sfai2)
      Cfai2 = Cfai**2
      S2fai = 2.0D0 * Sfai * Cfai
      C2fai = Cfai2 - Sfai2

      Sfaip  = rp(3) / rp1
      Sfaip2 = Sfaip**2
      Cfaip  = SQRT(1.0D0 - Sfaip2)
      Cfaip2 = Cfaip**2
      S2faip = 2.0D0 * Sfaip * Cfaip
      C2faip = Cfaip2 - Sfaip2

      rr = SQRT(r(1)**2 + r(2)**2)
      SLamda  = r(2) / rr
      CLamda  = r(1) / rr
      S2Lamda = 2.0D0 * SLamda * CLamda
      C2Lamda = CLamda**2 - SLamda**2

      rrp = SQRT(rp(1)**2 + rp(2)**2)
      WS = rp(2) / rrp
      WC = rp(1) / rrp
      WLamda = ATAN2(WS, WC)          ! 取代原 LATAN 调用
      WX = WLamda + tlag

      SLx  = SIN(WX)
      CLx  = COS(WX)
      S2Lx = 2.0D0 * SLx * CLx
      C2Lx = CLx**2 - SLx**2

      ! 连带勒让德相关量
      P2fai  = 1.5D0 * Sfai2 - 0.5D0
      P21fai = 1.5D0 * S2fai
      P22fai = 3.0D0 * Cfai2
      P2p    = 3.0D0 * Sfai
      P21p   = 3.0D0 * C2fai / Cfai
      P22p   = -6.0D0 * Sfai

      P2faip  = 1.5D0 * Sfaip2 - 0.5D0
      P21faip = 1.5D0 * S2faip
      P22faip = 3.0D0 * Cfaip2

      SLL  = SLamda * CLx - CLamda * SLx
      CLL  = CLamda * CLx + SLamda * SLx
      S2LL = 2.0D0 * SLL * CLL
      C2LL = CLL**2 - SLL**2

      ! 长周期项 (J2)
      FJ2X = -Tidek2 * P2faip * (3.0D0 * P2fai + Sfai * P2p) * r(1) / r5
      FJ2Y = -Tidek2 * P2faip * (3.0D0 * P2fai + Sfai * P2p) * r(2) / r5
      FJ2Z = -Tidek2 * P2faip * ((3.0D0 * P2fai + Sfai * P2p) * r(3) - r1 * P2p) / r5

      ! 日潮项 (J21)
      FJ21X = -Tidek2 * P21faip / 3.0D0 * &
                ( (3.0D0 * P21fai + Sfai * P21p) * r(1) * CLL / r5 &
                 - (P21fai * SLamda) * SLL / (rr * r3_local) )
      FJ21Y = -Tidek2 * P21faip / 3.0D0 * &
                ( (3.0D0 * P21fai + Sfai * P21p) * r(2) * CLL / r5 &
                 + (P21fai * CLamda) * SLL / (rr * r3_local) )
      FJ21Z = -Tidek2 * P21faip / 3.0D0 * &
                ( (3.0D0 * P21fai + Sfai * P21p) * r(3) - r1 * P21p ) * CLL / r5

      ! 半日潮项 (J22)
      FJ22X = -Tidek2 * P22faip / 12.0D0 * &
                ( (3.0D0 * P22fai + Sfai * P22p) * r(1) * C2LL / r5 &
                 - 2.0D0 * (P22fai * SLamda) * S2LL / (rr * r3_local) )
      FJ22Y = -Tidek2 * P22faip / 12.0D0 * &
                ( (3.0D0 * P22fai + Sfai * P22p) * r(2) * C2LL / r5 &
                 + 2.0D0 * (P22fai * CLamda) * S2LL / (rr * r3_local) )
      FJ22Z = -Tidek2 * P22faip / 12.0D0 * &
                ( (3.0D0 * P22fai + Sfai * P22p) * r(3) - r1 * P22p ) * C2LL / r5

      ! 总潮汐力
      F34t(1) = FJ2X + FJ21X + FJ22X
      F34t(2) = FJ2Y + FJ21Y + FJ22Y
      F34t(3) = FJ2Z + FJ21Z + FJ22Z

      RETURN
    END SUBROUTINE F34tX
    
    end module