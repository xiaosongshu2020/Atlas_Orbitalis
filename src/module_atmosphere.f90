
    module module_atmosphere
    use module_const
    use module_time
    use module_planets
    
    integer, parameter, public :: model1 = 1 ! Model 1: Ellipsoidal atmospheric model for Sunday variation
    integer, parameter, public :: model2 = 2 ! Model 2: NRLMSISE-00 model
    
    integer, public :: atmosphereModel = 1 ! Default atmospheric model
    logical, public :: atmospherePrepared = .false.
    
    double precision, save :: r_p0, rho_p0, h_p0

    contains
    

    ! Subroutine to calculate atmospheric drag acceleration
    ! Inputs:
    !   UTC     - UTC time in Modified Julian Date (MJD)
    !   rs      - Satellite position in J2000 coordinate system [km]
    !   rv      - Satellite velocity in J2000 coordinate system [km/s]
    !   Cd      - Drag coefficient
    !   amrDrag - Area-to-mass ratio [m^2/kg]
    !   r_sun   - Sun position in J2000 coordinate system [km]
    ! Output:
    !   F       - Atmospheric drag acceleration [m/s^2]
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine atmospherical_drag_perturbation(UTC, rs, rv, Cd, amrDrag, F107, Ap, F)
    implicit none
    double precision, intent(in) :: UTC, rs(3), rv(3), Cd, amrDrag
    double precision, intent(inOut) :: F107, Ap
    double precision, intent(out) :: F(3)
    double precision rho
    double precision Omega(3), va(3), delta_v(3)
    double precision TT, solarPV(6), r_sun(3), h_perigee
    character*100 fn_atmosphere
    integer n_atmos_line, i_line
    character*10 tmp1
    
    double precision rsu(3), rvu(3), amrDragu
    integer threadNum
    !-------------------------------------------------------
    
    rsu = rs/length_unit
    rvu = rv/length_unit*time_unit
    amrDragu = amrDrag/(length_unit*1d3)**2*mass_unit
        
    select case (atmosphereModel)
    case(model1) ! �����������ָ������ģ�ͣ��ο������춯��ѧ���ۡ�
        
        !print*, "thread ", OMP_get_thread_num(), atmospherePrepared
        !pause
        
        if (.not. atmospherePrepared) then 
            call prepareAtmosphere(norm2(rs) - length_unit, F107, Ap)
        endif
        
        !get TT, solar position
        call UTC2TT(UTC, TT)
        call planetsPV(TT, nEarth, nSun, SolarPV)
        r_sun = solarPV(1:3)/length_unit
        
        call get_density_model1(rsu, r_sun, r_p0, rho_p0, h_p0, rho)
        !rho = rho_p0
        
    case(model2)
        ! NRLMSISE-00 model
        !call get_nrlmsise00_density(UTC, rsu, rho)
    case default
        write(*, *) "Invalid atmospherical model@subroutine atmospherical_drag_perturbation"
        write(*, *) "atmosphereModel = ", atmosphereModel
        F = 0d0
        return
    end select
    
    

    !get the velocity of atmosphere (va)
    Omega = (/-rsu(2), rsu(1), 0d0/)
    va = 3.141592653589793d0*2d0/(86400d0/time_unit)*Omega
    delta_v = rvu - va

    F = -0.5d0*Cd*amrDragu*rho*norm2(delta_v)*delta_v
    
    F = F*acc_unit
    
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_density_model1(rs, r_sun, r_p0, rho_p0, h_p0, rho)
    implicit none
    double precision, intent(in) :: rs(3), r_sun(3), r_p0, rho_p0, h_p0
    double precision, intent(out) :: rho
    double precision F_star, epsi, lambda_m
    double precision rm(3)
    double precision rs_norm, cos_psi_star, sin_phi, sigma_p0, H_r
    !---------------------------------------------

    F_star = 0.3D0
    epsi = 1D0/297D0
    lambda_m = 5.235987755982988D-01

    rm(1) = 1d0/norm2(r_sun)*(r_sun(1)*dcos(lambda_m) - r_sun(2)*dsin(lambda_m))
    rm(2) = 1d0/norm2(r_sun)*(r_sun(2)*dcos(lambda_m) + r_sun(1)*dsin(lambda_m))
    rm(3) = 1d0/norm2(r_sun)*r_sun(3)

    rs_norm = norm2(rs)
    cos_psi_star = (rs(1)*rm(1) + rs(2)*rm(2) + rs(3)*rm(3))/rs_norm
    sin_phi = rs(3)/rs_norm
    sigma_p0 = r_p0*(1-epsi*sin_phi**2)/9.991582491582491D-1
    H_r = h_p0 + 0.05D0*(rs_norm - r_p0)
    rho = rho_p0*(1d0 + F_star*cos_psi_star)*exp(-(rs_norm - sigma_p0)/H_r)

    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine prepareAtmosphere(height, F107, AP)
    implicit none
    double precision, intent(in) :: height, F107, Ap
    double precision F107n, Apn
    double precision TT, solarPV(6), r_sun(3), h_perigee
    character*100 fn_atmosphere
    integer n_atmos_line, i_line
    character*10 tmp1
    
    if (atmospherePrepared) then
        return
    endif
    h_perigee = height
    if (.true.) then
        
        if (h_perigee >= 1500d0) then
            h_perigee = 1500d0
        endif
        
        ! Check Ap, F107
        if((AP /= 5d0) .and. (Ap /= 10d0) .and. (Ap /= 15d0)) then
            write(*, *) 'Ap value error, which should be 5, 10 or 15.'
            Apn = 15d0
        else
            Apn = Ap
        endif
        if(F107 < 60d0) then
            write(*, *) 'F10.7 value error, which should be between 60 and 240.'
            F107n = 70d0
        elseif (F107 > 240d0) then
            write(*, *) 'F10.7 value error, which should be between 60 and 240.'
            F107n = 230d0
        else
            F107n = dble(int(floor(F107/20d0))*20 + 10)
        endif    
        
        write(fn_atmosphere, 1001) int(Apn), int(F107n)
1001    format('data/AtlasData/nrlmsise00_cheng/Ap=',I2.2,'_F10.7=',I3.3,'.txt')
        
        h_perigee = dble(floor(h_perigee/10d0))*10d0
        n_atmos_line = int(floor(h_perigee/10d0)) - 3
        open(204, file = trim(fn_atmosphere))
        do i_line = 1, n_atmos_line
            read(204,*) tmp1
        enddo
        read(204, *) tmp1, rho_p0, H_p0
        close(204)
        r_p0 = h_perigee + length_unit

    endif
            
    r_p0 = r_p0/length_unit
    rho_p0 = rho_p0/mass_unit*(length_unit*1000d0)**3
    h_p0 = h_p0/length_unit
            
    atmospherePrepared = .true.
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !================================================================================
    ! ---------------------------------------------------------------------------
    ! Subroutine: tropospheric_refraction
    ! Purpose: 
    !   Corrects the observed elevation angle of a celestial object by accounting
    !   for tropospheric refraction effects.
    !
    ! Input Parameters:
    !   Elev  - Original observed elevation angle (in radians). Ensure conversion
    !           from degrees to radians if necessary.
    !   T     - Temperature (in Kelvin). Represents the atmospheric temperature.
    !   fh    - Relative humidity (range: 0.0 to 1.0).
    !
    ! Output Parameters:
    !   Elev1 - Corrected elevation angle (in radians), which includes the
    !           tropospheric refraction correction.
    !
    ! Calculation Steps:
    !   1. Compute water vapor pressure:
    !      eh = 6.10 * fh * exp(17.15 * TC / (234.7 + TC))
    !      where TC = T - T0 (temperature in Celsius).
    !   2. Compute refractivity:
    !      Ns = 77.64 * pa / T + 3.734e5 * eh / (T^2)
    !   3. Compute refraction correction:
    !      dElev = Ns * 1e-6 / tan(Elev)
    !
    ! Constants:
    !   T0 = 273.2 K   - Reference temperature for 0°C.
    !   pa = 1024.0 mb - Standard atmospheric pressure.
    !
    ! Notes:
    !   - Ensure the input Elev is in radians. Convert degrees to radians using:
    !     Elev_rad = Elev_deg * pi / 180.0
    !   - The corrected elevation angle (Elev1) is added to the original angle.
    !   - Typical correction values range from 0.01-0.03 radians (0.0002-0.0005 degrees).
    !
    ! Example Usage:
    !   call tropospheric_refraction(0.5236d0, 293.15d0, 0.7d0, Elev1)
    !   ! Corrects for 30 degrees elevation with given temperature and humidity.
    ! ---------------------------------------------------------------------------
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine tropospheric_refraction(Elev, T, fh, Elev1)
    implicit none
    double precision, intent(in)  :: Elev
    double precision, optional :: T, fh
    double precision, intent(out) :: Elev1
    double precision              :: T0, pa, TC, eh, Ns, dElev
    
    ! Fixed media data at ground site
    T0 = 273.2d0          ! Temperature at 0 deg C [K]
    pa = 1024.0d0         ! Partial pressure of dry air [mb]
    
    if(.not.present(T)) T = 293.15d0
    if(.not.present(fh)) fh = 0.7d0
    
    ! Calculate temperature in Celsius
    TC = T - T0
    
    ! Calculate partial water pressure
    eh = 6.10d0 * fh * dexp(17.15d0 * TC / (234.7d0 + TC))
    
    ! Calculate refractivity
    Ns = 77.64d0 * pa / T + 3.734d5 * eh / (T**2)
    
    ! Calculate elevation correction (ensure Elev is in radians)
    dElev = Ns * 1d-6 / dtan(Elev)
    
    ! Return corrected elevation
    Elev1 = Elev + dElev
    
    end subroutine tropospheric_refraction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endmodule