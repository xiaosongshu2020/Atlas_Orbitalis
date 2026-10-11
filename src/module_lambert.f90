module module_LambertSolution
    use module_const
    implicit none
    
    contains


!--------------------------------------------------------------------------
! Elements: Computes orbital elements from two given position vectors and
!           associated times
!
! Source: the algorithm and the documentation below follow O. Montenbruck &
!         E. Gill, "Satellite Orbits -- Models, Methods and Applications",
!         Springer, 2000 (the "Elements" routine of the book's code suite).
!         Kept here with attribution; not relicensed by this project.
!
! Inputs:
!   GM        Gravitational coefficient
!             (gravitational constant * mass of central body)
!   Mjd_a     Time t_a (Modified Julian Date)
!   Mjd_b     Time t_b (Modified Julian Date)
!   r_a       Position vector at time t_a
!   r_b       Position vector at time t_b
!
! Outputs:
!             Keplerian elements (a,ecc,i,ou,ol,M)
!               a      Semimajor axis
!               ecc    Eccentricity
!               i      Inclination [rad]
!               ou     Longitude of the ascending node [rad]
!               ol     Argument of pericenter  [rad]
!               M      Mean anomaly  [rad]
!             at time t_a
!
! Notes:
!   The function cannot be used with state vectors describing a circular
!   or non-inclined orbit.
!
!--------------------------------------------------------------------------
subroutine LambertSolution(GM, Mjd_a, Mjd_b, r_a, r_b, a, ecc, i, ou, ol, M)
    implicit none
    
    ! Input/Output parameters
    real(8), intent(in)  :: GM         ! Gravitational coefficient
    real(8), intent(in)  :: Mjd_a      ! Time t_a (Modified Julian Date)
    real(8), intent(in)  :: Mjd_b      ! Time t_b (Modified Julian Date)
    real(8), intent(in)  :: r_a(3)     ! Position vector at time t_a
    real(8), intent(in)  :: r_b(3)     ! Position vector at time t_b
    real(8), intent(out) :: a          ! Semimajor axis
    real(8), intent(out) :: ecc        ! Eccentricity
    real(8), intent(out) :: i          ! Inclination [rad]
    real(8), intent(out) :: ou         ! Longitude of the ascending node [rad]
    real(8), intent(out) :: ol         ! Argument of pericenter [rad]
    real(8), intent(out) :: M          ! Mean anomaly [rad]
    
    ! Local variables
    real(8) :: pi2, s_a, s_b, s_0, fac, tau, eta, p
    real(8) :: cos_dnu, sin_dnu, ecos_nu, esin_nu, nu, n
    real(8) :: E, sinhH
    real(8) :: e_a(3), e_0(3), r_0(3), W(3), u

    ! Constants
    pi2 = 2.0d0 * 3.14159265358979323846d0
    
    ! Calculate vector r_0 (fraction of r_b perpendicular to r_a) and the
    ! magnitudes of r_a, r_b and r_0
    s_a = sqrt(r_a(1)**2 + r_a(2)**2 + r_a(3)**2)
    e_a = r_a / s_a
    
    s_b = sqrt(r_b(1)**2 + r_b(2)**2 + r_b(3)**2)
    fac = dot_product(r_b, e_a)
    r_0 = r_b - fac * e_a
    
    s_0 = sqrt(r_0(1)**2 + r_0(2)**2 + r_0(3)**2)
    e_0 = r_0 / s_0
    
    ! Inclination and ascending node
    call cross_product(e_a, e_0, W)
    ou = atan2(W(1), -W(2))     ! Long. ascend. node
    ou = mod(ou, pi2)
    
    i = atan2(sqrt(W(1)**2 + W(2)**2), W(3))  ! Inclination
    
    if (i == 0.0d0) then
        u = atan2(r_a(2), r_a(1))
    else
        u = atan2(e_a(3), (-e_a(1)*W(2) + e_a(2)*W(1)))
    end if
    
    ! Semilatus rectum
    tau = sqrt(GM) * 86400.0d0 * abs(Mjd_b - Mjd_a)
    eta = FindEta(r_a, r_b, tau)
    p = (s_a * s_0 * eta / tau)**2
    
    ! Eccentricity, true anomaly and argument of perihelion
    cos_dnu = fac / s_b
    sin_dnu = s_0 / s_b
    ecos_nu = p / s_a - 1.0d0
    esin_nu = (ecos_nu * cos_dnu - (p / s_b - 1.0d0)) / sin_dnu
    
    ecc = sqrt(ecos_nu**2 + esin_nu**2)
    nu = atan2(esin_nu, ecos_nu)
    ol = mod(u - nu, pi2)
    
    ! Perihelion distance, semimajor axis and mean motion
    a = p / (1.0d0 - ecc**2)
    n = sqrt(GM / abs(a**3))
    
    ! Mean anomaly and time of perihelion passage
    if (ecc < 1.0d0) then
        E = atan2(sqrt((1.0d0-ecc)*(1.0d0+ecc)) * esin_nu, ecos_nu + ecc**2)
        M = mod(E - ecc * sin(E), pi2)
    else
        sinhH = sqrt((ecc-1.0d0)*(ecc+1.0d0)) * esin_nu / (ecc + ecc*ecos_nu)
        M = ecc * sinhH - log(sinhH + sqrt(1.0d0 + sinhH**2))
    end if
    
end subroutine LambertSolution

!--------------------------------------------------------------------------
! Helper function to calculate cross product
!--------------------------------------------------------------------------
subroutine cross_product(a, b, c)
    implicit none
    real(8), intent(in)  :: a(3), b(3)
    real(8), intent(out) :: c(3)
    
    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
end subroutine cross_product

!--------------------------------------------------------------------------
! FindEta: Computes the sector-triangle ratio from two position vectors and
!          the intermediate time
!
! Inputs:
!   r_a        Position at time t_a
!   r_b        Position at time t_b
!   tau        Normalized time (sqrt(GM)*(t_a-t_b))
!
! Output:
!   eta2       Sector-triangle ratio
!
!--------------------------------------------------------------------------
function FindEta(r_a, r_b, tau) result(eta2)
    implicit none
    
    ! Input parameters
    real(8), intent(in) :: r_a(3)     ! Position at time t_a
    real(8), intent(in) :: r_b(3)     ! Position at time t_b
    real(8), intent(in) :: tau        ! Normalized time
    real(8)             :: eta2       ! Output: Sector-triangle ratio
    
    ! Local variables
    integer, parameter :: maxit = 30
    real(8), parameter :: delta = 1.0d-16
    real(8) :: s_a, s_b, kappa, m, l, eta_min
    real(8) :: eta1, F1, F2, d_eta
    integer :: i
    
    ! Calculate magnitudes and dot product
    s_a = sqrt(r_a(1)**2 + r_a(2)**2 + r_a(3)**2)
    s_b = sqrt(r_b(1)**2 + r_b(2)**2 + r_b(3)**2)
    
    kappa = sqrt(2.0d0 * (s_a*s_b + dot_product(r_a, r_b)))
    
    m = (tau**2) / (kappa**3)
    l = (s_a + s_b) / (2.0d0 * kappa) - 0.5d0
    
    eta_min = sqrt(m / (l + 1.0d0))
    
    ! Start with Hansen's approximation
    eta2 = (12.0d0 + 10.0d0 * sqrt(1.0d0 + (44.0d0/9.0d0) * m / (l + 5.0d0/6.0d0))) / 22.0d0
    eta1 = eta2 + 0.1d0
    
    ! Secant method
    F1 = F_function(eta1, m, l)
    F2 = F_function(eta2, m, l)
    
    i = 0
    do while (abs(F2 - F1) > delta)
        d_eta = -F2 * (eta2 - eta1) / (F2 - F1)
        
        eta1 = eta2
        F1 = F2
        
        do while ((eta2 + d_eta) <= eta_min)
            d_eta = d_eta / 2.0d0
        end do
        
        eta2 = eta2 + d_eta
        F2 = F_function(eta2, m, l)
        
        i = i + 1
        
        if (i == maxit) then
            write(*, *) 'Convergence problems in FindEta'
            exit
        end if
    end do
    
end function FindEta

!--------------------------------------------------------------------------
!    F : local function for use by FindEta()
!    F = 1 - eta +(m/eta**2)*wu(m/eta**2-l)
!
!--------------------------------------------------------------------------
function F_function(eta, m, l) result(ou)
    implicit none
    real(8), intent(in) :: eta, m, l
    real(8) :: ou
    
    ! Local variables
    real(8) :: epsilon, wl, wu, a, g
    integer :: n
    
    ! Set epsilon approximately equivalent to 100*eps in MATLAB
    epsilon = 100.0d0 * 2.2204460492503131d-16  ! 100 * machine epsilon for double precision
    
    wl = m/(eta*eta) - l
    
    if (abs(wl) < 0.1d0) then
        ! Series expansion
        wu = 4.0d0/3.0d0
        a = 4.0d0/3.0d0
        n = 0
        
        do
            n = n + 1
            a = a * (wl * (n+2.0d0) / (n+1.5d0))
            wu = wu + a
            if (abs(a) < epsilon) then
                exit
            end if
        end do
        
    else if (wl > 0.0d0) then
        g = 2.0d0 * asin(sqrt(wl))
        wu = (2.0d0*g - sin(2.0d0*g)) / (sin(g)**3)
        
    else  ! wl <= 0
        ! g = 2.0*log(sqrt(-wl) + sqrt(1.0-wl))  ! =2.0*arsinh(sqrt(-wl))
        g = 2.0d0 * asinh(sqrt(-wl))
        wu = (sinh(2.0d0*g) - 2.0d0*g) / (sinh(g)**3)
    end if
    
    ou = (1.0d0 - eta + (wl+l)*wu)
    
end function F_function

end module module_LambertSolution

! ! test program
! !--------------------------------------------------------------------------
! program test_LambertSolution
!     implicit none

!     ! Test variables
!     real(8) :: GM, Mjd_a, Mjd_b
!     real(8) :: r_a(3), r_b(3)
!     real(8) :: a, e, i, ou, ol, M

!     ! Initialize test inputs
!     GM = 398600.4418d0  ! Gravitational parameter for Earth [km^3/s^2]
!     Mjd_a = 58000.0d0   ! Example Modified Julian Date for time t_a
!     Mjd_b = 58000.25d0   ! Example Modified Julian Date for time t_b
!     r_a = [42000.0d0, 0.0d0, 0.0d0]  ! Position vector at t_a [km]
!     r_b = [0.0d0, 42000.0d0, 0.0d0]  ! Position vector at t_b [km]

!     ! Call the subroutine
!     call LambertSolution(GM, Mjd_a, Mjd_b, r_a, r_b, a, e, i, ou, ol, M)

!     ! Print results
!     print *, "Test Results for LambertSolution:"
!     print *, "Semimajor axis (a): ", a
!     print *, "Eccentricity (e): ", e
!     print *, "Inclination (i): ", i
!     print *, "Longitude of ascending node (ou): ", ou
!     print *, "Argument of pericenter (ol): ", ol
!     print *, "Mean anomaly (M): ", M

! end program test_LambertSolution