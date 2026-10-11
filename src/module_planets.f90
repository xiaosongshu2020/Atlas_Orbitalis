!================================================================================
! MODULE: module_planets
!
! DESCRIPTION:
! Planetary ephemeris provider: position and velocity of solar system bodies
! (Sun, Moon, planets) with respect to a chosen origin, in the J2000 frame.
!
! The ephemeris source is selected by the module variable planetsEph:
!   useJplEph          JPL DE ephemeris via pleph (high accuracy)
!   useAnaEph          simple analytical models (low accuracy; only the
!                      Earth-centered Moon and Sun are supported)
!   useEarthMoonCRTBP  circular lunar orbit for the Earth-Moon CRTBP model
!
! SAFETY NET: when useJplEph is selected and the JPL ephemeris cannot be read
! (file missing, unreadable or corrupt, epoch outside the file's time span),
! planetsPV does NOT abort.  The Moon and the Sun fall back to the analytical
! models of the useAnaEph branch -- still in km/km/s, with a warning printed
! once per body -- while any other body stops with an explanatory message,
! since no analytical model exists for it and its perturbation must not be
! dropped silently.
!
! MAIN ENTRY: planetsPV(TT, origin, object, PV)
!   TT     epoch in MJD (TT time scale)
!   origin reference body number
!   object target body number
!   PV     position and velocity of object w.r.t. origin, J2000 frame
!          The JPL and CRTBP paths return SI units (km, km/s); the
!          analytical paths return normalized units (Earth radius)
!================================================================================

    module module_planets
    use module_const
    use module_elements
    implicit none

    ! ephemeris source options for planetsEph
    integer, parameter  ::  useJplEph = 1
    integer, parameter  ::  useAnaEph = 2
    integer, parameter  ::  useEarthMoonCRTBP = 3

    ! body numbers (same convention as the JPL ephemeris)
    integer, parameter  ::  nMercury = 1
    integer, parameter  ::  nVenus = 2
    integer, parameter  ::  nEarth = 3
    integer, parameter  ::  nMars = 4
    integer, parameter  ::  nJupiter = 5
    integer, parameter  ::  nSaturn = 6
    integer, parameter  ::  nUranus = 7
    integer, parameter  ::  nPluto = 9
    integer, parameter  ::  nMoon = 10
    integer, parameter  ::  nSun = 11
    integer, parameter  ::  nSolarSystemBarycenter = 12
    integer, parameter  ::  nEarthMoonBarycenter = 13
    integer, parameter  ::  nNutations = 14
    integer, parameter  ::  nLibrations = 15

    ! planet/Earth mass ratios (Earth mass = 1.0)
    double precision, parameter  ::  moonEarthMassRatio = 1.2300037226627106D-2      ! Moon
    double precision, parameter  ::  sunEarthMassRatio = 332946.03864607110D0        ! Sun
    double precision, parameter  ::  mercuryEarthMassRatio = 0.0553d0          ! Mercury
    double precision, parameter  ::  venusEarthMassRatio = 0.815d0            ! Venus
    double precision, parameter  ::  earthEarthMassRatio = 1.0d0              ! Earth
    double precision, parameter  ::  marsEarthMassRatio = 0.107d0             ! Mars
    double precision, parameter  ::  jupiterEarthMassRatio = 317.8d0          ! Jupiter
    double precision, parameter  ::  saturnEarthMassRatio = 95.16d0           ! Saturn
    double precision, parameter  ::  uranusEarthMassRatio = 14.54d0           ! Uranus
    double precision, parameter  ::  neptuneEarthMassRatio = 17.15d0          ! Neptune
    double precision, parameter  ::  plutoEarthMassRatio = 0.00218d0          ! Pluto

    ! mean (volumetric) radii of the Sun, planets and Moon (km)
    double precision, parameter   ::  radius_sun = 696000d0      ! Sun
    double precision, parameter   ::  radius_mercury = 2439.7d0  ! Mercury
    double precision, parameter   ::  radius_venus = 6051.8d0    ! Venus
    double precision, parameter   ::  radius_earth = 6371.0d0    ! Earth
    double precision, parameter   ::  radius_moon = 1736.4d0     ! Moon
    double precision, parameter   ::  radius_mars = 3389.5d0     ! Mars
    double precision, parameter   ::  radius_jupiter = 69911d0   ! Jupiter
    double precision, parameter   ::  radius_saturn = 58232d0   ! Saturn
    double precision, parameter   ::  radius_uranus = 25362d0    ! Uranus
    double precision, parameter   ::  radius_neptune = 24622d0   ! Neptune
    double precision, parameter   ::  radius_pluto = 1188.3d0    ! Pluto

    ! ephemeris source used by planetsPV (see the options above), default: JPL
    integer, save   ::  planetsEph = 1

    ! ---- JPL fallback support (see planetsPV / jplFallback) ----
    ! status reported by pleph_status (src/PLEPH.for)
    integer, parameter  ::  jplOk = 0           ! ephemeris read successfully
    integer, parameter  ::  jplInitFailed = 1   ! file missing or unreadable
    integer, parameter  ::  jplOutOfRange = 2   ! epoch outside the file's span
    integer, parameter  ::  jplReadError = 3    ! record read error

    ! central-difference step for the fallback velocities (days).  The
    ! harmonics of the analytical series are slow, so 0.01 day keeps the
    ! truncation error around 0.3 mm/s -- orders below the model error --
    ! while staying far from the round-off cliff.
    double precision, parameter  ::  fallbackFdStep = 0.01d0

    ! the fallback warning sits on the integration hot path (RKF78 evaluates
    ! the forces 13 times per step), so warn once per body
    logical, save   ::  warnedMoonFallback = .false.
    logical, save   ::  warnedSunFallback = .false.

    contains

    !============================================================================
    ! Position and velocity of `object` with respect to `origin` at epoch TT.
    !   input : TT     epoch in MJD (TT time scale)
    !           origin reference body number
    !           object target body number
    !   output: PV     position and velocity of object w.r.t. origin, J2000 frame
    !           the JPL and CRTBP paths return SI units (km, km/s);
    !           the analytical paths return normalized units (Earth radius)
    ! Calls with identical arguments return the cached result.
    !============================================================================
    subroutine planetsPV(TT, origin, object, PV)
    implicit none
    double precision, intent(in)  ::  TT
    integer, intent(in) ::  origin, object
    double precision, intent(out) ::  PV(6)
    double precision, save :: lastTT = 0d0
    double precision, save :: lastPV(6) = 0d0
    integer, save :: lastOrigin = 0
    integer, save :: lastObject = 0
    integer  ::  iJplStatus

    ! return the cached result when called again with the same arguments
    if (dabs(TT - lastTT) < 1d-8 .and. origin == lastOrigin .and. object == lastObject) then
        PV = lastPV
        return
    endif

    select case (planetsEph)
    case (useJplEph)
        ! pleph works in JD and returns AU and AU/day
        call pleph(TT + 2400000.5d0, object, origin, PV)
        call pleph_status(iJplStatus)
        if (iJplStatus == jplOk) then
            PV(1:3) = PV(1:3)*AU_TO_EARTH_RADIUS                   ! AU -> Earth radius
            PV(4:6) = PV(4:6)*AU_TO_EARTH_RADIUS/86400d0*time_unit ! AU/day -> Earth radius/time unit
        else
            ! JPL ephemeris unavailable: use the analytical models instead of
            ! aborting.  They fill PV in Earth radii (and Earth radii/time
            ! unit), so the common conversion below still hands the caller
            ! the same km/km/s contract as the JPL path.
            call jplFallback(TT, origin, object, iJplStatus, PV)
        endif
    case (useAnaEph)
        ! the analytical models return before the SI conversion below,
        ! i.e. their output stays in normalized units
        if (origin == nEarth .and. object == nMoon) then
            call fit_lunar_position(TT, PV)
            return
        endif
        if (origin == nEarth .and. object == nSun) then
            call mean_solar_position(TT, PV)
            return
        endif
        write(*, *) "Error:analytical planets' trajectory doesn't ", &
            &   "support this object or origin! (subroutine planetsPV)"
        STOP
    case (useEarthMoonCRTBP)
        ! circular Earth-Moon orbit, for CRTBP computations
        if (origin == nEarth .and. object == nMoon) then
            call circular_lunar_position(TT, PV)
        else
            write(*, *) "Error:analytical planets' trajectory doesn't ", &
            &   "support this object or origin! (subroutine planetsPV)"
            STOP
        endif
    case default
        write(*, *) "Error: unaccepted planets trajectory type! (subroutine planetsPV)"
        STOP
    endselect

    ! convert to SI units: Earth radius -> km, Earth radius/time unit -> km/s
    ! (the JPL branch, the CRTBP branch and the fallback all enter here with
    ! PV already in Earth radii and Earth radii/time unit)
    PV(1:3) = PV(1:3)*length_unit
    PV(4:6) = PV(4:6)*length_unit/time_unit

    ! cache the result for repeated calls
    lastOrigin = origin
    lastObject = object
    lastTT = TT
    lastPV = PV

    endsubroutine

    !===== JPL fallback =======================================================
    ! Analytical replacement for the JPL ephemeris, used when PLEPH reports a
    ! failure.  Fills PV in Earth radii and Earth radii/time unit, so that the
    ! common conversion in planetsPV applies unchanged and the caller gets the
    ! same km/km/s contract as the JPL path.
    ! Only the Earth-centred Moon and Sun are supported; any other body stops
    ! with an explanation rather than silently dropping its perturbation.
    !==========================================================================
    subroutine jplFallback(TT, origin, object, status, PV)
    implicit none
    double precision, intent(in)  ::  TT
    integer, intent(in)           ::  origin, object, status
    double precision, intent(out) ::  PV(6)
    double precision  ::  r0(3), rPlus(3), rMinus(3)

    if (origin /= nEarth .or. (object /= nMoon .and. object /= nSun)) then
        call fatalNoAnalyticalModel(object, status)
    endif

    call analyticalPosition(TT, object, r0)
    call analyticalPosition(TT + fallbackFdStep, object, rPlus)
    call analyticalPosition(TT - fallbackFdStep, object, rMinus)

    PV(1:3) = r0
    ! velocity by central difference: the analytic velocities available in
    ! this module are unusable, because mean_solar_position inherits mu = 1
    ! from kpl2cts while its elements are heliocentric (see that routine)
    PV(4:6) = (rPlus - rMinus)/(2d0*fallbackFdStep)*time_unit/86400d0

    call warnFallbackOnce(object, status)

    endsubroutine

    ! Position of `object` with respect to the Earth at epoch TT (MJD, TT),
    ! in Earth radii, J2000 frame.
    subroutine analyticalPosition(TT, object, r)
    implicit none
    double precision, intent(in)  ::  TT
    integer, intent(in)           ::  object
    double precision, intent(out) ::  r(3)
    double precision  ::  pvSun(6)

    select case (object)
    case (nMoon)
        call fit_lunar_position(TT, r)
    case (nSun)
        call mean_solar_position(TT, pvSun)
        r = pvSun(1:3)
    endselect

    endsubroutine

    ! Announce the fallback.  Warn once per body: this code sits on the
    ! integrator hot path (RKF78 evaluates the forces 13 times per step).
    subroutine warnFallbackOnce(object, status)
    implicit none
    integer, intent(in) :: object, status

    select case (object)
    case (nMoon)
        if (warnedMoonFallback) return
        warnedMoonFallback = .true.
    case (nSun)
        if (warnedSunFallback) return
        warnedSunFallback = .true.
    endselect

    ! diagnostics go to stderr (unit 0), never to standard output
    write(0, *) "Warning @planetsPV: JPL DE ephemeris unavailable (" // &
        &   trim(jplStatusText(status)) // ")."
    select case (object)
    case (nMoon)
        write(0, *) "  Using the analytical Moon ephemeris instead " // &
            &   "(truncated ELP series, velocity by central difference, " // &
            &   "step 0.01 day). Accuracy is degraded."
    case (nSun)
        write(0, *) "  Using the analytical Sun ephemeris instead " // &
            &   "(mean Keplerian orbit, velocity by central difference, " // &
            &   "step 0.01 day). Accuracy is degraded."
    endselect

    endsubroutine

    ! No analytical model exists for this body: stop with an explanation
    ! instead of silently dropping its perturbation.
    subroutine fatalNoAnalyticalModel(object, status)
    implicit none
    integer, intent(in) :: object, status
    character(len=80)   :: ephFile

    call pleph_filename(ephFile)

    write(0, *) "Error @planetsPV: JPL DE ephemeris unavailable (" // &
        &   trim(jplStatusText(status)) // ")."
    write(0, '(3a,i0,a)') "  Body ", trim(bodyName(object)), " (id ", &
        &   object, ") has no analytical model; its perturbation is not " // &
        &   "dropped silently."
    write(0, *) "  Ephemeris file: ", trim(ephFile)
    write(0, *) "  Options:"
    write(0, *) "    1) run from the repository root (all data paths are " // &
        &   "relative) or restore the file above;"
    write(0, *) "    2) choose an epoch inside the ephemeris time span;"
    write(0, *) "    3) switch this body off on the satellite's opMod."

    ! flush first: otherwise the runtime prints the ERROR STOP backtrace
    ! above this explanation
    flush(0)
    error stop

    endsubroutine

    ! Reason text for a code returned by pleph_status
    function jplStatusText(status) result(txt)
    implicit none
    integer, intent(in) :: status
    character(len=48)   :: txt

    select case (status)
    case (jplInitFailed)
        txt = "cannot open or read the ephemeris file"
    case (jplOutOfRange)
        txt = "epoch outside the ephemeris time span"
    case (jplReadError)
        txt = "error reading an ephemeris record"
    case default
        txt = "unknown error"
    endselect

    endfunction

    ! Human-readable name of a JPL body number
    function bodyName(object) result(nam)
    implicit none
    integer, intent(in) :: object
    character(len=8)    :: nam

    select case (object)
    case (nMercury)
        nam = "Mercury"
    case (nVenus)
        nam = "Venus"
    case (nEarth)
        nam = "Earth"
    case (nMars)
        nam = "Mars"
    case (nJupiter)
        nam = "Jupiter"
    case (nSaturn)
        nam = "Saturn"
    case (nUranus)
        nam = "Uranus"
    case (nPluto)
        nam = "Pluto"
    case (nMoon)
        nam = "Moon"
    case (nSun)
        nam = "Sun"
    case default
        nam = "body"
    endselect

    endfunction

    ! Approximate geocentric solar position (Sun w.r.t. Earth) from the
    ! mean Keplerian elements of the Earth's heliocentric orbit.
    !   input : TT in MJD
    !   output: PV in normalized units (Earth radius, Earth radius/time unit)
    ! CAUTION: crude approximation -- the mean elements ignore planetary
    ! perturbations; the caller warns when this replaces the JPL ephemeris.
    subroutine mean_solar_position(TT, PV)
    implicit none
    double precision TT, PV(6)
    double precision d, T
    double precision elements(6)

    d = TT - 51544.50074287037d0 ! days since J2000.0
    T = d/36525d0                ! Julian centuries since J2000.0

    ! mean Keplerian elements of the Earth's orbit: (a, e, i, Omega, omega, M)
    elements(1) = 1.00000102d0*AU_TO_EARTH_RADIUS
    elements(2) = 0.016709d0
    elements(3) = 23.4393d0*deg2rad
    elements(4) = 0d0
    elements(5) = 282.9373d0*deg2rad + 0.32d0*deg2rad*T
    elements(6) = 357.5291d0*deg2rad + 0.9856d0*deg2rad*d

    call kpl2cts(elements, PV)

    ! kpl2cts assumes mu = 1 (the Earth's GM in normalized units), but these
    ! are heliocentric elements: the Sun's GM is sunEarthMassRatio times
    ! larger, so the velocity it returns is too small by
    ! sqrt(sunEarthMassRatio) = 577 (0.05 km/s instead of 29.8 km/s).
    PV(4:6) = PV(4:6)*dsqrt(sunEarthMassRatio)

    end subroutine

    ! Mean lunar orbit from fixed Keplerian elements.
    !   input : TT in MJD
    !   output: PV in normalized units (Earth radius, Earth radius/time unit)
    ! CAUTION: low accuracy, the error can reach a few degrees
    subroutine mean_lunar_position(TT, PV)
    implicit none
    double precision TT, PV(6)
    double precision elements(6)
    double precision d, T

    d = TT - 51544.50074287037d0 ! days since J2000.0
    T = d/36525d0                ! Julian centuries since J2000.0

    ! mean Keplerian elements of the lunar orbit: (a, e, i, Omega, omega, M)
    elements(1) = 384747.981d0/length_unit
    elements(2) = 0.054880
    elements(3) = 5.1298d0*deg2rad
    elements(4) = 125.0446d0*deg2rad - 1934.14d0*deg2rad*T
    elements(5) = 318.3087d0*deg2rad + 6003.15d0*deg2rad*T
    elements(6) = 134.9634d0*deg2rad + 13.0650d0*deg2rad*d

    call kpl2cts(elements, PV)

    end subroutine

    ! Circular lunar orbit approximation (e = 0), suitable for the
    ! Earth-Moon CRTBP model.
    !   input : TT in MJD
    !   output: PV in normalized units (Earth radius, Earth radius/time unit)
    subroutine circular_lunar_position(TT, PV)
    implicit none
    double precision TT, PV(1:6)
    double precision elements(6)
    double precision d, T

    d = TT - 51544.5d0   ! days since J2000.0
    T = d/36525d0        ! Julian centuries since J2000.0

    ! mean Keplerian elements with the eccentricity set to zero
    elements(1) = 384400d0/length_unit
    elements(2) = 0d0
    elements(3) = 5.1298d0*deg2rad
    elements(4) = 125.0446d0*deg2rad - 1934.14d0*deg2rad*T
    elements(5) = 318.3087d0*deg2rad + 6003.15d0*deg2rad*T
    elements(6) = 134.9634d0*deg2rad + 13.0650d0*deg2rad*d

    call kpl2cts(elements, PV)

    end subroutine

    ! Analytical lunar position from a low-order series expansion
    ! (position only, no velocity).
    !   input : TT in MJD
    !   output: r_moon  lunar position in the J2000 frame, Earth radius
    !
    ! Method: series for the ecliptic longitude, latitude and parallax with
    ! the fundamental arguments alpha(1:10); geocentric distance from the
    ! parallax (1/sin(pi)); rotate ecliptic -> equator of date (obliquity
    ! ~23.44 deg); finally precess to the J2000 frame.
    subroutine fit_lunar_position(TT, r_moon)
    implicit none
    double precision TT, r_moon(3)
    double precision t
    double precision pi_const, deg2rad
    double precision lambda, beta, pi
    double precision K(14), alpha(10)
    double precision R(3), R_hat(3), R_norm
    double precision PR(3, 3), PR_T(3, 3)
    integer i

    pi_const = 3.141592653589793d0
    deg2rad = pi_const/180d0

    ! series amplitude coefficients K (rad)
    K(1:4) = (/6.29d0, -1.27d0, 0.66d0, 0.21d0/)*deg2rad
    K(5:8) = (/-0.19d0, -0.11d0, 5.13d0, 0.28d0/)*deg2rad
    K(9:12) = (/-0.28d0, 0.17d0, 0.0518d0, 0.0095d0/)*deg2rad
    K(13:14) = (/0.0078d0, 0.0028d0/)*deg2rad

    t = (TT - 51544.5d0)/36525d0  ! Julian centuries since J2000.0

    ! fundamental arguments of the series (deg -> rad)
    alpha = (/134.9d0 + 477198.85d0*t, 259.2d0 - 413335.38d0*t, &
            & 235.7d0 + 890534.23d0*t, 269.9d0 + 954397.70d0*t, &
            & 357.5d0 + 35999.05d0*t, 186.6d0 + 966404.05d0*t, &
            & 93.3d0 + 483202.13d0*t, 228.2d0 + 960400.87d0*t, &
            & 318.3d0 + 6003.18d0*t, 217.6d0 - 407332.20d0*t/)*deg2rad

    lambda = 0d0;
    beta = 0d0;
    pi = 0d0;

    ! ecliptic longitude: base term + 6 periodic terms
    lambda = (218.32d0 + 481267.883d0*t)*deg2rad
    do i = 1, 6
        lambda = lambda + K(i)*dsin(alpha(i))
    end do

    ! ecliptic latitude: 4 periodic terms
    do i = 7,10
        beta = beta + K(i)*dsin(alpha(i))
    end do

    ! lunar parallax: base term + 4 periodic terms
    pi = 0.9508d0*deg2rad
    do i = 11, 14
        pi = pi + K(i)*dcos(alpha(i - 10))
    end do

    ! Earth-Moon distance in Earth radii
    R_norm = 1d0/dsin(pi)

    ! unit vector in the equator-of-date frame (ecliptic coordinates
    ! rotated by the obliquity, ~23.44 deg)
    R_hat(1) = dcos(beta)*dcos(lambda)
    R_hat(2) = 0.9175d0*dcos(beta)*dsin(lambda) - 0.3978d0*dsin(beta)
    R_hat(3) = 0.3978d0*dcos(beta)*dsin(lambda) + 0.9175*dsin(beta)

    ! lunar position in the equator-of-date frame, Earth radii
    R = R_hat*R_norm

    ! IAU1976 precession matrix J2000 -> of date; its transpose rotates
    ! the of-date vector back to J2000
    call PR_get_1314(t, PR)
    PR_T = transpose(PR)
    r_moon = matmul(PR_T, R)

    r_moon = r_moon

    end subroutine

    ! IAU1976 precession matrix from J2000 to the mean equator/equinox
    ! of date: PR = Rz(-zA) * Ry(thetaA) * Rz(-zetaA).
    !   input : t2  Julian centuries since J2000
    !   output: PR  precession matrix (3x3)
    subroutine PR_get_1314(t2, PR)
    implicit none
    double precision zetaA, thetaA, zA, matrix_zA(3, 3), matrix_thetaA(3, 3), matrix_zetaA(3, 3), PR1(3, 3), PR(3, 3)
    double precision pi, t2
    parameter(pi=3.141592653589793d0)

    ! precession angles in arcsec, converted to rad
    zetaA=(2306.2181d0*t2+0.30188d0*t2**2)/648000d0*pi
    zA=(2306.2181d0*t2+1.09468d0*t2**2)/648000d0*pi
    thetaA=(2004.3109d0*t2-0.42665d0*t2**2)/648000d0*pi

    call Rzmatrix_1314(-zA,matrix_zA)
    call Rymatrix_1314(thetaA,matrix_thetaA)
    call Rzmatrix_1314(-zetaA,matrix_zetaA)
    PR1 = matmul(matrix_thetaA,matrix_zetaA)
    PR = matmul(matrix_zA,PR1)

    end subroutine

    ! elementary rotation matrix about the z axis
    subroutine Rzmatrix_1314(theta,a)
    implicit none
    double precision theta,a(3,3)
    a(1,1)=dcos(theta)
    a(1,2)=dsin(theta)
    a(1,3)=0
    a(2,1)=-dsin(theta)
    a(2,2)=dcos(theta)
    a(2,3)=0
    a(3,1)=0
    a(3,2)=0
    a(3,3)=1D0
    return
    end subroutine

    ! elementary rotation matrix about the y axis
    subroutine Rymatrix_1314(theta,a)
    implicit none
    double precision theta,a(3,3)
    a(1,1)=dcos(theta)
    a(1,2)=0
    a(1,3)=-dsin(theta)
    a(2,1)=0
    a(2,2)=1D0
    a(2,3)=0
    a(3,1)=dsin(theta)
    a(3,2)=0
    a(3,3)=dcos(theta)
    return
    end subroutine

    ! elementary rotation matrix about the x axis
    subroutine Rxmatrix_1314(theta,a)
    implicit none
    double precision theta,a(3,3)
    a(1,1)=1D0
    a(1,2)=0
    a(1,3)=0
    a(2,1)=0
    a(2,2)=dcos(theta)
    a(2,3)=dsin(theta)
    a(3,1)=0
    a(3,2)=-dsin(theta)
    a(3,3)=dcos(theta)
    return
    end subroutine

    endModule
