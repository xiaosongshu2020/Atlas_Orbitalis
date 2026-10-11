!! BeiDou MEO satellite orbit prediction
!! 2025-07-28
!!
!! This program demonstrates 1-day orbit prediction for a BeiDou MEO satellite
!! using the numerical propagator with a 21×21 gravity field and luni-solar
!! third-body perturbations.
!!
!! STK comparison notes:
!!   模型一：21×21 EGM2008 + 日月第三体引力
!!     → 与 STK 预报结果对比，1 天位置误差 < 1 cm
!!   模型二：模型一 + 固体潮汐（月球 + 太阳）
!!     → 与 STK 预报结果对比，1 天位置误差 < 1 cm
!!   模型三：模型二 + 木星、火星、土星、金星第三体引力
!!     → 与 STK 预报结果对比，1 天位置误差 ≈ 1 cm
!!  模型四：模型三 + PN 引力修正
!!     → 与 STK 预报结果对比，1 天位置误差 ≈ 2 cm
!!  模型五：模型四 + Earth Radiation Pressure (ERP)
!!     → Knocke et al. (1988) 地球辐射压模型
!!
!! Dynamics:
!!   - Earth gravity field: 21×21 (GGM03C spherical harmonics)
!!   - Lunar third-body perturbation
!!   - Solar third-body perturbation
!!   - Jupiter, Mars, Saturn, Venus third-body perturbations
!!   - Earth Radiation Pressure (ERP, Knocke 1988 model)
!!   - No drag, SRP, or relativistic effects
!!
!! Output:
!!   - data/output/BDS_MEO_ephemeris.txt : Cartesian position/velocity (SI units)


program BDS_orbitPrediction
    use module_satellite
    implicit none

    ! Satellite object
    type(type_satellite) :: sat

    ! Time variables
    type(type_t) :: t0
    type(type_t), allocatable :: epochs(:)

    ! Orbital elements
    double precision :: kepler0(6)

    ! Loop index and number of output points
    integer :: i, n_points
    integer, parameter :: step_minutes = 5    ! output every 5 minutes
    double precision :: dt_minutes

    ! Output file
    character(len=100) :: output_file

    output_file = "data/output/BDS_MEO_ephemeris.txt"

    !--------------------------------------------------------------------------
    ! 1. Set initial epoch: 2025-07-28 00:00:00 UTC
    !--------------------------------------------------------------------------
    call t0%init(2025.0d0, 7.0d0, 28.0d0, 0.0d0, 0.0d0, 0.0d0)

    !--------------------------------------------------------------------------
    ! 2. Define BeiDou MEO Keplerian elements
    !--------------------------------------------------------------------------
    ! BeiDou MEO nominal orbit:
    !   Altitude ≈ 21528 km  →  a ≈ 27906 km (normalized)
    !   Near-circular (e ≈ 0), inclination ≈ 55°
    !
    ! Semi-major axis: 21528 + 6378.1363 = 27906.1363 km
    kepler0(1) = 27906.1363d0 * u_km        ! a  [km → normalized]
    kepler0(2) = 0.001d0                     ! e  (near-circular)
    kepler0(3) = 55.0d0 * deg                ! i  [deg → rad]
    kepler0(4) = 120.0d0 * deg               ! Ω  (RAAN)
    kepler0(5) = 30.0d0 * deg                ! ω  (argument of perigee)
    kepler0(6) = 0.0d0 * deg                 ! M  (mean anomaly)

    !--------------------------------------------------------------------------
    ! 3. Initialize satellite with numerical propagator
    !--------------------------------------------------------------------------
    ! Force model switches:
    !   key_moonGravity = 1 → lunar third-body
    !   key_moonTide    = 1 → lunar tide
    !   key_sunGravity  = 1 → solar third-body
    !   key_sunTide     = 1 → solar tide
    !   key_srp         = 0 → SRP (off)
    !   key_drag        = 0 → drag (off)
    !   key_postNewton  = 1 → post-Newtonian

    call sat%init(t0, kepler0,                        &
                  number      = 1,                     &
                  name        = "BeiDou-MEO-01",       &
                  mass        = 1000.0d0,              &
                  eleType     = STATE_TYPE_KEPLERIAN,  &
                  propagator  = "numerical",           &
                  degree      = 21, order = 21)           ! 21x21 gravity field
     sat%opMod%key_moonGravity = 1
     sat%opMod%key_moonTide    = 1
     sat%opMod%key_sunGravity  = 1
     sat%opMod%key_sunTide     = 1
     sat%opMod%key_postNewton  = 1

     ! Add planetary third-body perturbations
     sat%opMod%Jupiter = 1
     sat%opMod%Mars    = 1
     sat%opMod%Saturn  = 1
     sat%opMod%Venus   = 1

     ! Add Earth Radiation Pressure (ERP)
     sat%opMod%key_ERP          = 1       ! enable ERP
     sat%opMod%reflectivityERP  = 0d0   ! surface reflectivity (0~1)
     sat%opMod%amrERP           = 0.02d0  ! area-to-mass ratio [m^2/kg]

    ! Display satellite information
    call sat%show()

    !--------------------------------------------------------------------------
    ! 4. Build output epoch array: every 5 minutes for 1 day
    !--------------------------------------------------------------------------
    n_points = (24 * 60) / step_minutes + 1    ! 289 points (0 min to 1440 min)
    allocate(epochs(n_points))

    do i = 1, n_points
        dt_minutes = dble((i - 1) * step_minutes)
        epochs(i) = timeAdd(t0, dt_minutes / 1440.0d0)   ! minutes → days
    end do

    print *, "Orbit prediction configuration:"
    print *, "  Satellite:      BeiDou-MEO-01"
    print *, "  Start epoch:    ", t0%toDate()
    print *, "  End epoch:      ", epochs(n_points)%toDate()
    print *, "  Step:           ", step_minutes, " min"
    print *, "  Output points:  ", n_points
    print *, "  Gravity order:  21 × 21"
    print *, "  Perturbations:  non-spherical + luni-solar + planets (J,M,S,V) + ERP"
    print *, ""

    !--------------------------------------------------------------------------
    ! 5. Generate trajectory
    !--------------------------------------------------------------------------
    print *, "Generating trajectory..."
    call sat%generateTrajectory(epochs)

    print *, "Done. Ephemeris points: ", sat%nEph

    !--------------------------------------------------------------------------
    ! 6. Save ephemeris in Cartesian format (position + velocity, SI units)
    !--------------------------------------------------------------------------
    print *, "Saving Cartesian ephemeris to: ", trim(output_file)
    call sat%savetrajectory(STATE_TYPE_CARTESIAN, output_file, freeEph=.false.)

    print *, "Output saved successfully."

    !--------------------------------------------------------------------------
    ! Cleanup
    !--------------------------------------------------------------------------
    if (allocated(epochs)) deallocate(epochs)
    call sat%freeTrajectory()

end program BDS_orbitPrediction
