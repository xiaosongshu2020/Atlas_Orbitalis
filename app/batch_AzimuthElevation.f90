!================================================================================
! PROGRAM: batch_AzimuthElevation
!
! Description:
!   Batch orbit determination using Azimuth/Elevation (Ah) observations.
!   This program initializes satellite, station, access, and ODor objects directly
!   in code (no command-line arguments or file inputs for initialization).
!
!   Steps:
!   1. Define time series for the scenario.
!   2. Initialize satellite scenario object.
!   3. Initialize satellite object with added error in initial orbit.
!   4. Initialize ground station object and add to scenario.
!   5. Initialize access link object with OBS_TYPE_AZIMUTH_ELEVATION.
!   6. Read observation data from file (obs_LEO_Beijing_Ah.txt).
!   7. Initialize orbit determinator (ODor) with direct parameter settings.
!   8. Perform batch orbit determination.
!   9. Output results.
!
! Author: Shawn Hu (modified by assistant)
! Date: 2026-01-06
!================================================================================

 !!================================================================================
 !Main program (optional - can be commented out if called from elsewhere)
 !!================================================================================
 program  batch_AzimuthElevation

    use module_Atlas
    
    implicit none

    ! Local variables
    type(type_satelliteScenario) :: ss
    type(type_ODor) :: ODor
    type(type_t) :: t0, te
    real*8 :: interval
    integer :: nStep
    character(len=1000) :: introduction

    ! Satellite, station, and access objects
    type(type_satellite) :: sat
    type(type_station) :: station
    type(type_access) :: access

    ! Temporary variables
    type(type_t) :: sat_t0
    real*8 :: sat_x0(6), sat_x0_with_error(6)
    real*8 :: nsg1(6)
    integer :: i

    type(type_varGroup) :: Groups(1)

    ! Step 1: Define time series for the scenario
    ! Start time: 2025-01-01 00:00:00 UTC (same as satelliteScenario.f90)
    call t0%init(2025.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0)
    ! End time: 2025-01-02 00:00:00 UTC (1 day later)
    call te%init(2025.0d0, 1.0d0, 2.0d0, 0.0d0, 0.0d0, 0.0d0)
    interval = 60.0d0          ! Time step in seconds
    !nStep = 1440               ! Number of steps (1 day with 1-minute steps)
    introduction = "Batch orbit determination with Azimuth/Elevation observations"

    ! Initialize satellite scenario
    call ss%init(t0, te, interval, eleType = STATE_TYPE_NONSINGULAR_1, &
                 introduction = trim(introduction))

    ! Step 2: Initialize satellite object with added error
    sat_t0 = t0
    ! Satellite initial Keplerian elements (similar to satelliteScenario.f90 LEO orbit)
    ! a: semi-major axis (km), e: eccentricity, i: inclination (deg),
    ! Omega: RAAN (deg), omega: argument of perigee (deg), M: mean anomaly (deg)
    sat_x0(1) = 7000.0d0 * u_km        ! Semi-major axis [km -> normalized]
    sat_x0(2) = 0.01d0                 ! Eccentricity [dimensionless]
    sat_x0(3) = 70.0d0 * DEG           ! Inclination [deg -> rad]
    sat_x0(4) = 30.0d0 * DEG           ! RAAN [deg -> rad]
    sat_x0(5) = 60.0d0 * DEG           ! Argument of perigee [deg -> rad]
    sat_x0(6) = 90.0d0 * DEG           ! Mean anomaly [deg -> rad]

    ! Add error to initial orbit (5% error in semi-major axis, 1 deg error in angles)
    sat_x0_with_error = sat_x0
    sat_x0_with_error(1) = sat_x0(1) * 1.00d0          ! 1% larger semi-major axis
    sat_x0_with_error(3) = sat_x0(3) + 1.0d0 * DEG     ! 1 deg error in inclination
    sat_x0_with_error(4) = sat_x0(4) + 1.0d0 * DEG     ! 1 deg error in RAAN
    sat_x0_with_error(5) = sat_x0(5) - 1.0d0 * DEG     ! 1 deg error in argument of perigee
    sat_x0_with_error(6) = sat_x0(6) + 1.0d0 * DEG     ! 2 deg error in mean anomaly
    
    call kpl2nsg1(sat_x0_with_error, nsg1)
    
    call sat%init(sat_t0, nsg1, &
                  number = 1, &
                  name = 'LEO-Sat-1', &
                  eleType = STATE_TYPE_NONSINGULAR_1, &
                  propagator = 'numerical', &
                  mass = 500.0d0, &
                  degree=8, order=8, & ! Gravity field degree/order
                  Cd=2.2d0, &          ! Drag coefficient
                  amrDrag=0.02d0, &    ! Area-to-mass ratio for drag
                  F107=130.0d0, &      ! Solar flux index
                  Ap=15.0d0)           ! Geomagnetic index
    sat%opMod%key_moonGravity = 1   ! Moon third-body
    sat%opMod%key_sunGravity  = 1   ! Sun third-body
    sat%opMod%key_srp         = 1   ! solar radiation pressure
    sat%opMod%key_drag        = 1   ! atmospheric drag

    ! Add satellite to scenario
    call ss%addSatellite(sat)

    ! Step 3: Initialize ground station object (Beijing station)
    call station%init(number = 1, &
                      name = 'Beijing', &
                      longitude = 116.4d0, &
                      latitude = 39.9d0, &
                      altitude = 50.0d0)

    ! Generate station's trajectory in J2000 frame for the scenario time series
    call station%generateTrajectory(ss%ts)

    ! Add ground station to scenario
    call ss%addStation(station)

    ! Step 4: Initialize access link object with Azimuth/Elevation observation type
    call access%init(number = 1, &
                     name = 'LEO-Beijing-Ah', &
                     fromType = OBJECT_TYPE_SATELLITE, &
                     fromNum = 1, &
                     toType = OBJECT_TYPE_STATION, &
                     toNum = 1, &
                     obsType = OBS_TYPE_AZIMUTH_ELEVATION, &
                     minEle = 10.0d0, &
                     solarDistance = 10.0d0, &
                     lunarDistance = 10.0d0)

    ! Set observation error standard deviation (0.1 deg ~ 6 arcmin)
    access%errStd(1:2) = 0.1d0 * DEG
    ! Enable light-time delay correction (optional)
    key_light_time_delay = .false.
    
    ! Read observation data from file generated by satelliteScenario.f90
    call access%readData('data/output/obs_LEO_Beijing_Ah.txt')

    ! Add access link to scenario
    call ss%addAccess(access)

    ! Report access relationships
    call ss%reportAccess(.true.)

    ! Step 5: Initialize orbit determinator (ODor)
    ! Set up variable group parameters
    Groups(1)%nVar = 6                     ! 6 state variables (Keplerian elements)
    Groups(1)%label = 1                    ! Label 1: states only
    Groups(1)%origin = 1                   ! Belongs to satellite 1
    ! Allocate and set indices
    allocate(Groups(1)%idx(6))
    Groups(1)%idx = (/1, 2, 3, 4, 5, 6/)   ! Indices for the 6 state variables
    
    ! Initialize ODor with batch method
    call ODor%init(ss, ss%ts%t0, 1, Groups, OD_METHOD_BATCH, &
                   intro = "Batch orbit determination with Azimuth/Elevation observations")
    
    ! Step 6: Perform batch orbit determination
    call ODor%solve(ss)

    ! Step 7: Output results
    call ODor%saveEph(ss, "data/output/preciseEph_LEO_Ah_OD.txt")

    ! Step 8: Clean up
    call ss%freeSS
    deallocate(Groups(1)%idx)

    print *, ""
    print *, "================================================================"
    print *, "BATCH ORBIT DETERMINATION WITH AZIMUTH/ELEVATION COMPLETED"
    print *, "================================================================"
    print *, ""

end program batch_AzimuthElevation
