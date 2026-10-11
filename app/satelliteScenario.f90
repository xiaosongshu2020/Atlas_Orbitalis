!! ============================================================================
!! DEMONSTRATION PROGRAM: satelliteScenario.f90
!!
!! DESCRIPTION:
!!   This example demonstrates the usage of the satellite scenario object
!!   (type_satelliteScenario) for managing satellite simulation scenarios.
!!   It shows how to create a scenario with multiple satellites and ground
!!   stations, define access links between them, perform observation calculations,
!!   and report the results.
!!
!! KEY FEATURES DEMONSTRATED:
!!   1. Satellite scenario initialization with time window
!!   2. Adding satellite objects to the scenario
!!   3. Adding ground station objects to the scenario
!!   4. Creating and adding access links with different observation types
!!   5. Performing observation calculations for all access links
!!   6. Reporting access link information and observation statistics
!!
!! DEPENDENCIES:
!!   - module_satelliteScenario
!!   - module_satellite
!!   - module_station
!!   - module_access
!!   - module_time
!!
!! AUTHOR: Shawn Hu
!! DATE: 2025-12-24
!!
!! LAST MODIFIED: 2025-12-24
!! ============================================================================


PROGRAM demo_satelliteScenario
  use module_satelliteScenario
  implicit none
  
  ! -------------------------------------------------------------------------
  ! DECLARATIONS
  ! -------------------------------------------------------------------------
  
  ! Time objects
  type(type_t) :: t0, te
  ! Satellite objects
  type(type_satellite) :: sat1, sat2
  ! Ground station objects
  type(type_station) :: station1, station2
  ! Access objects
  type(type_access) :: access1, access2, access3
  ! Satellite scenario object
  type(type_satelliteScenario) :: scenario
  ! Keplerian elements for satellites (a in km, angles in degrees)
  real*8 :: kepler1(6), kepler2(6)
  ! Output file paths
  character(len=100) :: output_file1, output_file2, output_file3
  ! Other variables
  integer :: i, n_points
  type(type_t), allocatable :: epochs(:)
  
  ! -------------------------------------------------------------------------
  ! INITIALIZATION
  ! -------------------------------------------------------------------------
  
  print *, ""
  print *, "================================================================"
  print *, "SATELLITE SCENARIO DEMONSTRATION"
  print *, "================================================================"
  print *, ""
  
  ! -------------------------------------------------------------------------
  ! EXAMPLE 1: SET UP TIME WINDOW
  ! -------------------------------------------------------------------------
  
  ! Create start time: 2025-01-01 00:00:00 UTC
  call t0%init(2025.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0)
  
  ! Create end time: 2025-01-02 00:00:00 UTC
  call te%init(2025.0d0, 1.0d0, 2.0d0, 0.0d0, 0.0d0, 0.0d0)
  
  print *, "Time window set:"
  print *, "  Start: ", t0%toDate()
  print *, "  End:   ", te%toDate()
  print *, ""
  
  ! -------------------------------------------------------------------------
  ! EXAMPLE 2: INITIALIZE SATELLITE SCENARIO
  ! -------------------------------------------------------------------------
  
  ! Initialize satellite scenario with 60-second time steps
  call scenario%init(t0, te, 60.0d0, &
                     introduction="Demo scenario with 2 satellites and 2 stations")
  
  ! -------------------------------------------------------------------------
  ! EXAMPLE 3: CREATE AND ADD SATELLITES
  ! -------------------------------------------------------------------------
  
  ! Define Keplerian elements for satellite 1 (LEO orbit)
  kepler1(1) = 7000.0d0 * u_km        ! Semi-major axis [km -> normalized]
  kepler1(2) = 0.01d0                 ! Eccentricity [dimensionless]
  kepler1(3) = 70.0d0 * deg           ! Inclination [deg -> rad]
  kepler1(4) = 30.0d0 * deg           ! RAAN [deg -> rad]
  kepler1(5) = 60.0d0 * deg           ! Argument of perigee [deg -> rad]
  kepler1(6) = 90.0d0 * deg           ! Mean anomaly [deg -> rad]
  
  ! Initialize satellite 1 with numerical propagator
  call sat1%init(t0, kepler1, &
                 number=1, &
                 name="LEO-Sat-1", &
                 mass=500.0d0, &
                 eleType=STATE_TYPE_KEPLERIAN, &
                 propagator="numerical", &
                 degree=8, order=8, & ! Gravity field degree/order
                 Cd=2.2d0, &          ! Drag coefficient
                 amrDrag=0.02d0, &    ! Area-to-mass ratio for drag
                 F107=130.0d0, &      ! Solar flux index
                 Ap=15.0d0)           ! Geomagnetic index
  sat1%opMod%key_moonGravity = 1   ! Moon third-body
  sat1%opMod%key_sunGravity  = 1   ! Sun third-body
  sat1%opMod%key_srp         = 1   ! solar radiation pressure
  sat1%opMod%key_drag        = 1   ! atmospheric drag
  
  ! Define Keplerian elements for satellite 2 (higher orbit)
  kepler2(1) = 26560.0d0 * u_km      ! Semi-major axis [km -> normalized] ~GPS orbit
  kepler2(2) = 0.001d0               ! Eccentricity [dimensionless]
  kepler2(3) = 55.0d0 * deg          ! Inclination [deg -> rad]
  kepler2(4) = 120.0d0 * deg         ! RAAN [deg -> rad]
  kepler2(5) = 30.0d0 * deg          ! Argument of perigee [deg -> rad]
  kepler2(6) = 180.0d0 * deg         ! Mean anomaly [deg -> rad]
  
  ! Initialize satellite 2 with two-body propagator
  call sat2%init(t0, kepler2, &
                 number=2, &
                 name="GPS-Sat-2", &
                 mass=1000.0d0, &
                 eleType=STATE_TYPE_KEPLERIAN, &
                 propagator="two-body")
  
  ! Add satellites to the scenario
  call scenario%addSatellite(sat1)
  call scenario%addSatellite(sat2)
  
  print *, "Satellites added to scenario:"
  print *, "  1. ", trim(sat1%name)
  print *, "  2. ", trim(sat2%name)
  print *, ""
  
  ! -------------------------------------------------------------------------
  ! EXAMPLE 4: CREATE AND ADD GROUND STATIONS
  ! -------------------------------------------------------------------------
  
  ! Initialize ground station 1 (Beijing)
  call station1%init(number=1, name="Beijing", &
                     longitude=116.4d0, latitude=39.9d0, altitude=50.0d0)
  
  ! Initialize ground station 2 (Tokyo)
  call station2%init(number=2, name="Tokyo", &
                     longitude=139.76d0, latitude=35.68d0, altitude=40.0d0)
  
  ! Add stations to the scenario
  call scenario%addStation(station1)
  call scenario%addStation(station2)
  
  print *, "Ground stations added to scenario:"
  print *, "  1. ", trim(station1%name)
  print *, "  2. ", trim(station2%name)
  print *, ""
  
  ! -------------------------------------------------------------------------
  ! EXAMPLE 5: CREATE AND ADD ACCESS LINKS
  ! -------------------------------------------------------------------------
  
  ! Access 1: Satellite 1 to Station 1 (Azimuth/Elevation observations)
  call access1%init(number=1, name="LEO-Beijing-Ah", &
                   fromType=OBJECT_TYPE_SATELLITE, fromNum=1, &
                   toType=OBJECT_TYPE_STATION, toNum=1, &
                   obsType=OBS_TYPE_AZIMUTH_ELEVATION, &
                   minEle=10.0d0, &          ! Minimum elevation angle
                   solarDistance=10.0d0, &   ! Minimum angular distance from Sun
                   lunarDistance=10.0d0)     ! Minimum angular distance from Moon
  
  ! Access 2: Satellite 1 to Station 2 (Range and Range-Rate observations)
  call access2%init(number=2, name="LEO-Tokyo-RRD", &
                   fromType=OBJECT_TYPE_SATELLITE, fromNum=1, &
                   toType=OBJECT_TYPE_STATION, toNum=2, &
                   obsType=OBS_TYPE_RANGE_DOPPLER, &
                   minEle=10.0d0, &          ! Minimum elevation angle
                   solarDistance=10.0d0, &   ! Minimum angular distance from Sun
                   lunarDistance=10.0d0)     ! Minimum angular distance from Moon
  
  ! Access 3: Satellite 2 to Station 1 (GNSS position/velocity observations)
  call access3%init(number=3, name="GPS-Beijing-GNSS", &
                   fromType=OBJECT_TYPE_SATELLITE, fromNum=2, &
                   toType=OBJECT_TYPE_STATION, toNum=1, &
                   obsType=OBS_TYPE_GNSS, &
                   minEle=5.0d0)             ! Minimum elevation angle
  
  ! Add access links to the scenario
  call scenario%addAccess(access1)
  call scenario%addAccess(access2)
  call scenario%addAccess(access3)
  
  print *, "Access links created:"
  print *, "  1. ", trim(access1%name), " (Azimuth/Elevation)"
  print *, "  2. ", trim(access2%name), " (Range/Range-Rate)"
  print *, "  3. ", trim(access3%name), " (GNSS position/velocity)"
  print *, ""
  
  ! -------------------------------------------------------------------------
  ! EXAMPLE 6: GENERATE TRAJECTORIES FOR OBSERVATION
  ! -------------------------------------------------------------------------
  
  ! Create array of output epochs for trajectory generation
  n_points = scenario%ts%n
  allocate(epochs(n_points))
  
  do i = 1, n_points
    epochs(i) = scenario%ts%time(i)
  end do
  
  ! Generate trajectories for satellites
  print *, "Generating satellite trajectories..."
  call scenario%sat(1)%generateTrajectory(epochs)
  call scenario%sat(2)%generateTrajectory(epochs)
  print *, "  LEO-Sat-1 trajectory points: ", scenario%sat(1)%nEph
  print *, "  GPS-Sat-2 trajectory points: ", scenario%sat(2)%nEph
  print *, ""
  
  ! Generate trajectories for stations
  print *, "Generating station trajectories..."
  call scenario%station(1)%generateTrajectory(scenario%ts)
  call scenario%station(2)%generateTrajectory(scenario%ts)
  print *, "  Station trajectories generated."
  print *, ""
  
  ! -------------------------------------------------------------------------
  ! EXAMPLE 7: PERFORM OBSERVATION CALCULATIONS
  ! -------------------------------------------------------------------------
  
  print *, "Performing observation calculations for all access links..."
  print *, "  This may take a moment..."
  call scenario%observation()
  print *, "  Observation calculations completed."
  print *, ""
  
  ! -------------------------------------------------------------------------
  ! EXAMPLE 8: REPORT ACCESS LINK INFORMATION
  ! -------------------------------------------------------------------------
  
  print *, "Access link report:"
  print *, "==================="
  call scenario%reportAccess()
  
  ! -------------------------------------------------------------------------
  ! EXAMPLE 9: SAVE OBSERVATION DATA TO FILES (OPTIONAL)
  ! -------------------------------------------------------------------------
  
  ! Set output file paths
  output_file1 = "data/output/obs_LEO_Beijing_Ah.txt"
  output_file2 = "data/output/obs_LEO_Tokyo_RRD.txt"
  output_file3 = "data/output/obs_GPS_Beijing_GNSS.txt"
  
  print *, "Saving observation data to files:"
  
  if (scenario%access(1)%nData > 0) then
    call scenario%access(1)%saveData(scenario%ts, output_file1)
    print *, "  Access 1 data saved to: ", trim(output_file1)
  end if
  
  if (scenario%access(2)%nData > 0) then
    call scenario%access(2)%saveData(scenario%ts, output_file2)
    print *, "  Access 2 data saved to: ", trim(output_file2)
  end if
  
  if (scenario%access(3)%nData > 0) then
    call scenario%access(3)%saveData(scenario%ts, output_file3)
    print *, "  Access 3 data saved to: ", trim(output_file3)
  end if
  
  print *, ""
  
  ! -------------------------------------------------------------------------
  ! CLEANUP
  ! -------------------------------------------------------------------------
  
  ! Free scenario memory (this will deallocate all satellite, station, and access arrays)
  call scenario%freeSS()
  
  ! Deallocate arrays
  if (allocated(epochs)) deallocate(epochs)
  
  ! -------------------------------------------------------------------------
  ! COMPLETION
  ! -------------------------------------------------------------------------
  
  print *, "================================================================"
  print *, "DEMONSTRATION COMPLETED SUCCESSFULLY"
  print *, "================================================================"
  print *, ""
  
end program demo_satelliteScenario
