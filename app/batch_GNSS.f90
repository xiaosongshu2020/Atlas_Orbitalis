!================================================================================
! PROGRAM: batch_GNSS
!
! Description:
!   Batch orbit determination using GNSS observations (position / velocity).
!   This program initializes satellite, station, access, and ODor objects directly
!   in code (no command-line arguments or file inputs for initialization).
!
!   Steps:
!   1. Define time series for the scenario.
!   2. Initialize satellite scenario object.
!   3. Initialize satellite object and add to scenario.
!   4. Initialize ground station object and add to scenario.
!   5. Initialize access link object and add to scenario.
!   6. Initialize orbit determinator (ODor) with direct parameter settings.
!   7. Read observation data from a file (optional, can be simulated).
!   8. Perform batch orbit determination.
!   9. Output results.
!
! Author: Shawn Hu
! Date: 2025-12-24
!================================================================================

!================================================================================
! Main program (optional - can be commented out if called from elsewhere)
!================================================================================
program  batch_GNSS

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
    real*8 :: sat_x0(6)
    integer :: i

    type(type_varGroup) :: Groups(1)

    ! Step 1: Define time series for the scenario
    ! Start time: 
    call t0%init(57445d0, 57624d0)
    ! End time: 2 days later
    call te%init(57447d0, 57624d0)
    interval = 60d0          ! Time step in seconds
    nStep = 2 * 1440    ! Number of steps
    introduction = "Batch orbit determination with angle observations (RA/DEC)"

    ! Initialize satellite scenario
    call ss%init(t0, te, interval, eleType = STATE_TYPE_NONSINGULAR_1, &
                 introduction = trim(introduction))

    ! Step 2: Initialize satellite object
    sat_t0 = t0
    ! Satellite initial state (Cartesian coordinates, km and km/s)
    sat_x0(1:3) = (/-779.0d0, 1394d0, -7176d0/) * u_km
    sat_x0(4:6) = (/2.54d0, 6.81d0, 1.04d0/) * u_km / u_sec
    call cts2nsg1(sat_x0, sat_x0) ! Convert Cartesian to nonsingular elements

    call sat%init(sat_t0, sat_x0, &
                  number = 1, &
                  name = 'HY-2A', &
                  eleType = STATE_TYPE_NONSINGULAR_1, &
                  propagator = 'numerical', &
                  Cd = 2.2d0, &
                  amrDrag = 0.01d0, &
                  F107 = 130d0, &
                  Ap = 15d0, &
                  reflectivity = 0.3d0, &
                  amrSRP = 0.02d0, &
                  degree = 10, order = 10)
    sat%opMod%key_moonGravity = 1   ! Moon third-body
    sat%opMod%key_moonTide    = 1   ! Moon solid tide
    sat%opMod%key_sunGravity  = 1   ! Sun third-body
    sat%opMod%key_sunTide     = 1   ! Sun solid tide
    sat%opMod%key_srp         = 1   ! solar radiation pressure
    sat%opMod%key_drag        = 1   ! atmospheric drag
    sat%opMod%key_postNewton  = 1   ! post-Newtonian
    sat%opMod%usePV = 0   ! Use Cartesian elements for propagation

    ! Add satellite to scenario
    call ss%addSatellite(sat)

    ! Step 3: Initialize ground station object
    call station%init(number = 1, &
                      name = 'Kashi', &
                      longitude = 75.9797d0, &
                      latitude = 39.4547d0, &
                      altitude = 0.0d0)

    ! Generate station's trajectory in J2000 frame for the scenario time series
    call station%generateTrajectory(ss%ts)

    ! Add ground station to scenario
    call ss%addStation(station)

    ! Step 4: Initialize access link object
    call access%init(number = 1, &
                     name = 'access_HY-2A_Kashi', &
                     fromType = OBJECT_TYPE_SATELLITE, &
                     fromNum = 1, &
                     toType = OBJECT_TYPE_STATION, &
                     toNum = 1, &
                     obsType = OBS_TYPE_GNSS, &
                     minEle = 10d0)

    ! Set observation error standard deviation (5e-6 rad ~ 1 arcsec)
    access%errStd = 5d-6
    ! Enable light-time delay correction
    key_light_time_delay = .false.
    call access%readData('data/input/eph_fit8h_refer.txt')

    ! Add access link to scenario
    call ss%addAccess(access)

    ! Report access relationships
    call ss%reportAccess(.true.)

    ! Step 5: Initialize orbit determinator (ODor)
    ! Set up variable group parameters
    Groups(1)%nVar = 6                     ! 6 state variables (position and velocity)
    Groups(1)%label = 1                    ! Label 1: states only
    Groups(1)%origin = 1                   ! Belongs to satellite 1
    ! Allocate and set indices
    allocate(Groups(1)%idx(6))
    Groups(1)%idx = (/1, 2, 3, 4, 5, 6/)   ! Indices for the 6 state variables
    
    ! Initialize ODor with batch method
    call ODor%init(ss, ss%ts%t0, 1, Groups, OD_METHOD_BATCH, &
                   intro = "Batch orbit determination with GNSS observations")
    
    ! Step 7: Perform batch orbit determination
    call ODor%solve(ss)

    ! Step 8: Output results
    call ODor%saveEph(ss, "data/output/preciseEph_HY-2A.txt")

    ! Step 9: Clean up
    call ss%freeSS
    deallocate(Groups(1)%idx)

end program batch_GNSS



