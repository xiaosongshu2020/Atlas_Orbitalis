!**************************************************************************
! Program Name: app_observation
! Purpose: Observation data generation demo - all observation types
! Function: Build one LEO satellite and one ground station scenario and
!           generate every observation type supported by module_access:
!             RA_DEC, AZIMUTH_ELEVATION, RANGING, RANGERATE,
!             POSITION, POSITION_VELOCITY
!           Each link is written to data/output/demo_observation/.
!
! All scenario settings are hard-coded below (modify and recompile to
! change them), following the batch-program convention. Run the program
! from the repository root so the relative data paths resolve.
!**************************************************************************
program app_observation
use module_const              ! Constants module
use module_time               ! Time processing module
use module_timeSeries
use module_satelliteScenario  ! Satellite scenario module
use module_satellite          ! Satellites module
use module_station            ! Station module
use module_access             ! Observation types and generation routines
IMPLICIT NONE

!---------------------- Scenario Setting Section (hard-coded) ----------------------
! Time settings
REAL(8), parameter :: step = 30d0              ! time step (seconds)
! Satellite settings: ~600 km altitude, inclination 60 deg, near-circular LEO
REAL(8), parameter :: altKm = 600d0            ! satellite altitude (km)
REAL(8), parameter :: inclination = 60d0       ! inclination (degrees)
REAL(8), parameter :: amrDrag = 0.02d0         ! drag area-to-mass ratio (m^2/kg)
! Ground station settings: Beijing
REAL(8), parameter :: longitude = 116.4d0      ! station longitude (degrees)
REAL(8), parameter :: latitude = 39.9d0        ! station latitude (degrees)
REAL(8), parameter :: altitude = 50d0          ! station altitude (meters)
! Observation settings
REAL(8), parameter :: minEle = 10d0            ! minimum elevation (degrees)
CHARACTER(len=*), parameter :: outDir = "data/output/demo_observation/"
CHARACTER(len=*), parameter :: satName = "LEO-600km"
! The six supported observation types, in the order they are demonstrated
INTEGER, parameter :: obsTypes(6) = (/ OBS_TYPE_RA_DEC, OBS_TYPE_AZIMUTH_ELEVATION, &
    &   OBS_TYPE_RANGING, OBS_TYPE_RANGERATE, OBS_TYPE_POSITION, &
    &   OBS_TYPE_POSITION_VELOCITY /)
CHARACTER(len=32), parameter :: obsNames(6) = (/ character(len=32) :: &
    &   "RA_DEC", "AZIMUTH_ELEVATION", "RANGING", "RANGERATE", &
    &   "POSITION", "POSITION_VELOCITY" /)

!---------------------- Variable Declaration Section ----------------------
type(type_t) :: t0, te
type(type_satelliteScenario) :: ss
type(type_satellite) :: sat
type(type_station) :: sta
type(type_access) :: ac
REAL(8) :: kpl(6)                         ! Keplerian elements (normalized units)
type(type_t), allocatable :: epochs(:)    ! trajectory output epochs
CHARACTER(len=160) :: fname
LOGICAL :: dirExists
INTEGER :: i, k

!---------------------- Time Setup Section ----------------------
call t0%init(2025d0, 1d0, 1d0, 0d0, 0d0, 0d0)   ! 2025-01-01 00:00:00 UTC
te = timeAdd(t0, 1d0)                            ! 1 day arc

!---------------------- Output Directory ----------------------
! saveAccessData opens the file directly and aborts if the directory is
! missing, so create it first. Skip the shell call when it already exists:
! mkdir on an existing directory prints a localized OS error.
INQUIRE(FILE="data/output/demo_observation", EXIST=dirExists)
IF (.NOT. dirExists) THEN
    call execute_command_line('mkdir "data/output/demo_observation"', wait=.true.)
END IF

!---------------------- Scenario Establishment Section ----------------------
WRITE(*,*) "Observation demo: all six observation types"
WRITE(*,*) "  Satellite: altitude=", altKm, "km, inclination=", inclination, "deg"
WRITE(*,*) "  Station:   longitude=", longitude, "deg, latitude=", latitude, "deg"
WRITE(*,*) "  Step:      ", step, "s, arc: 1 day"
WRITE(*,*) "  Output:    ", TRIM(outDir)

call ss%init(t0, te, step, STATE_TYPE_CARTESIAN, &
    & "Observation demo with all six observation types")

! Initialize satellite (Keplerian elements in normalized units)
kpl(1) = (LENGTH_UNIT + altKm) * u_km     ! semi-major axis
kpl(2) = 0.001d0                          ! eccentricity
kpl(3) = inclination * deg                ! inclination
kpl(4) = 0d0                              ! RAAN
kpl(5) = 0d0                              ! argument of perigee
kpl(6) = 0d0                              ! mean anomaly

call sat%init(t0, kpl, number=1, name=trim(satName), mass=500d0, &
    &   eleType=STATE_TYPE_KEPLERIAN, propagator="numerical", &
    &   degree=8, order=8, amrDrag=amrDrag)
sat%opMod%key_moonGravity = 1   ! Moon third-body
sat%opMod%key_sunGravity  = 1   ! Sun third-body
sat%opMod%key_drag        = 1   ! atmospheric drag
call ss%addSatellite(sat)

! Initialize ground station
call sta%init(number=1, name="Beijing", longitude=longitude, &
    &   latitude=latitude, altitude=altitude)
call ss%addStation(sta)

! One access per observation type. initAccess does not reset the constraint
! fields, so shadow is assigned explicitly on every pass: optical (RA_DEC,
! AZIMUTH_ELEVATION) skips epochs in the Earth shadow, radio link types do not.
do k = 1, 6
    call ac%init(number=k, name=trim(obsNames(k)), &
        &   fromType=OBJECT_TYPE_SATELLITE, fromNum=1, &
        &   toType=OBJECT_TYPE_STATION, toNum=1, &
        &   obsType=obsTypes(k), minEle=minEle)
    select case(obsTypes(k))
    case(OBS_TYPE_RA_DEC, OBS_TYPE_AZIMUTH_ELEVATION)
        ac%shadow = .true.
    case default
        ac%shadow = .false.
    end select
    call ss%addAccess(ac)
end do

!---------------------- Trajectory Generation Section ----------------------
allocate(epochs(ss%ts%n))
do i = 1, ss%ts%n
    epochs(i) = ss%ts%time(i)
end do

WRITE(*,*) "Generating satellite trajectory..."
call ss%sat(1)%generateTrajectory(epochs)
WRITE(*,*) "  Satellite trajectory points: ", ss%sat(1)%nEph

WRITE(*,*) "Generating station trajectory..."
call ss%station(1)%generateTrajectory(ss%ts)

!---------------------- Observation Calculation Section ----------------------
WRITE(*,*) "Calculating observations for all six types..."
call ss%observation
call ss%reportAccess

!---------------------- Save Section ----------------------
WRITE(*,*) "Saving observation data:"
do k = 1, ss%nAccess
    fname = TRIM(outDir)//"obs_"//TRIM(obsNames(k))//".txt"
    call ss%access(k)%saveData(ss%ts, TRIM(fname))
    WRITE(*, "(A,I0,A)") "  saved: "//TRIM(fname)//"   (epochs: ", &
        &   ss%access(k)%nData, ")"
end do

! Save the satellite ephemeris used to generate the observations
fname = TRIM(outDir)//"eph_"//TRIM(satName)//".txt"
call ss%sat(1)%savetrajectory(STATE_TYPE_CARTESIAN, TRIM(fname))
WRITE(*,*) "  saved: ", TRIM(fname)

!---------------------- Cleanup Section ----------------------
deallocate(epochs)
call ss%freeSS

!---------------------- Program End Section ----------------------
WRITE(*,*) "Processing complete!"

END PROGRAM app_observation
