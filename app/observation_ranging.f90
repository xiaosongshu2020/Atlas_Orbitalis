

!**************************************************************************
! Program Name: observation_ranging
! Purpose: Satellite ranging observation simulation
! Function: Build a scenario with one LEO satellite (~600 km altitude) and
!           one ground station, then generate range-only observation data
!           (obsType = RANGING, 13).
!
! All scenario settings are hard-coded below (modify and recompile to
! change them), following the batch-program convention.
!**************************************************************************
program observation_ranging
use module_const            ! Constants module
use module_time             ! Time processing module
use module_timeSeries
use module_satelliteScenario  ! Satellite scenario module
use module_satellite        ! Satellites module
use module_station          ! Station module
use module_access           ! Visibility calculation module
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
CHARACTER(len=*), parameter :: fnData = "data/output/obs_ranging.txt"
! Clock settings: the satellite (transmitter) carries a non-zero clock so the
! clock-error correction in calObsData1 is actually exercised. The ground
! station (receiver) keeps the default clock, i.e. zero bias and zero drift,
! so the observation epochs are left untouched and only the range carries the
! correction:  rho_obs = rho_true - clock_error(satellite)*c
REAL(8), parameter :: clkBias0 = 1.0d-4        ! satellite clock bias at t0 [s]
REAL(8), parameter :: clkDrift = 1.0d-11       ! satellite clock drift rate [s/s]

!---------------------- Variable Declaration Section ----------------------
! Scenario objects
type(type_t) :: t0, te
type(type_satelliteScenario) :: ss
type(type_satellite) :: sat
type(type_station) :: sta
type(type_access) :: ac
REAL(8) :: kpl(6)            ! Keplerian elements (normalized units)
type(type_t), allocatable :: epochs(:)  ! Trajectory output epochs

! Auxiliary variables
INTEGER :: i

!---------------------- Time Setup Section ----------------------
call t0%init(2025d0, 1d0, 1d0, 0d0, 0d0, 0d0)   ! 2025-01-01 00:00:00 UTC
te = timeAdd(t0, 1d0)                            ! 1 day arc

!---------------------- Scenario Establishment Section ----------------------
WRITE(*,*) "Ranging observation scenario:"
WRITE(*,*) "  Satellite: altitude=", altKm, "km, inclination=", inclination, "deg"
WRITE(*,*) "  Station: longitude=", longitude, "deg, latitude=", latitude, "deg, altitude=", altitude, "m"
WRITE(*,*) "  Minimum elevation: ", minEle, "deg"
WRITE(*,*) "  Satellite clock: bias0=", clkBias0, "s, drift=", clkDrift, "s/s"
WRITE(*,*) "  Station clock:   default (zero bias/drift)"
WRITE(*,*) "  Output file: ", TRIM(fnData)

call ss%init(t0, te, step, STATE_TYPE_CARTESIAN, &
    & "Ranging observation: one LEO satellite and one ground station")

! Initialize satellite (Keplerian elements in normalized units)
kpl(1) = (LENGTH_UNIT + altKm) * u_km     ! semi-major axis
kpl(2) = 0.001d0                          ! eccentricity
kpl(3) = inclination * deg                ! inclination
kpl(4) = 0d0                              ! RAAN
kpl(5) = 0d0                              ! argument of perigee
kpl(6) = 0d0                              ! mean anomaly

call sat%init(t0, kpl, number=1, name="LEO-600km", mass=500d0, &
    & eleType=STATE_TYPE_KEPLERIAN, propagator="numerical", &
    & degree=8, order=8, amrDrag=amrDrag)   ! gravity + luni-solar + drag
sat%opMod%key_moonGravity = 1   ! Moon third-body
sat%opMod%key_sunGravity  = 1   ! Sun third-body
sat%opMod%key_drag        = 1   ! atmospheric drag

! Satellite clock model. The bias and drift are referred to t0, so the bias
! at an epoch t is  bias0 + drift*(t - t0). Set it before addSatellite: the
! scenario stores a by-value copy of the whole satellite, clock included.
! Default calcMethod is the analytic formula, so no clock table is needed.
call sat%clock%init(t0, bias0=clkBias0, drift=clkDrift, drift2=0d0)

call ss%addSatellite(sat)

! Initialize ground station
call sta%init(number=1, name="station1", longitude=longitude, latitude=latitude, altitude=altitude)
call ss%addStation(sta)

! Initialize access relationship: satellite -> station, range only
call ac%init(number=1, name="ranging1", &
    & fromType=OBJECT_TYPE_SATELLITE, fromNum=1, toType=OBJECT_TYPE_STATION, toNum=1, &
    & obsType=OBS_TYPE_RANGING, minEle=minEle, solarDistance=0d0, lunarDistance=0d0)
call ss%addAccess(ac)

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
WRITE(*,*) "Calculating ranging observations..."
call ss%observation
call ss%reportAccess

! Save observation data to file
call ss%access(1)%saveData(ss%ts, trim(fnData))

!---------------------- Cleanup Section ----------------------
if (allocated(epochs)) deallocate(epochs)
call ss%freeSS

!---------------------- Program End Section ----------------------
WRITE(*,*) "Processing complete!"

END PROGRAM observation_ranging
