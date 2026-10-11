

!**************************************************************************
! Program Name: SatelliteObservation
! Purpose: Satellite observation scenario simulation program
! Function: Calculate satellite observation data based on given station position and parameters
!**************************************************************************
program SatelliteObservation
use module_const      ! Constants module
use module_time       ! Time processing module
use module_timeSeries
use module_satelliteScenario  ! Satellite scenario module
use module_satellite  ! Satellites module
use module_station    ! Station module
use module_access     ! Visibility calculation module
IMPLICIT NONE

!---------------------- Variable Declaration Section ----------------------
! Basic variables
INTEGER :: i, stat, nArgs
CHARACTER(len=256) :: arg

! Time parameters
REAL(8) :: ts, step, te  ! Start time(MJD), time step(seconds), end time(MJD)

! Observer station parameters
REAL(8) :: longitude, latitude, altitude  ! Longitude, latitude(degrees), altitude(meters)

! Observation condition parameters
INTEGER :: obsType                        ! Observation type
REAL(8) :: minEle                         ! Minimum elevation angle(degrees)
REAL(8) :: solarDistance                  ! Minimum angular distance from Sun(degrees)
REAL(8) :: lunarDistance                  ! Minimum angular distance from Moon(degrees)

! File parameters
CHARACTER(len=256) :: fnEph, fnData       ! trajectory file, data output file
LOGICAL :: hasDataFile                    ! Flag indicating if data file is provided

! Time series and scenario parameters
type(type_t) :: t0                        ! Initial time
integer :: nStep                          ! Total steps
type(type_timeSeries) :: timeS            ! Time series

type(type_satelliteScenario) :: ss

!---------------------- Initialization Section ----------------------
! Initialize optional parameter flags
hasDataFile = .FALSE.

! Get number of command line arguments
nArgs = COMMAND_ARGUMENT_COUNT()

!---------------------- Parameter Check Section ----------------------
! Check if sufficient arguments are provided
IF (nArgs < 11) THEN
    WRITE(*,*) "Usage: ./program ts step te longitude latitude altitude obsType minEle solarDistance lunarDistance fnEph [fnData]"
    WRITE(*,*) "Parameters:"
    WRITE(*,*) "  ts - Initial time (MJD)"
    WRITE(*,*) "  step - Time step (seconds)"
    WRITE(*,*) "  te - End time (MJD)"
    WRITE(*,*) "  longitude - Observer longitude (degrees)"
    WRITE(*,*) "  latitude - Observer latitude (degrees)"
    WRITE(*,*) "  altitude - Observer altitude (meters)"
    WRITE(*,*) "  obsType - Observation type (integer)"
    WRITE(*,*) "  minEle - Minimum elevation angle (degrees)"
    WRITE(*,*) "  solarDistance - Minimum angular distance from Sun (degrees)"
    WRITE(*,*) "  lunarDistance - Minimum angular distance from Moon (degrees)"
    WRITE(*,*) "  fnEph - Satellite trajectory filename"
    WRITE(*,*) "  fnData - Observation data output filename (optional)"
    WRITE(*,*) "Example: ./observation 60808 30 60809 108 30 20 16 10 0 0 ./satEph.txt ./obsData.txt"
    STOP
END IF

!---------------------- Parameter Reading Section ----------------------
! Read start time
CALL GET_COMMAND_ARGUMENT(1, arg)
READ(arg, *, IOSTAT=stat) ts
IF (stat /= 0) THEN
    WRITE(*,*) "Error: Unable to parse initial time parameter"
    STOP
END IF

! Read time step
CALL GET_COMMAND_ARGUMENT(2, arg)
READ(arg, *, IOSTAT=stat) step
IF (stat /= 0) THEN
    WRITE(*,*) "Error: Unable to parse time step parameter"
    STOP
END IF

! Read end time
CALL GET_COMMAND_ARGUMENT(3, arg)
READ(arg, *, IOSTAT=stat) te
IF (stat /= 0) THEN
    WRITE(*,*) "Error: Unable to parse end time parameter"
    STOP
END IF

! Read longitude
CALL GET_COMMAND_ARGUMENT(4, arg)
READ(arg, *, IOSTAT=stat) longitude
IF (stat /= 0) THEN
    WRITE(*,*) "Error: Unable to parse longitude parameter"
    STOP
END IF

! Read latitude
CALL GET_COMMAND_ARGUMENT(5, arg)
READ(arg, *, IOSTAT=stat) latitude
IF (stat /= 0) THEN
    WRITE(*,*) "Error: Unable to parse latitude parameter"
    STOP
END IF

! Read altitude
CALL GET_COMMAND_ARGUMENT(6, arg)
READ(arg, *, IOSTAT=stat) altitude
IF (stat /= 0) THEN
    WRITE(*,*) "Error: Unable to parse altitude parameter"
    STOP
END IF

! Read observation type
CALL GET_COMMAND_ARGUMENT(7, arg)
READ(arg, *, IOSTAT=stat) obsType
IF (stat /= 0) THEN
    WRITE(*,*) "Error: Unable to parse observation type parameter"
    STOP
END IF

! Read minimum elevation
CALL GET_COMMAND_ARGUMENT(8, arg)
READ(arg, *, IOSTAT=stat) minEle
IF (stat /= 0) THEN
    WRITE(*,*) "Error: Unable to parse minimum elevation parameter"
    STOP
END IF

! Read solar distance
CALL GET_COMMAND_ARGUMENT(9, arg)
READ(arg, *, IOSTAT=stat) solarDistance
IF (stat /= 0) THEN
    WRITE(*,*) "Error: Unable to parse solar distance parameter"
    STOP
END IF

! Read lunar distance
CALL GET_COMMAND_ARGUMENT(10, arg)
READ(arg, *, IOSTAT=stat) lunarDistance
IF (stat /= 0) THEN
    WRITE(*,*) "Error: Unable to parse lunar distance parameter"
    STOP
END IF

! Read trajectory filename
CALL GET_COMMAND_ARGUMENT(11, fnEph)

! Read optional data output filename
IF (nArgs >= 12) THEN
    CALL GET_COMMAND_ARGUMENT(12, fnData)
    hasDataFile = .TRUE.
END IF

!---------------------- Parameter Display Section ----------------------
! Display only the most critical parameter information
WRITE(*,*) "Time range: ", ts, "to", te, "MJD, Step:", step, "seconds"
WRITE(*,*) "Observer location: Longitude=", longitude, "deg, Latitude=", latitude, "deg, Altitude=", altitude, "m"
WRITE(*,*) "Observation conditions: Type=", obsType, ", Min. Elevation=", minEle, "deg"
WRITE(*,*) "trajectory file:", TRIM(fnEph)
IF (hasDataFile) THEN
    WRITE(*,*) "Output file:", TRIM(fnData)
END IF

!---------------------- File Check Section ----------------------
! Check if trajectory file exists
INQUIRE(FILE=TRIM(fnEph), EXIST=hasDataFile)
IF (.NOT. hasDataFile) THEN
    WRITE(*,*) "Error: trajectory file does not exist -", TRIM(fnEph)
    STOP
END IF

!---------------------- Calculation Preparation Section ----------------------
WRITE(*,*) "Processing started..."

! Initialize time
call t0%init(ts)
nStep = int((te - ts)*86400d0/step)  ! Calculate total steps
call timeS%init(t0, step, nStep)     ! Initialize time series

!---------------------- Scenario Establishment Section ----------------------
! Set time series
ss%ts = timeS
ss%eleType = STATE_TYPE_CARTESIAN
ss%nSat = 1

! Initialize satellite
allocate(ss%sat(1))
call ss%sat(1)%init(t0, (/0d0, 0d0, 0d0, 0d0, 0d0, 0d0/), name='sat1', number=1)
call ss%sat(1)%readtrajectory(STATE_TYPE_CARTESIAN, trim(fnEph))

! Initialize observation station
ss%nStation = 1
allocate(ss%station(1))
call ss%station(1)%init(number=1, name='station1', longitude=longitude, latitude=latitude, altitude=altitude)
call ss%station(1)%generateTrajectory(timeS)

! Initialize access relationship
ss%nAccess = 1
allocate(ss%access(1))
call ss%access(1)%init(number=1, name='ac1', &
&   fromType=102, fromNum=1, toType=101, toNum=1, &
&   obsType=obsType, minEle = minEle, solarDistance = solarDistance, &
&   lunarDistance = lunarDistance)

!---------------------- Calculation Execution Section ----------------------
! Calculate observation data
call ss%observation

! Save data to file
call ss%access(1)%saveData(ss%ts, trim(fnData))

!---------------------- Program End Section ----------------------
WRITE(*,*) "Processing complete!"

END PROGRAM SatelliteObservation