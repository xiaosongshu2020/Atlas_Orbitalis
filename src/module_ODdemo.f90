!================================================================================
! MODULE: module_ODdemo
!
! Shared driver for the six single-observation-type orbit-determination demos
! in app/od_*.f90.
!
! Scenario (identical for all six programs):
!   - one MEO satellite: a = 27906.1 km, e = 0.001, i = 55 deg
!   - one Beijing station (116.4 E, 39.9 N, 50 m), minimum elevation 5 deg
!   - 2-day arc, 60 s sampling
!   - truth force model: 12x12 gravity + Sun/Moon + cannonball SRP
!   - estimated parameters: Cartesian position/velocity (6) + amrSRP (1),
!     i.e. type_varGroup label = 3, nVar = 7
!
! The program generates the observations from the truth trajectory, saves
! them, perturbs the initial state and amrSRP, then runs the batch least
! squares solver. Outputs per type under data/output/od_demo/:
!   <TYPE>_obs.txt        simulated observations
!   <TYPE>_eph.txt        estimated precise ephemeris
!   <TYPE>_residuals.txt  final O-C residuals
!   <TYPE>_log.txt        iteration log (iteration number, RMS, x)
!
! Author: Shawn Hu
!================================================================================
module module_ODdemo

    use module_const
    use module_time, only : type_t, timeAdd
    use module_elements, only : kpl2cts, cts_normal2SI
    use module_satellite, only : type_satellite
    use module_station, only : type_station
    use module_access, only : type_access, OBS_TYPE_RA_DEC, &
        &   OBS_TYPE_AZIMUTH_ELEVATION, OBS_TYPE_RANGING, OBS_TYPE_RANGERATE, &
        &   OBS_TYPE_POSITION, OBS_TYPE_POSITION_VELOCITY
    use module_satelliteScenario, only : type_satelliteScenario
    use module_residual, only : res
    use module_ODor, only : type_ODor, type_varGroup, OD_METHOD_BATCH, resRMSstr

    implicit none
    private
    public :: runODdemo

contains

    !==========================================================================
    ! Run one observation-type OD demo.
    !   obsType - OBS_TYPE_* constant
    !   tag     - short upper-case name used in file names and the log header
    !==========================================================================
    subroutine runODdemo(obsType, tag)
    implicit none
    integer, intent(in) :: obsType
    character(len=*), intent(in) :: tag

    !---------------- scenario settings -------------------------------------
    real(8), parameter :: step        = 60d0          ! sampling [s]
    real(8), parameter :: arcDay      = 2d0           ! arc length [day]
    real(8), parameter :: aKm         = 27906.1d0     ! semi-major axis [km]
    real(8), parameter :: ecc         = 0.001d0       ! eccentricity
    real(8), parameter :: incDeg      = 55d0          ! inclination [deg]
    real(8), parameter :: satMass     = 1000d0        ! [kg]
    real(8), parameter :: reflectivity = 0.5d0        ! Cr = 1 + reflectivity
    real(8), parameter :: srpTruth    = 0.02d0        ! true amrSRP [m^2/kg]
    real(8), parameter :: srpGuess    = 0.025d0       ! initial amrSRP [m^2/kg]
    integer, parameter :: gravityDegree = 12
    integer, parameter :: gravityOrder  = 12

    ! initial-state perturbation (truth + perturbation -> estimate start)
    real(8), parameter :: posPert(3) = (/ 50d0, -80d0, 120d0 /)        ! [m]
    real(8), parameter :: velPert(3) = (/ 0.005d0, -0.008d0, 0.006d0 /) ! [m/s]

    ! station and visibility
    character(len=*), parameter :: staName = "Beijing"
    real(8), parameter :: staLon = 116.4d0, staLat = 39.9d0, staAlt = 50d0
    real(8), parameter :: minEle = 5d0

    character(len=*), parameter :: outDir = "data/output/od_demo/"

    !---------------- local variables ---------------------------------------
    type(type_t) :: t0, te
    type(type_satelliteScenario) :: ss
    type(type_satellite) :: sat, truthSat
    type(type_station)   :: sta
    type(type_access)    :: ac
    type(type_ODor)      :: ODor
    type(type_varGroup)  :: Groups(1)
    type(type_t), allocatable :: epochs(:)

    real(8) :: kpl(6), cts(6), truth0(6), pvEst(6), pvTru(6), dr(3)
    real(8) :: sum2, maxErr, rmsErr
    integer :: i, k
    character(len=240) :: fname, obsFile, ephFile, resFile, logFile
    character(len=200) :: intro
    character(len=80)  :: rmsStr

    !---------------- output directory --------------------------------------
    call execute_command_line('mkdir "data/output/od_demo"', wait=.true.)

    !---------------- time and scenario -------------------------------------
    call t0%init(2025d0, 1d0, 1d0, 0d0, 0d0, 0d0)
    te = timeAdd(t0, arcDay)

    write(fname, "(A,A,A)") "OD demo (", trim(tag), "): Cartesian state + amrSRP"
    intro = trim(fname)
    call ss%init(t0, te, step, STATE_TYPE_CARTESIAN, introduction = trim(intro))

    !---------------- truth satellite ---------------------------------------
    kpl(1) = aKm * u_km
    kpl(2) = ecc
    kpl(3) = incDeg * deg
    kpl(4) = 0d0
    kpl(5) = 0d0
    kpl(6) = 0d0
    call kpl2cts(kpl, cts)

    call sat%init(t0, cts, number = 1, name = "MEO", mass = satMass, &
        &   eleType = STATE_TYPE_CARTESIAN, propagator = 'numerical', &
        &   degree = gravityDegree, order = gravityOrder, &
        &   reflectivity = reflectivity, amrSRP = srpTruth)
    sat%opMod%usePV          = 1   ! integrate Cartesian position/velocity
    sat%opMod%key_moonGravity = 1
    sat%opMod%key_sunGravity  = 1
    sat%opMod%key_srp         = 1
    sat%opMod%SRP_Model       = 1   ! cannonball: amrSRP is the estimated parameter
    call ss%addSatellite(sat)

    !---------------- ground station ----------------------------------------
    call sta%init(number = 1, name = staName, longitude = staLon, &
        &   latitude = staLat, altitude = staAlt)
    call ss%addStation(sta)

    !---------------- observation link --------------------------------------
    call ac%init(number = 1, name = trim(tag), &
        &   fromType = OBJECT_TYPE_SATELLITE, fromNum = 1, &
        &   toType = OBJECT_TYPE_STATION, toNum = 1, &
        &   obsType = obsType, minEle = minEle)
    ac%shadow = .false.       ! this demo maximizes data; no eclipse filter
    ac%solarDistance = 0d0
    ac%lunarDistance = 0d0
    call ss%addAccess(ac)

    !---------------- truth trajectory and observations ---------------------
    allocate(epochs(ss%ts%n))
    do i = 1, ss%ts%n
        epochs(i) = ss%ts%time(i)
    end do

    call ss%sat(1)%generateTrajectory(epochs)
    call ss%station(1)%generateTrajectory(ss%ts)
    truthSat = ss%sat(1)          ! deep copy: truth ephemeris kept for comparison
    truth0   = ss%sat(1)%states%elements

    call ss%observation(1)        ! generate observations of this type
    if (ss%access(1)%nData <= 0) then
        print*, "ERROR @runODdemo: no observations generated for ", trim(tag)
        stop 1
    endif

    obsFile = trim(outDir)//trim(tag)//"_obs.txt"
    call ss%access(1)%saveData(ss%ts, trim(obsFile))
    write(*, "(A,A,A,I0,A)") " ", trim(tag), ": epochs = ", ss%access(1)%nData, &
        &   "  obs -> "//trim(obsFile)

    !---------------- perturb the initial guess ------------------------------
    ss%sat(1)%states%elements(1:3) = ss%sat(1)%states%elements(1:3) &
        &   + posPert/(LENGTH_UNIT*1000d0)
    ss%sat(1)%states%elements(4:6) = ss%sat(1)%states%elements(4:6) &
        &   + velPert/(LENGTH_UNIT*1000d0)*TIME_UNIT
    ss%sat(1)%opMod%amrSRP = srpGuess

    !---------------- estimation setup: 6 states + amrSRP --------------------
    Groups(1)%nVar   = 7
    Groups(1)%label  = 3          ! states + amrSRP
    Groups(1)%origin = 1
    allocate(Groups(1)%idx(7))
    Groups(1)%idx = (/ 1, 2, 3, 4, 5, 6, 7 /)

    call ODor%init(ss, ss%ts%t0, 1, Groups, OD_METHOD_BATCH, intro = trim(intro))
    logFile = trim(outDir)//trim(tag)//"_log.txt"
    ODor%logFile = trim(logFile)
    ODor%showProgress = .false.

    if (obsType == OBS_TYPE_AZIMUTH_ELEVATION) then
        write(*, "(A)") " note: module_ODor converts Az/El to RA/Dec for the batch estimator"
    endif

    !---------------- solve --------------------------------------------------
    call ODor%solve(ss)
    call ODor%computeResiduals(ss)

    resFile = trim(outDir)//trim(tag)//"_residuals.txt"
    call ODor%saveResiduals(trim(resFile))
    ephFile = trim(outDir)//trim(tag)//"_eph.txt"
    call ODor%saveEph(ss, trim(ephFile))

    !---------------- report -------------------------------------------------
    call resRMSstr(res%nAcc, rmsStr)
    call cts_normal2SI(ss%sat(1)%states%elements, cts)
    call cts_normal2SI(truth0, pvTru)

    write(*, "(A)") ""
    write(*, "(A)") "================ OD result: "//trim(tag)//" ================"
    write(*, "(A,I0)") " iterations            : ", ODor%nIter
    write(*, "(A,A)")  " final residual RMS    : ", trim(rmsStr)
    write(*, "(A,3F14.3)") " estimated r0 (km)     : ", cts(1:3)
    write(*, "(A,3F14.3)") " truth     r0 (km)     : ", pvTru(1:3)
    write(*, "(A,3ES14.5)") " error dr (m)          : ", (cts(1:3) - pvTru(1:3))*1000d0
    write(*, "(A,3ES14.5)") " error dv (m/s)        : ", (cts(4:6) - pvTru(4:6))*1000d0
    write(*, "(A,2ES14.5)") " amrSRP est / truth    : ", ss%sat(1)%opMod%amrSRP, srpTruth

    ! ephemeris error versus the truth trajectory
    call ss%sat(1)%generateTrajectory(epochs)
    sum2 = 0d0
    maxErr = 0d0
    do k = 1, ss%ts%n
        pvEst = ss%sat(1)%getPV(epochs(k))
        pvTru = truthSat%getPV(epochs(k))
        dr = (pvEst(1:3) - pvTru(1:3))*LENGTH_UNIT*1000d0
        sum2 = sum2 + dot_product(dr, dr)
        maxErr = max(maxErr, norm2(dr))
    end do
    rmsErr = dsqrt(sum2/dble(ss%ts%n))
    write(*, "(A,ES12.3,A,ES12.3,A)") " ephemeris position error: RMS = ", rmsErr, &
        &   " m, max = ", maxErr, " m"
    write(*, "(A)") " saved: "//trim(ephFile)
    write(*, "(A)") " saved: "//trim(resFile)
    write(*, "(A)") " saved: "//trim(logFile)
    write(*, "(A)") "================================================="

    !---------------- cleanup ------------------------------------------------
    deallocate(Groups(1)%idx)
    deallocate(epochs)
    call ss%freeSS

    end subroutine runODdemo

end module module_ODdemo
