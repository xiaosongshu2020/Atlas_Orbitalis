!================================================================================
! MODULE: module_satellite
!
! Author: Shawn Hu
! Last Modified: 2025-12-21
!
! DESCRIPTION:
! This module defines the 'type_satellite' object for orbital mechanics simulations.
! It provides a comprehensive framework for satellite orbit propagation, trajectory
! management, and force modeling.
!
! SATELLITE OBJECT PROPERTIES:
!   number      - Unique identifier (integer)
!   name        - Satellite name (character*100)
!   mass        - Mass in kg (double precision)
!   states      - Initial state (type_trajectory: time + 6 orbital elements)
!   opMod       - Propagation settings (type_opMod)
!   eleType     - Element type (integer: 1=Keplerian, 2=Cartesian, 3=NSG1, 4=NSG2)
!   trajectory  - Generated ephemeris array (allocatable type_trajectory)
!   nEph        - Number of trajectory points (integer)
!   propagator  - Propagation method (character*10: 'numerical', 'two-body', 'SGP4')
!   thruster    - Thruster model (type_thruster)
!   clock       - Satellite clock model (type_clock: bias/drift + clock table)
!
! SATELLITE METHODS:
!   init()                     - Initialize satellite with parameters
!   addThruster(aThruster)     - Attach thruster model
!   generateTrajectory(Epochs) - Generate trajectory at specified times
!   savetrajectory()           - Save trajectory to file
!   readtrajectory()           - Load trajectory from file
!   setTrajectory() - Load trajectory from matrix
!   getStates(t)           - Get orbital elements at time t via interpolation
!   getPV(t)                   - Get Cartesian position/velocity at time t
!   show()                     - Display satellite information
!   propagate(epoch)           - Propagate satellite state to new epoch
!
! FORCE MODEL PARAMETERS (type_opMod):
!   Cd, amrDrag, F107, Ap, reflectivity, amrSRP - Environmental parameters
!   degree, order - Earth gravity field degree n and order m (0-70): the field
!                   is summed over n = 2..degree, m = 0..min(n, order)
!   key_moonGravity, key_moonTide, key_sunGravity, key_sunTide,
!   key_srp, key_drag, key_postNewton - Force model switches (0=off, 1=on)
!   usePV      - 0=use NSG1 elements, 1=use Cartesian elements
!   showProgressBar - Display progress during propagation
!
! USAGE:
!   use module_satellite
!
!   type(type_satellite) :: sat
!   type(type_t) :: t0
!   double precision :: x0(6), epochs(10)
!   type(type_t) :: time_array(10)
!
!   ! Example 1: Basic initialization and propagation
!   call t0%init(2025d0, 1d0, 1d0, 0d0, 0d0, 0d0)  ! Jan 1, 2025, 00:00:00 UTC
!   x0 = [7000.0d0 * km, 0.01d0, 0.1d0 * deg, 1.0d0 * deg, 0.5d0 * deg, 0.2d0 * deg]  ! Keplerian: a,e,i,Ω,ω,M
!   call sat%init(t0, x0, number=1, name="TestSat", mass=500.0d0, &
!                 eleType=1, propagator="numerical", degree=8, order=8)
!   call sat%show()
!
!   ! Example 2: Generate and save trajectory
!   do i = 1, 10
!     time_array(i) = timeAdd(t0, dble(i-1)/24d0)  ! Hourly for 9 hours
!   end do
!   call sat%generateTrajectory(time_array)
!   call sat%savetrajectory(STATE_TYPE_CARTESIAN, "trajectory.txt", freeEph=.true.)  ! Save as Cartesian
!   print*, 'final states:', sat%trajectory(sat%nEph)%elements
!
!   ! Example 3: Load trajectory and interpolate
!   call sat%readtrajectory(STATE_TYPE_CARTESIAN, "existing_traj.txt")
!   call t0%init(2025d0, 1d0, 1d0, 1d0, 15d0, 0d0)  ! 1h15m
!   state = sat%getStates(t0)        ! Get orbital elements
!   pv = sat%getPV(t0)                   ! Get position/velocity
! NOTES:
!   - All internal calculations use normalized units (Earth radii, TU)
!   - File I/O converts to SI units (km, km/s) or degrees as appropriate
!   - Time system: UTC (converted internally as needed)
!   - Numerical propagator uses 8th-order Runge-Kutta with configurable forces
!   - Two-body propagator uses Keplerian motion (fast, no perturbations)
!================================================================================
module module_satellite
    use module_time
    use module_const
    use module_orbitPredictionNumerical
    use module_thruster
    use module_clock
    use CMD_ProgressBar
    use module_toolKits
    use module_matrix
    implicit none

    !--------------------------------------------------------------------------
    ! TYPE: type_opMod
    ! Defines dynamics model parameters for orbit propagation
    !--------------------------------------------------------------------------
    type :: type_opMod
        ! about drag
        double precision  :: Cd = 2.2d0            ! Drag coefficient
        double precision  :: amrDrag = 0.02d0      ! Area-to-mass ratio for drag [m^2/kg]
        double precision  :: F107 = 130d0          ! Solar flux index
        double precision  :: Ap = 15d0             ! Geomagnetic index
        ! about SRP
        integer :: SRP_Model = 1 ! 1 - canon ball; 2 - ECOM (9 paras)
        double precision  :: reflectivity = 0d0    ! Surface reflectivity (0-1)
        double precision  :: amrSRP = 0.02d0       ! Area-to-mass ratio for SRP [m^2/kg]
        ! ECOM (SRP_Model = 2) 的 9 个参数, 归一化加速度单位。读取顺序由
        ! module_perturbationSRP%solar_pressure_perturbation 定义, 别改这里:
        !   1 D0  2 Y0  3 B0  4 Bc  5 Bs  6 Dc  7 Ds  8 Yc  9 Ys
        !   F = (D0 + Dc*cos u + Ds*sin u)*e_D
        !     + (Y0 + Yc*cos u + Ys*sin u)*e_Y
        !     + (B0 + Bc*cos u + Bs*sin u)*e_B
        double precision  :: ECOM(9) = 0d0
        
        ! 地球引力场：degree = 阶数(n)，order = 次数(m)，二者独立设置，
        ! 实际取到 n = 2..degree、m = 0..min(n, order)。例：C(5,6) 这类项需要
        ! degree = 5、order = 6（m 恒 <= n，给大了也只到 n）。
        ! degree = 0 或 1 时整个地球引力场都不算(仅二体)。order = 0 即只算带谐项。
        ! 缺省 degree = order = 2，即 2x2 场（含 J2 与 C21/S21、C22/S22）。
        integer :: degree = 2            ! 引力场阶数 n (0-70)
        integer :: order = 2             ! 引力场次数 m (0-70)
        ! 摄动力开关 (0 = off, 1 = on)，默认全部关闭，即无日月三体/光压/阻力/
        ! 后牛顿项。地球引力场不在这组开关里，只由上面的 degree/order 决定
        integer :: key_moonGravity = 0   ! 月球三体引力
        integer :: key_moonTide = 0      ! 月球固体潮
        integer :: key_sunGravity = 0    ! 太阳三体引力
        integer :: key_sunTide = 0       ! 太阳固体潮
        integer :: key_srp = 0           ! 太阳光压
        integer :: key_drag = 0          ! 大气阻力
        integer :: key_postNewton = 0    ! 后牛顿效应

        integer :: orderMoon = 0        ! Moon gravity order, two-body model in default
        integer :: Mercury = 0          ! Include Mercury gravity (0=No, 1=Yes)
        integer :: Venus = 0            ! Include Venus gravity (0=No, 1=Yes)
        integer :: Mars = 0             ! Include Mars gravity (0=No, 1=Yes)
        integer :: Jupiter = 0          ! Include Jupiter gravity (0=No, 1=Yes)
        integer :: Saturn = 0           ! Include Saturn gravity (0=No, 1=Yes)
        integer :: Uranus = 0           ! Include Uranus gravity (0=No, 1=Yes)
        integer :: Pluto = 0            ! Include Pluto gravity (0=No, 1=Yes)
        integer :: usePV = 0             ! type of states used in propagation; 0 = use NSG1; 1 = use Cartesian
        
        ! about ERP
        integer :: key_ERP = 0              ! Earth Radiation Pressure switch (0=off, 1=on)
        double precision :: reflectivityERP = 0.0d0  ! ERP surface reflectivity (0~1)
        double precision :: amrERP = 0.02d0 ! Area-to-mass ratio for ERP [m^2/kg]
        
        LOGICAL :: showProgressBar = .False. ! Display progress bar during propagation
        double precision  :: step = 0d0            ! step length(time_unit), 0 means not fixed (adjust automatically).
    contains
        ! 内部桥接: 把 7 个命名开关组装为 NumOP 所需的 K(1:7) 数组。
        ! 声明为 private —— 面向 satellite 对象时不暴露任何 K 接口, 仅供本模块调用
        procedure, private :: forceSwitchArray => opMod_forceSwitchArray
        ! 关闭全部 7 个摄动力开关。只置这 7 个开关, 不影响 degree/order、
        ! key_ERP、orderMoon、各行星开关与 usePV 等
        procedure :: resetSwitches => opMod_resetSwitches
    end type type_opMod

    !--------------------------------------------------------------------------
    ! TYPE: type_trajectory
    ! Single epoch of satellite trajectory (time + orbital elements)
    !--------------------------------------------------------------------------
    type :: type_trajectory
        type(type_t) :: t               ! UTC time
        double precision       :: elements(6) = 0d0  ! Orbital elements (normalized units)
    end type type_trajectory
    
    !--------------------------------------------------------------------------
    ! TYPE: type_satellite
    ! Main satellite type containing all properties and methods
    !--------------------------------------------------------------------------
    type :: type_satellite
        integer              :: number = 0          ! Satellite ID number
        character*100        :: name = ""           ! Satellite name
        double precision     :: mass = 500d0        ! Mass [kg]
        type(type_trajectory) :: states             ! Initial state at epoch
        type(type_opMod)     :: opMod               ! Propagation settings
        integer              :: eleType = 3         ! Element type (default NSG1)
        type(type_trajectory), allocatable :: trajectory(:)  ! Generated ephemeris
        integer              :: nEph = 0            ! Number of trajectory points
        character*10         :: propagator = "numerical"  ! Propagation method
        type(type_thruster)  :: thruster            ! Thruster model (if any)
        type(type_clock)     :: clock               ! Satellite clock model (bias/drift + clock table)
    contains
        procedure :: init => sat_init
        procedure :: addThruster
        procedure :: generateTrajectory
        procedure :: generateTrajectoryTwobody
        procedure :: savetrajectory
        procedure :: readtrajectory
        procedure :: setTrajectory
        procedure :: freeTrajectory
        procedure :: getStates => sat_getStates
        procedure :: getPV => sat_getPV
        procedure :: show => showSatelliteInfor
        procedure :: propagate => sat_propagate
    end type type_satellite

contains

    !--------------------------------------------------------------------------
    ! SUBROUTINE: sat_init
    ! Initialize satellite object with initial state and parameters
    !
    ! INPUTS:
    !   t          - Initial epoch (UTC)
    !   x0(6)      - Initial orbital elements (type depends on eleType)
    !   number     - (optional) Satellite ID
    !   name       - (optional) Satellite name
    !   mass       - (optional) Mass [kg]
    !   eleType    - (optional) Element type (see module_const for constants)
    !   propagator - (optional) 'numerical', 'two-body', 'SGP4'
    !   Cd..degree - (optional) Force model parameters (see type_opMod)
    !
    ! NOTES:
    !   - 地球引力场用 degree/order 两个关键字直接给(阶数/次数);
    !     其余力模型开关不在本接口上, 初始化后直接对 sat%opMod%key_* 赋值, 例如
    !       sat%opMod%key_moonGravity = 1
    !   - 未给出的参数保持对象原有取值, 不会被清零
    !
    ! OUTPUTS:
    !   this - Initialized satellite object
    !--------------------------------------------------------------------------
    subroutine sat_init(this, t, x0, number, name, mass, eleType, propagator, &
                        Cd, amrDrag, F107, Ap, reflectivity, amrSRP, degree, order, &
                        reflectivityERP, amrERP, key_ERP)
        class(type_satellite), intent(inout) :: this
        type(type_t), intent(in) :: t
        double precision, intent(in) :: x0(6)
        integer, optional, intent(in) :: number
        character(*), optional, intent(in) :: name
        double precision, optional, intent(in) :: mass
        integer, optional, intent(in) :: eleType
        character(*), optional, intent(in) :: propagator
        double precision, optional, intent(in) :: Cd, amrDrag, F107, Ap, reflectivity, amrSRP
        double precision, optional, intent(in) :: reflectivityERP, amrERP
        integer, optional, intent(in) :: degree, order, key_ERP
        double precision :: elementsNormal(6), nsg1(6), cts(6)
        
        ! Set optional parameters if provided
        if (present(number)) this%number = number
        if (present(name)) this%name = trim(name)
        if (present(mass)) this%mass = mass
        if (present(eleType)) this%eleType = eleType
        if (present(propagator)) this%propagator = propagator

        ! Set force model parameters
        if (present(Cd)) this%opMod%Cd = Cd
        if (present(amrDrag)) this%opMod%amrDrag = amrDrag
        if (present(F107)) this%opMod%F107 = F107
        if (present(Ap)) this%opMod%Ap = Ap
        if (present(reflectivity)) this%opMod%reflectivity = reflectivity
        if (present(amrSRP)) this%opMod%amrSRP = amrSRP
        if (present(reflectivityERP)) this%opMod%reflectivityERP = reflectivityERP
        if (present(amrERP)) this%opMod%amrERP = amrERP
        if (present(key_ERP)) this%opMod%key_ERP = key_ERP
        if (present(degree)) this%opMod%degree = degree
        if (present(order)) this%opMod%order = order

        ! Store initial state
        this%states%t = t
        this%states%elements = x0

        ! Initialize trajectory array
        this%nEph = 0
        if (allocated(this%trajectory)) deallocate(this%trajectory)
    end subroutine sat_init

    !--------------------------------------------------------------------------
    ! FUNCTION: opMod_forceSwitchArray
    ! 把 7 个命名摄动力开关组装成 NumOP 所需的 K(1:7) 数组。
    !
    ! 这是 satellite 侧与 NumOP 之间唯一的转换点, 索引顺序只在此处定义:
    !   1 月球三体  2 月球固体潮  3 太阳三体
    !   4 太阳固体潮  5 太阳光压  6 大气阻力  7 后牛顿
    ! (地球引力场不在其中, 它由 degree/order 单独传给 NumOP)
    ! 该绑定声明为 private, 面向 satellite 对象时不暴露, 仅供本模块调用。
    !--------------------------------------------------------------------------
    function opMod_forceSwitchArray(this) result(K)
        class(type_opMod), intent(in) :: this
        integer :: K(7)

        K(1) = this%key_moonGravity
        K(2) = this%key_moonTide
        K(3) = this%key_sunGravity
        K(4) = this%key_sunTide
        K(5) = this%key_srp
        K(6) = this%key_drag
        K(7) = this%key_postNewton
    end function opMod_forceSwitchArray

    !--------------------------------------------------------------------------
    ! SUBROUTINE: opMod_resetSwitches
    ! 关闭全部 7 个摄动力开关(整组替换前的清零用)。
    !
    ! 只置这 7 个开关, 不触碰 degree/order、key_ERP、orderMoon、各行星开关与
    ! usePV —— 这些不属于摄动力开关组, 需由调用者自行决定。
    !--------------------------------------------------------------------------
    subroutine opMod_resetSwitches(this)
        class(type_opMod), intent(inout) :: this

        this%key_moonGravity = 0
        this%key_moonTide    = 0
        this%key_sunGravity  = 0
        this%key_sunTide     = 0
        this%key_srp         = 0
        this%key_drag        = 0
        this%key_postNewton  = 0
    end subroutine opMod_resetSwitches

    !--------------------------------------------------------------------------
    ! SUBROUTINE: addThruster
    ! Attach a thruster model to the satellite
    !
    ! INPUTS:
    !   aThruster - Thruster object containing force/torque model
    !--------------------------------------------------------------------------
    subroutine addThruster(this, aThruster)
        class(type_satellite), intent(inout) :: this
        type(type_thruster), intent(in) :: aThruster
        this%thruster = aThruster
    end subroutine addThruster

    !--------------------------------------------------------------------------
    ! SUBROUTINE: generateTrajectory
    ! Generate trajectory at specified epochs using selected propagator
    !
    ! INPUTS:
    !   Epochs(:) - Array of UTC times for trajectory output
    !--------------------------------------------------------------------------
    subroutine generateTrajectory(this, Epochs)
        class(type_satellite) :: this
        type(type_t), intent(in) :: Epochs(:)
        
        select case(trim(this%propagator))
        case("numerical") 
            call generateTrajectoryNumberical(this, Epochs)
        case("two-body")
            call generateTrajectoryTwobody(this, Epochs)
        case("SGP4")
            call generateTrajectorySGP4(this, Epochs)
        case default
            print*, "Error: unknown propagator '"//trim(this%propagator)//"'"
            stop
        end select
    end subroutine

    !--------------------------------------------------------------------------
    ! SUBROUTINE: generateTrajectoryNumberical
    ! Generate trajectory using numerical integration with full force models
    !
    ! INPUTS:
    !   Epochs(:) - Array of UTC output times
    !--------------------------------------------------------------------------
    subroutine generateTrajectoryNumberical(this, Epochs)
        class(type_satellite) :: this
        type(type_t), intent(in) :: Epochs(:)
        type(type_t) :: t0, t1
        integer :: nEpochs, i, nStep
        double precision :: x0(6), x1(6), utc0, dt
        logical :: isNewSat = .True.
        integer, save :: lastSat = -1
        double precision :: kpl(6), nsg1(6), cts(6), nsg2(6)
        integer :: threadNum
        integer :: eleType
        
        ! Check if this is a new satellite (for thread-local storage)
        if (this%number == lastSat) then
            isNewSat = .false.
        else
            isNewSat = .true.
            lastSat = this%number
        endif
        
        ! Publish the satellite properties that the force model reads. `mass` and
        ! `thisThruster` are the module variables of module_orbitPredictionNumerical
        ! (visible here through `use module_orbitPredictionNumerical`); they must NOT
        ! be declared locally in this routine, otherwise the assignments below only
        ! fill a shadow copy and the force model silently keeps its own defaults
        ! (is_active = .false. and mass = 500 kg).
        ! Note this makes the propagation state global: generating trajectories for
        ! several satellites concurrently is not safe.
        thisThruster = this%thruster
        mass = this%mass

        ! >>> Step 1: Get initial orbit elements in appropriate form
        nEpochs = size(Epochs)
        t0 = this%states%t
        
        if (this%opMod%usePV == 0) then
            ! Use non-singular 1 elements for propagation
            if (this%eleType == STATE_TYPE_NONSINGULAR_1) then 
                x0 = this%states%elements
            elseif (this%eleType == STATE_TYPE_KEPLERIAN) then
                kpl = this%states%elements
                call kpl2nsg1(kpl, x0)
            elseif (this%eleType == STATE_TYPE_CARTESIAN) then
                cts = this%states%elements
                call cts2nsg1(cts, x0)
            elseif (this%eleType == STATE_TYPE_NONSINGULAR_2) then
                nsg2 = this%states%elements
                call nsg22nsg1(nsg2, x0)
            else
                print*, "Error: unsupported element type @generateTrajectory."
                STOP
            endif
        else
            ! Use Cartesian elements for propagation
            if (this%eleType == STATE_TYPE_NONSINGULAR_1) then 
                nsg1 = this%states%elements
                call nsg12kpl(nsg1, kpl)
                call kpl2cts(kpl, x0)
            elseif (this%eleType == STATE_TYPE_KEPLERIAN) then
                kpl = this%states%elements
                call kpl2cts(kpl, x0)
            elseif (this%eleType == STATE_TYPE_CARTESIAN) then
                x0 = this%states%elements
            elseif (this%eleType == STATE_TYPE_NONSINGULAR_2) then
                nsg2 = this%states%elements
                call nsg22kpl(nsg2, kpl)
                call kpl2cts(kpl, x0)
            else
                print*, "Error: unsupported element type @generateTrajectory."
                STOP
            endif
        endif
        
        ! Allocate memory for trajectory storage
        this%nEph = nEpochs
        if (.not.allocated(this%trajectory)) then
            allocate(this%trajectory(nEpochs))
        else
            deallocate(this%trajectory)
            allocate(this%trajectory(nEpochs))
        endif
        
        !threadNum = omp_get_thread_num()  ! For OpenMP parallelization
        
        ! >>> Step 2: Propagate orbit to each output epoch
        progress%prefix = "prediction for "//this%name//':'

        ! 进度条刷新间隔: 逐历元刷新对 40 格的条没有任何意义, 却要多做 nEpochs 次
        ! 终端写 —— 输出重定向到文件时更会刷屏(每个历元留下半行)。每约 1% 刷一次,
        ! 末历元必刷(把条补满、补上换行)。
        nStep = max(1, nEpochs/100)
        do i = 1, nEpochs
            ! Show progress bar if requested
            if (this%opMod%showProgressBar) then
                if (mod(i, nStep) == 0 .or. i == nEpochs) &
                    &   call progress%put(nEpochs, i, .true.)
            endif
            
            ! Calculate time step and propagate
            utc0 = t0%day + t0%sec/86400d0
            dt = timeMinus(epochs(i), t0)
            
            if (this%opMod%usePV == 0) then
                eleType = STATE_TYPE_NONSINGULAR_1
            else
                eleType = STATE_TYPE_CARTESIAN
            endif
            
            call NumOP(this%opMod%forceSwitchArray(), this%opMod%degree, this%opMod%order, &
                    &   this%opMod%Cd, this%opMod%amrDrag, this%opMod%F107, this%opMod%Ap, &
                    &   this%opMod%reflectivity, this%opMod%amrSRP, &
                    &   this%opMod%SRP_model, this%opMod%ECOM, &
                    &   this%opMod%step, &
                    &   this%opMod%orderMoon, &
                    &   this%opMod%Mercury, &
                    &   this%opMod%Venus, &
                    &   this%opMod%Mars, &
                    &   this%opMod%Jupiter, &
                    &   this%opMod%Saturn, &
                    &   this%opMod%Uranus, &
                    &   this%opMod%Pluto, &
                    &   this%opMod%key_ERP, &
                    &   this%opMod%reflectivityERP, &
                    &   this%opMod%amrERP, &
                    &   eleType, &
                    &   this%number, &
                    &   x0, utc0, dt, x1)

            ! Settle the thruster ledger over the step that has just been
            ! propagated. The force model divides the thrust by
            ! (mass - total_fuel_used), and the ledger is caller-driven, so
            ! without this flush the burnt mass would stay frozen at its
            ! published value for the whole run and the engine would never run
            ! dry. `thisThruster` is the module copy the force model reads.
            if (thisThruster%is_active) call thisThruster%update_state(epochs(i))

            ! Store propagated state
            this%trajectory(i)%t = epochs(i)
            
            if (this%opMod%usePV == 0) then
                select case (this%eleType) 
                case (STATE_TYPE_KEPLERIAN)
                    call nsg12kpl(x1, kpl)
                    this%trajectory(i)%elements = kpl
                case (STATE_TYPE_CARTESIAN)
                    call nsg12cts(x1, cts)
                    this%trajectory(i)%elements = cts
                case (STATE_TYPE_NONSINGULAR_1)
                    this%trajectory(i)%elements = x1
                case (STATE_TYPE_NONSINGULAR_2)
                    call nsg12nsg2(x1, nsg2)
                    this%trajectory(i)%elements = nsg2
                end select
            else
                select case (this%eleType) 
                case (STATE_TYPE_KEPLERIAN)
                    call cts2kpl(x1, kpl)
                    this%trajectory(i)%elements = kpl
                case (STATE_TYPE_CARTESIAN)
                    this%trajectory(i)%elements = x1
                case (STATE_TYPE_NONSINGULAR_1)
                    call cts2kpl(x1, kpl)
                    call kpl2nsg1(kpl, nsg1)
                    this%trajectory(i)%elements = nsg1
                case (STATE_TYPE_NONSINGULAR_2)
                    call cts2kpl(x1, kpl)
                    call kpl2nsg2(kpl, nsg2)
                    this%trajectory(i)%elements = nsg2
                end select
            endif
            
            ! Prepare for next step (sequential propagation)
            t0 = epochs(i)
            x0 = x1
            isNewSat = .false.
        enddo
        
        isNewSat = .false.
    end subroutine

    !--------------------------------------------------------------------------
    ! SUBROUTINE: generateTrajectoryTwobody
    ! Generate trajectory using two-body Keplerian propagation
    !
    ! INPUTS:
    !   Epochs(:) - Array of UTC output times
    !--------------------------------------------------------------------------
    subroutine generateTrajectoryTwobody(this, epochs)
        class(type_satellite) :: this
        type(type_t), intent(in) :: Epochs(:)
        double precision :: kpl0(6), kpl1(6), cts0(6), cts1(6), dt
        integer :: nEpochs, i
        type(type_t) :: t0
        double precision :: n
        
        nEpochs = size(Epochs)
        t0 = this%states%t
        
        ! Allocate trajectory storage
        this%nEph = nEpochs
        if (.not.allocated(this%trajectory)) then
            allocate(this%trajectory(nEpochs))
        else
            deallocate(this%trajectory)
            allocate(this%trajectory(nEpochs))
        endif
        
        ! Propagate using two-body Keplerian motion
        select case(this%eleType)
        case(STATE_TYPE_CARTESIAN)
            cts0 = this%states%elements
            call cts2kpl(cts0, kpl0)
            do i = 1, nEpochs
                dt = timeMinus(Epochs(i), t0)*86400d0/time_unit
                kpl1 = kpl0
                n = kpl0(1)**(-1.5d0)  ! Mean motion
                kpl1(6) = kpl0(6) + n*dt  ! Update mean anomaly
                kpl1(6) = dmod(kpl1(6), 2d0*pi)  ! Wrap to [0, 2π)
                call kpl2cts(kpl1, cts1)
                this%trajectory(i)%elements = cts1
                this%trajectory(i)%t = epochs(i)
            enddo
        case default 
            kpl0 = this%states%elements
            do i = 1, nEpochs
                dt = timeMinus(Epochs(i), t0)*86400d0/time_unit
                kpl1 = kpl0
                n = kpl0(1)**(-1.5d0)
                kpl1(6) = kpl0(6) + n*dt
                kpl1(6) = dmod(kpl1(6), 2d0*pi)
                this%trajectory(i)%elements = kpl1
                this%trajectory(i)%t = epochs(i)
            enddo
        end select
    end subroutine

    !--------------------------------------------------------------------------
    ! SUBROUTINE: generateTrajectorySGP4
    ! Placeholder for SGP4 propagation (not implemented)
    !
    ! INPUTS:
    !   Epochs(:) - Array of UTC output times
    !--------------------------------------------------------------------------
    subroutine generateTrajectorySGP4(this, epochs)
        class(type_satellite) :: this
        type(type_t), intent(in) :: Epochs(:)
        print*, "SGP4 propagator not implemented @generateTrajectorySGP4"
    end subroutine

    !--------------------------------------------------------------------------
    ! SUBROUTINE: savetrajectory
    ! Save satellite trajectory to text file with comprehensive header
    !
    ! INPUTS:
    !   statesType - Element type for output file
    !   fileName   - Output filename
    !   freeEph    - (optional) Deallocate trajectory after saving
    !
    ! FILE FORMAT:
    !   A metadata header whose every line starts with '#' (so a reader skips
    !   it without knowing its length), followed by
    !   Trajectory data (MJD_day, MJD_sec, state vector).
    !   The header names the frame of the state vector (J2000, ECI) and the
    !   time system of the epochs (UTC); the angles of the Keplerian/NSG output
    !   forms (i, RAAN, argp, M) are referred to the same J2000 frame.
    !--------------------------------------------------------------------------
    subroutine savetrajectory(this, statesType, fileName, freeEph)
        class(type_satellite), intent(inOut) :: this
        character(len = *), intent(in) :: fileName
        integer, intent(in) :: statesType
        logical, intent(in), optional :: freeEph
        double precision  :: kpl(6), cts(6), nsg1(6)
        type(type_t) :: t
        integer :: i
        double precision :: day, sec
        integer :: number
        integer :: eleType
        integer :: file_unit
        character(len=100) :: state_type_str, ele_type_str, prop_str
        double precision :: start_day, start_sec, end_day, end_sec
        character(len=50) :: start_time_str, end_time_str
        
        ! Check if trajectory exists
        if (this%nEph == 0) then
            print*, "WARNING: failed to save trajectory - no ephemeris points."
            return
        endif
        if (.not. allocated(this%trajectory)) then
            print*, "WARNING: failed to save trajectory - not allocated."
            return
        endif
        
        number = this%number
        eleType = this%eleType
        
        ! Get time range
        start_day = this%trajectory(1)%t%day
        start_sec = this%trajectory(1)%t%sec
        end_day = this%trajectory(this%nEph)%t%day
        end_sec = this%trajectory(this%nEph)%t%sec
        
        ! Convert element types to strings
        select case(statesType)
        case(STATE_TYPE_CARTESIAN)
            state_type_str = "Cartesian (position/velocity)"
        case(STATE_TYPE_KEPLERIAN)
            state_type_str = "Keplerian (a,e,i,Ω,ω,M)"
        case(STATE_TYPE_NONSINGULAR_1)
            state_type_str = "Non-singular 1"
        case(STATE_TYPE_NONSINGULAR_2)
            state_type_str = "Non-singular 2"
        case default
            state_type_str = "Unknown"
        end select
        
        select case(eleType)
        case(STATE_TYPE_CARTESIAN)
            ele_type_str = "Cartesian"
        case(STATE_TYPE_KEPLERIAN)
            ele_type_str = "Keplerian"
        case(STATE_TYPE_NONSINGULAR_1)
            ele_type_str = "Non-singular 1"
        case(STATE_TYPE_NONSINGULAR_2)
            ele_type_str = "Non-singular 2"
        case default
            ele_type_str = "Unknown"
        end select
        
        prop_str = trim(this%propagator)
        
        ! Open file with new unit to avoid hardcoded file numbers
        open(newunit=file_unit, file=trim(fileName), status='replace')
        
        ! Write the header. Every line starts with '#', so readers can skip it
        ! without knowing how many lines it holds (see readtrajectory).
        write(file_unit, '(A)') "# ============================================================================"
        write(file_unit, '(A)') "# SATELLITE EPHEMERIS FILE"
        write(file_unit, '(A)') "# ============================================================================"
        write(file_unit, '(A, I6)')   "# Satellite ID:", this%number
        write(file_unit, '(A, A)')    "# Name:        ", trim(this%name)
        write(file_unit, '(A, F12.3)') "# Mass (kg):   ", this%mass
        write(file_unit, '(A, A)')    "# Propagator:  ", trim(prop_str)
        write(file_unit, '(A, A)')    "# Input type:  ", trim(ele_type_str)
        write(file_unit, '(A, A)')    "# Output type: ", trim(state_type_str)
        write(file_unit, '(A)')       "# Frame:       J2000 (ECI, inertial; Keplerian/NSG angles refer to it)"
        write(file_unit, '(A)')       "# Time system: UTC"
        write(file_unit, '(A)') "# ============================================================================"
        write(file_unit, '(A)') "# DATA FORMAT (Fortran style)"
        write(file_unit, '(A)') "# Format: (F8.1, 3x, F18.10, 3(3x, F20.6), 3(3x, F18.7))"
        write(file_unit, '(A)') "# Columns:"
        write(file_unit, '(A)') "#   1: MJD_day     - Modified Julian Date (day part)"
        write(file_unit, '(A)') "#   2: MJD_sec     - Modified Julian Date (seconds part)"
        write(file_unit, '(A)') "#   3: X (km)      - Position X component"
        write(file_unit, '(A)') "#   4: Y (km)      - Position Y component"
        write(file_unit, '(A)') "#   5: Z (km)      - Position Z component"
        write(file_unit, '(A)') "#   6: VX (km/s)   - Velocity X component"
        write(file_unit, '(A)') "#   7: VY (km/s)   - Velocity Y component"
        write(file_unit, '(A)') "#   8: VZ (km/s)   - Velocity Z component"
        write(file_unit, '(A)') "# ============================================================================"
        write(file_unit, '(A, I8)')    "# Number of points: ", this%nEph
        write(file_unit, '(A, F12.1, A, F12.6)') "# Time range (MJD): ", start_day, " + ", start_sec/86400d0
        write(file_unit, '(A, F12.1, A, F12.6)') "#                   ", end_day, " + ", end_sec/86400d0
        write(file_unit, '(A)') "# Units: km, km/s (Cartesian) or km, deg (Keplerian/NSG); frame J2000"
        write(file_unit, '(A)') "# ============================================================================"
        
        ! Write trajectory data in requested format
        if (eleType == STATE_TYPE_NONSINGULAR_1) then
            if (statesType == STATE_TYPE_CARTESIAN) then
                do i = 1, this%nEph
                    day = this%trajectory(i)%t%day
                    sec = this%trajectory(i)%t%sec
                    nsg1 = this%trajectory(i)%elements
                    call nsg12kpl(nsg1, kpl)
                    call kpl2cts(kpl, cts)
                    call cts_normal2SI(cts, cts)
                    write(file_unit, "(F8.1,3x,F18.10,3(3x,F20.6),3(3x,F18.7))") day, sec, cts
                enddo
            elseif (statesType == STATE_TYPE_NONSINGULAR_1) then
                do i = 1, this%nEph
                    day = this%trajectory(i)%t%day
                    sec = this%trajectory(i)%t%sec
                    nsg1 = this%trajectory(i)%elements
                    nsg1(1) = nsg1(1)*length_unit
                    nsg1(2:3) = nsg1(2:3)*rad2deg
                    nsg1(6) = nsg1(6)*rad2deg
                    write(file_unit, "(F8.1,3x,F18.10,3x,F20.7,5(3x,F20.10))") day, sec, nsg1
                enddo
            elseif (statesType == STATE_TYPE_KEPLERIAN) then
                do i = 1, this%nEph
                    day = this%trajectory(i)%t%day
                    sec = this%trajectory(i)%t%sec
                    nsg1 = this%trajectory(i)%elements
                    call nsg12kpl(nsg1, kpl)
                    call kpl_normal2SI(kpl, kpl)
                    write(file_unit, "(F8.1,3x,F18.10,3x,F18.3,3x,F13.10,4(3x,F18.12))") day, sec, kpl
                enddo
            else
                print*, "Unsupported statesType for NSG1 input."
            endif
        elseif (eleType == STATE_TYPE_KEPLERIAN) then
            if (statesType == STATE_TYPE_CARTESIAN) then
                do i = 1, this%nEph
                    day = this%trajectory(i)%t%day
                    sec = this%trajectory(i)%t%sec
                    kpl = this%trajectory(i)%elements
                    call kpl2cts(kpl, cts)
                    call cts_normal2SI(cts, cts)
                    write(file_unit, "(F8.1,3x,F18.10,3(3x,F20.6),3(3x,F18.7))") day, sec, cts
                enddo
            elseif (statesType == STATE_TYPE_KEPLERIAN) then
                do i = 1, this%nEph
                    day = this%trajectory(i)%t%day
                    sec = this%trajectory(i)%t%sec
                    kpl = this%trajectory(i)%elements
                    call kpl_normal2SI(kpl, kpl)
                    write(file_unit, "(F8.1,3x,F18.10,3x,F18.3,3x,F13.10,4(3x,F18.12))") day, sec, kpl
                enddo
            else
                print*, "Unsupported statesType for Keplerian input."
            endif
        elseif (eleType == STATE_TYPE_CARTESIAN) then
            if (statesType == STATE_TYPE_CARTESIAN) then
                do i = 1, this%nEph
                    day = this%trajectory(i)%t%day
                    sec = this%trajectory(i)%t%sec
                    cts = this%trajectory(i)%elements
                    call cts_normal2SI(cts, cts)
                    write(file_unit, "(F8.1,3x,F18.10,3(3x,F20.6),3(3x,F18.7))") day, sec, cts
                enddo
            elseif (statesType == STATE_TYPE_NONSINGULAR_1) then
                do i = 1, this%nEph
                    day = this%trajectory(i)%t%day
                    sec = this%trajectory(i)%t%sec
                    cts = this%trajectory(i)%elements
                    call cts2kpl(cts, kpl)
                    call kpl2nsg1(kpl, nsg1)
                    nsg1(1) = nsg1(1)*length_unit
                    nsg1(2:3) = nsg1(2:3)*rad2deg
                    nsg1(6) = nsg1(6)*rad2deg
                    write(file_unit, "(F8.1,3x,F18.10,3x,F20.7,5(3x,F20.10))") day, sec, nsg1
                enddo
            elseif (statesType == STATE_TYPE_KEPLERIAN) then
                do i = 1, this%nEph
                    day = this%trajectory(i)%t%day
                    sec = this%trajectory(i)%t%sec
                    cts = this%trajectory(i)%elements
                    call cts2kpl(cts, kpl)
                    call kpl_normal2SI(kpl, kpl)
                    write(file_unit, "(F8.1,3x,F18.10,3x,F18.3,3x,F13.10,4(3x,F18.12))") day, sec, kpl
                enddo
            else
                print*, "Unsupported statesType for Cartesian input."
            endif
        else
            print*, "Unsupported element type in satellite object."
        endif
        
        close(file_unit)
        
        ! Free memory if requested
        if(present(freeEph)) then
            if (freeEph) then
                deallocate(this%trajectory)
                this%nEph = 0
            endif
        else
            deallocate(this%trajectory)
            this%nEph = 0
        endif
    end subroutine

    !--------------------------------------------------------------------------
    ! SUBROUTINE: setTrajectory
    ! Load trajectory from a pre-loaded data matrix
    !
    ! INPUTS:
    !   eleType - Element type in the matrix
    !   nData   - Number of data points
    !   mat     - Matrix with columns: day, sec, element1..element6
    !--------------------------------------------------------------------------
    subroutine setTrajectory(this, eleType, nData, mat)
        class(type_satellite) :: this
        integer, intent(in) :: eleType, nData
        double precision, intent(in) :: mat(nData, 8)
        integer :: i
        type(type_t) :: ti
        
        ! Reallocate trajectory memory
        if (allocated(this%trajectory)) deallocate(this%trajectory)
        allocate(this%trajectory(nData))
        
        this%eleType = eleType
        this%nEph = nData
        
        ! Load each epoch
        do i = 1, nData
            ti%day = mat(i, 1)
            ti%sec = mat(i, 2)
            ti%val = ti%day + ti%sec/86400d0
            this%trajectory(i)%t = ti
            this%trajectory(i)%elements = mat(i, 3:8)
        enddo
    end subroutine

    !--------------------------------------------------------------------------
    ! SUBROUTINE: readtrajectory
    ! Load trajectory from text file with header support
    !
    ! INPUTS:
    !   statesType - Element type in file
    !   fn         - Filename
    !
    ! FILE FORMAT:
    !   Header: any number of lines starting with a comment marker
    !           ('#', '*', '!', '！'); blank lines are skipped as well
    !   Data:   one epoch per line. Two layouts are recognised from the first
    !           data line:
    !             8 columns: MJD_day, MJD_sec, x, y, z, vx, vy, vz
    !             7 columns: MJD,       x, y, z, vx, vy, vz
    !
    ! NOTES:
    !   - The header length is NOT assumed: comment lines are skipped wherever
    !     they appear in the file (see isCommentLine in module_toolKits)
    !   - Currently assumes Cartesian format for data
    !--------------------------------------------------------------------------
    subroutine readtrajectory(this, statesType, fn)
        class(type_satellite) :: this
        integer, intent(in) :: statesType
        character(len=*), intent(in) :: fn
        integer :: nEph, nCol, i, nRead
        double precision, allocatable :: mat(:, :)
        double precision :: mjd
        character(len=512) :: first_line
        logical :: found


        ! Number of data lines: blank and comment lines are skipped, so the
        ! header may have any number of lines
        call getDataLinesNum(trim(fn), nEph)

        if (nEph <= 0) then
            print*, "ERROR: No trajectory data found in file:", trim(fn)
            return
        endif

        ! Number of columns, from the first data line: 8 = MJD_day + MJD_sec +
        ! state, 7 = MJD + state
        call getFirstDataLine(trim(fn), first_line, found)
        nCol = countFields(first_line)

        if (nCol /= 7 .and. nCol /= 8) then
            print*, "ERROR @readtrajectory: unsupported number of columns:", nCol
            print*, trim(fn)
            return
        endif

        print*, "Reading trajectory file,", nEph, "data points,", nCol, "columns"

        ! Allocate matrix and read data; malformed lines are reported by
        ! readDataMatrix and dropped, so only nRead rows carry data
        allocate(mat(nEph, 8))

        if (nCol == 8) then
            call readDataMatrix(trim(fn), nEph, 8, mat, nRead)
            nEph = nRead
            if (nEph == 0) then
                print*, "ERROR @readtrajectory: no usable data line in", trim(fn)
                return
            endif
        else
            ! MJD + state: expand to MJD_day, MJD_sec, state
            call readDataMatrix(trim(fn), nEph, 7, mat(:, 1:7), nRead)
            nEph = nRead
            if (nEph == 0) then
                print*, "ERROR @readtrajectory: no usable data line in", trim(fn)
                return
            endif
            mat(:, 8) = mat(:, 7)
            mat(:, 7) = mat(:, 6)
            mat(:, 6) = mat(:, 5)
            mat(:, 5) = mat(:, 4)
            mat(:, 4) = mat(:, 3)
            mat(:, 3) = mat(:, 2)
            ! split the MJD into its day and second parts: the MJD is left in
            ! column 1 (it is not moved by the shift above), the shifted state
            ! occupies columns 3:8, so column 2 is free for the seconds
            do i = 1, nEph
                mjd = mat(i, 1)
                mat(i, 1) = floor(mjd)
                mat(i, 2) = (mjd - mat(i, 1))*86400d0
            end do
        endif
        
        ! Convert from SI to normalized units
        ! Note: This assumes Cartesian format (km, km/s)
        ! For other formats, additional conversion may be needed
        mat(1:nEph, 3:5) = mat(1:nEph, 3:5)/length_unit
        mat(1:nEph, 6:8) = mat(1:nEph, 6:8)/length_unit*time_unit
        
        ! Store trajectory in satellite object
        ! Note: Currently assumes Cartesian format; for other formats,
        ! additional parameter may be needed to specify input type
        call this%setTrajectory(statesType, nEph, mat)
        
        deallocate(mat)
        
        ! Header lines (any number, marked with '#', '*', '!' or '！') are
        ! skipped during reading; the trajectory itself is loaded above.
        print*, "Trajectory loaded successfully"
    end subroutine

    !--------------------------------------------------------------------------
    ! FUNCTION: sat_getStates
    ! Get orbital elements at specified time via interpolation
    !
    ! INPUTS:
    !   t - Requested time (must be within trajectory bounds)
    !
    ! RETURNS:
    !   yt(6) - Interpolated orbital elements
    !--------------------------------------------------------------------------
    function sat_getStates(this, t) result(yt)
        class(type_satellite) :: this
        type(type_t), intent(in) :: t
        double precision :: yt(6)
        double precision :: ntDble, step
        integer :: nEph, n
        type(type_t) :: t0
        integer :: nt, left1, left2, right1, right2, j
        double precision :: x(4), y(4)
        
        ! Check trajectory availability
        if (.not.allocated(this%trajectory)) then
            print*, "Error @getStates: trajectory not allocated."
            print*, "Execute orbit prediction first."
            stop
        endif
        if (timeMinus(t, this%trajectory(1)%t) < 0d0 .or. &
            &   timeMinus(t, this%trajectory(this%nEph)%t) > 0d0) then
            print*, "Error @getStates: time outside trajectory bounds."
            print*, "Requested t = ", t%val
            print*, "Trajectory range: ", this%trajectory(1)%t%val, " to ", &
                   this%trajectory(this%nEph)%t%val
            stop
        endif
            
        t0 = this%trajectory(1)%t
        step = timeMinus(this%trajectory(2)%t, this%trajectory(1)%t)
        nEph = this%nEph
        
        ntDble = timeMinus(t, t0)/step

        ! Exact match to stored epoch
        if (dabs(ntDble - dble(nint(ntDble))) <= 1d-12) then
            n = nint(ntDble) + 1
            yt = this%trajectory(n)%elements
        else
            ! Interpolation using 4-point Lagrange
            n = floor(ntDble) + 1
            if (n == 1 .or. n == 2) then
                left1 = 1; left2 = 2; right1 = 3; right2 = 4
            elseif (n == nEph .or. n == nEph - 1) then
                left1 = nEph - 3; left2 = nEph - 2
                right1 = nEph - 1; right2 = nEph
            else
                left1 = n - 1; left2 = n; right1 = n + 1; right2 = n + 2
            endif
            
            ! Interpolate each element
            x(1:4) = this%trajectory([left1, left2, right1, right2])%t%val
            do j = 1, 6
                y(1:4) = this%trajectory([left1, left2, right1, right2])%elements(j)
                call Interpolation(4, x(1:4), y(1:4), 1, (/t%val/), yt(j))
            enddo
        endif
    end function

    !--------------------------------------------------------------------------
    ! FUNCTION: sat_getPV
    ! Get Cartesian position and velocity at specified time
    !
    ! INPUTS:
    !   t - Requested time
    !
    ! RETURNS:
    !   PV(6) - Cartesian state [position; velocity] in normalized units
    !--------------------------------------------------------------------------
    function sat_getPV(this, t) result(PV)
        class(type_satellite) :: this
        type(type_t), intent(in) :: t
        double precision :: PV(6)
        double precision :: yt(6)
        
        yt = this%getStates(t)
        
        ! Convert to Cartesian
        select case (this%eleType) 
        case(STATE_TYPE_CARTESIAN)
            PV = yt
        case(STATE_TYPE_KEPLERIAN)
            call kpl2cts(yt, PV)
        case(STATE_TYPE_NONSINGULAR_1)
            call nsg12cts(yt, PV)
        case(STATE_TYPE_NONSINGULAR_2)
            call nsg22cts(yt, PV)
        end select
    end function

    !--------------------------------------------------------------------------
    ! SUBROUTINE: showSatelliteInfor
    ! Display satellite information to screen
    !--------------------------------------------------------------------------
    subroutine showSatelliteInfor(this)
        class(type_satellite) :: this
        double precision :: kpl(6), nsg1(6), cts(6)
        
        print*, "--------------------------------------------------------------"
        write(*, *) "Satellite number: ", this%number
        write(*, *) "Name: ", trim(this%name)
        write(*, *) "Mass: ", this%mass, " kg"
        write(*, *) "Initial epoch: ", this%states%t%toDate()
        
        ! Show initial state in human-readable form
        if (this%eleType == STATE_TYPE_NONSINGULAR_1) then
            nsg1 = this%states%elements
            call nsg12kpl(nsg1, kpl)
            call kpl2cts(kpl, cts)
            call cts_normal2SI(cts, cts)
            write(*, "(A,F12.3,1x,F12.3,1x,F12.3)") " Initial position (km, J2000): ", cts(1:3)
            write(*, "(A,F10.6,1x,F10.6,1x,F10.6)") " Initial velocity (km/s, J2000): ", cts(4:6)
        endif
        
        write(*, *) "Propagator: ", trim(this%propagator)
        print*, "------ Dynamic model ------"
        write(*, *) "Gravity degree/order = ", this%opMod%degree, this%opMod%order

        ! Show active force models
        if (this%opMod%degree >= 2) then
            if (this%opMod%order >= 1) then
                write(*, "(A,I0,A,I0,A,I0,A)") "  - Earth gravity field (", &
                    &   this%opMod%degree, " x ", this%opMod%order, ")"
            else
                write(*, "(A,I0,A)") "  - Earth gravity field (zonal only, degree ", &
                    &   this%opMod%degree, ")"
            endif
        endif
        if (this%opMod%key_moonGravity == 1) print*, "  - Moon gravity"
        if (this%opMod%key_sunGravity == 1) print*, "  - Sun gravity"
        if (this%opMod%key_moonTide == 1) print*, "  - Solid Earth tides (Moon)"
        if (this%opMod%key_sunTide == 1) print*, "  - Solid Earth tides (Sun)"
        if (this%opMod%key_drag == 1) print*, "  - Atmospheric drag"
        if (this%opMod%key_srp == 1) print*, "  - Solar radiation pressure"
        if (this%opMod%key_ERP == 1) print*, "  - Earth radiation pressure (ERP)"
        if (this%opMod%key_postNewton == 1) print*, "  - Post-Newtonian effect"
        
        print*, "---------------------------"
        print*, "--------------------------------------------------------------"
    end subroutine

    !--------------------------------------------------------------------------
    ! SUBROUTINE: propagate
    ! Propagate satellite state to new epoch (updates this%states)
    !
    ! INPUTS:
    !   epoch - Target epoch for propagation
    !--------------------------------------------------------------------------
    subroutine sat_propagate(this, epoch)
        class(type_satellite), intent(inout) :: this
        type(type_t), intent(in) :: epoch
        double precision :: x0(6), x1(6), utc0, dt, kpl(6), nsg1(6), nsg2(6), cts(6)
        type(type_t) :: t0
        integer :: opEleType
        
        t0 = this%states%t
        utc0 = t0%day + t0%sec/86400d0
        dt = timeMinus(epoch, t0)
        x0 = this%states%elements
        
        ! Convert to appropriate element type for propagation
        if (this%opMod%usePV == 1) then
            opEleType = STATE_TYPE_CARTESIAN
        else
            opEleType = STATE_TYPE_NONSINGULAR_1
            if (this%eleType == STATE_TYPE_NONSINGULAR_1) then 
                x0 = this%states%elements
            elseif (this%eleType == STATE_TYPE_KEPLERIAN) then
                kpl = this%states%elements
                call kpl2nsg1(kpl, x0)
            elseif (this%eleType == STATE_TYPE_CARTESIAN) then
                cts = this%states%elements
                call cts2nsg1(cts, x0)
            elseif (this%eleType == STATE_TYPE_NONSINGULAR_2) then
                nsg2 = this%states%elements
                call nsg22nsg1(nsg2, x0)
            else
                print*, "Error: unsupported element type @propagate."
                STOP
            endif
        endif
        
        ! Perform propagation based on selected method
        if (this%propagator == 'numerical') then
            call NumOP(this%opMod%forceSwitchArray(), this%opMod%degree, this%opMod%order, &
                &   this%opMod%Cd, this%opMod%amrDrag, this%opMod%F107, this%opMod%Ap, &
                &   this%opMod%reflectivity, this%opMod%amrSRP, &
                &   this%opMod%SRP_model, this%opMod%ECOM, &
                &   this%opMod%step, &
                &   this%opMod%orderMoon, &
                &   this%opMod%Mercury, &
                &   this%opMod%Venus, &
                &   this%opMod%Mars, &
                &   this%opMod%Jupiter, &
                &   this%opMod%Saturn, &
                &   this%opMod%Uranus, &
                &   this%opMod%Pluto, &
                &   this%opMod%key_ERP, &
                &   this%opMod%reflectivityERP, &
                &   this%opMod%amrERP, &
                &   opEleType, &
                &   this%number, &
                &   x0, utc0, dt, x1)
            
            ! Update satellite state
            this%states%t = epoch
            if (this%opMod%usePV == 0) then
                select case (this%eleType) 
                case (STATE_TYPE_KEPLERIAN)
                    call nsg12kpl(x1, kpl)
                    this%states%elements = kpl
                case (STATE_TYPE_CARTESIAN)
                    call nsg12cts(x1, cts)
                    this%states%elements = cts
                case (STATE_TYPE_NONSINGULAR_1)
                    this%states%elements = x1
                case (STATE_TYPE_NONSINGULAR_2)
                    call nsg12nsg2(x1, nsg2)
                    this%states%elements = nsg2
                end select
            else
                select case (this%eleType) 
                case (STATE_TYPE_KEPLERIAN)
                    call cts2kpl(x1, kpl)
                    this%states%elements = kpl
                case (STATE_TYPE_CARTESIAN)
                    this%states%elements = x1
                case (STATE_TYPE_NONSINGULAR_1)
                    call cts2kpl(x1, kpl)
                    call kpl2nsg1(kpl, nsg1)
                    this%states%elements = nsg1
                case (STATE_TYPE_NONSINGULAR_2)
                    call cts2kpl(x1, kpl)
                    call kpl2nsg2(kpl, nsg2)
                    this%states%elements = nsg2
                end select
            endif

        elseif (this%propagator == 'two-body') then
            print*, "Two-body propagator not implemented @propagate."
        elseif (this%propagator == 'SGP4') then
            print*, "SGP4 propagator not implemented @propagate."
        else
            print*, "Error @propagate: unrecognized propagator = ", trim(this%propagator)
            stop
        endif
    end subroutine

    !--------------------------------------------------------------------------
    ! SUBROUTINE: freeTrajectory
    ! Release trajectory memory and reset nEph to 0
    !
    ! DESCRIPTION:
    !   This method deallocates the trajectory array if it is allocated,
    !   and sets nEph to 0. It provides a clean way to free memory when
    !   trajectory data is no longer needed.
    !
    ! USAGE:
    !   call sat%freeTrajectory()
    !--------------------------------------------------------------------------
    subroutine freeTrajectory(this)
        class(type_satellite), intent(inout) :: this
        
        ! Check if trajectory is allocated
        if (allocated(this%trajectory)) then
            deallocate(this%trajectory)
            print*, "Trajectory memory released for satellite: ", trim(this%name)
        else
            print*, "Trajectory not allocated for satellite: ", trim(this%name)
        endif
        
        ! Reset nEph to 0
        this%nEph = 0
    end subroutine freeTrajectory

    !--------------------------------------------------------------------------
    ! NOTE: the local file readers readMatrix1/readMatrix2 that used to live in
    ! this module (and in module_access) skipped a fixed number of header lines.
    ! They are replaced by readDataMatrix in module_toolKits, which skips every
    ! line starting with a comment marker ('#', '*', '!', '！') and therefore
    ! does not need to know the header length.
    !--------------------------------------------------------------------------

end module module_satellite
