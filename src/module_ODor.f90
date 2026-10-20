    
    module module_ODor
    use module_satelliteScenario
    use module_SLE
    use module_residual
    use module_perturbationSRP
    use module_fastPrediction
    use module_access         ! OBS_TYPE_RANGING, key_clock_error
    use module_clock          ! type_clock: 钟差组的播种与写回
    use CMD_ProgressBar       ! 迭代进度条(全局单例 Progress)

    integer, parameter :: OD_METHOD_BATCH = 1
    integer, parameter :: OD_METHOD_EKF = 2
    integer, parameter :: OD_METHOD_SQUENTIAL = 3
    
    type :: type_varGroup
        integer :: nVar ! number of variables in this varGroup
        integer, allocatable :: idx(:) ! index in x0
        integer :: label ! type of variable group, 
        ! label = 
        ! 1 -> states; 
        ! 2 -> states + amr of drag; 
        ! 3 -> states + amr of SRP; 
        ! 4 -> states + amr of drag + amr of SRP; 
        ! 5 -> states + SRP ECOM m paras; m = nVar - 6;
        ! 11, 12, 13 -> clock of ONE object (satellite or station): the first nVar
        !               coefficients (nVar = 1..3) of
        !               bias(t) = bias0 + drift*dt + drift2*dt^2, dt = (t - this%t) [s]
        ! 1000 -> access bias (未接线)
        integer :: origin ! 组所属对象的编号: 卫星/地面站下标(label 1~5/11~13), access 编号(label 1000)
        integer :: objType = OBJECT_TYPE_SATELLITE ! 仅 label 11~13 使用: OBJECT_TYPE_SATELLITE / _STATION
    endtype

    !---------------------------------------------------------------------------
    ! TYPE: type_STM —— 状态转移/灵敏度矩阵的时间序列
    !
    ! 定轨需要 dX(t)/dp, p 是一颗卫星的状态量加上该变量组的参数(面质比 / ECOM)。
    !
    ! 为什么不现算: label 5 的数值 STM(numericalSTM_ECOM)是有限差分, 一个历元
    ! 要跑 15 次完整传播(~0.3 s); 若在"逐观测历元"的循环里现算, 1 天 30 s 采样
    ! 要 1.7 万次、单次迭代以小时计。
    !
    ! 做法: 先在一条比观测粗得多的网格上把 STM 整条算出来存下, 定轨时按观测
    ! 历元做 4 点 Lagrange 插值取用 —— 各 label 统一走这条路, calB 里不再分支
    ! (label 1~4 用解析 STM, label 5 用数值 STM, 由 stm_series_init 分派)。
    ! STM 随时间的变化主要是轨道周期尺度(~12.9 h), 远慢于观测采样, 网格取
    ! 600 s 时插值误差可忽略(实测 7e-4 相对)。
    !
    ! 存的内容(对某个变量组):
    !   stm(k, 1:6, 1:nVar)     = dX(t_k)/dp        (状态量对全部参数的偏导)
    !   stm(k, 7:nVar, 7:nVar)  = I                 (参数本身是常数)
    !   其余元素为 0。calB 只用到 stm(:, 1:6, :)。
    !---------------------------------------------------------------------------
    type :: type_STM
        logical :: valid = .false.  ! 序列是否已生成
        integer :: nEpoch = 0       ! 网格点数
        integer :: nVar = 0         ! 该组的待估参数个数
        type(type_t) :: t0          ! 网格起始历元
        double precision :: stepSec = 0d0 ! 网格间隔 [s]
        double precision, allocatable :: stm(:, :, :) ! (nEpoch, nVar, nVar)
    contains
        procedure :: init => stm_series_init
        procedure :: get => stm_series_get
        procedure :: free => stm_series_free
    end type type_STM
    
    type :: type_ODor
        character(len=200) :: introduction
        integer :: nVar ! number of the variables to be determined
        type(type_t) :: t ! epoch of x
        double precision, allocatable :: x(:, :) ! values of variables
        double precision, allocatable :: P(:, :) ! covariance matrix of x0
        integer :: nGroup ! number of the varGroups
        type(type_varGroup), allocatable :: group(:) ! groups of variables
        integer :: nSat
        integer :: satIdx(100) = 0
        integer :: nStation
        integer :: stationIdx(100) = 0
        integer :: nAccess
        integer :: accessIdx(100) = 0
        integer :: nFile(100) ! number of observational data for each access
        character*100 :: dataFN(100, 100) ! filenames of obs data for each access
        integer :: method = OD_METHOD_BATCH
        ! 各变量组的 STM 时间序列(每个 label 都预计算; 每次迭代重建一次)
        type(type_STM), allocatable :: STMseries(:)
        double precision :: stmStepSec = 600d0 ! STM 序列的网格间隔 [s]
        ! 迭代进度显示: 轨迹阶段复用卫星自带的 showProgressBar(单次迭代里最贵的一段),
        ! 观测遍历阶段用本模块的进度条(以"已处理历元数"计)。默认关闭, 与
        ! type_opMod%showProgressBar 的默认一致。
        logical :: showProgress = .false.
    contains
        procedure :: init => initialize_ODor ! initialize ODor
        procedure :: solve => orbitDetermination_batch ! execute orbit determination (batch process)
        procedure :: invoke ! update orbit with one data using the method of EKF
        procedure :: saveEph
        
        procedure :: SS2ODor ! transfer states from satellite scenario to orbit determinator
        procedure :: ODor2SS ! transfer states from orbit determinator to satellite scenario
        procedure, private :: ODor_propagate ! propagate x and covariance
        procedure, private :: buildSTMseries ! (re)build the STM time series of each group
    endtype

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initialize ODor
    ! inputs:
    !   ss          type_satelliteScenario     satellites, stations and access inside
    !   t0          type_t                     epoch of x0
    !   nGroup      integer                    number of variable groups
    !   Groups      type_varGroup(nGroup)     variable groups
    !   method      integer                    OD method
    !   P           double precision(nVar, nVar)       initial covariance matrix(optional), is needed in EKF and sequential method
    !   intro       character(len=*)        introduction of this ODor scenario(optional)
    subroutine initialize_ODor(this, ss, t0, nGroup, Groups, method, P, intro)
    implicit none
    class(type_ODor) :: this
    type(type_satelliteScenario),intent(inOut) :: ss
    type(type_t), intent(in) :: t0
    integer, intent(in) :: nGroup
    type(type_varGroup), intent(in) :: Groups(nGroup)
    integer, intent(in) :: method
    double precision, intent(in), optional, allocatable :: P(:, :)
    character(len=*), intent(in), optional :: intro
    integer :: i, stationNO
    !----------------------------------------------
    this%nSat = ss%nSat
    this%satIdx(1:this%nSat) = (/ (i, i = 1, ss%nSat) /)
    this%nStation = ss%nStation
    this%stationIdx(1:this%nStation) = (/ (i, i = 1, ss%nStation) /)
    this%nAccess = ss%nAccess
    this%accessIdx(1:this%nAccess) = (/ (i, i = 1, ss%nAccess) /)

    this%t = t0

    this%nGroup = nGroup
    allocate(this%group(this%nGroup))
    this%group = Groups
    if (allocated(this%STMseries)) deallocate(this%STMseries)
    allocate(this%STMseries(this%nGroup))
    this%nVar = 0
    do i = 1, this%nGroup
        this%nVar = this%nVar + this%group(i)%nVar
    enddo

    this%method = method

    ! 校验待估参数组(含钟差组的额外约定)，不合法立即停下
    call validateVarGroups(this, ss)

    allocate(this%x(this%nVar, 1))
    this%x = 0d0
    if (present(P)) then
        this%P = P
    else
        allocate(this%P(this%nVar, this%nVar))
        this%P = 0d0
    endif
    
    ! 
    if (present(intro)) then
        this%introduction = intro
    else
        this%introduction = "Orbit Determination"
    endif

    ! copy initial states from satellite scenario to OD scenario
    call this%SS2ODor(ss)

    ! if there is OBD_TYPE_AZIMUTH_ELEVATION access, 
    ! transfer it to OBS_TYPE_RA_DEC
    do i = 1, this%nAccess
        if (ss%access(i)%obsType == OBS_TYPE_AZIMUTH_ELEVATION &
        .and. ss%access(i)%toType == OBJECT_TYPE_STATION) then
            stationNO = ss%access(i)%toNum
            call AzimuthElevation2RaDec(ss%station(stationNO), ss%access(i))
        endif
    enddo

    endsubroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SUBROUTINE: validateVarGroups
    ! 校验待估参数组，不合法就报错停下(用本模块既有的 print + stop 风格)。
    !
    ! 为什么必须在这里前置检查，而不是交给求解器兜底：
    !   1) idx 越界或重复会静默写坏内存(SS2ODor/ODor2SS/calB 都按 idx 直接寻址)；
    !   2) 钟差组的列若恒为 0(所属对象不参与任何测距链路)，法方程是"对角非零但
    !      奇异"—— LUDCMP 只在整行元素都小于 1.5e-16 时才返回奇异标志，这种情况
    !      下它把零主元替换成 TINY，返回放大约 1e15 的逆，得到的是跑飞的垃圾解，
    !      而不是一个可诊断的报错。
    !
    ! 输入:
    !   ss - 卫星场景(判断钟差组所属对象是否参与测距链路)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine validateVarGroups(this, ss)
    implicit none
    class(type_ODor), intent(in) :: this
    type(type_satelliteScenario), intent(in) :: ss
    type(type_clock) :: clk
    integer :: i, j, k, label, nVar, origin, objType, nObj
    integer, allocatable :: nUsed(:)
    logical :: hasRanging

    ! ---- 通用: idx 必须落在 1..nVar 内，且各组的 idx 互不重复 ----
    allocate(nUsed(this%nVar))
    nUsed = 0
    do i = 1, this%nGroup
        nVar = this%group(i)%nVar
        if (.not. allocated(this%group(i)%idx)) then
            print*, "error @validateVarGroups: group ", i, " has no idx allocated"
            stop
        endif
        if (size(this%group(i)%idx) /= nVar) then
            print*, "error @validateVarGroups: group ", i, &
                &   " idx size (", size(this%group(i)%idx), ") /= nVar (", nVar, ")"
            stop
        endif
        do j = 1, nVar
            k = this%group(i)%idx(j)
            if (k < 1 .or. k > this%nVar) then
                print*, "error @validateVarGroups: group ", i, " idx out of range: ", k
                stop
            endif
            if (nUsed(k) /= 0) then
                print*, "error @validateVarGroups: idx ", k, " is used by group ", &
                    &   nUsed(k), " and group ", i
                stop
            endif
            nUsed(k) = i
        enddo
    enddo
    deallocate(nUsed)

    ! ---- 钟差组(label 11~13)的额外约定 ----
    do i = 1, this%nGroup
        label   = this%group(i)%label
        nVar    = this%group(i)%nVar
        origin  = this%group(i)%origin
        objType = this%group(i)%objType
        if (label < 11 .or. label > 13) cycle

        ! 只支持批处理: EKF/序贯走 invoke/numericalSTM，根本不过 calB
        if (this%method /= OD_METHOD_BATCH) then
            print*, "error @validateVarGroups: clock group ", i, " (label ", label, &
                &   ") is only supported by OD_METHOD_BATCH; method = ", this%method
            stop
        endif

        if (nVar < 1 .or. nVar > 3) then
            print*, "error @validateVarGroups: clock group ", i, &
                &   " needs nVar = 1..3 (bias0/drift/drift2), got ", nVar
            stop
        endif

        if (objType == OBJECT_TYPE_SATELLITE) then
            nObj = this%nSat
        elseif (objType == OBJECT_TYPE_STATION) then
            nObj = this%nStation
        else
            print*, "error @validateVarGroups: clock group ", i, &
                &   " has invalid objType = ", objType
            stop
        endif
        if (origin < 1 .or. origin > nObj) then
            print*, "error @validateVarGroups: clock group ", i, &
                &   " origin out of range: ", origin, " (objType ", objType, ")"
            stop
        endif

        ! 同一对象至多一个钟差组
        do j = 1, i - 1
            if (this%group(j)%label >= 11 .and. this%group(j)%label <= 13) then
                if (this%group(j)%objType == objType .and. this%group(j)%origin == origin) then
                    print*, "error @validateVarGroups: object (objType ", objType, &
                        &   ", origin ", origin, ") has more than one clock group"
                    stop
                endif
            endif
        enddo

        ! 该对象必须至少参与一条 RANGING 链路，否则钟差列恒为 0
        hasRanging = .false.
        do j = 1, this%nAccess
            k = this%accessIdx(j)
            if (ss%access(k)%obsType /= OBS_TYPE_RANGING) cycle
            if (ss%access(k)%fromType == objType .and. ss%access(k)%fromNum == origin) hasRanging = .true.
            if (ss%access(k)%toType   == objType .and. ss%access(k)%toNum   == origin) hasRanging = .true.
        enddo
        if (.not. hasRanging) then
            print*, "error @validateVarGroups: clock group ", i, " -- object (objType ", &
                &   objType, ", origin ", origin, ") takes part in no RANGING access."
            print*, "   (only ranging observations carry a clock term;"// &
                &   " any other observable gives an all-zero column)"
            stop
        endif

        ! 插值方法的钟: 解析系数与状态表是两个不同的函数，播种会不自洽
        if (objType == OBJECT_TYPE_SATELLITE) then
            clk = ss%sat(origin)%clock
        else
            clk = ss%station(origin)%clock
        endif
        if (clk%calcMethod == CLOCK_METHOD_INTERPOLATION) then
            print*, "error @validateVarGroups: clock group ", i, &
                &   " -- the object's clock uses CLOCK_METHOD_INTERPOLATION, but clock"
            print*, "   estimation works on the analytic (FORMULA) model"
            stop
        endif
    enddo

    end subroutine validateVarGroups

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! copy initial values of the OD variables from satelliteScenario to ODorcenario
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine SS2ODor(this, ss)
    use module_atmosphere
    use module_time
    implicit none
    class(type_ODor) :: this
    type(type_satelliteScenario), intent(in) :: ss
    integer :: label, origin, i, nVar, nECOM
    double precision :: beta, amrDrag
    double precision :: kappa
    double precision :: c(3)
    !-------------------------------------------------------------------

    ! 历元同步: 先对整个组的集合做一遍, 再逐组播种。
    ! 原来这一步在播种循环内部进行, 会在遍历中途改 this%t, 使后面的组看到的基准
    ! 历元依赖组序; 钟差组的参数基准也是 this%t, 不能留这种隐式依赖。
    do i = 1, this%nGroup
        label  = this%group(i)%label
        origin = this%group(i)%origin
        if (label < 1 .or. label > 3) cycle
        if (ss%eleType /= ss%sat(origin)%eleType) then
            print*, "error @SS2ODor"
            print*, "eleTypes are different!"
            print*, "ss%eleType = ", ss%eleType
            print*, "sat%eleType = ", ss%sat(origin)%eleType
        endif
        if (dabs(timeMinus(this%t, ss%sat(origin)%states%t)) > 1d-15) then
            this%t = ss%sat(origin)%states%t
        endif
    enddo

    ! copy states
    do i = 1, this%nGroup
        nVar = this%group(i)%nVar
        label = this%group(i)%label
        origin = this%group(i)%origin

        select case (label)
        case(1) ! elements
            this%x(this%group(i)%idx(1:6), 1) = ss%sat(origin)%states%elements
        case(2) ! elements + amrDrag
            this%x(this%group(i)%idx(1:6), 1) = ss%sat(origin)%states%elements
            amrDrag = ss%sat(origin)%opMod%amrDrag
            beta = amrDrag*mass_unit/(length_unit*1d3)**2*ss%sat(origin)%opMod%Cd*rho_p0
            this%x(this%group(i)%idx(7), 1) = beta
        case(3) ! elements + amrSRP
            this%x(this%group(i)%idx(1:6), 1) = ss%sat(origin)%states%elements
            kappa = 1d0 + ss%sat(origin)%opMod%reflectivity
            this%x(this%group(i)%idx(7), 1) = ss%sat(origin)%opMod%amrSRP &
                &   *mass_unit/(length_unit*1d3)**2 *kappa *rho_s ! m^2/kg -> normalized acceleration
        case(4) ! elements + amrDrag + amrSRP
            this%x(this%group(i)%idx(1:6), 1) = ss%sat(origin)%states%elements
            amrDrag = ss%sat(origin)%opMod%amrDrag
            beta = amrDrag*mass_unit/(length_unit*1d3)**2*ss%sat(origin)%opMod%Cd*rho_p0
            this%x(this%group(i)%idx(7), 1) = beta
            kappa = 1d0 + ss%sat(origin)%opMod%reflectivity
            this%x(this%group(i)%idx(8), 1) = ss%sat(origin)%opMod%amrSRP &
                &   *mass_unit/(length_unit*1d3)**2 *kappa *rho_s ! m^2/kg -> normalized acceleration
        case(5) ! elements + ECOMs
            nECOM = nVar - 6
            this%x(this%group(i)%idx(1:6), 1) = ss%sat(origin)%states%elements
            this%x(this%group(i)%idx(7:nVar), 1) = ss%sat(origin)%opMod%ECOM(1:nECOM)
        case(11, 12, 13) ! clock coefficients of one object (satellite or station)
            ! 参数基准历元是 this%t, 所以要把钟模型精确重基准化到 this%t —— 直接取
            ! clk%bias0/drift/drift2 只在 clk%t0 == this%t 时才等价。
            ! 这里只读、不改对象的钟: 重基准化不改变钟函数, 模型观测与播种值自洽;
            ! 对象的 clk%t0 留到 ODor2SS 写回时统一改成 this%t。
            if (this%group(i)%objType == OBJECT_TYPE_STATION) then
                call ss%station(origin)%clock%getModelCoeffs(this%t, c)
            else
                call ss%sat(origin)%clock%getModelCoeffs(this%t, c)
            endif
            this%x(this%group(i)%idx(1:nVar), 1) = c(1:nVar)
        case(1000) ! access bias
            this%x(this%group(i)%idx(1:nVar), 1) = ss%access(origin)%bias(1:nVar)
        case default
            print*, "error @SS2ODor"
            print*, "invalid label, label = ", label
            stop
        endselect
    enddo
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! copy initial values from ODorcenario to satelliteScenario
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ODor2SS(this, ss)
    use module_atmosphere
    implicit none
    class(type_ODor) :: this
    type(type_satelliteScenario), intent(inout) :: ss
    integer :: label, origin, i, nVar, nECOM
    double precision :: amrDrag, beta, kappa
    double precision :: c(3)
    !-------------------------------------------------------------------
    
    ! copy states
    do i = 1, this%nGroup
        nVar = this%group(i)%nVar
        label = this%group(i)%label
        origin = this%group(i)%origin
        
        ! if origin is satellite
        if (label > 0 .and. label < 4) then
            
            if (ss%eleType /= ss%sat(origin)%eleType) then
                print*, "error @ODor2SS"
                print*, "eleTypes are different!"
                print*, "ss%eleType = ", ss%eleType
                print*, "sat%eleType = ", ss%sat(origin)%eleType
            endif
        
        endif
        
        
        select case (label)
        case(1)
            ss%sat(origin)%states%elements = this%x(this%group(i)%idx(1:6), 1)
        case(2)
            ss%sat(origin)%states%elements = this%x(this%group(i)%idx(1:6), 1)
            beta = this%x(this%group(i)%idx(7), 1)
            amrDrag = beta/ss%sat(origin)%opMod%Cd/rho_p0/mass_unit*(length_unit*1d3)**2
            ss%sat(origin)%opMod%amrDrag = amrDrag
        case(3)
            ss%sat(origin)%states%elements = this%x(this%group(i)%idx(1:6), 1)
            kappa = 1d0 + ss%sat(origin)%opMod%reflectivity
            ss%sat(origin)%opMod%amrSRP = this%x(this%group(i)%idx(7), 1) / &
            &   ( mass_unit/(length_unit*1d3)**2 *kappa *rho_s ) ! normalized acc -> m^2/kg
        case(4)
            ss%sat(origin)%states%elements = this%x(this%group(i)%idx(1:6), 1)
            beta = this%x(this%group(i)%idx(7), 1)
            amrDrag = beta/ss%sat(origin)%opMod%Cd/rho_p0/mass_unit*(length_unit*1d3)**2
            ss%sat(origin)%opMod%amrDrag = amrDrag
            
            kappa = 1d0 + ss%sat(origin)%opMod%reflectivity
            ss%sat(origin)%opMod%amrSRP = this%x(this%group(i)%idx(8), 1) / &
            &   ( mass_unit/(length_unit*1d3)**2 *kappa *rho_s ) ! normalized acc -> m^2/kg
        case(5)
            nECOM = nVar - 6
            ss%sat(origin)%states%elements = this%x(this%group(i)%idx(1:6), 1)
            ss%sat(origin)%opMod%ECOM(1:nECOM) = this%x(this%group(i)%idx(7:nVar), 1)
        case(11, 12, 13) ! clock coefficients of one object (satellite or station)
            ! 组参数只是钟模型的前 nVar 个系数, 其余系数保持"重基准化到 this%t 之后"
            ! 的原值 —— 否则换个基准历元就等于悄悄改了钟函数。
            if (this%group(i)%objType == OBJECT_TYPE_STATION) then
                call clockInstall(ss%station(origin)%clock, this%t, &
                    &   this%x(this%group(i)%idx(1:nVar), 1), nVar)
            else
                call clockInstall(ss%sat(origin)%clock, this%t, &
                    &   this%x(this%group(i)%idx(1:nVar), 1), nVar)
            endif
        case(1000) ! access bias
            ss%access(origin)%bias(1:nVar) = this%x(this%group(i)%idx(1:nVar), 1)
        case default
            print*, "error @ODor2SS"
            print*, "invalid label, label = ", label
            stop
        endselect
    enddo

    endsubroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SUBROUTINE: clockInstall
    ! 把钟差组的估计值装进对象钟: 读出重基准化到 tRef 的三系数, 用 xVal 覆盖前
    ! nVar 个, 再整组装回(setModel 会同时把 calcMethod 置为公式方法)。
    !
    ! 输入:
    !   clk   - 目标钟 (inout)
    !   tRef  - 参数基准历元
    !   xVal  - 本组估计出的前 nVar 个系数
    !   nVar  - 系数个数 (1..3)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine clockInstall(clk, tRef, xVal, nVar)
    implicit none
    type(type_clock), intent(inout) :: clk
    type(type_t), intent(in) :: tRef
    double precision, intent(in) :: xVal(:)
    integer, intent(in) :: nVar
    double precision :: c(3)

    call clk%getModelCoeffs(tRef, c)
    c(1:nVar) = xVal(1:nVar)
    call clk%setModel(tRef, c(1), c(2), c(3))
    end subroutine clockInstall
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! FUNCTION: resUnitCode
    ! 观测残差的显示单位。残差在内部是归一化量, 直接打印是 1e-9 这种没法读的数,
    ! 屏幕/日志上按观测量类型换算成直观单位(换算在 resRMSstr 里做):
    !   1 = 角度 -> arcsec      2 = 长度 -> m
    !   3 = 速度 -> m/s         4 = 无量纲/未知 -> 保持归一化
    !
    ! 为什么按分量而不是按链路: 同一条链路的各分量可能不同量纲(POSITION_VELOCITY
    ! = 位置 + 速度), 把它们平方相加得到的那个"RMS"没有物理意义, 所以按单位分组,
    ! 每个单位组各给一个数。
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer function resUnitCode(obsType, jc)
    implicit none
    integer, intent(in) :: obsType ! 观测量类型 (OBS_TYPE_*)
    integer, intent(in) :: jc      ! 该观测量第几个分量 (1..npe)

    select case (obsType)
    case (OBS_TYPE_RA_DEC, OBS_TYPE_AZIMUTH_ELEVATION)
        resUnitCode = 1                        ! RA/Dec, Az/El (弧度)
    case (OBS_TYPE_RANGING, OBS_TYPE_POSITION)
        resUnitCode = 2                        ! Range, (x, y, z)
    case (OBS_TYPE_RANGERATE)
        resUnitCode = 3                        ! Range rate
    case (OBS_TYPE_POSITION_VELOCITY)
        if (jc <= 3) then
            resUnitCode = 2                    ! Position
        else
            resUnitCode = 3                    ! Velocity
        endif
    case default
        resUnitCode = 4
    end select
    end function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SUBROUTINE: resRMSstr
    ! 把已算过残差的前 nAcc 条链路(全局 res%acc)按显示单位分组, 每组统计一个 RMS
    ! 并拼成一行文字, 如 "1.234E-02 m", 或 "5.60E-06 arcsec, 1.20E-05 m/s"。
    ! 收敛判据用的仍是 res%calRMS 给出的归一化总 RMS, 这里只影响显示。
    ! 注意: 字符串一律 ASCII —— 中文在 GBK 终端/日志里会乱码。
    !
    ! 输入:
    !   nAcc - 参与统计的链路条数(本轮已处理到第几条)
    ! 输出:
    !   str  - 形如 "3.21E-02 m" 的文字; 无有效数据时给提示语
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine resRMSstr(nAcc, str)
    implicit none
    integer, intent(in) :: nAcc
    character(len=*), intent(out) :: str
    ! 归一化量 -> 显示单位的换算系数(顺序与 resUnitCode 的编码一致)
    double precision, parameter :: sc(4) = (/ RAD2ARCSEC, 1d0/u_meter, &
        &   1d0/u_meter/TIME_UNIT, 1d0 /)
    character(len=6), parameter :: unit(4) = (/ "arcsec", "m     ", "m/s   ", "n.u.  "/)
    double precision :: sumR(4), val
    integer :: ng(4), ic, i, j, k, jc
    character(len=32) :: buf

    sumR = 0d0
    ng = 0
    do i = 1, nAcc
        do j = 1, res%acc(i)%nValid
            k = res%acc(i)%validIdx(j)
            do jc = 1, res%acc(i)%npe
                ic = resUnitCode(res%acc(i)%obsType, jc)
                sumR(ic) = sumR(ic) + res%acc(i)%res(k, jc)**2
                ng(ic) = ng(ic) + 1
            enddo
        enddo
    enddo

    str = ""
    do ic = 1, 4
        if (ng(ic) == 0) cycle
        val = dsqrt(sumR(ic)/dble(ng(ic)))*sc(ic)
        write(buf, "(ES9.3,1x,a)") val, trim(unit(ic))
        if (len_trim(str) == 0) then
            str = trim(buf)
        elseif (len_trim(str) + len_trim(buf) + 2 <= len(str)) then
            str = trim(str)//", "//trim(buf)
        endif
    enddo
    if (len_trim(str) == 0) str = "(no valid residual)"
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine orbitDetermination_batch(this, ss)
    implicit none
    class(type_ODor), intent(inOut) :: this
    type(type_satelliteScenario), intent(inOut) :: ss
    type(type_SLE) :: SLE
    integer :: nSat, nVar, iIter, iac, it, npe
    integer, allocatable :: list(:)
    type(type_TimeSeries) :: tso
    type(type_access) :: accCal
    double precision, allocatable :: W(:, :), B(:, :), Y(:, :)
    integer :: aci, nData, tIdx
    integer :: lastNpe = 0
    double precision :: RMS, RMSprev
    double precision :: cts(6), kpl(6)
    
    double precision :: cts0(6), cts1(6), cts1p(6)
    double precision :: dkpl(6), kpl0(6), kpl1(6), kpl1p(6)
    double precision :: Bn(6, 6)
    integer :: i, j
    type(type_t) :: t0, t1
    type(type_satellite) :: sat1, sat1p
    ! 进度显示
    integer :: nObsTot, nObsDone, nStride, tc0, tc1, tcAll0, tcRate
    double precision :: tElap, tIter
    character(len=160) :: chr, chs
    character(len=80) :: rmsStr
    logical, allocatable :: barSaved(:)
    
    !------------------------------------------------------------------
    
    nSat = this%nSat
    if (nSat <= 0) then
        print*, "error @orbit_determination_batch"
        print*, "this%nSat = ", this%nSat
        stop
    endif
    allocate(list(nSat))
    list(1:nSat) = this%satIdx(1:nSat)
    nVar = this%nVar
    
    ! 1 get initial values
    call SLE%initialize(this%nVar)
    SLE%x = 0d0
    SLE%P = this%P
    
    ! initialize residual
    call res%allocRes(this%nAccess)
    
    open(901, file = "data/output/log_batchOD.txt")
    
    print*, ""
    write(*, "(' o_o iteration ',I2)") 0
    write(901, "('> iteration ',I2)") 0
    write(901, *) "x = ", this%x

    ! ---- 进度显示(可选) ----
    ! 观测总数: 各链路历元数之和, 循环前即可知(个别历元可能被判为粗差, 差异可忽略)。
    ! nStride: 每约 1% 刷新一次 —— 每个历元都刷新会白白拖慢程序(每次都是一次终端写)。
    nObsTot = 0
    do iac = 1, this%nAccess
        nObsTot = nObsTot + size(ss%access(this%accessIdx(iac))%timeTable)
    enddo
    nStride = max(1, nObsTot/100)
    if (this%showProgress) then
        allocate(barSaved(nSat))
        do i = 1, nSat
            barSaved(i) = ss%sat(list(i))%opMod%showProgressBar
            ss%sat(list(i))%opMod%showProgressBar = .true.
        enddo
        call system_clock(tcAll0, tcRate)     ! 整个 solve 的累计计时
    endif

    RMSprev = 1d30
    do iIter = 1, 15

        nObsDone = 0
        if (this%showProgress) call system_clock(tc0, tcRate)
        
        ! generate trajectory
        !call ss%generateSatellitetrajectory(nSat, list(1:nSat), .false.)
        do i = 1, nSat
            call ss%sat(list(i))%generateTrajectory(ss%ts%time)
        enddo

        ! 预计算各组的 STM 时间序列。必须在遍历观测数据之前做:
        ! 原先是每个观测历元现算一次有限差分 STM, 现在整条算完再插值取用。
        ! STM 沿参考轨迹线性化, 每次迭代参考轨迹都变, 所以放在迭代开头重建。
        call this%buildSTMseries(ss)

        ! 2 traverse the observational data
        do iac = 1, this%nAccess
            
            aci = this%accessIdx(iac)
            
            ! allocate matrix
            npe = ss%access(aci)%npe
            if (npe /= lastNpe) then
                if (allocated(B)) deallocate(B)
                if (allocated(Y)) deallocate(Y)
                if (allocated(W)) deallocate(W)
                allocate(B(npe, nVar))
                allocate(Y(npe, 1))
                allocate(W(npe, npe))
                lastNpe = npe
            endif
    
            ! create time Series of obs data, tso
            tso%time = ss%access(aci)%timeTable
            tso%n = size(ss%access(aci)%timeTable)
            
            ! calculate observational data
            accCal = ss%access(aci) ! use accCal to save the observational data
            call observationAccess(ss, accCal, tso, .true.)
            
            ! calculate residuals
            call Res%acc(iac)%calRes(ss%access(aci), accCal)
            call Res%calValidIdx(aci) ! except outliers
            call Res%calRMS ! calculate rms of residuals
            RMS = Res%RMS
            ! 显示用的 RMS: 按观测量类型换算成 m / m/s / arcsec(数值与归一化的
            ! RMS 一一对应, 只是换了单位)。此处含前 iac 条链路(累计值)。
            call resRMSstr(iac, rmsStr)

            ! 逐链路残差总是写日志。打开进度显示时屏幕上由"每轮一行汇总"代替
            ! (否则 6 条链路 × 15 轮就是 90 行噪音); 未打开时保持原来的屏幕输出,
            ! 其它定轨程序(多半只有一条链路)的行为不变。
            if (.not. this%showProgress) write(*, "(A,A)") "RMS = ", trim(rmsStr)
            write(901, "(A,A)") "RMS = ", trim(rmsStr)
    
            ! calculate W matrix
            call calW(ss%access(aci), npe, W)
            
            nData = tso%n
            
            do it = 1, Res%acc(iac)%nValid
                
                tIdx = Res%acc(iac)%validIdx(it)
                ! calculate the left and right term in the system of linear equations(SLE)
                ! B*x = Y
                call calB(this, ss, ss%access(aci), npe, nVar, tIdx, B)
                Y = reshape(res%acc(iac)%res(tIdx, 1:npe), (/npe, 1/))
    
                ! add data into SLE
                call SLE%addData(npe, nVar, B, Y, W)

                ! 进度: 以"已处理历元数"计, 每 nStride 个刷新一次
                nObsDone = nObsDone + 1
                if (this%showProgress .and. mod(nObsDone, nStride) == 0) then
                    call system_clock(tc1)
                    ! 条尾给"本次定轨累计已用时间": 观测遍历只占单轮的一小部分
                    ! (实测约 0.15 s / 1.88 s), 阶段级 ETA 没有意义且分辨率都不够;
                    ! 累计时间配合每轮汇总里的"本轮用时"就足以估计还要跑多久。
                    tElap = dble(tc1 - tcAll0)/dble(tcRate)
                    ! 终端输出一律 ASCII: 中文在 GBK 终端/重定向的日志里会乱码
                    ! 整行控制在 80 列以内: 条占 40+2+8 = 50 列, 前缀 + 条尾不超过
                    ! 30 列。超宽会让终端折行, 此时 char(13) 只能回到折行后的行首,
                    ! 上一行的残迹留在屏幕上。
                    write(chr, "(A,I2,A)") " iter ", iIter, " obs:"
                    Progress%prefix = trim(chr)
                    write(chs, "(I0,A,I0,A,F5.1,A)") nObsDone, "/", nObsTot, " ", &
                        &   tElap, "s"
                    call Progress%put(nObsTot, nObsDone, .true., trim(chs))
                endif
            enddo
        enddo

        ! 收尾: 把条补满并补一个换行, 收敛提前退出时条不会停在半格
        if (this%showProgress) then
            call Progress%put(nObsTot, nObsTot, .true.)
            call system_clock(tc1)
            tIter = dble(tc1 - tc0)/dble(tcRate)
        endif
    
        print*, ""
        write(*, "(' o_o iteration ',I2)") iIter
        write(901, *) ""
        write(901, "('> iteration ',I2)") iIter

        ! 本轮汇总(屏幕一行): 处理的星数/历元数 + 残差(直观单位) + 用时
        if (this%showProgress) then
            call resRMSstr(this%nAccess, rmsStr)   ! 全部链路的整体 RMS
            write(chs, "(A,I2,A,I0,A,I0,A,A,A,F6.2,A)") " iter ", iIter, ": ", nSat, &
                &   " sats, ", nObsTot, " obs, RMS = ", trim(rmsStr), ", ", tIter, " s"
            print "(A)", trim(chs)
        endif

        ! 4 solve SLE
        call SLE%solve
    
        ! 5 update satellite scenario, inlcuding change the initial states and ephemerides of satellites
        this%x = this%x + SLE%x
        this%P = SLE%P
        call this%ODor2SS(ss)
        write(901, *) "x = ", this%x
    
        !print*, "amrDrag = ", ss%sat(1)%opMod%amrDrag
        
        ! 6 exit iteration or not
        ! 判据用"本次迭代与上次迭代的整体残差 RMS"(calRMS 对全部 access 求和)。
        ! 原来的 lastRMS 是在 access 循环里被逐条覆盖的, 两侧比较的是两个不同
        ! 链路的残差, 与实际收敛与否无关, 会在残差还很大时提前退出。
        if (res%rms < 1d-9 .or. iIter > 1 .and. dabs(RMS - RMSprev) < 1d-12) then
            exit
        endif
        RMSprev = RMS
        !if (iIter == 10) exit
    
    enddo
    
    close(901)

    ! 还原卫星自己的进度条开关(OD 只是临时打开它来显示轨迹阶段)
    if (this%showProgress) then
        do i = 1, nSat
            ss%sat(list(i))%opMod%showProgressBar = barSaved(i)
        enddo
        deallocate(barSaved)
    endif

    ! free memory
    if (allocated(list)) deallocate(list)
    if (allocated(B)) deallocate(B)
    if (allocated(Y)) deallocate(Y)
    if (allocated(W)) deallocate(W)    
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine invoke(this, ss, epoch, iAcc, npe, nx, obsData, res)
    implicit none
    class(type_ODor) :: this
    type(type_satelliteScenario) :: ss
    type(type_t), intent(in) :: epoch
    integer, intent(in) :: iAcc, npe, nx
    double precision, intent(in) :: obsData(npe)
    double precision, intent(out) :: res(npe)
    
    double precision :: STM(nx, nx), DYDR(npe, 6), DYDX(npe, nx), DRDX(6, nx)
    type(type_timeSeries) :: tso
    integer :: i, fromNum, toNum, obsType
    double precision :: dY(npe, 1), x(6), PVs(6), PVo(6), W(npe, npe)
    double precision :: x_post(nx, 1), P_post(nx, nx)
    integer :: fromType, toType
    double precision :: t, longi, lati, obsData_cal(npe)
    double precision :: IntPNTM(nx, nx)
    double precision :: PVt(6), PVr(6)
    double precision :: dt
    !double precision :: sigma_ax = 1d-6
    !double precision :: sigma_ay = 1d-6
    !double precision :: sigma_az = 1d-6
    double precision :: sigma_ax = 1d-9
    double precision :: sigma_ay = 1d-9
    double precision :: sigma_az = 1d-9
    !-------------------------------------------------------------------
    
    ! calculate IntPNTM
    dt = (epoch%day - this%t%day)*86400d0 + epoch%sec - this%t%sec ! sec
    dt = dt/time_unit
    IntPNTM = 0d0
    call processNoise_cartesian(dt, sigma_ax, sigma_ay, sigma_az, IntPNTM(1:6, 1:6))

    ! >>>>> 1 time update
    call this%ODor_propagate(ss, epoch)

    t = epoch%val
    fromNum = ss%access(iAcc)%fromNum
    toNum = ss%access(iAcc)%toNum
    fromType = ss%access(iAcc)%fromType
    toType = ss%access(iAcc)%toType
    obsType = ss%access(iAcc)%obsType
    if (toType == OBJECT_TYPE_STATION) then
        longi = ss%station(toNum)%longitude
        lati = ss%station(toNum)%latitude
    endif
    
    if (fromType == OBJECT_TYPE_SATELLITE .and. toType == OBJECT_TYPE_STATION) then ! station-satellite observation
        ! >>>>>2 calculate observational data

        PVs = ss%sat(toNum)%states%elements ! object satellite
        PVo = ss%station(fromNum)%getPV(epoch) ! station satellite
        
        call observation_oneEpoch(t, longi, lati, PVs, PVo, obsType, npe, obsData_cal)
        obsData_cal(1:npe) = obsData_cal(1:npe) + ss%access(iAcc)%bias(1:npe)
        do i = 1, npe
            dY(i, 1) = obsData(i) - obsData_cal(i)
            res(i) = dY(i, 1) ! output
        enddo
        
        ! >>>>>3 matrix mapping states to observational data
        ! 3.1 DRDX
        call calDRDX(STATE_TYPE_CARTESIAN, PVs, DRDX)
        ! 3.2 DYDR
        call calDYDR(obsType, PVs, PVo, npe, DYDR)
        DYDX = matmul(DYDR, DRDX)
        
    elseif (fromType == OBJECT_TYPE_SATELLITE .and. toType == OBJECT_TYPE_SATELLITE) then ! inter-satellite observation
        
        PVr = ss%sat(toNum)%states%elements ! object satellite
        PVt = ss%sat(fromNum)%getPV(epoch)
        obsType = ss%access(iAcc)%obsType
        
        call observation_oneEpoch(t, longi, lati, PVt, PVr, obsType, npe, obsData_cal)
        
        obsData_cal(1:npe) = obsData_cal(1:npe) + ss%access(iAcc)%bias(1:npe) ! add bias 
        
        do i = 1, npe
            dY(i, 1) = obsData(i) - obsData_cal(i)
            res(i) = dY(i, 1) ! output
        enddo
        
        ! 3.2 DYDR
        call calDYDR(obsType, PVr, PVt, npe, DYDR)
        DYDX = 0d0
        DYDX(1:npe, 1:6) = DYDR
        DYDX(1, 6 + iAcc) = 1d0
        
    else
        print*, "unfinished part."
        stop
    endif
    
    ! >>>>> 4 kalman filter, update states and covariance
    ! define W
    W = 0d0
    do i = 1, npe
        W(i, i) = ss%access(iAcc)%errStd(i)**2
    enddo
    call KalmanFilter(nx, npe, dY, this%P, IntPNTM, DYDX, W, x_post, P_post)
    ! 4.2 transfer data to this
    this%x = this%x + x_post
    this%P = P_post
    
    ! 4.3 transfer data to ss
    call this%ODor2ss(ss)
    
    endsubroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! read observational data from files
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !subroutine readObsData(this, ss)
    !implicit none
    !class(type_ODor) :: this
    !type(type_satelliteScenario), intent(inout) :: ss
    !integer :: i, j, acIdx
    !character*200 :: fn
    !do i = 1, this%nAccess
    !    acIdx = this%accessIdx(i)
    !    if (this%nFile(i) > 0) then
    !        do j = 1, this%nFile(i)
    !            fn = trim(this%dataFN(i, j))
    !            call ss%access(acIdx)%readData(trim(fn))
    !        enddo
    !    endif
    !enddo
    !endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calB(this, ss, acO, m, n, it, B)
    implicit none
    type(type_ODor), intent(in) :: this
    type(type_satelliteScenario), intent(in) :: ss
    type(type_access), intent(in) :: acO
    integer, intent(in) :: m, n, it
    double precision, intent(out) :: B(m, n)
    type(type_t) :: t1
    integer :: fromType, toType, fromNum, toNum, eleType
    integer :: i, nVar, j
    logical :: keyDrag, keySRP
    logical :: isFromSat, isToSat
    double precision :: x0(6), x1(6), PVo(6), PVs(6)
    double precision :: STM(20, 20), dt, DRDX(6, 6), DYDR(m, 6)
    integer :: obsType, satIdx, idx0
    !-------------------------------------------------------------------
    
    t1 = acO%timeTable(it)
    ! B 是 intent(out): 只有参与本链路的那些组的列会被填。其余列必须显式清零 ——
    ! 否则 SLE%addData 会把整行(含上一轮 access 残留的、别的卫星的偏导)当成有效
    ! 观测累加进法方程, 相当于给不参与的卫星加上虚假方程。
    B = 0d0
    fromType = acO%fromType
    toType = acO%toType
    fromNum = acO%fromNum
    toNum = acO%toNum
    eleType = ss%eleType
    obsType = acO%obsType
    
    do i = 1, this%nGroup
        ! 钟差组必须最先分流并 cycle：钟差参数不经过轨道，而且下面的"本组卫星是否
        ! 参与本链路"判断会把地面站钟差组误判成参与 —— 例如 MEO1->Beijing 这条
        ! 链路的 fromNum = toNum = 1，地面站钟差组(origin=1)会满足 isFromSat，
        ! 随后把卫星 1 的轨道灵敏度写进地面站钟的列里，数值上完全看不出来。
        if (this%group(i)%label >= 11 .and. this%group(i)%label <= 13) then
            call calB_clock(this, acO, it, i, m, B)
            cycle
        endif

        satIdx = this%group(i)%origin
        ! 判断本组卫星是否参与这条链路。必须同时判断对象类型：卫星编号与地面站
        ! 编号都从 1 开始，只看数值会把"MEO2->Beijing(toNum=1)"误判成卫星 1 是
        ! 接收端，从而给卫星 1 的列塞进一份它不参与的偏导（且取 ss%sat(1) 当接收
        ! 机），污染法方程。
        isFromSat = (fromType == OBJECT_TYPE_SATELLITE .and. satIdx == fromNum)
        isToSat   = (toType == OBJECT_TYPE_SATELLITE .and. satIdx == toNum)
        if (isFromSat .or. isToSat) then
            nVar = this%group(i)%nVar
            
            ! calculate STM
            idx0 = this%group(i)%idx(1)
            x0(1:6) = this%x(idx0:idx0+5, 1)
            dt = timeMinus(t1, this%t)*86400d0/time_unit
            ! 统一从预计算的 STM 序列(见 type_STM)按历元插值取用 —— 不管 label
            ! 是什么都走这一条路。序列由 buildSTMseries 在遍历观测数据之前生成。
            if (this%STMseries(i)%valid) then
                call this%STMseries(i)%get(t1, STM(1:nVar, 1:nVar))
            else
                ! 退路: 序列未生成时现算(慢), 仅供不经 buildSTMseries 的调用
                select case(this%group(i)%label)
                case(1)
                    call calSTM(eleType, x0, dt, .false., .false., STM(1:8, 1:8))
                case(2)
                    call calSTM(eleType, x0, dt, .true., .false., STM(1:8, 1:8))
                case(3)
                    call calSTM(eleType, x0, dt, .false., .true., STM(1:8, 1:8))
                case(4)
                    call calSTM(eleType, x0, dt, .true., .true., STM(1:8, 1:8))
                case(5)
                    call numericalSTM_ECOM(ss%sat(satIdx), dt, STM(1:15, 1:15))
                case default ! error
                    print*, "error @cal_B_Y"
                    print*, "invalid label. label = ", this%group(i)%label
                endselect
            endif

            ! calculate dR/dX
            x1 = ss%sat(satIdx)%getStates(t1)
            call calDRDX(eleType, x1, DRDX)
                
            ! calculate dY/dR
            ! calDYDR(obsType, PVs, PVo, m, DYDR) 给的是对"第一个参数"的偏导。
            ! 本组卫星是发射端时直接用它; 是接收端时(星间链路才有这种情况)把
            ! 收发两者对调即可 —— 观测量只依赖相对几何, ρ = |r_T - r_R| 对 r_R
            ! 的偏导与把 (r_T, r_R) 换成 (r_R, r_T) 后对第一个参数的偏导相同。
            if (isFromSat) then ! 本组卫星是发射端
                ! get position and velocity of the transmitter satellite
                PVs = ss%sat(fromNum)%getPV(t1)
                if (toType == OBJECT_TYPE_STATION) then
                    ! get position and velocity of the recieve station
                    PVo = ss%station(toNum)%getPV(t1)
                elseif (toType == OBJECT_TYPE_SATELLITE) then
                    ! get position and velocity of the recieve satellite
                    PVo = ss%sat(toNum)%getPV(t1)
                endif
                call calDYDR(obsType, PVs, PVo, m, DYDR)

            else ! 本组卫星是接收端: 对调收发
                PVs = ss%sat(toNum)%getPV(t1)      ! 本组卫星(接收端)
                if (fromType == OBJECT_TYPE_STATION) then
                    PVo = ss%station(fromNum)%getPV(t1)  ! 实际发射站
                else
                    PVo = ss%sat(fromNum)%getPV(t1)      ! 实际发射星
                endif
                call calDYDR(obsType, PVs, PVo, m, DYDR)

            endif
            
            ! Bp -> B
            B(1:m, this%group(i)%idx) = matmul(DYDR, matmul(DRDX, STM(1:6, 1:nVar)))
        endif
    enddo

    endsubroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SUBROUTINE: calB_clock
    ! 钟差组(label 11~13)在观测方程里的一行偏导。
    !
    ! 钟差以"直接项"进入观测，不经过轨道，所以这里不用 STM/DRDX/DYDR，也不碰
    ! 任何状态量。测距模型(module_access 的 calObsData1/2)为
    !   rho_obs = rho + (clk_R - clk_T)*c_norm
    ! 于是对某个钟差系数 p_j(k 阶项, k = j-1):
    !   接收端: drho/dp_j = +c_norm*dt^k ;  发射端: -c_norm*dt^k
    ! dt 以秒计, 基准历元是 this%t(组的参数基准), 不是 clk%t0。
    !
    ! 输入:
    !   acO  - 观测链路
    !   it   - 该链路中的历元下标
    !   igrp - 本组在 this%group 中的下标
    ! 输入输出:
    !   B    - 已由 calB 清零的雅可比行, 本过程只填本组的列
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calB_clock(this, acO, it, igrp, m, B)
    implicit none
    type(type_ODor), intent(in) :: this
    type(type_access), intent(in) :: acO
    integer, intent(in) :: it, igrp, m
    double precision, intent(inout) :: B(m, this%nVar)
    integer :: j, nVar, origin, objType
    double precision :: dtSec, col
    logical :: isFrom, isTo
    !-------------------------------------------------------------------

    ! 只有测距观测的模型里含钟差项；其它观测量(角度/距离变率/位置速度)给全 0 列，
    ! 会让法方程奇异，validateVarGroups 已提前拦掉这种情况。
    if (acO%obsType /= OBS_TYPE_RANGING) return
    if (.not. key_clock_error) return

    origin  = this%group(igrp)%origin
    objType = this%group(igrp)%objType
    ! 参与判定必须"对象类型 + 编号"一起比，理由同 calB 里卫星组的判定
    isFrom = (objType == acO%fromType .and. origin == acO%fromNum)
    isTo   = (objType == acO%toType   .and. origin == acO%toNum)
    if (.not. (isFrom .or. isTo)) return

    nVar  = this%group(igrp)%nVar
    dtSec = timeMinus(acO%timeTable(it), this%t)*86400d0     ! [s]
    do j = 1, nVar
        col = SPEED_OF_LIGHT*u_meter*dtSec**(j - 1)
        if (isFrom) col = -col        ! 发射端的钟差在观测模型里是减号
        B(1, this%group(igrp)%idx(j)) = col
    enddo
    endsubroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calSTM(eleType, x0, dt, keyDrag, keySRP, STM)
    implicit none
    integer, intent(in) :: eleType
    double precision, intent(in) :: x0(6), dt
    logical, intent(in) :: keyDrag, keySRP
    double precision, intent(out) :: STM(8, 8)
    
    STM = 0d0
    
    if (eleType == STATE_TYPE_NONSINGULAR_1) then
        call calSTM_nsg1(x0, dt, .true., keyDrag, STM(1:7, 1:7))
    endif
    if (eleType == STATE_TYPE_KEPLERIAN) then
        call calSTM_kpl(x0, dt, .true., keyDrag, STM(1:7, 1:7))
    endif
    if (eleType == STATE_TYPE_NONSINGULAR_2) then
        call calSTM_nsg2(x0, dt, .true., keyDrag, STM(1:7, 1:7))
    endif
    
    if (keySRP) then
        ! ״̬���Թ�ѹ��Ч���ʱȵĵ��� ��ʱû�б�д
        print*, "unfinished part!"
        STM(8, 8) = 1d0
    endif
    
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate state transition matrix for the first kind of non-singular elements
    ! input:
    !   nsg1        double precision(6)       nsg1 elements at t0
    !   dt          double precision          t1 - t0, unit: day
    !   keyJ2       logical         consider J2 or not
    !   keyDrag     logical         consider atmosphere or not
    ! output:
    !   STM         double precision(7, 7)    state transition matrix
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calSTM_nsg1(nsg1, dt, keyJ2, keyDrag, STM)
    use module_atmosphere
    use module_const
    implicit none
    double precision, intent(in) :: nsg1(6)
    double precision dt
    logical, intent(in) :: keyJ2, keyDrag
    double precision, intent(out) :: STM(7, 7)
    double precision phi0(7, 7), phi1(7, 7)
    double precision a, e, i, p, n, ou, xi, eta, lambda
    integer ii, k
    double precision a61
    double precision J2, omega_up1, omega_low1, M1, lambda1
    double precision b17, b67
    double precision b31, b32, b33, b34, b35
    double precision b41, b42, b43, b44, b45
    double precision b51, b52, b53, b54, b55
    double precision b61, b62, b63, b64, b65
    double precision z, x_k, I0_z
    
    !--------------------------------------------------------------------------
    STM = 0d0
    phi0 = 0d0
    phi1 = 0d0
    
    ! if dt = 0
    if (dabs(dt) < 1d-8) then
        do ii = 1, 7
            STM(ii, ii) = 1d0
        end do
        return
    endif
    
    
    a = nsg1(1)
    e = dsqrt(nsg1(4)**2 + nsg1(5)**2)
    i = nsg1(2)
    ou = nsg1(3) ! uppercase omega
    xi = nsg1(4)
    eta = nsg1(5)
    lambda = nsg1(6)
    p = a*(1d0 - e**2)
    n = a**(-1.5d0)
    
    ! phi0
    do ii = 1, 7
        phi0(ii, ii) = 1d0
    end do
    a61 = -1.5d0*n/a*dt
    phi0(6, 1) = a61
    
    ! phi1
    ! J2 part
    if (keyJ2) then
        J2 = 1.082629989d-3
        omega_up1 = -1.5d0*J2/p**2*n*dcos(i)
        omega_low1 = 1.5d0*J2/p**2*n*(2d0 - 2.5d0*(dsin(i))**2)
        M1 = 1.5d0*J2/p**2*n*(1d0 - 1.5d0*(dsin(i))**2)*dsqrt(1d0 - e**2)
        lambda1 = omega_low1 + M1
        
        b31 = -3.5d0/a*omega_up1*dt
        b32 = -dtan(i)*omega_up1*dt
        b33 = 0d0
        b34 = 4d0*xi/(1d0 - e**2)*omega_up1*dt
        b35 = 4d0*eta/(1d0 - e**2)*omega_up1*dt
        
        b41 = -eta*(-3.5d0/a*omega_low1*dt)
        b42 = -eta*(-5d0*dsin(2d0*i)/(4d0 - 5d0*dsin(i)**2)*omega_low1*dt)
        b43 = 0d0
        b44 = -eta*(4d0*xi/(1d0 - e**2)*omega_low1*dt)
        b45 = -(1d0 + 4d0/(1d0 - e**2)*eta**2)*omega_low1*dt
        
        b51 = xi*(-3.5d0/a*omega_low1*dt)
        b52 = xi*(-5d0*dsin(2d0*i)/(4d0 - 5d0*dsin(i)**2)*omega_low1*dt)
        b53 = 0d0
        b54 = (1d0 + 4d0/(1d0 - e**2)*xi**2)*omega_low1*dt
        b55 = xi*4d0*eta/(1d0 - e**2)*omega_low1*dt
        
        b61 = -3.5d0/a*lambda1*dt
        b62 = -(3d0/(2d0 - 3d0*dsin(i)**2)*M1 + 5d0/(4d0 - 5d0*dsin(i)**2)*omega_low1)*dsin(2d0*i)*dt
        b63 = 0d0
        b64 = xi/(1d0 - e**2)*(3d0*M1 + 4d0*omega_low1)*dt
        b65 = eta/(1d0 - e**2)*(3d0*M1 + 4d0*omega_low1)*dt
    else
        b31 = 0d0
        b32 = 0d0
        b33 = 0d0
        b34 = 0d0
        b35 = 0d0
        b41 = 0d0
        b42 = 0d0
        b43 = 0d0
        b44 = 0d0
        b45 = 0d0
        b51 = 0d0
        b52 = 0d0
        b53 = 0d0
        b54 = 0d0
        b55 = 0d0
        b61 = 0d0
        b62 = 0d0
        b63 = 0d0
        b64 = 0d0
        b65 = 0d0
    endif
    ! atmosphere part
    if (keyDrag) then
    
        z = a*e/H_p0  
        !calculate the I0(z)
        k = 0
        x_k = 1d0
        I0_z = 0d0
        do while(.true.)
            I0_z = I0_z + x_k
            x_k = x_k*(z/2d0)**2/(k + 1d0)**2
            k = k + 1
            if (x_k/I0_z < 1d-13) then
                exit
            end if
        end do
        b17 = - a**2*I0_z*dexp(-(a - r_p0)/H_p0)*n*dt
        b67 = 1.5d0*a*I0_z*dexp(-(a - r_p0)/H_p0)*0.5d0*(n*dt)**2
    else
        b17 = 0d0
        b67 = 0d0
    endif
    
    phi1(1, 7) = b17
    phi1(6, 7) = b67
    phi1(3, 1:5) = (/b31, b32, b33, b34, b35/)
    phi1(4, 1:5) = (/b41, b42, b43, b44, b45/)
    phi1(5, 1:5) = (/b51, b52, b53, b54, b55/)
    phi1(6, 1:5) = (/b61, b62, b63, b64, b65/) 
    
    STM = phi0 + phi1
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calSTM_nsg2(nsg2, dt, keyJ2, keyDrag, STM)
    use module_const
    implicit none
    double precision, intent(in) :: nsg2(6)
    double precision dt
    logical, intent(in) :: keyJ2, keyDrag
    double precision, intent(out) :: STM(7, 7)
    double precision phi0(7, 7), phi1(7, 7)
    double precision a, e, i, p, n, xi, eta, lambda, h, k
    integer ii, kk
    double precision a61
    double precision J2, omega_up1, omega_low1, M1, omega_lt1, lambda1
    double precision b17, b67
    double precision b21, b22, b23, b24, b25
    double precision b31, b32, b33, b34, b35
    double precision b41, b42, b43, b44, b45
    double precision b51, b52, b53, b54, b55
    double precision b61, b62, b63, b64, b65
    double precision r_p0, H_p0, z, x_k, I0_z
    double precision i2, cosi, sini
    !--------------------------------------------------------------------------
    
    STM = 0d0
    phi0 = 0d0
    phi1 = 0d0
    
    ! if dt = 0
    if (dabs(dt) < 1d-8) then
        do ii = 1, 7
            STM(ii, ii) = 1d0
        end do
        return
    endif
    
    a = nsg2(1)
    xi = nsg2(2)
    eta = nsg2(3)
    h = nsg2(4)
    k = nsg2(5)
    lambda = nsg2(6)
    e = dsqrt(nsg2(2)**2 + nsg2(3)**2)
    ! get cos(i) and sin(i)
    i2 = h**2 + k**2 ! i2 is (sin(i/2))^2
    cosi = 1d0 - 2d0*i2 ! cosi is cos(i)
    sini = dsqrt(1d0 - cosi**2) ! sini is sin(i)
    p = a*(1d0 - e**2d0)
    n = a**(-1.5d0)
        
    ! phi0
    do ii = 1, 7
        phi0(ii, ii) = 1d0
    end do
    a61 = -1.5d0*n/a*dt
    phi0(6, 1) = a61
    
    ! phi1
    ! J2 part
    if (keyJ2) then
        J2 = 1.082629989d-3
        omega_up1 = -1.5d0*J2/p**2*n*dcos(i)
        omega_low1 = 1.5d0*J2/p**2*n*(2d0 - 2.5d0*(dsin(i))**2)
        omega_lt1 = omega_low1 + omega_up1
        M1 = 1.5d0*J2/p**2*n*(1d0 - 1.5d0*(dsin(i))**2)*dsqrt(1d0 - e**2)
        lambda1 = M1 + omega_lt1
        
        b21 = -eta*(-3.5d0/a*omega_lt1*dt)
        b22 = -eta*(-4d0*xi/(1d0 - e**2)*omega_lt1*dt)
        b23 = -(1d0 + 4d0/(1d0 - e**2)*eta**2)*omega_lt1*dt
        b24 = -eta*(-2.5d0*8d0*h*(1d0 - 2d0*(i2))+4d0*h)*(3d0*J2/2d0/p**2)*n*dt
        b25 = -eta*(-2.5d0*8d0*k*(1d0 - 2d0*(i2))+4d0*k)*(3d0*J2/2d0/p**2)*n*dt
        
        b31 = -xi*3.5d0/a*omega_lt1*dt
        b32 = (1d0 + 4d0/(1d0 - e**2)*xi**2)*omega_lt1*dt
        b33 = xi*(-4d0*eta/(1d0 - e**2)*omega_lt1*dt)
        b34 = xi*(-2.5d0*8d0*h*(1d0 - 2d0*(i2))+4d0*h)*(3d0*J2/2d0/p**2)*n*dt
        b35 = xi*(-2.5d0*8d0*k*(1d0 - 2d0*(i2))+4d0*k)*(3d0*J2/2d0/p**2)*n*dt
        
        b41 = -k*(-3.5d0/a*omega_up1*dt)
        b42 = -k*(4d0*xi/(1d0 - e**2)*omega_up1*dt)
        b43 = -k*(4d0*eta/(1d0 - e**2)*omega_up1*dt)
        b44 = -k*(4d0*h*(3d0*J2/2d0/p**2)*n*dt)
        b45 = -(omega_up1+4d0*k**2*(3d0*J2/2d0/p**2)*n)*dt
        
        b51 = h*(-3.5d0/a*omega_up1*dt)
        b52 = h*(4d0*xi/(1d0-e**2)*omega_up1*dt)
        b53 = h*(4d0*eta/(1d0-e**2)*omega_up1*dt)
        b54 = (omega_up1+4d0*h**2*(3d0*J2/2d0/p**2)*n)*dt
        b55 = h*(4d0*k*(3d0*J2/2d0/p**2)*n*dt)
        
        b61 = -3.5d0/a*lambda1*dt
        b62 = xi/(1d0-e**2)*(3d0*M1+4d0*omega_lt1)*dt
        b63 = eta/(1d0-e**2)*(3d0*M1+4d0*omega_lt1)*dt
        b64 = ((-1.5d0*dsqrt(1d0-e**2)-2.5d0)*8d0*h*(1d0-2d0*(i2))+4d0*h)*(3d0*J2/2d0/p**2)*n*dt
        b65 = ((-1.5d0*dsqrt(1d0-e**2)-2.5d0)*8d0*k*(1d0-2d0*(i2))+4d0*k)*(3d0*J2/2d0/p**2)*n*dt
    else
        b31 = 0d0
        b32 = 0d0
        b33 = 0d0
        b34 = 0d0
        b35 = 0d0
        b41 = 0d0
        b42 = 0d0
        b43 = 0d0
        b44 = 0d0
        b45 = 0d0
        b51 = 0d0
        b52 = 0d0
        b53 = 0d0
        b54 = 0d0
        b55 = 0d0
        b61 = 0d0
        b62 = 0d0
        b63 = 0d0
        b64 = 0d0
        b65 = 0d0
    endif
    ! atmosphere part
    if (keyDrag) then
    
        z = a*e/H_p0  
        !calculate the I0(z)
        kk = 0
        x_k = 1d0
        I0_z = 0d0
        do while(.true.)
            I0_z = I0_z + x_k
            x_k = x_k*(z/2d0)**2/(dble(kk) + 1d0)**2
            kk = kk + 1
            if (x_k/I0_z < 1d-13) then
                exit
            end if
        end do
        b17 = - a**2*I0_z*dexp(-(a - r_p0)/H_p0)*n*dt
        b67 = 1.5d0*a*I0_z*dexp(-(a - r_p0)/H_p0)*0.5d0*(n*dt)**2
    else
        b17 = 0d0
        b67 = 0d0
    endif
    
    phi1(1, 7) = b17
    phi1(6, 7) = b67
    phi1(3, 1:5) = (/b31, b32, b33, b34, b35/)
    phi1(4, 1:5) = (/b41, b42, b43, b44, b45/)
    phi1(5, 1:5) = (/b51, b52, b53, b54, b55/)
    phi1(6, 1:5) = (/b61, b62, b63, b64, b65/) 
    
    STM = phi0 + phi1
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate STM for keplerian using analytical solution of orbit
    ! input:
    !   elements    double precision(6)       keplerian elements at t0
    !   dt          double precision          t1 - t0, unit: tu, tu = 806.81099130673260d0 sec
    !   keyJ2       logical         consider J2 or not
    !   keyDrag     logical         consider atmosphere or not
    ! output:
    !   STM         double precision(7, 7)    state transition matrix
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! o_o verified
    subroutine calSTM_kpl(elements, dt, keyJ2, keyDrag, STM)
    use module_const
    use module_atmosphere
    implicit none
    double precision elements(6), STM(7, 7)
    double precision dt
    logical, intent(in) :: keyJ2, keyDrag
    double precision phi0(7, 7), phi1(7, 7)
    double precision a, e, i, p, n
    integer ii, k
    double precision a61
    double precision J2, omega_up1, omega_low1, M1
    double precision b17, b67
    double precision b41, b42, b43, b51, b52, b53, b61, b62, b63
    double precision z, x_k, I0_z
    
    !-----------------------------------------------------------------
    STM = 0d0
    phi0 = 0d0
    phi1 = 0d0
    !-----------------------------------------------------------------
    
    a = elements(1)
    e = elements(2)
    i = elements(3)
    p = a*(1d0 - e**2)
    n = a**(-1.5d0)
    
    !-----------------------------------------------------------------
        
    do ii = 1, 7
        phi0(ii, ii) = 1d0
    end do
    
    ! if dt = 0
    if (dabs(dt) < 1d-8) then
        STM = phi0
        return
    endif
    
    a61 = -1.5d0*n/a*dt
    phi0(6, 1) = a61
    
    !---------------------------------------------------------------
    if (keyJ2) then
        J2 = 1.082629989d-3
        omega_up1 = -1.5d0*J2/p**2*n*dcos(i)
        omega_low1 = 1.5d0*J2/p**2*n*(2d0 - 2.5d0*(dsin(i))**2)
        M1 = 1.5d0*J2/p**2*n*(1d0 - 1.5d0*(dsin(i))**2)*dsqrt(1d0 - e**2)
    
        b41 = -3.5d0/a*omega_up1*dt
        b42 = 4d0*e/(1d0 - e**2)*omega_up1*dt
        b43 = -dtan(i)*omega_up1*dt
        b51 = -3.5d0/a*omega_low1*dt
        b52 = 4d0*e/(1d0 - e**2)*omega_low1*dt
        b53 = -5d0*dsin(2d0*i)/(4d0 - 5d0*(dsin(i))**2)*omega_low1*dt
        b61 = -3.5d0/a*M1*dt
        b62 = 3d0*e/(1d0 - e**2)*M1*dt
        b63 = -3d0*dsin(2d0*i)/(2d0 - 3d0*(dsin(i))**2)*M1*dt
    else
        b41 = 0d0
        b42 = 0d0
        b43 = 0d0
        b51 = 0d0
        b52 = 0d0
        b53 = 0d0
        b61 = 0d0
        b62 = 0d0
        b63 = 0d0
    endif
    
    if (keyDrag) then
    
        !r_p0 = 1.069482963912503D0
        !H_p0 = 9.328742860180019D-03
        z = a*e/H_p0  
        !calculate the I0(z)
        k = 0
        x_k = 1d0
        I0_z = 0d0
        do while(.true.)
            I0_z = I0_z + x_k
            x_k = x_k*(z/2d0)**2/(k + 1d0)**2
            k = k + 1
            if (x_k/I0_z < 1d-13) then
                exit
            end if
        end do
        b17 = - a**2*I0_z*dexp(-(a - r_p0)/H_p0)*n*dt
        b67 = 1.5d0*a*I0_z*dexp(-(a - r_p0)/H_p0)*0.5d0*(n*dt)**2
    else
        b17 = 0d0
        b67 = 0d0
    endif
    
    phi1(1, 7) = b17
    phi1(6, 7) = b67
    phi1(4, 1:3) = (/b41, b42, b43/)
    phi1(5, 1:3) = (/b51, b52, b53/)
    phi1(6, 1:3) = (/b61, b62, b63/)
    
    STM = phi0 + phi1
    
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calSTM_cartesian(x0, t0, dt, STM)
    use module_planets
    implicit none
    double precision, intent(in) :: x0(7), t0, dt
    double precision, intent(out) :: STM(7, 7)
    double precision tau, x, y, z, vx, vy, vz, rv, r
    double precision r_moon(3), PV(6), deltaM, xm, ym, zm
    double precision r_sun(3), deltaS, xs, ys, zs
    double precision A, nu
    !--------------------------------------------------------------------------
    
    ! prepare variables
    x = x0(1)
    y = x0(2)
    z = x0(3)
    vx = x0(4)
    vy = x0(5)
    vz = x0(6)
    A = x0(7)
    rv = x*vx + y*vy + z*vz
    r = norm2(x0(1:3))
    STM = 0d0
    
    ! Calculate moon-related variables
    call planetsPV(t0, nEarth, nMoon, PV)
    r_moon = PV(1:3) * u_km
    deltaM = norm2(x0(1:3) - r_moon(1:3))
    xm = r_moon(1)
    ym = r_moon(2)
    zm = r_moon(3)
    ! Calculate sun-related variables
    call planetsPV(t0, nEarth, nsun, PV)
    r_sun = PV(1:3) * u_km
    deltaS = norm2(x0(1:3) - r_sun(1:3))
    xs = r_sun(1)
    ys = r_sun(2)
    zs = r_sun(3)
    call EarthShadowCoefficient(x0(1:3), r_sun, nu)
    
    STM(1, 1)= dt**2*(-(moonEarthMassRatio*(1/deltaM**3 - (3*(x - xm)*(x - xm))/deltaM**5))/2 - 1/(2*r**3) + &
        (A*nu)/(2*deltaS) + (3*x*x)/(2*r**5) - (A*nu*(x - xs)*(x - xs))/(2*deltaS**3)) + 1
    STM(1, 2)= -dt**2*(-(3*moonEarthMassRatio*(y - ym)*(x - xm))/(2*deltaM**5) - (3*x*y)/(2*r**5) + &
        (A*nu*(y - ys)*(x - xs))/(2*deltaS**3))
    STM(1, 3)= -dt**2*(-(3*moonEarthMassRatio*(z - zm)*(x - xm))/(2*deltaM**5) - (3*x*z)/(2*r**5) + &
        (A*nu*(z - zs)*(x - xs))/(2*deltaS**3))
    STM(1, 4)= dt
    STM(1, 5)= 0
    STM(1, 6)= 0
    STM(1, 7)= (dt**2*nu*(x - xs))/(2*deltaS)
    STM(2, 1)= -dt**2*(-(3*moonEarthMassRatio*(x - xm)*(y - ym))/(2*deltaM**5) - (3*y*x)/(2*r**5) + &
        (A*nu*(x - xs)*(y - ys))/(2*deltaS**3))
    STM(2, 2)= dt**2*(-(moonEarthMassRatio*(1/deltaM**3 - (3*(y - ym)*(y - ym))/deltaM**5))/2 - 1/(2*r**3) + &
        (A*nu)/(2*deltaS) + (3*y*y)/(2*r**5) - (A*nu*(y - ys)*(y - ys))/(2*deltaS**3)) + 1
    STM(2, 3)= -dt**2*(-(3*moonEarthMassRatio*(z - zm)*(y - ym))/(2*deltaM**5) - (3*y*z)/(2*r**5) + &
        (A*nu*(z - zs)*(y - ys))/(2*deltaS**3))
    STM(2, 4)= 0
    STM(2, 5)= dt
    STM(2, 6)= 0
    STM(2, 7)= (dt**2*nu*(y - ys))/(2*deltaS)
    STM(3, 1)= -dt**2*(-(3*moonEarthMassRatio*(x - xm)*(z - zm))/(2*deltaM**5) - (3*z*x)/(2*r**5) + &
        (A*nu*(x - xs)*(z - zs))/(2*deltaS**3))
    STM(3, 2)= -dt**2*(-(3*moonEarthMassRatio*(y - ym)*(z - zm))/(2*deltaM**5) - (3*z*y)/(2*r**5) + &
        (A*nu*(y - ys)*(z - zs))/(2*deltaS**3))
    STM(3, 3)= dt**2*(-(moonEarthMassRatio*(1/deltaM**3 - (3*(z - zm)*(z - zm))/deltaM**5))/2 - 1/(2*r**3) + &
        (A*nu)/(2*deltaS) + (3*z*z)/(2*r**5) - (A*nu*(z - zs)*(z - zs))/(2*deltaS**3)) + 1
    STM(3, 4)= 0
    STM(3, 5)= 0
    STM(3, 6)= dt
    STM(3, 7)= (dt**2*nu*(z - zs))/(2*deltaS)
    
    STM(4, 1)= dt*(-moonEarthMassRatio*(1d0/deltaM**3 - (3d0*(x - xm)*(x - xm))/deltaM**5) - 1d0/r**3 + &
        (A*nu)/deltaS + (3d0*x*x)/r**5 - (A*nu*(x - xs)*(x - xs))/deltaS**3) 

    STM(4, 2)= -dt*(-(3*moonEarthMassRatio*(y - ym)*(x - xm))/deltaM**5 - (3*x*y)/r**5 + &
        (A*nu*(y - ys)*(x - xs))/deltaS**3)
    STM(4, 3)= -dt*(-(3*moonEarthMassRatio*(z - zm)*(x - xm))/deltaM**5 - (3*x*z)/r**5 + &
        (A*nu*(z - zs)*(x - xs))/deltaS**3)
    STM(4, 4)= 1
    STM(4, 5)= 0
    STM(4, 6)= 0
    STM(4, 7)= (dt*nu*(x - xs))/deltaS
    STM(5, 1)= -dt*(-(3*moonEarthMassRatio*(x - xm)*(y - ym))/deltaM**5 - (3*y*x)/r**5 + &
        (A*nu*(x - xs)*(y - ys))/deltaS**3)
    STM(5, 2)= dt*(-moonEarthMassRatio*(1/deltaM**3 - (3*(y - ym)*(y - ym))/deltaM**5) - &
        1/r**3 + (A*nu)/deltaS + (3*y*y)/r**5 - (A*nu*(y - ys)*(y - ys))/deltaS**3)
    STM(5, 3)= -dt*(-(3*moonEarthMassRatio*(z - zm)*(y - ym))/deltaM**5 - (3*y*z)/r**5 + &
        (A*nu*(z - zs)*(y - ys))/deltaS**3)
    STM(5, 4)= 0
    STM(5, 5)= 1
    STM(5, 6)= 0
    STM(5, 7)= (dt*nu*(y - ys))/deltaS
    STM(6, 1)= -dt*(-(3*moonEarthMassRatio*(x - xm)*(z - zm))/deltaM**5 - (3*z*x)/r**5 + &
        (A*nu*(x - xs)*(z - zs))/deltaS**3)
    STM(6, 2)= -dt*(-(3*moonEarthMassRatio*(y - ym)*(z - zm))/deltaM**5 - (3*z*y)/r**5 + &
        (A*nu*(y - ys)*(z - zs))/deltaS**3)
    STM(6, 3)= dt*(-moonEarthMassRatio*(1/deltaM**3 - (3*(z - zm)*(z - zm))/deltaM**5) - &
        1/r**3 + (A*nu)/deltaS + (3*z*z)/r**5 - (A*nu*(z - zs)*(z - zs))/deltaS**3)
    STM(6, 4)= 0
    STM(6, 5)= 0
    STM(6, 6)= 1
    STM(6, 7)= (dt*nu*(z - zs))/deltaS
    STM(7, 1)= 0
    STM(7, 2)= 0
    STM(7, 3)= 0
    STM(7, 4)= 0
    STM(7, 5)= 0
    STM(7, 6)= 0
    STM(7, 7)= 1
    
    STM(1:3, 7) = STM(1:3, 7) + dt**3/6d0*nu*(x0(4:6)/deltaS - dot_product(x0(4:6), x0(1:3) - &
        r_sun(1:3))/deltaS**2*(x0(1:3) - r_sun(1:3)))
    STM(4:6, 7) = STM(4:6, 7) + dt**2/2d0*nu*(x0(4:6)/deltaS - dot_product(x0(4:6), x0(1:3) - &
        r_sun(1:3))/deltaS**2*(x0(1:3) - r_sun(1:3)))
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine EarthShadowCoefficient(r_sat, r_sun, nu)
    implicit none
    double precision nu
    double precision r_sat(3), r_sun(3), r_sun_unit(3), r_sat_unit(3)
    double precision sin_phi, cos_phi
    !---------------------------------------------------------------
    r_sun_unit = r_sun/norm2(r_sun)
    r_sat_unit = r_sat/norm2(r_sat)
    call dot(3, r_sun_unit, r_sat_unit, cos_phi)
    sin_phi = dsqrt(1d0 - cos_phi**2)
    if ((cos_phi < 0) .and. (sin_phi < 1d0/norm2(r_sat))) then
        nu = 0d0
    else
        nu = 1d0
    end if
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calW(ac, npe, W)
    type(type_access), intent(in) :: ac
    integer, intent(in) :: npe
    double precision, intent(out) :: W(npe, npe)
    integer i
    W = 0d0
    do i = 1, npe
        W(i, i) = 1d0
    enddo
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate matrix DR/DX, R is position and velocity, X is elements
    ! type of the elements can be keplerian, NSG1, NSG2 et al.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calDRDX(eleType, x, DRDX)
    implicit none
    integer, intent(in) :: eleType
    double precision, intent(in) :: x(6)
    double precision, intent(out) :: DRDX(6, 6)
    integer i
    select case (eleType)
    case(STATE_TYPE_NONSINGULAR_1)
        call calDRDX_NSG1(x, DRDX)
    case(STATE_TYPE_KEPLERIAN)
        call calDRDX_kpl(x, DRDX)
    case(STATE_TYPE_NONSINGULAR_2)
        call calDRDX_nsg2(x, DRDX)
    case(STATE_TYPE_CARTESIAN)
        DRDX = 0d0
        do i = 1, 6
            DRDX(i, i) = 1d0
        enddo
    endselect
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! verified
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calDRDX_NSG1(nsg1, DRDX)
    implicit none
    double precision, intent(in) :: nsg1(6)
    double precision, intent(out) :: DRDX(6, 6)
    double precision a, i, ou, xi, eta, lambda
    double precision e, p, n, r, u, ut, sinu, cosu
    double precision AA, B, C, D, Ap, Bp, Cp, Dp
    double precision kpl(6), cts(6), rs(3), vs(3), x, y, z, x_dot, y_dot, z_dot
    double precision omega_up, R_hat(3), omega_up_prime_vector(3)
    double precision I_prime_vector(3), omega_up_vector(3), I_vector(3)
    !-------------------------------------------------------------------
    a = nsg1(1)
    i = nsg1(2)
    ou = nsg1(3)
    xi = nsg1(4)
    eta = nsg1(5)
    lambda = nsg1(6)
    
    e = dsqrt(xi**2 + eta**2)
    p = a*(1d0 - e**2)
    n = a**(-1.5d0)
    
    call nsg12rut(nsg1, r, u, ut)
    
    sinu = dsin(u)
    cosu = dcos(u)
    
    AA = a/p*(-(cosu+xi) - r/p*(sinu+eta)*(xi*sinu-eta*cosu)) ! a is used by semi-major axis, use AA here
    B = a*r/dsqrt(p)*(sinu+a/r*dsqrt(1d0-e**2)/(1d0+dsqrt(1d0-e**2))*eta+r/p*(sinu+eta))
    C = -a/p*((sinu+eta)-r/p*(cosu+xi)*(xi*sinu-eta*cosu))
    D = -a*r/dsqrt(p)*(cosu+a/r*dsqrt(1-e**2)/(1d0+dsqrt(1d0-e**2))*xi+r/p*(cosu+xi))
    
    Ap = dsqrt(1d0/p)*(1d0/p*(sinu+eta)-1d0/r*(a/p*(sinu+eta)+a/r*sinu+ &
        (a/r)**2*dsqrt(1d0-e**2)/(1d0+dsqrt(1d0-e**2))*eta))
    Bp = r/p*cosu+a/p*xi
    Cp = -dsqrt(1d0/p)*(1d0/p*(cosu+xi)-1d0/r*(a/p*(cosu+xi)+a/r*cosu+ &
        (a/r)**2*dsqrt(1d0-e**2)/(1d0+dsqrt(1d0-e**2))*xi))
    Dp = r/p*sinu+a/p*eta
    
    call nsg12cts(nsg1, cts)
    rs = cts(1:3)
    vs = cts(4:6)
    
    x = rs(1)
    y = rs(2)
    z = rs(3)
    
    x_dot = vs(1)
    y_dot = vs(2)
    z_dot = vs(3)
    
    omega_up = ou
    R_hat = (/dsin(i)*dsin(omega_up), -dsin(i)*dcos(omega_up), dcos(i)/)
    omega_up_prime_vector = (/-y_dot, x_dot, 0d0/)
    
    I_prime_vector(1) = z_dot*dsin(omega_up)
    I_prime_vector(2) = -z_dot*dcos(omega_up)
    I_prime_vector(3) = -x_dot*dsin(omega_up) + y_dot*dcos(omega_up)
    
    omega_up_vector = (/-y, x, 0d0/)
    
    I_vector(1) = z*dsin(omega_up)
    I_vector(2) = -z*dcos(omega_up)
    I_vector(3) = -x*dsin(omega_up) + y*dcos(omega_up)
    
    DRDX(1:3, 1) = 1d0/a*rs
    DRDX(1:3, 2) = I_vector
    DRDX(1:3, 3) = omega_up_vector
    DRDX(1:3, 4) = AA*rs + B*vs
    DRDX(1:3, 5) = C*rs + D*vs
    DRDX(1:3, 6) = 1d0/n*vs
    
    DRDX(4:6, 1) = -1d0/2d0/a*vs
    DRDX(4:6, 2) = I_prime_vector
    DRDX(4:6, 3) = omega_up_prime_vector
    DRDX(4:6, 4) = Ap*rs + Bp*vs
    DRDX(4:6, 5) = Cp*rs + Dp*vs
    DRDX(4:6, 6) = -1d0/n*rs/r**3
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! o_o verified
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calDRDX_kpl(elements, AM)
    implicit none
    double precision elements(6), AM(6, 6)
    double precision rs(3), vs(3), cts(6)
    double precision a, e, i, omega_up, omega_low, M
    double precision x, y, z, x_dot, y_dot, z_dot
    double precision R_hat(3), omega_low_prime_vector(3), omega_up_prime_vector(3)
    double precision I_prime_vector(3)
    double precision r, p, n, EE
    double precision H_prime, K_prime, omega_low_vector(3), omega_up_vector(3)
    double precision I_vector(3), H, K
    !-----------------------------------------------------------------------
    
    call kpl2cts(elements, cts)
    rs = cts(1:3)
    vs = cts(4:6)
    
    a = elements(1)
    e = elements(2)
    i = elements(3)
    omega_up = elements(4)
    omega_low = elements(5)
    M = elements(6)
    
    x = rs(1)
    y = rs(2)
    z = rs(3)
    
    x_dot = vs(1)
    y_dot = vs(2)
    z_dot = vs(3)
    
    !-----------------------------------------------------
    
    R_hat = (/dsin(i)*dsin(omega_up), -dsin(i)*dcos(omega_up), dcos(i)/)
    call cross_3d(R_hat, vs, omega_low_prime_vector)
    omega_up_prime_vector = (/-y_dot, x_dot, 0d0/)
    
    I_prime_vector(1) = z_dot*dsin(omega_up)
    I_prime_vector(2) = -z_dot*dcos(omega_up)
    I_prime_vector(3) = -x_dot*dsin(omega_up) + y_dot*dcos(omega_up)
    
    call solve_kepler_equation(e, M, EE)
    
    r = norm2(rs)
    p = a*(1d0 - e**2)
    n = a**(-1.5d0)
    H_prime = dsqrt(a)/r/p*dsin(EE)*(1d0 - a/r*(1d0 + p/r))
    K_prime = a/p*dcos(EE)
    call cross_3d(R_hat, rs, omega_low_vector)
    omega_up_vector = (/-y, x, 0d0/)
    
    I_vector(1) = z*dsin(omega_up)
    I_vector(2) = -z*dcos(omega_up)
    I_vector(3) = -x*dsin(omega_up) + y*dcos(omega_up)
    
    H = -a/p*(dcos(EE) + e)
    K = dsin(EE)/n*(1d0 + r/p)
    
    !------------------------------------------------------------------------------
    
    AM(1:3, 1) = 1d0/a*rs
    AM(1:3, 2) = H*rs + K*vs
    AM(1:3, 3) = I_vector
    AM(1:3, 4) = omega_up_vector
    AM(1:3, 5) = omega_low_vector
    AM(1:3, 6) = 1d0/n*vs
    
    AM(4:6, 1) = -1d0/2d0/a*vs
    AM(4:6, 2) = H_prime*rs + K_prime*vs
    AM(4:6, 3) = I_prime_vector
    AM(4:6, 4) = omega_up_prime_vector
    AM(4:6, 5) = omega_low_prime_vector
    AM(4:6, 6) = -1d0/n/r**3*rs
    
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calDRDX_NSG2(nsg2, DRDX)
    implicit none
    double precision, intent(in) :: nsg2(6)
    double precision, intent(out) :: DRDX(6, 6)
    double precision a, xi, eta, h, k, lambda
    double precision e, p, n, r, u, ut, sinu, cosu
    double precision AA, B, C, D, Ap, Bp, Cp, Dp
    double precision cosi2, P_h_hat(3), Q_h_hat(3), P_k_hat(3), Q_k_hat(3)
    double precision kpl(6), rs(3), vs(3), cts(6)
    !--------------------------------------------------------------------------
    a = nsg2(1)
    xi = nsg2(2)
    eta = nsg2(3)
    h = nsg2(4)
    k = nsg2(5)
    lambda = nsg2(6)
    
    e = dsqrt(xi**2 + eta**2)
    p = a*(1d0 - e**2)
    n = a**(-1.5d0)
    
    call nsg22rut(nsg2, r, u, ut)
    
    sinu = dsin(u)
    cosu = dcos(u)
    
    AA = a/p*(-(cosu+xi) - r/p*(sinu+eta)*(xi*sinu-eta*cosu)) ! a is used by semi-major axis, use AA here
    B = a*r/dsqrt(p)*(sinu+a/r*dsqrt(1d0-e**2)/(1d0+dsqrt(1d0-e**2))*eta+r/p*(sinu+eta))
    C = -a/p*((sinu+eta)-r/p*(cosu+xi)*(xi*sinu-eta*cosu))
    D = -a*r/dsqrt(p)*(cosu + a/r*dsqrt(1-e**2)/(1d0+dsqrt(1d0-e**2))*xi+r/p*(cosu+xi))
    
    Ap = dsqrt(1d0/p)*(1d0/p*(sinu+eta)-1d0/r*(a/p*(sinu+eta)+a/r*sinu+ &
        (a/r)**2*dsqrt(1d0-e**2)/(1d0+dsqrt(1d0-e**2))*eta))
    Bp = r/p*cosu+a/p*xi
    Cp = -dsqrt(1d0/p)*(1d0/p*(cosu+xi)-1d0/r*(a/p*(cosu+xi)+a/r*cosu+ &
        (a/r)**2*dsqrt(1d0-e**2)/(1d0+dsqrt(1d0-e**2))*xi))
    Dp = r/p*sinu+a/p*eta
    
    cosi2 = dsqrt(1d0 - h**2 - k**2) ! cosi2 is cos(i/2)
    
    call nsg22kpl(nsg2, kpl)
    call kpl2cts(kpl, cts)
    rs = cts(1:3)
    vs = cts(4:6)
    
    P_h_hat = (/0d0, 2d0*k, 2d0*h*k/cosi2/)
    Q_h_hat = (/2d0*k, -4d0*h, 2d0*(1d0 - 2d0*h**2 - k**2)/cosi2/)
    P_k_hat = (/-4d0*k, 2d0*h, -2d0*(1d0 - h**2 - 2d0*k**2)/cosi2/)
    Q_k_hat = (/2d0*h, 0d0, -2d0*h*k/cosi2/)
    
    DRDX(1:3, 1) = 1d0/a*rs
    DRDX(1:3, 2) = AA*rs + B*vs
    DRDX(1:3, 3) = C*rs + D*vs
    DRDX(1:3, 4) = r*cosu*P_h_hat + r*sinu*Q_h_hat
    DRDX(1:3, 5) = r*cosu*P_k_hat + r*sinu*Q_k_hat
    DRDX(1:3, 6) = 1d0/n*vs
    
    DRDX(4:6, 1) = -1d0/2d0/a*vs
    DRDX(4:6, 2) = Ap*rs + Bp*vs
    DRDX(4:6, 3) = Cp*rs + Dp*vs
    DRDX(4:6, 4) = 1d0/dsqrt(p)*(-(sinu+eta)*P_h_hat + (cosu+xi)*Q_h_hat)
    DRDX(4:6, 5) = 1d0/dsqrt(p)*(-(sinu+eta)*P_k_hat + (cosu+xi)*Q_k_hat)
    DRDX(4:6, 6) = -1d0/n*rs/r**3
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calDXDR_nsg2(rv, dXdR)
    implicit none
    double precision, intent(in)  :: rv(6)      ! [x,y,z, vx,vy,vz]
    double precision, intent(out) :: dXdR(6,6)  ! d(nsg2)/d(r,v)

    ! Locals
    double precision :: r(3), v(3), mu
    double precision :: r_norm, v_norm, sigma
    double precision :: a, e, i, Omega_node, omega_peri, varpi
    double precision :: e_vec(3), P(3), h_vec(3), h, W(3), Q(3)
    double precision :: sin_i, sin_i_half, cos_i_half, cos_Omega, sin_Omega
    double precision :: C, S, E_angle, sinE, cosE, M

    ! Gradients
    double precision :: da_dr(3), da_dv(3)
    double precision :: de_vec_dr(3,3), de_vec_dv(3,3)
    double precision :: de_dr(3), de_dv(3)
    double precision :: identity(3,3), cross_v(3,3), cross_r(3,3), proj(3,3)
    double precision :: dW_dr(3,3), dW_dv(3,3)
    double precision :: di_dr(3), di_dv(3)
    double precision :: dOmega_node_dr(3), dOmega_node_dv(3)
    double precision :: dPx_dr(3), dPy_dr(3), dPz_dr(3)
    double precision :: dPx_dv(3), dPy_dv(3), dPz_dv(3)
    double precision :: dQz_dr(3), dQz_dv(3)
    double precision :: domega_peri_dr(3), domega_peri_dv(3)
    double precision :: dr_dr(3), dsigma_dr(3), dsigma_dv(3)
    double precision :: dC_dr(3), dC_dv(3), dS_dr(3), dS_dv(3)
    double precision :: dM_dr(3), dM_dv(3)

    double precision :: dkpl_dX(6,6), dnsg2_dkpl(6,6)
    integer :: j, k

    mu = 1.0d0

    r = rv(1:3)
    v = rv(4:6)

    r_norm = norm2(r)
    v_norm = norm2(v)
    sigma  = dot_product(r, v)

    ! ---- Semi-major axis a ----
    a = 1.0d0 / (2.0d0 / r_norm - v_norm**2 / mu)
    da_dr = 2.0d0 * a**2 / r_norm**3 * r
    da_dv = 2.0d0 * a**2 / mu * v

    ! ---- Eccentricity vector and eccentricity ----
    e_vec = (v_norm**2 - mu / r_norm) * r - sigma * v
    e     = norm2(e_vec)
    P     = e_vec / e

    ! ---- Angular momentum unit vector ----
    h_vec(1) = r(2)*v(3) - r(3)*v(2)
    h_vec(2) = r(3)*v(1) - r(1)*v(3)
    h_vec(3) = r(1)*v(2) - r(2)*v(1)
    h = norm2(h_vec)
    W = h_vec / h

    i           = acos(W(3))
    sin_i       = sin(i)
    Omega_node  = atan2(W(1), -W(2))   ! 升交点赤经

    ! ---- Q = W x P (transversal unit vector) ----
    Q(1) = W(2)*P(3) - W(3)*P(2)
    Q(2) = W(3)*P(1) - W(1)*P(3)
    Q(3) = W(1)*P(2) - W(2)*P(1)

    omega_peri = atan2(P(3), Q(3))     ! 近地点幅角
    varpi      = Omega_node + omega_peri

    ! ---- Eccentric and mean anomaly ----
    C      = 1.0d0 - r_norm / a
    S      = sigma / sqrt(mu * a)
    E_angle = atan2(S, C)
    sinE    = sin(E_angle)
    cosE    = cos(E_angle)
    M       = E_angle - e * sinE

    ! ============================================================
    !  Compute d(kpl)/d(r,v)   (6 x 6)
    ! ============================================================
    identity = 0.0d0
    do k = 1, 3
        identity(k,k) = 1.0d0
    end do

    ! --- de_vec / dr ,  de_vec / dv ---
    do j = 1, 3
        do k = 1, 3
            de_vec_dr(j,k) = (v_norm**2 - 1.0d0/r_norm) * identity(j,k) &
                           + (1.0d0 / r_norm**3) * r(j) * r(k)      &
                           - v(j) * v(k)

            de_vec_dv(j,k) = 2.0d0 * r(j) * v(k) - v(j) * r(k) &
                           - sigma * identity(j,k)
        end do
    end do

    ! de / d(r,v)
    do k = 1, 3
        de_dr(k) = sum(P * de_vec_dr(:,k))
        de_dv(k) = sum(P * de_vec_dv(:,k))
    end do

    ! --- Cross‑product matrices ---
    cross_v(1,1)=0.0d0; cross_v(1,2)=-v(3); cross_v(1,3)=v(2)
    cross_v(2,1)=v(3);  cross_v(2,2)=0.0d0; cross_v(2,3)=-v(1)
    cross_v(3,1)=-v(2); cross_v(3,2)=v(1);  cross_v(3,3)=0.0d0

    cross_r(1,1)=0.0d0; cross_r(1,2)=-r(3); cross_r(1,3)=r(2)
    cross_r(2,1)=r(3);  cross_r(2,2)=0.0d0; cross_r(2,3)=-r(1)
    cross_r(3,1)=-r(2); cross_r(3,2)=r(1);  cross_r(3,3)=0.0d0

    ! proj = I - W W^T   (替换原 A 矩阵，避免与半长轴 a 冲突)
    do j = 1, 3
        do k = 1, 3
            proj(j,k) = identity(j,k) - W(j)*W(k)
        end do
    end do

    ! dW/dr = - (1/h) * proj * [v×] ,   dW/dv = (1/h) * proj * [r×]
    dW_dr = - (1.0d0 / h) * matmul(proj, cross_v)
    dW_dv =   (1.0d0 / h) * matmul(proj, cross_r)

    ! --- Inclination i ---
    di_dr = - (1.0d0 / sin_i) * dW_dr(3,:)
    di_dv = - (1.0d0 / sin_i) * dW_dv(3,:)

    ! --- Right ascension of ascending node Omega_node ---
    dOmega_node_dr = ( -W(2)*dW_dr(1,:) + W(1)*dW_dr(2,:) ) / (sin_i**2)
    dOmega_node_dv = ( -W(2)*dW_dv(1,:) + W(1)*dW_dv(2,:) ) / (sin_i**2)

    ! --- Derivatives of P components (for omega) ---
    do k = 1, 3
        dPx_dr(k) = (1.0d0/e) * de_vec_dr(1,k) - (e_vec(1)/e**2) * de_dr(k)
        dPy_dr(k) = (1.0d0/e) * de_vec_dr(2,k) - (e_vec(2)/e**2) * de_dr(k)
        dPz_dr(k) = (1.0d0/e) * de_vec_dr(3,k) - (e_vec(3)/e**2) * de_dr(k)

        dPx_dv(k) = (1.0d0/e) * de_vec_dv(1,k) - (e_vec(1)/e**2) * de_dv(k)
        dPy_dv(k) = (1.0d0/e) * de_vec_dv(2,k) - (e_vec(2)/e**2) * de_dv(k)
        dPz_dv(k) = (1.0d0/e) * de_vec_dv(3,k) - (e_vec(3)/e**2) * de_dv(k)
    end do

    ! dQz / d(r,v) ,  Qz = Wx*Py - Wy*Px
    do k = 1, 3
        dQz_dr(k) = P(2) * dW_dr(1,k) + W(1) * dPy_dr(k) &
                  - P(1) * dW_dr(2,k) - W(2) * dPx_dr(k)
        dQz_dv(k) = P(2) * dW_dv(1,k) + W(1) * dPy_dv(k) &
                  - P(1) * dW_dv(2,k) - W(2) * dPx_dv(k)
    end do

    ! --- Argument of periapsis omega_peri ---
    domega_peri_dr = ( cos(omega_peri) * dPz_dr - sin(omega_peri) * dQz_dr ) / sin_i
    domega_peri_dv = ( cos(omega_peri) * dPz_dv - sin(omega_peri) * dQz_dv ) / sin_i

    ! --- Mean anomaly M ---
    dr_dr     = r / r_norm
    dsigma_dr = v
    dsigma_dv = r

    do k = 1, 3
        dC_dr(k) = - (1.0d0 / a) * dr_dr(k)   + (r_norm / a**2) * da_dr(k)
        dC_dv(k) =   (r_norm / a**2) * da_dv(k)

        dS_dr(k) = (1.0d0 / sqrt(a)) * v(k) &
                 - (sigma / (2.0d0 * a**1.5d0)) * da_dr(k)
        dS_dv(k) = (1.0d0 / sqrt(a)) * r(k) &
                 - (sigma / (2.0d0 * a**1.5d0)) * da_dv(k)
    end do

    do k = 1, 3
        dM_dr(k) = (- sinE / e) * dC_dr(k) + ((cosE - e) / e) * dS_dr(k)
        dM_dv(k) = (- sinE / e) * dC_dv(k) + ((cosE - e) / e) * dS_dv(k)
    end do

    ! Assemble d(kpl)/d(r,v)  (rows: a, e, i, Omega_node, omega_peri, M)
    dkpl_dX(1,1:3) = da_dr
    dkpl_dX(1,4:6) = da_dv
    dkpl_dX(2,1:3) = de_dr
    dkpl_dX(2,4:6) = de_dv
    dkpl_dX(3,1:3) = di_dr
    dkpl_dX(3,4:6) = di_dv
    dkpl_dX(4,1:3) = dOmega_node_dr
    dkpl_dX(4,4:6) = dOmega_node_dv
    dkpl_dX(5,1:3) = domega_peri_dr
    dkpl_dX(5,4:6) = domega_peri_dv
    dkpl_dX(6,1:3) = dM_dr
    dkpl_dX(6,4:6) = dM_dv

    ! ============================================================
    !  Matrix  d(nsg2) / d(kpl)   (6 x 6)
    ! ============================================================
    sin_i_half = sin(i / 2.0d0)
    cos_i_half = cos(i / 2.0d0)
    cos_Omega  = cos(Omega_node)
    sin_Omega  = sin(Omega_node)

    dnsg2_dkpl(1,:) = (/ 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    dnsg2_dkpl(2,:) = (/ 0.0d0, cos(varpi), 0.0d0, -e*sin(varpi), &
                         -e*sin(varpi), 0.0d0 /)
    dnsg2_dkpl(3,:) = (/ 0.0d0, sin(varpi), 0.0d0,  e*cos(varpi), &
                          e*cos(varpi), 0.0d0 /)
    dnsg2_dkpl(4,:) = (/ 0.0d0, 0.0d0, 0.5d0*cos_i_half*cos_Omega, &
                         -sin_i_half*sin_Omega, 0.0d0, 0.0d0 /)
    dnsg2_dkpl(5,:) = (/ 0.0d0, 0.0d0, 0.5d0*cos_i_half*sin_Omega, &
                          sin_i_half*cos_Omega, 0.0d0, 0.0d0 /)
    dnsg2_dkpl(6,:) = (/ 0.0d0, 0.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0 /)

    ! ============================================================
    !  Final Jacobian: d(nsg2)/d(r,v) = [d(nsg2)/d(kpl)] * [d(kpl)/d(r,v)]
    ! ============================================================
    dXdR = matmul(dnsg2_dkpl, dkpl_dX)

    end subroutine calDXDR_nsg2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calDYDR(obsType, PVs, PVo, m, DYDR)
    implicit none
    integer, intent(in) :: obsType, m
    double precision, intent(in) :: PVs(6), PVo(6)
    double precision, intent(out) :: DYDR(m, 6)
    integer i
    !-------------------------------------------------------------------
    select case(obsType)
    case(OBS_TYPE_RA_DEC)
        call calDYDR_alphaDelta(PVs, PVo, DYDR(1:2, 1:6))
    case(OBS_TYPE_RANGERATE)
        call calDYDR_velocity(PVs, PVo, DYDR(1:1, 1:6))
    case(OBS_TYPE_RANGING)
        call calDYDR_ranging(PVs, PVo, DYDR(1:1, 1:6))
    case(OBS_TYPE_POSITION, OBS_TYPE_POSITION_VELOCITY)
        DYDR = 0d0
        do i = 1, m
            DYDR(i, i) = 1d0
        enddo
    case default
        print*, "error@calDYDR"
        write(*, *) "unknown observational type: ", obsType
    endselect
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calDYDR_alphaDelta(rs, ro, M)
    implicit none
    double precision, intent(in) :: rs(6), ro(6)
    double precision, intent(out) :: M(2, 6)
    double precision dr(3), rho, cosAlpha, sinAlpha, cosDelta, sinDelta
    !--------------------------------------------------------------------------
    dr = rs(1:3) - ro(1:3)
    rho = norm2(dr)
    cosAlpha = dr(1)/norm2(dr(1:2))
    sinAlpha = dr(2)/norm2(dr(1:2))
    cosDelta = norm2(dr(1:2))/rho
    sinDelta = dr(3)/rho
    M(1, 1:6) = (/-sinAlpha, cosAlpha, 0d0, 0d0, 0d0, 0d0/) / rho / cosDelta
    M(2, 1:6) = (/-sinDelta*cosAlpha, -sinDelta*sinAlpha, cosDelta, 0d0, 0d0, 0d0/) / rho
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate partial derivative of range(rho) over satellite position and velocity
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calDYDR_ranging(rs, ro, DYDR)
    implicit none
    double precision, intent(in) :: rs(6), ro(6)
    double precision, intent(out) :: DYDR(1, 6)
    !--------------------------------------------------------------------------
    DYDR(1, 1:3) = 1d0/norm2(rs(1:3) - ro(1:3))*(rs(1:3) - ro(1:3))
    DYDR(1, 4:6) = 0d0
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate partial derivative of velocity(rho_dot) over satellite position and velocity
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calDYDR_velocity(rs, ro, DYDR)
    implicit none
    double precision, intent(in) :: rs(6), ro(6)
    double precision, intent(out) :: DYDR(1, 1:6)
    double precision rho, rhoDot
    !--------------------------------------------------------------------------
    rho = norm2(rs(1:3) - ro(1:3))
    rhoDot = 1d0/rho*dot_product(rs(1:3) - ro(1:3), rs(4:6) - ro(4:6))
    DYDR(1, 1:3) = 1d0/rho*(rs(4:6) - ro(4:6) - rhoDot/rho*(rs(1:3) - ro(1:3)))
    DYDR(1, 4:6) = 1d0/rho*(rs(1:3) - ro(1:3))
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initialize object res, allocate matrix
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine initializeRes(ODor, ss)
    implicit none
    type(type_ODor), intent(in) :: ODor
    type(type_SatelliteScenario), intent(in) :: ss
    integer nAc, i
    nAc = ODor%nAccess
    res%nAcc = nAc
    if (allocated(res%acc)) deallocate(res%acc)
    allocate(res%acc(nAc))
    do i = 1, nAc
        res%acc(i)%accIdx = ODor%accessIdx(i)
        res%acc(i)%obsType = ss%access(ODor%accessIdx(i))%obsType
        res%acc(i)%npe = ss%access(ODor%accessIdx(i))%npe
        res%acc(i)%nData = ss%access(ODor%accessIdx(i))%nData
        allocate(res%acc(i)%timeTable(res%acc(i)%nData))
        allocate(res%acc(i)%res(res%acc(i)%nData, res%acc(i)%npe))
    enddo
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine saveEph(this, ss, fn, ts_in)
    implicit none
    class(type_ODor) :: this
    type(type_satelliteScenario), intent(inOut) :: ss
    character(len=*), intent(in) :: fn
    type(type_timeSeries), intent(in), optional :: ts_in
    type(type_TimeSeries) :: ts
    integer :: nSat, list(100), i
    logical :: dirExists
    
    nSat = this%nSat
    list(1:nSat) = this%satIdx(1:nSat)

    if (present(ts_in)) then
        ts = ts_in
    else
        ts = ss%ts
    endif
    
    do i = 1, nSat
        call ss%sat(list(i))%generateTrajectory(ts%time)
        call ss%sat(list(i))%savetrajectory(STATE_TYPE_CARTESIAN, trim(fn), .false.)
        print"(A,A)", " satellite trajectory is saved into file : ", trim(fn)
    enddo
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ODor_propagate(this, ss, epoch)
    implicit none
    class(type_ODor) :: this
    type(type_satelliteScenario), intent(inOut) :: ss
    type(type_t), intent(in) :: epoch
    integer :: i, iSat
    double precision :: dt, STM(this%nVar, this%nVar)
    
    dt = (epoch%day - this%t%day)*86400d0 + epoch%sec - this%t%sec ! sec
    
    ! if this%t equals to epoch, return
    if (dabs(dt) < 1d-13) then
        return
    endif
    
    ! 2 get STM
    call numericalSTM(ss, dt, this%nVar, STM)
    
    ! 3 update covariance of states
    this%P = matmul(STM, matmul(this%P, transpose(STM)))
     
    ! 4 update satellites' states
    do i = 1, this%nSat
        iSat = this%satIdx(i)
        call ss%sat(iSat)%propagate(epoch)
    enddo
    call this%ss2ODor(ss)

    this%t = epoch
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine numericalSTM(ss, dt, nVar, STM)
    implicit none
    type(type_SatelliteScenario), intent(in) :: ss
    double precision, intent(in) :: dt ! sec
    integer, intent(in) :: nVar
    double precision, intent(out) :: STM(nVar, nVar)
    integer :: i
    type(type_satellite) :: sat0, sat1, sat1p
    type(type_t) :: epoch
    double precision :: x1(6), dx0(6), x1p(6)

    STM = 0d0
    sat0 = ss%sat(1)
    ! STM 用简化力模型: 地球引力场只到 2x2(含 J2), 加日月三体
    sat0%opMod%degree = 2
    sat0%opMod%order = 2
    ! 整组替换: 先全部关闭, 再只开日月三体引力
    call sat0%opMod%resetSwitches()
    sat0%opMod%key_moonGravity = 1
    sat0%opMod%key_sunGravity = 1
    sat0%opMod%usePV = 1
    epoch = timeAdd(sat0%states%t, 0d0, dt)
    
    sat1 = sat0
    call sat1%propagate(epoch)
    x1 = sat1%states%elements
    
    dx0 = (/1d-4, 1d-4, 1d-4, 1d-6, 1d-6, 1d-6/)
    
    do i = 1, 6
        sat1p = sat0
        sat1p%states%elements(i) = sat1p%states%elements(i) + dx0(i)
        call sat1p%propagate(epoch)
        x1p = sat1p%states%elements
        STM(1:6, i) = (x1p - x1)/dx0(i)
    enddo
    
    do i = 7, nVar
        STM(i, i) = 1d0
    enddo
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! subroutine numericalSTM_ECOM
    ! Calculate state transition matrix (STM) for orbit states + ECOM SRP parameters
    !
    ! The STM has dimension 15x15:
    !   - First 6 dimensions: orbital state variables (position/velocity or elements)
    !   - Last 9 dimensions: ECOM SRP parameters [D0, Dc, Ds, Y0, Yc, Ys, B0, Bc, Bs]
    !
    ! STM structure:
    !   STM = [dX/dX0  (6x6)    dX/dP0  (6x9) ]
    !         [dP/dX0  (9x6)    dP/dP0  (9x9) ]
    !
    ! where X = orbital states, P = ECOM parameters
    !
    ! Method: Finite difference (same approach as numericalSTM), but every
    !   propagation uses the fast predictor (module_fastPrediction): J2+C22/S22
    !   hard-coded, EOP = 0, analytical Sun/Moon, fixed 120 s steps.  The
    !   trajectories are therefore less accurate than NumOP's in absolute terms,
    !   but the two members of each finite difference share the same integrator,
    !   so their difference -- what the STM actually is -- keeps the accuracy,
    !   at a fraction of the cost.
    !
    ! Input:
    !   sat     type_satellite    Satellite object with initial state and force model
    !   dt      double precision  Time step (normalized unit, TU)
    ! Output:
    !   STM     double precision(15, 15)  State transition matrix
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine numericalSTM_ECOM(sat, dt, STM)
    implicit none
    type(type_satellite), intent(in) :: sat
    double precision, intent(in) :: dt ! 预报时长, 时间单位(TU)
    double precision, intent(out) :: STM(15, 15)

    type(type_satellite) :: sat0
    type(type_fastOption) :: opt, optp
    double precision :: x0(6), x1(6), xpert(6), xr(6), dx0(6), dECOM(9), dtSec, t0MJD
    integer :: i
    !--------------------------------------------------------------------------
    STM = 0d0

    sat0 = sat
    ! 快速预报的力模型: J2 + C22/S22 + 日月三体 + 太阳光压 —— 与原先 NumOP 的
    ! 简化模型一致(原先是 degree=2/order=2 的完整 2x2 场, 这里省掉量级 ~1e-10
    ! 的 C21/S21)。开关按输入的卫星设置, 见 satFastOption。
    opt = satFastOption(sat0)

    x0    = sat0%states%elements
    t0MJD = sat0%states%t%val        ! UTC
    dtSec = dt*time_unit             ! TU -> s

    ! 参考轨迹(快速预报, 定步长 2 min)
    call fastPrediction(x0, t0MJD, dtSec, opt, x1)

    !--------------------------------------------------------------------------
    ! Part 1: 状态量列 dX/dX0 与 ECOM 列 dX/dE(有限差分)
    !--------------------------------------------------------------------------
    dx0 = (/1d-4, 1d-4, 1d-4, 1d-4, 1d-4, 1d-4/)
    do i = 1, 6
        xpert = x0
        xpert(i) = xpert(i) + dx0(i)
        call fastPrediction(xpert, t0MJD, dtSec, opt, xr)
        STM(1:6, i) = (xr - x1)/dx0(i)
    enddo

    ! ECOM 参数量级 ~1e-8; 差分步长仍取 1e-6 —— 响应对参数是线性的, 步长不影响结果
    dECOM = 1d-6
    do i = 1, 9
        optp = opt
        optp%ECOM(i) = optp%ECOM(i) + dECOM(i)
        call fastPrediction(x0, t0MJD, dtSec, optp, xr)
        STM(1:6, 6 + i) = (xr - x1)/dECOM(i)
    enddo

    !--------------------------------------------------------------------------
    ! Part 2: ECOM parameters do not depend on initial orbital states (dP/dX0 = 0)
    !--------------------------------------------------------------------------
    ! STM(7:15, 1:6) = 0  (already zero from initialization)

    !--------------------------------------------------------------------------
    ! Part 3: ECOM parameters are constant (dP/dP0 = I)
    !--------------------------------------------------------------------------
    do i = 1, 9
        STM(6 + i, 6 + i) = 1d0
    enddo

    endsubroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! type_STM 的方法
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !--------------------------------------------------------------------------
    ! SUBROUTINE: buildSTMseries
    ! 重建各变量组的 STM 时间序列 —— 不管 label 是什么都先生成数组, 之后 calB
    ! 一律从数组里插值取用, 逻辑统一。label 1~4 用解析 STM(calSTM), label 5
    ! (状态量 + ECOM)用数值 STM(numericalSTM_ECOM, 很贵, 预计算收益最大)。
    !
    ! 输入:
    !   ss - 卫星场景, 提供各组所属卫星、场景单元类型与时间序列
    !--------------------------------------------------------------------------
    subroutine buildSTMseries(this, ss)
    implicit none
    class(type_ODor), intent(inout) :: this
    type(type_satelliteScenario), intent(in) :: ss
    integer :: i, origin

    ! 说明: 钟差组(label 11~13)有意不建 STM 序列 —— 钟差参数不经过轨道状态,
    ! 它的偏导是 calB_clock 里的直接项。钟差组的 STMseries(i)%valid 保持 .false.,
    ! calB 也在取 STM 之前就把它分流走了, 走不到 stm_series_get 的零返回路径。
    do i = 1, this%nGroup
        if (this%group(i)%label >= 1 .and. this%group(i)%label <= 5) then
            origin = this%group(i)%origin
            call this%STMseries(i)%init(ss%sat(origin), this%group(i)%label, &
                &   ss%eleType, ss%ts%t0, ss%ts%tEnd, &
                &   this%stmStepSec, this%group(i)%nVar)
        elseif (this%STMseries(i)%valid) then
            call this%STMseries(i)%free()
        endif
    enddo
    end subroutine buildSTMseries

    !--------------------------------------------------------------------------
    ! SUBROUTINE: stm_series_init
    ! 按等间隔网格生成一条 STM 序列: stm(k, :, :) = dX(t_k)/dp。
    !
    ! 数值列(笛卡尔状态列、ECOM 参数列)的做法: **每一列只做一次贯穿整个网格的
    ! 预报** —— 给初值(或参数)加一个偏差后逐节点推进, 在每个节点记下差商。原先
    ! 是按节点循环、每个节点都从基点重新积分一遍, 总积分量 O(n²); 现在按列循环,
    ! 总积分量 O(n), 省下约 n/2 倍(本场景 n = 145)。
    ! 网格等间隔, 相邻节点之间正好是整数个快速预报步(600 s = 5×120 s), 因此各列
    ! 与标称走的步点序列完全一致, 结果与旧写法在舍入级别上相同。
    !
    ! 输入:
    !   sat     - 参考卫星(沿它的状态线性化; 用 sat%states%t 作传播基点)
    !   t0      - 网格起始历元
    !   tEnd    - 网格结束历元(应覆盖全部观测历元)
    !   stepSec - 网格间隔 [s]
    !   nVar    - 该组参数个数(应为 6 + m, m <= 9)
    !--------------------------------------------------------------------------
    subroutine stm_series_init(this, sat, label, eleType, t0, tEnd, stepSec, nVar)
    implicit none
    class(type_STM), intent(inout) :: this
    type(type_satellite), intent(in) :: sat
    integer, intent(in) :: label, eleType
    type(type_t), intent(in) :: t0, tEnd
    double precision, intent(in) :: stepSec
    integer, intent(in) :: nVar
    integer :: k, n, i, j, nCol
    double precision :: dt, STM8(8, 8)
    double precision :: dSec, hCol, tNodeMJD
    double precision :: x0(6), xc(6), xr(6)
    double precision, allocatable :: xNom(:, :)
    type(type_fastOption) :: opt, optp
    type(type_t) :: tk

    call this%free()
    if (stepSec <= 0d0 .or. nVar < 6 .or. nVar > 15) return
    if (label < 1 .or. label > 5) return

    n = floor(timeMinus(tEnd, t0)*86400d0/stepSec) + 1
    if (n < 2) return

    this%nEpoch = n
    this%nVar = nVar
    this%t0 = t0
    this%stepSec = stepSec
    allocate(this%stm(n, nVar, nVar))
    this%stm = 0d0

    ! 参数本身是常数: dP/dP0 = I
    do k = 1, n
        do i = 7, nVar
            this%stm(k, i, i) = 1d0
        enddo
    enddo

    !------------------------------------------------------------------
    ! 情形一: 解析 STM(非笛卡尔状态的 label 1~4)。逐节点直接求值, 不涉及积分。
    ! dt 取时间单位(TU), 基点是 sat%states%t。
    !------------------------------------------------------------------
    if (label <= 4 .and. eleType /= STATE_TYPE_CARTESIAN) then
        nCol = min(nVar, 7)
        do k = 1, n
            tk = timeAdd(t0, dble(k-1)*stepSec/86400d0)
            dt = timeMinus(tk, sat%states%t)*86400d0/time_unit
            call calSTM(eleType, sat%states%elements, dt, &
                &   label == 2 .or. label == 4, label == 3 .or. label == 4, STM8)
            this%stm(k, 1:6, 1:nCol) = STM8(1:6, 1:nCol)
        enddo
        this%valid = .true.
        return
    endif

    !------------------------------------------------------------------
    ! 情形二: 数值 STM(笛卡尔状态的 label 1~4, 以及 label 5)。
    ! 状态列 6 个; label 5 再加前 (nVar-6) 个 ECOM 参数列。
    !------------------------------------------------------------------
    if (label > 1 .and. label < 5) then
        print*, "warning @stm_series_init: numerical STM gives state columns only;"// &
            &   " for a Cartesian state the amr column of label ", label, " is set to 0"
    endif

    nCol = 6
    if (label == 5) nCol = nVar         ! 6 状态列 + 前 (nVar-6) 个 ECOM 参数列
    x0 = sat%states%elements             ! 归一化笛卡尔状态(参数基点)
    opt = satFastOption(sat)
    dSec = timeMinus(t0, sat%states%t)*86400d0     ! 基点 -> 网格首点

    ! 标称轨迹: 基点 -> 首点, 再逐节点推进
    allocate(xNom(6, n))
    call fastPrediction(x0, sat%states%t%val, dSec, opt, xNom(1:6, 1))
    tNodeMJD = t0%val
    do k = 2, n
        call fastPrediction(xNom(1:6, k-1), tNodeMJD, stepSec, opt, xNom(1:6, k))
        tNodeMJD = tNodeMJD + stepSec/86400d0
    enddo

    ! 扰动列: 沿同一条节点路径走一遍, 逐节点记差商
    do j = 1, nCol
        xc = x0
        optp = opt
        if (j <= 6) then
            hCol = 1d-4                  ! 归一化状态步长(与 numericalSTM_ECOM 一致)
            xc(j) = xc(j) + hCol
        else
            hCol = 1d-6                  ! ECOM 参数步长(同上)
            optp%ECOM(j-6) = optp%ECOM(j-6) + hCol
        endif

        call fastPrediction(xc, sat%states%t%val, dSec, optp, xr)
        this%stm(1, 1:6, j) = (xr - xNom(1:6, 1))/hCol
        tNodeMJD = t0%val
        do k = 2, n
            xc = xr
            call fastPrediction(xc, tNodeMJD, stepSec, optp, xr)
            this%stm(k, 1:6, j) = (xr - xNom(1:6, k))/hCol
            tNodeMJD = tNodeMJD + stepSec/86400d0
        enddo
    enddo

    deallocate(xNom)
    this%valid = .true.
    end subroutine stm_series_init

    !---------------------------------------------------------------------------
    ! FUNCTION: satFastOption
    ! 由卫星的 opMod 组装快速预报的力模型选项(数值 STM 与状态列差分共用)
    !---------------------------------------------------------------------------
    function satFastOption(sat) result(opt)
    implicit none
    type(type_satellite), intent(in) :: sat
    type(type_fastOption) :: opt

    opt%useJ2   = .true.
    opt%useJ22  = .true.
    opt%useMoon = (sat%opMod%key_moonGravity == 1)
    opt%useSun  = (sat%opMod%key_sunGravity == 1)
    opt%useSRP  = (sat%opMod%key_srp == 1)
    opt%SRP_Model    = sat%opMod%SRP_Model
    opt%reflectivity = sat%opMod%reflectivity
    opt%amrSRP       = sat%opMod%amrSRP
    opt%ECOM         = sat%opMod%ECOM
    end function satFastOption

    !--------------------------------------------------------------------------
    ! SUBROUTINE: stm_series_get
    ! 取指定时刻的 STM: 网格上 4 点 Lagrange 插值(节点规则与 clock/edit 的
    ! 插值一致)。时刻落在网格外时夹到最近的端点, 不报错 —— 观测历元可能因
    ! 钟差标记而比网格起点早几个 ns。
    !
    ! 输入:
    !   t - 查询时刻
    ! 输出:
    !   M - STM(nVar, nVar); 序列未生成时返回 0
    !--------------------------------------------------------------------------
    subroutine stm_series_get(this, t, M)
    implicit none
    class(type_STM), intent(in) :: this
    type(type_t), intent(in) :: t
    double precision, intent(out) :: M(:, :)
    integer :: k, j, ii, nodes(4), nE, nv
    double precision :: xt, x(4), w(4)

    M = 0d0
    if (.not. this%valid) return

    nv = this%nVar
    nE = this%nEpoch
    if (size(M, 1) /= nv .or. size(M, 2) /= nv) return

    xt = timeMinus(t, this%t0)*86400d0/this%stepSec + 1d0   ! 1 起的连续网格下标
    if (xt <= 1d0) then
        M = this%stm(1, 1:nv, 1:nv)
        return
    elseif (xt >= dble(nE)) then
        M = this%stm(nE, 1:nv, 1:nv)
        return
    endif
    k = int(xt)

    if (nE >= 4) then
        if (k <= 1) then
            nodes = (/1, 2, 3, 4/)
        elseif (k >= nE - 2) then
            nodes = (/nE-3, nE-2, nE-1, nE/)
        else
            nodes = (/k-1, k, k+1, k+2/)
        endif
        do j = 1, 4
            x(j) = dble(nodes(j) - 1)*this%stepSec        ! 相对 t0 的秒数
        enddo
        xt = timeMinus(t, this%t0)*86400d0
        do j = 1, 4
            w(j) = 1d0
            do ii = 1, 4
                if (ii /= j) w(j) = w(j)*(xt - x(ii))/(x(j) - x(ii))
            enddo
            M = M + w(j)*this%stm(nodes(j), 1:nv, 1:nv)
        enddo
    else
        M = this%stm(k, 1:nv, 1:nv)
    endif
    end subroutine stm_series_get

    !--------------------------------------------------------------------------
    ! SUBROUTINE: stm_series_free
    ! 释放序列
    !--------------------------------------------------------------------------
    subroutine stm_series_free(this)
    implicit none
    class(type_STM), intent(inout) :: this

    if (allocated(this%stm)) deallocate(this%stm)
    this%valid = .false.
    this%nEpoch = 0
    end subroutine stm_series_free

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine processNoise_cartesian(dt, sigma_ax, sigma_ay, sigma_az, Q)

    ! Computes process noise covariance matrix for 3D tracking
    ! State vector: [x, y, z, vx, vy, vz]
    ! Process noise from random accelerations [ax, ay, az]
    
    implicit none
    double precision, intent(in) :: dt, sigma_ax, sigma_ay, sigma_az
    double precision, intent(out) :: Q(6,6)
    
    double precision :: dt2, dt3, dt4
    double precision :: sax2, say2, saz2
    
    ! Initialize Q to zero
    Q = 0.0d0
    
    ! Precompute values
    dt2 = dt * dt
    dt3 = dt2 * dt
    dt4 = dt3 * dt
    
    sax2 = sigma_ax * sigma_ax
    say2 = sigma_ay * sigma_ay
    saz2 = sigma_az * sigma_az
    
    ! Fill Q matrix
    Q(1,1) = 0.25d0 * sax2 * dt4
    Q(2,2) = 0.25d0 * say2 * dt4
    Q(3,3) = 0.25d0 * saz2 * dt4
    
    Q(1,4) = 0.5d0 * sax2 * dt3
    Q(2,5) = 0.5d0 * say2 * dt3
    Q(3,6) = 0.5d0 * saz2 * dt3
    
    Q(4,1) = Q(1,4)
    Q(5,2) = Q(2,5)
    Q(6,3) = Q(3,6)
    
    Q(4,4) = sax2 * dt2
    Q(5,5) = say2 * dt2
    Q(6,6) = saz2 * dt2
    
    end subroutine processNoise_cartesian

    end module