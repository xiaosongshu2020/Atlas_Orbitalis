!================================================================================
! MODULE: module_clock
!
! 说明：
!   本模块定义导航卫星的钟模型对象 type_clock，用于描述卫星钟相对于参考时间
!   系统的钟差及其变化规律，支持由钟模型生成(或直接读入)时钟状态表，再从表中
!   查询任意时刻的钟差与钟漂移率。
!
!   钟差符号约定：bias(t) = 卫星钟面时间 - 参考时间，即钟差为正表示卫星钟走快。
!   钟差、钟噪声的单位为秒(s)，钟漂移率为 s/s，二阶钟漂移率为 s/s^2。
!
! 属性：
!   name        - 钟所属对象名(如卫星名)，仅用于文件头与显示，不参与计算
!   t0          - 钟模型参考时刻 (type_t, UTC)
!   bias0       - 参考时刻的钟差 [s]
!   drift       - 钟漂移率(频偏) [s/s]
!   drift2      - 二阶钟漂移率(频漂) [s/s^2]
!   noise       - 钟噪声：钟差白噪声标准差 [s]
!   calcMethod  - 计算方法：CLOCK_METHOD_FORMULA (1，公式方法，默认)
!                           CLOCK_METHOD_INTERPOLATION (2，插值方法)
!   table_clock - 时钟状态表，每条记录含时刻(type_t)与钟差(s)
!
! 钟模型：
!   bias(t) = bias0 + drift*dt + drift2*dt^2 ,  dt = t - t0 [s]
!   (drift2 是 dt^2 项的系数，不除以 2；钟漂移率 v(t) = drift + 2*drift2*dt)
!
! 方法：
!   init(t0, bias0, drift, drift2, noise, calcMethod) - 初始化钟模型
!                          (整体复位，清除状态表)
!   setModel(t0, bias0, drift, drift2) - 直接设置解析钟模型参数并强制公式方法；
!                          不清表、不复位 name/noise，供定轨迭代中更新钟参数使用
!   getModelCoeffs(tRef, c(3)) - 把钟模型精确重基准化到 tRef，返回新的
!                          (bias0, drift, drift2)。换参考历元必须走这里，
!                          否则等于悄悄改变了钟函数
!   (内部函数 clock_get_model_state：由解析模型算钟差，供状态表生成与公式方法
!    调用；不是类型方法，模块外不可见)
!   generateTable(epochs) - 按给定历元生成状态表(模型钟差 + 白噪声)
!   setTable(n, times, bias) - 直接设置状态表(如读入实测/解算得到的钟差序列)
!   get_clock_state(t)   - 查询时刻 t 的钟差。公式方法(calcMethod = 1)用解析
!                          模型；插值方法(calcMethod = 2)用状态表插值(4点
!                          Lagrange，与 type_satellite%getStates 的星历插值
!                          同法)。选了插值但状态表未定义时，报警告并退回公式方法
!   get_clock_drift(t)   - 查询时刻 t 的钟漂移率。公式方法由解析钟模型求导，
!                          v(t) = drift + 2*drift2*dt；插值方法取 t 所在区间的
!                          两个数据点做线性差分 v = (bias2 - bias1)/(t2 - t1)。
!                          状态表未定义时同样退回公式方法
!   saveData(fn, name)   - 把钟模型与状态表写成文本文件(格式见下)
!   readData(fn)         - 从文本文件恢复钟模型与状态表；文件头声明的模型参数、
!                          名称与计算方法覆盖对象当前值，即文件为准
!   free()               - 释放状态表
!   show()               - 显示钟模型信息
!
! 文件格式(saveData / readData，与 module_access 的观测数据文件同一约定)：
!   文件头每行以 '#' 开头，自描述：生成者、对象名、数据类型、计算方法、钟模型
!   参数(bias0/drift/drift2/noise/t0)、历元数、时间跨度、数据格式；其后每行一个
!   历元："MJD日  MJD秒  钟差/s"。读取时跳过文件头，并把文件头声明的量恢复到对象
!   上；格式不对的数据行报警告并跳过(不补零)；文件头未声明计算方法而文件中有数据
!   行时，默认改用插值方法(否则读进来的表不会被使用)。
!
! 依赖：
!   module_time     - type_t 与时间运算(timeMinus/timeAdd)
!   module_toolKits - Interpolation(Lagrange 插值)
!   module_const    - 物理常数
!   module_randum   - normal(mean, sigma) 高斯随机数
!
! 用法：
!   type(type_clock) :: clk
!   type(type_t) :: t0, epochs(100)
!   call t0%init(2025d0, 1d0, 1d0, 0d0, 0d0, 0d0)             ! 2025-01-01 00:00:00 UTC
!   call clk%init(t0, bias0=1d-4, drift=1d-11, drift2=0d0, noise=1d-10)
!   do i = 1, 100
!       epochs(i) = timeAdd(t0, 30d0*dble(i-1)/86400d0)       ! 30 s 间隔
!   enddo
!   call clk%generateTable(epochs)
!   clk%calcMethod = CLOCK_METHOD_INTERPOLATION   ! 改用状态表插值(缺省是公式方法)
!   bias  = clk%get_clock_state(t)     ! 钟差 [s]
!   drift = clk%get_clock_drift(t)     ! 钟漂移率 [s/s]
!   call clk%show()
!   call clk%saveData("data/output/clock_LEO.txt", name="LEO-600km")   ! 写文件
!   call clk%readData("data/output/clock_LEO.txt")                    ! 再读回来
!
! 注：
!   - 状态表须按时间递增排列(查询按此前提做二分查找)；
!   - 查询时刻超出状态表范围时报错并 stop，与 type_satellite%getStates 一致；
!   - 所有历元应属于同一时间系统(建议 UTC)。
!
! 作者：Shawn Hu
! 最后修改：2026-10-08
!================================================================================
module module_clock

    use module_time
    use module_toolKits
    use module_const
    use module_randum

    implicit none
    private

    public :: type_clock, type_clockRecord
    public :: CLOCK_METHOD_FORMULA, CLOCK_METHOD_INTERPOLATION

    ! 钟差的计算方法
    integer, parameter :: CLOCK_METHOD_FORMULA       = 1  ! 公式方法(用解析钟模型)
    integer, parameter :: CLOCK_METHOD_INTERPOLATION = 2  ! 插值方法(查时钟状态表)

    !---------------------------------------------------------------------------
    ! TYPE: type_clockRecord
    ! 时钟状态表中的一条记录：历元 + 钟差
    !---------------------------------------------------------------------------
    type :: type_clockRecord
        type(type_t)     :: t              ! 历元 (UTC)
        double precision :: bias = 0d0     ! 钟差 [s]
    end type type_clockRecord

    !---------------------------------------------------------------------------
    ! TYPE: type_clock
    ! 卫星钟模型：解析模型参数 + 时钟状态表
    !---------------------------------------------------------------------------
    type :: type_clock
        ! 钟所属对象名(如卫星名)：仅用于文件头与显示，不参与任何计算
        character(len=100) :: name = ''
        ! 钟模型参数(解析模型)
        type(type_t)     :: t0             ! 参考时刻
        double precision :: bias0 = 0d0    ! 参考时刻钟差 [s]
        double precision :: drift = 0d0    ! 钟漂移率 [s/s]
        double precision :: drift2 = 0d0   ! 二阶钟漂移率(dt^2 项系数) [s/s^2]
        double precision :: noise = 0d0    ! 钟噪声(钟差白噪声标准差) [s]

        ! 钟差的计算方法：1 = 公式方法(默认)，2 = 插值方法
        ! 选 2 但状态表未定义时，查询会报警告并退回公式方法
        integer :: calcMethod = CLOCK_METHOD_FORMULA

        ! 时钟状态表
        type(type_clockRecord), allocatable :: table_clock(:)
    contains
        procedure :: init => clock_init
        procedure :: setModel => clock_set_model
        procedure :: getModelCoeffs => clock_get_model_coeffs
        procedure :: get_clock_state => clock_get_state
        procedure :: get_clock_drift => clock_get_drift
        procedure :: generateTable => clock_generate_table
        procedure :: setTable => clock_set_table
        procedure :: saveData => clock_save_data
        procedure :: readData => clock_read_data
        procedure :: free => clock_free
        procedure :: show => clock_show
    end type type_clock

contains

    !---------------------------------------------------------------------------
    ! SUBROUTINE: clock_init
    ! 初始化钟模型。整体复位：重置全部模型参数，并清除已有的状态表。
    !
    ! 输入：
    !   t0         - 参考时刻 (type_t)
    !   bias0      - (可选) 参考时刻钟差 [s]，缺省 0
    !   drift      - (可选) 钟漂移率 [s/s]，缺省 0
    !   drift2     - (可选) 二阶钟漂移率 [s/s^2]，缺省 0
    !   noise      - (可选) 钟噪声标准差 [s]，缺省 0(不加噪声)
    !   calcMethod - (可选) 计算方法，缺省 CLOCK_METHOD_FORMULA(公式方法)
    !   name       - (可选) 钟所属对象名(如卫星名)；不给时保留原值(与 module_access
    !                       的 initAccess 一致，对象名不随模型复位而丢失)
    !---------------------------------------------------------------------------
    subroutine clock_init(this, t0, bias0, drift, drift2, noise, calcMethod, name)
    implicit none
    class(type_clock), intent(inout) :: this
    type(type_t), intent(in) :: t0
    double precision, optional, intent(in) :: bias0, drift, drift2, noise
    integer, optional, intent(in) :: calcMethod
    character(len=*), optional, intent(in) :: name

    this%t0 = t0
    this%bias0 = 0d0
    this%drift = 0d0
    this%drift2 = 0d0
    this%noise = 0d0
    if (present(bias0))  this%bias0  = bias0
    if (present(drift))  this%drift  = drift
    if (present(drift2)) this%drift2 = drift2
    if (present(noise))  this%noise  = noise
    if (present(name))   this%name   = name
    this%calcMethod = CLOCK_METHOD_FORMULA
    if (present(calcMethod)) this%calcMethod = calcMethod

    ! 清除已有的状态表
    if (allocated(this%table_clock)) deallocate(this%table_clock)
    end subroutine clock_init

    !---------------------------------------------------------------------------
    ! SUBROUTINE: clock_set_model
    ! 直接设置解析钟模型的参数，并强制 calcMethod = FORMULA(公式方法)。
    !
    ! 与 init 的区别：不清除状态表，也不重置 name/noise —— 供定轨迭代中更新钟参数
    ! 使用(init 会清表，不能在迭代中调用)。强制公式方法的原因：估计出来的钟参数
    ! 只有在公式方法下才会被 get_clock_state 采用；插值方法读的是状态表，而且
    ! 查询超出表范围时会直接 stop。
    !
    ! 输入：
    !   t0, bias0, drift, drift2 - 钟模型参数
    !---------------------------------------------------------------------------
    subroutine clock_set_model(this, t0, bias0, drift, drift2)
    implicit none
    class(type_clock), intent(inout) :: this
    type(type_t), intent(in) :: t0
    double precision, intent(in) :: bias0, drift, drift2

    this%t0     = t0
    this%bias0  = bias0
    this%drift  = drift
    this%drift2 = drift2
    this%calcMethod = CLOCK_METHOD_FORMULA
    end subroutine clock_set_model

    !---------------------------------------------------------------------------
    ! SUBROUTINE: clock_get_model_coeffs
    ! 把解析钟模型精确重基准化到参考时刻 tRef，返回三个系数：
    !   bias(t) = c(1) + c(2)*(t - tRef) + c(3)*(t - tRef)^2
    ! 即同一个钟函数在 tRef 处的展开。换参考历元必须走这里 —— 只改 t0 而不变换
    ! 系数会悄悄改变钟函数(例如只估 bias0 时，要保留的 drift 是 drift + 2*drift2*Δ)。
    !
    ! 输入：
    !   tRef - 新的参考时刻 (type_t)
    ! 输出：
    !   c(3) - 重基准化后的 (bias0, drift, drift2)
    !---------------------------------------------------------------------------
    subroutine clock_get_model_coeffs(this, tRef, c)
    implicit none
    class(type_clock), intent(in) :: this
    type(type_t), intent(in) :: tRef
    double precision, intent(out) :: c(3)
    double precision :: dt

    dt = timeMinus(tRef, this%t0)*86400d0
    c(1) = this%bias0 + this%drift*dt + this%drift2*dt*dt
    c(2) = this%drift + 2d0*this%drift2*dt
    c(3) = this%drift2
    end subroutine clock_get_model_coeffs

    !---------------------------------------------------------------------------
    ! FUNCTION: clock_get_model_state
    ! 由解析钟模型计算指定时刻的钟差(确定性部分，不含噪声)
    !
    ! 模块内部函数，不是 type_clock 的方法，不对外暴露 —— 外部查询钟差请用
    ! get_clock_state / get_clock_drift。此处仅供状态表生成与公式方法调用。
    !
    !   bias(t) = bias0 + drift*dt + drift2*dt^2 ,  dt = t - t0 [s]
    !   (drift2 是 dt^2 项的系数，不除以 2)
    !
    ! 输入：
    !   t      - 查询时刻 (type_t)
    ! 返回：
    !   bias   - 模型钟差 [s]
    !---------------------------------------------------------------------------
    function clock_get_model_state(this, t) result(bias)
    implicit none
    class(type_clock), intent(in) :: this
    type(type_t), intent(in) :: t
    double precision :: bias
    double precision :: dt

    dt = timeMinus(t, this%t0)*86400d0
    bias = this%bias0 + this%drift*dt + this%drift2*dt*dt
    end function clock_get_model_state

    !---------------------------------------------------------------------------
    ! SUBROUTINE: clock_generate_table
    ! 按给定历元生成时钟状态表：bias_i = 模型钟差(epochs(i)) + 白噪声
    ! (noise > 0 时叠加)。会覆盖已有状态表。
    !
    ! 输入：
    !   epochs - 历元数组(type_t)，须按时间递增排列
    !---------------------------------------------------------------------------
    subroutine clock_generate_table(this, epochs)
    implicit none
    class(type_clock), intent(inout) :: this
    type(type_t), intent(in) :: epochs(:)
    integer :: i, n

    n = size(epochs)
    if (n <= 0) then
        write(0, *) "Error @generateTable: empty epoch array."
        return
    endif

    if (allocated(this%table_clock)) deallocate(this%table_clock)
    allocate(this%table_clock(n))
    do i = 1, n
        this%table_clock(i)%t = epochs(i)
        this%table_clock(i)%bias = clock_get_model_state(this, epochs(i))
        ! 叠加钟噪声(钟差白噪声)
        if (this%noise > 0d0) &
            this%table_clock(i)%bias = this%table_clock(i)%bias + normal(0d0, this%noise)
    enddo
    end subroutine clock_generate_table

    !---------------------------------------------------------------------------
    ! SUBROUTINE: clock_set_table
    ! 直接设置时钟状态表(用于读入实测/解算得到的钟差序列)。
    ! 仅设置状态表，不改变钟模型参数。
    !
    ! 输入：
    !   n     - 记录数
    !   times - 各记录时刻 (n)，须按时间递增排列
    !   bias  - 各记录钟差 (n) [s]
    !---------------------------------------------------------------------------
    subroutine clock_set_table(this, n, times, bias)
    implicit none
    class(type_clock), intent(inout) :: this
    integer, intent(in) :: n
    type(type_t), intent(in) :: times(n)
    double precision, intent(in) :: bias(n)
    integer :: i

    if (allocated(this%table_clock)) deallocate(this%table_clock)
    if (n > 0) then
        allocate(this%table_clock(n))
        do i = 1, n
            this%table_clock(i)%t = times(i)
            this%table_clock(i)%bias = bias(i)
        enddo
    endif
    end subroutine clock_set_table

    !---------------------------------------------------------------------------
    ! SUBROUTINE: clock_save_data
    ! 把钟模型与时钟状态表写成文本文件。
    !
    ! 文件头每行以 '#' 开头，自描述(生成者、对象名、数据类型、计算方法、钟模型
    ! 参数、历元数、时间跨度、数据格式)；其后每行一个历元：
    !   MJD日  MJD秒  钟差/s
    ! 状态表未生成时只写文件头，文件即"钟模型"的描述(读回时只恢复模型)。
    !
    ! 输入：
    !   fn   - 输出文件名
    !   name - (可选) 钟所属对象名(如卫星名)；给出时同时存入本对象的 name
    !---------------------------------------------------------------------------
    subroutine clock_save_data(this, fn, name)
    implicit none
    class(type_clock), intent(inout) :: this
    character(len=*), intent(in) :: fn
    character(len=*), optional, intent(in) :: name
    integer :: unitNum, i, n
    type(type_t) :: epoch

    if (present(name)) this%name = name

    if (allocated(this%table_clock)) then
        n = size(this%table_clock)
    else
        n = 0
    endif

    open(newunit = unitNum, file = trim(fn))
    call clock_write_header(unitNum, this, n)
    if (n > 0) then
        write(unitNum, "(A,2x,A5,2x,A13,2x,A18)") "#", "Day", "Sec", "Bias/s"
        write(unitNum, "(A,2x,5('-'),2x,13('-'),2x,18('-'))") "#"
        do i = 1, n
            epoch = this%table_clock(i)%t
            write(unitNum, "(2x,I5,2x,F13.6,2x,ES18.10)") &
                nint(epoch%day), epoch%sec, this%table_clock(i)%bias
        enddo
    endif
    close(unitNum)
    end subroutine clock_save_data

    !---------------------------------------------------------------------------
    ! SUBROUTINE: clock_read_data
    ! 从文本文件恢复钟模型与时钟状态表(clock_save_data 写出的格式)。
    !
    ! 先读文件头(见 clock_read_header)：文件头声明的模型参数、对象名与计算方法
    ! 覆盖对象当前值 —— 文件为准。再读数据行，覆盖已有状态表：
    !   - 文件头未声明计算方法、而文件中确有数据行时，改用插值方法(否则读进来的
    !     钟差表不会被使用，见 clock_effective_method)；
    !   - 格式不对的数据行由 readDataMatrix 报警告并跳过(不补零)，记录数按实际
    !     读到的行数计；
    !   - 文件中没有数据行时只恢复钟模型(允许只有模型、没有状态表的文件)。
    !
    ! 输入：
    !   fn - 输入文件名
    !---------------------------------------------------------------------------
    subroutine clock_read_data(this, fn)
    implicit none
    class(type_clock), intent(inout) :: this
    character(len=*), intent(in) :: fn
    double precision, allocatable :: x(:, :)
    integer :: nData, nRead, i
    logical :: methodDeclared

    ! 文件头：恢复钟模型参数、对象名与计算方法
    call clock_read_header(trim(fn), this, methodDeclared)

    call getDataLinesNum(trim(fn), nData)
    if (nData == 0) then
        print*, "warning @clock%readData: no data line in ", trim(fn)
        print*, "  (only the clock model, as declared in the header, was restored)"
        return
    endif

    allocate(x(nData, 3))
    call readDataMatrix(trim(fn), nData, 3, x, nRead)
    if (nRead == 0) then
        print*, "warning @clock%readData: no usable data line in ", trim(fn)
        deallocate(x)
        return
    endif
    nData = nRead

    if (allocated(this%table_clock)) deallocate(this%table_clock)
    allocate(this%table_clock(nData))
    do i = 1, nData
        call this%table_clock(i)%t%init(x(i, 1), x(i, 2))
        this%table_clock(i)%bias = x(i, 3)
    enddo

    ! 文件头没有声明计算方法：用插值方法，否则刚读入的钟差表不会被使用
    if (.not. methodDeclared) this%calcMethod = CLOCK_METHOD_INTERPOLATION

    deallocate(x)
    end subroutine clock_read_data

    !---------------------------------------------------------------------------
    ! SUBROUTINE: clock_write_header
    ! 写出 self-describing 文件头：每行以 '#' 开头，"# 关键字 : 值"。
    ! 关键字左对齐到 12 列，与 module_access 写观测数据文件头的写法一致。
    !
    ! 输入：
    !   unitNum - 已打开的输出单元
    !   n       - 状态表记录数(0 表示只有钟模型，没有状态表)
    !---------------------------------------------------------------------------
    subroutine clock_write_header(unitNum, this, n)
    implicit none
    integer, intent(in) :: unitNum, n
    class(type_clock), intent(in) :: this
    character(len=32) :: numStr, numStr2
    character(len=16) :: methodName
    double precision :: mjdFirst, mjdLast

    write(unitNum, "(A)") repeat("#", 78)
    write(unitNum, "(A)") "# Atlas Orbitalis - satellite clock offset data"
    call write_header_line(unitNum, "generator", "Atlas (Atlas Orbitalis)")
    if (len_trim(this%name) > 0) call write_header_line(unitNum, "object", trim(this%name))
    call write_header_line(unitNum, "data type", "CLOCK_OFFSET (seconds)")

    select case (this%calcMethod)
    case (CLOCK_METHOD_INTERPOLATION)
        methodName = "INTERPOLATION"
    case (CLOCK_METHOD_FORMULA)
        methodName = "FORMULA"
    case default
        methodName = "UNKNOWN"
    end select
    write(numStr, "(I6)") this%calcMethod
    call write_header_line(unitNum, "calc method", &
        trim(methodName)//" ("//trim(adjustl(numStr))//")")

    ! 解析钟模型：bias(t) = bias0 + drift*dt + drift2*dt^2, dt = t - t0 [s]
    write(numStr, "(ES16.8)") this%bias0
    call write_header_line(unitNum, "bias0", trim(adjustl(numStr))//" s")
    write(numStr, "(ES16.8)") this%drift
    call write_header_line(unitNum, "drift", trim(adjustl(numStr))//" s/s")
    write(numStr, "(ES16.8)") this%drift2
    call write_header_line(unitNum, "drift2", trim(adjustl(numStr))//" s/s^2")
    write(numStr, "(ES16.8)") this%noise
    call write_header_line(unitNum, "noise", trim(adjustl(numStr))//" s")

    write(numStr, "(F18.12)") this%t0%day + this%t0%sec/86400d0
    call write_header_line(unitNum, "t0", trim(adjustl(numStr))//" MJD (UTC)  = "//trim(this%t0%toDate()))

    write(numStr, "(I9)") n
    call write_header_line(unitNum, "epochs", trim(adjustl(numStr)))
    if (n > 0) then
        mjdFirst = this%table_clock(1)%t%day + this%table_clock(1)%t%sec/86400d0
        mjdLast  = this%table_clock(n)%t%day + this%table_clock(n)%t%sec/86400d0
        write(numStr,  "(F15.9)") mjdFirst
        write(numStr2, "(F15.9)") mjdLast
        call write_header_line(unitNum, "time span", &
            trim(adjustl(numStr))//" ... "//trim(adjustl(numStr2))//" MJD (UTC)")
    endif
    call write_header_line(unitNum, "data format", "MJD_day  MJD_sec  Bias/s")
    call write_header_line(unitNum, "note", &
        "header lines start with '#', blank lines are skipped on reading")
    write(unitNum, "(A)") repeat("#", 78)
    end subroutine clock_write_header

    !---------------------------------------------------------------------------
    ! SUBROUTINE: clock_read_header
    ! 读文件头(每行以 '#' 开头)，把声明的量恢复到对象上：对象名、钟模型参数、
    ! 参考时刻与计算方法。文件头中没有的项保持对象当前值不变(与外来文件兼容)。
    !
    ! 输出：
    !   methodDeclared - 文件头是否声明了计算方法
    !---------------------------------------------------------------------------
    subroutine clock_read_header(fn, this, methodDeclared)
    implicit none
    character(len=*), intent(in) :: fn
    class(type_clock), intent(inout) :: this
    logical, intent(out) :: methodDeclared
    character(len=512) :: line
    character(len=64)  :: key
    character(len=400) :: value
    integer :: unitNum, io_stat
    double precision :: x
    logical :: found

    methodDeclared = .false.

    open(newunit = unitNum, file = trim(fn), status = 'old', iostat = io_stat)
    if (io_stat /= 0) return          ! 打不开：由 clock_read_data 报错

    do
        read(unitNum, "(A)", iostat = io_stat) line
        if (io_stat /= 0) exit                     ! 文件结束
        line = adjustl(line)
        if (len_trim(line) == 0) cycle             ! 空行
        if (line(1:1) /= '#') exit                 ! 数据行：文件头结束

        call split_header_line(line, key, value)
        select case(trim(key))
        case ("object")
            this%name = trim(value)
        case ("bias0")
            call first_number(value, x, found)
            if (found) this%bias0 = x
        case ("drift")
            call first_number(value, x, found)
            if (found) this%drift = x
        case ("drift2")
            call first_number(value, x, found)
            if (found) this%drift2 = x
        case ("noise")
            call first_number(value, x, found)
            if (found) this%noise = x
        case ("t0")
            call first_number(value, x, found)
            if (found) call this%t0%init(x)
        case ("calc method")
            call first_number(value, x, found)
            if (found) then
                this%calcMethod = nint(x)
                methodDeclared = .true.
                if (this%calcMethod /= CLOCK_METHOD_FORMULA .and. &
                    this%calcMethod /= CLOCK_METHOD_INTERPOLATION) then
                    print*, "warning @clock%readData: unknown calc method in file header:", &
                            this%calcMethod, "-> formula method used"
                    this%calcMethod = CLOCK_METHOD_FORMULA
                endif
            endif
        end select
    end do
    close(unitNum)
    end subroutine clock_read_header

    !---------------------------------------------------------------------------
    ! SUBROUTINE: write_header_line
    ! 写一行文件头："# 关键字 : 值"，关键字左对齐到 12 列
    !---------------------------------------------------------------------------
    subroutine write_header_line(unitNum, key, value)
    implicit none
    integer, intent(in) :: unitNum
    character(len=*), intent(in) :: key, value
    character(len=12)  :: keyPad
    character(len=200) :: line

    ! 关键字先赋给定长 12 的变量(赋值右侧补空格 = 左对齐)，不用 A12 编辑符：
    ! A 编辑符在字段比字符串宽时是左补空格(右对齐)，位置就对不上了
    keyPad = trim(key)
    write(line, "(A,A,A,A)") "# ", keyPad, ": ", trim(value)
    write(unitNum, "(A)") trim(line)
    end subroutine write_header_line

    !---------------------------------------------------------------------------
    ! SUBROUTINE: split_header_line
    ! 把一行以 '#' 开头的文件头行拆成关键字与值："# 关键字 : 值"。
    ! 取第一个冒号处分割，故值里可以再出现冒号(如日期时间)。
    !---------------------------------------------------------------------------
    subroutine split_header_line(line, key, value)
    implicit none
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: key, value
    integer :: p

    key = ''
    value = ''
    p = index(line, ":")
    if (p <= 1) return
    key   = adjustl(line(2:p-1))
    value = adjustl(line(p+1:))
    end subroutine split_header_line

    !---------------------------------------------------------------------------
    ! SUBROUTINE: first_number
    ! 取一段文本中的第一个数值，如 "1.00000000E-04 s" -> 1.0d-4、
    ! "60676.000000000000 MJD (UTC) = 2025-01-01 ..." -> 60676.0d0。
    ! 先把数值记号(数字、正负号、小数点、指数)整体抠出来再按列表输入读，
    ! 这样数值后面跟着什么字符都不影响解析(列表输入遇到 ')' 会报错)。
    !---------------------------------------------------------------------------
    subroutine first_number(text, value, found)
    implicit none
    character(len=*), intent(in) :: text
    double precision, intent(out) :: value
    logical, intent(out) :: found
    character(len=32) :: token
    integer :: i, i1, i2, io_stat

    value = 0d0
    found = .false.
    i1 = 0
    i2 = 0
    do i = 1, len(text)
        if (is_number_char(text(i:i))) then
            if (i1 == 0) i1 = i
            i2 = i
        elseif (i1 > 0) then
            exit                              ! 数值记号到此结束
        endif
    end do
    if (i1 == 0) return
    token = text(i1:i2)
    read(token, *, iostat = io_stat) value
    found = (io_stat == 0)
    end subroutine first_number

    !---------------------------------------------------------------------------
    ! FUNCTION: is_number_char
    ! 是否为数值记号允许的字符
    !---------------------------------------------------------------------------
    logical function is_number_char(c)
    implicit none
    character, intent(in) :: c

    is_number_char = (c >= '0' .and. c <= '9') .or. c == '.' .or. &
                     c == '+' .or. c == '-' .or. &
                     c == 'e' .or. c == 'E' .or. c == 'd' .or. c == 'D'
    end function is_number_char

    !---------------------------------------------------------------------------
    ! FUNCTION: clock_effective_method
    ! 本次查询实际采用的计算方法。
    !
    ! 指定了插值方法但钟差表未定义(或记录数不足 2 条)时，给出警告并退回公式方法。
    !
    ! 注意：诊断信息一律写 stderr(unit 0)，不用 print*。查询函数经常被写在 I/O
    ! 语句里(如 print*, clk%get_clock_state(t))，此时若用 print* 输出到 unit 6，
    ! 就成了对同一文件的递归 I/O —— gfortran 运行库会在同样的单元锁上死锁。
    !---------------------------------------------------------------------------
    function clock_effective_method(this) result(method)
    implicit none
    class(type_clock), intent(in) :: this
    integer :: method

    method = this%calcMethod
    if (method == CLOCK_METHOD_INTERPOLATION) then
        if (.not. allocated(this%table_clock)) then
            write(0, *) "Warning @clock: interpolation selected but no clock table"// &
                &   " is defined; falling back to the formula method."
            write(0, *) "            call generateTable or setTable first."
            method = CLOCK_METHOD_FORMULA
        elseif (size(this%table_clock) < 2) then
            write(0, *) "Warning @clock: clock table has fewer than 2 records;"// &
                &   " falling back to the formula method."
            method = CLOCK_METHOD_FORMULA
        endif
    endif
    end function clock_effective_method

    !---------------------------------------------------------------------------
    ! SUBROUTINE: clock_locate
    ! 查询前的检查与定位：确认状态表存在、t 在表范围内，并定位 t 所在区间。
    !
    ! 输出：
    !   idx - 满足 table_clock(idx)%t <= t <= table_clock(idx+1)%t 的下标
    !   n   - 状态表记录数
    !---------------------------------------------------------------------------
    subroutine clock_locate(this, t, idx, n)
    implicit none
    class(type_clock), intent(in) :: this
    type(type_t), intent(in) :: t
    integer, intent(out) :: idx, n
    integer :: lo, hi, mid

    if (.not. allocated(this%table_clock)) then
        write(0, *) "Error @clock: the clock table has not been generated."
        write(0, *) "Call generateTable or setTable first."
        stop
    endif

    n = size(this%table_clock)
    if (n < 2) then
        write(0, *) "Error @clock: clock table has too few records (at least 2 needed)."
        stop
    endif

    if (timeMinus(t, this%table_clock(1)%t) < 0d0 .or. &
        timeMinus(t, this%table_clock(n)%t) > 0d0) then
        write(0, *) "Error @clock: requested epoch is outside the clock table range."
        write(0, *) "Requested t = ", t%val
        write(0, *) "Table range: ", this%table_clock(1)%t%val, " to ", this%table_clock(n)%t%val
        stop
    endif

    ! 二分查找(状态表按时间递增排列)
    lo = 1
    hi = n
    do while (hi - lo > 1)
        mid = (lo + hi)/2
        if (this%table_clock(mid)%t%val <= t%val) then
            lo = mid
        else
            hi = mid
        endif
    enddo
    idx = lo
    end subroutine clock_locate

    !---------------------------------------------------------------------------
    ! FUNCTION: clock_get_state
    ! 查询指定时刻的钟差
    !
    ! 按 calcMethod 选择算法：
    !   - 公式方法：直接由解析钟模型计算，任何时刻都可用；
    !   - 插值方法：查时钟状态表，与 type_satellite%getStates 的星历插值一致 ——
    !     命中表中历元时直接取该点的钟差；否则取 t 附近 4 个数据点做 4 点
    !     Lagrange 插值，表中记录少于 4 个时退化为 t 所在区间上的线性插值。
    !     插值自变量取"相对表中首历元的秒数"，避免 MJD 大数相减的精度损失。
    ! 选择插值方法但钟差表未定义时，报警告并退回公式方法。
    !
    ! 输入：
    !   t    - 查询时刻 (type_t)
    ! 返回：
    !   bias - 钟差 [s]
    !---------------------------------------------------------------------------
    function clock_get_state(this, t) result(bias)
    implicit none
    class(type_clock), intent(in) :: this
    type(type_t), intent(in) :: t
    double precision :: bias
    integer :: idx, n, j, nodes(4)
    double precision :: x(4), y(4), yy(1), dt

    ! 公式方法：解析钟模型直接给出钟差
    if (clock_effective_method(this) == CLOCK_METHOD_FORMULA) then
        bias = clock_get_model_state(this, t)
        return
    endif

    call clock_locate(this, t, idx, n)

    ! 命中表中历元，直接取该点值(1d-12 天 ≈ 0.09 微秒)
    if (dabs(timeMinus(t, this%table_clock(idx)%t)) <= 1d-12) then
        bias = this%table_clock(idx)%bias
        return
    endif
    if (dabs(timeMinus(t, this%table_clock(idx + 1)%t)) <= 1d-12) then
        bias = this%table_clock(idx + 1)%bias
        return
    endif

    ! 数据点不足 4 个，退化为线性插值
    if (n < 4) then
        dt = timeMinus(this%table_clock(idx + 1)%t, this%table_clock(idx)%t)
        bias = this%table_clock(idx)%bias + &
               (this%table_clock(idx + 1)%bias - this%table_clock(idx)%bias) * &
               (timeMinus(t, this%table_clock(idx)%t)/dt)
        return
    endif

    ! 取 4 个插值节点(边界处理与 sat_getStates 一致)
    if (idx <= 1) then
        nodes = [1, 2, 3, 4]
    elseif (idx >= n - 2) then
        nodes = [n - 3, n - 2, n - 1, n]
    else
        nodes = [idx - 1, idx, idx + 1, idx + 2]
    endif

    ! 插值自变量：相对表中首历元的秒数
    do j = 1, 4
        x(j) = timeMinus(this%table_clock(nodes(j))%t, this%table_clock(1)%t)*86400d0
        y(j) = this%table_clock(nodes(j))%bias
    enddo
    call Interpolation(4, x, y, 1, &
                       (/timeMinus(t, this%table_clock(1)%t)*86400d0/), yy)
    bias = yy(1)
    end function clock_get_state

    !---------------------------------------------------------------------------
    ! FUNCTION: clock_get_drift
    ! 查询指定时刻的钟漂移率
    !
    ! 按 calcMethod 选择算法：
    !   - 公式方法：对解析钟模型求导，v(t) = drift + 2*drift2*dt；
    !   - 插值方法：取 t 所在区间的两个数据点做线性差分
    !     v = (bias2 - bias1)/(t2 - t1)；t 恰好落在表中历元时，取该点与其后一点
    !     组成的区间(前向差分)。钟差表未定义时退回公式方法
    !
    ! 输入：
    !   t     - 查询时刻 (type_t)
    ! 返回：
    !   drift - 钟漂移率 [s/s]
    !---------------------------------------------------------------------------
    function clock_get_drift(this, t) result(drift)
    implicit none
    class(type_clock), intent(in) :: this
    type(type_t), intent(in) :: t
    double precision :: drift
    integer :: idx, n
    double precision :: dt

    ! 公式方法：对解析钟模型求导，v(t) = drift + 2*drift2*dt
    if (clock_effective_method(this) == CLOCK_METHOD_FORMULA) then
        dt = timeMinus(t, this%t0)*86400d0
        drift = this%drift + 2d0*this%drift2*dt
        return
    endif

    call clock_locate(this, t, idx, n)

    ! 时间差换算为秒
    dt = timeMinus(this%table_clock(idx + 1)%t, this%table_clock(idx)%t)*86400d0
    if (dabs(dt) < 1d-12) then
        write(0, *) "Warning @get_clock_drift: two adjacent records share the same"// &
            &   " epoch; drift rate set to 0."
        drift = 0d0
        return
    endif

    drift = (this%table_clock(idx + 1)%bias - this%table_clock(idx)%bias)/dt
    end function clock_get_drift

    !---------------------------------------------------------------------------
    ! SUBROUTINE: clock_free
    ! 释放状态表(钟模型参数保持不变)
    !---------------------------------------------------------------------------
    subroutine clock_free(this)
    implicit none
    class(type_clock), intent(inout) :: this

    if (allocated(this%table_clock)) deallocate(this%table_clock)
    end subroutine clock_free

    !---------------------------------------------------------------------------
    ! SUBROUTINE: clock_show
    ! 显示钟模型信息
    !---------------------------------------------------------------------------
    subroutine clock_show(this)
    implicit none
    class(type_clock), intent(in) :: this
    integer :: n

    print*, "--------------------------------------------------------------"
    print*, "Clock model:"
    if (len_trim(this%name) > 0) write(*, "(A,A)") "  Object          : ", trim(this%name)
    write(*, "(A,A)") "  Reference epoch : ", trim(this%t0%toDate())
    write(*, "(A,ES16.8,A)") "  Initial bias    : ", this%bias0, " s"
    write(*, "(A,ES16.8,A)") "  Drift rate      : ", this%drift, " s/s"
    write(*, "(A,ES16.8,A)") "  Drift rate 2nd  : ", this%drift2, " s/s^2"
    write(*, "(A,ES16.8,A)") "  Clock noise     : ", this%noise, " s"

    if (this%calcMethod == CLOCK_METHOD_INTERPOLATION) then
        write(*, "(A,I0,A)") "  Calc method     : ", this%calcMethod, " (interpolation)"
    else
        write(*, "(A,I0,A)") "  Calc method     : ", this%calcMethod, " (formula)"
    endif

    if (allocated(this%table_clock)) then
        n = size(this%table_clock)
        write(*, "(A,I0)") "  Table records   : ", n
        write(*, "(A,A)") "  Table start     : ", trim(this%table_clock(1)%t%toDate())
        write(*, "(A,A)") "  Table end       : ", trim(this%table_clock(n)%t%toDate())
    else
        print*, "  Table records   : 0 (not generated)"
    endif
    print*, "--------------------------------------------------------------"
    end subroutine clock_show

end module module_clock
