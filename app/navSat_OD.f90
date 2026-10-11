!================================================================================
! PROGRAM: navSat_OD
!
! 说明：
!   近地与地月空间导航卫星时空基准 —— 精密定轨程序（批处理最小二乘）。
!   读入 navSat_simulation 生成的 6 条测距观测（data/navSat/obs_*.txt），
!   联合估计 3 颗 MEO 卫星的
!     轨道状态量 6 个（笛卡尔位置速度）+ ECOM 光压参数 5 个（D0, Y0, B0, Dc, Ds）
!     + 卫星钟差 3 个（bias0, drift, drift2）
!   每颗星 14 个参数，共 42 个（type_varGroup: 轨道组 label = 5, nVar = 11；
!   钟差组 label = 13, nVar = 3）。
!
!   钟差：本次**建模并估计**。仿真里卫星钟差 |bias|*c ≈ 1 m，过去只能当作 0，
!   这 1 m 系统差被吸收进轨道（位置误差从零钟的 6 mm 恶化到 5~6 m）。
!   秩亏说明：测距只观测到 clk_R - clk_T 之差，所有钟同加常数不可观测，所以
!   必须留一个绝对时间基准 —— 本程序把**地面站钟固定为 0**（不建钟差组），
!   只估 3 颗卫星的钟。这是应用层的选择，ODor 本身不加任何约束。
!
!   力模型、星座、测站与链路与仿真程序 navSat_simulation 完全一致，只有初轨、
!   ECOM 参数与钟差系数取初值（真值 + 偏差 / 真值的一半 / clkInitScale×真值）。
!
! 输出：
!   data/navSat/odEph_MEO?.txt   定轨星历（估计值，笛卡尔 km / km/s）
!   data/navSat/odClk_MEO?.txt   估计出的卫星钟差模型（与 clock_MEO?.txt 同格式）
!   data/output/log_batchOD.txt  ODor%solve 的迭代日志（每次迭代的 RMS 与 x）
!   屏幕：估计的初轨、ECOM、钟差系数 vs 真值，以及与真值星历的位置差
!
! 注意：
!   - 必须从仓库根目录运行（数据路径都是相对路径）。
!   - 观测文件里的历元是接收设备的钟面时间；本场景钟差只有 ~1e-9 s，
!     落盘精度（微秒）下与理想格点无法区分，故不构成历元越界问题。
!
! 作者：Shawn Hu
! 日期：2026-10-10
!================================================================================

program navSat_OD

    use module_Atlas
    use module_const
    use module_time
    use module_timeSeries
    use module_elements

    implicit none

    !---------------------- 场景设置（与 navSat_simulation 保持一致） --------------
    real(8), parameter :: step   = 30d0          ! 采样间隔 [s]
    real(8), parameter :: arcDay = 1d0           ! 弧长 [天]

    integer, parameter :: nSat = 3
    real(8), parameter :: aKm    = 27906.1d0     ! 半长轴 [km]
    real(8), parameter :: ecc    = 0.001d0       ! 偏心率
    real(8), parameter :: incDeg = 55d0          ! 倾角 [deg]
    real(8), parameter :: raanDeg(nSat) = (/0d0, 120d0, 240d0/)
    real(8), parameter :: maDeg(nSat)   = (/0d0, 120d0, 240d0/)
    character(len=8), parameter :: satName(nSat) = (/"MEO1", "MEO2", "MEO3"/)
    real(8), parameter :: satMass      = 1000d0  ! [kg]
    real(8), parameter :: amrSRP       = 0.02d0  ! [m^2/kg]（ECOM 模式下不参与计算）
    real(8), parameter :: reflectivity = 0.5d0   ! 仅 cannonball 模式使用
    integer, parameter :: gravityDegree = 12     ! 引力场阶数
    integer, parameter :: gravityOrder  = 12     ! 引力场次数

    ! ECOM 真值（与仿真一致），顺序 [D0, Y0, B0, Bc, Bs, Dc, Ds, Yc, Ys]
    real(8), parameter :: ecomTrue(9) = (/ 1.4d-8,  3.0d-9, -1.5d-9, &
        &                                   0d0,     0d0,           &
        &                                   2.0d-9, -2.0d-9,        &
        &                                   0d0,     0d0 /)
    integer, parameter :: nECOM = 5              ! 待估的 ECOM 参数个数（前 5 个）
    ! 是否把 ECOM 光压参数一并作为待估参数(label 5)。
    ! .false. 时只估轨道状态量(label 1, 每星 6 个, 共 18 个), ECOM 固定不估,
    ! 取 ecomScale*真值 —— 相当于"光压模型已知(或已知其误差有多大)"的情形。
    logical, parameter :: estEcom = .true.
    real(8), parameter :: ecomScale = 1.0d0      ! ECOM 固定值 = 真值 * ecomScale

    ! 卫星钟真值（与 navSat_simulation 一致）：bias(t) = bias0 + drift*dt + drift2*dt^2
    real(8), parameter :: clkTrue0(nSat)  = (/ 1.2d-9,  -2.8d-9,  1.0d-9  /)  ! [s]
    real(8), parameter :: clkTrue1(nSat)  = (/ 1.5d-14,  8.0d-15, 2.5d-14 /)  ! [s/s]
    real(8), parameter :: clkTrue2(nSat)  = (/ 6.0d-20, -1.0d-19, 2.0d-20 /)  ! [s/s^2]
    ! 是否把卫星钟差作为待估参数（label 13, nVar = 3，每星 3 个共 9 个）。
    ! 地面站钟固定为 0 作时间基准 —— 测距只观测钟差之差，必须留一个绝对基准。
    logical, parameter :: estClock = .true.
    ! 钟差初值 = clkInitScale × 真值。0 = "对钟差一无所知"（默认，最难）；
    ! 1 = 从真值出发（用于核对解是不动点）。
    real(8), parameter :: clkInitScale = 0d0

    ! 初轨偏差：真值 + 偏差（笛卡尔，SI 单位）
    real(8), parameter :: pertScale = 1d0    ! 偏差缩放：1 = 正常；0 = 从真值出发
    logical, parameter :: checkDiag = .false. ! 是否跑 STM / calB 校核诊断
    logical, parameter :: useISL = .true.   ! 是否用星间链路（诊断用；.false. 只留星地）
    logical, parameter :: showProg = .true. ! 迭代进度条（ODor%showProgress）
    real(8), parameter :: posPert(3) = (/ 10d0, -10d0,  10d0 /)   ! [m]
    real(8), parameter :: velPert(3) = (/ 0.002d0, -0.002d0, 0.002d0 /) ! [m/s]

    character(len=16), parameter :: staName = "Beijing"
    real(8), parameter :: staLon = 116.4d0, staLat = 39.9d0, staAlt = 50d0
    real(8), parameter :: minEle = 10d0

    character(len=*), parameter :: dir = "data/navSat/"

    !---------------------- 变量声明 ----------------------
    type(type_t) :: t0, te
    type(type_satelliteScenario) :: ss
    type(type_satellite) :: sat
    type(type_station)   :: sta
    type(type_ODor)      :: ODor
    type(type_varGroup)  :: Groups(2 * nSat)   ! 前 nSat 个轨道组, 后 nSat 个钟差组
    type(type_t), allocatable :: epochs(:)

    type(type_satellite) :: truth(nSat)   ! 真值星历（读 eph_MEO?.txt）
    real(8) :: kpl(6), cts(6), pvEst(6), pvTru(6), dr(3)
    real(8) :: truth0(6, nSat)   ! 各星真值初轨（加偏差前），仅用于打印对比
    real(8) :: sum2, maxErr, rms
    integer :: i, j, k, nAccess, nVarPerSat, nGroup, nClk
    character(len=200) :: fname, intro

    nVarPerSat = 6
    if (estEcom) nVarPerSat = 6 + nECOM

    nClk = 0
    if (estClock) nClk = 3                      ! bias0 + drift + drift2
    nGroup = nSat
    if (estClock) nGroup = nSat + nSat          ! 轨道组 + 钟差组

    !====================================================================
    ! 1. 时间与场景
    !====================================================================
    call t0%init(2025d0, 1d0, 1d0, 0d0, 0d0, 0d0)
    te = timeAdd(t0, arcDay)

    if (estClock) then
        intro = "navSat POD: 3 MEO satellites, states + ECOM + satellite clocks"
    elseif (estEcom) then
        intro = "navSat POD: 3 MEO satellites, states + 5 ECOM parameters, clocks fixed at 0"
    else
        intro = "navSat POD: 3 MEO satellites, states only, clocks fixed at 0"
    endif
    call ss%init(t0, te, step, eleType = STATE_TYPE_CARTESIAN, introduction = trim(intro))

    allocate(epochs(ss%ts%n))
    do i = 1, ss%ts%n
        epochs(i) = ss%ts%time(i)
    enddo

    !====================================================================
    ! 2. 3 颗 MEO 卫星：力模型与仿真一致，初轨与 ECOM 取初值
    !    钟模型不设置 —— 缺省的钟 bias/drift 全为 0，即"钟差未知、当作 0"
    !====================================================================
    do i = 1, nSat
        kpl(1) = aKm * u_km
        kpl(2) = ecc
        kpl(3) = incDeg*deg
        kpl(4) = raanDeg(i)*deg
        kpl(5) = 0d0
        kpl(6) = maDeg(i)*deg
        call kpl2cts(kpl, cts)
        truth0(:, i) = cts              ! 记下本星真值初轨
        ! 加上初值偏差（SI -> 归一化：位置 /(LU*1000)，速度 /(LU*1000)*TU）
        cts(1:3) = cts(1:3) + pertScale*posPert/(LENGTH_UNIT*1000d0)
        cts(4:6) = cts(4:6) + pertScale*velPert/(LENGTH_UNIT*1000d0)*TIME_UNIT

        call sat%init(t0, cts, &
            &   number = i, &
            &   name = trim(satName(i)), &
            &   mass = satMass, &
            &   eleType = STATE_TYPE_CARTESIAN, &
            &   propagator = 'numerical', &
            &   degree = gravityDegree, &
            &   order = gravityOrder, &
            &   reflectivity = reflectivity, &
            &   amrSRP = amrSRP)
        sat%opMod%usePV          = 1   ! 积分用笛卡尔位置速度
        sat%opMod%key_moonGravity = 1
        sat%opMod%key_sunGravity  = 1
        sat%opMod%key_srp         = 1
        sat%opMod%SRP_Model       = 2   ! ECOM
        sat%opMod%ECOM            = ecomScale*ecomTrue   ! 初值（前 5 个待估，其余恒 0）

        ! 卫星钟：初值取 clkInitScale×真值（默认 0 = 对钟差一无所知）。
        ! 基准历元必须显式给出 —— 缺省构造的 type_clock 的 t0 是 MJD 0，
        ! 直接赋系数而不设 t0 会让 drift 项乘上 5.24e9 s。
        call sat%clock%setModel(t0, clkInitScale*clkTrue0(i), &
            &   clkInitScale*clkTrue1(i), clkInitScale*clkTrue2(i))

        call ss%addSatellite(sat)
    enddo

    !====================================================================
    ! 3. 地面站：站钟固定为 0，作为测距的绝对时间基准（不建钟差组）
    !====================================================================
    call sta%init(number = 1, name = staName, &
        &   longitude = staLon, latitude = staLat, altitude = staAlt)
    call sta%clock%setModel(t0, 0d0, 0d0, 0d0)
    call ss%addStation(sta)

    !====================================================================
    ! 4. 观测链路：与仿真同样的 6 条，数据从 obs_*.txt 读入
    !    （readData 按文件头恢复 obsType=13 与收发设备名）
    !====================================================================
    nAccess = 0
    if (useISL) then
        do i = 1, nSat - 1
            do j = i + 1, nSat
                nAccess = nAccess + 1
                call addReadAccess(ss, nAccess, trim(satName(i))//"-"//trim(satName(j)), &
                    &   i, j, .true.)
            enddo
        enddo
    endif
    do i = 1, nSat
        nAccess = nAccess + 1
        call addReadAccess(ss, nAccess, trim(satName(i))//"-"//trim(staName), i, 1, .false.)
    enddo
    call ss%reportAccess(.true.)

    ! 测站星历：od 求解时要用 station%getPV 算几何距离
    call ss%station(1)%generateTrajectory(ss%ts)

    !====================================================================
    ! 5. 待估参数分组
    !    前 nSat 组：轨道，每组 = 6 状态量 + nECOM 个 ECOM 参数
    !    后 nSat 组：卫星钟差（label 13, nVar 3），地面站钟不估（作时间基准）
    !    idx 是各组在全局 x 向量里的下标
    !====================================================================
    do i = 1, nSat
        Groups(i)%nVar   = nVarPerSat
        Groups(i)%label  = 1               ! 1 = 只估状态量; 见 estEcom
        if (estEcom) Groups(i)%label = 5
        Groups(i)%origin = i
        allocate(Groups(i)%idx(nVarPerSat))
        Groups(i)%idx = (/ (j, j = (i-1)*nVarPerSat + 1, i*nVarPerSat) /)
    enddo

    if (estClock) then
        do i = 1, nSat
            j = nSat + i
            Groups(j)%nVar   = nClk        ! 3: bias0, drift, drift2
            Groups(j)%label  = 10 + nClk   ! 13
            Groups(j)%origin = i
            Groups(j)%objType = OBJECT_TYPE_SATELLITE
            allocate(Groups(j)%idx(nClk))
            Groups(j)%idx = (/ (k, k = nSat*nVarPerSat + (i-1)*nClk + 1, &
                &                   nSat*nVarPerSat + i*nClk) /)
        enddo
    endif

    !====================================================================
    ! 6. 定轨
    !====================================================================
    call ODor%init(ss, ss%ts%t0, nGroup, Groups(1:nGroup), OD_METHOD_BATCH, intro = trim(intro))
    ODor%showProgress = showProg      ! 迭代进度条: 轨迹阶段 + 观测遍历阶段
    call ODor%solve(ss)

    !====================================================================
    ! 7. 输出：估计星历 + 与真值对比
    !====================================================================
    do i = 1, nSat
        ! 用估计出的状态重新生成星历后再存（solve 最后一次迭代里的星历还是
        ! 更新前那一步的）
        call ss%sat(i)%generateTrajectory(epochs)
        write(fname, "(A,A,A)") trim(dir)//"odEph_", trim(satName(i)), ".txt"
        call ss%sat(i)%savetrajectory(STATE_TYPE_CARTESIAN, trim(fname), freeEph = .false.)
        print "(A)", "  saved: "//trim(fname)
    enddo

    ! 估计参数 vs 真值
    print *, ""
    print *, "================ estimated vs truth ================"
    do i = 1, nSat
        call cts_normal2SI(ss%sat(i)%states%elements, cts)
        call cts_normal2SI(truth0(:, i), pvTru)   ! 本星真值初轨 -> km, km/s
        print "(A)", " satellite "//trim(satName(i))
        print "(A,3F16.3)", "   estimated r0 (km): ", cts(1:3)
        print "(A,3F16.3)", "   truth     r0 (km): ", pvTru(1:3)
        ! 误差用米: 定轨精度到 mm 级时, 以 km 记的 6 位小数只会印出 0.000000
        print "(A,3ES15.6)", "   error     dr (m):  ", (cts(1:3) - pvTru(1:3))*1000d0
        if (estEcom) then
            print "(A,5(ES13.5,1x))", "   estimated ECOM(1:5) = ", ss%sat(i)%opMod%ECOM(1:nECOM)
        else
            print "(A,5(ES13.5,1x))", "   fixed     ECOM(1:5) = ", ss%sat(i)%opMod%ECOM(1:nECOM)
        endif
        print "(A,5(ES13.5,1x))", "   truth     ECOM(1:5) = ", ecomTrue(1:nECOM)
        if (estClock) then
            print "(A,3ES15.6)", "   estimated clk = ", ss%sat(i)%clock%bias0, &
                &   ss%sat(i)%clock%drift, ss%sat(i)%clock%drift2
            print "(A,3ES15.6)", "   truth     clk = ", clkTrue0(i), clkTrue1(i), clkTrue2(i)
            print "(A,3ES15.6)", "   error     clk = ", ss%sat(i)%clock%bias0 - clkTrue0(i), &
                &   ss%sat(i)%clock%drift - clkTrue1(i), ss%sat(i)%clock%drift2 - clkTrue2(i)
        endif
    enddo

    ! 估计出的钟差模型落盘（与仿真 clock_MEO?.txt 同格式，便于对照）
    if (estClock) then
        print *, ""
        do i = 1, nSat
            write(fname, "(A,A,A)") trim(dir)//"odClk_", trim(satName(i)), ".txt"
            call ss%sat(i)%clock%saveData(trim(fname), name = trim(satName(i)))
            print "(A)", "  saved: "//trim(fname)
        enddo
    endif

    ! 与真值星历逐历元比较位置差
    do i = 1, nSat
        write(fname, "(A,A,A)") trim(dir)//"eph_", trim(satName(i)), ".txt"
        call truth(i)%readtrajectory(STATE_TYPE_CARTESIAN, trim(fname))
    enddo

    print *, ""
    print *, "================ position error vs truth ephemeris ================"
    do i = 1, nSat
        sum2 = 0d0; maxErr = 0d0
        do k = 1, ss%ts%n
            pvEst = ss%sat(i)%getPV(epochs(k))
            pvTru = truth(i)%getPV(epochs(k))
            dr = (pvEst(1:3) - pvTru(1:3))*LENGTH_UNIT*1000d0   ! 归一化 -> m
            sum2 = sum2 + dot_product(dr, dr)
            maxErr = max(maxErr, norm2(dr))
        enddo
        rms = dsqrt(sum2/dble(ss%ts%n))
        print "(A,F10.3,A,F10.3,A)", " "//trim(satName(i))//":  RMS = ", rms, " m   max = ", &
            &   maxErr, " m"
    enddo

    !====================================================================
    ! 7.5 诊断：校核定轨所用的 STM 与 calB（临时）
    !====================================================================
    if (checkDiag) then
        call checkSTM(ss, ODor)
        call checkB(ss, ODor, 1)     ! 星间链路 MEO1-MEO2
        call checkB(ss, ODor, 5)     ! 星地链路 MEO2-Beijing (toNum=1, 与卫星 1 编号相撞)
    endif

    !====================================================================
    ! 8. 清理
    !====================================================================
    do i = 1, nGroup
        deallocate(Groups(i)%idx)
    enddo
    deallocate(epochs)
    call ss%freeSS

contains

!==============================================================================
! 诊断（临时）：把定轨实际使用的 STM 序列与"全模型有限差分"逐列对比。
! 差分用同一个卫星的完整力模型（sat%propagate），因此若两者不一致，
! 差异只能来自 STM 序列本身（快速预报/插值），而不是力模型。
!==============================================================================
subroutine checkSTM(ss, ODor)
    type(type_satelliteScenario), intent(in) :: ss
    type(type_ODor), intent(in) :: ODor
    type(type_satellite) :: satA, satB
    type(type_t) :: t1
    double precision :: x0(6), x1(6), xb(6), FD(6, 6), M(6, 6), dlt
    double precision :: xOne(6), xSeq(6), xr(6), tSeq
    type(type_fastOption) :: opx
    integer :: i, j

    t1 = timeAdd(t0, 0.25d0)                 ! +6 h
    x0 = ss%sat(1)%states%elements
    print *, ""
    print *, "==== check: STM series vs full-model finite difference (sat 1, t0 -> +6h) ===="
    write(*, "(A,SP,3ES14.5)") "  x0(1:3) = ", x0(1:3)
    write(*, "(A,SP,3ES14.5)") "  x0(4:6) = ", x0(4:6)

    satA = ss%sat(1)
    call satA%propagate(t1)
    x1 = satA%states%elements

    call ODor%STMseries(1)%get(t1, M)

    dlt = 1d-6                               ! 归一化: 位置 6.4 m, 速度 0.79 m/s
    do j = 1, 6
        satB = ss%sat(1)
        xb = x0
        xb(j) = xb(j) + dlt
        satB%states%elements = xb
        call satB%propagate(t1)
        FD(1:6, j) = (satB%states%elements - x1)/dlt
        write(*, "(A,I2,A,ES11.3,A,ES11.3,A,ES10.2,A,ES10.2)") "  col ", j, &
            &   "  |FD|=", maxval(dabs(FD(:, j))), "  |STM|=", maxval(dabs(M(:, j))), &
            &   "  |FD|/|STM|=", maxval(dabs(FD(:, j)))/maxval(dabs(M(:, j)))
    enddo

    do j = 1, 6
        write(*, "(A,I2)") "  ---- col ", j, " (rows 1..6: FD / STM) ----"
        do i = 1, 6
            if (dabs(M(i, j)) > 1d-20) then
                write(*, "(I4,2X,2ES16.6,A,F10.4)") i, FD(i, j), M(i, j), &
                    &   "  ratio=", FD(i, j)/M(i, j)
            else
                write(*, "(I4,2X,2ES16.6,A)") i, FD(i, j), M(i, j), "  ratio= --"
            endif
        enddo
    enddo

    ! STM 序列现在的做法是"逐节点顺序推进"(每步 stepSec)。这里直接核对它是否
    ! 等价于"从基点一次积分到同一时刻" —— 两种走法若一致, 新写法的轨迹就是精确的。
    opx%useJ2 = .true.
    opx%useJ22 = .true.
    opx%useMoon = (ss%sat(1)%opMod%key_moonGravity == 1)
    opx%useSun = (ss%sat(1)%opMod%key_sunGravity == 1)
    opx%useSRP = (ss%sat(1)%opMod%key_srp == 1)
    opx%SRP_Model = ss%sat(1)%opMod%SRP_Model
    opx%reflectivity = ss%sat(1)%opMod%reflectivity
    opx%amrSRP = ss%sat(1)%opMod%amrSRP
    opx%ECOM = ss%sat(1)%opMod%ECOM

    xOne = x0
    call fastPrediction(x0, ss%sat(1)%states%t%val, &
        &   timeMinus(t1, ss%sat(1)%states%t)*86400d0, opx, xOne)
    xSeq = x0
    tSeq = ss%sat(1)%states%t%val
    do i = 1, 36                                   ! 36 × 600 s = 21600 s = t1
        xr = xSeq
        call fastPrediction(xr, tSeq, 600d0, opx, xSeq)
        tSeq = tSeq + 600d0/86400d0
    enddo
    write(*, "(A,ES12.4,A,ES12.4)") "  step-by-step vs one-shot: |dr| = ", &
        &   norm2(xSeq(1:3) - xOne(1:3))*LENGTH_UNIT*1000d0, " m,  |dv| = ", &
        &   norm2(xSeq(4:6) - xOne(4:6))*LENGTH_UNIT*1000d0/TIME_UNIT, " m/s"
end subroutine checkSTM

!==============================================================================
! 诊断（临时）：端到端校核 calB —— 对第 1 条链路的第 1 个历元，把 calB 给出的
! 一行偏导与"扰动卫星 1 的初值、重新生成轨迹、重算模型观测"得到的差分对比。
! 这一行 B 覆盖了 calB 的全部环节: STM 序列 + DRDX + DYDR + 列索引。
!==============================================================================
subroutine checkB(ss, ODor, iac)
    type(type_satelliteScenario), intent(inout) :: ss
    type(type_ODor), intent(inout) :: ODor
    integer, intent(in) :: iac
    type(type_access) :: accC
    type(type_timeSeries) :: tso
    double precision, allocatable :: B(:, :)
    double precision :: c0, c1, col(6), dlt
    integer :: it, j, isat, j0, iend, iendType

    allocate(B(1, ODor%nVar))          ! 按当前参数个数分配(ECOM 模式下不止 18)
    tso%time = ss%access(iac)%timeTable
    tso%n = size(ss%access(iac)%timeTable)
    it = tso%n/2 + 1                              ! 取弧段中间的历元: 起点处 ∂r/∂v0 = 0,
                                                  ! 速度列的差分只是舍入噪声

    print *, ""
    write(*, "(A,I2,A,2I3,A,F12.3,A)") " ==== check: calB vs full-model FD (access ", iac, &
        &   "  from/to: ", ss%access(iac)%fromNum, ss%access(iac)%toNum, &
        &   "  epoch sec = ", tso%time(it)%sec, " )"

    call calB(ODor, ss, ss%access(iac), 1, ODor%nVar, it, B)

    ! 链路两端分别校核(星间链路时两端都是卫星)
    do iend = 1, 2
        if (iend == 1) then
            iendType = ss%access(iac)%fromType
            isat = ss%access(iac)%fromNum
        else
            iendType = ss%access(iac)%toType
            isat = ss%access(iac)%toNum
        endif
        if (iendType /= OBJECT_TYPE_SATELLITE) cycle
        j0 = (isat - 1)*nVarPerSat                ! 该星状态在全局 x 中的起始列
                                                  ! (每星 nVarPerSat 个: 6 状态 + ECOM)

        accC = ss%access(iac)
        call observationAccess(ss, accC, tso, .true.)
        c0 = accC%dataTable(it, 1)

        dlt = 1d-6
        do j = 1, 6
            ss%sat(isat)%states%elements(j) = ss%sat(isat)%states%elements(j) + dlt
            call ss%sat(isat)%generateTrajectory(ss%ts%time)
            accC = ss%access(iac)
            call observationAccess(ss, accC, tso, .true.)
            c1 = accC%dataTable(it, 1)
            col(j) = (c1 - c0)/dlt
            ss%sat(isat)%states%elements(j) = ss%sat(isat)%states%elements(j) - dlt
        enddo
        call ss%sat(isat)%generateTrajectory(ss%ts%time)      ! 还原

        if (iend == 1) then
            write(*, "(A,I2,A)") "  -- transmitter sat ", isat, " --"
        else
            write(*, "(A,I2,A)") "  -- receiver sat ", isat, " --"
        endif
        do j = 1, 6
            write(*, "(A,I2,A,ES14.6,A,ES14.6,A,F9.4)") "    col ", j, &
                &   "  FD = ", col(j), "  calB = ", B(1, j0 + j), &
                &   "  ratio = ", col(j)/B(1, j0 + j)
        enddo

        ! ---- 钟差列: 直接项, 不经过轨道 ----
        ! 观测模型对钟差系数是严格线性的, 所以差分与解析值应当一位不差。
        ! 这一步同时验证符号(收/发)、dt 基准(this%t)、幂次与列映射。
        if (estClock) call checkB_clock(ss, ODor, iac, it, isat, iend, B, tso)
    enddo
    deallocate(B)
end subroutine checkB

!==============================================================================
! 诊断（临时）：钟差列的有限差分校核。
! 扰动该星钟的第 j 个系数 -> 重算模型观测 -> 与 calB 的钟差列比较。
!==============================================================================
subroutine checkB_clock(ss, ODor, iac, it, isat, iend, B, tso)
    type(type_satelliteScenario), intent(inout) :: ss
    type(type_ODor), intent(inout) :: ODor
    type(type_access) :: accC
    type(type_timeSeries), intent(in) :: tso
    double precision, intent(in) :: B(:, :)
    integer, intent(in) :: iac, it, isat, iend
    double precision :: c0, c1, col, dlt, save(3)
    double precision, parameter :: dltArr(3) = (/ 1d-11, 1d-15, 1d-19 /)
    integer :: j, jc, iclk

    ! 找该星的钟差组, 记下它的全局列下标
    iclk = 0
    do j = nSat + 1, nGroup
        if (ODor%group(j)%objType == OBJECT_TYPE_SATELLITE .and. &
            &   ODor%group(j)%origin == isat) iclk = j
    enddo
    if (iclk == 0) return

    accC = ss%access(iac)
    call observationAccess(ss, accC, tso, .true.)
    c0 = accC%dataTable(it, 1)

    save(1) = ss%sat(isat)%clock%bias0
    save(2) = ss%sat(isat)%clock%drift
    save(3) = ss%sat(isat)%clock%drift2

    if (iend == 1) then
        write(*, "(A,I2,A)") "  -- clock column (transmitter sat ", isat, ") --"
    else
        write(*, "(A,I2,A)") "  -- clock column (receiver sat ", isat, ") --"
    endif
    do j = 1, ODor%group(iclk)%nVar
        dlt = dltArr(j)                           ! 1e-11 / 1e-15 / 1e-19
        ! 三种步长给出的距离变化都约 4e-10(归一化), 远高于舍入, 便于逐列比对
        select case (j)
        case (1); ss%sat(isat)%clock%bias0  = save(1) + dlt
        case (2); ss%sat(isat)%clock%drift  = save(2) + dlt
        case (3); ss%sat(isat)%clock%drift2 = save(3) + dlt
        endselect
        accC = ss%access(iac)
        call observationAccess(ss, accC, tso, .true.)
        c1 = accC%dataTable(it, 1)
        col = (c1 - c0)/dlt
        ss%sat(isat)%clock%bias0  = save(1)
        ss%sat(isat)%clock%drift  = save(2)
        ss%sat(isat)%clock%drift2 = save(3)

        jc = ODor%group(iclk)%idx(j)
        write(*, "(A,I2,A,ES14.6,A,ES14.6,A,F9.4)") "    clk ", j, &
            &   "  FD = ", col, "  calB = ", B(1, jc), "  ratio = ", col/B(1, jc)
    enddo
end subroutine checkB_clock


!==============================================================================
! 建立一条"从文件读数"的测距 access 并加入场景。
! access 在过程内声明，每次调用都从类型默认值开始；约束按链路类型设置，
! 再用 readData 读入 obs_<链路名>.txt（文件头会覆盖 obsType 与设备名）。
!==============================================================================
subroutine addReadAccess(ss, num, aname, fromNum, toNum, isSatSat)
    type(type_satelliteScenario), intent(inout) :: ss
    integer, intent(in) :: num, fromNum, toNum
    character(len=*), intent(in) :: aname
    logical, intent(in) :: isSatSat
    type(type_access) :: ac
    character(len=200) :: fn

    if (isSatSat) then
        call ac%init(number = num, name = aname, &
            &   fromType = OBJECT_TYPE_SATELLITE, fromNum = fromNum, &
            &   toType = OBJECT_TYPE_SATELLITE, toNum = toNum, &
            &   obsType = OBS_TYPE_RANGING)
        ac%EarthBlockage = .true.
    else
        call ac%init(number = num, name = aname, &
            &   fromType = OBJECT_TYPE_SATELLITE, fromNum = fromNum, &
            &   toType = OBJECT_TYPE_STATION, toNum = toNum, &
            &   obsType = OBS_TYPE_RANGING, minEle = minEle)
        ac%shadow = .false.    ! 与仿真一致：测距不看地影
    endif

    write(fn, "(A,A,A)") trim(dir)//"obs_", trim(aname), ".txt"
    call ac%readData(trim(fn))
    print "(A,A,A,I0,A)", "  read ", trim(fn), "   epochs = ", ac%nData, ""

    call ss%addAccess(ac)
end subroutine addReadAccess

end program navSat_OD
