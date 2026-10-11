!================================================================================
! PROGRAM: navSat_simulation
!
! 说明：
!   近地与地月空间导航卫星时空基准 —— 数据仿真程序。
!   建立"3 颗 MEO 卫星 + 1 个地面站"的场景，两两建立 access（3 条星间链路
!   + 3 条星地链路，共 6 条），生成测距(RANGING)观测数据。卫星钟差非零，
!   三颗星各用一组不同的钟模型参数；地面站钟取零，作时间基准。
!
!   卫星状态量全部用位置速度（笛卡尔）：eleType = STATE_TYPE_CARTESIAN，
!   积分也在笛卡尔形式下进行（opMod%usePV = 1）。
!
!   本程序只做数据仿真，不涉及定轨。输出全部在 data/navSat/ 下：
!     eph_MEO?.txt          卫星星历（真值，笛卡尔 km / km/s）
!     obs_MEO?-MEO?.txt     星间测距观测
!     obs_MEO?-Beijing.txt  星地测距观测
!     clock_MEO?.txt        卫星钟解析模型（真值，公式法，无钟差表）
!
!   观测模型（与 module_access 中 calObsData1 / calObsData2 的 RANGING 一致）：
!     rho_obs = rho_true + (clk_R - clk_T)*c
!     t_obs   = t_true + clk_R
!   即"接收端计时"：测距值含收发钟差之差，观测历元按接收设备的钟面时间标记。
!   地面站钟为零，故星地观测的历元不偏移，测距值只含卫星钟差项；星间观测的
!   历元按接收星的钟面时间标记，偏移量即该星在收发时刻的钟差。
!
! 注意：
!   - 必须从仓库根目录运行（数据路径都是相对路径）。
!   - data/navSat/ 目录需事先存在（本程序不建目录）。
!
! 作者：Shawn Hu
! 日期：2026-10-10
!================================================================================

program navSat_simulation

    use module_Atlas        ! 公共 API（module_ODor / module_IODor / module_sp3）
    use module_const
    use module_time
    use module_timeSeries
    use module_elements
    use module_clock

    implicit none

    !---------------------- 场景设置（硬编码，改后重新编译） ----------------------
    ! 时间
    real(8), parameter :: step   = 30d0          ! 采样间隔 [s]
    real(8), parameter :: arcDay = 1d0           ! 仿真弧长 [天]

    ! 星座：Walker delta 3/3/1，MEO 近圆轨道（a、i 取北斗 MEO 量级）
    integer, parameter :: nSat = 3
    real(8), parameter :: aKm    = 27906.1d0     ! 半长轴 [km]（高度约 21528 km）
    real(8), parameter :: ecc    = 0.001d0       ! 偏心率
    real(8), parameter :: incDeg = 55d0          ! 轨道倾角 [deg]
    real(8), parameter :: raanDeg(nSat) = (/0d0, 120d0, 240d0/)  ! 升交点赤经 [deg]
    real(8), parameter :: maDeg(nSat)   = (/0d0, 120d0, 240d0/)  ! 平近点角 [deg]
    character(len=8), parameter :: satName(nSat) = (/"MEO1", "MEO2", "MEO3"/)
    real(8), parameter :: satMass      = 1000d0  ! 卫星质量 [kg]
    real(8), parameter :: amrSRP       = 0.02d0  ! 光压面积质量比 [m^2/kg]（ECOM 模式下不参与计算）
    real(8), parameter :: reflectivity = 0.5d0   ! 表面反射系数（仅 cannonball 模式使用）
    integer, parameter :: gravityDegree = 12     ! 地球引力场阶数（n）
    integer, parameter :: gravityOrder  = 12     ! 地球引力场次数（m），此处取满阶

    ! 太阳光压：ECOM 模型（SRP_Model = 2），9 参数框架里只启用 5 个，其余为 0。
    ! 读取顺序见 type_opMod%ECOM（由 module_perturbationSRP 定义）：
    !   1 D0  2 Y0  3 B0  4 Bc  5 Bs  6 Dc  7 Ds  8 Yc  9 Ys
    !   F = (D0 + Dc*cos u + Ds*sin u)*e_D
    !     + (Y0 + Yc*cos u + Ys*sin u)*e_Y
    !     + (B0 + Bc*cos u + Bs*sin u)*e_B
    ! 单位：归一化加速度（1 = ACC_UNIT = 9.79828762253515 m/s^2）。
    ! e_D 由太阳指向卫星，D0 > 0 即把卫星推离太阳（光压的物理方向）。
    ! 量级参照同参数的 cannonball：(1+reflectivity)*amrSRP*amr2k*日地距修正
    ! ≈ 1.44e-8（归一化）≈ 1.41e-7 m/s^2，故 D0 取同量级，周期项与 Y0/B0 取 10%~30%。
    ! 本组参数下 |F| 在一个轨道周期内在 1.16e-8 ~ 1.71e-8 之间变化，均值≈D0。
    real(8), parameter :: ecomParams(9) = (/ 1.4d-8,  3.0d-9, -1.5d-9, &   ! D0, Y0, B0
        &                                     0d0,     0d0,           &   ! Bc, Bs
        &                                     2.0d-9, -2.0d-9,        &   ! Dc, Ds
        &                                     0d0,     0d0 /)          ! Yc, Ys

    ! 钟模型：3 颗星各不相同（初始钟差、钟漂、二阶钟漂）
    ! 真值钟差由解析模型直接计算（公式法，calcMethod = CLOCK_METHOD_FORMULA）：
    !   bias(t) = bias0 + drift*dt + drift2*dt^2 ,  dt = t - t0 [s]
    ! 量级刻意压到 |bias|*c <= ~1 m：定轨程序不估计钟差（钟差未知，只能当 0 用），
    ! 这 1 m 就是钟差给测距带来的系统性误差。按 1 m / c = 3.34e-9 s 分配三个
    ! 系数在一整天末尾的贡献：bias0 ~ 1e-9 s、drift*86400 ~ 2e-9 s、
    ! drift2*86400^2 ~ 1e-9 s（=> drift2 ~ 1e-19 s/s^2），三者之和 < 3.3e-9 s。
    ! 钟差总开关（验证用）：1d0 = 按下面的钟模型产生钟差；0d0 = 钟差全零。
    ! 置 0 后观测里不再含收发钟差，用来检验"这 ~1 m 钟差是否就是定轨误差的来源"。
    real(8), parameter :: clkScale = 1d0
    real(8), parameter :: clkBias0(nSat)  = (/ 1.2d-9, -2.8d-9,  1.0d-9 /)   ! [s]
    real(8), parameter :: clkDrift(nSat)  = (/ 1.5d-14, 8.0d-15, 2.5d-14 /)  ! [s/s]
    real(8), parameter :: clkDrift2(nSat) = (/ 6.0d-20, -1.0d-19, 2.0d-20 /) ! [s/s^2]

    ! 地面站：北京。站钟为零，作时间基准
    character(len=16), parameter :: staName = "Beijing"
    real(8), parameter :: staLon = 116.4d0       ! 经度 [deg]
    real(8), parameter :: staLat = 39.9d0        ! 纬度 [deg]
    real(8), parameter :: staAlt = 50d0          ! 高程 [m]

    ! 链路约束
    real(8), parameter :: minEle = 10d0          ! 星地链路最低仰角 [deg]

    ! 输出目录（须已存在）
    character(len=*), parameter :: outDir = "data/navSat/"

    !---------------------- 变量声明 ----------------------
    type(type_t) :: t0, te
    type(type_satelliteScenario) :: ss
    type(type_satellite) :: sat
    type(type_station)   :: sta
    type(type_t), allocatable :: epochs(:)

    real(8) :: kpl(6), cts(6)
    integer :: i, j, k, nAccess
    character(len=200) :: fname, intro

    !====================================================================
    ! 1. 时间与场景
    !====================================================================
    call t0%init(2025d0, 1d0, 1d0, 0d0, 0d0, 0d0)   ! 2025-01-01 00:00:00 UTC
    te = timeAdd(t0, arcDay)

    intro = "navSat: 3 MEO satellites + 1 ground station, pairwise ranging links"
    call ss%init(t0, te, step, eleType = STATE_TYPE_CARTESIAN, introduction = trim(intro))

    ! 场景历元数组：卫星星历与观测量都在这些历元上取值
    allocate(epochs(ss%ts%n))
    do i = 1, ss%ts%n
        epochs(i) = ss%ts%time(i)
    enddo

    !====================================================================
    ! 2. 3 颗 MEO 卫星（Walker 3/3/1）
    !    状态量：位置速度；积分：笛卡尔；力模型：非球形引力(12 阶次) +
    !    月球三体 + 太阳三体 + 太阳光压(ECOM 模型，启用 5 个参数)
    !    （MEO 高度大气阻力可忽略，不开 key_drag）
    !====================================================================
    do i = 1, nSat
        kpl(1) = aKm * u_km         ! 半长轴（归一化：地球半径）
        kpl(2) = ecc                ! 偏心率
        kpl(3) = incDeg*deg         ! 轨道倾角
        kpl(4) = raanDeg(i)*deg     ! 升交点赤经
        kpl(5) = 0d0                ! 近地点幅角
        kpl(6) = maDeg(i)*deg       ! 平近点角
        call kpl2cts(kpl, cts)      ! 归一化开普勒根数 -> 归一化位置速度

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
        sat%opMod%key_moonGravity = 1  ! 月球三体引力
        sat%opMod%key_sunGravity  = 1  ! 太阳三体引力
        sat%opMod%key_srp         = 1  ! 太阳光压
        sat%opMod%SRP_Model       = 2  !   ECOM 模型（9 参数框架，见 ecomParams）
        sat%opMod%ECOM            = ecomParams

        ! 卫星钟：观测用解析模型直接求值，calcMethod 保持缺省的公式法（1）——
        ! 任意时刻（含场景时间段外的查询，如星间观测历元 t + clk_R）都能求值。
        call sat%clock%init(t0, bias0 = clkScale*clkBias0(i), &
            &   drift = clkScale*clkDrift(i), &
            &   drift2 = clkScale*clkDrift2(i), name = trim(satName(i)))

        ! 另按场景历元生成一份钟差表：仿真本身用不到它（公式法不查表），
        ! 但 saveData 靠它把钟差序列写进 clock_*.txt。表的生成不改变
        ! calcMethod，也不会影响观测里钟差的计算。noise = 0，表值即模型值。
        call sat%clock%generateTable(epochs)

        call ss%addSatellite(sat)   ! addSatellite 按值复制，含钟模型
    enddo

    !====================================================================
    ! 3. 地面站（站钟为零，作时间基准）
    !====================================================================
    call sta%init(number = 1, name = staName, &
        &   longitude = staLon, latitude = staLat, altitude = staAlt)
    call sta%clock%init(t0)      ! 显式零钟：bias0 = drift = drift2 = 0
    call ss%addStation(sta)

    !====================================================================
    ! 4. 观测链路：两两建立 access（3 条星间 + 3 条星地，全部测距）
    !====================================================================
    nAccess = 0
    ! 星间：MEO1-MEO2, MEO1-MEO3, MEO2-MEO3
    do i = 1, nSat - 1
        do j = i + 1, nSat
            nAccess = nAccess + 1
            call addRangingAccess(ss, nAccess, trim(satName(i))//"-"//trim(satName(j)), &
                &   i, j, .true.)
        enddo
    enddo
    ! 星地：每颗星到北京站
    do i = 1, nSat
        nAccess = nAccess + 1
        call addRangingAccess(ss, nAccess, trim(satName(i))//"-"//trim(staName), &
            &   i, 1, .false.)
    enddo

    !====================================================================
    ! 5. 生成星历（真值）
    !    注意：必须在 ss%observation() 之前生成——观测要用 getPV 插值星历
    !====================================================================
    do i = 1, nSat
        write(*, "(A)") "Propagating "//trim(satName(i))//" ..."
        call ss%sat(i)%generateTrajectory(epochs)
        write(*, "(A,I0,A)") "  ephemeris points: ", ss%sat(i)%nEph, &
            &   " ("//trim(ss%ts%t0%toDate())//" ~ "//trim(ss%ts%tEnd%toDate())//" UTC)"
    enddo
    call ss%station(1)%generateTrajectory(ss%ts)

    !====================================================================
    ! 6. 生成观测数据（含收发钟差，见程序头说明）
    !====================================================================
    write(*, "(A)") "Generating ranging observations ..."
    call ss%observation
    call ss%reportAccess(.true.)

    !====================================================================
    ! 7. 输出
    !    注意：savetrajectory 默认会释放内存中的星历，故放在观测之后
    !====================================================================
    ! 7.1 卫星星历（真值）
    do i = 1, nSat
        write(fname, "(A,A,A)") trim(outDir)//"eph_", trim(satName(i)), ".txt"
        call ss%sat(i)%savetrajectory(STATE_TYPE_CARTESIAN, trim(fname))
        write(*, "(A)") "  saved: "//trim(fname)
    enddo

    ! 7.2 观测数据（每条第 k 条链路一个文件，文件名即链路名）
    do k = 1, ss%nAccess
        write(fname, "(A,A,A)") trim(outDir)//"obs_", trim(ss%access(k)%name), ".txt"
        call ss%access(k)%saveData(ss%ts, trim(fname))
        write(*, "(A,I0,A)") "  saved: "//trim(fname)//"   (epochs: ", &
            &   ss%access(k)%nData, ")"
    enddo

    ! 7.3 卫星钟（解析模型；公式法无状态表，文件即模型描述）
    do i = 1, nSat
        write(fname, "(A,A,A)") trim(outDir)//"clock_", trim(satName(i)), ".txt"
        call ss%sat(i)%clock%saveData(trim(fname), name = trim(satName(i)))
        write(*, "(A)") "  saved: "//trim(fname)
    enddo

    !====================================================================
    ! 8. 清理
    !====================================================================
    deallocate(epochs)
    call ss%freeSS

    write(*, "(A)") ""
    write(*, "(A)") "navSat simulation complete. Output: "//trim(outDir)

contains

!==============================================================================
! 建立一条测距 access 并加入场景。
! access 对象在过程内声明，每次调用都从类型默认值开始，避免上一条链路的
! 约束残留（initAccess 不重置约束字段）。
!   isSatSat = .true.  星间链路：地球遮挡生效；仰角/地影约束无意义
!   isSatSat = .false. 星地链路：仰角门限 minEle 生效
!==============================================================================
subroutine addRangingAccess(ss, num, aname, fromNum, toNum, isSatSat)
    type(type_satelliteScenario), intent(inout) :: ss
    integer, intent(in) :: num, fromNum, toNum
    character(len=*), intent(in) :: aname
    logical, intent(in) :: isSatSat
    type(type_access) :: ac

    if (isSatSat) then
        call ac%init(number = num, name = aname, &
            &   fromType = OBJECT_TYPE_SATELLITE, fromNum = fromNum, &
            &   toType = OBJECT_TYPE_SATELLITE, toNum = toNum, &
            &   obsType = OBS_TYPE_RANGING)
        ac%EarthBlockage = .true.    ! 视线被地球遮挡的历元不建链
        ac%MoonBlockage  = .false.
    else
        call ac%init(number = num, name = aname, &
            &   fromType = OBJECT_TYPE_SATELLITE, fromNum = fromNum, &
            &   toType = OBJECT_TYPE_STATION, toNum = toNum, &
            &   obsType = OBS_TYPE_RANGING, minEle = minEle)
        ! 测距是无线电链路，卫星进出地影不影响观测，故不用 shadow 约束
        ! （shadow 面向光学观测，类型缺省值为 .true.）
        ac%shadow = .false.
    endif

    call ss%addAccess(ac)
end subroutine addRangingAccess

end program navSat_simulation
