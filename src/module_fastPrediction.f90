!================================================================================
! MODULE: module_fastPrediction
!
! 说明：
!   快速轨道预报 —— 供"要算很多次、只要求近似"的场合，典型用途是数值状态转移
!   矩阵(numericalSTM_ECOM)里的有限差分：那里每算一次 STM 要跑十几次整弧传播，
!   用正规预报(NumOP)时单次迭代以小时计。
!
!   相比正规预报的简化(每一条都是为了摆脱文件 I/O):
!     1. 地球非球形引力只取 J2 与 C22/S22 两项，系数直接硬编码(取自
!        data/AtlasData/EGM2008.txt 的 2 阶归一化系数)，不读引力场文件；
!     2. EOP 参数取 0：不用极移与 DUT1，地球转角用解析 GMST(IAU1982)；
!     3. 太阳、月亮位置用解析平均轨道根数求，不读 JPL 行星历表；
!     4. 定步长积分，步长 FAST_STEP_SEC = 120 s，最后一步取剩余时长。
!
!   保留的力：中心引力 + J2 + C22/S22 + 日月三体 + 太阳光压(cannonball/ECOM，
!   与 module_perturbationSRP 同一套公式与系数约定)。力的开关由 type_fastOption
!   给出，调用者自行组装(数值 STM 里是按卫星的 opMod 组装)。
!
!   单位与坐标系：全程归一化(地球半径 / 时间单位)，J2000 惯性系，与 Atlas 其余
!   部分一致 —— 注意模块里所有力都返回归一化加速度(1 = acc_unit = 9.798 m/s^2)，
!   不像 module_perturbationSRP 那样返回 SI。时刻用 UTC 的 MJD，内部按
!   UTC + 32.184 s + ΔLS 换算到 TT 供日月解析式使用(ΔLS 取 module_time 内置表)。
!
!   积分器：经典四阶 Runge-Kutta，定步长、无误差控制。对 STM 而言，有限差分里
!   的两个传播用同一积分器，系统误差在差分中大部分相消，2 min 步长足够；绝对
!   轨迹精度则明显低于 NumOP(见速度/精度对比)。
!
! 用法：
!   type(type_fastOption) :: opt
!   opt%useJ2 = .true.;  opt%useJ22 = .true.
!   opt%useMoon = .true.; opt%useSun = .true.
!   opt%useSRP = .true.; opt%SRP_Model = 2; opt%ECOM = ecom
!   call fastPrediction(x0, tMJD_utc, dtSec, opt, x1)   ! x0/x1: 归一化位置速度
!
! 依赖：module_const, module_time, module_elements,
!       module_planets(质量比与天体编号), module_perturbationSRP(rho_s/地影)
!
! 作者：Shawn Hu
! 日期：2026-10-10
!================================================================================

module module_fastPrediction

    use module_const
    use module_time
    use module_elements
    use module_planets, only : moonEarthMassRatio, sunEarthMassRatio, nEarth, nMoon, nSun
    use module_perturbationSRP, only : rho_s, in_shadow

    implicit none
    private

    public :: type_fastOption, fastPrediction, fastAcceleration, FAST_STEP_SEC

    ! 定步长 [s]
    double precision, parameter :: FAST_STEP_SEC = 120d0

    ! EGM2008 的 2 阶归一化系数(与 data/AtlasData/EGM2008.txt 逐位一致)
    double precision, parameter :: C20_bar = -0.484165143790815d-03
    double precision, parameter :: C22_bar =  0.243938357328313d-05
    double precision, parameter :: S22_bar = -0.140027370385934d-05

    ! 由归一化系数换算到封闭式公式里用的形式：引力位里 P̄20 = sqrt(5)·P20、
    ! P̄22 = (sqrt(15)/6)·P22，故 J2 = -sqrt(5)·C̄20，
    ! C22(非归一化) = sqrt(15)/6·C̄22，S22 同理。
    double precision, parameter :: J2   = -dsqrt(5d0)*C20_bar
    double precision, parameter :: C22u = dsqrt(15d0)/6d0*C22_bar
    double precision, parameter :: S22u = dsqrt(15d0)/6d0*S22_bar

    ! amrSRP [m^2/kg] -> 归一化乘子(与 module_perturbationSRP 的 amrSRPu 同)
    double precision, parameter :: amrUnit = mass_unit/(length_unit*1d3)**2

    !---------------------------------------------------------------------------
    ! TYPE: type_fastOption
    ! 快速预报的力模型开关与参数(调用者按需组装)
    !---------------------------------------------------------------------------
    type :: type_fastOption
        logical :: useJ2   = .true.   ! 带谐 J2 项
        logical :: useJ22  = .true.   ! 扇谐 C22/S22 项
        logical :: useMoon = .false.  ! 月球三体
        logical :: useSun  = .false.  ! 太阳三体
        logical :: useSRP  = .false.  ! 太阳光压
        integer :: SRP_Model = 1      ! 1 - cannonball; 2 - ECOM(9 参数)
        double precision :: reflectivity = 0d0
        double precision :: amrSRP = 0.02d0
        double precision :: ECOM(9) = 0d0
    end type type_fastOption

contains

    !---------------------------------------------------------------------------
    ! SUBROUTINE: fastPrediction
    ! 从 x0 预报 dtSec 秒后的状态，定步长 RK4。
    !
    ! 输入：
    !   x0    - 初始状态(6)，归一化位置速度，J2000
    !   t0MJD - 初始时刻(UTC, MJD)
    !   dtSec - 预报时长 [s]，可为负(负号表示往回预报)；分步按 FAST_STEP_SEC
    !   opt   - 力模型选项
    ! 输出：
    !   x1    - 预报得到的状态(6)
    !---------------------------------------------------------------------------
    subroutine fastPrediction(x0, t0MJD, dtSec, opt, x1)
    implicit none
    double precision, intent(in) :: x0(6), t0MJD, dtSec
    type(type_fastOption), intent(in) :: opt
    double precision, intent(out) :: x1(6)
    integer :: nStep, k
    double precision :: hDay, hStepDay, hTU, tDay, y(6), k1(6), k2(6), k3(6), k4(6)

    x1 = x0
    if (dtSec == 0d0) return

    ! 步长用"天"推进时刻(力模型里的日月位置、GMST 都要绝对时刻),
    ! 但状态量的推进必须用归一化时间单位(TU) —— 导数 dx/dt 里 dr/dt = v 是
    ! ER/TU, 两者相差 86400/TIME_UNIT ~ 107 倍, 不能混用
    hStepDay = sign(FAST_STEP_SEC, dtSec)/86400d0
    nStep = int(abs(dtSec)/FAST_STEP_SEC)
    tDay = t0MJD
    y = x0

    do k = 1, nStep + 1
        hDay = hStepDay
        ! 最后一步取剩余时长(不足一个整步)
        if (k == nStep + 1) then
            hDay = dtSec/86400d0 - dble(k-1)*hStepDay
            if (hDay == 0d0) exit
        endif
        hTU = hDay*86400d0/time_unit

        call fastDeriv(tDay, y, opt, k1)
        call fastDeriv(tDay + 0.5d0*hDay, y + 0.5d0*hTU*k1, opt, k2)
        call fastDeriv(tDay + 0.5d0*hDay, y + 0.5d0*hTU*k2, opt, k3)
        call fastDeriv(tDay + hDay,       y + hTU*k3,        opt, k4)
        y = y + (hTU/6d0)*(k1 + 2d0*k2 + 2d0*k3 + k4)
        tDay = tDay + hDay
    enddo

    x1 = y
    end subroutine fastPrediction

    !---------------------------------------------------------------------------
    ! SUBROUTINE: fastDeriv
    ! 运动方程右端：dr/dt = v, dv/dt = a(t, r, v)
    !---------------------------------------------------------------------------
    subroutine fastDeriv(tMJD, x, opt, dx)
    implicit none
    double precision, intent(in) :: tMJD, x(6)
    type(type_fastOption), intent(in) :: opt
    double precision, intent(out) :: dx(6)

    dx(1:3) = x(4:6)
    call fastAcceleration(tMJD, x, opt, dx(4:6))
    end subroutine fastDeriv

    !---------------------------------------------------------------------------
    ! SUBROUTINE: fastAcceleration
    ! 归一化加速度 = 中心引力 + J2 + C22/S22 + 日月三体 + 太阳光压
    !
    !   中心引力在 J2000 惯性系里算；J2 与 C22/S22 在 ECEF 里算(ECI->ECF 只用
    !   解析 GMST，即 EOP 取 0)，算完再转回惯性系。速度为光压的纬度幅角所需，
    !   与 module_perturbationSRP 一样从位置速度反算(加速度对位置的依赖不受影响)。
    !
    ! 输入：tMJD - 时刻(UTC, MJD)；x - 归一化状态(6, J2000)；opt - 力模型
    ! 输出：acc  - 归一化加速度(ER/TU^2, J2000)
    !---------------------------------------------------------------------------
    subroutine fastAcceleration(tMJD, x, opt, acc)
    implicit none
    double precision, intent(in) :: tMJD, x(6)
    type(type_fastOption), intent(in) :: opt
    double precision, intent(out) :: acc(3)
    double precision :: r(3), rn, rn2, rn3, rn5, rn7
    double precision :: theta, ct, st, xe, ye, ze, re2, ae(3), f
    double precision :: rSun(6), rMoon(6), d(3), dn, rho1, kappa, amrSRPu, ksrp, du(3)
    double precision :: spv(6), uLat
    logical :: flag_shadow

    r = x(1:3)
    rn2 = dot_product(r, r)
    rn  = dsqrt(rn2)
    rn3 = rn2*rn

    ! --- 中心引力 ---
    acc = -r/rn3

    ! --- J2 与 C22/S22 ---
    if (opt%useJ2 .or. opt%useJ22) then
        ! ECI -> ECEF: 只用解析 GMST(EOP 取 0)
        theta = fastGMST(tMJD)
        ct = dcos(theta); st = dsin(theta)
        xe =  ct*r(1) + st*r(2)
        ye = -st*r(1) + ct*r(2)
        ze =  r(3)
        re2 = xe*xe + ye*ye + ze*ze

        ! 2 阶位的标量部分 f = C20(3z^2-r^2)/2 + 3 C22(x^2-y^2) + 6 S22 x y
        ! (C20 = -J2; 用非归一化系数; 参考半径 = length_unit = 1 ER)
        f = 0d0
        if (opt%useJ2)  f = f - J2*(3d0*ze*ze - re2)/2d0
        if (opt%useJ22) f = f + 3d0*C22u*(xe*xe - ye*ye) + 6d0*S22u*xe*ye

        ! df/dr
        ae = 0d0
        if (opt%useJ2) then
            ae(1) = J2*xe
            ae(2) = J2*ye
            ae(3) = -2d0*J2*ze
        endif
        if (opt%useJ22) then
            ae(1) = ae(1) + 6d0*C22u*xe + 6d0*S22u*ye
            ae(2) = ae(2) - 6d0*C22u*ye + 6d0*S22u*xe
        endif

        ! a_ecef = mu Re^2 [ df/dr /r^5 - 5 r f /r^7 ]
        rn5 = rn3*rn2
        rn7 = rn5*rn2
        ae = ae/rn5 - 5d0*(/xe, ye, ze/)*f/rn7

        ! ECEF -> ECI
        acc(1) = acc(1) + ct*ae(1) - st*ae(2)
        acc(2) = acc(2) + st*ae(1) + ct*ae(2)
        acc(3) = acc(3) + ae(3)
    endif

    ! --- 日月三体(解析位置, 归一化 ER; mu_earth = 1) ---
    if (opt%useSun .or. opt%useSRP) call fastSunPosition(tMJD, rSun)
    if (opt%useMoon) then
        call fastMoonPosition(tMJD, rMoon)
        d = rMoon(1:3) - r
        dn = dsqrt(dot_product(d, d))
        rho1 = dsqrt(dot_product(rMoon(1:3), rMoon(1:3)))
        acc = acc + moonEarthMassRatio*(d/dn**3 - rMoon(1:3)/rho1**3)
    endif
    if (opt%useSun) then
        d = rSun(1:3) - r
        dn = dsqrt(dot_product(d, d))
        rho1 = dsqrt(dot_product(rSun(1:3), rSun(1:3)))
        acc = acc + sunEarthMassRatio*(d/dn**3 - rSun(1:3)/rho1**3)
    endif

    ! --- 太阳光压 ---
    if (opt%useSRP) then
        du = r - rSun(1:3)
        dn = dsqrt(dot_product(du, du))
        du = du/dn
        call in_shadow(r, rSun(1:3), flag_shadow)
        if (.not. flag_shadow) then
            if (opt%SRP_Model == 1) then
                ! cannonball(与 module_perturbationSRP 同一式; rho_s 为归一化加速度)
                kappa = 1d0 + opt%reflectivity
                amrSRPu = opt%amrSRP*amrUnit
                ksrp = kappa*amrSRPu*rho_s*(1.49591492d8/6378.1363d0/dn)**2
                acc = acc + ksrp*du
            elseif (opt%SRP_Model == 2) then
                ! ECOM(9 参数, 顺序 D0 Y0 B0 Bc Bs Dc Ds Yc Ys)
                call latitudeArgument(x, uLat)
                spv = 0d0
                spv(1:3) = du
                call ecomAcceleration(spv, rSun(1:3), x(1:3), uLat, opt%ECOM, acc)
            endif
        endif
    endif
    end subroutine fastAcceleration

    !---------------------------------------------------------------------------
    ! SUBROUTINE: ecomAcceleration
    ! ECOM 光压(SRP_Model = 2)，参数顺序与 module_perturbationSRP 一致：
    !   1 D0  2 Y0  3 B0  4 Bc  5 Bs  6 Dc  7 Ds  8 Yc  9 Ys
    !   F = (D0 + Dc*cos u + Ds*sin u)*e_D + (Y0 + Yc*cos u + Ys*sin u)*e_Y
    !     + (B0 + Bc*cos u + Bs*sin u)*e_B
    ! 输入：eD 由太阳指向卫星的单位向量；rSun/r 归一化位置；u 纬度幅角
    ! 输出：acc 累加(归一化加速度；ECOM 参数本身就是归一化加速度单位)
    !---------------------------------------------------------------------------
    subroutine ecomAcceleration(eD, rSun, r, u, ECOM, acc)
    implicit none
    double precision, intent(in) :: eD(3), rSun(3), r(3), u, ECOM(9)
    double precision, intent(inout) :: acc(3)
    double precision :: eY(3), eB(3)

    ! e_Y = r_sun x r_sat / |...|; e_B = e_D x e_Y (与 module_perturbationSRP 同)
    eY(1) = rSun(2)*r(3) - rSun(3)*r(2)
    eY(2) = rSun(3)*r(1) - rSun(1)*r(3)
    eY(3) = rSun(1)*r(2) - rSun(2)*r(1)
    eY = eY/dsqrt(dot_product(eY, eY))
    eB(1) = eD(2)*eY(3) - eD(3)*eY(2)
    eB(2) = eD(3)*eY(1) - eD(1)*eY(3)
    eB(3) = eD(1)*eY(2) - eD(2)*eY(1)

    acc = acc + (ECOM(1) + ECOM(6)*dcos(u) + ECOM(7)*dsin(u))*eD &
        &     + (ECOM(2) + ECOM(8)*dcos(u) + ECOM(9)*dsin(u))*eY &
        &     + (ECOM(3) + ECOM(4)*dcos(u) + ECOM(5)*dsin(u))*eB
    end subroutine ecomAcceleration

    !---------------------------------------------------------------------------
    ! SUBROUTINE: latitudeArgument
    ! 纬度幅角 u = omega + f：由轨道面法向 h = r x v 定升交点方向 n = z x h，
    ! u 即位置在轨道面内相对 n 的夹角(朝向运动方向为正)。
    !---------------------------------------------------------------------------
    subroutine latitudeArgument(x, u)
    implicit none
    double precision, intent(in) :: x(6)
    double precision, intent(out) :: u
    double precision :: h(3), n(3), nh(3), hh(3), hn, un

    h(1) = x(2)*x(6) - x(3)*x(5)
    h(2) = x(3)*x(4) - x(1)*x(6)
    h(3) = x(1)*x(5) - x(2)*x(4)
    hn = dsqrt(dot_product(h, h))
    if (hn < 1d-14) then
        u = 0d0
        return
    endif
    n = (/-h(2), h(1), 0d0/)
    un = dsqrt(dot_product(n, n))
    if (un < 1d-14) then
        nh = (/1d0, 0d0, 0d0/)      ! 近赤道轨道: 以 x 轴作参考
    else
        nh = n/un
    endif
    hh = h/hn
    ! m_hat = h_hat x n_hat (轨道面内超前 n_hat 90 度)
    u = atan2( x(1)*(hh(2)*nh(3) - hh(3)*nh(2)) &
        &    + x(2)*(hh(3)*nh(1) - hh(1)*nh(3)) &
        &    + x(3)*(hh(1)*nh(2) - hh(2)*nh(1)), &
        &      dot_product(x(1:3), nh) )
    end subroutine latitudeArgument

    !---------------------------------------------------------------------------
    ! FUNCTION: fastGMST
    ! 格林尼治平恒星时(IAU1982)，输入 UTC 的 MJD，返回弧度
    !---------------------------------------------------------------------------
    function fastGMST(tMJD) result(theta)
    implicit none
    double precision, intent(in) :: tMJD
    double precision :: theta, d

    d = (tMJD + 2400000.5d0) - 2451545d0
    theta = modulo(280.46061837d0 + 360.98564736629d0*d, 360d0)*deg2rad
    end function fastGMST

    !---------------------------------------------------------------------------
    ! SUBROUTINE: fastSunPosition
    ! 解析太阳位置(地心, J2000, 归一化 ER)。根数与 module_planets 的
    ! mean_solar_position 相同，但不直接调用它 —— 那个子程序每调用一次都打印
    ! 一行警告，放在快速预报的热路径上会淹没输出。
    !---------------------------------------------------------------------------
    subroutine fastSunPosition(tMJD, PV)
    implicit none
    double precision, intent(in) :: tMJD
    double precision, intent(out) :: PV(6)
    double precision :: d, T, elements(6)

    d = TT_days(tMJD) - 51544.50074287037d0
    T = d/36525d0
    elements(1) = 1.00000102d0*AU_TO_EARTH_RADIUS
    elements(2) = 0.016709d0
    elements(3) = 23.4393d0*deg2rad
    elements(4) = 0d0
    elements(5) = 282.9373d0*deg2rad + 0.32d0*deg2rad*T
    elements(6) = 357.5291d0*deg2rad + 0.9856d0*deg2rad*d
    call kpl2cts(elements, PV)
    end subroutine fastSunPosition

    !---------------------------------------------------------------------------
    ! SUBROUTINE: fastMoonPosition
    ! 解析月亮位置(地心, J2000, 归一化 ER)，根数与 module_planets 的
    ! mean_lunar_position 相同
    !---------------------------------------------------------------------------
    subroutine fastMoonPosition(tMJD, PV)
    implicit none
    double precision, intent(in) :: tMJD
    double precision, intent(out) :: PV(6)
    double precision :: d, T, elements(6)

    d = TT_days(tMJD) - 51544.50074287037d0
    T = d/36525d0
    elements(1) = 384747.981d0/length_unit
    elements(2) = 0.054880d0
    elements(3) = 5.1298d0*deg2rad
    elements(4) = 125.0446d0*deg2rad - 1934.14d0*deg2rad*T
    elements(5) = 318.3087d0*deg2rad + 6003.15d0*deg2rad*T
    elements(6) = 134.9634d0*deg2rad + 13.0650d0*deg2rad*d
    call kpl2cts(elements, PV)
    end subroutine fastMoonPosition

    !---------------------------------------------------------------------------
    ! FUNCTION: TT_days
    ! UTC 的 MJD -> TT 的 MJD（用 module_time 内置的跳秒表，不读文件）
    !---------------------------------------------------------------------------
    function TT_days(tMJD_utc) result(tMJD_tt)
    implicit none
    double precision, intent(in) :: tMJD_utc
    double precision :: tMJD_tt
    double precision :: day, sec, dayTT, secTT

    day = dble(floor(tMJD_utc))
    sec = (tMJD_utc - day)*86400d0
    call UTC2TT_daySec(day, sec, dayTT, secTT)
    tMJD_tt = dayTT + secTT/86400d0
    end function TT_days

end module module_fastPrediction
