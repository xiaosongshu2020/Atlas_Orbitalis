!================================================================================
! 测试：module_fastPrediction（快速轨道预报）
!
! 与正规预报 NumOP 对比，共四组检查：
!   [1] 引力场：快速模块的 J2+C22/S22 与 GravityField(degree=2,order=2) 对比
!       —— 校验硬编码系数与公式符号（用模块内部同一个 GMST 把 ECEF 转到 ECI）
!   [2] 速度：1 天传播耗时，NumOP vs 快速预报（MEO / LEO）
!   [3] 精度：1 天后位置速度差（MEO / LEO），力模型两边取一致
!   [4] STM  ：有限差分 STM 两法对比（数值 STM 是快速预报的实际用途）
!
! 结论看打印：速度应快 1~2 个量级，绝对精度是"米级以下"，而 STM 的相对差异
! 应远小于绝对轨迹差异（差分里系统误差相消）。
!================================================================================

program test_fastPrediction

    use module_Atlas
    use module_const
    use module_time
    use module_elements
    use module_Gravity
    use module_fastPrediction

    implicit none

    type(type_satellite) :: sat
    type(type_t) :: t0, ep
    type(type_fastOption) :: opt
    type(GravityField) :: gf2
    double precision :: kpl(6), x0(6), xf(6), xr(6)
    double precision :: ct, st, theta, rg(3), ag(3), af(3), ae(3), aecf(3), d
    real(8) :: c0, c1, tNumOP, tFast, tStmF, tStmN, tStmE
    double precision :: StmF(6, 11), StmN(6, 11), x1n(6), S15(15, 15)
    real(8), parameter :: ecom(9) = (/ 1.4d-8, 3.0d-9, -1.5d-9, 0d0, 0d0, &
        &                             2.0d-9, -2.0d-9, 0d0, 0d0 /)
    integer :: i, iear

    call t0%init(2025d0, 1d0, 1d0, 0d0, 0d0, 0d0)
    ep = timeAdd(t0, 1d0)

    opt%useJ2   = .true.
    opt%useJ22  = .true.
    opt%useMoon = .true.
    opt%useSun  = .true.
    opt%useSRP  = .true.
    opt%SRP_Model = 2
    opt%ECOM = ecom

    print *, ""
    print *, "================ [1] gravity field: J2+C22/S22 vs GravityField(2x2) ================"
    call gf2%init("EarthGravity", "Earth", 2, 2)
    gf2%GM = EARTH_GM
    gf2%referenceRadius = LENGTH_UNIT*1000d0
    call gf2%read_gravity_coefficients_4col("data/AtlasData/EGM2008.txt")

    ! 取 MEO 初始位置(ECI, 归一化)，按模块内部同一个 GMST 转到 ECEF 做比较
    kpl = (/ 27906.1d0*u_km, 0.001d0, 55d0*deg, 0d0, 0d0, 0d0 /)
    call kpl2cts(kpl, x0)
    ! 两边必须用同一个地球转角：模块内部用解析 GMST，故这里复刻同一公式
    theta = fastGMST_test(t0%val)
    ct = dcos(theta); st = dsin(theta)
    rg(1) =  ct*x0(1) + st*x0(2)          ! ECI -> ECEF
    rg(2) = -st*x0(1) + ct*x0(2)
    rg(3) =  x0(3)
    ! GravityField 要米、返回 m/s^2
    call gf2%acceleration(rg*LENGTH_UNIT*1000d0, 2, 2, ag)
    ! 快速模块：只留引力(关掉日月与光压)，取总加速度再减中心项
    opt%useMoon = .false.; opt%useSun = .false.; opt%useSRP = .false.
    call fastAcceleration(t0%val, x0, opt, af)          ! 归一化, ECI
    ae = af + x0(1:3)/norm2(x0(1:3))**3                  ! 去掉中心项
    aecf(1) = ct*ae(1) + st*ae(2)                        ! ECI -> ECEF, 同一方向
    aecf(2) = -st*ae(1) + ct*ae(2)
    aecf(3) = ae(3)
    d = maxval(abs(aecf - ag/acc_unit))                  ! m/s^2 -> 归一化 后作差
    write(*, "(A,ES12.3,A)") " max|fast - GravityField(2x2)| = ", d, " (normalized accel)"
    print *, " (expected ~1e-10: the only difference is the unmodelled C21/S21)"

    opt%useMoon = .true.; opt%useSun = .true.; opt%useSRP = .true.

    print *, ""
    print *, "================ [2][3] speed and accuracy: NumOP vs fast prediction ================"
    do iear = 1, 2
        if (iear == 1) then
            kpl = (/ 27906.1d0*u_km, 0.001d0, 55d0*deg, 0d0, 0d0, 0d0 /)
            write(*, "(A)") " --- MEO (a = 27906 km, i = 55 deg) ---"
        else
            kpl = (/ (6378.1363d0+600d0)*u_km, 0.001d0, 60d0*deg, 0d0, 0d0, 0d0 /)
            write(*, "(A)") " --- LEO (a = 6978 km, i = 60 deg) ---"
        endif
        call kpl2cts(kpl, x0)

        ! (a) 正规预报: 与快速预报取同一力模型(degree=2/order=2 + 日月 + ECOM 光压)
        call sat%init(t0, x0, eleType=STATE_TYPE_CARTESIAN, propagator='numerical', &
            &   degree=2, order=2)
        sat%opMod%usePV=1; sat%opMod%key_moonGravity=1; sat%opMod%key_sunGravity=1
        sat%opMod%key_srp=1; sat%opMod%SRP_Model=2; sat%opMod%ECOM=ecom
        call sat%propagate(ep)                   ! 预热
        call cpu_time(c0)
        do i = 1, 20
            sat%states%t = t0            ! 每轮都要把初值整组复位
            sat%states%elements = x0     ! (只复位 t 会让 states 停在 ep 上)
            call sat%propagate(ep)
        enddo
        call cpu_time(c1)
        tNumOP = (c1-c0)/20d0*1000d0
        x1n = sat%states%elements

        ! (b) 快速预报
        call fastPrediction(x0, t0%val, 86400d0, opt, xf)     ! 预热
        call cpu_time(c0)
        do i = 1, 200
            call fastPrediction(x0, t0%val, 86400d0, opt, xf)
        enddo
        call cpu_time(c1)
        tFast = (c1-c0)/200d0*1000d0

        write(*, "(A,F10.3,A,F10.3,A,F8.1,A)") " 1-day propagation: NumOP ", tNumOP, &
            &   " ms   fast ", tFast, " ms   speedup ", tNumOP/max(tFast,1d-6), " x"
        write(*, "(A,ES12.4,A)") " position difference after 1 day |dr| = ", &
            &   norm2(xf(1:3)-x1n(1:3))*LENGTH_UNIT*1000d0, " m"
        write(*, "(A,ES12.4,A)") " velocity difference after 1 day |dv| = ", &
            &   norm2(xf(4:6)-x1n(4:6))*LENGTH_UNIT*1000d0/TIME_UNIT, " m/s"

        ! (c) STM: 有限差分(步长同 numericalSTM_ECOM), 快速预报 vs 正规预报
        call cpu_time(c0)
        call FDstm_fast(x0, t0%val, opt, StmF)
        call cpu_time(c1)
        tStmF = (c1-c0)*1000d0
        call cpu_time(c0)
        call FDstm_numop(x0, t0, ep, StmN)
        call cpu_time(c1)
        tStmN = (c1-c0)*1000d0
        call cpu_time(c0)
        call numericalSTM_ECOM(sat, 86400d0/time_unit, S15)   ! 真正接线的那个
        call cpu_time(c1)
        tStmE = (c1-c0)*1000d0
        write(*, "(A,F9.2,A,F9.2,A,F8.1,A)") " one STM call: fast ", tStmF, &
            &   " ms   regular ", tStmN, " ms   speedup ", tStmN/max(tStmF,1d-6), " x"
        write(*, "(A,F9.2,A)") " (the two above use the same FD scheme, 11 columns); "// &
            &   "ODor%numericalSTM_ECOM (15 columns, wired in) = ", tStmE, " ms"
        write(*, "(A,ES12.3)") " FD STM relative difference (vs NumOP) = ", &
            &   maxval(abs(StmF-StmN))/maxval(abs(StmN))
    enddo

    print *, ""
    print *, "test_fastPrediction: done"

contains

    ! 测试里复刻模块内部的 GMST(IAU1982)，保证两边用同一个地球转角
    function fastGMST_test(tMJD) result(th)
        implicit none
        double precision, intent(in) :: tMJD
        double precision :: th, dd
        dd = (tMJD + 2400000.5d0) - 2451545d0
        th = modulo(280.46061837d0 + 360.98564736629d0*dd, 360d0)*deg2rad
    end function fastGMST_test

    ! 快速预报版本的有限差分 STM(与 numericalSTM_ECOM 同一套步长)
    subroutine FDstm_fast(x0, tMJD, op, S)
        implicit none
        double precision, intent(in) :: x0(6), tMJD
        type(type_fastOption), intent(in) :: op
        double precision, intent(out) :: S(6, 11)
        type(type_fastOption) :: opp
        double precision :: x1(6), xpert(6), xr(6), h, dtSec
        integer :: k

        dtSec = 86400d0
        call fastPrediction(x0, tMJD, dtSec, op, x1)
        h = 1d-4
        do k = 1, 6
            xpert = x0; xpert(k) = xpert(k) + h
            call fastPrediction(xpert, tMJD, dtSec, op, xr)
            S(1:6, k) = (xr - x1)/h
        enddo
        h = 1d-6
        do k = 1, 5
            xpert = x0
            opp = op
            opp%ECOM(k) = opp%ECOM(k) + h
            call fastPrediction(xpert, tMJD, dtSec, opp, xr)
            S(1:6, 6+k) = (xr - x1)/h
        enddo
    end subroutine FDstm_fast

    ! 正规预报版本的有限差分 STM(用 sat%propagate)
    subroutine FDstm_numop(x0, t0, ep, S)
        implicit none
        double precision, intent(in) :: x0(6)
        type(type_t), intent(in) :: t0, ep
        double precision, intent(out) :: S(6, 11)
        type(type_satellite) :: satn, satp
        double precision :: x1(6), xr(6), h
        integer :: k

        call satn%init(t0, x0, eleType=STATE_TYPE_CARTESIAN, propagator='numerical', &
            &   degree=2, order=2)
        satn%opMod%usePV=1; satn%opMod%key_moonGravity=1; satn%opMod%key_sunGravity=1
        satn%opMod%key_srp=1; satn%opMod%SRP_Model=2; satn%opMod%ECOM=ecom
        call satn%propagate(ep)
        x1 = satn%states%elements

        h = 1d-4
        do k = 1, 6
            satp = satn
            satp%states%t = t0
            satp%states%elements = x0
            satp%states%elements(k) = satp%states%elements(k) + h
            call satp%propagate(ep)
            xr = satp%states%elements
            S(1:6, k) = (xr - x1)/h
        enddo
        h = 1d-6
        do k = 1, 5
            satp = satn
            satp%states%t = t0
            satp%states%elements = x0
            satp%opMod%ECOM(k) = satp%opMod%ECOM(k) + h
            call satp%propagate(ep)
            xr = satp%states%elements
            S(1:6, 6+k) = (xr - x1)/h
        enddo
    end subroutine FDstm_numop

end program test_fastPrediction
