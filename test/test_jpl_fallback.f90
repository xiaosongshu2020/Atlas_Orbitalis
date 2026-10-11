!================================================================================
! 测试：JPL 星历读取失败时 planetsPV 的解析回退（module_planets）
!
! 本程序不读任何数据文件，可以任意目录运行 —— 运行目录决定了星历能否找到，
! 从而切换两条路径：
!
!   [1] JPL 路径：从仓库根目录运行（能找到 data/AtlasData/JPLEPH_DE440），
!       pleph 状态应为 0，返回 JPL 结果。
!   [2] 回退路径：从没有 data/ 的目录运行（如 /tmp），星历打不开，应当看到
!       每体一行 "Warning @planetsPV"，仍返回解析结果且不崩溃。
!
! 每次运行都会检查（两条路径都查）：
!   - 速度模：月球 ≈ 1.02 km/s、太阳 ≈ 29.78 km/s。这一项专门用来抓单位错误
!     （AU / 地球半径 / km 混用，以及 mean_solar_position 的 mu=1 缺陷 ——
!     若不修正，太阳速度只有 0.05 km/s）；
!   - 告警只打一次：用 20 个互不相同的历元调用（绕开 1e-8 天的结果缓存），
!     "Warning @planetsPV" 的行数应为 2（日、月各一），而不是 40。
!
! 用法（从仓库根目录）：
!   ./bin/test_jpl_fallback.exe            ! JPL 路径 + 日月检查
!   cd /tmp && /d/home/projects/atlas/bin/test_jpl_fallback.exe    ! 回退路径
!   ./bin/test_jpl_fallback.exe planet     ! 额外请求木星（无解析模型）
!                                          ! 回退时会打印原因并终止
!================================================================================
program test_jpl_fallback
    use module_planets
    implicit none
    double precision :: TT, PV(6), vMoon, vSun
    integer :: i, iStat, nArg
    character(len=64) :: arg
    double precision, parameter :: testEpoch = 60700d0   ! MJD, TT (2025-01-24)

    print *, "=============================================================="
    print *, " JPL ephemeris fallback test (module_planets)"
    print *, "=============================================================="
    print *, " working directory decides whether the ephemeris is found:"
    print *, "   repository root -> JPL path,  other dir -> fallback path"
    print *, ""

    !>>> [1] 日月各一次：看状态码、位置、速度模
    TT = testEpoch

    call planetsPV(TT, nEarth, nMoon, PV)
    call pleph_status(iStat)
    vMoon = dsqrt(PV(4)**2 + PV(5)**2 + PV(6)**2)
    write (*, '(a,i3)')     " [1] pleph status       : ", iStat
    write (*, '(a,3f16.3)') "     Moon position  km  : ", PV(1:3)
    write (*, '(a,3f16.6)') "     Moon velocity km/s : ", PV(4:6)
    write (*, '(a,f12.6)')  "     |v| Moon      km/s : ", vMoon

    call planetsPV(TT, nEarth, nSun, PV)
    call pleph_status(iStat)
    vSun = dsqrt(PV(4)**2 + PV(5)**2 + PV(6)**2)
    write (*, '(a,i3)')     "     pleph status       : ", iStat
    write (*, '(a,3f16.3)') "     Sun  position  km  : ", PV(1:3)
    write (*, '(a,3f16.6)') "     Sun  velocity km/s : ", PV(4:6)
    write (*, '(a,f12.6)')  "     |v| Sun       km/s : ", vSun
    print *, ""

    !>>> 速度模判断（两条路径都必须通过）
    print *, " expect |v| Moon ~ 1.02 km/s, |v| Sun ~ 29.78 km/s"
    if (vMoon < 0.5d0 .or. vMoon > 2.0d0) then
        print *, " FAIL: |v| Moon is outside 0.5 - 2.0 km/s  (unit error?)"
    else
        print *, " OK  : |v| Moon has the right magnitude"
    endif
    if (vSun < 20d0 .or. vSun > 40d0) then
        print *, " FAIL: |v| Sun is outside 20 - 40 km/s  (unit error?)"
    else
        print *, " OK  : |v| Sun has the right magnitude"
    endif
    print *, ""

    !>>> [2] 告警只打一次：20 个不同历元（1e-6 天的间隔绕开结果缓存）
    print *, " [2] calling planetsPV 20 times at distinct epochs..."
    do i = 1, 20
        TT = testEpoch + i*1d-6
        call planetsPV(TT, nEarth, nMoon, PV)
        call planetsPV(TT, nEarth, nSun, PV)
    enddo
    print *, "     done.  On the fallback path the 'Warning @planetsPV'"
    print *, "     lines above should number exactly 2 (one per body), not 40."
    print *, ""

    !>>> [3] 超出星历时间范围：随仓库分发的 DE440 只覆盖 1949-12-14 ~ 2150-01-22
    !        （MJD 33264 ~ 106352），这里用 MJD 120000 触发状态码 2。
    !        注意这与 [2] 不同：文件是好的，只是这一次的历元取不到。
    TT = 120000d0
    print *, " [3] epoch outside the ephemeris span (MJD 120000):"
    call planetsPV(TT, nEarth, nMoon, PV)
    call pleph_status(iStat)
    write (*, '(a,i3)')     "     pleph status       : ", iStat
    write (*, '(a,3f16.3)') "     Moon position  km  : ", PV(1:3)
    print *, "     (status 2 = epoch out of range, not a broken file;"
    print *, "      the warning above appears only once per body per run)"
    print *, ""

    !>>> [4] 可选：请求没有解析模型的天体（回退时应终止并说明原因）
    nArg = command_argument_count()
    if (nArg >= 1) then
        call get_command_argument(1, arg)
        if (trim(arg) == 'planet') then
            print *, " [4] planetary request (no analytical model exists):"
            call planetsPV(testEpoch, nEarth, nJupiter, PV)
            write (*, '(a,3f16.3)') "     Jupiter position km: ", PV(1:3)
            print *, ""
        endif
    endif

    print *, "=== test finished without abort ==="

end program
