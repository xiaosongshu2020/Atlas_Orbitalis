program test_sp3
    !---------------------------------------------------------------------
    ! 测试程序：使用 type_sp3 对象读取 SP3c 格式精密星历文件
    ! 测试文件：data/input/isr10620.sp3
    !---------------------------------------------------------------------
    use module_sp3
    use ErrorHandler
    implicit none

    ! 参数
    character(len=*), parameter :: SP3_FILE = "data/input/isr10620.sp3"

    ! 变量
    type(type_sp3) :: sp3
    double precision :: pos(3)
    integer :: i, j
    integer :: nPrintEpoch, nPrintSat

    print *, "============================================"
    print *, "  SP3c File Reader Test (type_sp3 object)"
    print *, "============================================"
    print *, "File: ", trim(SP3_FILE)
    print *, ""

    ! ===== 测试 load 方法 =====
    print *, ">>> Testing load() method ..."
    call sp3%load(SP3_FILE)
    if (errCode > 0) then
        call PrintErrorInfo
        stop
    endif
    print *, "    Load successful!"
    print *, ""

    ! ===== 输出头信息 =====
    print *, "--- Header Information ---"
    print *, "  Version:              ", sp3%header%ver
    print *, "  Position/velocity:    ", sp3%header%pvType
    print *, "  File type:            ", sp3%header%fType
    print *, "  Time system:          ", sp3%header%tSys
    print *, "  Agency:               ", sp3%header%meta%Agency
    print *, "  Coordinate frame:     ", sp3%header%meta%Frame
    print *, "  Orbit type:           ", sp3%header%meta%oType
    print *, "  Data used:            ", sp3%header%meta%uData
    print *, ""

    print *, "  Start epoch (Gregorian):"
    print '("    ", I4.4, "-", I2.2, "-", I2.2, " ", I2.2, ":", I2.2, ":", F11.8)', &
        sp3%header%t0%cal%yr, sp3%header%t0%cal%mo, sp3%header%t0%cal%day, &
        sp3%header%t0%cal%hr, sp3%header%t0%cal%mi, sp3%header%t0%cal%sec
    print '("  Start epoch (MJD):       day=", I5, ", frac=", F16.13)', &
        sp3%header%t0%mjd%day, sp3%header%t0%mjd%frac
    print '("  Start epoch (GPS week):  week=", I4, ", sec=", F15.8)', &
        sp3%header%t0%week%week, sp3%header%t0%week%sec
    print *, ""

    print *, "  Number of epochs:     ", sp3%header%nEpoch
    print *, "  Number of satellites: ", sp3%header%nSat
    print '("  Time interval (dt):    ", F14.8, " seconds")', sp3%header%dt
    print '("  Position accuracy base:", F10.7, " km")', sp3%header%pvAccBase
    print '("  Clock accuracy base:   ", F12.9, " us")', sp3%header%clkAccBase
    print *, ""

    ! ===== 输出卫星列表（前 20 颗） =====
    print *, "--- Satellite List (first 20) ---"
    do i = 1, min(20, sp3%header%nSat)
        print '("  [", I3, "] ", A1, I2.2)', i, sp3%header%sat(i)%sys, sp3%header%sat(i)%SatNum
    enddo
    if (sp3%header%nSat > 20) then
        print '("  ... and ", I0, " more satellites")', sp3%header%nSat - 20
    endif
    print *, ""

    ! ===== 输出位置数据（前 3 个历元，前 5 颗卫星） =====
    nPrintEpoch = min(3, sp3%header%nEpoch)
    nPrintSat   = min(5, sp3%header%nSat)

    print *, "--- Position Data (first ", nPrintEpoch, " epochs, first ", nPrintSat, " satellites) ---"
    do i = 1, nPrintEpoch
        print '("  Epoch ", I0, ": ", I4.4, "-", I2.2, "-", I2.2, " ", I2.2, ":", I2.2, ":", F11.8)', &
            i, sp3%data(i,1)%t%yr, sp3%data(i,1)%t%mo, sp3%data(i,1)%t%day, &
            sp3%data(i,1)%t%hr, sp3%data(i,1)%t%mi, sp3%data(i,1)%t%sec
        do j = 1, nPrintSat
            print '("    ", A1, I2.2, ":  X=", F14.6, "  Y=", F14.6, "  Z=", F14.6, "  Clk=", F14.6)', &
                sp3%data(i,j)%sat%sys, sp3%data(i,j)%sat%SatNum, &
                sp3%data(i,j)%pos(1), sp3%data(i,j)%pos(2), sp3%data(i,j)%pos(3), &
                sp3%data(i,j)%clk
        enddo
        print *, ""
    enddo

    ! ===== 测试 getPosition 方法 =====
    print *, "--- Testing getPosition() method ---"

    ! 测试 1: 查询起始历元 G01 的位置
    print *, "  Test 1: Query G01 at start epoch (2026-05-10 00:00:00)"
    pos = sp3%getPosition(2026, 5, 10, 0, 0, 0.0d0, "G01")
    if (errCode == 0) then
        print '("    G01:  X=", F14.6, "  Y=", F14.6, "  Z=", F14.6, " km")', pos(1), pos(2), pos(3)
    else
        print *, "    ERROR: ", trim(errMsg)
        call PrintErrorInfo
    endif

    ! 测试 2: 查询第 2 个历元 G02 的位置（起始时间 + 900 秒）
    print *, "  Test 2: Query G02 at epoch 2 (2026-05-10 00:15:00)"
    pos = sp3%getPosition(2026, 5, 10, 0, 15, 0.0d0, "G02")
    if (errCode == 0) then
        print '("    G02:  X=", F14.6, "  Y=", F14.6, "  Z=", F14.6, " km")', pos(1), pos(2), pos(3)
    else
        print *, "    ERROR: ", trim(errMsg)
        call PrintErrorInfo
    endif

    ! 测试 3: 查询 GLONASS 卫星 R02
    print *, "  Test 3: Query R02 at start epoch"
    pos = sp3%getPosition(2026, 5, 10, 0, 0, 0.0d0, "R02")
    if (errCode == 0) then
        print '("    R02:  X=", F14.6, "  Y=", F14.6, "  Z=", F14.6, " km")', pos(1), pos(2), pos(3)
    else
        print *, "    ERROR: ", trim(errMsg)
        call PrintErrorInfo
    endif

    ! 测试 4: 查询 BeiDou 卫星 C06
    print *, "  Test 4: Query C06 at start epoch"
    pos = sp3%getPosition(2026, 5, 10, 0, 0, 0.0d0, "C06")
    if (errCode == 0) then
        print '("    C06:  X=", F14.6, "  Y=", F14.6, "  Z=", F14.6, " km")', pos(1), pos(2), pos(3)
    else
        print *, "    ERROR: ", trim(errMsg)
        call PrintErrorInfo
    endif

    ! 测试 5: 查询 Galileo 卫星 E02
    print *, "  Test 5: Query E02 at start epoch"
    pos = sp3%getPosition(2026, 5, 10, 0, 0, 0.0d0, "E02")
    if (errCode == 0) then
        print '("    E02:  X=", F14.6, "  Y=", F14.6, "  Z=", F14.6, " km")', pos(1), pos(2), pos(3)
    else
        print *, "    ERROR: ", trim(errMsg)
        call PrintErrorInfo
    endif

    ! 测试 6: 查询最后一个历元 G01 的位置
    print *, "  Test 6: Query G01 at last epoch"
    pos = sp3%getPosition(2026, 5, 11, 23, 45, 0.0d0, "G01")
    if (errCode == 0) then
        print '("    G01:  X=", F14.6, "  Y=", F14.6, "  Z=", F14.6, " km")', pos(1), pos(2), pos(3)
    else
        print *, "    ERROR: ", trim(errMsg)
        call PrintErrorInfo
    endif

    print *, ""

    ! ===== 输出统计信息 =====
    print *, "--- Statistics ---"
    print '("  Total epochs:      ", I0)', sp3%header%nEpoch
    print '("  Total satellites:  ", I0)', sp3%header%nSat
    print '("  GPS satellites:    ", I0)', count(sp3%header%sat(1:sp3%header%nSat)%sys == 'G')
    print '("  GLONASS satellites:", I0)', count(sp3%header%sat(1:sp3%header%nSat)%sys == 'R')
    print '("  BeiDou satellites: ", I0)', count(sp3%header%sat(1:sp3%header%nSat)%sys == 'C')
    print '("  Galileo satellites:", I0)', count(sp3%header%sat(1:sp3%header%nSat)%sys == 'E')
    print '("  Time span:         ", F10.3, " hours")', (sp3%header%nEpoch - 1) * sp3%header%dt / 3600.0d0
    print *, ""

    ! ===== 测试 freeData 方法 =====
    print *, ">>> Testing freeData() method ..."
    call sp3%freeData()
    if (allocated(sp3%data)) then
        print *, "    ERROR: data should have been deallocated!"
    else
        print *, "    Data successfully freed."
    endif
    print *, ""

    print *, "============================================"
    print *, "  Test completed successfully!"
    print *, "============================================"


end program test_sp3
