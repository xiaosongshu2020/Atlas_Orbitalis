module module_sp3

    use ErrorHandler
    implicit none

    !! >> ------------------------------------------------- DEFINE PARAMETERS AND ENUMERATORS
    integer,parameter:: maxSatellite = 500
    integer,parameter:: maxEpoch = 300

    !! >> ------------------------------------------------- DEFINE ELEMENTARY TYPES
    type sp3_Gregorian
        integer:: yr,mo,day,hr,mi
        double precision:: sec
    end type sp3_Gregorian

    type sp3_mjd
        integer:: day
        double precision:: frac
    end type sp3_mjd

    type sp3_WeekSec
        integer:: week
        double precision:: sec
    end type sp3_WeekSec

    type sp3_t0
        type(sp3_Gregorian):: cal
        type(sp3_mjd):: mjd
        type(sp3_WeekSec):: week
    end type sp3_t0

    type sp3_descriptor
        character(len=3):: oType
        character(len=4):: Agency
        character(len=5):: uData,Frame
    end type sp3_descriptor

    type sp3_SatID
        integer:: SatNum=0
        character(len=1):: sys=' '
    end type sp3_SatID

    type sp3_pData
        integer:: posAccExp(3),clkAccExp
!        integer:: Flag(4)=0
        character(len=1):: Flag(4)=' '
        double precision:: pos(3)=0d0, clk=999999d0
        logical:: isClkEvent = .false., isClkPred = .false., isManeuver = .false., isOrbitPred = .false.
        type(sp3_SatID):: sat
        type(sp3_Gregorian):: t
    end type sp3_pData

    !! >> ------------------------------------------------- DEFINE MAJOR TYPES
    type sp3c_Header
        integer:: nEpoch,nSat,AccExp(maxSatellite)=0
        double precision:: dt,pvAccBase,clkAccBase
        character(len=1):: ver,pvType
        character(len=2):: fType
        character(len=3):: tSys
        type(sp3_t0):: t0
        type(sp3_descriptor):: meta
        type(sp3_SatID):: sat(maxSatellite)
    end type sp3c_Header

    !! >> ------------------------------------------------- DEFINE THE SP3 OBJECT (CLASS)
    type type_sp3
        type(sp3c_Header) :: header
        type(sp3_pData), allocatable :: data(:,:)
    contains
        procedure :: load        => load_sp3_file
        procedure :: freeData    => free_sp3_data
        procedure :: getPosition => get_satellite_position
    end type type_sp3

    !! >> ------------------------------------------------- DEFINE OPERATIONS ON SP3 HEADER AND DATA

    contains
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        type(sp3c_Header) function read_sp3cHeader(sp3cUnit)
            use String
            implicit none
            ! in/out
            integer,intent(in):: sp3cUnit
            ! local
            integer:: i,j
            character(len=4):: subName='r3ch'
            character(len=80):: rec
            type(sp3c_Header):: header

            ! code
            call checkin(subName)

            ! start reading sp3c file line by line
            ! > Lines 1 & 2
            read(sp3cUnit,'(a80)',err=99) rec
            read(rec,1,err=98) header%ver,header%pvType,header%t0%cal%yr,header%t0%cal%mo,header%t0%cal%day,header%t0%cal%hr, &
                            &  header%t0%cal%mi,header%t0%cal%sec,header%nEpoch,header%meta%uData,header%meta%Frame,          &
                            &  header%meta%oType,header%meta%Agency
            read(sp3cUnit,'(a80)',err=99) rec
            read(rec,2,err=98) header%t0%week%week,header%t0%week%sec,header%dt,header%t0%mjd%day,header%t0%mjd%frac

            if (lower(header%ver).ne.'c') then
                call RaiseError(62003,subName)
                return
            endif
            if (lower(Header%pvType).eq.'v') then
                call RaiseWarning('Velocity and clock rate-of-change are not yet supported.',subName)
            endif

            ! > Lines 3 - 12
            ! ȡ͵һPRN
            read(sp3cUnit,'(a80)',err=99) rec
            read(rec,3,err=98) header%nSat,(header%sat(j)%sys,header%sat(j)%SatNum,j=1,17)
            do i = 2,10
                read(sp3cUnit,'(a80)',err=99) rec
                read(rec,4,err=98) (header%sat(j)%sys,header%sat(j)%SatNum,j=1+(i-1)*17,17*i)
            enddo

            ! > Lines 13 - 22
            do i = 1,10
                read(sp3cUnit,'(a80)',err=99) rec
                read(rec,8,err=98) (header%AccExp(j),j=1+(i-1)*17,17*i)
            enddo

            ! > Line 23
            read(sp3cUnit,'(a80)',err=99) rec
            read(rec,13,err=98) header%fType,header%tSys

            ! > Line 24, skip
            read(sp3cUnit,'(a80)',err=99) rec

            ! > Line 25
            read(sp3cUnit,'(a80)',err=99) rec
            read(rec,15,err=98) header%pvAccBase,header%clkAccBase

            ! > Lines 16 - 22, skip
            do i = 26,32
                read(sp3cUnit,'(a80)',err=99) rec
            enddo

            ! normal exit
            read_sp3cHeader = header
            call checkout(subName)
            return

            ! restricted area
 99          call RaiseError(62001,subName)
            read_sp3cHeader = header
            return
 98          call RaiseError(62002,subName)
            read_sp3cHeader = header
            return

            ! format strings
  1         format(1x,1a,1a,i4,4(1x,i2),1x,f11.8,1x,i7,1x,a5,1x,a5,1x,a3,1x,a4)
  2         format(3x,i4,1x,f15.8,1x,f14.8,1x,i5,1x,f16.13)
  3         format(3x,i3,3x,17(a1,i2.2))    !!!!!!! revised by Shawn Hu
  4         format(9x,17(a1,i2.2))
  8         format(9x,17(i3))
 13         format(3x,a2,4x,a3)
 15         format(3x,f10.7,1x,f12.9)
        end function read_sp3cHeader

        ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function read_sp3cPData(sp3cUnit,nEpoch,nSat) result (pData)
            use String
            implicit none
            ! in/out
            integer,intent(in):: sp3cUnit,nEpoch,nSat
            type(sp3_pData):: pData(nEpoch,nSat)
            ! local
            integer:: i,j
            character(len=4):: subName = 'r3cp'
            character(len=80):: card

            ! code
            call checkin(subName)

            do i = 1,nEpoch
                read(sp3cUnit,'(a80)',err=99) card
                read(card,'(3x,i4,4(1x,i2),1x,f11.8)',err=99) pData(i,1)%t%yr,pData(i,1)%t%mo,pData(i,1)%t%day, &
                                                            & pData(i,1)%t%hr,pData(i,1)%t%mi,pData(i,1)%t%sec
                do j = 1,nSat
                    pData(i, j)%t = pData(i, 1)%t
                    read(sp3cUnit,'(a80)',err=99) card
                    read(card,24,err=98) pData(i,j)%sat%sys,pData(i,j)%sat%SatNum,pData(i,j)%pos(1:3),pData(i,j)%clk, &
                                      &  pData(i,j)%posAccExp(1:3),pData(i,j)%clkAccExp,pData(i,j)%Flag(1:4)
                    pData(i,j)%isClkEvent  = pData(i,j)%Flag(1).eq.'E'
                    pData(i,j)%isClkPred   = pData(i,j)%Flag(2).eq.'P'
                    pData(i,j)%isManeuver  = pData(i,j)%Flag(3).eq.'M'
                    pData(i,j)%isOrbitPred = pData(i,j)%Flag(4).eq.'P'
                enddo
            enddo

            ! normal exit
            call checkout(subName)
            return
            ! restricted area
 99         call RaiseError(62001,subName)
            return
 98         call RaiseError(62002,subName)
            return
            ! format strings
 24         format(1x,a1,i2.2,4(f14.6),3(1x,i2),1x,i3,1x,2(a1),2x,2(a1))

        end function read_sp3cPData
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! return satellite's NO. according its PRN
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function PRN2Num(header, PRN) result(satNO)
        implicit none
        type(sp3c_Header), intent(in) :: header
        character*3, intent(in) :: PRN
        integer :: satNO
        character*1 sys
        integer num, i
        !-------------------------------------------------------------------------
        read(PRN, "(A1,I2.2)") sys, num
        do i = 1, header%nSat
            if (header%sat(i)%sys == sys .and. header%sat(i)%SatNum == num) then
                satNO = i
                return
            endif
        enddo
        print*, "ERROR: failed to find satellite NO. @subroutine PRN2Num. PRN =", PRN
        STOP
        endfunction

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! >> ------------------------------------------------- TYPE_SP3 METHODS
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! load_sp3_file: 打开 SP3c 文件，读取头信息和位置数据
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine load_sp3_file(this, filename)
            implicit none
            ! in/out
            class(type_sp3), intent(inout) :: this
            character(len=*), intent(in) :: filename
            ! local
            integer :: unit, ioStat
            character(len=4) :: subName = 'load'

            call checkin(subName)

            ! 获取空闲文件单元
            unit = 101

            ! 打开文件
            open(unit=unit, file=filename, status='old', action='read', iostat=ioStat)
            if (ioStat /= 0) then
                call RaiseError(4005, subName)
                call checkout(subName)
                return
            endif

            ! 读取文件头
            this%header = read_sp3cHeader(unit)
            if (errCode > 0) then
                close(unit)
                call checkout(subName)
                return
            endif

            ! 分配数据数组
            if (allocated(this%data)) deallocate(this%data)
            allocate(this%data(this%header%nEpoch, this%header%nSat))

            ! 读取位置数据
            this%data = read_sp3cPData(unit, this%header%nEpoch, this%header%nSat)
            if (errCode > 0) then
                close(unit)
                deallocate(this%data)
                call checkout(subName)
                return
            endif

            ! 关闭文件
            close(unit)

            call checkout(subName)
            return
        end subroutine load_sp3_file

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! free_sp3_data: 释放 SP3 数据内存
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine free_sp3_data(this)
            implicit none
            ! in/out
            class(type_sp3), intent(inout) :: this
            character(len=4) :: subName = 'free'

            call checkin(subName)

            if (allocated(this%data)) then
                deallocate(this%data)
            endif

            call checkout(subName)
            return
        end subroutine free_sp3_data

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! get_satellite_position: 获取指定时刻指定卫星的位置
        ! 输入: 时间(格里高利历) + PRN号(如"G01")
        ! 返回: pos(3) 位置向量 (km)
        ! 说明: 通过起始时间和步长快速定位历元索引
        !       如果时间不在数据点上，报警告并返回附近数据
        !       如果时间超出数据范围，报错
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function get_satellite_position(this, year, month, day, hour, min, sec, PRN) result(pos)
            implicit none
            ! in/out
            class(type_sp3), intent(in) :: this
            integer, intent(in) :: year, month, day, hour, min
            double precision, intent(in) :: sec
            character(len=*), intent(in) :: PRN
            double precision :: pos(3)
            ! local
            integer :: satIdx, epochIdx
            double precision :: targetSec, epochSec, dt_sec, diff
            character(len=4) :: subName = 'gpos'

            call checkin(subName)

            ! 初始化返回值为0
            pos = 0.0d0

            ! 检查数据是否已加载
            if (.not. allocated(this%data)) then
                call RaiseError(0016, subName)
                call checkout(subName)
                return
            endif

            ! 查找卫星索引
            satIdx = PRN2Num(this%header, PRN)

            ! 计算目标时间相对于起始历元的秒偏移量
            call time_diff(this%header%t0%cal, year, month, day, hour, min, sec, targetSec)

            ! 通过步长快速定位历元索引
            dt_sec = this%header%dt
            epochIdx = nint(targetSec / dt_sec) + 1

            ! 检查是否超出数据范围
            if (epochIdx < 1 .or. epochIdx > this%header%nEpoch) then
                call RaiseError(0015, subName)
                call checkout(subName)
                return
            endif

            ! 检查时间是否在数据点上（允许1秒的容差）
            epochSec = (epochIdx - 1) * dt_sec
            diff = abs(targetSec - epochSec)
            if (diff > 1.0d0) then
                call RaiseWarning('Requested time is not exactly on an epoch. Returning nearest data.', subName)
            endif

            ! 返回位置
            pos(1:3) = this%data(epochIdx, satIdx)%pos(1:3)

            call checkout(subName)
            return
        end function get_satellite_position

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! time_diff: 计算两个格里高利时间之间的秒差 (t2 - t1)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine time_diff(t1, y2, m2, d2, h2, mi2, s2, diff_sec)
            implicit none
            ! in
            type(sp3_Gregorian), intent(in) :: t1
            integer, intent(in) :: y2, m2, d2, h2, mi2
            double precision, intent(in) :: s2
            ! out
            double precision, intent(out) :: diff_sec
            ! local
            integer :: jd1, jd2
            double precision :: sec1, sec2

            ! 将 t1 转换为儒略日和秒
            jd1 = greg2jd(t1%yr, t1%mo, t1%day)
            sec1 = t1%hr * 3600.0d0 + t1%mi * 60.0d0 + t1%sec

            ! 将 t2 转换为儒略日和秒
            jd2 = greg2jd(y2, m2, d2)
            sec2 = h2 * 3600.0d0 + mi2 * 60.0d0 + s2

            ! 计算秒差
            diff_sec = (jd2 - jd1) * 86400.0d0 + (sec2 - sec1)

            return
        end subroutine time_diff

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! greg2jd: 格里高利历转儒略日 (只计算日期部分，不含时间)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer function greg2jd(y, m, d)
            implicit none
            integer, intent(in) :: y, m, d
            integer :: a, b, yy, mm

            if (m <= 2) then
                yy = y - 1
                mm = m + 12
            else
                yy = y
                mm = m
            endif

            a = yy / 100
            b = 2 - a + a / 4

            greg2jd = int(365.25d0 * (yy + 4716)) + int(30.6001d0 * (mm + 1)) + d + b - 1524

            return
        end function greg2jd

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module module_sp3
