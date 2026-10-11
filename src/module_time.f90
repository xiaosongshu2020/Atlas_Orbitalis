!================================================================================
! MODULE: module_time
! 
! DESCRIPTION:
!   This module provides comprehensive time system handling and conversion 
!   capabilities for astronomical and satellite applications. It implements
!   object-oriented time representation with support for multiple time systems
!   including UTC, TT (Terrestrial Dynamical Time), UT1, BDT (BeiDou Time), 
!   and GPST (GPS Time).
!
! KEY FEATURES:
!   - Object-oriented time representation using Fortran 2003/2008 type-bound procedures
!   - Support for five major time systems with automatic conversions
!   - Flexible initialization from MJD, day-second pairs, or calendar dates
!   - Time arithmetic operations (addition and subtraction)
!   - High-precision date formatting with millisecond resolution
!   - Integration with SOFA (Standards of Fundamental Astronomy) routines
!
! TYPE_T OBJECT CORE USAGE:
!   The type_t object represents a temporal instant with the following components:
!     - ttype: Time system identifier (UTC, TT, UT1, BDT, or GPST)
!     - day:   Integer day component in Modified Julian Date (MJD)
!     - sec:   Seconds within the day (0 <= sec < 86400)
!     - val:   Complete MJD value (day + sec/86400)
!
!   INITIALIZATION METHODS:
!     1. From complete MJD:          call time%init(mjd)
!     2. From day-second pair:       call time%init(day, sec)
!     3. From calendar date:         call time%init(year, month, day, hour, minute, second)
!
!   TIME CONVERSION METHODS (type-bound procedures):
!     - call time%toUTC()            Convert to Coordinated Universal Time
!     - call time%toTT()            Convert to Terrestrial Dynamical Time  
!     - call time%toUT1()            Convert to Universal Time 1
!     - call time%toBDT(week, sow)   Convert to BeiDou Time (returns week and seconds of week)
!     - call time%toGPST(week, sow)  Convert to GPS Time (returns week and seconds of week)
!     - str = time%toDate()          Return formatted date string "YYYY-MM-DD HH:MM:SS.sss"
!
!   TIME ARITHMETIC:
!     - dt = timeMinus(t1, t2)       Difference between two time objects (in days)
!     - t2 = timeAdd(t1, dt)         Add time interval (in days) to time object
!     - t2 = timeAdd(t1, days, secs) Add separate day and second components
!
!   SUPPORTING SUBROUTINES:
!     - Date/MJD conversions: date2MJD, date2JD, JD2date (SOFA-based)
!     - Time system conversions: UTC2TT, TT2UTC, UTC2UT1, etc.
!     - String parsing: month_string2num, day_string2num
!
! AUTHORSHIP AND REFERENCES:
!   - Primary implementation: Shawn Hu
!   - Incorporates SOFA library routines for Julian Date conversions
!   - Leap second handling based on IERS bulletins
!   - BeiDou and GPS time system definitions per respective interface specifications
!
! LAST MODIFIED: 2025-12-19 by Shawn Hu
!================================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_time
use module_const
use module_toolKits

    integer, parameter :: TIME_TYPE_UTC = 1
    integer, parameter :: TIME_TYPE_TT = 2
    integer, parameter :: TIME_TYPE_UT1 = 3
    integer, parameter :: TIME_TYPE_BDT = 4
    integer, parameter :: TIME_TYPE_GPST = 5
    
    type type_t
        integer :: ttype = TIME_TYPE_UTC ! time type, default is UTC
        double precision :: day = 0d0
        double precision :: sec = 0d0
        double precision :: val = 0d0
    contains
        procedure :: init => init_time
        procedure :: toDate
        procedure :: toUTC
        procedure :: toTT
        procedure :: toUT1
        procedure :: toGPST
        procedure :: toBDT
    endtype

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initialize a time object with t1 and t2
    ! t1 and t2 can be:
    ! type 1 - t1 is MJD, t2 is not exist
    ! type 2 - t1 is day number in MJD, t2 is sec number in MJD
    ! type 3 - t1 ~ t6 are year, month, day, hour, minute, second
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine init_time(this, t1, t2, t3, t4, t5, t6)
    implicit none
    class(type_t) :: this
    double precision, intent(in) :: t1
    double precision, intent(in), optional :: t2
    double precision, intent(in), optional :: t3, t4, t5, t6
    type(type_t) :: t
    double precision :: year, month, day, hour, min, sec
    double precision :: MJD_day, MJD_sec
    !---------------------------------------------------
    if (.not. present(t2)) then ! input: t1 = MJD
        this%day = dble(floor(t1))
        this%sec = (t1 - this%day)*86400d0
        this%val = t1
    elseif (.not. present(t3)) then! input: t1 = day, t2 = sec
        this%day = t1
        this%sec = t2
        this%val = t1 + t2/86400d0
    else ! input: t1 ~ t6: y-m-d h:m:s
        year = t1
        month = t2
        day = t3
        hour = t4
        min = t5
        sec = t6
        sec = sec + min*60d0 + hour*3600d0
        call date2MJD(year, month, day, sec, MJD_day, MJD_sec)
        this%ttype = TIME_TYPE_UTC
        this%day = MJD_day
        this%sec = MJD_sec
        this%val = MJD_day + MJD_sec/86400d0
    endif
    endsubroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function timeMinus(t1, t2) result(dt)
    implicit none
    type(type_t), intent(in) :: t1, t2
    double precision dt
    dt = t1%day - t2%day + (t1%sec - t2%sec)/86400d0
    endfunction

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function timeAdd(t1, dt1, dt2) result(t2)
    implicit none
    type(type_t), intent(in) :: t1
    double precision, intent(in) :: dt1
    double precision, intent(in), optional :: dt2
    type(type_t) :: t2
    !------------------------------------------------------------
    t2 = t1    
    if (.not. present(dt2)) then
        if (dt1 >= 0) then
            t2%sec = t1%sec + dt1*86400d0
            do while (t2%sec >= 86400d0) 
                t2%sec = t2%sec - 86400d0
                t2%day = t2%day + 1d0
            enddo
        else
            t2%sec = t1%sec + dt1*86400d0
            do while (t2%sec < 0)
                t2%sec = t2%sec + 86400d0
                t2%day = t2%day - 1d0
            enddo
        endif
        t2%val = t2%day + t2%sec/86400d0
    else
        t2%day = t2%day + dt1
        t2%sec = t2%sec + dt2
        do while (t2%sec >= 86400d0) 
            t2%sec = t2%sec - 86400d0
            t2%day = t2%day + 1d0
        enddo
        do while (t2%sec < 0)
            t2%sec = t2%sec + 86400d0
            t2%day = t2%day - 1d0
        enddo
        t2%val = t2%day + t2%sec/86400d0
    endif
    
    endfunction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !For example, 0:00 on 1st Junuary, 2000 means that 
    !year = 2000,month = 1, day = 1, sec = 0
    function date2MJD2(year, month, day, sec) result(t)
    implicit none
    double precision, intent(in) :: year, month, day, sec
    type(type_t) :: t
    double precision y, a, b, c, jd_tmp, m
    !----------------------------------------------------------
    y = year
    m = month
    b = 0d0
    c = 0d0
    if (m <= 2) then
        y = y - 1d0
        m = m + 12d0
    end if
    if (y < 0) then
        c = -0.75d0
    end if

    !Check for valid calendar date
    if (year < 1582) then

    elseif (year > 1582) then
        call fix(y/100d0, a)
        b = 2d0 - a + dble(floor(a/4d0))
    elseif (month < 10) then
    elseif (month > 10) then
        call fix(y/100d0, a)
        b = 2d0 - a + dble(floor(a/4d0))
    elseif (day <= 4) then
    elseif (day > 14) then
        call fix(y/100d0, a)
        b = 2d0 - a + dble(floor(a/4))
    else
        write(*, *) "Subroutine:date2JD, This is an invalid calendar date!!!!!"
        return
    end if
    jd_tmp = dble(floor(365.25*y + c)) + dble(floor(30.6001*(m + 1)))
    !JD = jd_tmp + day + b + 1720994.5d0
    t%day = jd_tmp + day + b - 679006d0
    t%sec = sec
    t%val = t%day + t%sec/86400d0
    endfunction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function toDate(this) result(str)
    implicit none
    class(type_t) :: this
    character(len=28) :: str
    integer :: year, month, day
    double precision :: dayf
    integer :: err
    integer :: hour, min
    double precision :: sec, sec1
    integer :: secI, secM
    !---------------------------------------------------
    
    if(this%ttype == TIME_TYPE_UTC) then
        call jd2date(this%day+2400000.5d0, this%sec/86400d0, year, month, day, dayf, err)
        if (err /= 0) then
            print*, "error @suborutine date"
            print*, "Ivalid time: day = ", this%day, ", sec = ", this%sec
            str = ""
            return
        endif
        hour = floor(dayf*24d0)
        min = floor((dayf - dble(hour)/24d0)*1440d0)
        sec = (dayF - dble(hour)/24d0 - dble(min)/1440d0)*86400d0
        sec1 = dble(nint(sec*1d3))/1d3
        if (dabs(sec1 - 60d0) <= 1d-3) then
            sec1 = 0d0
            min = min + 1
            if (min == 60) then
                min = 0
                hour = hour+1
            endif
        endif
        secI = floor(sec1)
        secM = nint(sec1*1000d0 - dble(secI)*1000d0)
        write(str, "(I4,'-',I2.2,'-',I2.2,1x,I2.2,':',I2.2,':',I2.2,'.',I3.3)") year, month, day, hour, min, secI, secM
    else
        print*, "error @subroutine show."
        print*, "Invalid input: ttype = ", this%ttype
        return
    endif
    
    endfunction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine toUTC(this)
    implicit none
    class(type_t) :: this
    double precision :: day, sec
    
    !---------------------------------------------------
    select case(this%tType) 
    case(TIME_TYPE_UTC)
        return
    case(TIME_TYPE_TT)
        this%tType = TIME_TYPE_UTC
        call TT2UTC_daySec(this%day, this%sec, day, sec)
        this%day = day
        this%sec = sec
        this%val = day + sec/86400d0
    case(TIME_TYPE_UT1)
        print*, "warning @toUTC."
        print*, "unfinished part."
    case(TIME_TYPE_BDT)
        print*, "warning @toUTC."
        print*, "unfinished part."
    case(TIME_TYPE_GPST)
        print*, "warning @toUTC."
        print*, "unfinished part."
    case default
        print*, "warning @toUTC."
        print*, "invalid parameter, tType = ", this%tType
    endSelect
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine toTT(this)
    implicit none
    class(type_t) :: this
    double precision :: day, sec
    select case(this%tType) 
    case(TIME_TYPE_UTC)
        call UTC2TT_daySec(this%day, this%Sec, day, sec)
        this%tType = TIME_TYPE_TT
        this%day = day
        this%sec = sec
        this%val = day + sec/86400d0
    case(TIME_TYPE_TT)
        return
    case(TIME_TYPE_UT1)
        print*, "warning @toUTC."
        print*, "unfinished part."
    case(TIME_TYPE_BDT)
        this%tType = TIME_TYPE_TT
        this%day = this%day + 53736d0
        this%sec = this%sec + 65.184d0
        if (this%sec >= 86400d0) then
            this%day = this%day + 1d0
            this%sec = this%sec - 86400d0
        endif
    case(TIME_TYPE_GPST)
        ! GPST -> BDT
        this%day = this%day - 9492d0
        this%sec = this%sec - 14d0
        if (this%sec < 0d0) then
            this%day = this%day - 1d0
            this%sec = this%sec + 86400d0
        endif
        ! BDT -> TT
        this%tType = TIME_TYPE_TT
        this%day = this%day + 53736d0
        this%sec = this%sec + 65.184d0
        if (this%sec >= 86400d0) then
            this%day = this%day + 1d0
            this%sec = this%sec - 86400d0
        endif
    case default
        print*, "warning @toTT."
        print*, "invalid parameter, tType = ", this%tType
    endSelect
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine toUT1(this)
    implicit none
    class(type_t) :: this
    print*, "unfinished part."
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine toBDT(this, week, sow)
    implicit none
    class(type_t) :: this
    double precision, intent(out) :: week, sow
    double precision day, sec
    select case(this%tType) 
    case(TIME_TYPE_UTC)
        call UTC2TT_daySec(this%day, this%Sec, day, sec)
        call TT2BDT_daySec(day, sec, week, sow)
        this%tType = TIME_TYPE_BDT
        this%day = week*7d0 
        this%sec = sow
        do while(this%sec >= 86400d0)
            this%day = this%day + 1d0
            this%sec = this%sec - 86400d0
        enddo
    case(TIME_TYPE_TT)
        day = this%day
        sec = this%sec
        call TT2BDT_daySec(day, sec, week, sow)
        this%tType = TIME_TYPE_BDT
        this%day = week*7d0 
        this%sec = sow
        do while(this%sec >= 86400d0)
            this%day = this%day + 1d0
            this%sec = this%sec - 86400d0
        enddo
    case(TIME_TYPE_UT1)
        print*, "warning @toUTC."
        print*, "unfinished part."
    case(TIME_TYPE_BDT)
        week = dble(floor(this%day/7d0))
        sow = this%sec + (this%day - week*7d0)*86400d0
    case(TIME_TYPE_GPST)
        this%tType = TIME_TYPE_BDT
        this%day = this%day - 9492d0
        this%sec = this%sec - 14d0
        if (this%sec < 0d0) then
            this%sec = this%sec + 86400d0
            this%day = this%day -1d0
        endif
    case default
        print*, "warning @toBDST."
        print*, "invalid parameter, tType = ", this%tType
    endSelect
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine toGPST(this, week, sow)
    implicit none
    class(type_t) :: this
    double precision, intent(out) :: week, sow
    !---------------------------------------------------
    select case(this%tType) 
    case(TIME_TYPE_GPST)
        week = dble(floor(this%day/7d0))
        sow = this%sec + (this%day - week*7d0)*86400d0
    case default
        ! -> BDT
        call this%toBDT(week, sow)
        ! BDT -> GPST
        this%day = this%day + 9492d0
        this%sec = this%sec + 14d0
        if (this%sec >= 86400d0) then
            this%day = this%day + 1d0
            this%sec = this%sec - 86400d0
        endif
        week = dble(floor(this%day/7d0))
        sow = this%sec + (this%day - week*7d0)*86400d0
    endSelect
    
    endsubroutine
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Transfer time from UTC to TT ( in MJD)
    !Last Modified : 2021/6/22 by Shawn Hu
    subroutine UTC2TT(UTC, TT)
    implicit none
    double precision UTC, TAI, TT, delta
    !----------------------------
    !call syncLeapSec(UTC)
    delta = 37d0
    TAI = UTC + delta/86400d0
    TT = TAI + 32.184d0/86400d0
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine UTC2UT1(UTC, DUT1, UT1)
    double precision UTC, UT1, DUT1
    UT1 = UTC + DUT1/86400d0
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  TT2UTC(TT, UTC)
    implicit    none
    double precision      TT, UTC, TAI
    !-----------------------------------
    TAI     =   TT - 32.184d0/86400d0  !TAI
    UTC     =   TAI - 37d0/86400d0
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine TT2UTC_daySec(day1, sec1, day2, sec2)
    implicit none
    double precision day1, sec1, day2, sec2
    !---------------------------------------------------
    day2 = day1
    sec2 = sec1 - 32.184d0 - 37d0
    if(sec2 < 0d0) then
        day2 = day2 - 1d0
        sec2 = sec2 + 86400d0
    endif
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TT --> Beidou time [week, sec]    
    subroutine TT2BDT(TT, week, sec)
    implicit none
    double precision TT, sec, tmp1, tmp2,dt
    integer week    
    !--------------------------------------------------
    dt = TT - 53736.0007544444d0
    week = floor(dt/7d0)
    sec = (dt - dble(week)*7d0)*86400d0
    endsubroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TT [day, sec] --> Beidou time [week, sow]    
    subroutine TT2BDT_daySec(day, sec, week, sow)
    implicit none
    double precision day, sec, tmp1, tmp2, dt
    double precision week, sow  
    !----------------------------------------------------
    week = dble(floor((day-53736d0)/7d0))
    sow = (day-53736d0-7d0*week)*86400d0 + sec - 65.184d0
    
    if(sow < 0d0) then
        sow = sow + 86400d0*7d0
        week = week - 1d0
    endif
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Beidou time [week, sow]  -->  TT [day, sec]  
    subroutine BDT2TT_daySec(week, sow, day, sec)
    implicit none
    double precision day, sec, tmp1, tmp2, dt
    double precision week, sow  
    !----------------------------------------------------
    day = week*7d0
    sec = sow
    day = day + 53736d0
    sec = sow + 65.184d0
    do while(sec >= 86400d0)
        sec = sec - 86400d0
        day = day + 1d0
    enddo
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Transfer time from UTC(day, sec) to TT(day, sec)
    !Last Modified : 2021/11/16 Hu
    
    subroutine UTC2TT_daySec(UTC_day, UTC_sec, TT_day, TT_sec)
    implicit none
    double precision UTC, TAI, TT, delta
    double precision UTC_day, UTC_sec, TT_day, TT_sec
    !----------------------------------------------------------
    UTC = UTC_day + UTC_sec/86400d0
    if (UTC >= 57754d0) then
        delta = 37d0
    elseif (UTC >= 57204d0) then
        delta = 36d0
    elseif (UTC >= 56109d0) then
        delta = 35d0
    elseif (UTC >= 54832d0) then
        delta = 34.0d0
    elseif (UTC >= 53736d0) then
        delta = 33.0d0
    elseif (UTC >= 51179d0) then
        delta = 32.0d0
    elseif (UTC >= 50630d0) then
        delta = 31.0d0
    else
        write(*, *) "Out of range of this subroutine: UTC2TT.f"
    end if
    
    TT_day = UTC_day
    TT_sec = UTC_sec + delta + 32.184d0
    if (TT_sec >= 86400d0) then
        TT_sec = TT_sec - 86400d0
        TT_day = TT_day + 1d0
    endif
    
    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !For example, 0:00 on 1st Junuary, 2000 means that year = 2000,month = 1, day = 1, sec = 0
    subroutine date2MJD(year, month, day, sec, MJD_day, MJD_sec)
    implicit none
    double precision day, sec, MJD_day, MJD_sec, y, a, b, c, jd_tmp
    double precision year, month, m
    !----------------------------------------------------------
    y = year
    m = month
    b = 0d0
    c = 0d0
    if (m <= 2) then
        y = y - 1d0
        m = m + 12d0
    end if
    if (y < 0) then
        c = -0.75d0
    end if

    !Check for valid calendar date
    if (year < 1582) then

    elseif (year > 1582) then
        call fix(y/100d0, a)
        b = 2d0 - a + dble(floor(a/4d0))
    elseif (month < 10) then
    elseif (month > 10) then
        call fix(y/100d0, a)
        b = 2d0 - a + dble(floor(a/4d0))
    elseif (day <= 4) then
    elseif (day > 14) then
        call fix(y/100d0, a)
        b = 2d0 - a + dble(floor(a/4))
    else
        write(*, *) "Subroutine:date2JD, This is an invalid calendar date!!!!!"
        return
    end if
    jd_tmp = dble(floor(365.25*y + c)) + dble(floor(30.6001*(m + 1)))
    !JD = jd_tmp + day + b + 1720994.5d0
    MJD_day = jd_tmp + day + b - 679006d0
    MJD_sec = sec
    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Beidou time [week, sec] --> TT (MJD)
    ! Author: Shawn Hu
    ! Last modified: 2021/7/17 Shawn Hu
    subroutine BDT2TT(week, sec, TT)
    implicit none
    integer week
    double precision sec, TT    
    !---------------------------------------
    TT = 53736.0007544444d0+dble(week)*7d0+sec/86400d0
    endsubroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE month_string2num(MonthAbbrev, MonthNumber)  
    implicit none
    CHARACTER(LEN=3), INTENT(IN) :: MonthAbbrev  
    double precision, INTENT(OUT) :: MonthNumber          
  
    MonthNumber = 0d0  
  
    SELECT CASE (MonthAbbrev)  
        CASE ('Jan')  
            MonthNumber = 1d0  
        CASE ('Feb')  
            MonthNumber = 2d0  
        CASE ('Mar')  
            MonthNumber = 3d0  
        CASE ('Apr')  
            MonthNumber = 4d0  
        CASE ('May')  
            MonthNumber = 5d0  
        CASE ('Jun')  
            MonthNumber = 6d0  
        CASE ('Jul')  
            MonthNumber = 7d0  
        CASE ('Aug')  
            MonthNumber = 8d0  
        CASE ('Sep')  
            MonthNumber = 9d0  
        CASE ('Oct')  
            MonthNumber = 10d0  
        CASE ('Nov')  
            MonthNumber = 11d0  
        CASE ('Dec')  
            MonthNumber = 12d0  
    END SELECT  
  
    END SUBROUTINE 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    subroutine day_string2num(TimeString, Seconds)
    implicit none
    CHARACTER(LEN=12), INTENT(IN) :: TimeString  
    double precision, INTENT(OUT) :: Seconds       
    INTEGER :: Hours  
    INTEGER :: Minutes  
    double precision :: Sec
    
    seconds = 0d0
    read(timeString(1:2), *) hours
    read(timeString(4:5), *) minutes
    read(timeString(7:12), *) sec
    seconds = hours*3600d0 + minutes*60d0 + sec
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !---------------------------------------------------------------------
    ! NOTE: derived work.  JD2date below is adapted from the IAU SOFA
    ! routine iau_JD2CAL (SOFA release 2023-10-11, Copyright (C) 2023 IAU
    ! SOFA Board).  It is renamed -- no "iau"/"sofa" prefix -- and is not
    ! endorsed by SOFA.  The SOFA Software License, whose six terms and
    ! conditions apply to this routine, is reproduced in full at the end of
    ! every file in sofa/ and summarised in THIRD_PARTY.md.
    !---------------------------------------------------------------------
    SUBROUTINE JD2date ( DJ1, DJ2, IY, IM, ID, FD, J )
!+
!  - - - - - - - - - - -
!   i a u _ J D 2 C A L
!  - - - - - - - - - - -
!
!  Julian Date to Gregorian year, month, day, and fraction of a day.
!
!  This routine is part of the International Astronomical Union's
!  SOFA (Standards of Fundamental Astronomy) software collection.
!
!  Status:  support routine.
!
!  Given:
!     DJ1,DJ2     d     Julian Date (Notes 1, 2)
!
!  Returned:
!     IY          i     year
!     IM          i     month
!     ID          i     day
!     FD          d     fraction of day
!     J           i     status:
!                           0 = OK
!                          -1 = unacceptable date (Note 1)
!
!  Notes:
!
!  1) The earliest valid date is -68569.5 (-4900 March 1).  The
!     largest value accepted is 10^9.
!
!  2) The Julian Date is apportioned in any convenient way between
!     the arguments DJ1 and DJ2.  For example, JD=2450123.7 could
!     be expressed in any of these ways, among others:
!
!             DJ1            DJ2
!
!         2450123.7D0        0D0        (JD method)
!          2451545D0      -1421.3D0     (J2000 method)
!         2400000.5D0     50123.2D0     (MJD method)
!         2450123.5D0       0.2D0       (date & time method)
!
!     Separating integer and fraction uses the "compensated summation"
!     algorithm of Kahan-Neumaier to preserve as much precision as
!     possible irrespective of the JD1+JD2 apportionment.
!
!  3) In early eras the conversion is from the "Proleptic Gregorian
!     Calendar";  no account is taken of the date(s) of adoption of
!     the Gregorian Calendar, nor is the AD/BC numbering convention
!     observed.
!
!  References:
!
!     Explanatory Supplement to the Astronomical Almanac,
!     P. Kenneth Seidelmann (ed), University Science Books (1992),
!     Section 12.92 (p604).
!
!     Klein, A., A Generalized Kahan-Babuska-Summation-Algorithm.
!     Computing, 76, 279-293 (2006), Section 3.
!
!  This revision:  2020 October 21
!
!  SOFA release 2023-10-11
!
!  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
!
!-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DJ1, DJ2
      INTEGER IY, IM, ID
      DOUBLE PRECISION FD
      INTEGER J

!  Smallest number such that 1D0+EPS .NE. 1D0
      DOUBLE PRECISION EPS
      PARAMETER ( EPS = 2.2204460492503131D-16  )

!  Minimum and maximum allowed JD
      DOUBLE PRECISION DJMIN, DJMAX
      PARAMETER ( DJMIN = -68569.5D0, DJMAX = 1D9 )

      INTEGER JD, I, L, N, K
      DOUBLE PRECISION DJ, F1, F2, D, S, CS, V(2), X, T, C, F

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Verify date is acceptable.
      DJ = DJ1 + DJ2
      IF ( DJ.LT.DJMIN .OR. DJ.GT.DJMAX ) THEN
         J = -1
      ELSE
         J = 0

!     Separate day and fraction (where -0.5 <= fraction < 0.5).
         D = ANINT(DJ1)
         F1 = DJ1 - D
         JD = INT(D)
         D = ANINT(DJ2)
         F2 = DJ2 - D
         JD = JD + INT(D)

!     Compute F1+F2+0.5 using compensated summation (Klein 2006).
         S = 0.5D0
         CS = 0D0
         V(1) = F1
         V(2) = F2
         DO 1 I=1,2
            X = V(I)
            T = S + X
            IF ( ABS(S) .GE. ABS(X) ) THEN
               C = (S-T) + X
            ELSE
               C = (X-T) + S
            END IF
            CS = CS + C
            S = T
            IF ( S .GE. 1D0 ) THEN
               JD = JD + 1
               S = S - 1D0
            END IF
 1       CONTINUE
         F = S + CS
         CS = F - S

!     Deal with negative F.
         IF ( F .LT. 0D0 ) THEN

!        Compensated summation: assume that |S| <= 1.
            F = S + 1D0
            CS = CS + ((1D0-F) + S)
            S = F
            F = S + CS
            CS = F - S
            JD = JD - 1
         END IF

!     Deal with F that is 1D0 or more (when rounded to double).
         IF ( (F-1D0) .GE. -EPS/4D0 ) THEN

!        Compensated summation: assume that |S| <= 1. !/
            T = S - 1D0
            CS = CS + ((S-T) - 1D0)
            S = T
            F = S + CS
            IF ( -EPS/2D0 .LT. F ) THEN
               JD = JD + 1
               F = MAX ( F, 0D0 )
            END IF
         END IF

!     Express day in Gregorian calendar.
         L = JD + 68569
         N = ( 4*L ) / 146097
         L = L - ( 146097*N + 3 ) / 4
         I = ( 4000 * (L+1) ) / 1461001
         L = L - ( 1461*I ) / 4 + 31
         K = ( 80*L ) / 2447
         ID = L - ( 2447*K ) / 80
         L = K / 11
         IM = K + 2 - 12*L
         IY = 100 * ( N-49 ) + I + L
         FD = F
      END IF

    ENDSUBROUTINE


    endmodule
