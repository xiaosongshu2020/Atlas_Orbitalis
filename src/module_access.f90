    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! module_access - visibility links and observation data
    !
    ! type_access represents a link between two objects:
    !   satellite <-> ground station, or satellite <-> satellite,
    ! together with the observation data generated over that link.
    !
    ! Whether a link is established at a given epoch is decided by
    ! the constraints stored in the access object:
    !   receiver constraints     : minEle (maxEle reserved)
    !   link constraints         : solarDistance, lunarDistance,
    !                              shadow, EarthBlockage, MoonBlockage
    !   transmitter constraints  : pointTo, t_coverage, beamHalfAngle
    ! The `force` argument of the generation routines bypasses all
    ! visibility constraints; POSITION_VELOCITY/POSITION data are always exempt.
    !
    ! Observation generation routines:
    !   calObsData1 : satellite <-> ground station
    !   calObsData2 : satellite <-> satellite
    !   calObsData3 : satellite only (state at observation epochs)
    !
    ! Module switches:
    !   key_light_time_delay         : light-time delay correction
    !   key_atmospherical_refraction : tropospheric refraction correction
    !
    ! Unit convention: access%dataTable is stored in NORMALIZED units
    ! (angles in radians); file read/save converts to/from SI
    ! (degrees, km, km/s) at the boundary.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    module module_access
    use module_station
    use module_satellite
    use module_atmosphere
    use module_toolKits
    use module_matrix

    !===========================================================================
    ! OBSERVATION DATA TYPES
    !===========================================================================
    
    integer, parameter :: OBS_TYPE_RA_DEC            = 11  ! Right ascension and declination
    integer, parameter :: OBS_TYPE_AZIMUTH_ELEVATION = 12  ! Azimuth and elevation
    integer, parameter :: OBS_TYPE_RANGING           = 13  ! Range-only measurements
    integer, parameter :: OBS_TYPE_RANGERATE         = 14  ! Range-rate-only measurements
    integer, parameter :: OBS_TYPE_POSITION          = 15  ! Position-only measurements
    integer, parameter :: OBS_TYPE_POSITION_VELOCITY = 16  ! Position and velocity

    type :: type_access
        
        !---------------- link identity -------------------------------------
        integer :: number = 0      ! access ID
        character*100 :: name = '' ! access name
        integer :: fromType = 0    ! source object type: OBJECT_TYPE_SATELLITE or OBJECT_TYPE_STATION
        integer :: fromNum = 0     ! source object number
        character*100 :: fromName = '' ! source object name (the transmitter of the link)
        integer :: toType = 0      ! target object type: OBJECT_TYPE_SATELLITE or OBJECT_TYPE_STATION
        integer :: toNum = 0       ! target object number
        character*100 :: toName = ''   ! target object name (the receiver of the link)
        integer :: obsType = 0     ! observation type: OBS_TYPE_* constant (11-16); determines npe
        integer :: nData = 0       ! number of observation epochs actually generated
        integer :: npe = 2         ! data dimension per epoch; auto-set by setNpe() from obsType
        
        !---------------- generated observation data ------------------------
        type(type_t), allocatable :: timeTable(:)          ! epochs (UTC) of the generated data, size nData
        double precision, allocatable :: dataTable(:, :)   ! observation data (nData x npe), normalized units
        
        !---------------- observation errors ---------------------------------
        double precision :: bias(10) = 0d0   ! systematic bias per observable column of dataTable
        double precision :: errStd(10) = 0d0 ! noise std-dev per observable column of dataTable
        
        !---------------- receiver constraints (satellite-station links) ----
        double precision :: minEle = 10d0 ! minimum elevation of the target above the station horizon (deg)
        double precision :: maxEle = 90d0 ! maximum elevation (deg); reserved, not used in the visibility filter yet
        
        !---------------- link constraints -----------------------------------
        double precision :: solarDistance = 0d0 ! min angular separation between target and Sun (deg); <=0 disables
        double precision :: lunarDistance = 0d0 ! min angular separation between target and Moon (deg); <=0 disables
        logical :: shadow = .True.              ! skip epochs when the satellite is inside the Earth shadow
        logical :: EarthBlockage = .True.       ! skip epochs when the Earth blocks the line of sight (sat-sat links)
        LOGICAL :: MoonBlockage = .False.       ! skip epochs when the Moon blocks the line of sight (sat-sat links)

        !---------------- transmitter constraints (satellite-satellite links) -
        integer :: pointTo = 1
        ! reference pointing of the transmitter antenna:
        !   1: toward the Earth center
        !   2: toward the receiver satellite
        logical :: t_coverage = .False.            ! enable the beam coverage constraint?
        double precision :: beamHalfAngle = 180d0  ! half-angle of the transmitter beam cone (deg);
                                                   !   the receiver direction must stay within this
                                                   !   angle of the reference pointing; 180 = no limit

        
    contains
        procedure :: init => initAccess
        procedure :: readData => readAccessData
        procedure :: saveData => saveAccessData
        procedure :: show => showAccess ! display access information
    endtype
    
    logical :: key_light_time_delay = .true. ! light-time delay correction
    logical :: key_atmospherical_refraction = .false. ! atmospherical refraction correction
    logical :: key_clock_error = .true.
    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize an access object with basic parameters
    ! 
    ! This procedure sets the primary properties of an access relationship
    ! and determines the appropriate number of parameters per epoch (npe)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine initAccess(this, number, name, fromType, fromNum, &
        &   toType, toNum, obsType, minEle, solarDistance, lunarDistance, &
        &   fromName, toName)
    implicit none
    class(type_access), intent(inout) :: this
    integer, intent(in), optional :: number               ! Access ID number
    character(len=*), intent(in), optional :: name        ! Name of the access
    integer, intent(in), optional :: fromType             ! Source type
    integer, intent(in), optional :: fromNum              ! Source ID number
    integer, intent(in), optional :: toType               ! Target type
    integer, intent(in), optional :: toNum                ! Target ID number
    character(len=*), intent(in), optional :: fromName    ! Name of the transmitter
    character(len=*), intent(in), optional :: toName      ! Name of the receiver
    integer, intent(in), optional :: obsType              ! Observation type
    double precision, intent(in), optional :: minEle                ! Minimum elevation angle (deg)
    double precision, intent(in), optional :: solarDistance
    double precision, intent(in), optional :: lunarDistance
        
    ! Set access properties if provided
    if (present(number)) this%number = number
    if (present(name)) this%name = name
    if (present(fromType)) this%fromType = fromType
    if (present(fromNum)) this%fromNum = fromNum
    if (present(toType)) this%toType = toType
    if (present(toNum)) this%toNum = toNum
    if (present(fromName)) this%fromName = fromName
    if (present(toName)) this%toName = toName
        
    ! Set observation type and determine appropriate npe
    if (present(obsType)) then
        this%obsType = obsType
        call setNpe(this)  ! Update npe based on observation type
    endif
        
    ! Set minimum elevation angle if provided
    if (present(minEle)) this%minEle = minEle
        
    if (present(solarDistance)) this%solarDistance = solarDistance
    if (present(lunarDistance)) this%lunarDistance = lunarDistance
        
    ! Free any previously allocated memory
    if (allocated(this%timeTable)) deallocate(this%timeTable)
    if (allocated(this%dataTable)) deallocate(this%dataTable)
        
    ! Reset data count
    this%nData = 0
        
    end subroutine initAccess
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! read observational data from file
    !
    ! The header may have any number of lines: blank lines and lines starting
    ! with a comment marker ('#', '*', '!', '！') are skipped, so the header
    ! length is not an input any more (see isCommentLine in module_toolKits).
    ! A header written by saveAccessData is read back as well (see
    ! readAccessHeader): it declares the observation type, which fixes npe,
    ! and the names of the transmitter and of the receiver.
    ! Data lines hold: MJD_day, MJD_sec, obs1, ..., obsN (SI units).
    ! A data line whose format does not match - non-numeric text, missing
    ! columns - is reported and skipped, so one damaged line does not spoil
    ! the rest of the file.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine readAccessData(ac, fn)
    implicit none
    class(type_access), intent(inOut) :: ac
    character(len=*), intent(in) :: fn
    double precision, allocatable :: x(:, :)
    integer :: nData, npe, i, nRead
    type(type_t), allocatable :: tTmp(:)
    double precision, allocatable :: dataTmp(:, :)

    !------------------------------------------------------

    ! Adopt the metadata of the file header (observation type, names of the
    ! two devices) before npe is fixed from the observation type
    call readAccessHeader(trim(fn), ac)

    ! Number of data lines: comment and blank lines do not count
    call getDataLinesNum(trim(fn), nData)
    if (nData == 0) then
        print*, "warning @readAccessData"
        print*, "empty file:", trim(fn)
        return
    endif

    ! Get number of parameters per epoch
    npe = ac%npe

    ! Read data and store in matrix; unusable lines are dropped
    allocate(x(nData, 2+npe))
    call readDataMatrix(trim(fn), nData, 2+npe, x, nRead)
    if (nRead == 0) then
        print*, "warning @readAccessData"
        print*, "no usable data line in:", trim(fn)
        deallocate(x)
        return
    endif
    nData = nRead

    ! Free memory if allocated, delete
    if (allocated(ac%timeTable)) deallocate(ac%timeTable)
    if (allocated(ac%dataTable)) deallocate(ac%dataTable)
    ac%nData = 0
    ! Allocate memory
    allocate(ac%timeTable(nData))
    allocate(ac%dataTable(nData, npe))
    ! Initialize data
    do i = 1, nData
        call ac%timeTable(i)%init(x(i, 1), x(i, 2))
    enddo
    ac%dataTable(1:nData, 1:npe) = x(1:nData, 3:npe+2)
    ac%nData = nData
    ! Storage unit -> normalized unit
    call data_SI2normal(ac)

    deallocate(x)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Save the observation data of an access to a text file.
    !
    ! Every header line starts with '#' (see writeAccessHeader) and carries
    ! the identification of the link - generator, transmitter, receiver,
    ! data type - together with the layout of the data lines, so that the
    ! file describes itself. The data lines follow, one epoch per line, in
    ! SI units (MJD day, seconds of day, km, km/s, deg).
    !
    ! 秒列用 13 位、小数点后 6 位（F13.6，微秒分辨率）：RANGING 的观测历元是
    ! 接收设备的钟面时间 t + clk_R，钟差量级为 1e-4 s（0.1 ms）。若秒列只有
    ! 毫秒分辨率，这个偏移会被四舍五入抹掉，读回来的历元就退化成理想历元了。
    ! 字段宽 13 是为了放下 86399.999999。
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine saveAccessData(ac, ts, fn)
    implicit none
    class(type_access), intent(inOut) :: ac
    type(type_timeSeries), intent(in) :: ts
    character(len=*), intent(in) :: fn
    integer i
    double precision alpha, delta
    type(type_t) epoch
    double precision Azimuth, elevation, rho, rhoDot, PV(6)
    !------------------------------------------------------

    if (ac%nData == 0) then
        print*, "warning @saveAccessData. None data observed."
        return
    endif

    open(511, file = trim(fn))
    ! Common header: every line starts with '#', so that readAccessData can
    ! skip the header without knowing how many lines it has
    call writeAccessHeader(511, ac)
    select case(ac%obsType)
    case(OBS_TYPE_RA_DEC)
        write(511, "(A,2x,A5,2x,A13,2x,A10,2x,A10)") "#", "Day", "Sec", "RA/deg", "DEC/deg"
        write(511, "(A,2x,5('-'),2x,13('-'),2x,10('-'),2x,10('-'))") "#"
        do i = 1, ac%nData
            alpha = ac%dataTable(i, 1)*rad2deg
            delta = ac%dataTable(i, 2)*rad2deg
            epoch = ac%timeTable(i)
            write(511, "(2x,I5,2x,F13.6,2x,F10.5,2x,F10.5)") nint(epoch%day), epoch%sec, alpha, delta
        enddo
    case(OBS_TYPE_RANGERATE)
        write(511, "(A,2x,A5,2x,A13,2x,A18)") "#", "Day", "Sec", "Rd/(km/s)"
        write(511, "(A,2x,5('-'),2x,13('-'),2x,18('-'))") "#"
        do i = 1, ac%nData
            rhoDot = ac%dataTable(i, 1)*length_unit/time_unit
            epoch = ac%timeTable(i)
            write(511, "(2x,I5,2x,F13.6,2x,F18.6)") nint(epoch%day), epoch%sec, rhoDot
        enddo
    case(OBS_TYPE_RANGING)
        write(511, "(A,2x,A5,2x,A13,2x,A18)") "#", "Day", "Sec", "R/km"
        write(511, "(A,2x,5('-'),2x,13('-'),2x,18('-'))") "#"
        do i = 1, ac%nData
            epoch = ac%timeTable(i)
            rho = ac%dataTable(i, 1)*length_unit
            write(511, "(2x,I5,2x,F13.6,2x,F18.6)") nint(epoch%day), epoch%sec, rho
        enddo
    case(OBS_TYPE_POSITION_VELOCITY)
        write(511, "(A,2x,A5,2x,A13,6(2x,A18))") "#", "Day", "Sec", &
            &   "X/km", "Y/km", "Z/km", "XDot/(km/s)", "YDot/(km/s)", "ZDot/(km/s)"
        write(511, "(A,2x,5('-'),2x,13('-'),6(2x,'------------------'))") "#"
        do i = 1, ac%nData
            epoch = ac%timeTable(i)
            PV(1:3) = ac%dataTable(i, 1:3)*length_unit
            PV(4:6) = ac%dataTable(i, 4:6)*length_unit/time_unit
            write(511, "(2x,I5,2x,F13.6,3(2x,F18.3),3(2x,F18.6))") nint(epoch%day), epoch%sec, PV(1:6)
        enddo

    case(OBS_TYPE_POSITION)
        write(511, "(A,2x,A5,2x,A13,3(2x,A18))") "#", "Day", "Sec", &
            &   "X/km", "Y/km", "Z/km"
        write(511, "(A,2x,5('-'),2x,13('-'),3(2x,'------------------'))") "#"
        do i = 1, ac%nData
            epoch = ac%timeTable(i)
            PV(1:3) = ac%dataTable(i, 1:3)*length_unit
            write(511, "(2x,I5,2x,F13.6,3(2x,F18.3))") nint(epoch%day), epoch%sec, PV(1:3)
        enddo
    case(OBS_TYPE_AZIMUTH_ELEVATION)
        write(511, "(A,2x,A5,2x,A13,2x,A10,2x,A10)") "#", "Day", "Sec", "Azi/deg", "Ele/deg"
        write(511, "(A,2x,5('-'),2x,13('-'),2x,10('-'),2x,10('-'))") "#"
        do i = 1, ac%nData
            azimuth = ac%dataTable(i, 1)*rad2deg
            elevation = ac%dataTable(i, 2)*rad2deg
            epoch = ac%timeTable(i)
            write(511, "(2x,I5,2x,F13.6,2x,F10.5,2x,F10.5)") nint(epoch%day), epoch%sec, azimuth, elevation
        enddo
    case default
        print*, "warning @subroutine saveData."
        print*, "invalid obsType, obsType = ", ac%obsType
    end select
    
    close(511)
    endsubroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Write the header of an observation data file to unit unitNum.
    !
    ! Every header line starts with '#', so that the reader can skip the
    ! header without knowing how long it is. The header makes the file
    ! self-describing: who generated it, which link produced the data
    ! (transmitter, receiver), which kind of observations they are, and how
    ! the data lines below are laid out.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine writeAccessHeader(unitNum, ac)
    implicit none
    integer, intent(in) :: unitNum
    class(type_access), intent(in) :: ac
    character(len=128) :: typeName, formatStr
    character(len=160) :: fromDev, toDev
    character(len=200) :: line
    character(len=16) :: numStr
    double precision :: mjdFirst, mjdLast
    !------------------------------------------------------
    call obsTypeDescription(ac%obsType, typeName, formatStr)
    call deviceLabel(ac%fromType, ac%fromNum, ac%fromName, fromDev)
    call deviceLabel(ac%toType, ac%toNum, ac%toName, toDev)
    mjdFirst = ac%timeTable(1)%day + ac%timeTable(1)%sec/86400d0
    mjdLast  = ac%timeTable(ac%nData)%day + ac%timeTable(ac%nData)%sec/86400d0

    write(unitNum, "(A)") repeat("#", 78)
    write(unitNum, "(A)") "# Atlas Orbitalis - simulated observation data"
    write(unitNum, "(A)") "# generator   : Atlas (Atlas Orbitalis)"
    write(numStr, "(I6)") ac%number
    line = "# access      : "//trim(adjustl(numStr))//" ("//trim(ac%name)//")"
    write(unitNum, "(A)") trim(line)
    write(unitNum, "(A)") "# transmitter : "//trim(fromDev)
    write(unitNum, "(A)") "# receiver    : "//trim(toDev)
    write(numStr, "(I3)") ac%obsType
    line = "# data type   : "//trim(typeName)//" (obsType = "//trim(adjustl(numStr))//")"
    write(unitNum, "(A)") trim(line)
    write(numStr, "(I9)") ac%nData
    line = "# epochs      : "//trim(adjustl(numStr))
    write(unitNum, "(A)") trim(line)
    write(line, "(A,F15.9,A,F15.9,A)") "# time span   : ", mjdFirst, " ... ", mjdLast, " MJD (UTC)"
    write(unitNum, "(A)") trim(line)
    write(unitNum, "(A)") "# data format : "//trim(formatStr)
    write(unitNum, "(A)") "# note        : header lines start with '#', blank lines are skipped on reading"
    write(unitNum, "(A)") repeat("#", 78)
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Scan the '#' header of an observation data file and adopt the metadata
    ! written by writeAccessHeader: the observation type (which fixes npe,
    ! the number of observable columns per epoch) and the names of the two
    ! devices. A file without such a header, e.g. one produced by another
    ! program, leaves the access object untouched.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine readAccessHeader(fn, ac)
    implicit none
    character(len=*), intent(in) :: fn
    class(type_access), intent(inout) :: ac
    character(len=512) :: line
    integer :: unitNum, io_stat, pos, eq, iType
    logical :: found
    !------------------------------------------------------
    open(newunit = unitNum, file = trim(fn), status = 'old', iostat = io_stat)
    if (io_stat /= 0) return      ! not readable: readAccessData reports it
    do
        read(unitNum, "(A)", iostat = io_stat) line
        if (io_stat /= 0) exit                    ! end of file
        line = adjustl(line)
        if (len_trim(line) == 0) cycle            ! blank line
        if (line(1:1) /= '#') exit                ! data line: the header is over

        ! observation type: "# data type   : RANGING (obsType = 13)"
        pos = index(line, "obsType")
        if (pos > 0) then
            eq = index(line(pos:), "=")
            if (eq > 0) then
                call firstInteger(line(pos+eq:), iType, found)
                if (found) then
                    if (.not. knownObsType(iType)) then
                        print*, "warning @readAccessHeader: unsupported obsType in file header:", iType
                    elseif (.not. headerTypeNameMatches(line, iType)) then
                        ! The 2026 type redesign reused integers 12/14/16 with new
                        ! meanings. A header whose name and integer disagree was
                        ! written by the old table; trusting the integer would
                        ! silently reinterpret its columns, so ignore it.
                        print*, "warning @readAccessHeader: type name in file header does not"
                        print*, "match obsType =", iType, "; header type ignored, access object kept"
                    else
                        if (ac%obsType /= 0 .and. ac%obsType /= iType) then
                            print*, "warning @readAccessHeader: the file header declares obsType =", iType
                            print*, "access object holds obsType =", ac%obsType, ", the file wins"
                        endif
                        ac%obsType = iType
                        call setNpe(ac)
                    endif
                endif
            endif
        endif

        ! device names: "# transmitter : satellite #1 (LEO-600km)"
        if (index(line, "transmitter") > 0) call textInParentheses(line, ac%fromName)
        if (index(line, "receiver") > 0) call textInParentheses(line, ac%toName)
    end do
    close(unitNum)
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Whether the type name in a '#' header agrees with the new observation
    ! type table, e.g. "RANGING" for obsType = 13. Headers without a type
    ! name are accepted (nothing to check). This guards against old files
    ! whose integer code now means something else: the redesign reused
    ! 12/14/16 for AZIMUTH_ELEVATION/RANGERATE/POSITION_VELOCITY.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical function headerTypeNameMatches(line, iType)
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(in) :: iType
    character(len=128) :: nameInFile, typeName, formatStr
    integer :: pos, colon, paren
    !------------------------------------------------------
    headerTypeNameMatches = .true.       ! no type name: nothing to check
    pos = index(line, "data type")
    if (pos == 0) return
    colon = index(line(pos:), ":")
    if (colon == 0) return
    paren = index(line(pos+colon:), "(")
    if (paren == 0) return
    nameInFile = adjustl(line(pos+colon : pos+colon+paren-2))
    if (len_trim(nameInFile) == 0) return   ! empty name: nothing to check
    call obsTypeDescription(iType, typeName, formatStr)
    headerTypeNameMatches = (trim(nameInFile) == trim(typeName))
    end function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Name and column layout of an observation data type,
    ! e.g. "RANGING" and "MJD_day  MJD_sec  Range/km"
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine obsTypeDescription(obsType, typeName, formatStr)
    implicit none
    integer, intent(in) :: obsType
    character(len=*), intent(out) :: typeName
    character(len=*), intent(out) :: formatStr
    !------------------------------------------------------
    select case(obsType)
    case(OBS_TYPE_RA_DEC)
        typeName  = "RA_DEC"
        formatStr = "MJD_day  MJD_sec  RA/deg  DEC/deg"
    case(OBS_TYPE_AZIMUTH_ELEVATION)
        typeName  = "AZIMUTH_ELEVATION"
        formatStr = "MJD_day  MJD_sec  Az/deg  El/deg"
    case(OBS_TYPE_RANGING)
        typeName  = "RANGING"
        formatStr = "MJD_day  MJD_sec  Range/km"
    case(OBS_TYPE_RANGERATE)
        typeName  = "RANGERATE"
        formatStr = "MJD_day  MJD_sec  RangeRate/(km/s)"
    case(OBS_TYPE_POSITION)
        typeName  = "POSITION"
        formatStr = "MJD_day  MJD_sec  X/km  Y/km  Z/km"
    case(OBS_TYPE_POSITION_VELOCITY)
        typeName  = "POSITION_VELOCITY"
        formatStr = "MJD_day  MJD_sec  X/km  Y/km  Z/km  Vx/(km/s)  Vy/(km/s)  Vz/(km/s)"
    case default
        typeName  = "UNKNOWN"
        formatStr = "unknown"
    end select
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Whether npe (and hence the column layout) is defined for an
    ! observation data type
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical function knownObsType(obsType)
    implicit none
    integer, intent(in) :: obsType
    !------------------------------------------------------
    select case(obsType)
    case(OBS_TYPE_RA_DEC, OBS_TYPE_AZIMUTH_ELEVATION, OBS_TYPE_RANGING, &
        & OBS_TYPE_RANGERATE, OBS_TYPE_POSITION, OBS_TYPE_POSITION_VELOCITY)
        knownObsType = .true.
    case default
        knownObsType = .false.
    end select
    end function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Text label of one end of a link, e.g. "satellite #1 (LEO-600km)",
    ! used in the data file header and in showAccess
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine deviceLabel(objType, objNum, objName, label)
    implicit none
    integer, intent(in) :: objType, objNum
    character(len=*), intent(in) :: objName
    character(len=*), intent(out) :: label
    character(len=16) :: typeStr, numStr
    !------------------------------------------------------
    if (objType == OBJECT_TYPE_SATELLITE) then
        typeStr = "satellite"
    elseif (objType == OBJECT_TYPE_STATION) then
        typeStr = "station"
    else
        typeStr = "unknown"
    endif
    write(numStr, "(I3)") objNum
    label = trim(typeStr)//" #"//trim(adjustl(numStr))
    if (len_trim(objName) > 0) label = trim(label)//" ("//trim(objName)//")"
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! First integer of a piece of text, e.g. 13 in "= 13)", where the digits
    ! may be followed by text that a list-directed read would reject
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine firstInteger(text, value, found)
    implicit none
    character(len=*), intent(in) :: text
    integer, intent(out) :: value
    logical, intent(out) :: found
    integer :: i, i1, i2, io_stat
    !------------------------------------------------------
    value = 0
    found = .false.
    i1 = 0
    i2 = 0
    do i = 1, len(text)
        if (text(i:i) >= '0' .and. text(i:i) <= '9') then
            if (i1 == 0) i1 = i
            i2 = i
        elseif (i1 > 0) then
            exit                              ! the digits end here
        endif
    end do
    if (i1 == 0) return
    read(text(i1:i2), *, iostat = io_stat) value
    found = (io_stat == 0)
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Text between the first pair of parentheses of a line, e.g. "LEO-600km"
    ! in "# transmitter : satellite #1 (LEO-600km)"
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine textInParentheses(line, text)
    implicit none
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: text
    integer :: p1, p2
    !------------------------------------------------------
    text = ''
    p1 = index(line, "(")
    if (p1 == 0) return
    p2 = index(line(p1+1:), ")")
    if (p2 == 0) return
    text = line(p1+1 : p1+p2-1)
    end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine setNpe(this)
    implicit none
    class(type_access) :: this
    select case(this%obsType)
    case(OBS_TYPE_RA_DEC)
        this%npe = 2
    case(OBS_TYPE_AZIMUTH_ELEVATION)
        this%npe = 2
    case(OBS_TYPE_RANGING)
        this%npe = 1
    case(OBS_TYPE_RANGERATE)
        this%npe = 1
    case(OBS_TYPE_POSITION)
        this%npe = 3
    case(OBS_TYPE_POSITION_VELOCITY)
        this%npe = 6
    case default
        print*, "warning @suborutine setNpe."
        print*, "invalid obsType, obsType = ", this%obsType
    end select
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine data_SI2normal(this)
    implicit none
    class(type_access) :: this
    integer :: nData
    nData = this%nData
    if (nData == 0) then
        print*, "warning @data_SI2normal"
        write(*, *) "nData = 0 !"
    endif
    select case(this%obsType)
    case(OBS_TYPE_RA_DEC)
        this%dataTable(1:nData, 1:2) = this%dataTable(1:nData, 1:2)*deg2rad
    case(OBS_TYPE_AZIMUTH_ELEVATION)
        this%dataTable(1:nData, 1:2) = this%dataTable(1:nData, 1:2)*deg2rad
    case(OBS_TYPE_RANGING)
        this%dataTable(1:nData, 1) = this%dataTable(1:nData, 1)/length_unit
    case(OBS_TYPE_RANGERATE)
        this%dataTable(1:nData, 1) = this%dataTable(1:nData, 1)/length_unit*time_unit
    case(OBS_TYPE_POSITION)
        this%dataTable(1:nData, 1:3) = this%dataTable(1:nData, 1:3)/length_unit
    case(OBS_TYPE_POSITION_VELOCITY)
        this%dataTable(1:nData, 1:3) = this%dataTable(1:nData, 1:3)/length_unit
        this%dataTable(1:nData, 4:6) = this%dataTable(1:nData, 4:6)/length_unit*time_unit
    endSelect
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine data_normal2SI(this)
    implicit none
    class(type_access) :: this
    integer nData
    nData = this%nData
    select case(this%obsType)
    case(OBS_TYPE_RA_DEC)
        this%dataTable(1:nData, 1:2) = this%dataTable(1:nData, 1:2)/deg2rad
    case(OBS_TYPE_AZIMUTH_ELEVATION)
        this%dataTable(1:nData, 1:2) = this%dataTable(1:nData, 1:2)/deg2rad
    case(OBS_TYPE_RANGING)
        this%dataTable(1:nData, 1) = this%dataTable(1:nData, 1)*length_unit
    case(OBS_TYPE_RANGERATE)
        this%dataTable(1:nData, 1) = this%dataTable(1:nData, 1)*length_unit/time_unit
    case(OBS_TYPE_POSITION)
        this%dataTable(1:nData, 1:3) = this%dataTable(1:nData, 1:3)*length_unit
    case(OBS_TYPE_POSITION_VELOCITY)
        this%dataTable(1:nData, 1:3) = this%dataTable(1:nData, 1:3)*length_unit
        this%dataTable(1:nData, 4:6) = this%dataTable(1:nData, 4:6)*length_unit/time_unit
    endSelect
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! cal obs data for sat - station access
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calObsData1(ts, sat, station, access, force)
    use module_perturbationSRP
    use module_planets
    implicit none
    type(type_timeSeries), intent(in) :: ts
    type(type_satellite), intent(inOut) :: sat
    type(type_station), intent(inOut) :: station
    type(type_access), intent(inOut) :: access
    logical, intent(in) :: force
    integer :: nt, npe
    type(type_t), allocatable :: timeTable(:)
    double precision, allocatable :: dataTable(:, :)
    integer :: ndata, it
    double precision :: utc, pvt_J2000(6), pvr_J2000(6), pvt_fixed(6), pvr_fixed(6)
    double precision :: drHorizon(3), A, h, dr(3), alpha, delta, rho, rhoDot, h1
    type(type_t) :: t, TT
    double precision :: pvSun(6), pvMoon(6), v_station_sat(3), v_station_sun(3), v_station_moon(3)
    double precision :: clkErrT, clkErrR       ! 发射端(卫星)/接收端(测站)钟差 [s]
    logical :: flag_shadow = .false.
    !------------------------------------------------------
    if (allocated(access%timeTable)) deallocate(access%timeTable)
    if (allocated(access%dataTable)) deallocate(access%dataTable)
    
    ! ������ʱ�������Ȱѹ۲����ݷ�����ʱ������
    nt = ts%n ! ʱ����
    npe = access%npe! npe��ÿ��ʱ�̵Ĺ۲�������
    allocate(timeTable(nt))
    allocate(dataTable(nt, npe))
    
    dataTable = 0d0
    nData = 0
    ! ����ʱ��
    do it = 1, nt
        
        ! epoch, position/velocity of transmitter/reciever
        utc = ts%time(it)%val
        t = ts%time(it)
        pvt_J2000 = sat%getPV(t) ! pv of the transmitter
        pvr_fixed = (/ station%position, 0d0, 0d0, 0d0 /) ! pv of the reciever in ECF
        pvr_J2000 = station%getPV(t) ! pv of the reciever in ECI
        
        ! light-time delay
        if (key_light_time_delay) then
            call light_time_delay(pvt_J2000, pvr_J2000)
        endif
        
        if (force .or. access%obsType == OBS_TYPE_POSITION_VELOCITY .or. access%obsType == OBS_TYPE_POSITION) then
            if (access%obsType == OBS_TYPE_AZIMUTH_ELEVATION) then 
                call ECI2ECF_IAU1980_positionVelocity(utc, pvt_J2000(1:3), &
                    &   pvt_J2000(4:6), pvt_fixed(1:3), pvt_fixed(4:6)) ! ���ǣ��ع�����

                call ECF2horizon_position(station%longitude, &
                    &   station%latitude, pvt_fixed(1:3) - pvr_fixed(1:3), drHorizon)
                h = dasin(drHorizon(3)/norm2(drHorizon))                ! get h from sinh, -pi/2 <= h <= pi/2
                call get_angle(drHorizon(2)/norm2(drHorizon(1:2)), &
                    &   drHorizon(1)/norm2(drHorizon(1:2)), A)                  ! get A from sinA and cosA, 0 <= A < 2*pi
            endif
        else
            ! calculate azimuth/elevation
            call ECI2ECF_IAU1980_positionVelocity(utc, pvt_J2000(1:3), &
                &   pvt_J2000(4:6), pvt_fixed(1:3), pvt_fixed(4:6)) ! ���ǣ��ع�����

            call ECF2horizon_position(station%longitude, &
                &   station%latitude, pvt_fixed(1:3) - pvr_fixed(1:3), drHorizon)
            h = dasin(drHorizon(3)/norm2(drHorizon))                ! get h from sinh, -pi/2 <= h <= pi/2
            call get_angle(drHorizon(2)/norm2(drHorizon(1:2)), &
                &   drHorizon(1)/norm2(drHorizon(1:2)), A)                  ! get A from sinA and cosA, 0 <= A < 2*pi
                        
            ! height is too low
            if (h < access%minEle*deg2rad) then
                cycle
            endif

            ! the sun interrupt observation
            if (access%solarDistance > 0d0) then
                TT = t
                call TT%toTT
                call planetsPV(TT%val, nEarth, nSun, pvSun)
                pvSun(1:3) = pvSun(1:3)/length_unit
                v_station_sat = pvt_J2000(1:3) - pvr_J2000(1:3)
                v_station_sun = pvSun(1:3) - pvr_J2000(1:3)
                ! if angle between v_station_sat and v_station_sun is less than SolarDistance
                ! then the access is invalid
                if (dot_product(v_station_sat, v_station_sun)/norm2(v_station_sat)/norm2(v_station_sun) &
                & > dcos(access%solarDistance*deg2rad)) then
                    cycle
                endif
            endif
            
            ! the moon interrupt observation
            if (access%lunarDistance > 0d0) then
                TT = t
                call TT%toTT
                call planetsPV(TT%val, nEarth, nMoon, pvMoon)
                pvMoon(1:3) = pvMoon(1:3)/length_unit
                v_station_sat = pvt_J2000(1:3) - pvr_J2000(1:3)
                v_station_moon = pvMoon(1:3) - pvr_J2000(1:3)
                ! if angle between v_station_sat and v_station_moon is less than LunarDistance
                ! then the access is invalid
                if (dot_product(v_station_sat, v_station_moon)/norm2(v_station_sat)/norm2(v_station_moon) &
                & > dcos(access%lunarDistance*deg2rad)) then
                    cycle
                endif
            endif
            
            ! the satellite is in the shadow of Earth
            if (access%shadow) then
                TT = t
                call TT%toTT
                call planetsPV(TT%val, nEarth, nSun, pvSun)
                pvSun(1:3) = pvSun(1:3)/length_unit
                call in_shadow(pvt_J2000(1:3), pvSun(1:3), flag_shadow)
                if (flag_shadow) then
                    cycle
                endif
            endif
            
        endif
        
        ! consider atmospherical refraction
        if (key_atmospherical_refraction) then
            call tropospheric_refraction(h, 293.15d0, 0.7d0, h1)
            h = h1
        endif
        

        ! if access is constructed successfully
        nData = nData + 1
        timeTable(nData) = t
        select case(access%obsType)
        
        case(OBS_TYPE_RA_DEC)
            dr = pvt_J2000(1:3) - pvr_J2000(1:3)
            delta = dasin(dr(3)/norm2(dr)) ! get delta from sin(delta), -pi/2 <= delta <= pi/2
            call get_angle(dr(2)/norm2(dr(1:2)), &
                &   dr(1)/norm2(dr(1:2)), Alpha) ! get alpha from sinAlpha and cosAlpha, 0 <= Alpha < 2*pi
            dataTable(nData, 1:2) = (/alpha, delta/)

        case(OBS_TYPE_RANGERATE)
            rho = norm2(pvt_J2000(1:3) - pvr_J2000(1:3))
            rhoDot = dot_product(pvt_J2000(4:6) - pvr_J2000(4:6), pvt_J2000(1:3) - pvr_J2000(1:3))/rho
            dataTable(nData, 1) = rhoDot

        case(OBS_TYPE_RANGING)
            rho = norm2(pvt_J2000(1:3) - pvr_J2000(1:3))
            dataTable(nData, 1) = rho
            
            ! 这里，仿照sat-sat链路的做法，考虑收发设备的钟差修正
            ! 观测时刻timeTable(nData) = t + clock_error(接收端)
            ! dataTable(nData, 1) = rho + (clock_error(接收端) - clock_error(发射端))*c
            if (key_clock_error) then
                ! 钟差取自钟模型，单位 s，约定 bias = 钟面时间 - 参考时间
                ! (钟差为正表示钟走快)
                clkErrT = sat%clock%get_clock_state(t)      ! 发射端(卫星)钟差
                clkErrR = station%clock%get_clock_state(t)  ! 接收端(测站)钟差

                ! 收发钟差之差表现为距离偏差。SPEED_OF_LIGHT 为 m/s，u_meter 把米
                ! 换算成归一化长度，故乘得"归一化长度/秒"，与归一化的 rho 同量纲
                dataTable(nData, 1) = rho + (clkErrR - clkErrT)*SPEED_OF_LIGHT*u_meter

                ! 接收端是计时的一方，观测时刻记为其钟面读数
                ! (timeTable 在 select case 之前已按理想时刻赋值，这里改写)
                timeTable(nData) = timeAdd(t, clkErrR/86400d0)
            endif

        case(OBS_TYPE_POSITION_VELOCITY)
            dataTable(nData, 1:6) = pvt_J2000(1:6)

        case(OBS_TYPE_POSITION)
            dataTable(nData, 1:3) = pvt_J2000(1:3)

        case(OBS_TYPE_AZIMUTH_ELEVATION)
            dataTable(nData, 1:2) = (/A, h/)
            
        case default
            print*, "error @subroutine calObsData1."
            print*, "invalid obsType, obsType = ", access%obsType
        end select

        
                
    enddo
    
    ! transfer data from temporary variable to access
    access%nData = nData
    if (nData >= 1) then
        allocate(access%timeTable(nData))
        allocate(access%dataTable(nData, npe))
        access%timeTable(1:nData) = timeTable(1:nData)
        access%dataTable(1:ndata, 1:npe) = dataTable(1:nData, 1:npe)
    endif
    
    ! free RAM
    deallocate(timeTable)
    deallocate(dataTable)
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate observational data for sat-sat access
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calObsData2(ts, sat1, sat2, access, force)
    use module_planets
    implicit none
    type(type_timeSeries), intent(in) :: ts
    type(type_satellite), intent(inOut) :: sat1
    type(type_satellite), intent(inOut) :: sat2
    type(type_access), intent(inOut) :: access
    logical, intent(in) :: force
    integer :: nt, npe
    type(type_t), allocatable :: timeTable(:)
    double precision, allocatable :: dataTable(:, :)
    integer :: ndata, it
    double precision :: utc, pvt_J2000(6), pvr_J2000(6), pvt_fixed(6), pvr_fixed(6)
    double precision :: drHorizon(3), A, h, dr(3), alpha, delta, rho, rhoDot
    double precision :: r1r2c(3), d1
    double precision :: clkErrT, clkErrR          ! 发射星/接收星钟差 [s]
    type(type_t) :: t, TT
    double precision :: pvMoon(6), rt_moon(3), rr_moon(3), dm
    double precision :: transmitter_point(3), t2r(3), cosine
    !------------------------------------------------------
    
    if (access%obsType == OBS_TYPE_AZIMUTH_ELEVATION) then
        print*, "warning @calObsData2"
        print*, "the Ah-type data is invalid in sat-sat access!"
        return
    endif
    
    ! 两个数组都要先释放: access 可能是上一次调用留下的(定轨每迭代一次就重算
    ! 一遍观测量), 漏掉 timeTable 会导致第二次调用"重复分配"直接崩掉
    if (allocated(access%timeTable)) deallocate(access%timeTable)
    if (allocated(access%dataTable)) deallocate(access%dataTable)

    ! allocate RAM
    nt = ts%n !
    npe = access%npe
    allocate(timeTable(nt))
    allocate(dataTable(nt, npe))

    dataTable = 0d0
    nData = 0
    ! iteration for all epochs
    do it = 1, nt
        ! epoch
        utc = ts%time(it)%val
        t = ts%time(it)
        pvt_J2000 = sat1%getPV(t) ! pv for transmitter
        pvr_J2000 = sat2%getPV(t) ! pv for reciever

        ! Earth blockage
        if (access%EarthBlockage) then
            call cross_3d(pvt_J2000(1:3), pvr_J2000(1:3), r1r2c)
            d1 = norm2(r1r2c)/norm2(pvt_J2000(1:3) - pvr_J2000(1:3)) ! d1�ǵ��ĵ����������ߵľ���
            if (d1 < 1d0 + 80d0/length_unit) then ! ������ڵ�
                cycle
            endif
        endif
        
        ! Moon blockage
        if (access%MoonBlockage) then
            ! get pv of the Moon
            TT = t
            call TT%toTT
            call planetsPV(TT%val, nEarth, nMoon, pvMoon)
            pvMoon(1:3) = pvMoon(1:3)/length_unit
            rt_moon = pvt_J2000(1:3) - pvMoon(1:3)
            rr_moon = pvr_J2000(1:3) - pvMoon(1:3)
            call cross_3d(rt_moon, rr_moon, r1r2c)
            dm = NORM2(r1r2c)/NORM2(rr_moon - rt_moon)
            if (dm < radius_moon/length_unit) then
                cycle
            endif
        endif
        
        ! transmitter coverage constraint
        if (access%t_coverage) then ! if there is a coverage constraint for transmitter
            ! calculate the angle between satellite position vector and transmitter-to-reciever vector
            if (access%pointTo == 1) then
                ! point to the center of the Earth
                transmitter_point = -pvt_J2000(1:3)/norm2(pvt_J2000(1:3))
            elseif (access%pointTo == 2) then
                ! point to the reciever
                transmitter_point = (pvr_J2000(1:3) - pvt_J2000(1:3))/norm2(pvr_J2000(1:3) - pvt_J2000(1:3))
            else
                print*, "warning @subroutine calObsData2."
                print*, "invalid pointTo, pointTo = ", access%pointTo
                stop
            endif
            t2r = (pvr_J2000(1:3) - pvt_J2000(1:3))/norm2(pvr_J2000(1:3) - pvt_J2000(1:3))
            cosine = dot_product(transmitter_point, t2r)
            ! the off-boresight angle of the receiver must not exceed the beam half angle
            if (cosine < dcos(access%beamHalfAngle*deg2rad)) then
                cycle
            endif
        endif

        ! light-time delay
        if (key_light_time_delay) then
            call light_time_delay(pvt_J2000, pvr_J2000)
        endif
        

        ! record data at this epoch
        nData = nData + 1
        select case(access%obsType)
            
        case(OBS_TYPE_RA_DEC)
            dr = pvt_J2000(1:3) - pvr_J2000(1:3)
            delta = dasin(dr(3)/norm2(dr)) ! get delta from sin(delta), -pi/2 <= delta <= pi/2
            call get_angle(dr(2)/norm2(dr(1:2)), &
                &   dr(1)/norm2(dr(1:2)), Alpha) ! get alpha from sinAlpha and cosAlpha, 0 <= Alpha < 2*pi
            dataTable(nData, 1:2) = (/alpha, delta/)
            timeTable(nData) = t
        case(OBS_TYPE_RANGING)
            rho = norm2(pvt_J2000(1:3) - pvr_J2000(1:3))
            dataTable(nData, 1) = rho
            timeTable(nData) = t

            ! 考虑钟差导致的观测数据和观测时刻偏差
            ! 观测时刻timeTable(nData) = t + clock_error(接收星)
            ! dataTable(nData, 1) = rho + (clock_error(接受星) - clock_error(发射星))*c
            if (key_clock_error) then
                ! 钟差取自卫星钟模型，单位 s，约定 bias = 卫星钟面时间 - 参考时间
                ! (钟差为正表示卫星钟走快)
                clkErrT = sat1%clock%get_clock_state(t)   ! 发射星钟差
                clkErrR = sat2%clock%get_clock_state(t)   ! 接收星钟差

                ! 收发钟差之差表现为距离偏差。SPEED_OF_LIGHT 为 m/s，u_meter 把米
                ! 换算成归一化长度，故乘得"归一化长度/秒"，与归一化的 rho 同量纲
                dataTable(nData, 1) = rho + (clkErrR - clkErrT)*SPEED_OF_LIGHT*u_meter

                ! 接收星是计时的一方，观测时刻记为其钟面读数
                timeTable(nData) = timeAdd(t, clkErrR/86400d0)
            endif

        case(OBS_TYPE_POSITION_VELOCITY)
            dataTable(nData, 1:6) = pvt_J2000(1:6)
            timeTable(nData) = t
        case(OBS_TYPE_POSITION)
            dataTable(nData, 1:3) = pvt_J2000(1:3)
            timeTable(nData) = t
        case(OBS_TYPE_RANGERATE)
            rho = norm2(pvt_J2000(1:3) - pvr_J2000(1:3))
            rhoDot = dot_product(pvt_J2000(4:6) - pvr_J2000(4:6), pvt_J2000(1:3) - pvr_J2000(1:3))/rho
            dataTable(nData, 1) = rhoDot
            timeTable(nData) = t
        case default
            print*, "error @subroutine calObsData2."
            print*, "invalid obsType, obsType = ", access%obsType
        end select
                
    enddo
    
    access%nData = nData
    if (nData >= 1) then
        allocate(access%timeTable(nData))
        allocate(access%dataTable(nData, npe))
        access%timeTable(1:nData) = timeTable(1:nData)
        access%dataTable(1:ndata, 1:npe) = dataTable(1:nData, 1:npe)
    endif
    
    deallocate(dataTable)
    deallocate(timeTable)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calObsData3(ts, sat, access, force)
    implicit none
    type(type_timeSeries), intent(in) :: ts
    type(type_satellite), intent(inOut) :: sat
    type(type_access), intent(inOut) :: access
    logical, intent(in) :: force
    integer :: nt, npe
    type(type_t), allocatable :: timeTable(:)
    double precision, allocatable :: dataTable(:, :)
    integer :: ndata, it
    double precision :: utc, pvt_J2000(6), pvr_J2000(6), pvt_fixed(6), pvr_fixed(6)
    double precision :: drHorizon(3), A, h, dr(3), alpha, delta, rho, rhoDot
    type(type_t) :: t
    !------------------------------------------------------
    if (allocated(access%timeTable)) deallocate(access%timeTable)
    if (allocated(access%dataTable)) deallocate(access%dataTable)
    
    ! ������ʱ�������Ȱѹ۲����ݷ�����ʱ������
    nt = ts%n ! ʱ����
    npe = access%npe! npe��ÿ��ʱ�̵Ĺ۲�������
    allocate(timeTable(nt))
    allocate(dataTable(nt, npe))
    
    dataTable = 0d0
    nData = 0
    ! ����ʱ��
    do it = 1, nt
        ! ��ô˿̵ķ����豸λ�úͽ����豸λ��
        utc = ts%time(it)%val
        t = ts%time(it)
        pvt_J2000 = sat%getPV(t) ! ���ǣ�J2000����

        ! ����۲�����ֵ
        nData = nData + 1
        select case(access%obsType)

        case(OBS_TYPE_POSITION_VELOCITY)
            dataTable(nData, 1:6) = pvt_J2000(1:6)
            timeTable(nData) = t
        case(OBS_TYPE_POSITION)
            dataTable(nData, 1:3) = pvt_J2000(1:3)
            timeTable(nData) = t
        case default
            print*, "error @subroutine calObsData3."
            print*, "invalid obsType, obsType = ", access%obsType
        end select
                
    enddo
    
    ! �ѹ۲����ݴ���ʱ������ת�Ƶ�access������
    access%nData = nData
    allocate(access%dataTable(nData, npe))
    allocate(access%timeTable(nData))
    access%dataTable(1:ndata, 1:npe) = dataTable(1:nData, 1:npe)
    access%timeTable(1:nData) = timeTable
    ! �ͷ���ʱ�������ڴ�
    deallocate(timeTable)
    deallocate(dataTable)
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine showAccess(this)
    implicit none
    class(type_access) :: this
    character(len=160) :: fromDev, toDev
    character(len=64)  :: typeName, formatStr
    character(len=16)  :: numStr
    print*, "--------------------------------------------------------------"
    print*, "access number: ", this%number
    print"(A,A)", " name: ", this%name
    call deviceLabel(this%fromType, this%fromNum, this%fromName, fromDev)
    call deviceLabel(this%toType, this%toNum, this%toName, toDev)
    print"(A,A,A,A)", " from ", trim(fromDev), " to ", trim(toDev)
    call obsTypeDescription(this%obsType, typeName, formatStr)
    write(numStr, "(I3)") this%obsType
    print"(A,A,A,A,A)", " data type: ", trim(typeName), " (obsType = ", trim(adjustl(numStr)), ")"
    print*, "--------------------------------------------------------------"
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to calculate light time delay correction
    ! Takes satellite position/velocity and observer position vectors
    ! Iteratively corrects for light travel time
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine light_time_delay(rs, ro)
    use module_const
    implicit none
    
    ! Input/Output variables
    real(8), intent(inout) :: rs(6)  ! Satellite state vector (position and velocity)
    real(8), intent(in)    :: ro(6)  ! Observer position vector
    
    ! Local variables
    integer :: i
    real(8) :: dt, dist
    real(8) :: rs_corrected(3)
    
    ! Iteratively correct for light time delay
    ! Typically 3 iterations are sufficient for convergence
    rs_corrected(1:3) = rs(1:3)
    do i = 1, 3
                
        ! Calculate distance between satellite and observer
        dist = norm2(rs_corrected(1:3) - ro(1:3))
        
        ! Calculate light travel time
        dt = dist/SPEED_OF_LIGHT_UNIT
        
        ! Apply correction to satellite position based on its velocity
        rs_corrected(1:3) = rs(1:3) - rs(4:6)*dt
        
    enddo
    
    ! Update satellite position
    rs(1:3) = rs_corrected(1:3)
    
    end subroutine light_time_delay
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the observational data at one epoch
    ! pvt_J2000: satellite position and velocity in J2000 ECI frame
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine observation_oneEpoch(t, longi, lati, pvt, pvr, obsType, npe, data1)
    implicit none
    double precision, intent(in) :: t, longi, lati, pvt(6), pvr(6)
    integer, intent(in) :: obsType, npe
    double precision, intent(out) :: data1(npe)
    double precision :: pvt_J2000(6), pvr_J2000(6)
    double precision :: dr(3), alpha, delta, rho, rhoDot
    double precision :: drHorizon(3), A, h
    double precision :: pvt_fixed(6), pvr_fixed(6)
    
    
    pvt_J2000 = pvt
    pvr_J2000 = pvr
    if (key_light_time_delay) then
        call light_time_delay(pvt_J2000, pvr_J2000)
    endif
    ! Calculate relative position/velocity
    dr = pvt_J2000(1:3) - pvr_J2000(1:3)
    
    select case(obsType)
    case(OBS_TYPE_RA_DEC)
        ! Right ascension and declination
        delta = dasin(dr(3)/norm2(dr)) ! declination
        call get_angle(dr(2)/norm2(dr(1:2)), dr(1)/norm2(dr(1:2)), alpha) ! right ascension
        data1(1:2) = (/alpha, delta/)

    case(OBS_TYPE_RANGERATE)
        ! Range-rate only
        rho = norm2(dr)
        rhoDot = dot_product(pvt_J2000(4:6) - pvr_J2000(4:6), dr)/rho
        data1(1) = rhoDot

    case(OBS_TYPE_RANGING)
        ! Range only
        data1(1) = norm2(dr)

    case(OBS_TYPE_POSITION_VELOCITY)
        ! Position and velocity
        data1(1:6) = pvt_J2000(1:6)

    case(OBS_TYPE_POSITION)
        ! Position only
        data1(1:3) = pvt_J2000(1:3)

    case(OBS_TYPE_AZIMUTH_ELEVATION)
        ! Azimuth and elevation
        call ECI2ECF_IAU1980_positionVelocity(t, pvt_J2000(1:3), pvt_J2000(4:6), pvt_fixed(1:3), pvt_fixed(4:6))
        call ECI2ECF_IAU1980_positionVelocity(t, pvr_J2000(1:3), pvr_J2000(4:6), pvr_fixed(1:3), pvr_fixed(4:6))
        call ECF2horizon_position(longi, lati, pvt_fixed(1:3) - pvr_fixed(1:3), drHorizon)
        h = dasin(drHorizon(3)/norm2(drHorizon))
        call get_angle(drHorizon(2)/norm2(drHorizon(1:2)), drHorizon(1)/norm2(drHorizon(1:2)), A)
        data1(1:2) = (/A, h/)
        
    case default
        print*, "warning @observation_oneEpoch: invalid obsType = ", obsType
    end select

    endsubroutine observation_oneEpoch

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine AzimuthElevation2RaDec(station, access)
    ! Purpose: Convert azimuth/elevation observations to right ascension/declination
    ! Input: 
    !   station - ground station object (type_station)
    !   access  - access object with azimuth/elevation data (type_access)
    ! Output:
    !   access%obsType changed to OBS_TYPE_RA_DEC
    !   access%dataTable converted from (Azi, Ele) to (RA, DEC)
    !
    ! Conversion steps:
    ! 1. (Azi, Ele) -> cartesian in horizon frame
    ! 2. horizon -> ITRS (Earth-fixed)
    ! 3. ITRS -> ICRS (J2000 inertial)
    ! 4. cartesian -> spherical (RA, DEC)
    
    use module_space
    use module_const
    implicit none
    type(type_station), intent(in) :: station
    type(type_access), intent(inout) :: access
    
    integer :: i, nData
    double precision :: azi, ele, utc
    double precision :: r_horizon(3), r_ITRS(3), r_ICRS(3)
    double precision :: v_ITRS(3), v_ICRS(3)  ! dummy velocity vectors (zero)
    double precision :: ra, dec, r_norm
    
    ! Check input validity
    if (access%obsType /= OBS_TYPE_AZIMUTH_ELEVATION) then
        print*, "Warning @AzimuthElevation2RaDec: "
        print*, "  Expected obsType = OBS_TYPE_AZIMUTH_ELEVATION (", OBS_TYPE_AZIMUTH_ELEVATION, ")"
        print*, "  But got obsType = ", access%obsType
        return
    endif
    
    if (access%npe /= 2) then
        print*, "Warning @AzimuthElevation2RaDec: "
        print*, "  Expected npe = 2 for azimuth/elevation data"
        print*, "  But got npe = ", access%npe
        return
    endif
    
    nData = access%nData
    if (nData <= 0) then
        print*, "Warning @AzimuthElevation2RaDec: "
        print*, "  No data to convert (nData = ", nData, ")"
        return
    endif
    
    ! Initialize velocity vectors to zero (not used in position-only transformation)
    v_ITRS = 0d0
    v_ICRS = 0d0
    
    ! Convert each observation point
    do i = 1, nData
        ! Get azimuth and elevation (already in radians from data_SI2normal)
        azi = access%dataTable(i, 1)
        ele = access%dataTable(i, 2)
        
        ! Step 1: Convert (Azi, Ele) to cartesian in horizon frame
        ! Horizon frame: X points North, Y points East, Z points Up
        ! Standard spherical to cartesian conversion:
        ! x = cos(ele) * cos(azi)
        ! y = cos(ele) * sin(azi) 
        ! z = sin(ele)
        r_horizon(1) = dcos(ele) * dcos(azi)
        r_horizon(2) = dcos(ele) * dsin(azi)
        r_horizon(3) = dsin(ele)
        
        ! Step 2: Convert horizon frame to ITRS (Earth-fixed)
        ! Note: horizon2ECF_position expects position relative to station
        ! but r_horizon is a unit vector in local horizon frame
        ! We need to convert this unit vector to ECF coordinates
        call horizon2ECF_position(station%longitude, station%latitude, &
                                 r_horizon, r_ITRS)
        
        ! Step 3: Convert ITRS to ICRS (J2000 inertial)
        ! Get UTC time for this observation
        utc = access%timeTable(i)%val
        
        ! Convert position from ITRS to ICRS
        call ECF2ECI_IAU1980_positionVelocity(utc, r_ITRS, v_ITRS, &
                                             r_ICRS, v_ICRS)
        
        ! Step 4: Convert cartesian to spherical (RA, DEC)
        ! r_ICRS is a unit vector in ICRS frame
        r_norm = norm2(r_ICRS)
        if (r_norm > 0d0) then
            ! Declination: asin(z/r), range [-pi/2, pi/2]
            dec = dasin(r_ICRS(3)/r_norm)
            
            ! Right ascension: atan2(y, x), range [0, 2*pi)
            ra = datan2(r_ICRS(2), r_ICRS(1))
            if (ra < 0d0) ra = ra + 2d0*pi
        else
            ! Should not happen for unit vector, but handle gracefully
            ra = 0d0
            dec = 0d0
        endif
        
        ! Store converted values
        access%dataTable(i, 1) = ra
        access%dataTable(i, 2) = dec
    enddo
    
    ! Update observation type
    access%obsType = OBS_TYPE_RA_DEC
    
    ! Note: npe remains 2 (RA and DEC), so no need to change
    
    endsubroutine

    endmodule
