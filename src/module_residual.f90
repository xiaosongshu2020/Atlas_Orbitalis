
    module module_residual
    use module_access
    use module_const
    
    type :: type_residual
        integer :: nAcc = 0
        type(type_acc), allocatable :: acc(:)
        double precision :: outlier = 1d2
        double precision :: RMS = 1d2
    contains
        procedure :: calValidIdx
        procedure :: allocRes
        procedure :: freeRes
        procedure :: calRMS
    endtype
    
    type :: type_acc
        integer :: accIdx
        integer :: obsType
        integer :: npe
        integer :: nData
        type(type_t), allocatable :: timeTable(:)
        double precision, allocatable :: res(:, :)
        integer :: nValid = 0
        integer, allocatable :: validIdx(:)
    contains
        procedure :: calRes
        procedure :: saveRes
    endtype
    
    type(type_residual) :: res
    
    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine allocRes(this, n)
    implicit none
    class(type_residual) :: this
    integer, intent(in) :: n
    if (allocated(this%acc)) deallocate(this%acc)
    allocate(this%acc(n))
    this%nAcc = n
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calValidIdx(this, aci)
    implicit none
    class(type_residual) :: this
    integer, intent(in) :: aci
    integer i, nData
    nData = size(this%acc(aci)%timeTable)
    if (allocated(this%acc(aci)%validIdx)) deallocate(this%acc(aci)%validIdx)
    allocate(this%acc(aci)%validIdx(nData))
    this%acc(aci)%nValid = 0
    do i = 1, nData
        if (norm2(this%acc(aci)%res(i, :)) < this%outlier) then
            this%acc(aci)%nValid = this%acc(aci)%nValid + 1
            this%acc(aci)%validIdx(this%acc(aci)%nValid) = i
        endif
    enddo
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine freeRes(this)
    class(type_residual) :: this
    this%nAcc = 0
    deallocate(this%acc)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the resdiuals, which are the difference between accO and accC
    ! residuals are saved into this%res
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calRes(this, accO, accC)
    class(type_acc) :: this
    type(type_access), intent(in) :: accO, accC
    !--------------------------------------------------------
    this%accIdx = accO%number
    this%obsType = accO%obsType
    this%npe = accO%npe
    this%nData = accO%nData
    this%timeTable = accO%timeTable
    this%res = accO%dataTable - accC%dataTable
    ! Wrap the angular residual that is periodic. RA/Dec and Az/El both put
    ! the first component in [0, 2*pi); without this, an observation crossing
    ! the 0/2*pi seam produces a 2*pi outlier and wrecks the least-squares.
    select case(this%obsType)
    case(OBS_TYPE_RA_DEC, OBS_TYPE_AZIMUTH_ELEVATION)
        this%res(:, 1) = modulo(this%res(:, 1) + pi, 2d0*pi) - pi
    end select
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write residual data into files
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine saveRes(this, fn)
    implicit none
    class(type_acc) :: this
    character(len=*), intent(in), optional :: fn
    character(len=200) :: fn1
    character(len=64) :: typeName, formatStr
    character(len=80) :: unitStr
    character(len=80) :: fmt
    double precision :: scale(10)
    integer :: i, j, unit, ios, nOut
    integer, save :: NO = 1
    !--------------------------------------------------------
    if (.not. allocated(this%res)) return
    if (present(fn)) then
        fn1 = trim(fn)
    else
        write(fn1, "('data/output/res_',I3.3,'.txt')") NO
        NO = NO + 1
    endif

    ! Residuals are stored normalized; write them in display units.
    scale = 1d0
    unitStr = "normalized"
    select case(this%obsType)
    case(OBS_TYPE_RA_DEC, OBS_TYPE_AZIMUTH_ELEVATION)
        scale(1:2) = RAD2ARCSEC
        unitStr = "arcsec  arcsec"
    case(OBS_TYPE_RANGING)
        scale(1) = 1d0/u_meter
        unitStr = "m"
    case(OBS_TYPE_RANGERATE)
        scale(1) = 1d0/u_meter/TIME_UNIT
        unitStr = "m/s"
    case(OBS_TYPE_POSITION)
        scale(1:3) = 1d0/u_meter
        unitStr = "m  m  m"
    case(OBS_TYPE_POSITION_VELOCITY)
        scale(1:3) = 1d0/u_meter
        scale(4:6) = 1d0/u_meter/TIME_UNIT
        unitStr = "m  m  m  m/s  m/s  m/s"
    end select

    call obsTypeDescription(this%obsType, typeName, formatStr)

    open(newunit = unit, file = trim(fn1), status = 'replace', iostat = ios)
    if (ios /= 0) then
        print*, "warning @saveRes: cannot open ", trim(fn1)
        return
    endif

    write(unit, "(A)") repeat("#", 78)
    write(unit, "(A)") "# Atlas Orbitalis - OD residuals (observed - computed)"
    write(unit, "(A,I0,A,A,A)") "# observation type : ", this%obsType, " (", trim(typeName), ")"
    write(unit, "(A,I0)") "# components       : ", this%npe
    write(unit, "(A,I0)") "# epochs           : ", this%nData
    write(unit, "(A,A)") "# units            : ", trim(unitStr)
    write(unit, "(A)") "# data format      : MJD_day  MJD_sec  res1 ... resN"
    write(unit, "(A)") repeat("#", 78)

    if (this%npe >= 1 .and. this%npe <= 10) then
        write(fmt, "('(I6,2x,F13.6,',I0,'(2x,ES16.8))')") this%npe
        nOut = min(this%nData, size(this%res, 1))
        if (allocated(this%timeTable)) nOut = min(nOut, size(this%timeTable))
        do i = 1, nOut
            write(unit, fmt) nint(this%timeTable(i)%day), this%timeTable(i)%sec, &
                &   (this%res(i, j)*scale(j), j = 1, this%npe)
        enddo
    endif
    close(unit)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calRMS(this)
    implicit none
    class(type_residual) :: this
    double precision sum
    integer nData, i, j, idx
    sum = 0d0
    nData = 0
    do i = 1, this%nAcc
        do j = 1, this%acc(i)%nValid
            idx = this%acc(i)%validIdx(j)
            sum = sum + (norm2(this%acc(i)%res(idx, :)))**2
            nData = nData + 1
        enddo
    enddo
    this%RMS = dsqrt(sum / nData)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endmodule