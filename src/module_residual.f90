
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
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write residual data into files
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine saveRes(this, fn)
    implicit none
    class(type_acc) :: this
    character(len=*), intent(in), optional :: fn
    character*100 :: fn1
    integer, save :: NO = 1
    !--------------------------------------------------------
    if (present(fn)) then
        open(821, file = trim(fn))
    else
        write(fn1, "('data/output/res_',I3.3,'.txt')") NO
        open(821, file = trim(fn1))
        NO = NO + 1
    endif
    
    close(821)
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