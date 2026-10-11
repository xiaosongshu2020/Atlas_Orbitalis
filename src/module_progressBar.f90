Module CMD_ProgressBar
    Implicit None
  
    Type , public :: CLS_CMD_Progress
        Integer , private :: N
        Integer , private :: lens = 40
        Integer , private :: i
        Character :: M = "#" , O = "_"
        Character(len=64) :: Prefix = ''
    Contains
        Procedure :: Set
        Procedure :: Put
    End Type CLS_CMD_Progress
  
    Logical , parameter , public :: CMD_PROGRESS_ABSOLUTE = .true.
    type(CLS_CMD_Progress), public :: Progress
  
    contains
    
    ! setting, it is optional
    Subroutine Set( this , N , L , Prefix)
        Class( CLS_CMD_Progress ) :: this
        Integer , Intent( IN ) :: N , L
        Character(len=*), Intent( IN ) :: Prefix
        this % N    = N
        this % lens = L
        this % i = 0
        this % Prefix = trim(Prefix)
    End Subroutine Set
  
    subroutine Put( this , N , K , bAbsol , Suffix )
        Class( CLS_CMD_Progress ) :: this
        Integer , Intent( IN ) :: N, K
        Logical , optional :: bAbsol
        ! Suffix: 条尾附带的文字(如"  5392/17383  剩余 1.4 s")。可选，不给则与旧行为一致
        Character(len=*), optional , Intent( IN ) :: Suffix
        Character(len=1) :: br
        integer :: jm
        this % N = N
        this % i = this % i + K
        if ( present( bAbsol ) ) then
            if ( bAbsol ) this % i = K
        else
            this % i = K
        end if
        if ( this % i > this % n ) this % i = this % n
        jm = Nint( real( this%i * this%lens ) / real( this%N ) )
        if ( this%i < this%n ) then
          br = char(13)
        else
          br = char(10)
        end if
        write( * , '(5a,f6.2,a)',advance="no") trim(this%Prefix) , '[' , &
        ! write( * , '(5a,f6.2,2a\)') trim(this%Prefix) , '[' , &
        repeat(this%M , jm ) , repeat( this%O , this%lens-jm ) , '] ' , this%i*100.0/this%N , "%"
        if ( present( Suffix ) ) then
            write( * , '(a)', advance="no") trim(Suffix)
        end if
        write( * , '(a)', advance="no") br
    End subroutine Put
  
    End Module CMD_ProgressBar

    ! subroutine example_progressBar
    ! use CMD_progressBar
    ! integer i, n
    ! n = 1000
    ! progress%prefix = "word:"
    ! do i = 1, n
    !     call progress%put(n, i, .true.)
    ! enddo
    ! ! 可选的第 5 个参数是在条尾追加文字(如已处理个数、累计用时):
    ! ! call progress%put(n, i, .true., "  1200/1000  累计 3.4 s")
    ! endsubroutine