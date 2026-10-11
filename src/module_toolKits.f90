    
    module module_toolKits

    implicit none
    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! count lines number of a file
    subroutine get_lines_num(filename, n_lines)
    implicit none
    character(len=*), intent(in)    ::  fileName
    character*100 c_tmp
    integer is_end, n_lines
    logical :: file_exists
    !-----------------------------------------------------------
    inquire(file=trim(filename), exist=file_exists)
    if (.not. file_exists) then
        write(*, *) 'Warning @subroutine get_lines_num'
        write(*,*) 'File "', trim(filename), '" does not exist.'
        n_lines = 0
        return
    end if

    open(372, file = trim(filename))
    is_end = 0
    n_lines = 0 !number of lines in this file
    do while (is_end == 0)
        read(372, *, iostat = is_end) c_tmp
        n_lines = n_lines + 1
    end do
    n_lines = n_lines - 1
    close(372)
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine fix(a, b)
    double precision a, b
    if (a == 0) then
        b = 0
    elseif (a > 0) then
        b = dble(floor(a))
    else
        b = dble(floor(a)) + 1
    end if
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This subroutine is written to get the value of the angle while 
    ! we know its sine and cosine.
    !The returned value of function 'dacos' is between 0 and pi, 
    ! however the value of the angle
    ! is between 0 and 2*pi
    subroutine get_angle(sin_phi, cos_phi, phi)
    implicit none
    double precision sin_phi, cos_phi, phi
    double precision sin_phi2, cos_phi2
    !----------------------------------------------
    ! Normalize sin_phi and cos_phi
    double precision norm
    norm = sqrt(sin_phi**2 + cos_phi**2)
    if (norm > 0d0) then
        sin_phi2 = sin_phi / norm
        cos_phi2 = cos_phi / norm
    end if
    phi = dacos(cos_phi2)
    if (sin_phi2 < 0) then
        phi = 2d0*3.141592653589793d0 - phi
    end if
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine demo_LagrangianInterpolation
    Implicit none
    Integer :: i, fileid 
    Integer, parameter :: n1 = 5!// ��ֵ���ڵ���
    Real(kind=8) :: x(n1), y(n1)
    integer, parameter :: n2 = 10
    real(kind=8) :: xx(n2), yy(n2)
  
    !Open ( newunit = fileid, file = 'origin.dat' )
    !Do i = 1, n1
    !  x(i) = t0 + ( i-1 ) * ( t1-t0 ) / ( n1-1 )
    !  y(i) = exp( x(i) )
    !  Write ( fileid,'(2es)' ) x(i), y(i)
    !End do 
    !Close ( fileid )
    x = [1.0, 2.0, 3.0, 4.0, 5.0]
    y = [10.0, 5.0, 8.0, 2.0, 6.0]
    xx = (/1.0d0, 1.5d0, 2.0d0, 2.5d0, 3d0, 3.5d0, 4d0, 4.5d0, 5d0, 5.5d0/)
    call Interpolation(n1, x, y, n2, xx, yy)
    write(*, *) xx
    write(*, *) yy
  
    Endsubroutine 
    
    subroutine demo2_LagrangianInterpolation
    Implicit none
    Integer :: i, fileid 
    Integer, parameter :: n1 = 1000!// ��ֵ���ڵ���
    Real(kind=8) :: x(n1), y(n1)
    integer, parameter :: n2 = 1000
    real(kind=8) :: xx(n2), yy(n2)
    double precision, parameter :: pi = 3.141592653589793d0
    
    !// ��ȡ��ֵ���ڵ�����
    do i = 1, n1
        x(i) = dble(i)
    enddo
    do i = 1, n1
        y(i) = 10d0 + 2d0*x(i) + 5d0*dcos(x(i)*2d0*pi/100d0)
    enddo
    
    do i = 1, n2
        xx(i) = dble(i) + 0.5d0
    enddo
    call Interpolation(n1, x, y, n2, xx, yy)
    
    open(261, file = 'data/data.txt')
    do i = 1, n1
        write(261, *) x(i), y(i), xx(i), yy(i)
    enddo
    close(261)
    
    Endsubroutine 
    
    Subroutine Interpolation( n1, x, y, n2, xx, yy)
    Implicit none
    Integer :: i, j, k, fileid
    Integer, intent(in) :: n1
    Integer, intent(in) :: n2
    Real(kind=8), intent(in) :: x(n1), y(n1)
    Real(kind=8), intent(in) :: xx(n2)
    real(kind=8), intent(out) :: yy(n2)
    real(kind=8)  L1(n1), L2(n1), tmp
    
    Do k = 1, n2
        Do i = 1, n1
        L1(i) = 1.d0
        Do j = 1, n1
            If ( j /= i ) then
            L1(i) = ( x(i) - x(j) ) * L1(i)  !// ���ĸ
            End if 
        End do
        
        L2(i) = 1.d0
        Do j = 1, n1
            If ( j /= i ) then
            L2(i) = ( xx(k) - x(j) ) * L2(i)  !// �����
            End if
        End do
        End do
        
        tmp = 0.d0
        Do i = 1, n1
        tmp = y(i) * L2(i) / L1(i) + tmp  !// ���ֵ��
        End do
        yy(k) = tmp
    End do
    End subroutine Interpolation
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine cross_3d ( x, y, z )
    implicit none
    double precision x(3), y(3), z(3)
    !----------------------------------------
    z(1) = y(3)*x(2) - y(2)*x(3)
    z(2) = y(1)*x(3) - y(3)*x(1)
    z(3) = y(2)*x(1) - y(1)*x(2)
    end subroutine cross_3d
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dot(n, x, y, outcome)
    implicit none
    integer n, ii
    double precision x(n), y(n), outcome
    !-----------------------------------
    outcome = 0d0
    do ii = 1, n
        outcome = outcome + x(ii)*y(ii)
    end do
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Whether a line of a data file is a comment (or blank) line, in which
    ! case the data readers skip it. The first non-blank character decides:
    !   '#'   '*'   '!'     or any non-ASCII character (the full-width '！'
    ! is a multi-byte character, so its encoding - UTF-8, GBK - is not assumed)
    ! Data lines always start with a sign, a digit or a decimal point, so this
    ! test never rejects a data line.
    logical function isCommentLine(line)
    implicit none
    character(len=*), intent(in) :: line
    integer :: i
    !-----------------------------------------------------------
    ! blank (or empty) line -> skip
    isCommentLine = .true.
    do i = 1, len(line)
        if (line(i:i) /= ' ' .and. line(i:i) /= achar(9)) exit
    end do
    if (i > len(line)) return
    ! comment markers
    if (line(i:i) == '#' .or. line(i:i) == '*' .or. line(i:i) == '!') return
    if (line(i:i) > achar(127)) return      ! full-width '！' and other non-ASCII
    ! it is a data line
    isCommentLine = .false.
    end function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Whether a line carries data: it is neither blank nor a comment line (see
    ! isCommentLine) and its first field parses as a number.
    !
    ! The numeric test makes the readers work with the older files as well,
    ! whose header lines carry no marker at all; such lines are then reported
    ! by the readers instead of being taken for data.
    logical function isDataLine(line)
    implicit none
    character(len=*), intent(in) :: line
    double precision :: first_field
    integer :: io_stat
    !-----------------------------------------------------------
    isDataLine = .false.
    if (isCommentLine(line)) return
    read(line, *, iostat = io_stat) first_field
    isDataLine = (io_stat == 0)
    end function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Number of blank-separated fields of a line. Used to recognise the layout
    ! of a data file (7 or 8 columns) without assuming a format.
    integer function countFields(line)
    implicit none
    character(len=*), intent(in) :: line
    integer :: i
    logical :: in_field
    !-----------------------------------------------------------
    countFields = 0
    in_field = .false.
    do i = 1, len(line)
        if (line(i:i) /= ' ' .and. line(i:i) /= achar(9)) then
            if (.not. in_field) then
                countFields = countFields + 1
                in_field = .true.
            end if
        else
            in_field = .false.
        end if
    end do
    end function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! First data line of a file (see isDataLine), returned in line with found
    ! set to .false. when the file holds no data line at all.
    subroutine getFirstDataLine(filename, line, found)
    implicit none
    character(len=*), intent(in)    ::  fileName
    character(len=*), intent(out)   ::  line
    logical,          intent(out)   ::  found
    character(len=512) tmp
    integer is_end, file_unit
    !-----------------------------------------------------------
    line  = ''
    found = .false.
    open(newunit = file_unit, file = trim(filename), status = 'old', iostat = is_end)
    if (is_end /= 0) return
    do
        read(file_unit, '(A)', iostat = is_end) tmp
        if (is_end /= 0) exit
        if (.not. isDataLine(tmp)) cycle
        line  = tmp
        found = .true.
        exit
    end do
    close(file_unit)
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Number of data lines of a text file (see isDataLine). The file header may
    ! therefore have any number of lines: it does not have to be counted.
    subroutine getDataLinesNum(filename, n_data)
    implicit none
    character(len=*), intent(in)    ::  fileName
    integer,          intent(out)   ::  n_data
    character(len=512) line
    integer is_end
    logical :: file_exists
    !-----------------------------------------------------------
    n_data = 0
    inquire(file=trim(filename), exist=file_exists)
    if (.not. file_exists) then
        write(*, *) 'Warning @subroutine getDataLinesNum'
        write(*,*) 'File "', trim(filename), '" does not exist.'
        return
    end if

    open(373, file = trim(filename))
    do
        read(373, '(A)', iostat = is_end) line
        if (is_end /= 0) exit
        if (.not. isDataLine(line)) cycle
        n_data = n_data + 1
    end do
    close(373)
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read the first n data lines of a file into mat(n, m), skipping blank
    ! lines, comment lines (see isCommentLine) and any other line that does not
    ! start with a number (see isDataLine).
    !
    ! A line that starts with a number but cannot be read as m numbers -
    ! non-numeric text, too few columns, ... - is reported and skipped, so that
    ! a single damaged line does not spoil the rest of the file; the lines that
    ! follow it move up.
    !
    ! The number of lines actually stored is returned in nRead.
    subroutine readDataMatrix(fn, n, m, mat, nRead)
    implicit none
    character(len=*), intent(in)    :: fn
    integer,          intent(in)    :: n, m
    double precision, intent(out)   :: mat(n, m)
    integer,          intent(out), optional :: nRead
    character(len=512) line
    integer file_unit, io_stat, nData, nBad, nSkip
    !-----------------------------------------------------------
    mat    = 0d0
    nData  = 0
    nBad   = 0
    nSkip  = 0

    open(newunit = file_unit, file = trim(fn), status = 'old', iostat = io_stat)
    if (io_stat /= 0) then
        write(*, *) 'Error @subroutine readDataMatrix: cannot open file'
        write(*, *) trim(fn)
        if (present(nRead)) nRead = 0
        return
    end if

    do
        read(file_unit, '(A)', iostat = io_stat) line
        if (io_stat /= 0) exit                  ! end of file

        if (isCommentLine(line)) cycle          ! blank line or marked comment

        if (.not. isDataLine(line)) then
            ! text that is neither marked as a comment nor a set of numbers,
            ! e.g. an old style header line: skip it, but not silently
            nSkip = nSkip + 1
            if (nSkip <= 3) then
                write(*, *) 'Warning @subroutine readDataMatrix: line skipped, not data', &
                            ' and not marked with #, *, ! or ！:'
                write(*, *) '   ', trim(line)
            end if
            cycle
        end if

        if (nData >= n) exit                    ! keep only the first n data lines

        ! Read tentatively into the next row and keep it only when the line
        ! holds all m numbers; a malformed line is dropped, not zero-filled
        read(line, *, iostat = io_stat) mat(nData + 1, :)
        if (io_stat /= 0) then
            nBad = nBad + 1
            mat(nData + 1, :) = 0d0
            write(*, *) 'Warning @subroutine readDataMatrix: malformed data line skipped', &
                        '(does not hold the expected', m, 'numbers):'
            write(*, *) trim(line)
            cycle
        end if
        nData = nData + 1
    end do
    close(file_unit)

    if (nSkip > 0) then
        write(*, *) 'Warning @subroutine readDataMatrix: ', nSkip, ' line(s) skipped in'
        write(*, *) trim(fn)
        write(*, *) '  (header lines should start with #, *, ! or ！)'
    end if
    if (nBad > 0) then
        write(*, *) 'Warning @subroutine readDataMatrix: ', nBad, ' malformed data line(s) skipped in'
        write(*, *) trim(fn)
    end if

    if (present(nRead)) nRead = nData
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endmodule