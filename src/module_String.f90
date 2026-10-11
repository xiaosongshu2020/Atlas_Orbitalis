module String
    use ErrorHandler
    implicit none

    contains
    
        !function int2str(int) result(str)
        !    !
        !    ! 2017-10-31
        !    !
        !    implicit none
        !    ! input
        !    integer,intent(in):: int
        !    ! output
        !    character(:),allocatable:: str
        !    ! local
        !    integer:: n
        !    if (int.eq.0) then
        !        allocate(character(len=1)::str)
        !        write(str,'("0")')
        !    elseif (int.gt.0) then
        !        n = idint(dlog10(dble(int)))+1
        !        allocate(character(len=n)::str)
        !        write(str,'(i)') int
        !    else    ! negative
        !        n = idint(dlog10(dble(-int)))+1
        !        allocate(character(len=n+1)::str)
        !        write(str,'("-",i)') -int
        !    endif
        !    return
        !end function int2str
    
        function int2str(int) result(str)
            implicit none
            integer, intent(in) :: int
            character(:), allocatable :: str
            character(32) :: buffer  ! 临时缓冲区

            if (int == 0) then
                str = '0'
            else
                ! 使用I0格式自动处理宽度
                write(buffer, '(I0)') int
                str = trim(buffer)
            end if
        end function int2str

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 2023-11-27 Shawn Hu 
        ! inputs:
        !   int         integer             input integer
        !   fmt         character(len=*)    output string's format
        !   n           integer             output string's length
        ! outputs:
        !   str         character           output string
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function int2str2(int, fmt, n) result(str)
        implicit none
        integer, intent(in) :: int, n
        character(len=*), intent(in) :: fmt
        character(len = n) :: str
        !----------------------------------------------------
        write(str, fmt) int
        endfunction
    !
        function double2str(x,prec) result(str)
            !
            ! 2018-12-29
            !
            implicit none
            ! input
            double precision,intent(in):: x
            integer,intent(in),optional:: prec
            ! output
            character(:),allocatable:: str
            character(64) :: buffer
            ! local
            integer:: n,p
            if (present(prec)) then
                p = prec
            else
                p = 6
            endif
            if (dabs(x).gt.1d-15) then
                n = idint(dlog10(x))+1+p
            else
                n = 1+p
            endif
            allocate(character(len=n+1)::str)
            ! 使用缓冲区写入，然后修剪空格
            write(buffer,'(f'//trim(int2str(n+1))//'.'//trim(int2str(p))//')') x
            str = trim(adjustl(buffer))
            return
        end function double2str
    !
        integer function strFindInArray(str,strArr,priority)
            !
            ! Find index of the designated string in an array
            ! 2017-10-31
            !
            implicit none
            ! input
            character(len=*),intent(in):: str,strArr(:)
            character(len=*),intent(in),optional:: priority
            ! local
            integer:: nStr,i,idx
            logical:: retFirst
            character(len=4):: subName='fsia'
            ! save
            logical,save:: warnDup = .true.
            ! code
            call checkin(subName)
            nStr = size(strArr)
            idx = 0
            retFirst = .not.present(priority)
            if (present(priority)) retFirst = strcmpit(priority,'first')
            do i = 1,nStr
                if (strcmpit(str,strArr(i))) idx = i
                if ((idx.gt.0).and.retFirst) exit
            enddo
            strFindInArray = idx
            call checkout(subName)
            return
        end function strFindInArray

    !
        logical function strcmpit(str1,str2)
            implicit none
            ! input
            character(len=*),intent(in):: str1,str2
            ! local
            character(len=4):: subName='scit'
            call checkin(subName)
            strcmpit = strcmpi(trim(str1),trim(str2))
            if (errCode.gt.0) return
            call checkout(subName)
            return
        end function strcmpit
    !
        logical function strcmpi(str1,str2)
            implicit none
            ! input
            character(len=*),intent(in):: str1,str2
            ! local
            integer:: i,strlen,U1,L1,U2,L2
            character(len=4):: subName='scpi'
            character(len=26):: upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
            character(len=26):: lower='abcdefghijklmnopqrstuvwxyz'
            ! code
            call checkin(subName)
            strcmpi = .false.
            if (len(str1).ne.len(str2)) goto 99
            strlen = len(str1)
            do i = 1, strlen
                if (str1(i:i).ne.str2(i:i)) then
                    U1 = index(upper,str1(i:i))
                    L1 = index(lower,str1(i:i))
                    U2 = index(upper,str2(i:i))
                    L2 = index(lower,str2(i:i))
                    if ((U1*U1+L1*L1.eq.0).or.(U2*U2+L2*L2.eq.0)) goto 99
                    if (.not.((U1.eq.L2).or.(L1.eq.U2))) goto 99
                endif
            enddo
            strcmpi = .true.
 99         call checkout(subName)
            return
        end function strcmpi
    !
        subroutine UpperCase(instr,outstr)
            implicit none
            integer i,j,k
            character(len=*):: instr,outstr
            character(len=26):: lower,upper
            character(len=4):: subname='loca'
            call checkin(subname)
            data lower/'abcdefghijklmnopqrstuvwxyz'/
            data upper/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

            i=len(instr)
            outstr = ''
            outstr(1:i) = instr(1:i)
            do j = 1,i
                k = index(lower,instr(j:j))
                if (k.gt.0) outstr(j:j) = upper(k:k)
            end do
            call checkout(subname)
        end subroutine UpperCase
    !
        subroutine LowerCase(instr,outstr)
        !   subroutine: converts a string to lower case.
            use ErrorHandler
            implicit none
            integer i,j,k
            character(len=*):: instr,outstr
            character(len=26):: lower,upper
            character(len=4):: subname='loca'
            call checkin(subname)
            data lower/'abcdefghijklmnopqrstuvwxyz'/
            data upper/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

            i=len(instr)
            outstr = ''
            outstr(1:i) = instr(1:i)
            do j = 1,i
                k = index(upper,instr(j:j))
                if (k.gt.0) outstr(j:j) = lower(k:k)
            end do
            call checkout(subname)
        end subroutine LowerCase

        function Lower(instr)
            ! function: converts a string to lower case
            use ErrorHandler
            implicit none
            integer j,k
            character(len=*),intent(in):: instr
            character(len=:),allocatable:: Lower
            character(len=26):: loweralphabet,upperalphabet
            character(len=4):: subname='f2lo'

            call checkin(subname)
            data loweralphabet/'abcdefghijklmnopqrstuvwxyz'/
            data upperalphabet/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

            Lower = instr
            do j = 1,len(instr)
                k = index(upperalphabet,instr(j:j))
                if (k.gt.0) Lower(j:j) = loweralphabet(k:k)
            end do
            call checkout(subname)
            return
        end function Lower

        function Upper(instr)
            ! function: converts a string to upper case
            use ErrorHandler
            implicit none
            integer j,k
            character(len=*),intent(in):: instr
            character(len=:),allocatable:: Upper
            character(len=26):: loweralphabet,upperalphabet
            character(len=4):: subname='f2up'

            call checkin(subname)
            data loweralphabet/'abcdefghijklmnopqrstuvwxyz'/
            data upperalphabet/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

            Upper = instr
            do j = 1,len(instr)
                k = index(loweralphabet,instr(j:j))
                if (k.gt.0) Upper(j:j) = upperalphabet(k:k)
            end do
            call checkout(subname)
            return
        end function Upper

        function trim2(str)
            ! function: trim both leading and trailing blanks of a string
            use ErrorHandler
            implicit none
            character(len=*),intent(in):: str
            character(len=:),allocatable:: trim2
            trim2 = trim(adjustl(str))
            return
        end function trim2

        integer function countWord(inStr)
            implicit none
            character(len=*),intent(in):: inStr
            character(len=1000):: tmpStr
            integer:: iSpace,nWords
            tmpStr = adjustl(inStr)
            nWords = 0
            do while (len_trim(tmpStr).ne.0)
                iSpace = index(tmpStr,' ')
                tmpStr = trim(adjustl(tmpStr(iSpace:len_trim(tmpStr))))
                nWords = nWords + 1
            enddo
            countWord = nWords
            return
        end function countWord

        subroutine getWord(inStr,nWords,strList)
            implicit none
            ! input
            character(len=*),intent(in):: inStr
            ! output
            integer,intent(out),optional:: nWords
            character(len=*),allocatable:: strList(:)
            ! local
            integer:: n
            character(len=4):: subName='gtwd'
            ! code
            call checkin(subName)
            n = countWord(inStr)
            if (allocated(strList)) deallocate(strList)
            allocate(strList(n))
            read(inStr,*,err=99) strList(1:n)
            if (present(nWords)) nWords = n
            call checkout(subName)
            return
 99         call RaiseError(0008,subName)
            return
        end subroutine getWord

        function strRep(inStr,oldStr,newStr) result(outStr)
            implicit none
            ! input
            character(len=*),intent(in):: inStr,oldStr,newStr
            ! output
            character(len=:),allocatable:: outStr
            ! local
            integer:: subIdx
            logical:: isFound
            character(len=4):: subName='srep'
            ! code
            call checkin(subName)
            outStr = trim(inStr)
            if (trim(inStr).eq.'') then
                call RaiseWarning('Input string is empty.',subName)
                return
            endif
            subIdx = index(outStr,oldStr)
            if (subIdx.eq.0) return

            do while (subIdx.gt.0)
                outStr = outStr(1:subIdx-1)//trim(newStr)//trim(outStr(subIdx+len(oldStr):len_trim(outStr)))
                subIdx = index(outStr,oldStr)
            enddo

            call checkout(subName)
            return
        end function strRep

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 判断字符串str2是否是str1的一部分，不区分大小写
        ! 例如str1 = '0 CHINASAT 32 (BEIDOU 1)', str2 = 'BEIDOU', 返回.true.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function stringContains(str1, str2) result(con)
        character(len=*), intent(in) :: str1, str2
        logical :: con
        integer :: i, j, len1, len2
        character(len=:),allocatable :: string1, string2
        
        con = .false.
        
        string1 = lower(str1)
        string2 = lower(str2)
        
        len1 = len_trim(string1)
        len2 = len_trim(string2)
        
        if (len1 > len2) then
            return
        end if
    
        do i = 1, len2 - len1 + 1
            if (string1(1:len1) == string2(i:i+len1-1)) then
                con = .true.
                return
            end if
        end do
        end function stringContains
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end
