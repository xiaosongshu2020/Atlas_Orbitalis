module FileDef
    !----------------------------------------------------------------------------------------------------------------------------------
    !   Instruction:
    !       This module defines the file struct and the routines to open and close an instance.
    !   Date:
    !       2013-09-26
    !   Contains:
    !       OpenFile    sub     Opens a file with given file type
    !       CloseFile   sub     Closes a file with given file type
    !   Note:
    !       When subroutines in this module encounters an error, we do not promptly return, but jump to the end and
    !       check out before quitting. This is because file manipulations are some fundamental operations, it would be
    !       meaningless if you check a callstack and always find OpenFile/CloseFile at the end.
    !----------------------------------------------------------------------------------------------------------------------------------
    implicit none

    integer,parameter:: maxFilePerTab=100
    integer,parameter:: DataFileUnitBase=100, TempFileUnitBase=300, UnitCeiling=9999

    type File
        integer:: unit = -1
        character(len=200):: name = ''
        logical:: newfile=.false.
        logical:: temporary=.false.
        contains
            procedure:: open => OpenFile
    end type File
    type FileTable
        integer:: nFile
        type(File):: Files(maxFilePerTab)
    end type FileTable

    type(FileTable),protected:: OpenedFile

    contains
        integer function RequestNewFileunit(base)
        !----------------------------------------------------------------------------------------------------------------------------------
        !   Instruction:
        !       This function returns a unit for a new file, if possible.
        !   Date:
        !       2014-05-11
        !   Inputs:
        !       base    integer     unit base from which a new available unit is searched
        !   Outputs:
        !       The unit for the new file is returned.
        !----------------------------------------------------------------------------------------------------------------------------------
            use ErrorHandler
            implicit none
            ! Declare input parameters:
            integer,intent(in):: base
            ! Declare local variables:
            integer:: i,reqUnit
            character(len=4):: subName='rqfu'

            call checkin(subName)
            if (OpenedFile%nFile.ge.maxFilePerTab) then
                call RaiseError(4018,subName)
                return
            endif

            reqUnit = base
 99         continue
            do i = 1, OpenedFile%nFile
                if (OpenedFile%Files(i)%unit.eq.reqUnit) then
                    if (reqUnit.ge.UnitCeiling) then
                        call RaiseError(4019,subName)
                        return
                    endif
                    reqUnit = reqUnit+1
                    goto 99
                endif
            enddo

!            OpenedFile%nFile = OpenedFile%nFile+1
            RequestNewFileUnit = reqUnit
            call checkout(subName)
            return

        end function RequestNewFileunit

        type(File) function OpenFile2(Fname,isNew)
            ! A transparent wrapper function of OpenFile.
            ! Do not (need to) check in/out or detect error.
            implicit none
            logical,intent(in),optional:: isNew
            character(len=*), intent(in):: Fname
            OpenFile2%name = trim(Fname)
            if (present(isNew)) OpenFile2%newfile = isNew
            call OpenFile(OpenFile2)
            return
        end function OpenFile2

        subroutine OpenFile(F)
        !----------------------------------------------------------------------------------------------------------------------------------
        !   Instruction:
        !       This subroutine opens a file with given file struct.
        !   Date:
        !       2013-09-26
        !   Inputs:
        !       F       type(file)      the file struct to open
        !   Outputs:
        !       No explicit output. If it fails to open the file, an error is raised.
        !----------------------------------------------------------------------------------------------------------------------------------
            use ErrorHandler
            implicit none
            ! Declare input parameters:
            class(File):: F
            ! Declare the local variables:
            logical:: isExist,isOpen
            character(len=4):: subName='opfi'
            character(len=10):: fstatus

            call checkin(subName)
            inquire(file=trim(F%name),exist=isExist,opened=isOpen)
            if (isOpen) then
                call RaiseWarning('Opening an opened file is attempted and will be ignored.',subName)
                goto 98
            endif
            if ((.not.isExist).and.(.not.F%newfile)) then
                call RaiseError(4005,subName)
                goto 98
            endif

            if (F%newfile) then
                fstatus = 'unknown'
            else
                fstatus = 'old'
            endif

            if (F%temporary) then
                F%unit = RequestNewFileunit(TempFileUnitBase)
            else
                F%unit = RequestNewFileunit(DataFileUnitBase)
            endif
            if (errCode.gt.0) return
            OpenedFile%nFile = OpenedFile%nFile+1
            OpenedFile%Files(OpenedFile%nFile) = F

            open(F%unit,file=trim(F%name),status=trim(fstatus),err=99)
            goto 98
 99         call RaiseError(4003,subName)

 98         call checkout(subName)
            return
        end subroutine OpenFile

        subroutine CloseFile(FUnit)
        !----------------------------------------------------------------------------------------------------------------------------------
        !   Instruction:
        !       This subroutine closes a file with given file struct.
        !   Date:
        !       2013-09-26
        !   Inputs:
        !       FUnit      integer      the unit of the file to be closed
        !   Outputs:
        !       No explicit output. If it fails to close the file, an error is raised.
        !----------------------------------------------------------------------------------------------------------------------------------
            use ErrorHandler
            implicit none
            ! Declare input parameters:
            integer,intent(in):: FUnit
            ! Declare local variables:
            integer:: i,j
            logical:: isOpen
            type(File):: F
            character(len=4):: subName='clfi'
            character(len=10):: fstatus

            call checkin(subName)
            inquire(FUnit,opened=isOpen)
            if (.not.isOpen) then
                call RaiseWarning('Closing an unopened file is attempted and will be ignored.',subName)
                goto 98
            endif

            F%unit = 0
            do i = 1,OpenedFile%nFile
                if (OpenedFile%Files(i)%unit.eq.FUnit) then
                    F = OpenedFile%Files(i)
                    do j = i,OpenedFile%nFile-1
                        OpenedFile%Files(j) = OpenedFile%Files(j+1)
                    enddo
                    OpenedFile%nFile = OpenedFile%nFile-1
                    exit
                endif
            enddo
            if (F%unit.eq.0) then
                call RaiseError(4020,subName)
                return
            endif

            if (F%temporary) then
                fstatus = 'delete'
            else
                fstatus = 'keep'
            endif

            close(F%unit,status=trim(fstatus),err=99)
            goto 98
 99         call RaiseError(4004,subName)

 98         call checkout(subName)
            return
        end subroutine CloseFile

end module FileDef
