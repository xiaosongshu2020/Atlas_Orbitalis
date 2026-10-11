module ErrorHandler
    !----------------------------------------------------------------------------------------------------------------------------------
    !   Instruction:
    !       This module generates error message
    !   Date:
    !       2013-08-28
    !   Update:
    !       2014-06-21
    !   Public variables:
    !       errCode         i               code of the latest error
    !       errMsg          str             error messge of the given code
    !       errHeap         str             heap of calling subs' function ID
    !   Contains:
    !       RaiseError      sub             process newly raised error
    !       RaiseWarning    sub             write warning message to stdout
    !       UpdateErrMsg    sub             update error message with given code
    !----------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer,parameter:: errHeapSize = 100
    integer,parameter:: errMsgLen  = 200
    integer,parameter:: errFileLen = 200
    integer,parameter:: maxErrCode = 999999
    integer,parameter:: maxUserErr = 999
    integer,protected:: errCode = 0
    integer,protected:: nWarning = 0
    integer,protected:: errHeapTop = 0
    integer,protected:: errLineNum = 0
    integer,private:: nUserErr = 0
    character(len=4),protected:: errHeap(errHeapSize) = ''
    character(len=errMsgLen),protected::  errMsg  = ''
    character(len=errMsgLen),protected::  suppErrMsg  = ''
    character(len=errFileLen),protected:: errFileName = ''
    character(len=errHeapSize*5),protected:: errHeapMsg = ''
    logical,protected:: errON(maxErrCode) = .true.
    integer,protected:: nErrHit(maxErrCode) = 0

    type aUserErr
        integer:: code
        character(len=errMsgLen):: msg
    end type aUserErr
    type(aUserErr),private:: userErr(maxUserErr)

    contains
        subroutine SetSuppErrMsg(suppMsg)
            ! 2019/01/06
            implicit none
            ! input
            character(len=*):: suppMsg
            suppErrMsg = ' '//trim(suppMsg)
            return
        end subroutine SetSuppErrMsg

        subroutine SwitchError(thisCode,isEnabled)
            ! 2017/09/24
            implicit none
            ! input/output
            integer,intent(in):: thisCode
            logical,optional,intent(in):: isEnabled
            ! local
            ! code
            if ((thisCode.gt.maxErrCode).or.(thisCode.lt.0)) return   ! do nothing in case of invalid errCode
            if (present(isEnabled)) then
                errON(thisCode) = isEnabled
            else
                errON(thisCode) = .not.errON(thisCode)
            endif
            return
        end subroutine SwitchError

        subroutine IncrementError(thisCode)
            ! 2017/09/24
            implicit none
            ! in
            integer,intent(in):: thisCode
            ! local
            ! code
            if ((thisCode.gt.maxErrCode).or.(thisCode.lt.0)) return   ! do nothing in case of invalid errCode
            nErrHit(thisCode) = nErrHit(thisCode)+1
            return
        end subroutine IncrementError


        subroutine RaiseWarning(warnMsg,warnSub)
            !--------------------------------------------------------------------------------------------------------------------------
            !   Instruction:
            !       This subroutine processes warnings, generally warning message is displayed to stdout and/or log immediately.
            !   Date:
            !       2013-09-25
            !   Inputs:
            !       warnMsg     str         warning message
            !       warnSub     str         the function ID in which this subroutine is called
            !--------------------------------------------------------------------------------------------------------------------------
            implicit none
            ! Declare input parameters:
            character(len=*) warnMsg,warnSub
            nWarning = nWarning+1
            write(*,1) trim(adjustl(warnSub)),trim(adjustl(warnMsg))
  1         format('WARNING[',a,']: ',a)
            return
        end subroutine RaiseWarning

        subroutine RaiseError(iErr,errSub,errLine,errFile)
            !--------------------------------------------------------------------------------------------------------------------------
            !   Instruction:
            !       This subroutine processes newly raised error - update error code and the error heap.
            !   Date:
            !       2013-09-24
            !   Update:
            !       J. Tang @ 2016-07-15    added optional parameters of errLine and errFile
            !   Inputs:
            !       iErr        i           error ID number
            !       errSub      str         the function ID in which this subroutine is called
            !--------------------------------------------------------------------------------------------------------------------------
            implicit none
            ! Declare input parameters:
            integer,intent(in):: iErr
            character(len=*),intent(in):: errSub
            integer,intent(in),optional:: errLine
            character(len=*),intent(in),optional:: errFile
            ! Declare local variables:
            integer:: currentLen

            ! Update the errCode only when iErr is positive, which represents the occurrance of a new error.
            errCode = iErr
            call UpdateErrMsg(errCode,errMsg)
            if (present(errLine)) errLineNum = errLine
            if (present(errFile)) errFileName = errFile

            return
        end subroutine RaiseError

        subroutine CheckIn(inSub)
            !--------------------------------------------------------------------------------------------------------------------------
            !   Instruction:
            !       This subroutine checks in a sub.
            !   Date:
            !       2013-11-18
            !   Inputs:
            !       inSub       str         the checking-in sub
            !--------------------------------------------------------------------------------------------------------------------------
            implicit none
            character*4:: inSub
            if (errCode.gt.0) then
                call PrintErrorInfo
                print*,  ' >> Emergent abort due to an unexpected error upon check-in!'
            endif
            if (errHeapTop.ge.errHeapSize) then
                call RaiseError(0003,'ckin')
                call PrintErrorInfo
                print*,  ' >> Abort due to an error !'
            else
                errHeapTop = errHeapTop+1
                errHeap(errHeapTop) = inSub
            endif
            return
        end subroutine CheckIn

        subroutine CheckOut(outSub)
            !--------------------------------------------------------------------------------------------------------------------------
            !   Instruction:
            !       This subroutine checks out a sub.
            !   Date:
            !       2013-11-18
            !   Inputs:
            !       outSub      str         the checking-out sub
            !--------------------------------------------------------------------------------------------------------------------------
            implicit none
            character*4:: outSub
            if (errHeapTop.le.0) then
                call PrintErrorInfo
                print*,  ' >> Emergent abort due to an unexpected error upon check-out !'
            elseif (outSub.ne.errHeap(errHeapTop)) then
                call RaiseError(0005,'ckot')
                call PrintErrorInfo
                print*,  ' >> Abort due to an error !'
            else
                errHeap(errHeapTop) = ''
                errHeapTop = errHeapTop-1
            endif
            return
        end subroutine CheckOut

        subroutine UpdateErrMsg(iCode, thisMsg)
            !--------------------------------------------------------------------------------------------------------------------------
            !   Instruction:
            !       This subroutine updates the error message according to given error code.
            !   Date:
            !       2013-09-24
            !   Inputs:
            !       iCode       i           error code
            !   Outputs:
            !       thisMsg      str         error message
            !--------------------------------------------------------------------------------------------------------------------------
            implicit none
            ! Declare input parameters:
            integer,intent(in):: iCode
            ! Declare output parameters:
            character(len=*),intent(out):: thisMsg
            ! local
            integer:: j

            ! - Suggest errors be catalogued by the first digit of the 4-/5-digit error code.
            ! - Recommend error strings to use full sentence, not a phrase.
            ! - Recommend error string to directly describe the error, not indication of what is expected.
            !   Prefer: X is 0/X is not non-zero. (direct description of error).
            !   Not prefer: X should be non-zero. (expectation, indirect)
            select case (iCode)
            ! 0 stands for internal error, for debugging/system purpose.
            ! 0001 is the test error.
                case (0001)
                    thisMsg = 'No error actually occurred. It is a test error.'
                case (0002)
                    thisMsg = 'An unimplemented function has been called.'
                case (0003)
                    thisMsg = 'Error-heap overflows.'
                case (0004)
                    thisMsg = 'Popping empty error-heap is attempted.'
                case (0005)
                    thisMsg = 'Subname to be popped does not match error-heap top.'
                case (0006)
                    thisMsg = 'A conflict in configuration is detected.'
                case (0007)
                    thisMsg = 'Not all compulsory arguments are available.'
                case (0008)
                    thisMsg = 'Fail to read internal/external file.'
                case (0009)
                    thisMsg = 'Requested variable is not found in the given VariableModel.'
                case (0010)
                    thisMsg = 'An unimplemented TypeID(1) is specified.'
                case (0011)
                    thisMsg = 'An unimplemented TypeID(2) is specified.'
                case (0012)
                    thisMsg = 'An unaccessible logic is detected.'
                case (0013)
                    thisMsg = 'Fail to allocate/deallocate memory.'
                case (0014)
                    thisMsg = 'An obsolete function has been called.'
                case (0015)
                    thisMsg = 'An unexpected end of file is encountered.'
                case (0016)
                    thisMsg = 'Subroutine/function arguments are invliad or in conflict.'
                case (0017)
                    thisMsg = 'A conflict in setup/configuration is detected.'
                case (0018)
                    thisMsg = 'Program command-line arguments are in conflict.'
                case (0019)
                    thisMsg = 'One or more program command-line arguments are invalid.'
                case (0020)
                    thisMsg = 'File to open does not exist.'
            ! 1 stands for common maths error.
                case (1001)
                    thisMsg = 'Zero is divided.'
                case (1002)
                    thisMsg = 'Iteration exceeds maximum before convergence.'
                case (1003)
                    thisMsg = 'Matrix is not positive definite.'
                case (1004)
                    thisMsg = 'Square root of a negative number is computed.'
                case (1005)
                    thisMsg = 'Matrix is not square.'
                case (1006)
                    thisMsg = 'A vector has all three components of naught.'
                case (1007)
                    thisMsg = 'The x/y vectors for interpolation do not have consistent dimension.'
                case (1008)
                    thisMsg = 'The nodes for interpolation have the same x values.'
                case (1009)
                    thisMsg = 'Array dimensions do not match.'
                case (1010)
                    thisMsg = 'Negative variance is found.'
                case (1011)
                    thisMsg = 'Fail to solve the large sparse linear equation, due to invalid dimension (0) specified.'
                case (1012)
                    thisMsg = 'Fail to solve the large sparse linear equation, due to a null row.'
                case (1013)
                    thisMsg = 'Fail to solve the large sparse linear equation, due to duplicate rows in the matrix.'
                case (1014)
                    thisMsg = 'Fail to solve the large sparse linear equation, due to a zero-pivot row.'
                case (1015)
                    thisMsg = 'Fail to solve the large sparse linear equation, due to storage exceeding limit.'
                case (1016)
                    thisMsg = 'Matrix is singular.'
            ! 2 stands for errors concerning state vectors.
                case (2001)
                    thisMsg = 'Semi-major axis is not positive.'
                case (2002)
                    thisMsg = 'Eccentricity is not between 0 (inclusive) and 1 (exclusive).'
                case (2003)
                    thisMsg = 'Specified state type is invalid.'
                case (2004)
                    thisMsg = 'Unknown operation is specified.'
                case (2005)
                    thisMsg = 'Unknown element is specified.'
                case (2006)
                    thisMsg = 'Solving for eccentricity failed, due to invalid xi and eta.'
            ! 3 stands for errors in time and reference system process.
                case (3001)
                    thisMsg = 'Unknown time system is specified.'
                case (3002)
                    thisMsg = 'Size of EOP table exceeds maximum.'
                case (3003)
                    thisMsg = 'Size of leap second table exceeds maximum.'
                case (3004)
                    thisMsg = 'Given UTC is earlier than the first entry of leap second table.'
                case (3005)
                    thisMsg = 'Given UTC is beyond the range in the EOP table.'
                case (3006)
                    thisMsg = 'Unknown EOP model is specified.'
                case (3007)
                    thisMsg = 'Fro/To frame is not ''inertial''/''fixed''.'
                case (3008)
                    thisMsg = 'Unknwon body ID is specified for frame conversion.'
                case (3009)
                    thisMsg = 'Failed in applying memory address for storing rotation series.'
                case (3010)
                    thisMsg = 'Invalid axis is specified for a rotation.'
                case (3011)
                    thisMsg = 'An unexpected end of rotation series is encountered.'
                case (3012)
                    thisMsg = 'Leaptab is called before being initialized.'
                case (3013)
                    thisMsg = 'Required degree of nutation series is smaller than that of given IAU1980 series.'
                case (3014)
                    thisMsg = 'Invalid ID is specified for a celestial body.'
                case (3015)
                    thisMsg = 'Invalid trajectory type is specified for current celestial body.'
                case (3016)
                    thisMsg = 'Cannot sync time since EOP/Leap is not initialized.'
                case (3017)
                    thisMsg = 'Cell frame conversions exceed limit.'
                case (3018)
                    thisMsg = 'Failed in applying memory address for storing frame conversion cache.'
                case (3019)
                    thisMsg = 'Required frame conversion cannot be implemented with given cell conversions.'
                case (3020)
                    thisMsg = 'Conversion has not be defined for given cell.'
                case (3021)
                    thisMsg = 'Failed in applying memory address for storing added cell conversion to the cell sequence.'
                case (3022)
                    thisMsg = 'Failed in applying memory address for storing cell conversion subroutine pointers.'
                case (3023)
                    thisMsg = 'Invalid type is specified when requesting for lunisolar parameters.'
                case (3024)
                    thisMsg = 'Cannot convert week/sec to day/sec without a valid starting epoch for an unknown system.'
            ! 4 stands for the errors when processing (read/write) external files.
                case (4001)
                    thisMsg = 'No EOP records are read.'
                case (4002)
                    thisMsg = 'No leap second records are read.'
                case (4003)
                    thisMsg = 'Cannot open file due to unspecific cause.'
                case (4004)
                    thisMsg = 'Cannot close file due to unspecific cause.'
                case (4005)
                    thisMsg = 'Opening an inexistent file is attempted.'
                case (4006)
                    thisMsg = 'Cannot read all the requied gravity coefficients in gravity model file.'
                case (4007)
                    thisMsg = 'Number of dependent files exceeds uplimit.'
                case (4008)
                    thisMsg = 'Cannot read one of lines in dependent_file_directories file.'
                case (4009)
                    thisMsg = 'The required maxmum order or degree need to read exceeds general setting.'
                case (4010)
                    thisMsg = 'File name for the JPLDE exceeds 80 characters.'
                case (4011)
                    thisMsg = 'Specified JPLDE does not exist.'
                case (4012)
                    thisMsg = 'Specified record length (RECL) for JPLDE is unacceptable.'
                case (4013)
                    thisMsg = 'Specified KSIZE for JPLDE is unacceptable.'
                case (4014)
                    thisMsg = 'Specifier FSIZER for JPLDE is unacceptable.'
                case (4015)
                    thisMsg = 'Neither KSIZE or FSIZER is specified for JPLDE.'
                case (4016)
                    thisMsg = 'KSIZE is not specified for JPLDE when FSIZER is not 3.'
                case (4017)
                    thisMsg = 'JPLDE is requested before being properly initialized.'
                case (4018)
                    thisMsg = 'Cannot request new FileUnit since maximum files are opened.'
                case (4019)
                    thisMsg = 'Failed to request a unit with given base'
                case (4020)
                    thisMsg = 'The unit to be closed has no match in the opened files.'
                case (4021)
                    thisMsg = 'Unsupported JPLDE is specified for SetDE2.'
                case (4022)
                    thisMsg = 'Cannot write file due to unspecific cause.'
                case (4023)
                    thisMsg = 'Cannot read file due to unspecific cause.'
                case (4024)
                    thisMsg = 'Specified JPL DE version does not match the trajectory file.'
                case (4025)
                    thisMsg = 'Cannot find the requested keyword in input.'
                case (4026)
                    thisMsg = 'Log file is not ready.'
            ! 5 stands for the error in propagation
                case (5001)   ! adapted and used
                    thisMsg = 'Unknown integrator is specified.'
                case (5002)   ! adapted and used
                    thisMsg = 'Initial guess of ACC/ATT interpolation location is not proper.'
                case (5003)   ! adapted and used
                    thisMsg = 'Cannot set the active OPModel.'
                case (5004)   ! adapted and used
                    thisMsg = 'The propagation model has not been properly initialized.'
                case (5005)   ! adapted and used
                    thisMsg = 'Variable dimension changes between epochs when multi-step numerical method is employed.'
                case (5006)   ! adapted and used
                    thisMsg = 'Non-uniform steps or non-consecutive epochs are found in multi-step numerical method.'
                case (5007)
                    thisMsg = 'There is no gravity model corresponds to the specified gravity body.'
                case (5008)
                    thisMsg = "Celestrail trajectory model haven't been asigned for non-central gravity body."
                case (5009)
                    thisMsg = 'There is no celestrial trajectory model corresponds to the specified gravity body.'
                case (5010)
                    thisMsg = 'Integrator has not been set when integrating procedure starts.'
                case (5011)
                    thisMsg = 'Initial stepsize for integrator must be positive.'
                case (5012)
                    thisMsg = "Force model hasn't been set before integration."
                case (5013)   ! tang's def below
                    thisMsg = 'Fail to find an optimal stepsize in a variable-stepsize integrator.'
                case (5014)
                    thisMsg = 'Truncation error exceeds given limit for fixed stepsize.'
                case (5015)
                    thisMsg = 'Fail to initialize YHC Variable.'
                case (5016)
                    thisMsg = 'Fail to initialize OD Variable.'
                case (5017)
                    thisMsg = 'Input VAR is not specified for YHC.'
                case (5018)
                    thisMsg = 'Invalid location is specified for trajectory storage.'
                case (5019)
                    thisMsg = 'Consecutive trajectory entries have inconsistent dimensions.'
                case (5020)
                    thisMsg = 'Unknown type of orbit analytical solution is specified.'
                case (5021)
                    thisMsg = 'Cannot find supported state variables in YHC Variable.'
                case (5022)
                    thisMsg = 'State type is not properly configured in the VAR.'
            ! 6 stands for the errors regarding various types of observations.
                case (6001)
                    thisMsg = 'Iteration exceeds limit when solving for TOA/range measurement.'
                case (6002)
                    thisMsg = 'The number of observations exceeds the limit.'
                case (6003)
                    thisMsg = 'Unsupported type of observations is found.'
                case (6004)
                    thisMsg = 'An unrecognizable satellite/station is related to the observations.'
                case (6005) ! new
                    thisMsg = 'Incompatible measurements are specified at the same time.'
                case (6006) ! new
                    thisMsg = 'An invalid time flag is specified.'
                case (6007) ! new
                    thisMsg = 'Available Cartesians of the target/tracker are not consistent in coordinate system.'
                case (6008) ! new
                    thisMsg = 'Requested measurement and given Cartesian are not compatible.'
                case (6009) ! new
                    thisMsg = 'An invalid equator type is specified.'
                case (6010) ! new
                    thisMsg = 'An invalid link type is specified.'
                case (6011) ! new
                    thisMsg = 'Observations are not ordered chronologically.'
                case (6012)
                    thisMsg = 'Optional input variables must appear in pears.'
            ! 7 stands for the errors regarding force models.
                case (7001)
                    thisMsg = 'The specified gravity field degree/order is invalid.'
                case (7002)
                    thisMsg = 'The specified gravity field degree is smaller than its order.'
                case (7003)
                    thisMsg = 'The density model is not correctly set before use.'
                case (7004)
                    thisMsg = 'A zero or negative altitude is requested for reference density.'
                case (7005)
                    thisMsg = 'Invalid degree/order is requested for solid Earth tide.'
                case (7006)
                    thisMsg = 'Invalid degree/order is requested for ocean tide.'
                case (7007)
                    thisMsg = 'Invalid tide type is specified when requesting for spharm correction.'
            ! 8 stands for the errors regarding satellite models
                case (8001)
                    thisMsg = 'The number of satellites exceeds the limit.'
            ! 9 stands for errors regarding orbit determination
                case (9001)
                    thisMsg = 'PODVar elements are inconsistent with allocated PODModel.'
                case (9002)
                    thisMsg = 'Invalid evalType is specified.'
                case (9003)
                    thisMsg = 'Observations fed to the sequential estimator have different epochs.'
                case (9004)
                    thisMsg = 'Unknown estimator is specified for POD.'
                case (9005)
                    thisMsg = 'No related device is found when collecting variables to propagate.'
                case (9006)
                    thisMsg = 'Not all observations are associated with proper rejection models.'
                case (9007)
                    thisMsg = 'Variables of PODVar have inconsistent epochs.'
            ! 10 stands for internal storage exploration and release.
                case (10001)
                    thisMsg = 'Cannot explore requied internal storage.'
                case (10002)
                    thisMsg = 'Cannot release explored internal storage.'
            ! 11 stands for the errors during process of GRACE data
                case (11002)
                    thisMsg = 'Satellite names in GRACE L1B are inconsistent.'
                case (11003)
                    thisMsg = 'Unexpected value is found in the GRACE L1B header.'
                case (11004)
                    thisMsg = 'Unrecognizable record is found in the GRACE L1B file.'
            ! 51 stands for the matlab plot
                case (51001)
                    thisMsg = 'Cannot get matlab engine open.'
                case (51002)
                    thisMsg = 'Cannot implant created array in matlab engine.'
                case (51003)
                    thisMsg = 'Cannot get current work directory.'
                case (51004)
                    thisMsg = 'Error in inporting and runing command into matlab engine.'
                case (51005)
                    thisMsg = 'Cannot get matlab engine close.'
                case (51006)
                    thisMsg = 'There is no line data and properties has been added.'
                case (51007)
                    thisMsg = 'Not all the added lines have been loaded.'
                case (51008)
                    thisMsg = 'Matlab package cannot run in non-Win32 system.'
            ! 61 stands for RINEX routines
                case (61001)
                    thisMsg = 'Failed to interpret a description.'
                case (61002)
                    thisMsg = 'Specific wavelength factor exceeds limit.'
                case (61003)
                    thisMsg = 'OBS type exceeds limit.'
                case (61004)
                    thisMsg = 'Actual # of OBSERV exceeds specified value.'
                case (61005)
                    thisMsg = '# of OBS TYPE is missing in its first entry.'
                case (61006)
                    thisMsg = '# of SATELLTES exceeds limit.'
                case (61007)
                    thisMsg = 'Actual # of SATELLTES exceeds specified value.'
                case (61008)
                    thisMsg = 'Actual # of OBS TYPE per PRN exceeds # of OBS TYPE.'
                case (61009)
                    thisMsg = 'The RINEX file does not start with RINEX VERSION / TYPE label.'
                case (61010)
                    thisMsg = 'Type of file is invalid.'
                case (61011)
                    thisMsg = 'Type of OBS is invalid.'
                case (61012)
                    thisMsg = 'Default Wavelength Fact is invalid.'
                case (61013)
                    thisMsg = '# of PRN per specific WLFact is invalid.'
                case (61014)
                    thisMsg = 'Satellite system in specific WLFact is invalid.'
                case (61015)
                    thisMsg = 'Specified # of PRN per specific WLFact exceeds format limit.'
                case (61016)
                    thisMsg = 'PRN number in specific WLFact is invalid or exceeds limit.'
                case (61017)
                    thisMsg = 'Satllite system for OBS file is invalid.'
                case (61018)
                    thisMsg = 'Explicit time system is missing for mixed satellite system.'
                case (61019)
                    thisMsg = 'Time system and satellite system do not match.'
                case (61020)
                    thisMsg = 'The RINEX file is not of a supported version.'
                case (61021)
                    thisMsg = 'The RINEX file is not of a supported type.'
                case (61022)
                    thisMsg = 'Failed to interpret an epoch in the data segment.'
                case (61023)
                    thisMsg = 'SatNum in OBS data exceeds SatNum in Epoch.'
            ! 62 stands for sp3 format routines
                case (62001)
                    thisMsg = 'Failed to read in the sp3c Header.'
                case (62002)
                    thisMsg = 'Failed to intepret the sp3c Header.'
                case (62003)
                    thisMsg = 'The sp3 file is not of a supported version.'
            ! unknown error
                case default
                    thisMsg = 'Unknown error.'
                    do j = 1,nUserErr
                        if (iCode.eq.userErr(j)%code) then
                            thisMsg = trim(userErr(j)%msg)
                            exit
                        endif
                    enddo
            end select
        end subroutine UpdateErrMsg

        subroutine UpdateErrHeapMsg
            !--------------------------------------------------------------------------------------------------------------------------
            !   Instruction:
            !       This subroutine updates the error-heap string.
            !   Date:
            !       2013-11-18
            !   Inputs:
            !       errHeap     str(4,errHeapSize)  error-heap, module variable
            !       errHeapTop  i                   error-heap top, module variable
            !   Outputs:
            !       errHeapMsg  str                 error-heap message, module variable
            !--------------------------------------------------------------------------------------------------------------------------
            implicit none
            integer:: i
            ! FILO for heap
            errHeapMsg = ''
            do i = 1, errHeapTop-1
                errHeapMsg((i-1)*5+1:i*5) = errHeap(i)//'>'
            enddo
            if (errHeapTop.gt.0) then
                errHeapMsg(errHeapTop*5-4:errHeapTop*5) = errHeap(i)//'.'
            else
                errHeapMsg = 'EMPTY!'
            endif
            return
        end subroutine UpdateErrHeapMsg

!        subroutine PrintErrorInfo
!            implicit none
!            integer:: j
!            call UpdateErrHeapMsg
!            if (errCode.gt.0) then
!                print 1,trim(adjustl(errHeapMsg))
!                print '(1x,a)','Errors detected:'
!                print 2,errCode
!                print 3,trim(errMsg)//trim(suppErrMsg)
!                if ((trim(errFileName).ne.'').or.(errLineNum.ne.0)) print 4
!                if (trim(errFileName).ne.'') print 5,trim(errFileName)
!                if (errLineNum.ne.0) print 6,errLineNum
!            else
!                if (count(nErrHit.gt.0).gt.0) then
!                    print '(1x,a)','Errors recorded:'
!                    do j = 1,maxErrCode
!                        if (nErrHit(j).gt.0) then
!                            print 7,j,nErrHit(j)
!                        endif
!                    enddo
!                else
!                    print '(1x,a)','No errors detected.'
!                endif
!            endif
! 1          format(' Current calling chain: ',a)
!! 2          format(' - Error code: ',i6)
! 2          format(' - Error code: ',i<idint(dlog10(dble(errCode)))+1>)
! 3          format(' - Error message: ',a)
! 4          format(' Approximate error location:')
! 5          format(' - File: ',a)
! 6          format(' - Line: ',i<idint(dlog10(dble(errLineNum)))+1>)
! 7          format(' - #',i<idint(dlog10(dble(j)))+1>,' is hit ',i<idint(dlog10(dble(nErrHit(j))))+1>,' times.')
!            return
!        end subroutine PrintErrorInfo
        
        subroutine PrintErrorInfo
        implicit none
        integer:: j
        character(20) :: fmt_str  ! ���ڶ�̬������ʽ�ַ���
    
        call UpdateErrHeapMsg
        if (errCode.gt.0) then
            print 1, trim(adjustl(errHeapMsg))
            print '(1x,a)', 'Errors detected:'
        
            ! ʹ��I0��ʽ�Զ�������������
            print '(" - Error code: ",I0)', errCode
        
            print 3, trim(errMsg)//trim(suppErrMsg)
            if ((trim(errFileName).ne.'').or.(errLineNum.ne.0)) print 4
            if (trim(errFileName).ne.'') print 5, trim(errFileName)
            if (errLineNum.ne.0) print '(" - Line: ",I0)', errLineNum
        else
            if (count(nErrHit.gt.0).gt.0) then
                print '(1x,a)', 'Errors recorded:'
                do j = 1, maxErrCode
                    if (nErrHit(j).gt.0) then
                        ! ʹ��I0��ʽ�Զ�������������
                        print '(" - #",I0," is hit ",I0," times.")', j, nErrHit(j)
                    endif
                enddo
            else
                print '(1x,a)', 'No errors detected.'
            endif
        endif

1       format(' Current calling chain: ',a)
3       format(' - Error message: ',a)
4       format(' Approximate error location:')
5       format(' - File: ',a)
    
        return
        end subroutine PrintErrorInfo
        

        !subroutine DefineUserError(code,msg)
        !    !--------------------------------------------------------------------------------------------------------------------------
        !    !   Instruction:
        !    !       This subroutine adds a new user-defined error.
        !    !   Date:
        !    !       2019-04-27
        !    !   Inputs:
        !    !       code        i           error code
        !    !       msg         str         error message
        !    !   Outputs:
        !    !       No explicit output, the private list is extended.
        !    !--------------------------------------------------------------------------------------------------------------------------
        !    implicit none
        !    integer:: code
        !    character(len=*):: msg
        !    if (nUserErr.ge.maxUserErr) then
        !        write(*,'(1x,"User error #",i<idint(dlog10(dble(code)))+1>," is not added. User-defined errors are full. Increase the hard-coded limit and try again.")') code
        !    elseif (code.lt.90000.or.code.gt.99999) then
        !        write(*,'(1x,"User error #",i<idint(dlog10(dble(code)))+1>," is not added. User-defined errors should be numbered between 90000 and 99999.")') code
        !    else
        !        nUserErr = nUserErr+1
        !        userErr(nUserErr)%code = code
        !        userErr(nUserErr)%msg  = trim(msg)
        !    endif
        !    return
        !end subroutine DefineUserError

        subroutine DefineUserError(code, msg)
            !--------------------------------------------------------------------------------------------------------------------------
            !   Instruction:
            !       This subroutine adds a new user-defined error.
            !   Date:
            !       2019-04-27
            !   Inputs:
            !       code        i           error code
            !       msg         str         error message
            !   Outputs:
            !       No explicit output, the private list is extended.
            !--------------------------------------------------------------------------------------------------------------------------
            implicit none
            integer:: code
            character(len=*):: msg
    
            if (nUserErr.ge.maxUserErr) then
                write(*, '(" User error #", I0, " is not added.")') code
            elseif (code.lt.90000 .or. code.gt.99999) then
                write(*, '(" User error #", I0, " is not added.")') code
            else
                nUserErr = nUserErr + 1
                userErr(nUserErr)%code = code
                userErr(nUserErr)%msg = trim(msg)
            endif
    
            return
        end subroutine DefineUserError


end module ErrorHandler
