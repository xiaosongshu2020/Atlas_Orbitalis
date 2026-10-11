module EOPParam
    !----------------------------------------------------------------------------------------------------------------------------------
    !   Instruction:
    !       This module defines the EOP parameters.
    !   Date:
    !       2013-09-25
    !   Contains:
    !       AddEOP      sub
    !       AddLeap     sub
    !       SyncEOP     sub
    !----------------------------------------------------------------------------------------------------------------------------------
    implicit none

    integer,parameter:: maxEOP = 100000
    integer,parameter:: maxLeap = 50

    type EOPTable
        double precision:: T
        double precision:: Xp,Yp,DUT1
    end type EOPTable
    type(EOPTable),protected:: eopTab(maxEOP)
    integer,protected:: nEOP=0

    type LeapTable
        integer:: yr,mo,day,leap
        double precision:: mjd
    end type LeapTable
    type(LeapTable),protected:: leapTab(maxLeap)
    integer,protected:: nLeap=0

    integer,protected:: EOPModel=1    ! 1 for IAU1976/1980; 2 for IAU2000
    integer,protected:: Leap,iLeap
    double precision,protected:: Xp, Yp, DUT1   ! xp,yp->radians, dut1->second
    double precision,protected:: Xp_rate, Yp_rate, DUT1_rate     ! xp_rate,yp_rate->rad/day, dut1_rate->sec/day
    logical,protected:: isEOPSet=.false., isLeapSet=.false.
    logical,protected:: isNoEOPaccept=.true. ! if the Default value of EOP and leap second data can be accepted when the data file cannot cover all the times needed.

    contains
        subroutine AddEOP(utc,xp,yp,dut1)
        !--------------------------------------------------------------------------------------------------------------------------
        !   Instruction:
        !       This subroutine adds a new entry in EOP table, supposed to be called once and once only during initialization.
        !   Date:
        !       2013-09-25
        !   Inputs:
        !       utc         d           utc of the entry
        !       xp          d           x component of polar motion
        !       yp          d           y component of polar motion
        !       dut1        d           dut1=ut1-utc
        !   Outputs:
        !       No explicit output parameters, yet the EOP table is updated.
        !--------------------------------------------------------------------------------------------------------------------------
            use ErrorHandler
            implicit none
            ! Declare the input parameters:
            double precision,intent(in):: utc,xp,yp,dut1
            ! Declare local variables:
            character(len=4):: subName='adeo'
            logical,save:: warnSize=.true.

            call checkin(subName)
            if (nEOP.ge.maxEOP) then
                if (warnSize) then
                    call RaiseWarning('Size of the EOP table exceeds maximum. No more data will be accepted.',subName)
                    warnSize = .false.
                endif
            else
                nEOP = nEOP+1
                eopTab(nEOP)%T  = utc
                eopTab(nEOP)%Xp = xp
                eopTab(nEOP)%Yp = yp
                eopTab(nEOP)%DUT1 = dut1
            endif
            call checkout(subName)

            return
        end subroutine AddEOP

        subroutine AddLeap(yr,mo,day,leap)
        !--------------------------------------------------------------------------------------------------------------------------
        !   Instruction:
        !       This subroutine adds a new entry in Leap table, supposed to be called once and once only during initialization.
        !   Date:
        !       2013-09-25
        !   Inputs:
        !       yr          i           year of the leap second
        !       mo          i           month of the leap second
        !       day         i           day of the leap second
        !       leap        i           accumulated leap second till then
        !   Outputs:
        !       No explicit output parameters, yet the Leap table is updated.
        !--------------------------------------------------------------------------------------------------------------------------
            use ErrorHandler
            implicit none
            ! Declare the input parameters:
            integer,intent(in):: yr,mo,day,leap
            ! Declare local variables:
            character(len=4):: subName='adls'
            logical,save:: warnSize=.true.
            ! Declare external functions:

            call checkin(subName)
            if (nLeap.ge.maxLeap) then
                if (warnSize) then
                    call RaiseWarning('Size of the leap second table exceeds maximum. No more data will be accepted.',subName)
                    warnSize = .false.
                endif
            else
                nLeap = nLeap+1
                leapTab(nLeap)%yr   = yr
                leapTab(nLeap)%mo   = mo
                leapTab(nLeap)%day  = day
                leapTab(nLeap)%leap = leap
                leapTab(nLeap)%mjd  = EncodeMJD(yr,mo,day,0,0,0d0)
            endif

            call checkout(subName)
            return
        end subroutine AddLeap

        subroutine SyncLeapSec(utc)
        !--------------------------------------------------------------------------------------------------------------------------
        !   Instruction:
        !       This subroutine synchronize leap second data.
        !   Date:
        !       2013-09-24
        !   Update:
        !       2015-01-27
        !       2017-10-15 bug fix for that iLeap is not set when the epoch is later than the last leap.
        !   Inputs:
        !       utc         d           current utc
        !   Outputs:
        !       No explicit outputs, leap data will be updated with respect of input time.
        !--------------------------------------------------------------------------------------------------------------------------
            use ErrorHandler
            implicit none
            double precision,intent(in)::utc
            integer i
            character(len=4)::subname='synl'
            call checkin(subname)
            
            !if LeapSec parameters are not been initialized, initialized it
            if ((.not. isEOPSet) .or. (.not. isLeapSet)) then
                call EOPLeapInit('data/AtlasData/EOP20170101_20230426.txt', 'data/AtlasData/leap.txt')
            endif
            
            ! find current leap second
            if (isLeapSet.and.(nLeap.gt.0)) then    ! find only when leaptab is filled
                if (utc.lt.leapTab(1)%mjd) then
                    call RaiseError(3004,subName)
                    return
                elseif (utc.ge.leapTab(nLeap)%mjd) then
                    Leap = leapTab(nLeap)%leap
                    iLeap = nLeap
                else
                    do i = 1,nLeap-1
                        if (utc.ge.leapTab(i)%mjd.and.utc.lt.leapTab(i+1)%mjd) then
                          Leap = leapTab(i)%leap
                          ileap=i
                          exit
                        endif
                    enddo
                endif
            else     ! Leap is used to compute TT, it will be required in any case. Raise error is leaptab is empty.
                call RaiseError(3012,subName)
                return
            endif

            call checkout(subname)
        end subroutine SyncLeapSec

        subroutine SyncEOP(utc)
        !--------------------------------------------------------------------------------------------------------------------------
        !   Instruction:
        !       This subroutine synchronize xp, yp and dUT1 from EOP table.
        !       The input time SHOULD BE UTC. But with the fact that the change of xp, yp, and dUT1 cannot exceed 0.005 (arcsec/second) per day,
        !       the input time COULD BE UT1, TT and TDB instead, with ERROR no lager than 5d-6 arcsec/second ( take the largest input time error as 0.001 day).
        !   Date:
        !       2013-09-24
        !   Update:
        !       2015-01-27
        !   Inputs:
        !       utc         d           current utc
        !   Outputs:
        !       No explicit output. xp, yp and dUT1 are updated with respect of input time.
        !--------------------------------------------------------------------------------------------------------------------------
            use ErrorHandler
            implicit none
            ! Declare input parameters:
            double precision,intent(in):: utc
            ! Declare local variables:
            integer:: i
            character(len=4):: subName='syne'
            ! Declare external functions:

            call checkin(subName)
            
            !if EOP parameters are not been initialized, initialized it
            if ((.not. isEOPSet) .or. (.not. isLeapSet)) then
                call EOPLeapInit('data/AtlasData/EOP.txt', 'data/AtlasData/leap.txt')
            endif
            
            ! find current xp, yp and dut1
            if (nEOP.gt.0) then    ! find when eopTab is not empty
                if (utc.lt.eopTab(1)%T.or.utc.gt.eopTab(nEOP)%T) then
                    if(isNoEOPaccept)then
                        xp   = 0d0
                        yp   = 0d0
                        dut1 = 0d0
                        call checkout(subname)
                    else
                        call RaiseError(3005,subName)
                    endif
                    return
                endif
                i = idint(utc-eopTab(1)%T) + 1
                if (utc.ge.eopTab(i)%T.and.utc.le.eopTab(i+1)%T) then    ! test to directly locate the interpolation
                    xp   = LinInt(eopTab(i)%T,eopTab(i)%xp,eopTab(i+1)%T,eopTab(i+1)%xp,utc)
                    if (errCode.gt.0) return
                    ! If no error occurs here, it is unlikely to happen in the following two cases.
                    yp   = LinInt(eopTab(i)%T,eopTab(i)%yp,eopTab(i+1)%T,eopTab(i+1)%yp,utc)
                    dut1 = LinInt(eopTab(i)%T,eopTab(i)%dut1,eopTab(i+1)%T,eopTab(i+1)%dut1,utc)
                    xp   = xp/206264.806246D0   ! mas to rad
                    yp   = yp/206264.806246D0   ! mas to rad
                    xp_rate   = (eopTab(i+1)%xp-eopTab(i)%xp)/206264.806246D0/(eopTab(i+1)%T-eopTab(i)%T)    ! change rate of xp, in rad/day
                    yp_rate   = (eopTab(i+1)%yp-eopTab(i)%yp)/206264.806246D0/(eopTab(i+1)%T-eopTab(i)%T)    ! change rate of yp, in rad/day
                    dut1_rate = (eopTab(i+1)%dut1-eopTab(i)%dut1)/(eopTab(i+1)%T-eopTab(i)%T)  ! change rate of dut1 in sec/day
                else                    ! cannot directly locate, sequential search
                    do i = 1,nEOP-1
                        if (utc.ge.eopTab(i)%T.and.utc.le.eopTab(i+1)%T) then
                            xp   = LinInt(eopTab(i)%T,eopTab(i)%xp,eopTab(i+1)%T,eopTab(i+1)%xp,utc)
                            if (errCode.gt.0) return
                            ! If no error occurs here, it is unlikely to happen in the following two cases.
                            yp   = LinInt(eopTab(i)%T,eopTab(i)%yp,eopTab(i+1)%T,eopTab(i+1)%yp,utc)
                            dut1 = LinInt(eopTab(i)%T,eopTab(i)%dut1,eopTab(i+1)%T,eopTab(i+1)%dut1,utc)
                            xp   = xp/206264.806246D0   ! mas to rad
                            yp   = yp/206264.806246D0   ! mas to rad
                            xp_rate   = (eopTab(i+1)%xp-eopTab(i)%xp)/206264.806246D0/(eopTab(i+1)%T-eopTab(i)%T)    ! change rate of xp, in rad/day
                            yp_rate   = (eopTab(i+1)%yp-eopTab(i)%yp)/206264.806246D0/(eopTab(i+1)%T-eopTab(i)%T)    ! change rate of yp, in rad/day
                            dut1_rate = (eopTab(i+1)%dut1-eopTab(i)%dut1)/(eopTab(i+1)%T-eopTab(i)%T)  ! change rate of dut1 in sec/day
                            exit
                        endif
                    enddo
                endif

            endif    ! Xp/Yp and UT1 are all about Earth. They may not be required for other bodies. Empty eopTab is acceptable.

            call checkout(subName)
            return
        end subroutine SyncEOP


        subroutine SetLeap(FName)
        !--------------------------------------------------------------------------------------------------------------------------
        !   Instruction:
        !       This subroutine reads the leap second series.
        !   Date:
        !       2014-05-11
        !   Inputs:
        !       FName       str         File name of Leapsecond series
        !   Outputs:
        !       No explicit output parameters, the Leap table is updated.
        !--------------------------------------------------------------------------------------------------------------------------
            use ErrorHandler
            use FileDef
            implicit none
            ! Declare input parameters:
            character(len=*),intent(in):: FName
            ! Declare local variables:
            integer:: yr,mo,day,ls
            type(File):: F
            character(len=4):: subName='slps'
            character(len=200):: image

            call checkin(subName)

            F%name = trim(FName)
            call OpenFile(F)
            if (errCode.gt.0) return
 98         read(F%unit,'(a200)',end=99) image
            read(image,1,err=97) yr,mo,day,ls
            call AddLeap(yr,mo,day,ls)
            if (errCode.gt.0) return
 97         goto 98
 99         continue

            if (nLeap.eq.0) then
                call RaiseError(4002,subName)
                return
            endif

            call CloseFile(F%unit)
            if (errCode.gt.0) return
            isLeapSet=.true.

 1          FORMAT(1x,i4,2(1x,i1),2x,i2)

            call checkout(subName)
            return
        end subroutine SetLeap

        subroutine SetEOP(FName)
        !--------------------------------------------------------------------------------------------------------------------------
        !   Instruction:
        !       This subroutine reads the EOP series.
        !   Date:
        !       2014-05-11
        !   Inputs:
        !       FName       str         File name of EOP series
        !   Outputs:
        !       No explicit output parameters, the EOP table is updated.
        !--------------------------------------------------------------------------------------------------------------------------
            use ErrorHandler
            use FileDef
            implicit none
            ! Declare input parameters:
            character(len=*),intent(in):: FName
            ! Declare local variables:
            integer:: t0,yr,mo,day
            double precision:: x0(4)
            double precision :: x1
            type(File):: F
            character(len=4):: subName='seop'
            character(len=200):: image,fmtstr
            character(len=1) :: c_tmp, c_tmp2
            integer :: ios = 0

            call checkin(subName)

            F%name = trim(FName)
            call OpenFile(F)
            if (errCode.gt.0) return
 98         read(F%unit,'(a200)',end=99) image
            if (trim(image).eq.'') goto 97
            if (index(image,'FORMAT(').ne.0) then
                fmtstr = image(index(image,'FORMAT(')+6:index(image,')',back=.true.))
                goto 97
            endif
            if (trim(fmtstr).eq.'') goto 97
!            read(image,1,err=97) t0,x0
            !read(image,trim(fmtstr),err=97) yr,mo,day,t0,x0
            read(image, "(I2,I2,I2,1x,I5,3x,1x,A1,1x,F9.6,F9.6,1x,F9.6,F9.6,2x,A1,F10.7)", iostat = ios) &
                yr, mo, day, t0, c_tmp, x0(1:4), c_tmp2, x1
            ! print*, ios, yr, mo, day, t0, c_tmp, x0(1:4), c_tmp2, x1
            if (ios == 0) then
                call AddEOP(dble(t0),x0(1),x0(3),x1)
                ! print*, nEOP
            endif
            if (errCode.gt.0) return
 97         goto 98
 99         continue

            if (nEOP.eq.0) then
                print*, "warning@setEOP, read empty EOPs."
                return
            else
                write(*, "('Read ', I5, ' EOPs')") nEOP
            endif

!            FORMAT(3(I4),I7,2(F11.6),2(F12.7),2(F11.6),2(F11.6),2(F11.7),2F12.6)
 1          FORMAT(12x,I7,2(F11.6),F12.7)

            call CloseFile(F%unit)
            if (errCode.gt.0) return
            isEOPSet=.true.

            call checkout(subName)
            return
        end subroutine SetEOP

        subroutine EOPLeapInit(dirEOPfile,dirLeapSecFile,isNoEOPLegal)
        !--------------------------------------------------------------------------------------------------------------------------
        !   Instruction:
        !       Initializing data for EOP and Leapsecond
        !   Date:
        !       2015-01-26
        !   Inputs:
        !       dirEOPfile(optional)             char*200        filename of EOP file
        !       dirLeapSecFile(optional)         char*200        filename of Leapsecond file
        !       isNoEOPLegal(optional)           logical         if the default value is accepted when EOP file is not avalible
        !   Outputs:
        !       No explicit output parameters, the directories of EOP file and leapsecond file and the logical varialble isDefautEOPLegal will be set.
        !--------------------------------------------------------------------------------------------------------------------------
            use ErrorHandler
            implicit none
            character(len=*),intent(in)::dirEOPfile,dirLeapSecFile
            logical,intent(in),optional:: isNoEOPLegal
            character(len=4):: subname='eint'
            call checkin(subname)
            ! Open and read EOP file
            call setEOP(dirEOPfile)
            if(ErrCode/=0)return
            ! Open and read leap second file
            call setLeap(dirLeapSecFile)
            if(ErrCode/=0)return
            if(present(isNoEOPLegal))isNoEOPaccept=isNoEOPLegal
            call checkout(subname)
        end subroutine EOPLeapInit

        subroutine NullifyPolarMotion
            implicit none
            Xp = 0d0
            Yp = 0d0
            Xp_rate = 0d0
            Yp_rate = 0d0
            return
        end subroutine NullifyPolarMotion

        subroutine ResetEOP(n,eop_list)
        !--------------------------------------------------------------------------------------------------------------------------
        !   Instruction:
        !       This subroutine copies the given EOP series to the shared memory or resets the memory (if n=0).
        !   Date:
        !       2018-11-29
        !   Inputs:
        !       n           integer         size of the list
        !       eop_list    double(n,4)     optional, eop list to be copied
        !   Outputs:
        !       No explicit output parameters, the EOP table is updated.
        !--------------------------------------------------------------------------------------------------------------------------
            use ErrorHandler
            implicit none
            ! Declare input parameters:
            integer,intent(in):: n
            double precision,intent(in),optional:: eop_list(n,4)
            ! Declare local variables:
            integer:: i
            character(len=4):: subName='reop'

            call checkin(subName)
            if (n.eq.0) then
                if (present(eop_list)) call RaiseWarning('The explicit EOP list will be ignored regardlessly when n=0.',subName)
                nEOP = 0
                isEOPSet = .false.
            else
                if (.not.present(eop_list)) then
                    !call RaiseError(0016,subName,__LINE__,__FILE__)
                    return
                endif
                do i = 1,n
                    call AddEOP(eop_list(i,1),eop_list(i,2),eop_list(i,3),eop_list(i,4))
                    if (errCode.gt.0) return
                enddo
                isEOPSet=.true.
            endif

            call checkout(subName)
            return
        end subroutine ResetEOP

        subroutine ResetLeap(n,leap_list)
        !--------------------------------------------------------------------------------------------------------------------------
        !   Instruction:
        !       This subroutine copies the given leap second series to the shared memory or resets the memory (if n=0).
        !   Date:
        !       2018-11-29
        !   Inputs:
        !       n           integer         size of the list
        !       leap_list   double(n,5)     optional,leap second list
        !   Outputs:
        !       No explicit output parameters, the Leap table is updated.
        !--------------------------------------------------------------------------------------------------------------------------
            use ErrorHandler
            implicit none
            ! Declare input parameters:
            integer,intent(in):: n
            double precision,intent(in),optional:: leap_list(n,5)
            ! Declare local variables:
            integer:: i
            character(len=4):: subName='rlps'

            call checkin(subName)

            if (n.eq.0) then
                if (present(leap_list)) call RaiseWarning('The explicit Leap list will be ignored regardlessly when n=0.',subName)
                nLeap = 0
                isLeapSet = .false.
            else
                if (.not.present(leap_list)) then
                    !call RaiseError(0016,subName,__LINE__,__FILE__)
                    return
                endif
                do i = 1,n
                call AddLeap(idint(leap_list(i,2)),idint(leap_list(i,3)),idint(leap_list(i,4)),idint(leap_list(i,5)))
                    if (errCode.gt.0) return
                enddo
                isLeapSet=.true.
            endif

            call checkout(subName)
            return
        end subroutine ResetLeap
        
        double precision function EncodeMJD(yr,mo,day,hr,min,sec)
        !------------------------------------------------------------
        !   Instruction:
        !       This function computes MJD of given date.
        !   Date:
        !       2013-09-25
        !   Inputs:
        !       yr,mo,day   i           year,month,day of the epoch
        !       hr,min      i           hour,minute of the epoch
        !       sec         d           second of the epoch
        !   Outputs:
        !       EncodeMJD   d           corresponding MJD
        !------------------------------------------------------------
        implicit none
        ! Declare input parameters:
        integer,intent(in):: yr,mo,day,hr,min
        double precision,intent(in):: sec
        ! Declare local variables:
        integer:: itimes(5),jd
        double precision:: mjd

        itimes = (/yr,mo,day,hr,min/)

        JD = ITIMES(3) - 32075 + &
             1461 * (ITIMES(1) + 4800 + (ITIMES(2)-14)/12) /4  + &
             367 * (ITIMES(2) - 2 - (ITIMES(2)-14)/12*12) /12  - &
             3 * ((ITIMES(1) + 4900 + (ITIMES(2)-14)/12)/100) /4

        MJD = DBLE(JD) - 0.5D0  -2400000.5D0 + &
              DBLE(ITIMES(4)) / 24.0D0  +      &
              DBLE(ITIMES(5)) / 1440.0D0  +    &
              DBLE(sec) / 86400.0D0

        EncodeMJD = MJD
        return
        end function EncodeMJD
        
        double precision function LinInt(x1,y1,x2,y2,x)
        !--------------------------------------------------------------------------------------------------------------------------
        !   Instruction:
        !       This function interpolates the given points with a linear function.
        !   Date:
        !       2013-09-25
        !   Inputs:
        !       x1,y1       d           first point
        !       x2,y2       d           second point
        !       x           d           x of the point whose y is to be evaluated
        !   Outputs:
        !       LinInt      d           evaluated y of the point
        !--------------------------------------------------------------------------------------------------------------------------
        use ErrorHandler
        implicit none
        double precision,parameter:: EPS = 1d-10
        ! Declare input parameters:
        double precision,intent(in):: x1,y1,x2,y2,x
        ! Declare local variables:
        character(len=4):: subName='dlin'
        call checkin(subName)
        if (dmax1(dabs(1d0-x2/x1),dabs(x1),dabs(x2)).lt.EPS) then
            call RaiseError(1008,subName)
            return
        endif
        if (x.lt.dmin1(x1,x2).or.x.gt.dmax1(x1,x2)) then
            call RaiseWarning('Given x to interpolate is not between the nodes. Going to extrapolate.',subName)
        endif
        LinInt = y1+(y2-y1)*(x-x1)/(x2-x1)
        call checkout(subName)
        return
        end function LinInt

end module EOPParam
