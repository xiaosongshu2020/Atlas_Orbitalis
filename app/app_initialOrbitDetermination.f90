
    program main

    use module_const
    use module_time
    use module_station
    use module_satellite
    use module_access
    use module_IODor
    implicit none
    
    ! �����в���
    character(len=256) :: arg
    real*8 :: lon_deg, lat_deg, altitude_m
    character(len=10) :: obs_type_str
    character(len=256) :: obs_filename
    
    ! �������
    type(type_station) :: station
    type(type_satellite) :: satellite
    type(type_access) :: access
    type(type_IODor) :: iodor
    type(type_t) :: t0
    integer :: obsType
    integer :: nargs, ios
    logical :: file_exists
    real*8 :: initial_kpl(6)
    
    ! ��ȡ�����в�������
    nargs = command_argument_count()
    
    ! ����������
    if (nargs /= 5) then
        call print_usage()
        stop
    endif
    
    ! ��ȡ�����в���
    ! 1. ��վ����
    call get_command_argument(1, arg)
    read(arg, *, iostat=ios) lon_deg
    if (ios /= 0) then
        print*, "Error: Invalid longitude value"
        call print_usage()
        stop
    endif
    
    ! 2. ��վγ��
    call get_command_argument(2, arg)
    read(arg, *, iostat=ios) lat_deg
    if (ios /= 0) then
        print*, "Error: Invalid latitude value"
        call print_usage()
        stop
    endif
    
    ! 3. ��վ����
    call get_command_argument(3, arg)
    read(arg, *, iostat=ios) altitude_m
    if (ios /= 0) then
        print*, "Error: Invalid altitude value"
        call print_usage()
        stop
    endif
    
    ! 4. �۲���������
    call get_command_argument(4, obs_type_str)
    call to_upper(obs_type_str)
    
    ! 5. �۲������ļ���
    call get_command_argument(5, obs_filename)
    
    ! ����ļ��Ƿ����
    inquire(file=trim(obs_filename), exist=file_exists)
    if (.not. file_exists) then
        print*, "Error: Observation file not found: ", trim(obs_filename)
        stop
    endif
    
    ! ת���۲������ַ�������������
    select case(trim(obs_type_str))
    case('AH')
        obsType = OBS_TYPE_AZIMUTH_ELEVATION
        print*, "Observation type: Azimuth/Elevation"
    case('AD')
        obsType = OBS_TYPE_RA_DEC
        print*, "Observation type: Right Ascension/Declination"
    case('RR')
        obsType = OBS_TYPE_RANGERATE
        print*, "Observation type: Range Rate"
        print*, "Warning: range-rate data is not supported by the General Laplace method"
    case default
        print*, "Error: Unknown observation type: ", trim(obs_type_str)
        print*, "Supported types: AH, AD, RR"
        stop
    end select
    
    !! ��ӡ�������
    !print*, "=============================================================="
    !print*, "             INITIAL ORBIT DETERMINATION PROGRAM              "
    !print*, "=============================================================="
    !print*, ""
    !print*, "Input Parameters:"
    !print*, "-----------------"
    !write(*, '(A,F10.5,A)') " Station Longitude:  ", lon_deg, " deg"
    !write(*, '(A,F10.5,A)') " Station Latitude:   ", lat_deg, " deg"
    !write(*, '(A,F10.2,A)') " Station Altitude:   ", altitude_m, " m"
    !write(*, '(A,A)')       " Observation Type:   ", trim(obs_type_str)
    !write(*, '(A,A)')       " Data File:          ", trim(obs_filename)
    !print*, ""
    
    ! ��ʼ����վ
    call station%init(1, "User Station", lon_deg, lat_deg, altitude_m)
    
    ! ��ʼ�����ǣ�ʹ��Ĭ�ϳ�ֵ������IOD������ǣ�
    initial_kpl = (/ 1.2d0, 0.001d0, 45d0*deg2rad, 0d0, 0d0, 0d0 /)
    t0%val = 60000d0  ! Ĭ����Ԫ�����ӹ۲������и���
    t0%day = int(t0%val)
    t0%sec = (t0%val - t0%day) * 86400d0
    
    satellite%number = 1
    satellite%name = "Target Satellite"
    call satellite%init(t0, initial_kpl, &
        number=1, &
        name="Target Satellite", &
        eleType=STATE_TYPE_KEPLERIAN, &
        propagator="two-body")
    
    ! ��ʼ���۲����ݷ��ʶ���
    call access%init(number=1, name="IOD Observations", &
        &   fromType=OBJECT_TYPE_STATION, fromNum=station%number, &
        &   toType=OBJECT_TYPE_SATELLITE, toNum=satellite%number, &
        &   obsType=obsType, minEle=5d0)   ! minimum elevation 5 deg
    
    ! ��ȡ�۲�����
    !print*, "Reading observation data..."
    call access%readData(obs_filename)
    
    if (access%nData <= 0) then
        print*, "Error: No valid observation data found"
        stop
    endif
    
    !print*, "Successfully read ", access%nData, " observations"
    !print*, ""
    
    ! �����ʼ��Ԫ��ʹ�ù۲����ݵĳ�ʼ��
    t0 = access%timeTable(1)
    
    ! ��ʼ��IODor����
    call iodor%init(t0%val, station, satellite, access)
    
    !! ִ�г���ȷ��
    !print*, "Performing initial orbit determination..."
    !print*, "Using General Laplace Method"
    !print*, ""
    
    call iodor%solve()
    
    ! ��ʾ���
    call iodor%show()
    
    ! ���������ļ�
    call iodor%saveIODor("iod_results.txt")
    
    !print*, ""
    !print*, "Program completed successfully!"
    
contains

    !===========================================================================
    ! ��ӡʹ��˵��
    !===========================================================================
    subroutine print_usage()
        print*, "Usage: main_iod <longitude> <latitude> <altitude> <obs_type> <data_file>"
        print*, ""
        print*, "Parameters:"
        print*, "  longitude  : Station longitude in degrees (-180 to 180)"
        print*, "  latitude   : Station latitude in degrees (-90 to 90)"
        print*, "  altitude   : Station altitude in meters"
        print*, "  obs_type   : Observation type (AH, AD, or RR)"
        print*, "               AH - Azimuth/Elevation"
        print*, "               AD - Right Ascension/Declination"
        print*, "               RR - Range rate (not supported by General Laplace)"
        print*, "  data_file  : Path to observation data file"
        print*, ""
        print*, "Example:"
        print*, "  ./main_iod 116.385 39.905 50.0 AH observations.dat"
    end subroutine print_usage
    
    !===========================================================================
    ! ת���ַ���Ϊ��д
    !===========================================================================
    subroutine to_upper(str)
        character(len=*), intent(inout) :: str
        integer :: i, ich
        
        do i = 1, len_trim(str)
            ich = ichar(str(i:i))
            if (ich >= ichar('a') .and. ich <= ichar('z')) then
                str(i:i) = char(ich - 32)
            endif
        enddo
    end subroutine to_upper

end program main