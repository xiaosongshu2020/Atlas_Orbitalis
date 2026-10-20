
    module module_IODor
    use module_const
    use module_planets
    use module_satellite
    use module_station
    use module_access

    ! 默认所有子程序为private
    !private

    type :: type_IODor
        ! 这里可以定义一些类型相关的变量
        double precision :: t0
        type(type_station) :: station
        type(type_satellite) :: satellite
        type(type_access) :: access
        double precision :: max_data_len = 800d0 ! 最大数据长度, sec
        contains
        procedure :: init => init_IODor
        procedure :: solve => solve_IODor
        procedure :: show => show_IODor
        procedure :: saveIODor => save_IODor
    end type type_IODor
    ! Public interface for the module
    contains

    subroutine init_IODor(this, t0, station, satellite, access)
        implicit none
        class(type_IODor), intent(inout) :: this
        double precision, intent(in) :: t0
        type(type_station), intent(in) :: station
        type(type_satellite), intent(in) :: satellite
        type(type_access), intent(in) :: access

        this%t0 = t0
        this%station = station
        this%satellite = satellite
        this%access = access
    end subroutine init_IODor

    subroutine solve_IODor(this)
        use module_time
        use module_space
        implicit none
        class(type_IODor), intent(inout) :: this
        integer :: nD, i, j, nValid
        double precision, allocatable :: t(:), alpha(:), delta(:), X(:), Y(:), Z(:)
        double precision, allocatable :: A(:), h(:), alpha_valid(:), delta_valid(:)
        double precision, allocatable :: t_valid(:), X_valid(:), Y_valid(:), Z_valid(:)
        double precision :: cts0(6), kpl0(6)
        double precision :: rStation_ECF(6), rStation_ECI(6)
        double precision :: lon_rad, lat_rad
        double precision :: dr_horizon(3), dr_ECI(3), dr_ECF(3)
        double precision :: sinA, cosA, sinh, cosh
        double precision :: utc, dt_total
        type(type_t) :: t_epoch
        type(type_timeSeries) :: ts
        logical :: valid_data

        ! 检查是否有观测数据
        if (this%access%nData <= 0) then
            print*, "Error @solve_IODor: No observation data available"
            return
        endif

        ! 获取数据点数量
        nD = this%access%nData
        
        ! 检查数据时间跨度
        dt_total = timeMinus(this%access%timeTable(nD), this%access%timeTable(1))
        if (dt_total > this%max_data_len/86400d0) then
            print*, "Warning @solve_IODor: Data span exceeds maximum allowed length"
            print*, "Data span: ", dt_total*86400d0, " seconds"
            print*, "Maximum allowed: ", this%max_data_len, " seconds"
            ! 截取数据到最大长度
            do i = 1, nD
                if (timeMinus(this%access%timeTable(i), this%access%timeTable(1)) > this%max_data_len/86400d0) then
                    nD = i - 1
                    exit
                endif
            enddo
        endif

        ! 分配内存
        allocate(t(nD), alpha(nD), delta(nD), X(nD), Y(nD), Z(nD))
        
        ! 获取测站地心固连坐标
        rStation_ECF = (/ this%station%position, 0d0, 0d0, 0d0 /)
        lon_rad = this%station%longitude
        lat_rad = this%station%latitude
        ts%time = this%access%timeTable(1:nD)
        call this%station%generateTrajectory(ts)
        do i = 1, nD
            t(i) = this%access%timeTable(i)%val
            X(i) = this%station%trajectory(i, 1)
            Y(i) = this%station%trajectory(i, 2)
            Z(i) = this%station%trajectory(i, 3)
        enddo

        ! 处理不同类型的观测数据
        select case(this%access%obsType)
        
        case(OBS_TYPE_RA_DEC)
            ! 直接从access中获取赤经赤纬数据
            nValid = 0
            do i = 1, nD
                
                ! 获取赤经赤纬（已经是弧度）
                alpha(i) = this%access%dataTable(i, 1)
                delta(i) = this%access%dataTable(i, 2)

                nValid = nValid + 1
            enddo
            
        case(OBS_TYPE_AZIMUTH_ELEVATION)
            
            if (allocated(A)) deallocate(A)
            if (allocated(h)) deallocate(h)
            allocate(A(nD))
            allocate(h(nD))
            
            ! 从方位角和高度角转换到赤经赤纬
            nValid = 0
            do i = 1, nD
                ! 获取时间
                
                utc = t(i)
                
                ! 获取方位角和高度角（已经是弧度）
                A(i) = this%access%dataTable(i, 1)
                h(i) = this%access%dataTable(i, 2)

                ! 转换到地平坐标系的方向向量
                sinh = dsin(h(i))
                cosh = dcos(h(i))
                sinA = dsin(A(i))
                cosA = dcos(A(i))
                
                dr_horizon(1) = cosh * cosA  ! 北向
                dr_horizon(2) = cosh * sinA  ! 东向
                dr_horizon(3) = sinh         ! 天向
                
                ! 从地平坐标系转换到地心固连坐标系
                call horizon2ECF_position(lon_rad, lat_rad, dr_horizon, dr_ECF)

                ! 从地心固连坐标系转换到J2000坐标系
                call ECF2ECI_IAU1980_positionVelocity(utc, dr_ECF, &
                    (/0d0, 0d0, 0d0/), dr_ECI, (/0d0, 0d0, 0d0/))
                
                ! ! 计算赤经赤纬
                delta(i) = dasin(dr_ECI(3))
                call get_angle(dr_ECI(2), dr_ECI(1), alpha(i))

                nValid = nValid + 1
            enddo
            
        case(OBS_TYPE_RANGERATE)
            print*, "Error @solve_IODor: Range-rate data not supported by General Laplace"
            deallocate(t, alpha, delta, X, Y, Z)
            return
            
        case default
            print*, "Error @solve_IODor: Unsupported observation type"
            deallocate(t, alpha, delta, X, Y, Z)
            return
        end select

        ! 检查是否有足够的有效数据
        if (nValid < 3) then
            print*, "Error @solve_IODor: Insufficient valid observations"
            print*, "Valid observations: ", nValid
            print*, "Minimum required: 3"
            deallocate(t, alpha, delta, X, Y, Z)
            return
        endif

        ! 如果有无效数据，重新整理数组
        if (nValid < nD) then
            allocate(t_valid(nValid), alpha_valid(nValid), delta_valid(nValid))
            allocate(X_valid(nValid), Y_valid(nValid), Z_valid(nValid))
            
            j = 0
            do i = 1, nD
                if (this%access%obsType == OBS_TYPE_AZIMUTH_ELEVATION) then
                    if (h(i) < this%access%minEle * deg2rad) cycle
                endif
                j = j + 1
                t_valid(j) = t(i)
                alpha_valid(j) = alpha(i)
                delta_valid(j) = delta(i)
                X_valid(j) = X(i)
                Y_valid(j) = Y(i)
                Z_valid(j) = Z(i)
            enddo
            
            ! 使用有效数据进行初轨确定
            call initialOrbitDetermination_generalLaplaceMethod(nValid, this%t0, &
                t_valid, alpha_valid, delta_valid, X_valid, Y_valid, Z_valid, cts0)
            
            deallocate(t_valid, alpha_valid, delta_valid, X_valid, Y_valid, Z_valid)
        else
            ! 使用所有数据进行初轨确定
            call initialOrbitDetermination_generalLaplaceMethod(nD, this%t0, &
                t, alpha, delta, X, Y, Z, cts0)
        endif

        ! 将结果存储到卫星对象中
        t_epoch%val = this%t0
        t_epoch%day = int(this%t0)
        t_epoch%sec = (this%t0 - int(this%t0)) * 86400d0
        
        ! 转换到开普勒轨道根数
        call cts2kpl(cts0, kpl0)
        
        call kpl_normal2SI(kpl0, kpl0)
        
        ! 更新卫星状态
        call this%satellite%init(t_epoch, kpl0, &
            number=this%satellite%number, &
            name=trim(this%satellite%name), &
            eleType=STATE_TYPE_KEPLERIAN)
        
        !! 输出结果
        !print*, "Initial orbit determination completed successfully!"
        !print*, "Epoch (MJD): ", this%t0
        !print*, "Position (km): ", cts0(1:3) * length_unit
        !print*, "Velocity (km/s): ", cts0(4:6) * length_unit / time_unit
        !print*, "Semi-major axis (km): ", kpl0(1) * length_unit
        !print*, "Eccentricity: ", kpl0(2)
        !print*, "Inclination (deg): ", kpl0(3) * rad2deg
        !print*, "RAAN (deg): ", kpl0(4) * rad2deg
        !print*, "Argument of perigee (deg): ", kpl0(5) * rad2deg
        !print*, "Mean anomaly (deg): ", kpl0(6) * rad2deg
        
        ! 释放内存
        deallocate(t, alpha, delta, X, Y, Z)
        if (allocated(A)) deallocate(A)
        if (allocated(h)) deallocate(h)
        
    end subroutine solve_IODor

    subroutine show_IODor(this)
        use module_time
        implicit none
        class(type_IODor), intent(in) :: this
        type(type_t) :: epoch
        double precision :: kpl(6), cts(6), nsg1(6)
        character(len=30) :: dateStr
        character(len=20) :: obsTypeStr
        integer :: i
        
        print*, "=============================================================="
        print*, "            INITIAL ORBIT DETERMINATION RESULTS               "
        print*, "=============================================================="
        
        ! 时间信息
        epoch = this%satellite%states%t
        !call epoch%toDateStr(dateStr)
        dateStr = epoch%toDate()
        print*, ""
        print*, "EPOCH INFORMATION:"
        print*, "------------------"
        write(*, '(A,F12.6)') " MJD (UTC):        ", epoch%val
        write(*, '(A,A)')     " Date/Time:        ", trim(dateStr)
        
        ! 测站信息
        print*, ""
        print*, "STATION INFORMATION:"
        print*, "--------------------"
        write(*, '(A,I6)')    " Station Number:   ", this%station%number
        write(*, '(A,A)')     " Station Name:     ", trim(this%station%name)
        write(*, '(A,F10.5,A)') " Longitude:        ", this%station%longitude*rad2deg, " deg"
        write(*, '(A,F10.5,A)') " Latitude:         ", this%station%latitude*rad2deg, " deg"
        write(*, '(A,F10.2,A)') " Altitude:         ", this%station%altitude / u_meter, " m"
        write(*, '(A,3F12.3,A)') " ECEF Position:    ", this%station%position*EARTH_RADIUS_EQUATORIAL, " km"
        
        ! 观测数据信息
        print*, ""
        print*, "OBSERVATION INFORMATION:"
        print*, "------------------------"
        write(*, '(A,I6)')    " Access Number:    ", this%access%number
        write(*, '(A,A)')     " Access Name:      ", trim(this%access%name)
        
        ! 观测类型
        select case(this%access%obsType)
        case(OBS_TYPE_RA_DEC)
            obsTypeStr = "RA/Dec"
        case(OBS_TYPE_AZIMUTH_ELEVATION)
            obsTypeStr = "Azimuth/Elevation"
        case(OBS_TYPE_RANGING)
            obsTypeStr = "Range Only"
        case(OBS_TYPE_RANGERATE)
            obsTypeStr = "Range Rate"
        case(OBS_TYPE_POSITION)
            obsTypeStr = "Position Only"
        case(OBS_TYPE_POSITION_VELOCITY)
            obsTypeStr = "Position/Velocity"
        case default
            obsTypeStr = "Unknown"
        end select
        
        write(*, '(A,A)')     " Observation Type: ", trim(obsTypeStr)
        write(*, '(A,I6)')    " Number of Data:   ", this%access%nData
        
        if (this%access%nData > 0) then
            write(*, '(A,F12.6)') " First Obs (MJD):  ", this%access%timeTable(1)%val
            write(*, '(A,F12.6)') " Last Obs (MJD):   ", this%access%timeTable(this%access%nData)%val
            write(*, '(A,F10.2,A)') " Data Span:        ", &
                timeMinus(this%access%timeTable(this%access%nData), this%access%timeTable(1))*86400d0, " sec"
        endif
        
        ! 轨道确定结果
        print*, ""
        print*, "ORBIT DETERMINATION RESULTS:"
        print*, "----------------------------"
        
        ! 根据卫星的轨道要素类型显示相应的结果
        if (this%satellite%eleType == STATE_TYPE_KEPLERIAN) then
            kpl = this%satellite%states%elements
            call kpl2cts(kpl, cts)
            
            print*, ""
            print*, "Keplerian Elements:"
            write(*, '(A,F12.3,A)')  "  Semi-major axis:        ", kpl(1)*length_unit, " km"
            write(*, '(A,F12.10)')   "  Eccentricity:           ", kpl(2)
            write(*, '(A,F12.6,A)')  "  Inclination:            ", kpl(3)*rad2deg, " deg"
            write(*, '(A,F12.6,A)')  "  RAAN:                   ", kpl(4)*rad2deg, " deg"
            write(*, '(A,F12.6,A)')  "  Argument of Perigee:    ", kpl(5)*rad2deg, " deg"
            write(*, '(A,F12.6,A)')  "  Mean Anomaly:           ", kpl(6)*rad2deg, " deg"
            
            print*, ""
            print*, "Cartesian State Vector (J2000):"
            write(*, '(A,3F15.6,A)') "  Position:         ", cts(1:3)*length_unit, " km"
            write(*, '(A,3F12.6,A)') "  Velocity:         ", cts(4:6)*length_unit/time_unit, " km/s"
            
        elseif (this%satellite%eleType == STATE_TYPE_CARTESIAN) then
            cts = this%satellite%states%elements
            call cts2kpl(cts, kpl)
            
            print*, ""
            print*, "Cartesian State Vector (J2000):"
            write(*, '(A,3F15.6,A)') "  Position:         ", cts(1:3)*length_unit, " km"
            write(*, '(A,3F12.6,A)') "  Velocity:         ", cts(4:6)*length_unit/time_unit, " km/s"
            
            print*, ""
            print*, "Keplerian Elements:"
            write(*, '(A,F12.3,A)')  "  Semi-major axis:        ", kpl(1)*length_unit, " km"
            write(*, '(A,F12.10)')   "  Eccentricity:           ", kpl(2)
            write(*, '(A,F12.6,A)')  "  Inclination:            ", kpl(3)*rad2deg, " deg"
            write(*, '(A,F12.6,A)')  "  RAAN:                   ", kpl(4)*rad2deg, " deg"
            write(*, '(A,F12.6,A)')  "  Argument of Perigee:    ", kpl(5)*rad2deg, " deg"
            write(*, '(A,F12.6,A)')  "  Mean Anomaly:           ", kpl(6)*rad2deg, " deg"
            
        elseif (this%satellite%eleType == STATE_TYPE_NONSINGULAR_1) then
            nsg1 = this%satellite%states%elements
            call nsg12kpl(nsg1, kpl)
            call nsg12cts(nsg1, cts)
            
            print*, ""
            print*, "Non-Singular Elements (Type 1):"
            write(*, '(A,F12.3,A)')  "  Semi-major axis:        ", nsg1(1)*length_unit, " km"
            write(*, '(A,F12.6,A)')  "  h = e*sin(omega+Omega): ", nsg1(2)*rad2deg, " deg"
            write(*, '(A,F12.6,A)')  "  k = e*cos(omega+Omega): ", nsg1(3)*rad2deg, " deg"
            write(*, '(A,F12.10)')   "  p = tan(i/2)*sin(Omega):", nsg1(4)
            write(*, '(A,F12.10)')   "  q = tan(i/2)*cos(Omega):", nsg1(5)
            write(*, '(A,F12.6,A)')  "  Mean Longitude:         ", nsg1(6)*rad2deg, " deg"
            
            print*, ""
            print*, "Equivalent Keplerian Elements:"
            write(*, '(A,F12.3,A)')  "  Semi-major axis:        ", kpl(1)*length_unit, " km"
            write(*, '(A,F12.10)')   "  Eccentricity:           ", kpl(2)
            write(*, '(A,F12.6,A)')  "  Inclination:            ", kpl(3)*rad2deg, " deg"
            write(*, '(A,F12.6,A)')  "  RAAN:                   ", kpl(4)*rad2deg, " deg"
            write(*, '(A,F12.6,A)')  "  Argument of Perigee:    ", kpl(5)*rad2deg, " deg"
            write(*, '(A,F12.6,A)')  "  Mean Anomaly:           ", kpl(6)*rad2deg, " deg"
        endif
        
        ! 轨道特征
        print*, ""
        print*, "Orbital Characteristics:"
        write(*, '(A,F12.3,A)')  "  Perigee Altitude:       ", (kpl(1)*(1d0-kpl(2))-1d0)*length_unit, " km"
        write(*, '(A,F12.3,A)')  "  Apogee Altitude:        ", (kpl(1)*(1d0+kpl(2))-1d0)*length_unit, " km"
        write(*, '(A,F12.6,A)')  "  Orbital Period:         ", 2d0*pi*sqrt(kpl(1)**3)*time_unit/60d0, " min"
        
        print*, ""
        print*, "=============================================================="
        
    end subroutine show_IODor
    
    
    subroutine save_IODor(this, filename)
        use module_time
        implicit none
        class(type_IODor), intent(in) :: this
        character(len=*), intent(in) :: filename
        type(type_t) :: epoch
        double precision :: kpl(6), cts(6), nsg1(6)
        character(len=30) :: dateStr
        character(len=20) :: obsTypeStr
        integer :: unit_num, ios
        character(len=8) :: date
        character(len=10) :: time
        
        ! 打开文件
        unit_num = 999
        open(unit=unit_num, file=trim(filename), status='replace', iostat=ios)
        if (ios /= 0) then
            print*, "Error: Cannot open file ", trim(filename)
            return
        endif
        
        ! 写入文件头
        write(unit_num, *) "=============================================================="
        write(unit_num, *) "            INITIAL ORBIT DETERMINATION RESULTS               "
        write(unit_num, *) "=============================================================="
        
        ! 写入计算时间
        call date_and_time(date, time)
        write(unit_num, '(A,A,A,A)') " Generated on: ", date(1:4)//"-"//date(5:6)//"-"//date(7:8), &
                                     " at ", time(1:2)//":"//time(3:4)//":"//time(5:6)
        write(unit_num, *) ""
        
        ! 时间信息
        epoch = this%satellite%states%t
        !call epoch%toDateStr(dateStr)
        dateStr = epoch%toDate()
        write(unit_num, *) "EPOCH INFORMATION:"
        write(unit_num, *) "------------------"
        write(unit_num, '(A,F12.6)') " MJD (UTC):        ", epoch%val
        write(unit_num, '(A,A)')     " Date/Time:        ", trim(dateStr)
        
        ! 测站信息
        write(unit_num, *) ""
        write(unit_num, *) "STATION INFORMATION:"
        write(unit_num, *) "--------------------"
        write(unit_num, '(A,I6)')    " Station Number:   ", this%station%number
        write(unit_num, '(A,A)')     " Station Name:     ", trim(this%station%name)
        write(unit_num, '(A,F10.5,A)') " Longitude:        ", this%station%longitude*rad2deg, " deg"
        write(unit_num, '(A,F10.5,A)') " Latitude:         ", this%station%latitude*rad2deg, " deg"
        write(unit_num, '(A,F10.2,A)') " Altitude:         ", this%station%altitude / u_meter, " m"
        write(unit_num, '(A,3F12.3,A)') " ECEF Position:    ", this%station%position*EARTH_RADIUS_EQUATORIAL, " km"
        
        ! 观测数据信息
        write(unit_num, *) ""
        write(unit_num, *) "OBSERVATION INFORMATION:"
        write(unit_num, *) "------------------------"
        write(unit_num, '(A,I6)')    " Access Number:    ", this%access%number
        write(unit_num, '(A,A)')     " Access Name:      ", trim(this%access%name)
        
        ! 观测类型
        select case(this%access%obsType)
        case(OBS_TYPE_RA_DEC)
            obsTypeStr = "RA/Dec"
        case(OBS_TYPE_AZIMUTH_ELEVATION)
            obsTypeStr = "Azimuth/Elevation"
        case(OBS_TYPE_RANGING)
            obsTypeStr = "Range Only"
        case(OBS_TYPE_RANGERATE)
            obsTypeStr = "Range Rate"
        case(OBS_TYPE_POSITION)
            obsTypeStr = "Position Only"
        case(OBS_TYPE_POSITION_VELOCITY)
            obsTypeStr = "Position/Velocity"
        case default
            obsTypeStr = "Unknown"
        end select
        
        write(unit_num, '(A,A)')     " Observation Type: ", trim(obsTypeStr)
        write(unit_num, '(A,I6)')    " Number of Data:   ", this%access%nData
        
        if (this%access%nData > 0) then
            write(unit_num, '(A,F12.6)') " First Obs (MJD):  ", this%access%timeTable(1)%val
            write(unit_num, '(A,F12.6)') " Last Obs (MJD):   ", this%access%timeTable(this%access%nData)%val
            write(unit_num, '(A,F10.2,A)') " Data Span:        ", &
                timeMinus(this%access%timeTable(this%access%nData), this%access%timeTable(1))*86400d0, " sec"
        endif
        
        ! 轨道确定结果
        write(unit_num, *) ""
        write(unit_num, *) "ORBIT DETERMINATION RESULTS:"
        write(unit_num, *) "----------------------------"
        
        ! 根据卫星的轨道要素类型显示相应的结果
        if (this%satellite%eleType == STATE_TYPE_KEPLERIAN) then
            kpl = this%satellite%states%elements
            call kpl2cts(kpl, cts)
            
            write(unit_num, *) ""
            write(unit_num, *) "Keplerian Elements:"
            write(unit_num, '(A,F12.3,A)')  "  Semi-major axis:        ", kpl(1)*length_unit, " km"
            write(unit_num, '(A,F12.10)')   "  Eccentricity:           ", kpl(2)
            write(unit_num, '(A,F12.6,A)')  "  Inclination:            ", kpl(3)*rad2deg, " deg"
            write(unit_num, '(A,F12.6,A)')  "  RAAN:                   ", kpl(4)*rad2deg, " deg"
            write(unit_num, '(A,F12.6,A)')  "  Argument of Perigee:    ", kpl(5)*rad2deg, " deg"
            write(unit_num, '(A,F12.6,A)')  "  Mean Anomaly:           ", kpl(6)*rad2deg, " deg"
            
            write(unit_num, *) ""
            write(unit_num, *) "Cartesian State Vector (J2000):"
            write(unit_num, '(A,3F15.6,A)') "  Position:         ", cts(1:3)*length_unit, " km"
            write(unit_num, '(A,3F12.6,A)') "  Velocity:         ", cts(4:6)*length_unit/time_unit, " km/s"
            
        elseif (this%satellite%eleType == STATE_TYPE_CARTESIAN) then
            cts = this%satellite%states%elements
            call cts2kpl(cts, kpl)
            
            write(unit_num, *) ""
            write(unit_num, *) "Cartesian State Vector (J2000):"
            write(unit_num, '(A,3F15.6,A)') "  Position:         ", cts(1:3)*length_unit, " km"
            write(unit_num, '(A,3F12.6,A)') "  Velocity:         ", cts(4:6)*length_unit/time_unit, " km/s"
            
            write(unit_num, *) ""
            write(unit_num, *) "Keplerian Elements:"
            write(unit_num, '(A,F12.3,A)')  "  Semi-major axis:        ", kpl(1)*length_unit, " km"
            write(unit_num, '(A,F12.10)')   "  Eccentricity:           ", kpl(2)
            write(unit_num, '(A,F12.6,A)')  "  Inclination:            ", kpl(3)*rad2deg, " deg"
            write(unit_num, '(A,F12.6,A)')  "  RAAN:                   ", kpl(4)*rad2deg, " deg"
            write(unit_num, '(A,F12.6,A)')  "  Argument of Perigee:    ", kpl(5)*rad2deg, " deg"
            write(unit_num, '(A,F12.6,A)')  "  Mean Anomaly:           ", kpl(6)*rad2deg, " deg"
            
        elseif (this%satellite%eleType == STATE_TYPE_NONSINGULAR_1) then
            nsg1 = this%satellite%states%elements
            call nsg12kpl(nsg1, kpl)
            call nsg12cts(nsg1, cts)
            
            write(unit_num, *) ""
            write(unit_num, *) "Non-Singular Elements (Type 1):"
            write(unit_num, '(A,F12.3,A)')  "  Semi-major axis:        ", nsg1(1)*length_unit, " km"
            write(unit_num, '(A,F12.6,A)')  "  h = e*sin(omega+Omega): ", nsg1(2)*rad2deg, " deg"
            write(unit_num, '(A,F12.6,A)')  "  k = e*cos(omega+Omega): ", nsg1(3)*rad2deg, " deg"
            write(unit_num, '(A,F12.10)')   "  p = tan(i/2)*sin(Omega):", nsg1(4)
            write(unit_num, '(A,F12.10)')   "  q = tan(i/2)*cos(Omega):", nsg1(5)
            write(unit_num, '(A,F12.6,A)')  "  Mean Longitude:         ", nsg1(6)*rad2deg, " deg"
            
            write(unit_num, *) ""
            write(unit_num, *) "Equivalent Keplerian Elements:"
            write(unit_num, '(A,F12.3,A)')  "  Semi-major axis:        ", kpl(1)*length_unit, " km"
            write(unit_num, '(A,F12.10)')   "  Eccentricity:           ", kpl(2)
            write(unit_num, '(A,F12.6,A)')  "  Inclination:            ", kpl(3)*rad2deg, " deg"
            write(unit_num, '(A,F12.6,A)')  "  RAAN:                   ", kpl(4)*rad2deg, " deg"
            write(unit_num, '(A,F12.6,A)')  "  Argument of Perigee:    ", kpl(5)*rad2deg, " deg"
            write(unit_num, '(A,F12.6,A)')  "  Mean Anomaly:           ", kpl(6)*rad2deg, " deg"
        endif
        
        ! 轨道特征
        write(unit_num, *) ""
        write(unit_num, *) "Orbital Characteristics:"
        write(unit_num, '(A,F12.3,A)')  "  Perigee Altitude:       ", (kpl(1)*(1d0-kpl(2))-1d0)*length_unit, " km"
        write(unit_num, '(A,F12.3,A)')  "  Apogee Altitude:        ", (kpl(1)*(1d0+kpl(2))-1d0)*length_unit, " km"
        write(unit_num, '(A,F12.6,A)')  "  Orbital Period:         ", 2d0*pi*sqrt(kpl(1)**3)*time_unit/60d0, " min"
        
        ! 额外的数据表格（适合后续处理）
        write(unit_num, *) ""
        write(unit_num, *) "=============================================================="
        write(unit_num, *) "COMPACT DATA FORMAT (for post-processing):"
        write(unit_num, *) "=============================================================="
        write(unit_num, '(A)') "# Epoch(MJD), a(km), e, i(deg), RAAN(deg), omega(deg), M(deg)"
        write(unit_num, '(F12.6,6(1X,F15.8))') epoch%val, &
            kpl(1)*length_unit, kpl(2), &
            kpl(3)*rad2deg, kpl(4)*rad2deg, &
            kpl(5)*rad2deg, kpl(6)*rad2deg
        
        write(unit_num, '(A)') "# Epoch(MJD), X(km), Y(km), Z(km), VX(km/s), VY(km/s), VZ(km/s)"
        write(unit_num, '(F12.6,6(1X,F15.8))') epoch%val, &
            cts(1:3)*length_unit, &
            cts(4:6)*length_unit/time_unit
        
        write(unit_num, *) ""
        write(unit_num, *) "=============================================================="
        write(unit_num, *) "End of Report"
        write(unit_num, *) "=============================================================="
        
        ! 关闭文件
        close(unit_num)
        
        print*, "Results saved to file: ", trim(filename)
        
    end subroutine save_IODor

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! General Laplace method for initial orbit determination    
    subroutine initialOrbitDetermination_generalLaplaceMethod(nD, t0, t, alpha, delta, X, Y, Z, cts0)
    use module_const
    implicit none
    integer, intent(in) :: nD
    double precision, intent(in) :: t0, t(nD), alpha(nD), delta(nD), X(nD), Y(nD), z(nD)
    double precision, intent(out) :: cts0(6)
    integer i, j
    double precision F(nD), G(nD), Fz(nD), Gz(nD), A(2*nD, 6), b(2*nD, 1), lastCts0(6)
    double precision tau(nD), lambda(nD), mu(nD), nu(nD)
    double precision Fx, Fy, Gx, Gy
    !--------------------------------------------------------------------------
    !>>>>>1 ��ֵtau, lambda, mu, nu
    tau = (t - t0)*86400d0/time_unit
    do i = 1, nD
        lambda(i) = dcos(delta(i))*dcos(alpha(i))
        mu(i) = dcos(delta(i))*dsin(alpha(i))
        nu(i) = dsin(delta(i))
    enddo
    
    !>>>>>2 ��ʼ�������������10��
    do i = 1, 50
        !>>>>>3 ��ֵF, G
        if (i == 1) then
            F = 1d0
            Fz = 1d0
            G = tau
            Gz = tau
        else
            do j = 1, nD
                !call calculateFG(t0, tau(j), cts0, F(j), G(j))
                call calculateFG2(t0, tau(j), cts0, .true., .true., .true., .false., 0d0, &
        &   Fx, Gx, Fy, Gy, Fz(j), Gz(j))
                F(j) = Fx
                G(j) = Gx
            enddo
        endif
        
        !>>>>>4 �༭����A��b��A*X = b
        do j = 1, nD
            A(2*j - 1, 1:3) = (/F(j)*nu(j), 0d0,    -Fz(j)*lambda(j)/)
            A(2*j - 1, 4:6) = (/G(j)*nu(j), 0d0,    -Gz(j)*lambda(j)/)
            A(2*j, 1:3)     = (/0d0,    F(j)*nu(j),     -Fz(j)*mu(j)/)
            A(2*j, 4:6)     = (/0d0,    G(j)*nu(j),     -Gz(j)*mu(j)/)
            b(2*j - 1, 1)   = nu(j)*X(j) - lambda(j)*Z(j)
            b(2*j, 1)       = nu(j)*Y(j) - mu(j)*Z(j)
        enddo
        
        !>>>>>5 ��С���˽�X
        lastCts0 = cts0 !�����ϴε����Ľ��
        call least_squares_solution(2*nD, 6, A, b, cts0)
        
        print"(6(2x,F10.5))", cts0(1:6)
        
        !>>>>>6 �ж��Ƿ�����ѭ��
        if (norm2(cts0 - lastCts0) <= 1d-7) exit
    enddo
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate F and G, r = F*r0 + G*v0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calculateFG(t0, tau, X, F, G)
    use module_planets
    implicit none
    double precision, intent(in) :: t0, tau
    double precision, dimension(6), intent(in) :: X
    double precision :: F, G
    double precision :: u1, u2, u3, u4, u5, u6, u7, u8, u9, sigma, v02
    double precision :: tau2, tau3, tau4, tau5, tau6
    double precision, dimension(3) :: r_moon
    double precision :: u3_moon, r0_moon
    double precision, dimension(6) :: PV
    ! Calculate u values
    u1 = 1.0d0 / dsqrt(X(1)**2 + X(2)**2 + X(3)**2)
    u2 = u1 * u1
    u3 = u2 * u1
    u4 = u3 * u1
    u5 = u4 * u1
    u6 = u5 * u1
    u7 = u6 * u1
    u8 = u7 * u1
    u9 = u8 * u1
    sigma = X(1) * X(4) + X(2) * X(5) + X(3) * X(6)
    v02 = X(4)**2 + X(5)**2 + X(6)**2
    
    ! Calculate tau powers
    tau2 = tau**2
    tau3 = tau2 * tau
    tau4 = tau3 * tau
    tau5 = tau4 * tau
    tau6 = tau5 * tau
    
    ! Calculate F and G with pure two-body terms
    ! F
    F = 1.0d0 + tau2 / 2.0d0 * (-u3) + tau3 / 6.0d0 * (3.0d0 * u5 * sigma) + &
        tau4 / 24.0d0 * (u5 * (3.0d0 * v02 - 2.0d0 * u1 - 15.0d0 * u2 * sigma**2)) + &
        tau5 / 120.0d0 * u7 * (15.0d0 * sigma * (-3.0d0 * v02 + 2.0d0 * u1 + 7.0d0 * u2 * sigma**2)) + &
        tau6 / 720.0d0 * u7 * (u2 * sigma**2 * (630.0d0 * v02 - 420.0d0 * u1 - 945.0d0 * u2 * sigma**2) - &
        (22.0d0 * u2 - 66.0d0 * u1 * v02 + 45.0d0 * v02**2))
    
    ! G
    G = tau + tau3 / 6.0d0 * (-u3) + tau4 / 24.0d0 * (6.0d0 * u5 * sigma) + &
        tau5 / 120.0d0 * u5 * (9.0d0 * v02 - 8.0d0 * u1 - 45.0d0 * u2 * sigma**2) + &
        tau6 / 720.0d0 * u7 * (30.0d0 * sigma * (-6.0d0 * v02 + 5.0d0 * u1 + 14.0d0 * u2 * sigma**2))
    
    ! Calculate moon-related terms
    call planetsPV(t0, nEarth, nMoon, PV)
    r_moon = PV(1:3)
    r0_moon = dsqrt(r_moon(1)**2 + r_moon(2)**2 + r_moon(3)**2)
    u3_moon = 1.0d0 / r0_moon**3
    
    ! Add moon terms
    F = F + tau2 / 2.0d0 * (-moonEarthMassRatio * u3_moon) + tau4 / 24.0d0 * (u3 * moonEarthMassRatio * u3_moon)
    G = G + tau3 / 6.0d0 * (-moonEarthMassRatio * u3_moon)
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate F and G function
    ! inputs:
    !   t0          utc in mjd
    !   tau         dt from t0
    !   X           position and velocity, normalized units
    !   useJ2       consider J2 term or not
    !   useMoon     consider third body gravity (the moon) or not
    !   useSun      consider third body gravity (the sun) or not
    !   useSRP      consider solar radiation pressure or not
    ! outputs:
    !   F, G, Fz, Gz    4 functions used to construct r and v of spacecraft, 
    !       [x; y; z] = [Fx, Fy, Fz]*[x0; y0; z0] + [Gx, Gy, Gz]*[vx0; vy0; vz0]
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calculateFG2(t0, tau, X, useJ2, useMoon, useSun, useSRP, A, &
        &   Fx, Gx, Fy, Gy, Fz, Gz)
    use module_planets
    implicit none
    double precision, intent(in) :: t0, tau
    double precision, dimension(6), intent(in) :: X
    logical, intent(in) :: useJ2, useMoon, useSun, useSRP
    double precision, intent(in) :: A
    double precision, intent(out) :: Fx, Gx, Fy, Gy, Fz, Gz
    double precision :: u1, u2, u3, u4, u5, u6, u7, u8, u9, sigma, v02
    double precision :: tau2, tau3, tau4, tau5, tau6
    double precision, dimension(3) :: r_moon, r_sun
    double precision :: u3_moon, r0_moon, u3_sun, r0_sun
    double precision, dimension(6) :: PV
    double precision x0, y0, z0, x0d, y0d, z0d, dF, dG, J2, F, G
    double precision delta(3), deltaU(3), nu
    
    !>>>>>1 prepare some inetrmediate values
    u1 = 1.0d0 / dsqrt(X(1)**2 + X(2)**2 + X(3)**2) ! u_n = 1/r0**n
    u2 = u1 * u1
    u3 = u2 * u1
    u4 = u3 * u1
    u5 = u4 * u1
    u6 = u5 * u1
    u7 = u6 * u1
    u8 = u7 * u1
    u9 = u8 * u1
    sigma = X(1) * X(4) + X(2) * X(5) + X(3) * X(6) ! sigma = r0'*v0
    v02 = X(4)**2 + X(5)**2 + X(6)**2 ! v02 = v0**2
    x0 = X(1)
    y0 = X(2)
    z0 = X(3)
    x0d = X(4)
    y0d = X(5)
    z0d = X(6)
    tau2 = tau**2
    tau3 = tau2 * tau
    tau4 = tau3 * tau
    tau5 = tau4 * tau
    tau6 = tau5 * tau
    
    !>>>>>2 Calculate F, G, Fz, Gz for pure two-body terms
    F = 1.0d0 + &
        tau2 / 2.0d0 * (-u3) + &
        tau3 / 6.0d0 * (3.0d0 * u5 * sigma) + &
        tau4 / 24.0d0 * (u5 * (3.0d0 * v02 - 2.0d0 * u1 - 15.0d0 * u2 * sigma**2)) + &
        tau5 / 120.0d0 * u7 * (15.0d0 * sigma * (-3.0d0 * v02 + 2.0d0 * u1 + 7.0d0 * u2 * sigma**2)) + &
        tau6 / 720.0d0 * u7 * (u2 * sigma**2 * (630.0d0 * v02 - 420.0d0 * u1 - 945.0d0 * u2 * sigma**2) & 
        &   - (22.0d0 * u2 - 66.0d0 * u1 * v02 + 45.0d0 * v02**2))
    
    G = tau + &
        tau3 / 6.0d0 * (-u3) + tau4 / 24.0d0 * (6.0d0 * u5 * sigma) + &
        tau5 / 120.0d0 * u5 * (9.0d0 * v02 - 8.0d0 * u1 - 45.0d0 * u2 * sigma**2) + &
        tau6 / 720.0d0 * u7 * (30.0d0 * sigma * (-6.0d0 * v02 + 5.0d0 * u1 + 14.0d0 * u2 * sigma**2))
    
    Fx = F
    
    Gx = G
    
    Fy = F
    
    Gy = G
    
    Fz = F
    
    Gz = G
    
    !>>>>>3 add J2 term
    if (useJ2) then
        J2 = 1.082629989d-3
        dF = tau2/2d0*(3d0*J2/2d0)*(5d0*u7*z0**2 - u5) + &
             tau3/6d0*(3d0*J2/2d0)*(5*(u7 - 7d0*u9*z0**2)*sigma + 10d0*u7*z0*z0d) + &
             tau4/24d0*(3d0*J2/2d0)*(6d0*u8*(4d0*u2*z0**2 - 1d0) - &
                5d0*u7*(7d0*u2*z0**2 - 1d0)*v02 + 10*u7*z0d**2 + &
                35d0*u9*(9*u2*z0**2 - 1d0)*sigma**2 - 140d0*u9*sigma*z0*z0d)
        
        dG = tau3/6d0*(3d0*J2/2d0)*(5*u7*z0**2 - u5) + &
             tau4/24d0*(3d0*J2/2d0)*(20d0*u7*z0*z0d - 10d0*u7*(7d0*u2*z0**2 - 1)*sigma)
        
        Fx = Fx + dF
        
        Gx = Gx + dG
        
        Fy = Fy + dF
        
        Gy = Gy + dG
        
        Fz = Fz + dF + (3d0*J2/2d0)*(tau2/2d0*(-2d0*u5) + tau3/6d0*(10d0*u7*sigma) + &
            tau4/24d0*u7*(10d0*v02 - 6d0*u1 - 70d0*u2*sigma**2))
        
        Gz = Gz + dG + (3d0*J2/2d0)*(tau3/2d0*(-2d0*u5) + tau4/24d0*(20d0*u7*sigma))
    endif
    
    !>>>>>4 add lunar gravity term
    if (useMoon) then
        call planetsPV(t0, nEarth, nMoon, PV)
        r_moon = PV(1:3)
        r0_moon = dsqrt(r_moon(1)**2 + r_moon(2)**2 + r_moon(3)**2)
        u3_moon = 1.0d0 / r0_moon**3
        
        dF = tau2 / 2.0d0 * (-moonEarthMassRatio * u3_moon) + &
             tau4 / 24.0d0 * (u3 * moonEarthMassRatio * u3_moon)
        
        dG = tau3 / 6.0d0 * (-moonEarthMassRatio * u3_moon)
        
        Fx = Fx + dF
        
        Gx = Gx + dG
        
        Fy = Fy + dF
        
        Gy = Gy + dF
        
        Fz = Fz + dF
        
        Gz = Gz + dG
    endif
    
    !>>>>>5 add solar gravity term
    if (useSun) then
        call planetsPV(t0, nEarth, nSun, PV)
        r_sun = PV(1:3)
        r0_sun = dsqrt(r_sun(1)**2 + r_sun(2)**2 + r_sun(3)**2)
        u3_sun = 1.0d0 / r0_sun**3
        
        dF = tau2 / 2.0d0 * (-sunEarthMassRatio * u3_sun) + &
             tau4 / 24.0d0 * (u3 * sunEarthMassRatio * u3_sun)
        
        dG = tau3 / 6.0d0 * (-sunEarthMassRatio * u3_sun)
        
        Fx = Fx + dF
        
        Gx = Gx + dG
        
        Fy = Fy + dF
        
        Gy = Gy + dG
        
        Fz = Fz + dF
        
        Gz = Gz + dG
    endif
    
    !>>>>>6 add SRP term
    if (useSRP) then
        !get solar position
        call planetsPV(t0, nEarth, nSun, PV)
        r_sun = PV(1:3)        
        delta = X(1:3) - r_sun
        deltaU = delta/norm2(delta)
        !calculate earth shadow coefficient nu
        call EarthShadowCoefficient(X(1:3), r_sun, nu)
        !refresh F G
        Fx = Fx + tau2/2d0*nu*A*deltaU(1)/x0
        Fy = Fy + tau2/2d0*nu*A*deltaU(2)/y0
        Fz = Fz + tau2/2d0*nu*A*deltaU(3)/z0
    endif
    
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine EarthShadowCoefficient(r_sat, r_sun, nu)
    implicit none
    double precision nu
    double precision r_sat(3), r_sun(3), r_sun_unit(3), r_sat_unit(3)
    double precision sin_phi, cos_phi
    !---------------------------------------------------------------
    r_sun_unit = r_sun/norm2(r_sun)
    r_sat_unit = r_sat/norm2(r_sat)
    cos_phi = dot_product(r_sun_unit, r_sat_unit)
    sin_phi = dsqrt(1d0 - cos_phi**2)
    if ((cos_phi < 0) .and. (sin_phi < 1d0/norm2(r_sat))) then
        nu = 0d0
    else
        nu = 1d0
    end if
    end subroutine
    
    endmodule