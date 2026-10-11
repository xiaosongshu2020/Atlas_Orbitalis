!==============================================================================
! Module: module_station
! Purpose: Defines ground station object for satellite observation and orbit determination
! Last Updated: 2025-12-22 by Shawn Hu
!==============================================================================
!
! Type: type_station
! ------------------
! Properties:
!   - number      : integer       - Station ID
!   - name        : character*100 - Station name
!   - longitude   : double precision        - Geographic longitude (radius, East positive)
!   - latitude    : double precision        - Geographic latitude (radius, North positive)
!   - altitude    : double precision        - Altitude above sea level (Earth radii)
!   - position(3) : double precision        - Position in Earth-Centered Fixed frame (Earth radii)
!   - ts          : type_timeSeries - Associated time series
!   - trajectory(:, :): double precision, allocatable - Trajectory in J2000 frame [x,y,z,vx,vy,vz]
!   - clock       : type_clock    - Station clock model (bias/drift + clock table)
!
! Methods:
!   - init      : Initialize station with basic parameters
!   - generateTrajectory : Generate station trajectory in J2000 frame
!   - getPV     : Get position/velocity(normalized unit) at specific time in J2000 frame
!   - show      : Display station information
!
! Usage Example:
! -------------
! type(type_station) :: station
! type(type_t) :: t0, t
! type(type_timeSeries) :: ts
! double precision :: pv(6)
! 
! call station%init(number=1, name="Beijing", &
!                  longitude=116.4d0, latitude=39.9d0, altitude=50.0d0)
! call t0%init(2025d0, 1d0, 1d0, 0d0, 0d0, 0d0) ! "2025-01-01 00:00:00"
! call ts%init(t0, 60d0, 1440) ! 1 day of data at 1-minute intervals
! call station%generateTrajectory(ts) ! Generate station's trajectory in J2000 frame
! call station%show()
!
! call t%init(2025d0, 1d0, 1d0, 0d0, 5d0, 10d0) ! "2025-01-01 00:05:10"
! pv = station%getPV(t) ! Get position/velocity at specified time
!
!==============================================================================

module module_station

    use module_time
    use module_timeSeries
    use module_clock

    private :: Geodetic, geo2pos, pos2geo
    
    type, public ::  type_station
        integer         ::  number = 0      ! Station number/ID
        character*100   ::  name = ""       ! Station name
        double precision          ::  longitude = 0d0 ! Geographic longitude (rad)
        double precision          ::  latitude = 0d0  ! Geographic latitude (rad)
        double precision          ::  altitude = 0d0  ! Altitude above sea level (Earth radius)
        double precision          ::  position(3) = 0d0 ! Position in Earth-Centered Fixed frame (units of Earth radius)
        type(type_timeSeries) :: ts         ! Time series associated with the station
        double precision, allocatable ::  trajectory(:, :)   ! trajectory data (position & velocity in J2000 frame)
        type(type_clock) :: clock           ! Station clock model (bias/drift + clock table)
    contains
        procedure :: init => initStation    ! Initialize the station object
        procedure :: generateTrajectory           ! Generate station trajectory in J2000 coordinate system
        procedure :: getPV => station_getPV         ! Get position/velocity at a specific time in J2000 system
        procedure :: show => showStation   ! Display station information
    end type type_station
    
    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize a station with basic parameters
    ! This procedure sets the primary properties of a station and calculates the position vector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine initStation(this, number, name, longitude, latitude, altitude)
        implicit none
        class(type_station), intent(inout) :: this
        integer, intent(in), optional :: number           ! Station number/ID
        character(len=*), intent(in), optional :: name    ! Station name
        double precision, intent(in), optional :: longitude         ! Geographic longitude (degrees)
        double precision, intent(in), optional :: latitude          ! Geographic latitude (degrees)
        double precision, intent(in), optional :: altitude          ! Altitude above sea level (meters)
        
        ! Set station properties if provided
        if (present(number)) this%number = number
        if (present(name)) this%name = name
        if (present(longitude)) this%longitude = longitude * deg2rad
        if (present(latitude)) this%latitude = latitude * deg2rad
        if (present(altitude)) this%altitude = altitude * u_meter
        
        ! Convert geographic coordinates to ECF position if all necessary params are provided
        if (present(longitude) .and. present(latitude) .and. present(altitude)) then
            call geo2pos(this)
        endif
        
    end subroutine initStation
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Generate trajectory data for the station in Earth-Centered Inertial (J2000) frame
    ! This transforms the station's fixed position in ECF to multiple positions in ECI
    ! accounting for Earth's rotation at each time point
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine generateTrajectory(this, ts)
        use module_space
        class (type_station) :: this
        type(type_timeSeries), intent(in) :: ts  ! Time series for which to generate trajectory
        integer i, nEpoch
        double precision rECF(6), rECI(6)  ! Position and velocity vectors in ECF and ECI frames
        !---------------------------------------------------------------------------
        
        ! Error check: ensure position is defined
        if (norm2(this%position) == 0d0) then
            print*, "Error @subroutine generateTrajectory"
            print*, "station's position is unknown."
            return
        endif
        
        this%ts = ts
        
        ! Allocate memory for trajectory data: position and velocity (6 values) at each epoch
        nEpoch = size(ts%time)
        if (allocated(this%trajectory)) deallocate(this%trajectory)
        allocate(this%trajectory(nEpoch, 6))
        this%trajectory = 0d0
        
        ! Calculate position and velocity in ECI frame for each time point
        do i = 1, nEpoch
            rECF = (/ this%position, 0d0, 0d0, 0d0 /)  ! Station has velocity=0 in ECF
            call ECF2ECI_IAU1980_positionVelocity(ts%time(i)%val, rECF(1:3), rECF(4:6), &
                &   rECI(1:3), rECI(4:6))
            this%trajectory(i, 1:6) = rECI
        enddo
        
    endsubroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Convert geographic coordinates (longitude, latitude, altitude) to 
    ! Earth-Centered Fixed (ECF) position vector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine geo2pos(this)
        use module_const
        implicit none
        class(type_station) :: this
        double precision      longitude, latitude, altitude
        double precision      position(3), NN, epsi
        !-----------------------------------------------------------------
        longitude = this%longitude
        latitude = this%latitude
        altitude = this%altitude
        
        epsi        =   EARTH_FLATTENING            ! Earth's flattening parameter
        ! Calculate radius of curvature in the prime vertical
        NN          =   ((dcos(latitude))**2d0 + (1d0 - epsi)**2d0*(dsin(latitude))**2d0)**(-0.5d0)
        
        ! Calculate cartesian coordinates in ECF frame
        position(1) =   (NN + altitude)*dcos(latitude)*dcos(longitude)  ! X coordinate
        position(2) =   (NN + altitude)*dcos(latitude)*dsin(longitude)  ! Y coordinate
        position(3) =   (NN*(1 - epsi)**2 + altitude)*dsin(latitude)    ! Z coordinate
        
        this%position = position
        
        
        
    endsubroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Convert Earth-Centered Fixed (ECF) position vector to 
    ! geographic coordinates (longitude, latitude, altitude)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pos2geo(this)
        use module_const
        implicit none
        class(type_station) :: this
        double precision :: r(3), lon, lat, h
        r = this%position*length_unit*1d3  ! Convert to meters
        call geodetic(r, lon, lat, h)      ! Call helper routine to do the conversion
        this%longitude = lon       
        this%latitude = lat        
        this%altitude = h / u_meter
    endsubroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Convert Earth-Centered Fixed (ECF) cartesian coordinates to geodetic coordinates
    ! Uses iterative algorithm to handle Earth's ellipticity
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Geodetic(r, lon, lat, h)
        use module_const
        implicit none
        double precision, intent(in)  :: r(3)   ! Position vector [m]
        double precision, intent(out) :: lon    ! Longitude [rad]
        double precision, intent(out) :: lat    ! Latitude [rad]
        double precision, intent(out) :: h      ! Height [m]
        
        double precision :: R_equ, f, epsRequ, e2
        double precision :: X, Y, Z, rho2, dZ, ZdZ, Nh, SinPhi, N, dZ_new
    
        ! Get Earth parameters from module
        R_equ = EARTH_RADIUS_EQUATORIAL*1d3      ! Equatorial radius in meters
        f = EARTH_FLATTENING              ! Earth's flattening parameter
    
        ! Calculate convergence threshold and squared eccentricity
        epsRequ = epsilon(1.d0) * R_equ  ! Convergence threshold based on machine precision
        e2 = f * (2.d0 - f)              ! First eccentricity squared
    
        ! Decompose position vector into components
        X = r(1)
        Y = r(2)
        Z = r(3)
        rho2 = X**2 + Y**2      ! Square of distance from Z-axis
    
        ! Check for null vector input
        if (X == 0.d0 .and. Y == 0.d0 .and. Z == 0.d0) then
            write(*,*) 'Invalid input in Geodetic constructor'
            lon = 0.d0
            lat = 0.d0
            h = -R_equ
            return
        endif
    
        ! Initialize height correction term
        dZ = e2 * Z
    
        ! Iterative algorithm to solve for latitude
        do
            ZdZ = Z + dZ
            Nh = sqrt(rho2 + ZdZ**2)
            SinPhi = ZdZ / Nh             ! Sine of geodetic latitude
            N = R_equ / sqrt(1.d0 - e2*SinPhi**2)  ! Radius of curvature in prime vertical
            dZ_new = N * e2 * SinPhi      ! New height correction term
            
            ! Check for convergence
            if (abs(dZ - dZ_new) < epsRequ) exit
            dZ = dZ_new
        end do
    
        ! Calculate final geodetic coordinates
        lon = atan2(Y, X)               ! Longitude (-pi, pi]
        lat = atan2(ZdZ, sqrt(rho2))    ! Latitude (-pi/2, pi/2)
        h = Nh - N                      ! Height above ellipsoid
    
    end subroutine Geodetic
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get position and velocity of the station at a specific time
    ! Uses interpolation if the time is not exactly at an trajectory point
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function station_getPV(this, t) result(PV)
        implicit none
        class(type_station) :: this
        type(type_t), intent(in) :: t   ! Time at which to get position/velocity
        double precision :: PV(6)                 ! Result: position (1:3) and velocity (4:6)
        double precision :: ntDble, step
        integer :: nEph, n
        type(type_t) :: t0
        integer nt, left1, left2, right1, right2, j
        double precision x(4), y(4), yt(6)
        !---------------------------------------------------------------------------
        ! Check if trajectory data is available
        if (.not.allocated(this%trajectory)) then
            print*, "error @function getPV."
            print*, "station trajectory is not been prepared."
            stop
        endif
        
        nEph = size(this%ts%time)
        
        ! Check if requested time is within the trajectory span
        if (timeMinus(t, this%ts%time(1)) < 0d0 .or. &
            &   timeMinus(t, this%ts%time(nEph)) > 0d0) then
            print*, "error @function getPV."
            print*, "trajectory is not been prepared at t = ", t%val
            stop
        endif
        
        t0 = this%ts%time(1)          ! Initial time of trajectory
        step = timeMinus(this%ts%time(2), this%ts%time(1))  ! Time step size
        
        ! Calculate the fractional index corresponding to time t
        ntDble = timeMInus(t, t0)/step
    
        if (dabs(ntDble - dble(nint(ntDble))) <= 1d-12) then
            ! If t matches exactly an trajectory time point, use that point directly
            n = nint(ntDble) + 1
            PV = this%trajectory(n, 1:6)
        else
            ! Otherwise, use interpolation between neighboring points
            n = floor(ntDble) + 1
            
            ! Select four points for interpolation, handling boundary cases
            if (n == 1 .or. n == 2) then
                left1 = 1
                left2 = 2
                right1 = 3
                right2 = 4
            elseif (n == nEph .or. n == nEph - 1) then
                left1 = nEph - 3
                left2 = nEph - 2
                right1 = nEph - 1
                right2 = nEph
            else
                left1 = n - 1
                left2 = n
                right1 = n + 1
                right2 = n + 2
            endif
            
            ! Interpolate each component of position and velocity
            do j = 1, 6
                x(1:4) = this%ts%time([left1, left2, right1, right2])%val
                y(1:4) = this%trajectory([left1, left2, right1, right2], j)
                call Interpolation(4, x(1:4), y(1:4), 1, (/t%val/), PV(j))
            enddo
        endif
        
    endfunction
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Display station information in a formatted way
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine showStation(this)
        implicit none
        class(type_station) :: this
        print*, "--------------------------------------------------------------"
        print*, "Station Number: ", this%number
        print*, "Station Name: ", trim(this%name)
        print*, "Longitude (deg): ", this%longitude * rad2deg
        print*, "Latitude (deg): ", this%latitude * rad2deg
        print*, "Altitude (m): ", this%altitude / u_meter
        print*, "Position in Earth-Fixed Frame (m): ", this%position*EARTH_RADIUS_EQUATORIAL*1d3
        print*, "--------------------------------------------------------------"
    endsubroutine
    
endmodule
