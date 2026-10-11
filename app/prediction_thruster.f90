
    program orbit_prediction_thrust
    use module_satellite
    use module_thruster
    use module_time
    use module_const
    implicit none
    
    ! Satellite object
    type(type_satellite) :: sat_with_thruster, sat_without_thrust
    
    ! Thruster object
    type(type_thruster) :: thruster
    
    ! Time parameters
    type(type_t) :: start_time, end_time
    type(type_t), allocatable :: epochs(:)
    
    ! Orbital parameters - 400km altitude near-circular orbit
    real*8 :: kepler_elements(6)
    real*8 :: nsg1_elements(6)
    
    ! Other variables
    integer :: i, num_steps
    real*8 :: step
    character(len=255) :: filename_with_thrust, filename_without_thrust
    
    ! Initialize time parameters
    print *, "=== Low-Thrust Orbit Prediction Demo Program ==="
    print *, "Scenario: Low-thrust orbit transfer at 400km altitude"
    
    ! Set initial time: January 1, 2025 00:00:00
    call start_time%init(2025d0, 1d0, 1d0, 0d0, 0d0, 0d0)
    
    ! Set end time: January 2, 2025 00:00:00 (24-hour prediction)
    call end_time%init(2025d0, 1d0, 2d0, 0d0, 0d0, 0d0)
    
    ! Set time step: 60 seconds
    step = 60.0d0
    
    ! Calculate number of steps
    num_steps = ceiling(abs(timeMinus(end_time, start_time) * 86400.0d0 / step)) + 1
    
    ! Allocate and generate time series
    allocate(epochs(num_steps))
    epochs(1) = start_time
    do i = 2, num_steps
        epochs(i)%day = start_time%day
        epochs(i)%sec = start_time%sec + (i-1) * step
        ! Time normalization
        do while (epochs(i)%sec >= 86400.0d0)
            epochs(i)%day = epochs(i)%day + 1
            epochs(i)%sec = epochs(i)%sec - 86400.0d0
        end do
        epochs(i)%val = epochs(i)%day + epochs(i)%sec/86400.0d0
    end do
    
    ! Set Keplerian elements for 400km altitude orbit
    ! Earth radius: 6378.1363 km, orbit altitude: 400 km
    ! Semi-major axis: 6778.1363 km, eccentricity: 0.001 (near-circular orbit)
    kepler_elements(1) = 6778.1363d0  ! Semi-major axis (km)
    kepler_elements(2) = 0.001d0      ! Eccentricity
    kepler_elements(3) = 45.0d0       ! Inclination (degrees)
    kepler_elements(4) = 0.0d0        ! Right ascension of ascending node (degrees)
    kepler_elements(5) = 0.0d0        ! Argument of perigee (degrees)
    kepler_elements(6) = 0.0d0        ! Mean anomaly (degrees)
    
    print *, "Initial orbital parameters:"
    print *, "  Semi-major axis: ", kepler_elements(1), " km"
    print *, "  Eccentricity: ", kepler_elements(2)
    print *, "  Inclination: ", kepler_elements(3), " degrees"
    print *, "  Orbit altitude: 400 km"
    
    ! Convert to normalized units
    call kpl_SI2normal(kepler_elements, kepler_elements)
    call kpl2nsg1(kepler_elements, nsg1_elements)
    
    ! =================================================================
    ! Case 1: Orbit prediction with thruster
    ! =================================================================
    print *, ""
    print *, "=== Case 1: Orbit prediction with low-thruster ==="
    
    ! Initialize satellite (with thruster)
    call sat_with_thruster%init(start_time, nsg1_elements, 1, "Satellite_with_Thrust", 1000.0d0, 3, "numerical", &
                             Cd=2.2d0, amrDrag=0.02d0, F107=130.0d0, Ap=15.0d0, &
                             reflectivity=0.0d0, amrSRP=0.02d0, degree=2, order=2)
    sat_with_thruster%opMod%key_drag = 1   ! Only atmospheric drag enabled
    
    ! Initialize thruster.
    ! 100 kg of propellant: a day at 20 N costs 19 N / (3000 s * 9.7983 m/s^2)
    ! * 86400 s = 55.8 kg, so the tank lasts the whole run and the mass the force
    ! model uses drops from 1000 kg to ~944 kg as the engine burns.
    call thruster%init(id=1, name="Small_Thruster", &
                      specific_impulse=3000.0d0, &    ! Specific impulse 3000 seconds
                      max_thrust=1000d0, &             ! Maximum thrust 1000N
                      coordinate_frame=3, &           ! NTH orbital coordinate system
                      thrust_direction=[0.0d0, 1.0d0, 0.0d0], &  ! Along velocity direction
                      efficiency=0.95d0, &            ! Efficiency 95%
                      min_thrust=0.001d0, &            ! Minimum thrust 0.001N
                      response_time=0d0, &            ! Response time 0 seconds
                      initial_fuel_mass=100.0d0)      ! Propellant load 100 kg

    ! Set thrust magnitude and activate thruster
    call thruster%set_current_thrust(20d0)  ! Set thrust to 10N
    call thruster%activate(start_time)       ! 引擎，启动！！！（记录开机时刻）
    
    ! Add thruster to satellite
    call sat_with_thruster%addThruster(thruster)
    
    ! Display thruster information
    print *, "Thruster configuration:"
    print *, "  Thrust magnitude: 10 N"
    print *, "  Thrust direction: Along velocity direction (NTH coordinate system [0, 1, 0])"
    print *, "  Specific impulse: 3000 seconds"
    print *, "  Application: Orbit altitude transfer"
    
    ! Display satellite information
    call sat_with_thruster%show()
    
    ! Execute orbit prediction
    print *, "Starting orbit prediction with thruster..."
    call sat_with_thruster%generateTrajectory(epochs)
    
    ! Save results
    filename_with_thrust = "data/output/orbit_with_thrust.txt"
    call sat_with_thruster%savetrajectory(STATE_TYPE_CARTESIAN, filename_with_thrust, .false.)
    print *, "Orbit prediction results with thruster saved to: ", trim(filename_with_thrust)
    
    ! =================================================================
    ! Case 2: Orbit prediction without thruster (comparison)
    ! =================================================================
    print *, ""
    print *, "=== Case 2: Orbit prediction without thruster (comparison) ==="
    
    ! Initialize satellite (without thruster)
    call sat_without_thrust%init(start_time, nsg1_elements, 2, "Satellite_without_Thrust", 1000.0d0, 3, "numerical", &
                                Cd=2.2d0, amrDrag=0.02d0, F107=130.0d0, Ap=15.0d0, &
                                reflectivity=0.0d0, amrSRP=0.02d0, degree=2, order=2)
    sat_without_thrust%opMod%key_drag = 1   ! Only atmospheric drag enabled
    
    ! Display satellite information
    call sat_without_thrust%show()
    
    ! Execute orbit prediction
    print *, "Starting orbit prediction without thruster..."
    call sat_without_thrust%generateTrajectory(epochs)
    
    ! Save results
    filename_without_thrust = "data/output/orbit_without_thrust.txt"
    call sat_without_thrust%savetrajectory(STATE_TYPE_CARTESIAN, filename_without_thrust, .false.)
    print *, "Orbit prediction results without thruster saved to: ", trim(filename_without_thrust)
    
    ! =================================================================
    ! Result analysis
    ! =================================================================
    print *, ""
    print *, "=== Result Analysis ==="
    
    ! Get initial and final states for comparison
    if (sat_with_thruster%nEph > 0 .and. sat_without_thrust%nEph > 0) then
        ! Calculate initial orbit altitude
        call nsg12kpl(sat_with_thruster%states%elements, kepler_elements)
        call kpl_normal2SI(kepler_elements, kepler_elements)
        print *, "Initial orbit altitude: ", (kepler_elements(1)*(1.0d0-kepler_elements(2)) - 6378.1363d0), " km"
        
        ! Calculate final orbit altitude (with thruster)
        call nsg12kpl(sat_with_thruster%trajectory(sat_with_thruster%nEph)%elements, kepler_elements)
        call kpl_normal2SI(kepler_elements, kepler_elements)
        print *, "Final orbit altitude (with thruster): ", (kepler_elements(1)*(1.0d0-kepler_elements(2)) - 6378.1363d0), " km"
        
        ! Calculate final orbit altitude (without thruster)
        call nsg12kpl(sat_without_thrust%trajectory(sat_without_thrust%nEph)%elements, kepler_elements)
        call kpl_normal2SI(kepler_elements, kepler_elements)
        print *, "Final orbit altitude (without thruster): ", (kepler_elements(1)*(1.0d0-kepler_elements(2)) - 6378.1363d0), " km"
        
        ! Display thruster operation statistics.
        ! The ledger is caller-driven: nothing inside the integrator accrues, so it
        ! must be settled up to the end of the burn before the totals are read.
        ! Work on the satellite's own copy of the thruster: `thruster` and
        ! sat_with_thruster%thruster are independent objects from the moment
        ! addThruster copied it.
        print *, ""
        print *, "Thruster operation statistics:"
        call sat_with_thruster%thruster%update_state(end_time)
        print *, "  Total operating time: ", &
                 sat_with_thruster%thruster%total_operation_time, " seconds"
        print *, "  Total impulse: ", &
                 sat_with_thruster%thruster%get_total_impulse(end_time), " N*s"
        print *, "  Total fuel consumption: ", &
                 sat_with_thruster%thruster%get_total_fuel_used(end_time), " kg"
        print *, "  Fuel remaining: ", &
                 sat_with_thruster%thruster%current_fuel_mass, " kg"
        print *, "  Satellite mass at the end: ", &
                 sat_with_thruster%mass - sat_with_thruster%thruster%total_fuel_used, " kg"
    end if
    
    print *, ""
    print *, "=== Low-Thrust Orbit Prediction Demo Completed ==="
    print *, "Result files:"
    print *, "  With thruster: ", trim(filename_with_thrust)
    print *, "  Without thruster: ", trim(filename_without_thrust)
    print *, "These files can be used to compare the effects of low-thrust on orbit."
    
    ! Clean up memory
    deallocate(epochs)
    
end program orbit_prediction_thrust
