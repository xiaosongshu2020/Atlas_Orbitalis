
    
    program orbit_prediction

    use module_Atlas
    
    implicit none
    
    ! Variables for command line argument parsing
    integer :: i, argc, status
    character(len=255) :: arg
    
    ! All parameters (no optional parameters)
    type(type_t) :: start_time, end_time
    real*8 :: kepler_elements(6), nsg1(6)
    real*8 :: step
    character(len=255) :: trajectory_filename
    integer :: ephType
    real*8 :: Cd
    real*8 :: amrDrag
    real*8 :: F107
    real*8 :: Ap
    real*8 :: reflectivity
    real*8 :: amrSRP
    integer :: degree
    integer :: order
    integer :: K_lunar, K_lunarTide, K_solor, K_solarTide, K_SRP, K_drag, K_PN
    
    ! Fixed parameters
    character(len=100) :: name = ""
    integer :: eleType = 3
    character(len=10) :: propagator = "numerical"
    
    ! Other variables
    type(type_satellite) :: sat
    integer :: num_steps
    type(type_t), allocatable :: epochs(:)
    
    ! Get number of command line arguments
    argc = command_argument_count()
    
    ! Check if we have exactly 21 parameters
    if (argc /= 21) then
        call print_usage()
        print *, "Error: Expected 21 parameters, got ", argc
        stop
    end if
    
    ! Parse all parameters in order
    ! Parameter 1: start_time
    call get_command_argument(1, arg)
    call parseDateTime(trim(arg), start_time)
    
    ! Parameter 2: kepler_elements
    call get_command_argument(2, arg)
    call parse_elements(trim(arg), kepler_elements)
    
    ! Parameter 3: end_time
    call get_command_argument(3, arg)
    call parseDateTime(trim(arg), end_time)
    
    ! Parameter 4: step
    call get_command_argument(4, arg)
    read(arg, *, iostat=status) step
    if (status /= 0) then
        print *, "Error: Invalid step value"
        call print_usage()
        stop
    end if
    
    ! Parameter 5: trajectory_filename
    call get_command_argument(5, trajectory_filename)
    
    ! Parameter 6: ephType
    call get_command_argument(6, arg)
    read(arg, *, iostat=status) ephType
    if (status /= 0) then
        print *, "Error: Invalid ephType value"
        call print_usage()
        stop
    end if
    
    ! Parameter 7: Cd
    call get_command_argument(7, arg)
    read(arg, *, iostat=status) Cd
    if (status /= 0) then
        print *, "Error: Invalid Cd value"
        call print_usage()
        stop
    end if
    
    ! Parameter 8: amrDrag
    call get_command_argument(8, arg)
    read(arg, *, iostat=status) amrDrag
    if (status /= 0) then
        print *, "Error: Invalid amrDrag value"
        call print_usage()
        stop
    end if
    
    ! Parameter 9: F107
    call get_command_argument(9, arg)
    read(arg, *, iostat=status) F107
    if (status /= 0) then
        print *, "Error: Invalid F107 value"
        call print_usage()
        stop
    end if
    
    ! Parameter 10: Ap
    call get_command_argument(10, arg)
    read(arg, *, iostat=status) Ap
    if (status /= 0) then
        print *, "Error: Invalid Ap value"
        call print_usage()
        stop
    end if
    
    ! Parameter 11: reflectivity
    call get_command_argument(11, arg)
    read(arg, *, iostat=status) reflectivity
    if (status /= 0) then
        print *, "Error: Invalid reflectivity value"
        call print_usage()
        stop
    end if
    
    ! Parameter 12: amrSRP
    call get_command_argument(12, arg)
    read(arg, *, iostat=status) amrSRP
    if (status /= 0) then
        print *, "Error: Invalid amrSRP value"
        call print_usage()
        stop
    end if
    
    ! Parameter 13: degree (引力场阶数 n)
    call get_command_argument(13, arg)
    read(arg, *, iostat=status) degree
    if (status /= 0) then
        print *, "Error: Invalid degree value"
        call print_usage()
        stop
    end if
    
    ! Parameter 14: order (引力场次数 m)
    call get_command_argument(14, arg)
    read(arg, *, iostat=status) order
    if (status /= 0) then
        print *, "Error: Invalid order value"
        call print_usage()
        stop
    end if
    
    call get_command_argument(15, arg)
    read(arg, *, iostat=status) K_lunar
    if (status /= 0) then
        print *, "Error: Invalid K_lunar value"
        call print_usage()
        stop
    end if
    
    call get_command_argument(16, arg)
    read(arg, *, iostat=status) K_lunarTide
    if (status /= 0) then
        print *, "Error: Invalid K_lunarTide value"
        call print_usage()
        stop
    end if
    
    call get_command_argument(17, arg)
    read(arg, *, iostat=status) K_solor
    if (status /= 0) then
        print *, "Error: Invalid K_solor value"
        call print_usage()
        stop
    end if
    
    call get_command_argument(18, arg)
    read(arg, *, iostat=status) K_solarTide
    if (status /= 0) then
        print *, "Error: Invalid K_solarTide value"
        call print_usage()
        stop
    end if
    
    call get_command_argument(19, arg)
    read(arg, *, iostat=status) K_SRP
    if (status /= 0) then
        print *, "Error: Invalid K_SRP value"
        call print_usage()
        stop
    end if
    
    call get_command_argument(20, arg)
    read(arg, *, iostat=status) K_drag
    if (status /= 0) then
        print *, "Error: Invalid K_drag value"
        call print_usage()
        stop
    end if
    
    call get_command_argument(21, arg)
    read(arg, *, iostat=status) K_PN
    if (status /= 0) then
        print *, "Error: Invalid K_PN value"
        call print_usage()
        stop
    end if
    
    ! Assign force model switches (argv 15-21, same order as before)
    sat%opMod%key_moonGravity = K_lunar
    sat%opMod%key_moonTide = K_lunarTide
    sat%opMod%key_sunGravity = K_solor
    sat%opMod%key_sunTide = K_solarTide
    sat%opMod%key_srp = K_SRP
    sat%opMod%key_drag = K_drag
    sat%opMod%key_postNewton = K_PN

    ! Calculate number of steps
    num_steps = ceiling(abs(timeMinus(end_time, start_time) * 86400.0d0 / step)) + 1
    
    ! Allocate and populate epochs array
    allocate(epochs(num_steps))
    epochs(1) = start_time
    do i = 2, num_steps
        epochs(i)%day = start_time%day
        epochs(i)%sec = start_time%sec + (i-1) * step
        ! Normalize if needed
        do while (epochs(i)%sec >= 86400.0d0)
            epochs(i)%day = epochs(i)%day + 1
            epochs(i)%sec = epochs(i)%sec - 86400.0d0
        end do
        epochs(i)%val = epochs(i)%day + epochs(i)%sec/86400.0d0
    end do

    ! Initialize satellite
    call kpl_SI2normal(kepler_elements, kepler_elements)  ! Convert to standard units
    call kpl2nsg1(kepler_elements, nsg1)
    call sat%init(start_time, nsg1, number = 1, name = name, eleType = eleType, propagator = propagator, &
                 Cd = Cd, amrDrag = amrDrag, F107 = F107, Ap = Ap, reflectivity = reflectivity, &
                amrSRP = amrSRP, degree = degree, order = order)
    
    ! Display satellite information
    call sat%show()
    
    ! Generate trajectory
    print *, "Generating trajectory..."
    call sat%generateTrajectory(epochs)
    
    ! Save trajectory to file
    print *, "Saving trajectory to file: ", trim(trajectory_filename)
    call sat%savetrajectory(ephType, trajectory_filename)
    
    print *, "Orbit prediction completed successfully."
    
    ! Clean up
    deallocate(epochs)

contains

    ! Function to parse date-time string in format YYYY MM DD HH MM SS.SSS
    subroutine parseDateTime(datetime_str, t)
        character(len=*), intent(in) :: datetime_str
        type(type_t), intent(out) :: t
        integer :: year, month, day, hour, minute
        real*8 :: second
        integer :: status
        
        ! Parse the date-time string
        read(datetime_str, '(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f6.3)', iostat=status) &
            year, month, day, hour, minute, second
            
        if (status /= 0) then
            print *, "Error parsing date-time: ", trim(datetime_str)
            print *, "Expected format: YYYY MM DD HH MM SS.SSS"
            stop
        end if
        
        call t%init(DBLE(year), DBLE(month), DBLE(day), DBLE(hour), DBLE(minute), second)
        
    end subroutine parseDateTime
    
    ! Function to parse Keplerian elements string (space-separated values)
    subroutine parse_elements(elements_str, elements)
        character(len=*), intent(in) :: elements_str
        real*8, intent(out) :: elements(6)
        integer :: i, pos1, pos2, status
        character(len=255) :: temp_str
        
        ! Parse the elements string
        pos1 = 1
        do i = 1, 6
            pos2 = index(elements_str(pos1:), ' ')
            if (pos2 == 0 .and. i < 6) then
                print *, "Error: Invalid Keplerian elements format. Expected 6 space-separated values."
                print *, "elements_str: ", elements_str
                stop
            end if
            if (pos2 == 0) then
                temp_str = elements_str(pos1:)
            else
                temp_str = elements_str(pos1:pos1+pos2-2)
            end if
            read(temp_str, *, iostat=status) elements(i)
            if (status /= 0) then
                print *, "Error: Invalid value in Keplerian elements at position", i
                stop
            end if
            pos1 = pos1 + pos2
        end do
    end subroutine parse_elements

    ! Print usage information
    subroutine print_usage()
        print *, "Usage: orbit_prediction <start_time> <kepler_elements>", &
        " <end_time> <step> <filename> <ephType> <Cd> <amrDrag> <F107> <Ap>", &
        " <reflectivity> <amrSRP> <degree> <order> <K_lunar> <K_lunarTide>", &
        " <K_solor> <K_solarTide> <K_SRP> <K_drag> <K_PN>"
        print *, ""
        print *, "Required parameters (21 total):"
        print *, "  1.  start_time      - Start time in format 'YYYY MM DD HH MM SS.SSS'"
        print *, "  2.  kepler_elements - Six Keplerian elements as space-separated values, km, deg"
        print *, "  3.  end_time        - End time in format 'YYYY MM DD HH MM SS.SSS'"
        print *, "  4.  step            - Time step in seconds"
        print *, "  5.  filename        - Output trajectory filename"
        print *, "  6.  ephType         - trajectory type (1: Keplerian elements, ", &
        "2: position/velocity, typically use 2)"
        print *, "  7.  Cd              - Drag coefficient (typically 2.2)"
        print *, "  8.  amrDrag         - Area-to-mass ratio for drag (typically 0.02)"
        print *, "  9.  F107            - Solar flux index (typically 130)"
        print *, "  10. Ap              - Geomagnetic index (typically 15)"
        print *, "  11. reflectivity    - Surface reflectivity (typically 0)"
        print *, "  12. amrSRP          - Area-to-mass ratio for solar radiation pressure ", &
        "(typically 0.02)"
        print *, "  13. degree          - Gravity field degree n (typically 2)"
        print *, "  14. order           - Gravity field order m (typically = degree; 0 = zonal only)"
        print *, "  15. K_lunar         - Lunar perturbation flag (typically 0)"
        print *, "  16. K_lunarTide     - Lunar tide flag (typically 0)"
        print *, "  17. K_solor         - Solar perturbation flag (typically 0)"
        print *, "  18. K_solarTide     - Solar tide flag (typically 0)"
        print *, "  19. K_SRP           - Solar radiation pressure flag (typically 0)"
        print *, "  20. K_drag          - Atmospheric drag flag (typically 0)"
        print *, "  21. K_PN            - Post-Newtonian flag (typically 0)"
        print *, ""
        print *, "Example:"
        print *, '  orbit_prediction "2023 01 01 12 00 00.000" "7000 0.01 45 30 60 90" ', &
        '"2023 01 02 12 00 00.000" 60 output.txt 2 2.2 0.02 130 15 0 0.02 2 2 0 0 0 0 0 0 0'
    end subroutine print_usage
    
    end program orbit_prediction