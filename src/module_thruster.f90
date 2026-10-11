module module_thruster

    use module_const
    use module_time
    use module_space

    implicit none

    !===========================================================================
    ! Thruster model with an epoch-driven ledger and a fuel budget
    !
    ! The counters are accrued LAZILY: the caller passes the current epoch as a
    ! type_t (UTC) whenever it switches the engine on/off or reads a total, and
    ! the ledger is flushed up to that epoch. Nothing inside the integrator
    ! accrues, so a burn is accounted as long as the caller flushes at the end
    ! of it.
    !
    ! Two ledgers are kept side by side:
    !
    !   current_*   the burn in progress: fuel, impulse and time since the
    !               engine was switched on. Zeroed when the engine goes off
    !               (deactivate / fault), so "off" always reads as zero.
    !   total_*     the lifetime budget: everything the thruster ever burned.
    !               Never reset by activate / deactivate; init() is the only
    !               reset (the lifetime counters are also what the force model
    !               reads for the mass).
    !
    ! Accounting rules:
    !   - activate(t)     : records t_on = t, starts a fresh burn (idempotent;
    !                       no-op if already on)
    !   - deactivate(t)   : flushes the ledger up to t, switches off, and zeroes
    !                       the current_* counters
    !   - update_state(t) : flush up to t
    !   - get_total_*(t)  : flush up to t, then return the lifetime total
    !
    ! Fuel: initial_fuel_mass is the propellant loaded at the start and
    ! current_fuel_mass what is left, decremented by each flush. The two are in
    ! kg (unlike the satellite mass, whose unit the caller owns). Because the
    ! ledger is lazy, both only move at a flush. initial_fuel_mass = 0 means the
    ! fuel budget is NOT modelled: nothing is depleted, get_thrust never flames
    ! out, and callers that do not care about fuel keep working unchanged.
    !
    ! The rate model is piecewise constant: a flush applies the CURRENT
    ! current_thrust over the whole [t_last, t] interval. Call update_state(t)
    ! before changing the thrust level (or pass t to set_current_thrust),
    ! otherwise the new level is applied retroactively.
    !
    ! All epochs must be UTC: timeMinus() compares day/sec without inspecting
    ! ttype, so mixing in a BDT/GPST epoch (whose day is week*7) would inject a
    ! spurious interval of tens of thousands of days.
    !===========================================================================

    ! Thruster type definition
    type :: type_thruster
        ! Basic properties
        integer         :: id = 0                    ! Thruster ID
        character*100   :: name = ""                 ! Thruster name
        double precision          :: efficiency = 1.0d0        ! Efficiency coefficient (0-1)
        double precision          :: response_time = 0d0       ! Response time (s); NOT implemented yet
        double precision          :: specific_impulse = 0d0    ! Specific impulse (s)
        double precision          :: min_thrust = 0d0          ! Minimum thrust (N)
        double precision          :: max_thrust = 1d0          ! Maximum thrust (N)

        ! Fuel budget (kg). initial_fuel_mass = 0: budget not modelled
        double precision          :: initial_fuel_mass = 0d0   ! Propellant loaded at the start (kg)
        double precision          :: current_fuel_mass = 0d0   ! Propellant still on board (kg)

        ! status
        logical         :: is_active = .false.       ! On/off status
        double precision          :: current_thrust = 0d0      ! Commanded thrust level (N), efficiency applied
        double precision          :: thrust_direction(3) = [0d0, 0d0, 0d0]  ! Thrust direction unit vector

        ! Ledger of the burn in progress; zeroed when the engine is off
        double precision          :: current_operation_time = 0d0  ! Time on since activation (s)
        double precision          :: current_impulse = 0d0         ! Impulse delivered since activation (N·s)
        double precision          :: current_fuel_used = 0d0       ! Fuel burned since activation (kg)

        ! Lifetime ledger; only init() resets these
        double precision          :: total_operation_time = 0d0    ! Cumulative operating time (s)
        double precision          :: total_impulse = 0d0           ! Total impulse delivered (N·s)
        double precision          :: total_fuel_used = 0d0         ! Fuel consumed (kg)

        ! Ledger epochs (UTC). Only meaningful once clock_started is true: an
        ! unset type_t defaults to MJD 0, which is indistinguishable from a
        ! legitimate epoch, hence the separate flag.
        type(type_t)    :: t_on                       ! Epoch of the last activation
        type(type_t)    :: t_last                     ! Epoch the ledger has been accrued up to
        logical         :: clock_started = .false.    ! Have t_on / t_last been assigned?

        ! Thrust direction related properties
        integer         :: coordinate_frame = 1      ! Coordinate system type
        ! 1 - J2000 (Earth-centered inertial)
        ! 2 - SBCS (Spacecraft body coordinate system)
        ! 3 - NTH (Orbital coordinate system: N-normal, T-velocity, H-angular momentum)

        ! Fault detection
        logical         :: is_faulty = .false.       ! Fault status
        character*200   :: fault_description = ""    ! Fault description

    contains
        ! Initialization methods
        procedure :: init => thrust_init
        procedure :: set_thrust_direction => set_thrust_direction
        procedure :: set_current_thrust => set_current_thrust

        ! Control methods (all take a UTC epoch)
        procedure :: activate => thrust_activate
        procedure :: deactivate => thrust_deactivate
        procedure :: toggle => thrust_toggle

        ! Getter methods
        procedure :: get_thrust => get_thrust
        procedure :: get_fuel_flow_rate => get_fuel_flow_rate
        procedure :: get_total_impulse => get_total_impulse
        procedure :: get_total_fuel_used => get_total_fuel_used
        procedure :: get_thrust_ECI => get_thrust_ECI

        ! Status update methods (take a UTC epoch)
        procedure :: update_state => update_thrust_state

        ! Fault handling
        procedure :: set_fault => set_thrust_fault
        procedure :: clear_fault => clear_thrust_fault

        ! Information display
        procedure :: show_info => show_thrust_info
    end type type_thruster

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize thruster
    !
    ! init() is a full reset: it clears both ledgers and reloads the tank.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine thrust_init(this, id, name, specific_impulse, max_thrust, &
                          coordinate_frame, thrust_direction, efficiency, &
                          min_thrust, response_time, initial_fuel_mass)
    class(type_thruster), intent(inout) :: this
    integer, optional, intent(in) :: id
    character(*), optional, intent(in) :: name
    double precision, optional, intent(in) :: specific_impulse
    double precision, optional, intent(in) :: max_thrust
    integer, optional, intent(in) :: coordinate_frame
    double precision, optional, intent(in) :: thrust_direction(3)
    double precision, optional, intent(in) :: efficiency
    double precision, optional, intent(in) :: min_thrust
    double precision, optional, intent(in) :: response_time
    double precision, optional, intent(in) :: initial_fuel_mass  ! kg; 0 = budget not modelled

    ! Set basic parameters
    if (present(id)) this%id = id
    if (present(name)) this%name = trim(name)
    if (present(specific_impulse)) this%specific_impulse = specific_impulse
    if (present(max_thrust)) this%max_thrust = max_thrust
    if (present(coordinate_frame)) this%coordinate_frame = coordinate_frame
    if (present(efficiency)) this%efficiency = efficiency
    if (present(min_thrust)) this%min_thrust = min_thrust
    if (present(response_time)) this%response_time = response_time
    if (present(initial_fuel_mass)) this%initial_fuel_mass = initial_fuel_mass

    ! Set thrust direction
    if (present(thrust_direction)) then
        this%thrust_direction = thrust_direction
        ! Normalize thrust direction vector
        call normalize_vector(this%thrust_direction)
    else
        ! Default thrust direction is Z-axis
        this%thrust_direction = [0d0, 0d0, 1d0]
    endif

    ! Initialize other parameters
    this%current_thrust = 0d0
    this%is_active = .false.
    this%is_faulty = .false.
    this%fault_description = ""

    ! Tank is full, both ledgers empty
    this%current_fuel_mass = this%initial_fuel_mass
    this%current_operation_time = 0d0
    this%current_impulse = 0d0
    this%current_fuel_used = 0d0
    this%total_operation_time = 0d0
    this%total_impulse = 0d0
    this%total_fuel_used = 0d0

    this%clock_started = .false.   ! no valid epoch yet; set by the first activate()

    end subroutine thrust_init

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set thrust direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine set_thrust_direction(this, direction)
    class(type_thruster), intent(inout) :: this
    double precision, intent(in) :: direction(3)

    if (this%is_faulty) then
        write(*, *) "Warning: Thruster is faulty, cannot set thrust direction"
        return
    endif

    this%thrust_direction = direction
    call normalize_vector(this%thrust_direction)

    end subroutine set_thrust_direction

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set thrust magnitude
    !
    ! The stored current_thrust is the commanded magnitude clamped to
    ! [min_thrust, max_thrust] and scaled by efficiency. The fuel flow rate is
    ! NOT stored any more: it follows from this level and the specific impulse
    ! (see fuel_flow_rate), so changing the level is enough to change the rate.
    !
    ! Because a flush applies the CURRENT rate, pass t (UTC) to settle the
    ! elapsed interval at the OLD rate before the level changes; otherwise the
    ! new level is applied retroactively to everything since t_last.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine set_current_thrust(this, magnitude, t)
    class(type_thruster), intent(inout) :: this
    double precision, intent(in) :: magnitude
    type(type_t), optional, intent(in) :: t   ! UTC epoch; if given, flush before changing level

    if (this%is_faulty) then
        write(*, *) "Warning: Thruster is faulty, cannot set thrust magnitude"
        return
    endif

    ! Settle the elapsed interval at the old rate before overwriting it
    if (present(t)) call this%update_state(t)

    ! Check if thrust is within allowed range
    if (magnitude < this%min_thrust) then
        write(*, *) "Warning: commanded thrust ", magnitude, &
                    " N is below min_thrust; clamped to ", this%min_thrust, " N"
        this%current_thrust = this%min_thrust * this%efficiency
    else if (magnitude > this%max_thrust) then
        write(*, *) "Warning: commanded thrust ", magnitude, &
                    " N is above max_thrust; clamped to ", this%max_thrust, " N"
        this%current_thrust = this%max_thrust * this%efficiency
    else
        this%current_thrust = magnitude * this%efficiency
    endif

    end subroutine set_current_thrust

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Activate thruster
    !
    ! Records the turn-on epoch and opens a fresh burn ledger. The lifetime
    ! totals are untouched -- switching off and on again continues them.
    ! Calling it on an already-running thruster is a no-op, so repeated
    ! activate() calls never wipe an in-progress burn.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine thrust_activate(this, t)
    class(type_thruster), intent(inout) :: this
    type(type_t), intent(in) :: t   ! UTC epoch of turn-on

    ! The fault test must come first: a faulty-but-active thruster must not be
    ! silently "re-activated" into staying on.
    if (this%is_faulty) then
        write(*, *) "Warning: Thruster is faulty, cannot activate"
        return
    endif

    ! Already running: idempotent, keep the current burn
    if (this%is_active) return

    ! Open a fresh burn. The counters should already be zero (deactivate and
    ! set_fault both clear them), but a burn cut short without going through
    ! either of those must not leak into the new one.
    this%current_operation_time = 0d0
    this%current_impulse = 0d0
    this%current_fuel_used = 0d0

    this%t_on = t
    this%t_last = t                ! the burn is measured from here
    this%clock_started = .true.
    this%is_active = .true.

    end subroutine thrust_activate

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deactivate thruster
    !
    ! Flushes the burn up to t, switches off, and zeroes the current_* counters:
    ! an engine that is off has burnt nothing "in the current burn", and the
    ! totals keep the history. Switching off is always permitted, including for
    ! a faulty thruster.
    !
    ! current_thrust is deliberately NOT zeroed: the off state is expressed by
    ! is_active (get_thrust / get_fuel_flow_rate both test it), and wiping the
    ! commanded level would make the next activation silently produce no thrust.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine thrust_deactivate(this, t)
    class(type_thruster), intent(inout) :: this
    type(type_t), intent(in) :: t   ! UTC epoch of switch-off

    ! Nothing was running: do not touch the configured thrust level
    if (.not. this%is_active) return

    call this%update_state(t)
    this%is_active = .false.

    this%current_operation_time = 0d0
    this%current_impulse = 0d0
    this%current_fuel_used = 0d0

    end subroutine thrust_deactivate

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Toggle thruster status
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine thrust_toggle(this, t)
    class(type_thruster), intent(inout) :: this
    type(type_t), intent(in) :: t   ! UTC epoch of the switch

    if (this%is_active) then
        call this%deactivate(t)
    else
        call this%activate(t)
    endif

    end subroutine thrust_toggle

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get thrust
    !
    ! Kept argument-free and intent(in): this is called from the RKF78 right-hand
    ! side at every stage, so it must stay free of ledger bookkeeping and I/O.
    ! A thruster with a modelled fuel budget stops thrusting once the tank is dry
    ! (the ledger itself is only settled by a flush, so the flame-out is seen
    ! here no earlier than the flush that emptied it).
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function get_thrust(this) result(thrust_vector)
    class(type_thruster), intent(in) :: this
    double precision :: thrust_vector(3)

    if (this%initial_fuel_mass > 0d0 .and. this%current_fuel_mass <= 0d0) then
        thrust_vector = [0d0, 0d0, 0d0]     ! flame-out
    elseif (this%is_active .and. .not. this%is_faulty) then
        thrust_vector = this%thrust_direction * this%current_thrust
    else
        thrust_vector = [0d0, 0d0, 0d0]
    endif

    end function get_thrust

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get fuel flow rate
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function get_fuel_flow_rate(this) result(flow_rate)
    class(type_thruster), intent(in) :: this
    double precision :: flow_rate

    if (this%is_active .and. .not. this%is_faulty) then
        flow_rate = fuel_flow_rate(this%current_thrust, this%specific_impulse)
    else
        flow_rate = 0d0
    endif

    end function get_fuel_flow_rate

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get total impulse accrued up to epoch t (N·s, whole life of the thruster)
    !
    ! intent(inout): the ledger is flushed to t before reading, so this cannot be
    ! called on an intent(in) dummy.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function get_total_impulse(this, t) result(total_impulse)
    class(type_thruster), intent(inout) :: this
    type(type_t), intent(in) :: t   ! UTC epoch to account up to
    double precision :: total_impulse

    call this%update_state(t)
    total_impulse = this%total_impulse

    end function get_total_impulse

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get total fuel consumed up to epoch t (kg, whole life of the thruster)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function get_total_fuel_used(this, t) result(total_fuel_used)
    class(type_thruster), intent(inout) :: this
    type(type_t), intent(in) :: t   ! UTC epoch to account up to
    double precision :: total_fuel_used

    call this%update_state(t)
    total_fuel_used = this%total_fuel_used

    end function get_total_fuel_used

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get thrust in ECI coordinate system
    !
    ! Kept argument-free and intent(in): this is called from the RKF78 right-hand
    ! side at every stage, so it must stay free of ledger bookkeeping and I/O.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function get_thrust_ECI(this, rs, vs) result(thrust_ECI)
    class(type_thruster), intent(in) :: this
    double precision, intent(in) :: rs(3)  ! Satellite position in ECI (km)
    double precision, intent(in) :: vs(3)  ! Satellite velocity in ECI (km/s)
    double precision :: thrust_ECI(3)      ! Thrust vector in ECI (N)
    double precision :: thrust_local(3)    ! Thrust vector in local coordinate system
    double precision :: M(3, 3)           ! Transformation matrix

    ! Get thrust in local coordinate system
    thrust_local = this%get_thrust()

    ! Convert thrust to ECI coordinate system based on coordinate frame
    select case(this%coordinate_frame)
    case(1)  ! J2000 (already in ECI)
        thrust_ECI = thrust_local

    case(2)  ! SBCS (Spacecraft body coordinate system)
        ! SBCS to ECI conversion requires attitude information
        ! Currently not implemented - return zero thrust
        write(*, *) "Warning: SBCS to ECI conversion not implemented"
        thrust_ECI = [0d0, 0d0, 0d0]

    case(3)  ! NTH (Orbital coordinate system)
        ! Get transformation matrix from NTH to ECI
        call NTH2ECI(rs, vs, M)
        ! Transform thrust vector: thrust_ECI = M * thrust_local
        thrust_ECI = matmul(M, thrust_local)

    case default
        write(*, *) "Warning: Unknown coordinate frame, returning zero thrust"
        thrust_ECI = [0d0, 0d0, 0d0]
    end select

    end function get_thrust_ECI

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Update thruster state: accrue the interval [t_last, t] onto both ledgers
    !
    ! INPUTS:
    !   t - UTC epoch to account up to
    !
    ! Notes:
    !   - Does nothing while the thruster is off, so reading a total after a burn
    !     never picks up the idle time in between.
    !   - t_last is advanced for every accounted interval even when the thruster
    !     is faulty or out of fuel (a dead engine burns no fuel, but the pointer
    !     must move or the span would be accrued after clear_fault / refuelling).
    !   - A non-positive interval is ignored without moving t_last backwards, so
    !     repeated flushes at the same epoch and out-of-order epochs cannot
    !     double-count. Out-of-order epochs are reported, since silently dropping
    !     ledger time is worse than double counting.
    !   - A burn that runs the tank dry inside the interval is truncated at the
    !     flame-out epoch, so time / impulse / fuel stay consistent with each
    !     other and mass never drops below the dry mass.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine update_thrust_state(this, t)
    class(type_thruster), intent(inout) :: this
    type(type_t), intent(in) :: t   ! UTC epoch to account up to
    double precision :: dt          ! Elapsed time since the last flush (s)
    double precision :: mdot        ! Fuel flow rate over the interval (kg/s)
    double precision :: dI, dm      ! Impulse (N·s) and fuel (kg) accrued by this flush

    ! Nothing to account for while the engine is off; leave the ledger untouched
    if (.not. this%is_active) return
    if (.not. this%clock_started) return

    if (t%ttype /= TIME_TYPE_UTC) then
        write(*, *) "Warning: update_state expects a UTC epoch, got ttype = ", t%ttype
    endif

    dt = timeMinus(t, this%t_last) * 86400d0   ! timeMinus returns days

    if (dt < 0d0) then
        write(*, *) "Warning: epoch is earlier than the last accounted epoch; " // &
                    "interval ignored (dt = ", dt, " s)"
    endif
    if (dt <= 0d0) return

    ! Advance the accounting pointer before the fault / fuel tests: a dead span
    ! must be consumed, not deferred.
    this%t_last = t
    if (this%is_faulty) return

    ! An empty tank delivers nothing (get_thrust agrees): consume the interval
    ! without accruing. Only the modelled-budget case can get here.
    if (this%initial_fuel_mass > 0d0 .and. this%current_fuel_mass <= 0d0) return

    ! Burn only as long as the tank allows. The `>` below keeps the warning to
    ! one per flame-out: an already-empty tank returned above, and the rate is
    ! positive here (current_fuel_mass > 0, so mdot*dt > it forces mdot > 0).
    mdot = fuel_flow_rate(this%current_thrust, this%specific_impulse)
    if (this%initial_fuel_mass > 0d0 .and. mdot*dt > this%current_fuel_mass) then
        dt = this%current_fuel_mass / mdot
        write(*, *) "Warning: thruster fuel exhausted; burn truncated to ", dt, " s"
    endif

    dI = this%current_thrust * dt
    dm = mdot * dt

    ! The burn in progress ...
    this%current_operation_time = this%current_operation_time + dt
    this%current_impulse = this%current_impulse + dI
    this%current_fuel_used = this%current_fuel_used + dm

    ! ... and the lifetime totals
    this%total_operation_time = this%total_operation_time + dt
    this%total_impulse = this%total_impulse + dI
    this%total_fuel_used = this%total_fuel_used + dm

    ! Tanks full of rounding error are still empty
    if (this%initial_fuel_mass > 0d0) then
        this%current_fuel_mass = max(0d0, this%current_fuel_mass - dm)
    endif

    end subroutine update_thrust_state

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set fault status
    !
    ! Pass t to settle the ledger before the fault is raised; without it the
    ! interval since the last flush cannot be accounted and is lost.
    ! The engine is cut, so the burn counters are cleared like in deactivate().
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine set_thrust_fault(this, description, t)
    class(type_thruster), intent(inout) :: this
    character(*), intent(in) :: description
    type(type_t), optional, intent(in) :: t   ! UTC epoch of the fault

    ! Flush while the thruster is still healthy, then flag the fault
    if (present(t)) then
        call this%deactivate(t)
    endif

    this%is_faulty = .true.
    this%fault_description = trim(description)
    this%is_active = .false.   ! engine cut; t_last stays frozen if no epoch was given

    this%current_operation_time = 0d0
    this%current_impulse = 0d0
    this%current_fuel_used = 0d0

    end subroutine set_thrust_fault

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Clear fault status
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine clear_thrust_fault(this)
    class(type_thruster), intent(inout) :: this

    this%is_faulty = .false.
    this%fault_description = ""

    end subroutine clear_thrust_fault

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Show thruster information
    !
    ! The current_* counters describe the burn in progress (zero while the engine
    ! is off, even just after a long burn); the total_* counters are the lifetime
    ! budget, and total_fuel_used / current_fuel_mass are what the force model
    ! uses to keep the mass up to date.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine show_thrust_info(this)
    class(type_thruster), intent(in) :: this
    character*20 :: coord_name
    character*28 :: on_string, last_string
    double precision :: thrust_vector(3)

    ! Determine coordinate system name
    select case(this%coordinate_frame)
    case(1)
        coord_name = "J2000"
    case(2)
        coord_name = "SBCS"
    case(3)
        coord_name = "NTH"
    case default
        coord_name = "Unknown"
    end select

    ! Format the ledger epochs; a never-activated thruster has no valid epoch
    if (this%clock_started) then
        on_string = this%t_on%toDate()
        last_string = this%t_last%toDate()
    else
        on_string = "<unset>"
        last_string = "<unset>"
    endif

    ! Get current thrust vector
    thrust_vector = this%get_thrust()

    print*, "=============================================================="
    write(*, *) "Thruster Information:"
    write(*, *) "ID: ", this%id
    write(*, *) "Name: ", trim(this%name)
    if (this%is_active) then
        write(*, *) "Status: Active"
    else
        write(*, *) "Status: Inactive"
    endif

    if (this%is_faulty) then
        write(*, *) "Fault: Yes"
    else
        write(*, *) "Fault: No"
    endif
    if (this%is_faulty) then
        write(*, *) "Fault Description: ", trim(this%fault_description)
    endif
    write(*, *) "Specific Impulse: ", this%specific_impulse, " s"
    write(*, *) "Maximum Thrust: ", this%max_thrust, " N"
    write(*, *) "Current Thrust: ", this%current_thrust, " N"
    write(*, *) "Coordinate System: ", trim(coord_name)
    write(*, "(A,3F10.6)") " Thrust Direction: ", this%thrust_direction
    write(*, "(A,3F10.3)") " Thrust Vector: ", thrust_vector
    write(*, *) "Efficiency: ", this%efficiency
    write(*, *) "Fuel Flow Rate: ", this%get_fuel_flow_rate(), " kg/s"
    print*, "--------------------------------------------------------------"
    write(*, *) "Initial Fuel Mass: ", this%initial_fuel_mass, " kg"
    write(*, *) "Current Fuel Mass: ", this%current_fuel_mass, " kg"
    print*, "--------------------------------------------------------------"
    write(*, *) "Turn-on Epoch: ", trim(on_string)
    write(*, *) "Ledger Settled To: ", trim(last_string)
    write(*, *) "Current Operating Time: ", this%current_operation_time, " s"
    write(*, *) "Current Impulse: ", this%current_impulse, " N·s"
    write(*, *) "Current Fuel Used: ", this%current_fuel_used, " kg"
    print*, "--------------------------------------------------------------"
    write(*, *) "Total Operating Time: ", this%total_operation_time, " s"
    write(*, *) "Total Impulse: ", this%total_impulse, " N·s"
    write(*, *) "Total Fuel Used: ", this%total_fuel_used, " kg"
    print*, "=============================================================="

    end subroutine show_thrust_info

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Helper function: Normalize vector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine normalize_vector(vec)
    double precision, intent(inout) :: vec(3)
    double precision :: norm

    norm = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
    if (norm > 0d0) then
        vec = vec / norm
    else
        vec = [0d0, 0d0, 1d0]  ! Default Z-axis direction
    endif

    end subroutine normalize_vector

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Helper function: fuel flow rate of a given (effective) thrust level, kg/s
    !
    ! mdot = F / (Isp * g0). acc_unit plays the role of g0 here, as it did when
    ! the rate was stored on set_current_thrust -- the module has always measured
    ! Isp against the normalized acceleration unit rather than 9.80665 m/s^2.
    ! The rate is derived on demand instead of being cached, so it can never go
    ! stale against a thrust level changed behind its back.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pure function fuel_flow_rate(thrust, isp) result(mdot)
    double precision, intent(in) :: thrust   ! Effective thrust (N)
    double precision, intent(in) :: isp      ! Specific impulse (s)
    double precision :: mdot                 ! Flow rate (kg/s)

    if (isp > 0d0) then
        mdot = thrust / (isp * acc_unit)
    else
        mdot = 0d0
    endif

    end function fuel_flow_rate

end module module_thruster
