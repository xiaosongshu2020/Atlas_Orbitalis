!==============================================================================
! test_thruster.f90
!
! Purpose: verify the epoch-driven ledger and the fuel budget of module_thruster
!          (type_thruster).
!
! The thruster accrues its counters lazily: activate / deactivate / update_state
! / get_total_* all take a UTC epoch and settle the interval [t_last, t]. This
! test drives that ledger with hand-made epochs and checks the arithmetic; it
! does not touch the integrator.
!
! Two ledgers are checked side by side: current_* covers the burn in progress
! and is zeroed when the engine goes off, total_* covers the whole life of the
! thruster and is only cleared by init().
!
! No thruster of the reference orbit is involved, so the test is independent of
! the force model and of the (currently broken) OpenMP build.
!
! Reference values used throughout:
!   thrust command   0.05 N, efficiency 0.95  ->  current_thrust = 0.0475 N
!   specific impulse 3000 s                   ->  mdot = 0.0475 / (3000 * ACC_UNIT)
!   one hour of burn  ->  171.0 N*s, 5.817e-3 kg of fuel
!==============================================================================
program test_thruster

    use module_const, only: ACC_UNIT
    use module_time
    use module_thruster

    implicit none

    integer :: nFail
    double precision :: thr_expected, mdot_expected

    nFail = 0

    thr_expected = 0.05d0 * 0.95d0                  ! commanded N after efficiency
    mdot_expected = thr_expected / (3000d0 * ACC_UNIT)

    print *, "=============================================================="
    print *, " test_thruster - epoch-driven ledger of type_thruster"
    print *, "--------------------------------------------------------------"
    print *, " effective thrust   : ", thr_expected, " N"
    print *, " fuel flow rate     : ", mdot_expected, " kg/s"
    print *, " one hour of burn   : ", thr_expected * 3600d0, " N*s"
    print *, " one hour of fuel   : ", mdot_expected * 3600d0, " kg"
    print *, "=============================================================="

    call test_continuous_burn(thr_expected, mdot_expected, nFail)
    call test_repeat_flush(thr_expected, nFail)
    call test_out_of_order_epoch(thr_expected, nFail)
    call test_idle_time_not_counted(thr_expected, mdot_expected, nFail)
    call test_deactivate_zeroes_current(thr_expected, mdot_expected, nFail)
    call test_totals_accumulate_across_burns(thr_expected, mdot_expected, nFail)
    call test_fault(thr_expected, mdot_expected, nFail)
    call test_faulted_span_not_deferred(thr_expected, nFail)
    call test_rate_change_mid_burn(thr_expected, nFail)
    call test_deactivate_without_activation(thr_expected, nFail)
    call test_toggle(thr_expected, mdot_expected, nFail)
    call test_fuel_budget(mdot_expected, nFail)
    call test_flame_out(thr_expected, mdot_expected, nFail)
    call test_unlimited_fuel(thr_expected, nFail)

    print *, "--------------------------------------------------------------"
    if (nFail == 0) then
        print *, " test_thruster: ALL CHECKS PASSED"
    else
        print *, " test_thruster: ", nFail, " CHECK(S) FAILED"
    endif
    print *, "=============================================================="

    if (nFail > 0) stop 1

contains

    !--------------------------------------------------------------------------
    ! A thruster configured like the prediction_thruster demo.
    ! fuel < 0 (the default) leaves the fuel budget unmodelled.
    !--------------------------------------------------------------------------
    function new_thruster(fuel) result(thr)
    double precision, optional, intent(in) :: fuel
    type(type_thruster) :: thr

    if (present(fuel)) then
        call thr%init(id=1, name="Test_Thruster", &
                      specific_impulse=3000d0, &
                      max_thrust=0.1d0, &
                      coordinate_frame=3, &
                      thrust_direction=[0.0d0, 1.0d0, 0.0d0], &
                      efficiency=0.95d0, &
                      min_thrust=0.01d0, &
                      initial_fuel_mass=fuel)
    else
        call thr%init(id=1, name="Test_Thruster", &
                      specific_impulse=3000d0, &
                      max_thrust=0.1d0, &
                      coordinate_frame=3, &
                      thrust_direction=[0.0d0, 1.0d0, 0.0d0], &
                      efficiency=0.95d0, &
                      min_thrust=0.01d0)
    endif
    call thr%set_current_thrust(0.05d0)   ! 0.0475 N after efficiency

    end function new_thruster

    !--------------------------------------------------------------------------
    ! UTC epoch on 2025-01-01
    !--------------------------------------------------------------------------
    subroutine mkTime(t, hh, mm)
    type(type_t), intent(out) :: t
    integer, intent(in) :: hh, mm

    call t%init(2025d0, 1d0, 1d0, dble(hh), dble(mm), 0d0)

    end subroutine mkTime

    !--------------------------------------------------------------------------
    ! Scalar check
    !--------------------------------------------------------------------------
    subroutine check(name, got, expect, tol, nFail)
    character(*), intent(in) :: name
    double precision, intent(in) :: got, expect, tol
    integer, intent(inout) :: nFail
    double precision :: err

    err = abs(got - expect)
    if (err <= tol) then
        write(*, "(A,A,A,F16.6)") "  PASS  ", name, "   = ", got
    else
        nFail = nFail + 1
        write(*, "(A,A)") "  FAIL  ", name
        write(*, "(A,F16.6,A,F16.6,A,ES10.2)") "          got ", got, &
            "  expected ", expect, "  err ", err
    endif

    end subroutine check

    !--------------------------------------------------------------------------
    ! A single uninterrupted burn accrues linearly with time, on both ledgers
    !--------------------------------------------------------------------------
    subroutine test_continuous_burn(thrVal, mdot, nFail)
    double precision, intent(in) :: thrVal, mdot
    integer, intent(inout) :: nFail
    type(type_thruster) :: thr
    type(type_t) :: t0, t1, t2

    print *, ""
    print *, "[1] continuous burn over two hours"

    thr = new_thruster()
    call mkTime(t0, 0, 0)
    call mkTime(t1, 1, 0)
    call mkTime(t2, 2, 0)

    call thr%activate(t0)
    call check("total impulse after 1 h", thr%get_total_impulse(t1), thrVal * 3600d0, 1d-9, nFail)
    call check("total impulse after 2 h", thr%get_total_impulse(t2), thrVal * 7200d0, 1d-9, nFail)
    call check("total fuel after 2 h", thr%get_total_fuel_used(t2), mdot * 7200d0, 1d-12, nFail)
    call check("total operating time after 2 h", thr%total_operation_time, 7200d0, 1d-9, nFail)
    ! The same burn, seen from the current_* side (engine still on)
    call check("current impulse after 2 h", thr%current_impulse, thrVal * 7200d0, 1d-9, nFail)
    call check("current fuel after 2 h", thr%current_fuel_used, mdot * 7200d0, 1d-12, nFail)
    call check("current operating time after 2 h", thr%current_operation_time, 7200d0, 1d-9, nFail)

    end subroutine test_continuous_burn

    !--------------------------------------------------------------------------
    ! Flushing twice at the same epoch must not double count
    !--------------------------------------------------------------------------
    subroutine test_repeat_flush(thrVal, nFail)
    double precision, intent(in) :: thrVal
    integer, intent(inout) :: nFail
    type(type_thruster) :: thr
    type(type_t) :: t0, t1
    double precision :: first, second

    print *, ""
    print *, "[2] repeated flush at the same epoch"

    thr = new_thruster()
    call mkTime(t0, 0, 0)
    call mkTime(t1, 1, 0)

    call thr%activate(t0)
    first = thr%get_total_impulse(t1)
    second = thr%get_total_impulse(t1)
    call check("second read is not larger", second, first, 1d-12, nFail)
    call check("value is still one hour", second, thrVal * 3600d0, 1d-9, nFail)

    end subroutine test_repeat_flush

    !--------------------------------------------------------------------------
    ! An epoch before t_last is ignored and leaves the ledger intact
    !--------------------------------------------------------------------------
    subroutine test_out_of_order_epoch(thrVal, nFail)
    double precision, intent(in) :: thrVal
    integer, intent(inout) :: nFail
    type(type_thruster) :: thr
    type(type_t) :: t0, t1, t2, t3
    double precision :: afterReverse

    print *, ""
    print *, "[3] out-of-order epoch is ignored"

    thr = new_thruster()
    call mkTime(t0, 0, 0)
    call mkTime(t1, 1, 0)     ! 1 h
    call mkTime(t2, 2, 0)     ! 2 h
    call mkTime(t3, 3, 0)     ! 3 h

    call thr%activate(t0)
    afterReverse = thr%get_total_impulse(t2)
    ! Rewinding below t_last must be a no-op ...
    call check("rewind to 1 h ignored", thr%get_total_impulse(t1), afterReverse, 1d-12, nFail)
    call thr%update_state(t0)   ! far in the past; must not corrupt t_last
    ! ... so the next forward flush still covers exactly [2 h, 3 h]
    call check("3 h total after rewind", thr%get_total_impulse(t3), thrVal * 10800d0, 1d-9, nFail)

    end subroutine test_out_of_order_epoch

    !--------------------------------------------------------------------------
    ! Idle time between burns is not accrued
    !--------------------------------------------------------------------------
    subroutine test_idle_time_not_counted(thrVal, mdot, nFail)
    double precision, intent(in) :: thrVal, mdot
    integer, intent(inout) :: nFail
    type(type_thruster) :: thr
    type(type_t) :: t0, t1, t5

    print *, ""
    print *, "[4] idle time after switch-off is not counted"

    thr = new_thruster()
    call mkTime(t0, 0, 0)
    call mkTime(t1, 1, 0)
    call mkTime(t5, 5, 0)

    call thr%activate(t0)
    call thr%deactivate(t1)
    call check("burn of 1 h", thr%get_total_impulse(t1), thrVal * 3600d0, 1d-9, nFail)
    call check("fuel of 1 h", thr%get_total_fuel_used(t5), mdot * 3600d0, 1d-12, nFail)
    call check("4 h later, unchanged", thr%get_total_impulse(t5), thrVal * 3600d0, 1d-9, nFail)
    call check("total operating time still 1 h", thr%total_operation_time, 3600d0, 1d-9, nFail)

    end subroutine test_idle_time_not_counted

    !--------------------------------------------------------------------------
    ! Switching the engine off zeroes the current_* counters but keeps the totals
    !--------------------------------------------------------------------------
    subroutine test_deactivate_zeroes_current(thrVal, mdot, nFail)
    double precision, intent(in) :: thrVal, mdot
    integer, intent(inout) :: nFail
    type(type_thruster) :: thr
    type(type_t) :: t0, t1

    print *, ""
    print *, "[5] deactivate zeroes the current_* ledger only"

    thr = new_thruster()
    call mkTime(t0, 0, 0)
    call mkTime(t1, 1, 0)

    call thr%activate(t0)
    call thr%deactivate(t1)

    call check("current operating time cleared", thr%current_operation_time, 0d0, 1d-12, nFail)
    call check("current impulse cleared", thr%current_impulse, 0d0, 1d-12, nFail)
    call check("current fuel cleared", thr%current_fuel_used, 0d0, 1d-12, nFail)
    ! ... while the finished burn survives in the totals
    call check("total operating time kept", thr%total_operation_time, 3600d0, 1d-9, nFail)
    call check("total impulse kept", thr%total_impulse, thrVal * 3600d0, 1d-9, nFail)
    call check("total fuel kept", thr%total_fuel_used, mdot * 3600d0, 1d-12, nFail)

    end subroutine test_deactivate_zeroes_current

    !--------------------------------------------------------------------------
    ! The totals span the life of the thruster: a second burn adds to them, while
    ! current_* describes only the new burn
    !--------------------------------------------------------------------------
    subroutine test_totals_accumulate_across_burns(thrVal, mdot, nFail)
    double precision, intent(in) :: thrVal, mdot
    integer, intent(inout) :: nFail
    type(type_thruster) :: thr
    type(type_t) :: t0, t1, t2, t3, t4

    print *, ""
    print *, "[6] totals accumulate across burns"

    thr = new_thruster()
    call mkTime(t0, 0, 0)
    call mkTime(t1, 1, 0)
    call mkTime(t2, 2, 0)
    call mkTime(t3, 3, 0)
    call mkTime(t4, 4, 0)

    call thr%activate(t0)
    call thr%deactivate(t1)
    call check("first burn total", thr%get_total_impulse(t2), thrVal * 3600d0, 1d-9, nFail)

    call thr%activate(t2)                   ! second burn of 1 h, same lifetime
    call check("current cleared on re-activation", thr%current_impulse, 0d0, 1d-12, nFail)
    ! Calling activate again while already running must NOT clear the burn
    call thr%activate(t3)
    call check("re-activate while running is a no-op", &
               thr%get_total_impulse(t3), thrVal * 7200d0, 1d-9, nFail)
    call check("second burn current", thr%current_impulse, thrVal * 3600d0, 1d-9, nFail)
    call thr%deactivate(t3)                 ! close the second burn

    call check("lifetime total after two burns", thr%get_total_impulse(t4), &
               thrVal * 7200d0, 1d-9, nFail)   ! and frozen while off
    call check("lifetime fuel after two burns", thr%total_fuel_used, mdot * 7200d0, 1d-12, nFail)
    call check("lifetime operating time", thr%total_operation_time, 7200d0, 1d-9, nFail)

    end subroutine test_totals_accumulate_across_burns

    !--------------------------------------------------------------------------
    ! A fault settles the ledger up to the fault epoch, then cuts the engine
    !--------------------------------------------------------------------------
    subroutine test_fault(thrVal, mdot, nFail)
    double precision, intent(in) :: thrVal, mdot
    integer, intent(inout) :: nFail
    type(type_thruster) :: thr
    type(type_t) :: t0, t1, t5

    print *, ""
    print *, "[7] fault with an epoch"

    thr = new_thruster()
    call mkTime(t0, 0, 0)
    call mkTime(t1, 1, 0)
    call mkTime(t5, 5, 0)

    call thr%activate(t0)
    call thr%set_fault("test fault", t1)
    call check("hour before fault retained", thr%get_total_impulse(t1), thrVal * 3600d0, 1d-9, nFail)
    call check("engine cut by fault", merge(1d0, 0d0, thr%is_active), 0d0, 1d-12, nFail)
    call check("faulted span not accrued", thr%get_total_impulse(t5), thrVal * 3600d0, 1d-9, nFail)
    call check("fuel before fault retained", thr%get_total_fuel_used(t5), mdot * 3600d0, 1d-12, nFail)
    ! A fault cuts the engine, so it also closes the burn ledger
    call check("burn ledger closed by fault", thr%current_impulse, 0d0, 1d-12, nFail)

    end subroutine test_fault

    !--------------------------------------------------------------------------
    ! A span consumed while active+faulty must not be accrued after clear_fault.
    ! Reaching that state needs a direct component write, since set_fault()
    ! always switches the engine off; the write is deliberate here.
    !--------------------------------------------------------------------------
    subroutine test_faulted_span_not_deferred(thrVal, nFail)
    double precision, intent(in) :: thrVal
    integer, intent(inout) :: nFail
    type(type_thruster) :: thr
    type(type_t) :: t0, t1, t2

    print *, ""
    print *, "[8] faulted span is consumed, not deferred"

    thr = new_thruster()
    call mkTime(t0, 0, 0)
    call mkTime(t1, 1, 0)
    call mkTime(t2, 2, 0)

    call thr%activate(t0)
    thr%is_faulty = .true.                 ! direct write: active but faulted
    call check("faulted hour accrues nothing", thr%get_total_impulse(t1), 0d0, 1d-12, nFail)
    thr%is_faulty = .false.                ! fault cleared
    ! Only [1 h, 2 h] may accrue; [0, 1 h] was consumed while faulted
    call check("no retroactive accrual", thr%get_total_impulse(t2), thrVal * 3600d0, 1d-9, nFail)

    end subroutine test_faulted_span_not_deferred

    !--------------------------------------------------------------------------
    ! Passing t to set_current_thrust settles the old rate before the change
    !--------------------------------------------------------------------------
    subroutine test_rate_change_mid_burn(thrVal, nFail)
    double precision, intent(in) :: thrVal
    integer, intent(inout) :: nFail
    type(type_thruster) :: thr
    type(type_t) :: t0, t1, t2
    double precision :: fastVal

    print *, ""
    print *, "[9] rate change mid-burn settles the old rate first"

    thr = new_thruster()
    call mkTime(t0, 0, 0)
    call mkTime(t1, 1, 0)
    call mkTime(t2, 2, 0)

    fastVal = 0.1d0 * 0.95d0                ! raised to max_thrust at 1 h

    call thr%activate(t0)
    call thr%set_current_thrust(0.1d0, t1)  ! flush at 0.0475 N, then raise to 0.095 N
    ! 1 h at the old rate + 1 h at the new rate
    call check("mixed-rate total", thr%get_total_impulse(t2), &
               (thrVal + fastVal) * 3600d0, 1d-9, nFail)

    end subroutine test_rate_change_mid_burn

    !--------------------------------------------------------------------------
    ! Deactivating a thruster that was never running keeps the thrust level
    !--------------------------------------------------------------------------
    subroutine test_deactivate_without_activation(thrVal, nFail)
    double precision, intent(in) :: thrVal
    integer, intent(inout) :: nFail
    type(type_thruster) :: thr
    type(type_t) :: t0, t1

    print *, ""
    print *, "[10] deactivate without activation keeps the commanded level"

    thr = new_thruster()
    call mkTime(t0, 0, 0)
    call mkTime(t1, 1, 0)

    call thr%deactivate(t0)                 ! never activated
    call check("commanded level preserved", thr%current_thrust, thrVal, 1d-12, nFail)

    call thr%activate(t0)
    call check("burn works after that", thr%get_total_impulse(t1), thrVal * 3600d0, 1d-9, nFail)

    end subroutine test_deactivate_without_activation

    !--------------------------------------------------------------------------
    ! toggle switches on and off, settling on the way down
    !--------------------------------------------------------------------------
    subroutine test_toggle(thrVal, mdot, nFail)
    double precision, intent(in) :: thrVal, mdot
    integer, intent(inout) :: nFail
    type(type_thruster) :: thr
    type(type_t) :: t0, t1, t3, t4

    print *, ""
    print *, "[11] toggle on / off / on"

    thr = new_thruster()
    call mkTime(t0, 0, 0)
    call mkTime(t1, 1, 0)
    call mkTime(t3, 3, 0)
    call mkTime(t4, 4, 0)

    call thr%toggle(t0)                     ! on
    call check("toggled on", merge(1d0, 0d0, thr%is_active), 1d0, 1d-12, nFail)
    call thr%toggle(t1)                     ! off, settling the burn
    call check("toggled off", merge(1d0, 0d0, thr%is_active), 0d0, 1d-12, nFail)
    call check("burn of 1 h", thr%get_total_impulse(t3), thrVal * 3600d0, 1d-9, nFail)
    call check("no thrust while off", norm2(thr%get_thrust()), 0d0, 1d-12, nFail)

    call thr%toggle(t3)                     ! on again, fresh burn
    call check("fresh burn after toggle on", thr%current_impulse, 0d0, 1d-12, nFail)
    call check("second burn", thr%get_total_impulse(t4), thrVal * 3600d0 * 2d0, 1d-9, nFail)
    call check("second burn fuel", thr%get_total_fuel_used(t4), mdot * 7200d0, 1d-12, nFail)

    end subroutine test_toggle

    !--------------------------------------------------------------------------
    ! With a fuel budget, each burn draws the tank down
    !--------------------------------------------------------------------------
    subroutine test_fuel_budget(mdot, nFail)
    double precision, intent(in) :: mdot
    integer, intent(inout) :: nFail
    type(type_thruster) :: thr
    type(type_t) :: t0, t1, t2, t3
    double precision :: fuel0

    print *, ""
    print *, "[12] fuel budget is drawn down by each burn"

    fuel0 = 5d0
    thr = new_thruster(fuel0)
    call mkTime(t0, 0, 0)
    call mkTime(t1, 1, 0)
    call mkTime(t2, 2, 0)
    call mkTime(t3, 3, 0)

    call check("tank starts full", thr%current_fuel_mass, fuel0, 1d-12, nFail)

    call thr%activate(t0)
    call thr%update_state(t1)
    call check("after 1 h of burn", thr%current_fuel_mass, fuel0 - mdot * 3600d0, 1d-12, nFail)
    call thr%deactivate(t1)
    call check("tank unchanged while off", thr%current_fuel_mass, fuel0 - mdot * 3600d0, 1d-12, nFail)

    call thr%activate(t2)                   ! 1 h more, after an idle hour
    call thr%update_state(t3)
    call check("after 2 h of burn", thr%current_fuel_mass, fuel0 - mdot * 7200d0, 1d-12, nFail)
    call check("used equals the drawdown", thr%total_fuel_used, &
               fuel0 - thr%current_fuel_mass, 1d-12, nFail)
    call check("fuel flow rate", thr%get_fuel_flow_rate(), mdot, 1d-15, nFail)

    end subroutine test_fuel_budget

    !--------------------------------------------------------------------------
    ! A burn longer than the fuel allows is truncated at the flame-out epoch,
    ! and the engine delivers nothing afterwards
    !--------------------------------------------------------------------------
    subroutine test_flame_out(thrVal, mdot, nFail)
    double precision, intent(in) :: thrVal, mdot
    integer, intent(inout) :: nFail
    type(type_thruster) :: thr
    type(type_t) :: t0, t1, t2, t3
    double precision :: fuel0, burnMax

    print *, ""
    print *, "[13] flame-out truncates the burn (expect one warning)"

    burnMax = 5400d0                        ! 1.5 h of propellant ...
    fuel0 = mdot * burnMax
    thr = new_thruster(fuel0)
    call mkTime(t0, 0, 0)
    call mkTime(t1, 1, 0)
    call mkTime(t2, 2, 0)                   ! ... asked to burn for 2 h
    call mkTime(t3, 3, 0)

    call thr%activate(t0)
    call check("impulse stops at flame-out", thr%get_total_impulse(t2), &
               thrVal * burnMax, 1d-9, nFail)
    call check("fuel stops at flame-out", thr%get_total_fuel_used(t2), fuel0, 1d-12, nFail)
    call check("operating time stops too", thr%total_operation_time, burnMax, 1d-9, nFail)
    call check("tank is empty", thr%current_fuel_mass, 0d0, 1d-12, nFail)
    call check("no thrust from an empty tank", norm2(thr%get_thrust()), 0d0, 1d-12, nFail)

    ! A dry tank accrues nothing, at any later epoch
    call check("nothing accrues after flame-out", thr%get_total_impulse(t3), &
               thrVal * burnMax, 1d-12, nFail)
    call check("still no thrust", norm2(thr%get_thrust()), 0d0, 1d-12, nFail)

    end subroutine test_flame_out

    !--------------------------------------------------------------------------
    ! initial_fuel_mass = 0 means "no fuel budget modelled": nothing flames out
    !--------------------------------------------------------------------------
    subroutine test_unlimited_fuel(thrVal, nFail)
    double precision, intent(in) :: thrVal
    integer, intent(inout) :: nFail
    type(type_thruster) :: thr
    type(type_t) :: t0, t9

    print *, ""
    print *, "[14] a zero initial fuel mass means no budget"

    thr = new_thruster()                    ! default: initial_fuel_mass = 0
    call mkTime(t0, 0, 0)
    call mkTime(t9, 9, 0)

    call thr%activate(t0)
    call check("9 h of burn, no flame-out", thr%get_total_impulse(t9), thrVal * 32400d0, 1d-9, nFail)
    call check("thrust still available", norm2(thr%get_thrust()), thrVal, 1d-12, nFail)

    end subroutine test_unlimited_fuel

end program test_thruster
