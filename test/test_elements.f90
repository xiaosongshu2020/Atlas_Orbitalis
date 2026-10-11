!===============================================================================
! test_elements.f90
!
! Purpose: Test program for type_orbit — OOP encapsulation of module_elements
! Author: Shawn Hu
! Date: 2026-07-11
!
! Test Cases:
!   1. init()          — initialize from SI units for Keplerian
!   2. toCartesian()   — convert Keplerian → Cartesian
!   3. toKeplerian()   — round-trip: Keplerian → Cartesian → Keplerian
!   4. toNSG1/toNSG2() — convert between non-singular element types
!   5. get_rEf()       — solve Kepler equation
!   6. get_rut()       — solve generalized Kepler equation
!   7. propagateTwoBody() — two-body (Keplerian) propagation
!   8. toSI/toNormalized() — unit conversion round-trip
!   9. show()          — display orbit info
!  10. chained usage   — full round-trip KPL → CTS → NSG2 → NSG1 → KPL
!  11. init() Cartesian — initialize from Cartesian SI input
!===============================================================================

program test_elements
    use module_const
    use module_elements
    implicit none

    ! -- orbit objects --
    type(type_orbit) :: orb, orb2

    ! -- reference Keplerian elements (SI units) --
    double precision :: kpl_si(6), kpl_norm(6), cts_si(6), cts_norm(6)
    double precision :: a_ref, e_ref

    ! -- Kepler equation outputs --
    double precision :: r_val, E_val, f_val, u_val, ut_val

    ! -- saved intermediate states for round-trip comparison --
    double precision :: nsg1_saved(6), norm_saved(6)

    ! -- two-body propagation --
    double precision :: n_mean, dt_prop, M0, M1_expected

    ! -- error tracking --
    double precision :: err_val, max_err, tol
    integer :: nPass, nFail, nTotal

    !===========================================================================
    ! Reference orbit: near-circular LEO
    !   Keplerian (SI): a=7000 km, e=0.01, i=45°, Omega=30°, omega=60°, M=90°
    !===========================================================================
    nPass = 0
    nFail = 0
    nTotal = 0
    tol   = 1.0d-9

    a_ref = 7000.0d0
    e_ref = 0.01d0

    kpl_si(1) = a_ref
    kpl_si(2) = e_ref
    kpl_si(3) = 45.0d0 * deg
    kpl_si(4) = 30.0d0 * deg
    kpl_si(5) = 60.0d0 * deg
    kpl_si(6) = 90.0d0 * deg

    ! Pre-compute normalized Keplerian for comparisons
    call kpl_SI2normal(kpl_si, kpl_norm)

    write(*, '(a)') ""
    write(*, '(a)') "============================================================"
    write(*, '(a)') "  test_elements — type_orbit OOP verification"
    write(*, '(a)') "============================================================"
    write(*, '(a)') ""
    write(*, '(a,f10.3)') "  Reference: a = ", a_ref
    write(*, '(a,f10.6)') "             e = ", e_ref
    write(*, '(a,f10.4)') "             i = ", kpl_si(3) / deg
    write(*, '(a,f10.4)') "         Omega = ", kpl_si(4) / deg
    write(*, '(a,f10.4)') "         omega = ", kpl_si(5) / deg
    write(*, '(a,f10.4)') "             M = ", kpl_si(6) / deg

    !===========================================================================
    ! TEST 1: init() — Keplerian from SI
    !===========================================================================
    nTotal = nTotal + 1
    write(*, '(a)') ""
    write(*, '(a)') "--- TEST 1: init() with Keplerian SI input ---"

    call orb%init(kpl_si, STATE_TYPE_KEPLERIAN, isNormalized=.false.)

    write(*, '(a)') "  After init (internal normalized units):"
    write(*, '(a,6es16.8)') "    elements = ", orb%elements
    write(*, '(a,i0)')     "    eleType  = ", orb%eleType

    err_val = abs(orb%elements(1) - a_ref / length_unit)
    if (err_val < tol) then
        write(*, '(a)') "  PASS: a converted correctly to normalized units."
        nPass = nPass + 1
    else
        write(*, '(a,es16.8,a,es16.8)') "  FAIL: a mismatch: ", orb%elements(1), &
                                          " vs ", a_ref / length_unit
        nFail = nFail + 1
    end if

    !===========================================================================
    ! TEST 2: toCartesian() — convert to Cartesian
    !===========================================================================
    nTotal = nTotal + 1
    write(*, '(a)') ""
    write(*, '(a)') "--- TEST 2: toCartesian() ---"

    call orb%toCartesian()

    write(*, '(a)') "  After toCartesian():"
    write(*, '(a,3es16.8)') "    r = ", orb%elements(1:3)
    write(*, '(a,3es16.8)') "    v = ", orb%elements(4:6)
    write(*, '(a,i0)')     "    eleType = ", orb%eleType

    if (orb%eleType == STATE_TYPE_CARTESIAN) then
        write(*, '(a)') "  PASS: eleType correctly set to CARTESIAN."
        nPass = nPass + 1
    else
        write(*, '(a)') "  FAIL: wrong eleType after toCartesian()."
        nFail = nFail + 1
    end if

    !===========================================================================
    ! TEST 3: toKeplerian() — round-trip back to Keplerian
    !===========================================================================
    nTotal = nTotal + 1
    write(*, '(a)') ""
    write(*, '(a)') "--- TEST 3: toKeplerian() round-trip ---"

    call orb%toKeplerian()

    max_err = maxval(abs(orb%elements - kpl_norm))
    write(*, '(a,es12.4)')   "  Max round-trip error (kpl→cts→kpl): ", max_err
    write(*, '(a,6es16.8)') "    recovered = ", orb%elements
    write(*, '(a,6es16.8)') "    reference = ", kpl_norm

    if (max_err < tol) then
        write(*, '(a)') "  PASS: round-trip KPL → CTS → KPL preserved."
        nPass = nPass + 1
    else
        write(*, '(a)') "  FAIL: round-trip error too large."
        nFail = nFail + 1
    end if

    !===========================================================================
    ! TEST 4: toNSG1() and toNSG2()
    !===========================================================================
    nTotal = nTotal + 1
    write(*, '(a)') ""
    write(*, '(a)') "--- TEST 4: toNSG1() and toNSG2() ---"

    call orb%init(kpl_si, STATE_TYPE_KEPLERIAN, isNormalized=.false.)

    ! Convert to NSG1
    call orb%toNSG1()
    if (orb%eleType == STATE_TYPE_NONSINGULAR_1) then
        write(*, '(a)') "  toNSG1() — eleType correct."
        write(*, '(a,6es16.8)') "    nsg1 = ", orb%elements
    else
        write(*, '(a)') "  FAIL: wrong eleType after toNSG1()."
        nFail = nFail + 1
    end if

    ! Save NSG1 for later comparison
    nsg1_saved = orb%elements

    ! Convert to NSG2
    call orb%toNSG2()
    if (orb%eleType == STATE_TYPE_NONSINGULAR_2) then
        write(*, '(a)') "  toNSG2() — eleType correct."
        write(*, '(a,6es16.8)') "    nsg2 = ", orb%elements
    else
        write(*, '(a)') "  FAIL: wrong eleType after toNSG2()."
        nFail = nFail + 1
    end if

    ! Round-trip: NSG2 → NSG1
    call orb%toNSG1()
    max_err = maxval(abs(orb%elements - nsg1_saved))
    write(*, '(a,es12.4)') "  nsg1→nsg2→nsg1 max error = ", max_err

    ! Round-trip all the way: NSG1 → Keplerian → compare
    call orb%toKeplerian()
    max_err = maxval(abs(orb%elements - kpl_norm))
    write(*, '(a,es12.4)') "  kpl→nsg1→nsg2→nsg1→kpl max error = ", max_err

    if (max_err < tol) then
        write(*, '(a)') "  PASS: full NSG round-trip preserved."
        nPass = nPass + 1
    else
        write(*, '(a)') "  FAIL: NSG round-trip error too large."
        nFail = nFail + 1
    end if

    !===========================================================================
    ! TEST 5: get_rEf() — solve Kepler equation
    !===========================================================================
    nTotal = nTotal + 1
    write(*, '(a)') ""
    write(*, '(a)') "--- TEST 5: get_rEf() — Kepler equation solver ---"

    call orb%init(kpl_si, STATE_TYPE_KEPLERIAN, isNormalized=.false.)
    call orb%get_rEf(r_val, E_val, f_val)

    write(*, '(a,f12.6)') "  r  (norm. units) = ", r_val
    write(*, '(a,f12.6)') "  E  (rad)         = ", E_val
    write(*, '(a,f12.6)') "  f  (rad)         = ", f_val

    ! Verify: r should equal a*(1 - e*cos(E))
    err_val = abs(r_val - orb%elements(1)*(1.0d0 - e_ref*cos(E_val)))
    write(*, '(a,es12.4)') "  r - a*(1-e*cosE)  = ", err_val

    if (err_val < tol) then
        write(*, '(a)') "  PASS: Kepler equation solved correctly."
        nPass = nPass + 1
    else
        write(*, '(a)') "  FAIL: Kepler equation error too large."
        nFail = nFail + 1
    end if

    ! Test get_rEf auto-converts non-Keplerian input
    nTotal = nTotal + 1
    call orb%toCartesian()
    call orb%get_rEf(r_val, E_val, f_val)   ! should auto-convert
    if (orb%eleType == STATE_TYPE_KEPLERIAN) then
        write(*, '(a)') "  PASS: get_rEf auto-converted Cartesian to Keplerian."
        nPass = nPass + 1
    else
        write(*, '(a)') "  FAIL: get_rEf did not auto-convert."
        nFail = nFail + 1
    end if

    !===========================================================================
    ! TEST 6: get_rut() — solve generalized Kepler equation
    !===========================================================================
    nTotal = nTotal + 1
    write(*, '(a)') ""
    write(*, '(a)') "--- TEST 6: get_rut() — generalized Kepler equation ---"

    call orb%init(kpl_si, STATE_TYPE_KEPLERIAN, isNormalized=.false.)
    call orb%toNSG1()
    call orb%get_rut(r_val, u_val, ut_val)

    write(*, '(a,f12.6)') "  r  (norm. units) = ", r_val
    write(*, '(a,f12.6)') "  u  (rad)         = ", u_val
    write(*, '(a,f12.6)') "  ut (rad)         = ", ut_val

    ! Verify: r = a*(1 - xi*cos(ut) - eta*sin(ut))
    err_val = abs(r_val - orb%elements(1) &
                  * (1.0d0 - orb%elements(4)*cos(ut_val) &
                            - orb%elements(5)*sin(ut_val)))
    write(*, '(a,es12.4)') "  r - a*(1-xi*cos_ut-eta*sin_ut) = ", err_val

    if (err_val < tol) then
        write(*, '(a)') "  PASS: generalized Kepler equation (NSG1) solved."
        nPass = nPass + 1
    else
        write(*, '(a)') "  FAIL: generalized Kepler equation error."
        nFail = nFail + 1
    end if

    ! Test get_rut on NSG2
    nTotal = nTotal + 1
    call orb%toNSG2()
    call orb%get_rut(r_val, u_val, ut_val)
    err_val = abs(r_val - orb%elements(1) &
                  * (1.0d0 - orb%elements(2)*cos(ut_val) &
                            - orb%elements(3)*sin(ut_val)))
    write(*, '(a,es12.4)') "  NSG2 — r error = ", err_val
    if (err_val < tol) then
        write(*, '(a)') "  PASS: get_rut works on NSG2."
        nPass = nPass + 1
    else
        write(*, '(a)') "  FAIL: get_rut on NSG2 error too large."
        nFail = nFail + 1
    end if

    !===========================================================================
    ! TEST 7: propagateTwoBody() — two-body propagation
    !===========================================================================
    nTotal = nTotal + 1
    write(*, '(a)') ""
    write(*, '(a)') "--- TEST 7: propagateTwoBody() ---"

    call orb%init(kpl_si, STATE_TYPE_KEPLERIAN, isNormalized=.false.)

    ! mean motion n = sqrt(mu/a^3), mu = 1 in normalized units
    n_mean  = orb%elements(1)**(-1.5d0)
    dt_prop = 100.0d0 / time_unit        ! 100 seconds in normalized time
    M0      = orb%elements(6)

    call orb%propagateTwoBody(dt_prop)

    M1_expected = dmod(M0 + n_mean*dt_prop, 3.141592653589793d0*2.0d0)
    err_val = abs(orb%elements(6) - M1_expected)
    write(*, '(a,es12.4)') "  M propagation error = ", err_val

    if (err_val < tol .and. orb%eleType == STATE_TYPE_KEPLERIAN) then
        write(*, '(a)') "  PASS: two-body propagation correct."
        nPass = nPass + 1
    else
        write(*, '(a)') "  FAIL: two-body propagation error."
        nFail = nFail + 1
    end if

    !===========================================================================
    ! TEST 8: toSI() / toNormalized() — unit conversion round-trip
    !===========================================================================
    nTotal = nTotal + 1
    write(*, '(a)') ""
    write(*, '(a)') "--- TEST 8: toSI() / toNormalized() round-trip ---"

    call orb%init(kpl_si, STATE_TYPE_KEPLERIAN, isNormalized=.false.)
    norm_saved = orb%elements

    ! Convert to SI and back
    call orb%toSI()
    write(*, '(a)') "  After toSI():"
    write(*, '(a,6es16.8)') "    elements = ", orb%elements
    write(*, '(a,i0)')     "    eleType  = ", orb%eleType

    max_err = maxval(abs(orb%elements - kpl_si))
    write(*, '(a,es12.4)') "  SI values vs original kpl_si, max err = ", max_err

    call orb%toNormalized()
    max_err = maxval(abs(orb%elements - norm_saved))
    write(*, '(a,es12.4)') "  Normalized round-trip error = ", max_err

    if (max_err < tol) then
        write(*, '(a)') "  PASS: unit conversion round-trip correct."
        nPass = nPass + 1
    else
        write(*, '(a)') "  FAIL: unit conversion error too large."
        nFail = nFail + 1
    end if

    ! Test toSI/toNormalized on Cartesian
    nTotal = nTotal + 1
    call orb%toCartesian()
    call orb%toSI()
    call orb%toNormalized()
    write(*, '(a)') "  toSI()/toNormalized() works on Cartesian."
    nPass = nPass + 1

    ! Test toSI/toNormalized on NSG1 (goes through Keplerian internally)
    nTotal = nTotal + 1
    call orb%toNSG1()
    nsg1_saved = orb%elements
    call orb%toSI()
    call orb%toNormalized()
    max_err = maxval(abs(orb%elements - nsg1_saved))
    write(*, '(a,es12.4)') "  NSG1 toSI/toNormalized round-trip error = ", max_err
    if (max_err < tol) then
        write(*, '(a)') "  PASS: NSG1 unit round-trip correct."
        nPass = nPass + 1
    else
        write(*, '(a)') "  FAIL: NSG1 unit round-trip error too large."
        nFail = nFail + 1
    end if

    !===========================================================================
    ! TEST 9: show() — display orbit
    !===========================================================================
    nTotal = nTotal + 1
    write(*, '(a)') ""
    write(*, '(a)') "--- TEST 9: show() ---"

    call orb%init(kpl_si, STATE_TYPE_KEPLERIAN, isNormalized=.false.)
    call orb%show()
    write(*, '(a)') "  PASS: show() ran without error."
    nPass = nPass + 1

    !===========================================================================
    ! TEST 10: Full chain — KPL → CTS → NSG2 → NSG1 → KPL
    !===========================================================================
    nTotal = nTotal + 1
    write(*, '(a)') ""
    write(*, '(a)') "--- TEST 10: Full chain KPL → CTS → NSG2 → NSG1 → KPL ---"

    call orb%init(kpl_si, STATE_TYPE_KEPLERIAN, isNormalized=.false.)

    call orb%toCartesian()
    call orb%toNSG2()
    call orb%toNSG1()
    call orb%toKeplerian()

    max_err = maxval(abs(orb%elements - kpl_norm))
    write(*, '(a,es12.4)') "  Full chain round-trip error = ", max_err

    if (max_err < tol) then
        write(*, '(a)') "  PASS: full chain preserved the orbit."
        nPass = nPass + 1
    else
        write(*, '(a)') "  FAIL: full chain error too large."
        nFail = nFail + 1
    end if

    !===========================================================================
    ! TEST 11: init() with Cartesian SI input
    !===========================================================================
    nTotal = nTotal + 1
    write(*, '(a)') ""
    write(*, '(a)') "--- TEST 11: init() with Cartesian SI input ---"

    ! Build Cartesian SI: kpl_norm → kpl2cts → Cartesian normalized → Cartesian SI
    call kpl2cts(kpl_norm, cts_norm)
    call cts_normal2SI(cts_norm, cts_si)

    call orb2%init(cts_si, STATE_TYPE_CARTESIAN, isNormalized=.false.)
    write(*, '(a)') "  After init with Cartesian SI:"
    write(*, '(a,3es16.8)') "    r (norm) = ", orb2%elements(1:3)
    write(*, '(a,3es16.8)') "    v (norm) = ", orb2%elements(4:6)

    ! Convert to Keplerian and compare with reference
    call orb2%toKeplerian()
    max_err = maxval(abs(orb2%elements - kpl_norm))
    write(*, '(a,es12.4)') "  Max error vs reference KPL = ", max_err

    if (max_err < tol) then
        write(*, '(a)') "  PASS: Cartesian SI init and conversion correct."
        nPass = nPass + 1
    else
        write(*, '(a)') "  FAIL: Cartesian init error."
        nFail = nFail + 1
    end if

    !===========================================================================
    ! SUMMARY
    !===========================================================================
    write(*, '(a)') ""
    write(*, '(a)') "============================================================"
    write(*, '(a)') "  TEST SUMMARY"
    write(*, '(a)') "============================================================"
    write(*, '(a,i0)') "  Passed: ", nPass
    write(*, '(a,i0)') "  Failed: ", nFail
    write(*, '(a,i0)') "  Total:  ", nTotal
    write(*, '(a)') ""

    if (nFail == 0) then
        write(*, '(a)') "  ALL TESTS PASSED."
    else
        write(*, '(a)') "  SOME TESTS FAILED!"
        stop 1
    end if

    write(*, '(a)') ""

end program test_elements
