! Test program for moon coordinate transformation
! Test LCF2LCI and LCI2LCF functions in module_space
! Author: Cline
! Date: 2025/10/22

program test_moon_coordinate_transformation
    use module_space
    implicit none
    
    ! Variable declarations
    real*8 :: t_utc, t_TT           ! TT time (JD)
    real*8 :: rf(3)           ! Position in moon-center-fixed
    real*8 :: ri(3)           ! Position in moon-center-inertial
    real*8 :: rf_back(3)      ! Position after back transformation
    real*8 :: error(3)        ! Error vector
    real*8 :: max_error       ! Maximum error
    
    ! Initialize test data
    t_utc = 60977d0
    call UTC2TT(t_utc, t_TT)
    t_TT = t_TT + 2400000.5d0
    
    ! Test position: point on moon surface (unit: km)
    rf = [-1232.102835d0, -1621.706463d0, 54.676440d0]  ! Point on moon equator
    
    write(*,*) "================================================"
    write(*,*) "Moon Coordinate Transformation Test Program"
    write(*,*) "================================================"
    write(*,*)
    write(*,*) "Test time (TT JD):", t_TT
    write(*,*) "Initial moon-center-fixed position (km):", rf
    write(*,*)
    
    ! Test LCF2LCI transformation
    write(*,*) "1. Test LCF2LCI transformation (Moon-fixed -> Moon-inertial)"
    call LCF2LCI(t_TT, rf, ri)
    write(*,*) "   Transformed moon-center-inertial position (km):", ri
    write(*,*)
    
    ! Test LCI2LCF transformation
    write(*,*) "2. Test LCI2LCF transformation (Moon-inertial -> Moon-fixed)"
    call LCI2LCF(t_TT, ri, rf_back)
    write(*,*) "   Back-transformed moon-center-fixed position (km):", rf_back
    write(*,*)
    
    ! Calculate error
    error = rf_back - rf
    max_error = maxval(abs(error))
    
    write(*,*) "3. Error analysis"
    write(*,*) "   Position error vector (km):", error
    write(*,*) "   Maximum absolute error (km):", max_error
    write(*,*)
    
    ! Check test result
    if (max_error < 1.0d-10) then
        write(*,*) "Test PASSED: Coordinate transformation accuracy meets requirement"
    else
        write(*,*) "Test WARNING: Transformation error detected"
    endif
    
    write(*,*)
    write(*,*) "================================================"
    write(*,*) "Test completed"
    write(*,*) "================================================"
    
end program test_moon_coordinate_transformation
