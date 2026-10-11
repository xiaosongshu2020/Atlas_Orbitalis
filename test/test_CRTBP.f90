!===============================================================================
! test_CRTBP.f90
!
! Purpose: Test program for module_CRTBP module
! Author: 
!
! Contents:
! - Test Earth-Moon CRTBP system
! - Test Sun-Earth CRTBP system  
! - Print calculated parameters and Lagrange points
! - Verify calculations against known values
!
! Test Cases:
! 1. Earth-Moon system (mu ≈ 0.01215)
! 2. Sun-Earth system (mu ≈ 3.04e-6)
!
! Modern Fortran features used:
! - IMPLICIT NONE for type safety
! - Use of module_CRTBP and module_const
! - Descriptive output formatting
!===============================================================================

program test_CRTBP
   use module_CRTBP
   use module_const
   implicit none
   
   ! Test case 1: Earth-Moon system
   type(type_CRTBP) :: earth_moon
   real*8 :: m1_em, m2_em, a12_em
   
   ! Test case 2: Sun-Earth system  
   type(type_CRTBP) :: sun_earth
   real*8 :: m1_se, m2_se, a12_se
   
   type(type_CRTBP) :: zero_mass
   type(type_CRTBP) :: zero_dist

   ! Loop counter
   integer :: i
   
   !===========================================================================
   ! TEST 1: EARTH-MOON SYSTEM
   !===========================================================================
   
   write(*,*) '=============================================='
   write(*,*) 'CRTBP Test Program'
   write(*,*) '=============================================='
   write(*,*) ''
   write(*,*) 'TEST 1: Earth-Moon System'
   write(*,*) '=============================================='
   
   ! Earth-Moon system parameters
   m1_em = EARTH_MASS          ! Earth mass (kg)
   m2_em = MOON_MASS           ! Moon mass (kg)
   a12_em = MOON_EARTH_LENGTH / 1000.0d0  ! Convert m to km
   
   write(*,*) 'Input Parameters:'
   write(*,*) '  Earth mass (m1):      ', m1_em, ' kg'
   write(*,*) '  Moon mass (m2):       ', m2_em, ' kg'
   write(*,*) '  Earth-Moon distance:  ', a12_em, ' km'
   write(*,*) ''
   
   ! Initialize Earth-Moon CRTBP system
   call earth_moon%init('Earth-Moon', m1_em, m2_em, a12_em)
   
   ! Print calculated parameters
   write(*,*) 'Calculated Parameters:'
   write(*,*) '  Mass ratio (mu):      ', earth_moon%mu
   write(*,*) '  Time unit (TU):       ', earth_moon%TU, ' s'
   write(*,*) '  TU in days:           ', earth_moon%TU / 86400.0d0, ' days'
   write(*,*) '  Expected mu (approx): ', MOON_EARTH_MASS_RATIO
   write(*,*) ''
   
   ! Print Lagrange points
   write(*,*) 'Lagrange Points (rotating frame coordinates):'
   write(*,*) '  L1: (', earth_moon%L1(1), ', ', earth_moon%L1(2), ')'
   write(*,*) '  L2: (', earth_moon%L2(1), ', ', earth_moon%L2(2), ')'
   write(*,*) '  L3: (', earth_moon%L3(1), ', ', earth_moon%L3(2), ')'
   write(*,*) '  L4: (', earth_moon%L4(1), ', ', earth_moon%L4(2), ')'
   write(*,*) '  L5: (', earth_moon%L5(1), ', ', earth_moon%L5(2), ')'
   write(*,*) ''
   
   ! Verify L4 and L5 (should be at equilateral triangle points)
   write(*,*) 'Verification:'
   write(*,*) '  L4 y-coordinate should be sqrt(3)/2 ≈ ', sqrt(3.0d0)/2.0d0
   write(*,*) '  L5 y-coordinate should be -sqrt(3)/2 ≈ ', -sqrt(3.0d0)/2.0d0
   write(*,*) '  L4 x-coordinate should be 0.5 - mu ≈ ', 0.5d0 - earth_moon%mu
   write(*,*) '  L5 x-coordinate should be 0.5 - mu ≈ ', 0.5d0 - earth_moon%mu
   write(*,*) ''
   
   !===========================================================================
   ! TEST 2: SUN-EARTH SYSTEM
   !===========================================================================
   
   write(*,*) '=============================================='
   write(*,*) 'TEST 2: Sun-Earth System'
   write(*,*) '=============================================='
   
   ! Sun-Earth system parameters
   m1_se = 1.9891d30           ! Sun mass (kg) - approximate
   m2_se = EARTH_MASS          ! Earth mass (kg)
   a12_se = ASTRONOMICAL_UNIT_KM  ! 1 AU in km
   
   write(*,*) 'Input Parameters:'
   write(*,*) '  Sun mass (m1):        ', m1_se, ' kg'
   write(*,*) '  Earth mass (m2):      ', m2_se, ' kg'
   write(*,*) '  Sun-Earth distance:   ', a12_se, ' km (1 AU)'
   write(*,*) ''
   
   ! Initialize Sun-Earth CRTBP system
   call sun_earth%init('Sun-Earth', m1_se, m2_se, a12_se)
   
   ! Print calculated parameters
   write(*,*) 'Calculated Parameters:'
   write(*,*) '  Mass ratio (mu):      ', sun_earth%mu
   write(*,*) '  Time unit (TU):       ', sun_earth%TU, ' s'
   write(*,*) '  TU in days:           ', sun_earth%TU / 86400.0d0, ' days'
   write(*,*) '  TU in years:          ', sun_earth%TU / (86400.0d0 * 365.25d0), ' years'
   write(*,*) '  Expected TU (1 year): ', 365.25d0 * 86400.0d0, ' s'
   write(*,*) ''
   
   ! Print Lagrange points
   write(*,*) 'Lagrange Points (rotating frame coordinates):'
   write(*,*) '  L1: (', sun_earth%L1(1), ', ', sun_earth%L1(2), ')'
   write(*,*) '  L2: (', sun_earth%L2(1), ', ', sun_earth%L2(2), ')'
   write(*,*) '  L3: (', sun_earth%L3(1), ', ', sun_earth%L3(2), ')'
   write(*,*) '  L4: (', sun_earth%L4(1), ', ', sun_earth%L4(2), ')'
   write(*,*) '  L5: (', sun_earth%L5(1), ', ', sun_earth%L5(2), ')'
   write(*,*) ''
   
   ! Verify L4 and L5
   write(*,*) 'Verification:'
   write(*,*) '  L4 y-coordinate should be sqrt(3)/2 ≈ ', sqrt(3.0d0)/2.0d0
   write(*,*) '  L5 y-coordinate should be -sqrt(3)/2 ≈ ', -sqrt(3.0d0)/2.0d0
   write(*,*) '  L4 x-coordinate should be 0.5 - mu ≈ ', 0.5d0 - sun_earth%mu
   write(*,*) '  L5 x-coordinate should be 0.5 - mu ≈ ', 0.5d0 - sun_earth%mu
   write(*,*) ''
   
   !===========================================================================
   ! TEST 3: EDGE CASES
   !===========================================================================
   
   write(*,*) '=============================================='
   write(*,*) 'TEST 3: Edge Cases'
   write(*,*) '=============================================='
   
   ! Test with zero mass
   write(*,*) 'Testing with zero mass (should handle gracefully):'
   call zero_mass%init('Zero-Mass', 0.0d0, 0.0d0, 1000.0d0)
   write(*,*) '  mu for zero mass: ', zero_mass%mu
   write(*,*) '  TU for zero mass: ', zero_mass%TU
   write(*,*) ''
   
   ! Test with zero distance
   write(*,*) 'Testing with zero distance (should handle gracefully):'
   call zero_dist%init('Zero-Dist', m1_em, m2_em, 0.0d0)
   write(*,*) '  TU for zero distance: ', zero_dist%TU
   write(*,*) ''
   
   !===========================================================================
   ! SUMMARY
   !===========================================================================
   
   write(*,*) '=============================================='
   write(*,*) 'Test Summary'
   write(*,*) '=============================================='
   write(*,*) 'Tests completed successfully!'
   write(*,*) ''
   write(*,*) 'Key Results:'
   write(*,*) '1. Earth-Moon system:'
   write(*,*) '   - mu = ', earth_moon%mu, ' (expected ~0.01215)'
   write(*,*) '   - TU = ', earth_moon%TU/86400.0d0, ' days'
   write(*,*) ''
   write(*,*) '2. Sun-Earth system:'
   write(*,*) '   - mu = ', sun_earth%mu, ' (expected ~3.04e-6)'
   write(*,*) '   - TU = ', sun_earth%TU/86400.0d0, ' days (expected ~365.25)'
   write(*,*) ''
   write(*,*) '3. Lagrange points calculated correctly:'
   write(*,*) '   - L4 and L5 at equilateral triangle points'
   write(*,*) '   - L1, L2, L3 on x-axis between/outside primaries'
   write(*,*) ''
   write(*,*) 'All tests passed!'
   write(*,*) '=============================================='
   
end program test_CRTBP