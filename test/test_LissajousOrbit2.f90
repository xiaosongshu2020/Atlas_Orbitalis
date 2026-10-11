!===============================================================================
! test_LissajousOrbit2.f90
!
! Purpose: Generate a Lissajous orbit around the Earth-Moon L2 point and
!          save its ephemeris (position/velocity) to a text file.
!
! Output file: data/output/Lissajous_ephemeris.txt
!   Columns : MJD  x(km)  y(km)  z(km)  vx(km/s)  vy(km/s)  vz(km/s)
!   Frame   : Earth-Moon synodic (rotating) frame
!     origin : Earth-Moon barycenter
!     x-axis : along the Earth -> Moon line
!     z-axis : normal to the Earth-Moon orbital plane (angular momentum)
!     y-axis : completes the right-handed system
!   Units   : km for position, km/s for velocity
!
! Note: type_LissajousOrbit%getPV returns the linearized Lissajous solution
!       in CRTBP normalized units (length = Earth-Moon distance, time = TU),
!       relative to the Lagrange point. This program converts the state to
!       the barycentric synodic frame and scales to km and km/s.
!===============================================================================

program test_LissajousOrbit2
   use module_CRTBP
   use module_LissajousOrbit
   use module_const
   use module_time
   implicit none

   ! CRTBP system and Lissajous orbit
   type(type_CRTBP) :: earth_moon
   type(type_LissajousOrbit) :: em_lissajous

   ! Time variables
   type(type_t) :: t0, t
   real*8 :: T_xy_days, span_days, step_days

   ! Lagrange point parameters
   real*8 :: x0, C0, d2, d3, k
   real*8 :: omega_xy, omega_z, k_coeff
   logical :: success

   ! Lissajous state (normalized, relative to Lagrange point)
   real*8 :: xi, eta, zeta, xi_dot, eta_dot, zeta_dot
   ! Synodic-frame state in km and km/s (barycentric)
   real*8 :: pos_syn(3), vel_syn(3)

   integer :: i, nSteps, file_unit
   integer :: io_stat
   character(len=100) :: output_file

   !===========================================================================
   ! 1. BUILD EARTH-MOON CRTBP SYSTEM
   !===========================================================================

   call earth_moon%init('Earth-Moon', EARTH_MASS, MOON_MASS, &
                        MOON_EARTH_LENGTH / 1000.0d0)  ! m -> km

   write(*,*) '=============================================='
   write(*,*) 'Lissajous Orbit Ephemeris Generation'
   write(*,*) '=============================================='
   write(*,*) 'CRTBP system : ', trim(earth_moon%name)
   write(*,*) 'Mass ratio mu: ', earth_moon%mu
   write(*,*) 'Length unit  : ', earth_moon%a12, ' km'
   write(*,*) 'Time unit TU : ', earth_moon%TU, ' s (', &
              earth_moon%TU / 86400.0d0, ' days)'
   write(*,*) ''

   !===========================================================================
   ! 2. INITIALIZE LISSAJOUS ORBIT AROUND L2
   !===========================================================================

   call t0%init(58000.0d0)  ! MJD 58000.0

   call em_lissajous%init('EM-L2-Lissajous', LAGRANGE_L2, &
                          0.03d0, 0.035d0, &  ! A_xy = 0.05, A_z = 0.02
                          0.0d0, 0.0d0, &    ! phi_xy = 0, phi_z = 0
                          t0, earth_moon)

   call earth_moon%getLagrangePointParams(LAGRANGE_L2, x0, C0, d2, d3, k, success)
   call em_lissajous%getOrbitParams(omega_xy, omega_z, k_coeff, success)

   if (.not. success) then
       write(*,*) 'Error: Failed to get Lagrange point parameters'
       stop
   endif

   ! In-plane period in days: T_xy = 2*pi / omega_xy (in TU)
   T_xy_days = 2.0d0*PI / omega_xy * earth_moon%TU / 86400.0d0

   write(*,*) 'Lissajous Orbit Parameters:'
   write(*,*) '  Lagrange point: L', em_lissajous%LID
   write(*,*) '  In-plane amplitude A_xy: ', em_lissajous%A_xy, &
              ' (', em_lissajous%A_xy * earth_moon%a12, ' km)'
   write(*,*) '  Out-of-plane amplitude A_z: ', em_lissajous%A_z, &
              ' (', em_lissajous%A_z * earth_moon%a12, ' km)'
   write(*,*) '  In-plane period T_xy: ', T_xy_days, ' days'
   write(*,*) '  Out-of-plane period T_z: ', &
              2.0d0*PI / omega_z * earth_moon%TU / 86400.0d0, ' days'
   write(*,*) ''

   !===========================================================================
   ! 3. GENERATE EPHEMERIS AND SAVE TO FILE
   !===========================================================================

   ! Time span: 
   span_days = 50.0d0 * T_xy_days
   step_days = 0.5d0
   nSteps = int(span_days / step_days)

   output_file = 'data/output/Lissajous_ephemeris.txt'

   open(newunit=file_unit, file=trim(output_file), status='replace', &
        action='write', iostat=io_stat)
   if (io_stat /= 0) then
       write(*,*) 'Error: Cannot open output file ', trim(output_file)
       stop
   endif

   ! Header (all lines start with '#')
   write(file_unit, '(A)') '# ============================================================================'
   write(file_unit, '(A)') '# LISSAJOUS ORBIT EPHEMERIS'
   write(file_unit, '(A)') '# ============================================================================'
   write(file_unit, '(A)') '# Orbit name: EM-L2-Lissajous (Earth-Moon L2 point)'
   write(file_unit, '(A, F12.6)') '# Mass ratio mu      : ', earth_moon%mu
   write(file_unit, '(A, F12.3)') '# Length unit (km)   : ', earth_moon%a12
   write(file_unit, '(A, F12.3)') '# Time unit TU (s)   : ', earth_moon%TU
   write(file_unit, '(A, F12.6)') '# A_xy (normalized)  : ', em_lissajous%A_xy
   write(file_unit, '(A, F12.6)') '# A_z  (normalized)  : ', em_lissajous%A_z
   write(file_unit, '(A, F12.4)') '# In-plane period (days): ', T_xy_days
   write(file_unit, '(A)') '# ============================================================================'
   write(file_unit, '(A)') '# Coordinate frame: Earth-Moon synodic (rotating) frame'
   write(file_unit, '(A)') '#   origin: Earth-Moon barycenter'
   write(file_unit, '(A)') '#   x-axis: along the Earth -> Moon line'
   write(file_unit, '(A)') '#   z-axis: normal to the Earth-Moon orbital plane'
   write(file_unit, '(A)') '#   y-axis: completes the right-handed system'
   write(file_unit, '(A)') '# Velocity is the time derivative of the rotating-frame position'
   write(file_unit, '(A)') '# ============================================================================'
   write(file_unit, '(A)') '# Units: km, km/s'
   write(file_unit, '(A)') '# Columns: MJD  x(km)  y(km)  z(km)  vx(km/s)  vy(km/s)  vz(km/s)'
   write(file_unit, '(A)') '# ============================================================================'

   write(*,*) 'Generating ephemeris...'
   write(*,*) '  Time span: ', span_days, ' days'
   write(*,*) '  Time step: ', step_days, ' day'
   write(*,*) '  Number of epochs: ', nSteps + 1
   write(*,*) ''

   do i = 0, nSteps
       call t%init(t0%val + i * step_days)

       ! Lissajous state relative to L2 (normalized units)
       call em_lissajous%getPV(t, xi, eta, zeta, xi_dot, eta_dot, zeta_dot)

       ! Convert to barycentric synodic frame, km and km/s
       !   position: (x0 + xi, eta, zeta) * a12            [km]
       !   velocity: (xi_dot, eta_dot, zeta_dot) * a12/TU  [km/s]
       pos_syn(1) = (x0 + xi)   * earth_moon%a12
       pos_syn(2) = eta         * earth_moon%a12
       pos_syn(3) = zeta        * earth_moon%a12
       vel_syn(1) = xi_dot      * earth_moon%a12 / earth_moon%TU
       vel_syn(2) = eta_dot     * earth_moon%a12 / earth_moon%TU
       vel_syn(3) = zeta_dot    * earth_moon%a12 / earth_moon%TU

       write(file_unit, '(F16.8, 6F16.6)') t%val, pos_syn, vel_syn

       ! Print first few and last epochs for inspection
       if (i <= 2 .or. i == nSteps) then
           write(*, '(F16.8, 6F16.6)') t%val, pos_syn, vel_syn
       endif
   end do

   close(file_unit)

   write(*,*) ''
   write(*,*) 'Ephemeris saved to: ', trim(output_file)
   write(*,*) '  Total epochs: ', nSteps + 1
   write(*,*) ''
   write(*,*) '=============================================='
   write(*,*) 'Test completed successfully!'
   write(*,*) '=============================================='

end program test_LissajousOrbit2
