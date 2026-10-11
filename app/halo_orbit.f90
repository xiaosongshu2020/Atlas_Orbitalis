!===============================================================================
! Program: halo_orbit
!
! Purpose:
!   Computes a third-order analytical Halo orbit about the Earth-Moon L2 point
!   with module_halo and writes its ephemeris to
!       data/output/halo/halo_L2_ephemeris.txt
!
! Output file:
!   Columns : MJD  x(km)  y(km)  z(km)  vx(km/s)  vy(km/s)  vz(km/s)
!   Frame   : Earth-Moon synodic (rotating) frame
!     origin : Earth-Moon barycenter
!     x-axis : along the Earth -> Moon line
!     z-axis : normal to the Earth-Moon orbital plane (angular momentum)
!     y-axis : completes the right-handed system
!   Units   : km for position, km/s for velocity
!
! Note:
!   type_halo%getPV returns the analytical solution in CRTBP normalized units
!   (length = Earth-Moon distance a12, time = TU), relative to the libration
!   point and scaled by gamma (the distance from L2 to the Moon), i.e. the
!   length unit of the analytical solution is r1 = gamma*a12, as in
!   Richardson (1980). This program converts the state to the barycentric
!   synodic frame and scales it to km and km/s.
!
!   The synodic frame rotates with the Earth-Moon line, so the velocities are
!   the time derivatives of the rotating-frame position, not inertial ones.
!===============================================================================

program halo_orbit
 use module_const
 use module_time
 use module_CRTBP
 use module_halo
 implicit none

 ! Objects
 type(type_CRTBP) :: earth_moon_system
 type(type_halo)  :: orbit
 type(type_t)     :: t0, t

 ! Halo orbit parameters
 real(8)          :: beta, phi0
 integer          :: LID

 ! State relative to the libration point (normalized by gamma)
 real(8)          :: xi, eta, zeta, xi_dot, eta_dot, zeta_dot
 ! Synodic-frame state in km and km/s (barycentric)
 real(8)          :: pos_syn(3), vel_syn(3)

 ! Libration point coordinates (canonical units)
 real(8)          :: xL, yL

 ! Time settings
 real(8)          :: period_days, step_days
 integer, parameter :: nPeriod = 1      ! number of revolutions to write
 integer, parameter :: nSteps  = 100    ! steps per revolution
 integer          :: i, nEpoch

 ! Closure check
 real(8)          :: xi0, eta0, zeta0, xi0d, eta0d, zeta0d
 real(8)          :: xiT, etaT, zetaT, xiTd, etaTd, zetaTd

 ! File output
 integer          :: file_unit, io_stat
 character(len=256) :: output_file
 logical          :: success, dir_exists

 ! Earth-Moon system parameters
 real(8) :: m_earth, m_moon, a12

 !===========================================================================
 ! 1. BUILD EARTH-MOON CRTBP SYSTEM
 !===========================================================================
 m_earth = 5.9722d24        ! Earth mass (kg)
 m_moon  = 7.342d22         ! Moon mass (kg)
 a12     = 384400.0d0       ! Earth-Moon distance (km)

 call earth_moon_system%init("Earth-Moon", m_earth, m_moon, a12)

 !===========================================================================
 ! 2. INITIALIZE THE HALO ORBIT
 !===========================================================================
 LID   = LAGRANGE_L2
 beta  = 0.2d0             ! out-of-plane amplitude, in units of gamma
 phi0  = 0.0d0              ! initial phase

 ! Reference epoch (the CRTBP is autonomous; t0 fixes the calendar date)
 call t0%init(2025d0, 1d0, 1d0, 0d0, 0d0, 0d0)

 call orbit%init(earth_moon_system, LID, beta, phi0, t0, &
                 "Earth-Moon L2 Halo", success)

 if (.not. success) then
   write(*,*) 'Error: Halo orbit could not be initialized'
   stop
 endif

 call earth_moon_system%getLagrangePointCoord(LID, xL, yL)

 !===========================================================================
 ! 3. SUMMARY
 !===========================================================================
 write(*,*) '=============================================='
 write(*,*) 'Halo Orbit Ephemeris Generation'
 write(*,*) '=============================================='
 write(*,*) 'CRTBP system : ', trim(earth_moon_system%name)
 write(*,*) 'Mass ratio mu: ', earth_moon_system%mu
 write(*,*) 'Length unit  : ', earth_moon_system%a12, ' km'
 write(*,*) 'Time unit TU : ', earth_moon_system%TU, ' s (', &
            earth_moon_system%TU / 86400.0d0, ' days)'
 write(*,*) ''

 call orbit%show()

 period_days = orbit%period * earth_moon_system%TU / 86400.0d0
 step_days   = period_days / dble(nSteps)
 nEpoch      = nPeriod * nSteps

 write(*,*) ''
 write(*,*) 'Orbit amplitudes in km:'
 write(*,*) '  in-plane   alpha*gamma*a12 = ', &
            orbit%alpha * orbit%gamma * earth_moon_system%a12, ' km'
 write(*,*) '  vertical   beta *gamma*a12 = ', &
            orbit%beta  * orbit%gamma * earth_moon_system%a12, ' km'
 write(*,*) 'Time span    : ', period_days * nPeriod, ' days'
 write(*,*) 'Time step    : ', step_days, ' days'
 write(*,*) 'Number of epochs: ', nEpoch + 1
 write(*,*) ''

 !===========================================================================
 ! 4. GENERATE EPHEMERIS AND SAVE TO FILE
 !===========================================================================
 output_file = 'data/output/halo/halo_L2_ephemeris.txt'

 inquire(file='data/output/halo/', exist=dir_exists)
 if (.not. dir_exists) then
   call execute_command_line('mkdir "data/output/halo"', wait=.true.)
 endif

 open(newunit=file_unit, file=trim(output_file), status='replace', &
      action='write', iostat=io_stat)
 if (io_stat /= 0) then
   write(*,*) 'Error: Cannot open output file ', trim(output_file)
   stop
 endif

 ! Header (all lines start with '#')
 write(file_unit, '(A)') '# ============================================================================'
 write(file_unit, '(A)') '# HALO ORBIT EPHEMERIS'
 write(file_unit, '(A)') '# ============================================================================'
 write(file_unit, '(A,A)')   '# Orbit name: ', trim(orbit%name)
 write(file_unit, '(A,I2)')  '# Libration point: L', LID
 write(file_unit, '(A,F14.10)') '# Mass ratio mu      : ', earth_moon_system%mu
 write(file_unit, '(A,F14.3)')  '# Length unit (km)   : ', earth_moon_system%a12
 write(file_unit, '(A,F14.3)')  '# Time unit TU (s)   : ', earth_moon_system%TU
 write(file_unit, '(A,F14.10)') '# gamma              : ', orbit%gamma
 write(file_unit, '(A,F14.8)')  '# beta  (vertical,   units of gamma): ', orbit%beta
 write(file_unit, '(A,F14.8)')  '# alpha (in-plane,   units of gamma): ', orbit%alpha
 write(file_unit, '(A,F14.8)')  '# Orbit frequency (1/TU): ', orbit%omega
 write(file_unit, '(A,F14.6)')  '# Period (days)      : ', period_days
 write(file_unit, '(A)') '# ============================================================================'
 write(file_unit, '(A)') '# Coordinate frame: Earth-Moon synodic (rotating) frame'
 write(file_unit, '(A)') '#   origin: Earth-Moon barycenter'
 write(file_unit, '(A)') '#   x-axis: along the Earth -> Moon line'
 write(file_unit, '(A)') '#   z-axis: normal to the Earth-Moon orbital plane'
 write(file_unit, '(A)') '#   y-axis: completes the right-handed system'
 write(file_unit, '(A)') '# Velocity is the time derivative of the rotating-frame position'
 write(file_unit, '(A)') '# Third-order analytical solution (Richardson 1980), CRTBP model:'
 write(file_unit, '(A)') '# the ephemeris is a reference orbit, not a perturbed propagation.'
 write(file_unit, '(A)') '# ============================================================================'
 write(file_unit, '(A)') '# Units: km, km/s'
 write(file_unit, '(A)') '# Columns: MJD  x(km)  y(km)  z(km)  vx(km/s)  vy(km/s)  vz(km/s)'
 write(file_unit, '(A)') '# ============================================================================'

 write(*,*) 'Generating ephemeris...'
 write(*,*) ''

 do i = 0, nPeriod * nSteps
    call t%init(t0%val + i * step_days)

    ! Halo state relative to L2, gamma-normalized (getPV handles the
    ! conversion of the type_t epoch to canonical time)
    call orbit%getPV(t, xi, eta, zeta, xi_dot, eta_dot, zeta_dot)

    ! Convert to the barycentric synodic frame, km and km/s.
    ! The analytical solution is scaled by gamma (unit length r1 = gamma*a12):
    !   position: (xL + gamma*xi, gamma*eta, gamma*zeta) * a12
    !   velocity: gamma*(xi_dot, eta_dot, zeta_dot) * a12/TU
    pos_syn(1) = (xL + orbit%gamma*xi) * earth_moon_system%a12
    pos_syn(2) = orbit%gamma * eta        * earth_moon_system%a12
    pos_syn(3) = orbit%gamma * zeta       * earth_moon_system%a12
    vel_syn(1) = orbit%gamma * xi_dot     * earth_moon_system%a12 / earth_moon_system%TU
    vel_syn(2) = orbit%gamma * eta_dot    * earth_moon_system%a12 / earth_moon_system%TU
    vel_syn(3) = orbit%gamma * zeta_dot   * earth_moon_system%a12 / earth_moon_system%TU

    write(file_unit, '(F16.8, 6F16.6)') t%val, pos_syn, vel_syn

    ! Print the first few and the last epochs for inspection
    if (i <= 2 .or. i == nPeriod * nSteps) then
       write(*, '(F16.8, 6F16.6)') t%val, pos_syn, vel_syn
    endif
 end do

 close(file_unit)

 !===========================================================================
 ! 5. CLOSURE CHECK (one period)
 !===========================================================================
 ! (a) in canonical time: validates that the analytical orbit closes
 call orbit%compute_state(0d0,          xi0, eta0, zeta0, xi0d, eta0d, zeta0d)
 call orbit%compute_state(orbit%period, xiT, etaT, zetaT, xiTd, etaTd, zetaTd)

 write(*,*) ''
 write(*,*) '=== Closure check (one period) ==='
 write(*,'(A,ES16.6)') ' canonical time, max |position| difference: ', &
     max(abs(xiT-xi0), abs(etaT-eta0), abs(zetaT-zeta0))
 write(*,'(A,ES16.6)') ' canonical time, max |velocity| difference: ', &
     max(abs(xiTd-xi0d), abs(etaTd-eta0d), abs(zetaTd-zeta0d))

 ! (b) through the MJD epochs actually written: adds the day <-> TU round trip
 call orbit%getPV(t0, xi0, eta0, zeta0, xi0d, eta0d, zeta0d)
 call orbit%getPV(timeAdd(t0, period_days), xiT, etaT, zetaT, xiTd, etaTd, zetaTd)

 write(*,'(A,ES16.6)') ' via MJD epochs,  max |position| difference: ', &
     max(abs(xiT-xi0), abs(etaT-eta0), abs(zetaT-zeta0))

 write(*,*) ''
 write(*,*) 'Ephemeris saved to: ', trim(output_file)
 write(*,*) '  Total epochs: ', nPeriod * nSteps + 1

end program halo_orbit
