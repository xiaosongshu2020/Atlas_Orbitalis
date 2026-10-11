program orbit_prediction_ECOM
   !=======================================================================
   ! Purpose: Test program for ECOM solar radiation pressure model
   !          Define a GPS-like satellite, use ECOM SRP model,
   !          predict orbit for 1 day, and save trajectory data.
   !
   ! Author: Shawn Hu
   !=======================================================================
   use module_satellite
   use module_const
   use module_elements
   use module_time
   implicit none

   !--------------------------------------------------------------------------
   ! Satellite and orbit parameters
   !--------------------------------------------------------------------------
   type(type_satellite) :: sat
   type(type_t) :: t0, t1
   double precision :: kep(6)           ! Keplerian elements: a(km), e, i(deg), Omega(deg), omega(deg), M(deg)
   double precision :: nsg1(6)          ! Non-singular elements type 1
   double precision :: ECOM_params(9)   ! ECOM 9 parameters (normalized acceleration units)
   integer :: i, num_steps
   double precision :: step_sec, dt_days
   type(type_t), allocatable :: epochs(:)
   character(len=255) :: output_file

   !--------------------------------------------------------------------------
   ! Step 1: Define initial orbit (GPS-like satellite)
   !--------------------------------------------------------------------------
   write(*, '(A)') '===================================================='
   write(*, '(A)') '  Orbit Prediction with ECOM SRP Model'
   write(*, '(A)') '===================================================='

   ! GPS satellite: semi-major axis ~26560 km, MEO orbit
   kep(1) = 26560.0d0       ! a  - semi-major axis [km]
   kep(2) = 0.01d0          ! e  - eccentricity
   kep(3) = 55.0d0          ! i  - inclination [deg]
   kep(4) = 30.0d0          ! Omega - RAAN [deg]
   kep(5) = 60.0d0          ! omega - argument of perigee [deg]
   kep(6) = 90.0d0          ! M  - mean anomaly [deg]

   write(*, '(A)') 'Initial Orbit (Keplerian Elements):'
   write(*, '(A, F10.3, A)') '  Semi-major axis:   ', kep(1), ' km'
   write(*, '(A, F10.5)')    '  Eccentricity:      ', kep(2)
   write(*, '(A, F10.3)')    '  Inclination:       ', kep(3), ' deg'
   write(*, '(A, F10.3)')    '  RAAN:              ', kep(4), ' deg'
   write(*, '(A, F10.3)')    '  Arg. of Perigee:   ', kep(5), ' deg'
   write(*, '(A, F10.3)')    '  Mean Anomaly:      ', kep(6), ' deg'

   !--------------------------------------------------------------------------
   ! Step 2: Define ECOM parameters
   !--------------------------------------------------------------------------
   ! ECOM 9-parameter model, in the order the force model reads them
   ! (module_perturbationSRP%solar_pressure_perturbation; see type_opMod%ECOM):
   !   1 D0  2 Y0  3 B0  4 Bc  5 Bs  6 Dc  7 Ds  8 Yc  9 Ys
   !   F = (D0 + Dc*cos u + Ds*sin u)*e_D
   !     + (Y0 + Yc*cos u + Ys*sin u)*e_Y
   !     + (B0 + Bc*cos u + Bs*sin u)*e_B
   ! Units: normalized acceleration, 1 = ACC_UNIT = 9.798 m/s^2.
   ! e_D points from the Sun to the satellite, so D0 > 0 pushes the satellite
   ! away from the Sun - that is the physical SRP direction. (This value was
   ! negative while the table was being read in a different order.)
   ! Illustrative magnitudes for a GPS-like satellite: D0 ~ 1e-7 (direct solar
   ! radiation pressure), Dc, Ds ~ 1e-8 (D-direction periodic terms),
   ! Y0 ~ 1e-8 (Y-direction bias, related to satellite panel design),
   ! B0 ~ 1e-8 (B-direction bias).
   ECOM_params(1) =  1.0d-7    ! D0
   ECOM_params(2) =  1.0d-8    ! Y0
   ECOM_params(3) =  1.0d-8    ! B0
   ECOM_params(4) =  0.0d0     ! Bc
   ECOM_params(5) =  0.0d0     ! Bs
   ECOM_params(6) =  1.0d-8    ! Dc
   ECOM_params(7) =  1.0d-8    ! Ds
   ECOM_params(8) =  0.0d0     ! Yc
   ECOM_params(9) =  0.0d0     ! Ys

   write(*, '(A)') 'ECOM SRP Model Parameters (normalized units):'
   write(*, '(A, 3ES12.3)') '  D: D0, Dc, Ds = ', ECOM_params(1), ECOM_params(6:7)
   write(*, '(A, 3ES12.3)') '  Y: Y0, Yc, Ys = ', ECOM_params(2), ECOM_params(8:9)
   write(*, '(A, 3ES12.3)') '  B: B0, Bc, Bs = ', ECOM_params(3:5)

   !--------------------------------------------------------------------------
   ! Step 3: Set time parameters
   !--------------------------------------------------------------------------
   ! Start time: 2025-01-01 00:00:00 UTC
   call t0%init(2025d0, 1d0, 1d0, 0d0, 0d0, 0d0)
   ! End time: 2025-01-02 00:00:00 UTC (1 day later)
   call t1%init(2025d0, 1d0, 2d0, 0d0, 0d0, 0d0)

   step_sec = 60.0d0  ! Output every 60 seconds
   dt_days = timeMinus(t1, t0)  ! Duration in days
   num_steps = int(dt_days * 86400.0d0 / step_sec) + 1

   write(*, '(A)') 'Time Parameters:'
   write(*, '(A, F10.3)') '  Start time: 2025-01-01 00:00:00 UTC (MJD: ', t0%val, ')'
   write(*, '(A, F10.3)') '  End time:   2025-01-02 00:00:00 UTC (MJD: ', t1%val, ')'
   write(*, '(A, F10.1, A)') '  Duration:    ', dt_days * 24.0d0, ' hours'
   write(*, '(A, F10.1, A)') '  Step:        ', step_sec, ' seconds'
   write(*, '(A, I6)')       '  Total steps: ', num_steps

   !--------------------------------------------------------------------------
   ! Step 4: Set force model switches
   !   (satellite 初始化后直接设置命名开关)
   !--------------------------------------------------------------------------

   !--------------------------------------------------------------------------
   ! Step 5: Initialize satellite
   !--------------------------------------------------------------------------
   ! Convert Keplerian elements to normalized units
   call kpl_SI2normal(kep, kep)
   ! Convert to non-singular elements type 1 (for numerical propagation)
   call kpl2nsg1(kep, nsg1)

   ! Initialize satellite with basic parameters
   call sat%init(t0, nsg1, &
                 number = 1, &
                 name = "GPS-ECOM-Test", &
                 mass = 500.0d0, &
                 eleType = STATE_TYPE_NONSINGULAR_1, &
                 propagator = "numerical", &
                 degree = 8, order = 8)
   sat%opMod%key_moonGravity = 1   ! Moon gravity (third-body)
   sat%opMod%key_sunGravity  = 1   ! Sun gravity (third-body)
   sat%opMod%key_srp         = 1   ! Solar radiation pressure (SRP) - ENABLED
   ! key_drag = 0       ! Atmospheric drag (disabled for MEO orbit)
   ! key_postNewton = 0 ! Post-Newtonian effect (disabled)

   write(*, '(A)') 'Force Model:'
   write(*, '(A, I2, A, I2)') '  Earth Gravity degree/order: ', sat%opMod%degree, ' x ', &
        &   sat%opMod%order
   write(*, '(A, I2)') '  key_moonGravity:       ', sat%opMod%key_moonGravity
   write(*, '(A, I2)') '  key_sunGravity:        ', sat%opMod%key_sunGravity
   write(*, '(A, I2)') '  key_srp (ECOM Model):  ', sat%opMod%key_srp

   ! Manually set ECOM SRP model parameters (sat_init does not support these directly)
   sat%opMod%SRP_Model = 2       ! 1 = cannon ball, 2 = ECOM (9 parameters)
   sat%opMod%ECOM = ECOM_params  ! Set ECOM 9 parameters
   sat%opMod%reflectivity = 0.0d0 ! Not used in ECOM model, but set to 0
   sat%opMod%amrSRP = 0.02d0     ! Area-to-mass ratio for SRP [m^2/kg]

   ! Display satellite information
   call sat%show()

   !--------------------------------------------------------------------------
   ! Step 6: Generate time epochs
   !--------------------------------------------------------------------------
   allocate(epochs(num_steps))
   epochs(1) = t0
   do i = 2, num_steps
      epochs(i)%day = t0%day
      epochs(i)%sec = t0%sec + (i - 1) * step_sec
      ! Normalize seconds to [0, 86400)
      do while (epochs(i)%sec >= 86400.0d0)
         epochs(i)%day = epochs(i)%day + 1
         epochs(i)%sec = epochs(i)%sec - 86400.0d0
      end do
      epochs(i)%val = epochs(i)%day + epochs(i)%sec / 86400.0d0
   end do

   !--------------------------------------------------------------------------
   ! Step 7: Generate trajectory
   !--------------------------------------------------------------------------
   write(*, '(A)') 'Generating trajectory with ECOM SRP model...'
   call sat%generateTrajectory(epochs)

   !--------------------------------------------------------------------------
   ! Step 8: Save trajectory to file
   !--------------------------------------------------------------------------
   output_file = 'data/output/orbit_ECOM_prediction.txt'
   write(*, '(A)') 'Saving trajectory to: ' // trim(output_file)
   ! ephType = 2: save as position/velocity (Cartesian) in km, km/s
   call sat%savetrajectory(STATE_TYPE_CARTESIAN, output_file)

   !--------------------------------------------------------------------------
   ! Step 9: Output summary
   !--------------------------------------------------------------------------
   write(*, '(A)') '===================================================='
   write(*, '(A)') '  Prediction Complete - Summary'
   write(*, '(A)') '===================================================='
   write(*, '(A, I6)')    '  Total epochs:     ', sat%nEph
   write(*, '(A, A)')     '  Output file:      ', trim(output_file)

   ! Show initial and final states
   if (sat%nEph > 0) then
      write(*, '(A)') '  Initial State (Cartesian):'
      write(*, '(A, 3F15.3)') '    Position (km):   ', sat%trajectory(1)%elements(1:3)
      write(*, '(A, 3F15.6)') '    Velocity (km/s): ', sat%trajectory(1)%elements(4:6)

      write(*, '(A)') '  Final State (Cartesian):'
      write(*, '(A, 3F15.3)') '    Position (km):   ', sat%trajectory(sat%nEph)%elements(1:3)
      write(*, '(A, 3F15.6)') '    Velocity (km/s): ', sat%trajectory(sat%nEph)%elements(4:6)
   endif

   write(*, '(A)') '===================================================='
   write(*, '(A)') '  Test Completed Successfully!'
   write(*, '(A)') '===================================================='

   ! Clean up
   deallocate(epochs)

end program orbit_prediction_ECOM
