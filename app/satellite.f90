!! Demonstration of module_satellite
!! 2025-12-21 by Cline (assistant)
!!
!! This example demonstrates the usage of the satellite object defined in
!! module_satellite.f90. It covers various functionalities including
!! initialization, trajectory generation, file I/O, state query, and propagation.
!!
!! Key features demonstrated:
!! 1. Satellite initialization with Keplerian elements
!! 2. Information display using show() method
!! 3. Trajectory generation with numerical propagator (including perturbations)
!! 4. Trajectory generation with two-body propagator
!! 5. Saving and loading trajectory files
!! 6. State interpolation at arbitrary times
!! 7. State propagation to new epochs
!! 8. Thruster model integration


program demo_satellite
 use module_satellite
 implicit none
 
 ! Declare satellite objects
 type(type_satellite) :: sat1, sat2, sat3
 type(type_t) :: t0, t1, t_interp
 double precision :: kepler0(6), kepler1(6), cartesian0(6)
 type(type_t), allocatable :: epochs(:)
 double precision :: interpolated_state(6), interpolated_pv(6)
 integer :: i, n_points
 character(len=100) :: output_file1, output_file2, output_file3
 
 ! Output file paths
 output_file1 = "data/output/satellite_traj_numerical.txt"
 output_file2 = "data/output/satellite_traj_twobody.txt"
 output_file3 = "data/output/satellite_traj_loaded.txt"

 !------------------------------------------------------------------------
 ! EXAMPLE 1: Basic initialization and information display
 !------------------------------------------------------------------------

 ! Create initial epoch: 2025-01-01 00:00:00 UTC
 call t0%init(2025.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0)
 
 ! Define Keplerian elements (a in km, angles in degrees)
 ! a = 7000 km, e = 0.01, i = 45°, Ω = 30°, ω = 60°, M = 90°
 kepler0(1) = 7000.0d0 * u_km          ! Semi-major axis [km -> normalized]
 kepler0(2) = 0.01d0                 ! Eccentricity [dimensionless]
 kepler0(3) = 45.0d0 * deg           ! Inclination [deg -> rad]
 kepler0(4) = 30.0d0 * deg           ! RAAN [deg -> rad]
 kepler0(5) = 60.0d0 * deg           ! Argument of perigee [deg -> rad]
 kepler0(6) = 90.0d0 * deg           ! Mean anomaly [deg -> rad]
 
 ! Initialize satellite with Keplerian elements
 call sat1%init(t0, kepler0, &
                number=1, &
                name="DemoSat-1", &
                mass=500.0d0, &
                eleType=STATE_TYPE_KEPLERIAN, &
                propagator="numerical", &
                degree=8, order=8, & ! Gravity field degree/order
                Cd=2.2d0, &          ! Drag coefficient
                amrDrag=0.02d0, &    ! Area-to-mass ratio for drag
                F107=130.0d0, &      ! Solar flux index
                Ap=15.0d0)           ! Geomagnetic index
 sat1%opMod%key_moonGravity = 1   ! Moon third-body
 sat1%opMod%key_sunGravity  = 1   ! Sun third-body
 sat1%opMod%key_srp         = 1   ! solar radiation pressure
 sat1%opMod%key_drag        = 1   ! atmospheric drag
 
 ! Display satellite information
 call sat1%show() ! print satellite information
 
 !------------------------------------------------------------------------
 ! EXAMPLE 2: Trajectory generation with numerical propagator
 !------------------------------------------------------------------------

 ! Create array of output epochs (every hour for 24 hours)
 n_points = 25  ! 0 to 24 hours inclusive
 allocate(epochs(n_points))
 
 do i = 1, n_points
   call epochs(i)%init(2025.0d0, 1.0d0, 1.0d0, real(i-1, 8), 0.0d0, 0.0d0)
 end do
 
 print*, "Generating trajectory with numerical propagator..."
 print*, "  Start time: ", epochs(1)%toDate()
 print*, "  End time:   ", epochs(n_points)%toDate()
 print*, "  Number of points: ", n_points
 
 ! Generate trajectory
 call sat1%generateTrajectory(epochs)
 
 print*, "Trajectory generated successfully."
 print*, "Number of ephemeris points: ", sat1%nEph
 
 ! Save trajectory to file (in Keplerian format)
 print*, "Saving trajectory to: ", trim(output_file1)
 call sat1%savetrajectory(STATE_TYPE_KEPLERIAN, output_file1, freeEph=.false.)
 
 ! Display final state
 print*, "Final Keplerian elements:"
 kepler1 = sat1%trajectory(sat1%nEph)%elements
 call kpl_normal2SI(kepler1, kepler1)  ! Convert to SI units for display
 print*, "  a = ", kepler1(1), " km"
 print*, "  e = ", kepler1(2)
 print*, "  i = ", kepler1(3), " deg"
 print*, "  Omega = ", kepler1(4), " deg"
 print*, "  omega = ", kepler1(5), " deg"
 print*, "  M = ", kepler1(6), " deg"
 
 !------------------------------------------------------------------------
 ! EXAMPLE 3: Trajectory generation with two-body propagator
 !------------------------------------------------------------------------

 ! Initialize a second satellite with two-body propagator
 call sat2%init(t0, kepler0, &
                number=2, &
                name="DemoSat-2", &
                mass=500.0d0, &
                eleType=STATE_TYPE_KEPLERIAN, &
                propagator="two-body")
   
 ! Generate trajectory with two-body propagator
 call sat2%generateTrajectory(epochs)
 
 print*, "Two-body trajectory generated."
 print*, "Number of ephemeris points: ", sat2%nEph
 
 ! Save trajectory to file
 print*, "Saving trajectory to: ", trim(output_file2)
 call sat2%savetrajectory(STATE_TYPE_KEPLERIAN, output_file2, freeEph=.false.)
 
 !------------------------------------------------------------------------
 ! EXAMPLE 4: Loading trajectory from file
 !------------------------------------------------------------------------

 ! Initialize a third satellite
 call sat3%init(t0, kepler0, &
                number=3, &
                name="DemoSat-3", &
                mass=500.0d0, &
                eleType=STATE_TYPE_KEPLERIAN, &
                propagator="numerical")
 
 ! Load trajectory from the file saved earlier
 print*, "Loading trajectory from: ", trim(output_file1)
 call sat3%readtrajectory(STATE_TYPE_KEPLERIAN, output_file1)
 
 print*, "Trajectory loaded successfully."
 print*, "Number of loaded ephemeris points: ", sat3%nEph
 
 ! Save the loaded trajectory to a new file for verification
 print*, "Saving loaded trajectory to: ", trim(output_file3)
 call sat3%savetrajectory(STATE_TYPE_KEPLERIAN, output_file3, freeEph=.false.)
 
 !------------------------------------------------------------------------
 ! EXAMPLE 5: State interpolation at arbitrary time
 !------------------------------------------------------------------------

 ! Create a time for interpolation (2025-01-01 12:30:00 UTC)
 call t_interp%init(2025.0d0, 1.0d0, 1.0d0, 12.0d0, 30.0d0, 0.0d0)
 
 ! Get orbital elements at interpolated time
 interpolated_state = sat1%getStates(t_interp)
 
 ! Convert to SI units for display
 call kpl_normal2SI(interpolated_state, interpolated_state)
 
 print*, "Interpolated state at ", t_interp%toDate()
 print*, "  a = ", interpolated_state(1), " km"
 print*, "  e = ", interpolated_state(2)
 print*, "  i = ", interpolated_state(3), " deg"
 print*, "  Omega = ", interpolated_state(4), " deg"
 print*, "  omega = ", interpolated_state(5), " deg"
 print*, "  M = ", interpolated_state(6), " deg"
 
 ! Get Cartesian position/velocity at interpolated time
 interpolated_pv = sat1%getPV(t_interp)
 call cts_normal2SI(interpolated_pv, interpolated_pv)
 
 print*, "Interpolated Cartesian state:"
 print*, "  Position (km): ", interpolated_pv(1:3)
 print*, "  Velocity (km/s): ", interpolated_pv(4:6)
 
 !------------------------------------------------------------------------
 ! EXAMPLE 6: State propagation to new epoch
 !------------------------------------------------------------------------

 ! Create a new epoch (1 day later)
 call t1%init(2025.0d0, 1.0d0, 2.0d0, 0.0d0, 0.0d0, 0.0d0)
 
 print*, "Propagating satellite state from:"
 print*, "  ", sat1%states%t%toDate()
 print*, "  to: ", t1%toDate()
 
 ! Propagate satellite state
 call sat1%propagate(t1)
 
 print*, "Propagation completed."
 print*, "New satellite state at: ", sat1%states%t%toDate()
 
 ! Display new orbital elements
 kepler1 = sat1%states%elements
 call kpl_normal2SI(kepler1, kepler1)
 print*, "  a = ", kepler1(1), " km"
 print*, "  e = ", kepler1(2)
 print*, "  i = ", kepler1(3), " deg"
 print*, "  Omega = ", kepler1(4), " deg"
 print*, "  omega = ", kepler1(5), " deg"
 print*, "  M = ", kepler1(6), " deg"

 !------------------------------------------------------------------------
 ! Cleanup and completion
 !------------------------------------------------------------------------
 
 ! Deallocate arrays
 if (allocated(epochs)) deallocate(epochs)
 
 ! Free trajectory memory
 call sat1%freeTrajectory()
 call sat2%freeTrajectory()
 call sat3%freeTrajectory()

end program demo_satellite
