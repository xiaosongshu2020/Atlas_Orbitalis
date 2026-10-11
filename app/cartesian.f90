! Demonstration of module_cartesian
! 2025-12-20 by Shawn Hu
program demo_cartesian
 use module_cartesian
 implicit none
 
 type(type_cartesian) :: cart
 type(type_t) :: time1
 real*8 :: pos(3), vel(3)
 real*8 :: longitude, latitude
 
 ! 1. Create time object
 call time1%init(51544.5d0)  ! MJD for 2000-01-01 12:00:00
 
 ! 2. Initialize cartesian object with position, velocity, reference frame and time
 pos = [7000.0d0, 0.0d0, 0.0d0] * u_km  ! km -> normalized unit
 vel = [0.0d0, 7.5d0, 0.0d0] * u_km / u_sec     ! km/s -> normalized unit
 call cart%init(pos, vel, COORD_ICRS, time1)
 
 ! 3. Display object properties
 print*, "Initial cartesian object:"
 print*, "  Position (km):", cart%position / u_km
 print*, "  Velocity (km/s):", cart%velocity / (u_km / u_sec)
 print*, "  Reference frame:", cart%getReferenceName()
 print*, "  Time:", cart%t%toDate()
 
 ! 4. Demonstrate coordinate transformations
 ! ICRS -> ITRS (uses internal time)
 call cart%ICRS2ITRS()
 print*, "After ICRS2ITRS:"
 print*, "  Position (km):", cart%position / u_km
 print*, "  Velocity (km/s):", cart%velocity / (u_km / u_sec)
 print*, "  Reference frame:", cart%getReferenceName()
 print*, "  Time:", cart%t%toDate()
 
 ! ITRS -> ICRS (uses internal time)
 call cart%ITRS2ICRS()

 ! ITRS -> Horizon (need to be in ITRS first)
 call cart%ICRS2ITRS()
 longitude = 116.4d0 * deg  ! Beijing longitude (deg -> radians)
 latitude = 39.9d0 * deg    ! Beijing latitude (deg -> radians)
 call cart%ITRS2Horizon(longitude, latitude)

 ! Horizon -> ITRS
 call cart%Horizon2ITRS(longitude, latitude)

 print*, "Demo completed successfully!"
 
end program demo_cartesian
