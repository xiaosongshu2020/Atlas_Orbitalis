! Test program for moon coordinate transformation - Earth position calculation
! Calculate Earth position in Moon-centered inertial frame and convert to lunar coordinates
! Author: Cline
! Date: 2025/10/29

program test_moon_coordinate_transformation2
   use module_space
   use module_const
   use module_time
   implicit none
   
   ! Variable declarations
   real*8 :: t_utc, t_TT           ! TT time (JD)
   real*8 :: t_start, t_end         ! Start and end times in MJD
   real*8 :: t_current              ! Current time in MJD
   real*8 :: time_step              ! Time step in hours
   real*8 :: earth_pos_helio(6)     ! Earth position in heliocentric frame (AU)
   real*8 :: moon_pos(6)      ! Moon position in heliocentric frame (AU)
   real*8 :: earth_pos_moon(6)      ! Earth position relative to Moon (AU)
   real*8 :: earth_pos_moon_km(3)   ! Earth position relative to Moon (km)
   real*8 :: earth_pos_mci(3)       ! Earth position in Moon-centered inertial frame (km)
   real*8 :: earth_pos_mcf(3)       ! Earth position in Moon-centered fixed frame (km)
   real*8 :: longitude, latitude, altitude  ! Lunar coordinates
   integer :: i, n_steps
   character(len=100) :: output_file
   integer :: file_unit
   integer :: iPlephStatus

   ! Initialize parameters
   time_step = 2.0d0  ! 2 hours
   output_file = 'data/output/EarthPosition.txt'
   
   ! Calculate start and end times in MJD
   ! 2025-01-01 00:00:00 to 2025-07-01 00:00:00
   call date2MJD(2010.0d0, 1.0d0, 1.0d0, 0.0d0, t_start, t_current)
   call date2MJD(2010.0d0, 7.0d0, 1.0d0, 0.0d0, t_end, t_current)
   
   ! Convert MJD to JD for PLEPH function
   t_start = t_start + 2400000.5d0
   t_end = t_end + 2400000.5d0
   
   ! Calculate number of steps
   n_steps = int((t_end - t_start) * 24.0d0 / time_step) + 1
   
   write(*,*) "================================================"
   write(*,*) "Earth Position in Moon Coordinate System"
   write(*,*) "================================================"
   write(*,*)
   write(*,*) "Time span: 2025-01-01 to 2025-07-01"
   write(*,*) "Time step:", time_step, "hours"
   write(*,*) "Number of steps:", n_steps
   write(*,*)
   
   ! Open output file
   open(newunit=file_unit, file=output_file, status='replace', action='write')
   write(file_unit, '(A)') "Time(JD)          X(km)            Y(km)            "//&
   "Z(km)            Longitude(deg)   Latitude(deg)   Altitude(km)"
   
   ! Main loop over time
   do i = 0, n_steps-1
       t_current = t_start + (i * time_step) / 24.0d0
       
       ! Convert UTC to TT
       t_utc = t_current - 2400000.5d0
       call UTC2TT(t_utc, t_TT)
       t_TT = t_TT + 2400000.5d0
       
       ! Get Moon position relative to Sun (heliocentric)
       ! NTARG=10 (Moon), NCENT=11 (Sun)
       call pleph(t_TT, 10, 3, moon_pos)
       ! pleph now returns zeros instead of aborting on a bad ephemeris --
       ! fail loudly here rather than writing meaningless coordinates
       call pleph_status(iPlephStatus)
       if (iPlephStatus /= 0) then
           write(0, *) "Error @test_moon_coordinate_transformation2: JPL ", &
               &   "ephemeris unavailable (pleph status ", iPlephStatus, ")."
           write(0, *) "  Run from the repository root, or restore the JPL file."
           error stop
       endif

       ! Calculate Earth position relative to Moon
       earth_pos_moon = - moon_pos
       
       ! Convert from AU to km
       earth_pos_moon_km = earth_pos_moon(1:3) * ASTRONOMICAL_UNIT_KM
       
       ! Transform to Moon-centered inertial frame
       ! Note: PLEPH already gives positions in J2000 inertial frame
       ! So earth_pos_moon_km is already in Moon-centered inertial frame
       
       ! For this test, we'll use the position directly as MCI
       earth_pos_mci = earth_pos_moon_km
       
       ! Transform from MCI to MCF (Moon-centered fixed)
       call LCI2LCF(t_TT, earth_pos_mci, earth_pos_mcf)
       
       ! Convert Cartesian coordinates to lunar longitude, latitude, altitude
       call cartesian2lunar(earth_pos_mcf, longitude, latitude, altitude)
       
       ! Convert from radians to degrees
       longitude = longitude * rad2deg
       latitude = latitude * rad2deg
       
       ! Write results to file
       write(file_unit, '(F16.8, 6F16.8)') t_current, earth_pos_mcf, longitude, latitude, altitude
       
       ! Progress indicator
       if (mod(i, 100) == 0) then
           write(*, '(A, I6, A, I6)') "Progress: ", i, "/", n_steps
       endif
   end do
   
   ! Close output file
   close(file_unit)
   
   write(*,*)
   write(*,*) "Calculation completed successfully!"
   write(*,*) "Results saved to: ", trim(output_file)
   write(*,*)
   write(*,*) "================================================"
   write(*,*) "Program completed"
   write(*,*) "================================================"
   
contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Convert Cartesian coordinates to lunar longitude, latitude, altitude
   ! Input: position in Moon-centered fixed frame (km)
   ! Output: longitude, latitude (radians), altitude (km)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine cartesian2lunar(position, longitude, latitude, altitude)
       implicit none
       real*8, intent(in) :: position(3)
       real*8, intent(out) :: longitude, latitude, altitude
       real*8 :: r, x, y, z
       real*8, parameter :: R_moon = 1737.4d0  ! Moon radius in km
       
       x = position(1)
       y = position(2)
       z = position(3)
       r = sqrt(x*x + y*y + z*z)
       
       ! Calculate longitude (0 to 2*pi)
       longitude = atan2(y, x)
       if (longitude < 0.0d0) then
           longitude = longitude + 2.0d0 * pi
       endif
       
       ! Calculate latitude (-pi/2 to pi/2)
       latitude = asin(z / r)
       
       ! Calculate altitude (distance from Moon surface)
       altitude = r - R_moon
       
   end subroutine cartesian2lunar

end program test_moon_coordinate_transformation2
