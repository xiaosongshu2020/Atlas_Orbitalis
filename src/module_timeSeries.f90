
!================================================================================
! MODULE: module_timeSeries
!
! DESCRIPTION:
!   This module defines the type_TimeSeries class for managing time series data
!   in satellite orbit analysis and other time-dependent simulations.
!
! AUTHOR: Shawn Hu
! LASTEST MODIFIED: 2025-12-22
!
! TYPE: type_TimeSeries
!   A container for time series data with evenly spaced time points.
!
! ATTRIBUTES:
!   n        : integer        - Number of time points
!   interval : double precision         - Time interval between points (seconds)
!   t0       : type(type_t)   - Initial time point
!   tEnd     : type(type_t)   - Final time point
!   duration : double precision         - Total duration of time series (days)
!   time     : type(type_t), allocatable(:) - Array of time points (UTC)
!
! METHODS:
!   init  : Initialize time series with start time, step size, and number of points
!   free  : Deallocate memory and reset object
!   show  : Display time series information to standard output
!
! DEPENDENCIES:
!   module_time - For time manipulation and type_t definition
!
! USAGE EXAMPLE:
!   type(type_TimeSeries) :: ts
!   type(type_t) :: start_time
!   call start_time%init(2025d0, 12d0, 22d0, 0d0, 0d0, 0d0)  ! Initialize start time
!   call ts%init(start_time, 60.0d0, 100)  ! 60-second intervals, 100 points
!   call ts%show() ! Display time series information
!   call ts%free() ! free memory
!
!================================================================================

module module_timeSeries
use module_time
    
type :: type_TimeSeries
    integer :: n = 0           ! Number of time points
    double precision :: interval         ! Time interval between points (seconds)
    type(type_t) :: t0         ! Initial time point
    type(type_t) :: tEnd       ! Final time point
    double precision :: duration         ! Total duration (days)
    type(type_t), allocatable :: time(:) ! Array of time points (UTC, MJD format)
contains 
    procedure :: init => timeSeries_initialize
    procedure :: free => timeSeries_free
    procedure :: show => timeSeries_show
end type type_TimeSeries

contains
    
    
!==============================================================================
! SUBROUTINE: timeSeries_initialize
!
! DESCRIPTION:
!   Initialize a timeSeries object with evenly spaced time points.
!
! INPUTS:
!   this  : type_timeSeries object (implicit)
!   t0    : type(type_t) - Starting time point
!   step  : double precision       - Time interval between points (seconds)
!   n     : integer      - Number of time points to generate
!
! OUTPUTS:
!   this  : Initialized timeSeries object with allocated time array
!
! ALGORITHM:
!   1. Deallocate existing time array if allocated
!   2. Allocate time array with size n
!   3. Generate time points: time(i) = t0 + (i-1)*step/86400.0 days
!   4. Set object attributes: n, interval, t0, tEnd, duration
!
! NOTES:
!   - Time points are stored in type_t format (UTC, MJD)
!   - Duration is calculated as (n-1)*step/86400.0 days
!   - Step is converted from seconds to days using 86400.0 seconds/day
!
!==============================================================================
subroutine timeSeries_initialize(this, t0, step, n)
implicit none
class (type_timeSeries) :: this
type(type_t), intent(in) :: t0
double precision, intent(in) :: step
integer, intent(in) :: n
integer i

! Deallocate existing time array if allocated
if (allocated(this%time)) deallocate(this%time) 

! Allocate and generate time points
allocate(this%time(n))
do i = 1, n
    ! Add (i-1)*step seconds to initial time, convert seconds to days
    this%time(i) = timeAdd(t0, step*dble(i-1)/86400d0)
enddo

! Set object attributes
this%n = n
this%interval = step
this%t0 = this%time(1)
this%tEnd = this%time(n)
this%duration = step*dble(n-1)/86400d0
endsubroutine

!==============================================================================
! SUBROUTINE: timeSeries_free
!
! DESCRIPTION:
!   Deallocate memory and reset a timeSeries object.
!
! INPUTS:
!   this : type_timeSeries object (implicit)
!
! OUTPUTS:
!   this : Cleaned timeSeries object with n=0 and deallocated time array
!
! NOTES:
!   - Always check allocation status before deallocating
!   - Resets n to 0 to indicate empty object
!
!==============================================================================
subroutine timeSeries_free(this)
class(type_timeSeries), intent(in out) :: this
if (allocated(this%time)) deallocate(this%time)
this%n = 0
endsubroutine
    
!==============================================================================
! SUBROUTINE: timeSeries_show
!
! DESCRIPTION:
!   Display time series information to standard output.
!
! INPUTS:
!   this : type_timeSeries object (implicit)
!
! OUTPUT FORMAT:
!   Time Series Info:
!   Number of times: [n]
!   Time interval (sec): [interval]
!   Start time: [t0]
!   End time: [tEnd]
!   Duration (days): [duration]
!
!==============================================================================
subroutine timeSeries_show(this)
implicit none
class(type_timeSeries), intent(in) :: this
print *, "Time Series Info:"
print *, "Number of times: ", this%n
print *, "Time interval (sec): ", this%interval
print *, "Start time: ", this%t0
print *, "End time: ", this%tEnd
print *, "Duration (days): ", this%duration
endsubroutine
    
endmodule
