!**************************************************************************
! Program Name: od_MIXED
! Purpose: Batch orbit determination from MIXED observations:
!          RA/Dec (angle) + RANGING (range) + RANGERATE (range-rate),
!          all from the same MEO satellite and Beijing station.
!          Cartesian state + amrSRP (7 parameters).
! The shared scenario lives in od_demo_common.f90.
!**************************************************************************
program od_mixed
use od_demo_common
implicit none

call runMixedODdemo()

end program od_mixed
