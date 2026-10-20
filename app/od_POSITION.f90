!**************************************************************************
! Program Name: od_POSITION
! Purpose: Batch orbit determination from position-only data
!          (OBS_TYPE_POSITION). Cartesian state + amrSRP, 7 parameters.
!**************************************************************************
program od_position
use module_access, only : OBS_TYPE_POSITION
use module_ODdemo
implicit none

call runODdemo(OBS_TYPE_POSITION, "POSITION")

end program od_position
