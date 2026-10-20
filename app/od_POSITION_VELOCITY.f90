!**************************************************************************
! Program Name: od_POSITION_VELOCITY
! Purpose: Batch orbit determination from position/velocity data
!          (OBS_TYPE_POSITION_VELOCITY). Cartesian state + amrSRP, 7 params.
!**************************************************************************
program od_position_velocity
use module_access, only : OBS_TYPE_POSITION_VELOCITY
use od_demo_common
implicit none

call runODdemo(OBS_TYPE_POSITION_VELOCITY, "POSITION_VELOCITY")

end program od_position_velocity
