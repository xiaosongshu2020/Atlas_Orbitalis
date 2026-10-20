!**************************************************************************
! Program Name: od_RANGERATE
! Purpose: Batch orbit determination from range-rate-only data
!          (OBS_TYPE_RANGERATE). Cartesian state + amrSRP, 7 parameters.
!**************************************************************************
program od_rangerate
use module_access, only : OBS_TYPE_RANGERATE
use od_demo_common
implicit none

call runODdemo(OBS_TYPE_RANGERATE, "RANGERATE")

end program od_rangerate
