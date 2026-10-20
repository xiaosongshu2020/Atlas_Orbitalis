!**************************************************************************
! Program Name: od_RANGING
! Purpose: Batch orbit determination from range-only data
!          (OBS_TYPE_RANGING). Cartesian state + amrSRP, 7 parameters.
!**************************************************************************
program od_ranging
use module_access, only : OBS_TYPE_RANGING
use module_ODdemo
implicit none

call runODdemo(OBS_TYPE_RANGING, "RANGING")

end program od_ranging
