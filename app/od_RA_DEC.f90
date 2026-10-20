!**************************************************************************
! Program Name: od_RA_DEC
! Purpose: Batch orbit determination from right-ascension/declination data
!          (OBS_TYPE_RA_DEC). Cartesian state + amrSRP, 7 parameters.
! The shared scenario and all output handling live in module_ODdemo.
!**************************************************************************
program od_ra_dec
use module_access, only : OBS_TYPE_RA_DEC
use module_ODdemo
implicit none

call runODdemo(OBS_TYPE_RA_DEC, "RA_DEC")

end program od_ra_dec
