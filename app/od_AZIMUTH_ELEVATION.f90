!**************************************************************************
! Program Name: od_AZIMUTH_ELEVATION
! Purpose: Batch orbit determination from azimuth/elevation data
!          (OBS_TYPE_AZIMUTH_ELEVATION). Cartesian state + amrSRP, 7 params.
! Note: module_ODor converts Az/El to RA/Dec for the batch estimator; the
!       input observations and the saved *_obs.txt file are Az/El.
!**************************************************************************
program od_azimuth_elevation
use module_access, only : OBS_TYPE_AZIMUTH_ELEVATION
use module_ODdemo
implicit none

call runODdemo(OBS_TYPE_AZIMUTH_ELEVATION, "AZIMUTH_ELEVATION")

end program od_azimuth_elevation
