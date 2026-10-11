

 program  station_satellite_angle_ekf

   use module_Atlas

   implicit none
   
   ! Define variables
   integer :: i, j, nEpochs
   type(type_t) :: currentEpoch, t0, te, sat_t0
   character(len=200) :: outputDir
   logical :: dirExists
   real*8 :: interval, sat_x0(6)
   integer :: nStep
   real*8 :: errStd_ele(6), errStd_amrDrag, errStd_amrSRP
   type(type_access) :: acc(10)
   character*100 :: fn_obsData
   real*8 :: residuals(2)
   
   type(type_satellite) :: sat
   type(type_station) :: station
   type(type_access) :: access
   type(type_satelliteScenario) :: ss
   type(type_ODor) :: ODor
   
   type(type_varGroup) :: groups(1)
   real*8, allocatable :: P(:, :)
   
   character(len=200) :: introduction
   
   print*, " o_o Station-Satellite Angle EKF Orbit Determination"
   
   ! 1. Initialize satellite scenario

   ! Initialize time series
   call t0%init(60306d0, 3600d0*4d0)  ! Initial time, 2023-12-28 04:00:00
   call te%init(60316d0, 3600d0*4d0)
   interval = 30d0  ! Time interval, seconds
   nStep = 2*1440*10     ! Number of time series points
   call ss%init(t0, te, interval, eleType = STATE_TYPE_CARTESIAN)
   

   
   ! Initialize satellite object
   sat_t0 = t0
   ! Satellite initial state (Cartesian coordinates, km and km/s)
   sat_x0 = (/20.000000d0, -42000.000000d0, 20d0, 3.08d0, 0.0d0, 0.10d0/)
   ! Convert to normalized units
   call cts_SI2normal(sat_x0, sat_x0)
   ! Initialize satellite
   call sat%init(sat_t0, sat_x0, &
                      number=1, &
                      name='HY-2A', &
                      eleType=STATE_TYPE_CARTESIAN, &
                      propagator='numerical', &
                      Cd=2.2d0, &
                      amrDrag=0.02d0, &
                      F107=130d0, &
                      Ap=15d0, &
                      reflectivity=0.0d0, &
                      amrSRP=0.02d0, &
                      degree=5, order=5)
   sat%opMod%key_moonGravity = 1   ! Moon third-body
   sat%opMod%key_sunGravity  = 1   ! Sun third-body
   sat%opMod%key_srp         = 1   ! solar radiation pressure
   sat%opMod%usePV = 1
    ! Add satellite
   call ss%addSatellite(sat)
   
   
   
   ! Initialize ground station object
   call station%init(number=1, &
                          name='Sanya', &
                          longitude=108d0, &
                          latitude=18d0, &
                          altitude=7d0)
   ! generate station's PV eph in J2000
   call station%generateTrajectory(ss%ts)
   ! Add ground station
   call ss%addStation(station)
   

   
   ! Initialize access link object
   call access%init(number=1, &
                         name='access_sat1_station1', &
                         fromType=OBJECT_TYPE_SATELLITE, &
                         fromNum=1, &
                         toType=OBJECT_TYPE_STATION, &
                         toNum=1, &
                         obsType=OBS_TYPE_ALPHA_DELTA, &
                         minEle=10d0)
   access%errStd = 5d-6
   key_light_time_delay = .true.
   ! Add access link
   call ss%addAccess(access)
   
   ! Report access relationships
   call ss%reportAccess(.true.)

   ! 6. Initialize orbit determinator

   ! Set observation data file
   fn_obsData = './data/input/obsData_alphaDelta1.txt'

   ! Set variable group parameters (state)
   Groups(1)%nVar = 6
   Groups(1)%label = 1  ! states only
   Groups(1)%origin = 1 ! satellite 1
   allocate(Groups(1)%idx(6))
   Groups(1)%idx = (/1, 2, 3, 4, 5, 6/)

   ! Set initial covariance (override the zero matrix from init)
   errStd_ele = (/1d-4, 1d-4, 1d-4, 1d-6, 1d-6, 1d-6/)*100d0
   errStd_amrDrag = 0.1d0
   errStd_amrSRP = 0.1d0
   ! covariance matrix
   allocate(P(6, 6))
   P = 0d0
   P(1,1) = errStd_ele(1)**2
   P(2,2) = errStd_ele(2)**2
   P(3,3) = errStd_ele(3)**2
   P(4,4) = errStd_ele(4)**2
   P(5,5) = errStd_ele(5)**2
   P(6,6) = errStd_ele(6)**2
   !P(7,7) = errStd_amrDrag**2
   !P(8,8) = errStd_amrSRP**2
   
   ! Set introduction and output trajectory parameters
   introduction = 'EKF Orbit Determination for Station-Satellite Angle Observation'
   
   ! Initialize ODor using init procedure
   call ODor%init(ss, ss%ts%t0, 1, Groups, OD_METHOD_EKF, P = P, intro = trim(introduction))
   
   ! 7. Read observation data
   acc(1) = ss%access(1)
   call acc(1)%readData(fn_obsData)
   
   ! 8. EKF orbit filtering
   
   ! Get time series of observation data
   nEpochs = acc(1)%nData
   print*, "Number of observation epochs: ", nEpochs
   
   ! Create output directory
   outputDir = './data/output/'

   ! Open result file
   open(unit=200, file=trim(outputDir)//'ekf_filter_results.txt', status='replace')
   write(200, "(A)") "EKF Filter Results - Station-Satellite Angle Observation"
   write(200, "(A)") "Epoch(UTC MJD)    X(km)       Y(km)       Z(km)       VX(km/s)    VY(km/s)    VZ(km/s)"
   write(200, "(A)") "--------------------------------------------------------------------------------------"
     
   ! Open residual file
   open(unit=201, file=trim(outputDir)//'ekf_residuals.txt', status='replace')
   write(201, "(A)") "EKF Residuals - Station-Satellite Angle Observation"
   write(201, "(A)") "Epoch(UTC MJD)    Residual_Alpha(deg)    Residual_Delta(deg)"
   write(201, "(A)") "----------------------------------------------------------------"


   ! EKF filtering loop
   do i = 1, nEpochs
       
       currentEpoch = acc(1)%timeTable(i)
              
       ! Use EKF for orbit update
       call ODor%invoke(ss, currentEpoch, 1, acc(1)%npe, ODor%nVar, &
                       acc(1)%dataTable(i, 1:acc(1)%npe), residuals)
       
       ! Output current state
       write(*, "('num: ',I5.5,'/',I5.5,', utc = ',F12.6,', Res = ',F12.5)") &
           &    i, nEpochs, currentEpoch%val, norm2(residuals(:))
       
       ! Output to file
       write(200, '(I5,3(2x,F18.8),3(2x,F18.8))') nint(currentEpoch%day), currentEpoch%sec, &
           ss%sat(1)%states%elements(1:3)*length_unit, &
           ss%sat(1)%states%elements(4:6)*length_unit/time_unit
           
       ! Save residual data (convert to deg)
       write(201, '(I5,2x,F18.8,2x,F18.8,2x,F18.8)') nint(currentEpoch%day), currentEpoch%sec, &
           residuals(1)*rad2deg, residuals(2)*rad2deg
                       
   enddo
   
   close(200)
   close(201)

   deallocate(P)
   
   ! 11. Summary
   print*, " o_o EKF Orbit Determination Completed Successfully!"
   print*, ""
   
end program station_satellite_angle_ekf
