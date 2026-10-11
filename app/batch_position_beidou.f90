!================================================================================
! PROGRAM: batch_position_beidou
!
! Description:
!   Batch orbit determination using SP3 position observations for BeiDou C06.
!   This program reads SP3 file, extracts C06 satellite position data,
!   converts from ECF to ECI frame, estimates initial velocity, and performs
!   batch orbit determination.
!
!   Steps:
!   1. Read SP3 file and extract C06 position data.
!   2. Convert position from ECF to GCRS (ECI) frame.
!   3. Estimate initial velocity using finite difference.
!   4. Initialize satellite, satelliteScenario, access, and ODor objects.
!   5. Execute orbit determination.
!
! Author: Shawn Hu
! Date: 2026-05-15
!================================================================================

program batch_position_beidou

    use module_Atlas

    implicit none

    !--------------------------------------------------------------------------
    ! Variables
    !--------------------------------------------------------------------------
    
    ! SP3 data
    type(type_sp3) :: sp3
    double precision :: pos_ECF(3)
    
    ! Time
    type(type_t) :: t0, te
    type(type_t), allocatable :: epoch(:)
    double precision :: utc_val
    
    ! Position and velocity arrays (ECF and ECI)
    double precision, allocatable :: pos_ECF_all(:, :)
    double precision, allocatable :: pos_ECI_all(:, :)
    double precision, allocatable :: vel_ECI_all(:, :)
    
    ! Initial state
    double precision :: x0_cts(6)   ! Cartesian [x, y, z, vx, vy, vz] (km, km/s)
    double precision :: x0_nsg1(6)  ! NSG1 elements
    
    ! Satellite scenario
    type(type_satelliteScenario) :: ss
    
    ! Satellite
    type(type_satellite) :: sat
    
    ! Access
    type(type_access) :: access
    
    ! ODor
    type(type_ODor) :: ODor
    type(type_varGroup) :: Groups(1)
    
    ! Prediction time series
    type(type_timeSeries) :: pred_ts
    
    ! Loop variables
    integer :: i, nEpoch
    integer :: year, month, day, hour, min
    double precision :: sec
    
    ! Satellite PRN
    character(len=3), parameter :: SAT_PRN = "C20"
    
    ! File I/O
    character(len=*), parameter :: SP3_FILE      = "data/input/isr10620.sp3"
    character(len=*), parameter :: OBS_FILE      = "data/output/obs_" // SAT_PRN // "_ECI.txt"
    character(len=*), parameter :: EPH_FILE      = "data/output/preciseEph_" // SAT_PRN // ".txt"
    character(len=*), parameter :: PRED_EPH_FILE = "data/output/predictedEph_" // SAT_PRN // ".txt"
    character(len=*), parameter :: REF_SP3_FILE  = "data/input/isr10621.sp3"
    character(len=*), parameter :: ERR_FILE      = "data/output/positionError_" // SAT_PRN // ".txt"
    
    !--------------------------------------------------------------------------
    ! Step 1: Read SP3 file and extract C06 position data
    !--------------------------------------------------------------------------
    
    print *, "============================================"
    print '("  Batch OD for BeiDou ", A, " using SP3 data")', SAT_PRN
    print *, "============================================"
    print *, ""
    
    print *, ">>> Step 1: Loading SP3 file ..."
    call sp3%load(SP3_FILE)
    if (errCode > 0) then
        call PrintErrorInfo
        stop
    endif
    print *, "    Load successful!"
    print '("    Number of epochs: ", I0)', sp3%header%nEpoch
    print '("    Time interval:    ", F14.8, " seconds")', sp3%header%dt
    print *, ""
    
    ! Get number of epochs
    nEpoch = sp3%header%nEpoch
    
    ! Allocate arrays
    allocate(pos_ECF_all(nEpoch, 3))
    allocate(pos_ECI_all(nEpoch, 3))
    allocate(vel_ECI_all(nEpoch, 3))
    allocate(epoch(nEpoch))
    
    ! Extract satellite position data for all epochs
    print '(">>> Step 2: Extracting ", A, " position data ...")', SAT_PRN
    do i = 1, nEpoch
        ! Get time from SP3 data
        year = sp3%data(i, 1)%t%yr
        month = sp3%data(i, 1)%t%mo
        day = sp3%data(i, 1)%t%day
        hour = sp3%data(i, 1)%t%hr
        min = sp3%data(i, 1)%t%mi
        sec = sp3%data(i, 1)%t%sec
        
        ! Create time object
        call epoch(i)%init(dble(year), dble(month), dble(day), dble(hour), dble(min), sec)
        
        ! Get satellite position in ECF (km)
        pos_ECF = sp3%getPosition(year, month, day, hour, min, sec, SAT_PRN)
        if (errCode > 0) then
            call PrintErrorInfo
            stop
        endif
        pos_ECF_all(i, :) = pos_ECF
    enddo
    print '("    Extracted ", I0, " epochs of ", A, " data")', nEpoch, SAT_PRN
    print *, ""
    
    !--------------------------------------------------------------------------
    ! Step 2: Convert ECF positions to ECI (GCRS) using IAU2006A
    !--------------------------------------------------------------------------
    
    print *, ">>> Step 3: Converting ECF to ECI (IAU2006A) ..."
    
    do i = 1, nEpoch
        ! Get UTC value (MJD)
        utc_val = epoch(i)%val
        
        ! Convert ECF to ECI (position only, IAU2006A)
        call ECF2ECI_IAU2006A_position(utc_val, pos_ECF_all(i, :), pos_ECI_all(i, :))
    enddo
    
    print *, "    Conversion complete!"
    print *, ""
    
    !--------------------------------------------------------------------------
    ! Step 3: Estimate initial velocity using first two epochs
    !--------------------------------------------------------------------------
    
    print *, ">>> Step 4: Estimating initial velocity ..."
    
    ! Use the difference between epoch 2 and epoch 1 divided by time interval
    vel_ECI_all(1, :) = (pos_ECI_all(2, :) - pos_ECI_all(1, :)) / sp3%header%dt
    
    print '("    dt = ", F14.8, " seconds")', sp3%header%dt
    print '("    Initial velocity (km/s): ", 3F14.8)', vel_ECI_all(1, :)
    print *, ""
    
    !--------------------------------------------------------------------------
    ! Step 4: Write observation data to temporary file
    !--------------------------------------------------------------------------
    
    print *, ">>> Step 5: Writing observation data to file ..."
    
    open(unit=201, file=OBS_FILE)
    write(201, '(A)') SAT_PRN // " position observations in ECI frame (km)"
    write(201, '(A)') "Format: Day Sec X/km Y/km Z/km"
    do i = 1, nEpoch
        write(201, '(I5, F10.3, 3F18.6)') nint(epoch(i)%day), epoch(i)%sec, &
            pos_ECI_all(i, 1), pos_ECI_all(i, 2), pos_ECI_all(i, 3)
    enddo
    close(201)
    print '("    Written to: ", A)', trim(OBS_FILE)
    print *, ""
    
    !--------------------------------------------------------------------------
    ! Step 5: Initialize satellite scenario
    !--------------------------------------------------------------------------
    
    print *, ">>> Step 6: Initializing satellite scenario ..."
    
    ! Use the first epoch time as t0, last epoch as te
    t0 = epoch(1)
    te = epoch(nEpoch)
    
    ! Initialize satellite scenario
    call ss%init(t0, te, sp3%header%dt, &
                 introduction="Batch orbit determination for BeiDou " // SAT_PRN // " using SP3 position data")
    
    !--------------------------------------------------------------------------
    ! Step 6: Initialize satellite object
    !--------------------------------------------------------------------------
    
    print *, ">>> Step 7: Initializing satellite object ..."
    
    ! Initial state in ECI (Cartesian, km and km/s)
    x0_cts(1:3) = pos_ECI_all(1, :)
    x0_cts(4:6) = vel_ECI_all(1, :)
    call cts_SI2normal(x0_cts, x0_cts)
    
    ! Convert Cartesian to NSG1 elements (for satellite initialization)
    call cts2nsg1(x0_cts, x0_nsg1)
    
    ! Initialize satellite
    call sat%init(t0, x0_nsg1, &
                  number = 1, &
                  name = SAT_PRN, &
                  eleType = STATE_TYPE_NONSINGULAR_1, &
                  propagator = 'numerical', &
                  Cd = 2.2d0, &
                  amrDrag = 0.02d0, &
                  F107 = 130d0, &
                  Ap = 15d0, &
                  reflectivity = 0.3d0, &
                  amrSRP = 0.02d0, &
                  degree = 20, order = 20)
    sat%opMod%key_moonGravity = 1   ! Moon third-body
    sat%opMod%key_moonTide    = 1   ! Moon solid tide
    sat%opMod%key_sunGravity  = 1   ! Sun third-body
    sat%opMod%key_sunTide     = 1   ! Sun solid tide
    sat%opMod%key_srp         = 1   ! solar radiation pressure
    sat%opMod%key_postNewton  = 1   ! post-Newtonian
    sat%opMod%SRP_Model = 2
    sat%opMod%ECOM = 0d0
    sat%opMod%key_ERP = 1
    
    ! Add satellite to scenario
    call ss%addSatellite(sat)
    print '("    Satellite ", A, " added to scenario.")', trim(SAT_PRN)
    print *, ""
    
    !--------------------------------------------------------------------------
    ! Step 7: Initialize access object
    !--------------------------------------------------------------------------
    
    print *, ">>> Step 8: Initializing access object ..."
    
    ! Initialize access with OBS_TYPE_POSITION (15)
    ! fromType = OBJECT_TYPE_SATELLITE, toType = 0 (no receiver station)
    call access%init(number = 1, &
                     name = 'access_' // SAT_PRN, &
                     fromType = OBJECT_TYPE_SATELLITE, &
                     fromNum = 1, &
                     toType = 0, &
                     toNum = 0, &
                     obsType = OBS_TYPE_POSITION)
    
    ! Set observation error standard deviation
    ! Position accuracy: ~1m = 0.001 km, normalized: 0.001 / length_unit
    access%errStd(1:3) = 1.0d-3 / LENGTH_UNIT
    
    ! Read observation data from file (skip 2 header lines)
    call access%readData(OBS_FILE)
    
    ! Add access to scenario
    call ss%addAccess(access)
    
    ! Report access relationships
    call ss%reportAccess(.true.)
    
    !--------------------------------------------------------------------------
    ! Step 8: Initialize orbit determinator (ODor)
    !--------------------------------------------------------------------------
    
    print *, ">>> Step 9: Initializing orbit determinator ..."
    
    !! Set up variable group parameters
    !Groups(1)%nVar = 6                     ! 6 state variables (position and velocity)
    !Groups(1)%label = 1                    ! states + ECOM 9
    !Groups(1)%origin = 1                   ! Belongs to satellite 1
    !! Allocate and set indices
    !allocate(Groups(1)%idx(6))
    !Groups(1)%idx = (/ (i, i = 1, 6) /)   ! Indices for the 6 state variables
    
    ! Set up variable group parameters
    Groups(1)%nVar = 15                     ! 6 state variables (position and velocity)
    Groups(1)%label = 5                    ! states + ECOM 9
    Groups(1)%origin = 1                   ! Belongs to satellite 1
    ! Allocate and set indices
    allocate(Groups(1)%idx(15))
    Groups(1)%idx = (/ (i, i = 1, 15) /)   ! Indices for the 6 state variables
    
    !! Set up variable group parameters
    !Groups(1)%nVar = 11                     ! 6 state variables (position and velocity)
    !Groups(1)%label = 5                    ! states + ECOM 9
    !Groups(1)%origin = 1                   ! Belongs to satellite 1
    !! Allocate and set indices
    !allocate(Groups(1)%idx(11))
    !Groups(1)%idx = (/ (i, i = 1, 11) /)   ! Indices for the 6 state variables
    
    ! Initialize ODor with batch method
    call ODor%init(ss, ss%ts%t0, 1, Groups, OD_METHOD_BATCH, &
                   intro = "Batch orbit determination for BeiDou " // SAT_PRN // " using SP3 position data")
    
    print *, "    ODor initialized."
    print *, ""
    
    !--------------------------------------------------------------------------
    ! Step 9: Perform batch orbit determination
    !--------------------------------------------------------------------------
    
    print *, ">>> Step 10: Performing batch orbit determination ..."
    print *, ""
    
    call ODor%solve(ss)
    
    print *, ""
    print *, "    Batch orbit determination complete!"
    print *, ""
    
    !--------------------------------------------------------------------------
    ! Step 10: Output results
    !--------------------------------------------------------------------------
    
    print *, ">>> Step 11: Saving results ..."
    
    call ODor%saveEph(ss, EPH_FILE)
    print '("    Precise ephemeris saved to: ", A)', trim(EPH_FILE)
    print *, ""
    
    !--------------------------------------------------------------------------
    ! Step 11: Generate prediction for the next day
    !--------------------------------------------------------------------------
    
    print *, ">>> Step 12: Generating prediction for next day ..."
    
    ! Create prediction time series
    ! Prediction period: from te+900s (2026-05-11 00:00:00) to te + 1 day
    ! Step size: same as SP3 data (900 seconds), 96 points for 1 day
    call pred_ts%init(timeAdd(te, 0d0, 900d0), sp3%header%dt, 96)
    
    ! Generate trajectory for the satellite at prediction times
    call ss%sat(1)%generateTrajectory(pred_ts%time)
    
    ! Save predicted ephemeris in Cartesian format
    call ss%sat(1)%savetrajectory(STATE_TYPE_CARTESIAN, PRED_EPH_FILE, .false.)
    print '("    Predicted ephemeris saved to: ", A)', trim(PRED_EPH_FILE)
    print *, ""
    
    !--------------------------------------------------------------------------
    ! Step 12: Compare predicted ephemeris with SP3 reference data
    !--------------------------------------------------------------------------
    
    print *, ">>> Step 13: Computing position errors against SP3 reference ..."
    
    call calcPositionError(PRED_EPH_FILE, &
                           REF_SP3_FILE, &
                           SAT_PRN, &
                           ERR_FILE)
    
    print *, ""
    
    !--------------------------------------------------------------------------
    ! Step 13: Clean up
    !--------------------------------------------------------------------------
    
    print *, ">>> Cleaning up ..."
    
    call pred_ts%free()
    call ss%freeSS
    deallocate(Groups(1)%idx)
    call sp3%freeData()
    
    if (allocated(pos_ECF_all)) deallocate(pos_ECF_all)
    if (allocated(pos_ECI_all)) deallocate(pos_ECI_all)
    if (allocated(vel_ECI_all)) deallocate(vel_ECI_all)
    if (allocated(epoch)) deallocate(epoch)
    
    print *, ""
    print *, "============================================"
    print *, "  Program completed successfully!"
    print *, "============================================"
    
contains

!==============================================================================
! SUBROUTINE: calcPositionError
!
! Description:
!   Compare predicted ephemeris with SP3 reference data and compute
!   position differences in meters.
!
!   The predicted ephemeris file is in Cartesian format (ECI frame) with
!   30 header lines. The SP3 file contains positions in ECF frame.
!   SP3 positions are converted to ECI before comparison.
!
! Inputs:
!   pred_file   - Predicted ephemeris file (Cartesian format)
!   sp3_file    - Reference SP3 file
!   sat_id      - Satellite ID (e.g., "C06")
!   output_file - Output file for error results
!==============================================================================
subroutine calcPositionError(pred_file, sp3_file, sat_id, output_file)
    character(len=*), intent(in) :: pred_file, sp3_file, sat_id, output_file
    
    ! Local variables
    type(type_sp3) :: ref_sp3
    double precision :: pos_ECF(3), pos_ECI(3)
    double precision :: pred_pos(3), ref_pos(3)
    double precision :: diff(3), diff_norm
    double precision :: pred_day, pred_sec, pred_x, pred_y, pred_z
    double precision :: utc_val, jd_val
    integer :: i, nRef, iostat
    integer :: year, month, day, hour, min
    double precision :: sec, day_frac
    integer :: file_unit, out_unit
    integer :: j_status
    
    print *, "    Loading reference SP3 file: ", trim(sp3_file)
    
    ! Load reference SP3 file
    call ref_sp3%load(sp3_file)
    if (errCode > 0) then
        call PrintErrorInfo
        return
    endif
    
    nRef = ref_sp3%header%nEpoch
    print '("    Reference SP3 has ", I0, " epochs")', nRef
    
    ! Open predicted ephemeris file (skip 30 header lines)
    open(newunit=file_unit, file=pred_file, status='old')
    do i = 1, 30
        read(file_unit, *)
    enddo
    
    ! Open output file
    open(newunit=out_unit, file=output_file, status='replace')
    write(out_unit, '(A)') "# Position error between predicted ephemeris and SP3 reference"
    write(out_unit, '(A)') "# Reference SP3: " // trim(sp3_file)
    write(out_unit, '(A)') "# Satellite: " // trim(sat_id)
    write(out_unit, '(A)') "# Units: meters"
    write(out_unit, '(A)') "# Format: Day, Sec, dx(m), dy(m), dz(m), |dr|(m)"
    write(out_unit, '(A)') "# ============================================================"
    
    ! Read predicted ephemeris and compare with SP3 reference
    iostat = 0
    do while (iostat == 0)
        ! Read one line from predicted ephemeris
        read(file_unit, *, iostat=iostat) pred_day, pred_sec, pred_x, pred_y, pred_z
        if (iostat /= 0) exit
        
        ! Convert MJD (pred_day, pred_sec) to calendar date using JD2date
        ! JD = MJD + 2400000.5, where MJD = pred_day + pred_sec/86400
        jd_val = pred_day + 2400000.5d0 + pred_sec / 86400.0d0
        call JD2date(jd_val, 0d0, year, month, day, day_frac, j_status)
        if (j_status /= 0) then
            print *, "Warning: JD2date conversion failed at day=", pred_day, " sec=", pred_sec
            cycle
        endif
        ! Convert day fraction to hour, minute, second
        hour = floor(day_frac * 24d0)
        min = floor((day_frac - dble(hour)/24d0) * 1440d0)
        sec = (day_frac - dble(hour)/24d0 - dble(min)/1440d0) * 86400d0
        
        ! Get reference position from SP3 (ECF)
        pos_ECF = ref_sp3%getPosition(year, month, day, hour, min, sec, sat_id)
        if (errCode > 0) then
            call PrintErrorInfo
            cycle
        endif
        
        ! Convert reference position from ECF to ECI
        ! MJD value = pred_day + pred_sec / 86400
        utc_val = pred_day + pred_sec / 86400.0d0
        call ECF2ECI_IAU2006A_position(utc_val, pos_ECF, pos_ECI)
        
        ! Compute position difference (km -> m)
        pred_pos = (/pred_x, pred_y, pred_z/)
        ref_pos = pos_ECI
        diff = (pred_pos - ref_pos) * 1000.0d0  ! Convert km to m
        diff_norm = sqrt(diff(1)**2 + diff(2)**2 + diff(3)**2)
        
        ! Write to output file
        write(out_unit, '(F8.1, 3x, F10.3, 4(3x, F12.4))') &
            pred_day, pred_sec, diff(1), diff(2), diff(3), diff_norm
    enddo
    
    close(file_unit)
    close(out_unit)
    
    call ref_sp3%freeData()
    
    print '("    Position errors saved to: ", A)', trim(output_file)
    
end subroutine calcPositionError

end program batch_position_beidou
