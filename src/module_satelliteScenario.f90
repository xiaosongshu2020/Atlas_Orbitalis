
! ============================================================================
! MODULE: module_satelliteScenario
!
! DESCRIPTION:
!   This module defines the satellite scenario object, which serves as a 
!   central container for managing satellite simulation scenarios. It 
!   integrates multiple components including time series, satellites, 
!   ground stations, and access calculations into a unified framework.
!
!   The satellite scenario object provides capabilities to:
!   - Initialize simulation time windows and parameters
!   - Add and manage satellite and ground station objects
!   - Define and calculate access links between objects
!   - Perform observation calculations for various object types
!   - Generate reports on access links and scenario status
!
!   The object uses Fortran's object-oriented features with type-bound
!   procedures to encapsulate related functionality and maintain clean
!   interfaces between components.
!
! COMPONENTS:
!   - Time series (ts): Defines the simulation time window and steps
!   - Satellite array (sat): Collection of satellite objects
!   - Station array (station): Collection of ground station objects
!   - Access array (access): Collection of access/link definitions
!   - Scenario metadata: Element type, introduction text
!
! KEY FEATURES:
!   - Flexible time series configuration
!   - Dynamic memory allocation for scenario objects
!   - Type-safe object classification and handling
!   - Comprehensive error checking and reporting
!   - Multiple observation modes (single access or batch)
!
! DEPENDENCIES:
!   - module_station: Ground station definitions
!   - module_satellite: Satellite definitions
!   - module_access: Access calculation routines
!   - module_attitude: Attitude dynamics
!   - module_atmosphere: Atmospheric models
!   - module_const: Physical and mathematical constants
!
! USAGE EXAMPLE:
!   type(type_satelliteScenario) :: scenario
!   call scenario%init(t0, te, 60.0d0, introduction="My simulation")
!   call scenario%addSatellite(mySat)
!   call scenario%addStation(myStation)
!   call scenario%observation()
!   call scenario%reportAccess()
!
! AUTHOR: Shawn Hu
!
! LAST MODIFIED: 2025-12-24 by Shawn Hu
!
! ============================================================================
module module_satelliteScenario
    use string
    use module_station
    use module_satellite
    use module_attitude
    use module_access
    use module_atmosphere
    use module_const
    
    implicit none
    
    
    ! Module constants
    integer, parameter :: MAX_DISPLAY_ACCESSES = 20
    integer, parameter :: UNKNOWN_OBJECT_TYPE = 0
    character(len=20), parameter :: STR_UNKNOWN = 'unknown'
    character(len=10), parameter :: STR_SATELLITE = 'satellite'
    character(len=10), parameter :: STR_STATION = 'station'
    
    ! Error codes
    integer, parameter :: ERROR_NONE = 0
    integer, parameter :: ERROR_INVALID_INPUT = 1
    integer, parameter :: ERROR_INDEX_OUT_OF_BOUNDS = 2
    integer, parameter :: ERROR_MEMORY_ALLOCATION = 3
    
    type :: type_satelliteScenario
        type(type_timeSeries) :: ts
        integer :: eleType = 3
        integer :: nSat = 0
        type(type_satellite), allocatable :: sat(:)
        integer :: nStation = 0
        type(type_station), allocatable :: station(:)
        integer :: nAccess = 0
        type(type_access), allocatable :: access(:)
        character(len=1000) :: introduction = ''
        
    contains
        procedure :: init => initializeSS
        procedure :: addSatellite
        procedure :: addStation
        procedure :: addAccess
        procedure :: observation
        procedure :: freeSS
        procedure :: reportAccess
    end type type_satelliteScenario
    
contains

! ============================================================================
! @brief Initialize a satellite scenario
! @param[inout] this   Satellite scenario object
! @param[in]    t0     Start time
! @param[in]    te     End time
! @param[in]    step   Time step in seconds
! @param[in]    eleType Element type (optional)
! @param[in]    introduction Description (optional)
! ============================================================================
subroutine initializeSS(this, t0, te, step, eleType, introduction)
    implicit none
    class(type_satelliteScenario) :: this
    type(type_t), intent(in) :: t0
    type(type_t), intent(in) :: te
    real(kind=8), intent(in) :: step
    integer, intent(in), optional :: eleType
    character(len=*), intent(in), optional :: introduction
    
    real(kind=8) :: dt_days, step_days
    integer :: n
    
    ! Check input validity
    if (step <= 0.0d0) then
        print *, "ERROR: Time step must be positive."
        return
    end if
    
    dt_days = timeMinus(te, t0)
    if (dt_days < 0.0d0) then
        print *, "ERROR: End time must be after start time."
        return
    end if
    
    step_days = step / 86400.0d0
    n = floor(dt_days / step_days) + 1
    
    ! Initialize time series
    call this%ts%init(t0, step, n)
    
    ! Set optional parameters
    if (present(eleType)) then
        this%eleType = eleType
    end if
    
    if (present(introduction)) then
        this%introduction = trim(introduction)
    end if
    
    ! Print scenario information
    print *, ""
    print *, "Created a satellite scenario"
    write(*, *) "Start time (UTC): ", t0%toDate()
    write(*, *) "End time (UTC):   ", te%toDate()
    write(*, *) "Time step (sec):  ", step
    write(*, *) "Number of epochs: ", n
    write(*, *) "Element type:     ", this%eleType
    
    if (len_trim(this%introduction) > 0) then
        write(*, *) "Description:      ", trim(this%introduction)
    end if
    
    print *, ""
    
end subroutine initializeSS

! ============================================================================
! @brief Add a satellite to the scenario
! @param[inout] this Satellite scenario object
! @param[in]    sat  Satellite object to add
! ============================================================================
subroutine addSatellite(this, sat)
    implicit none
    class(type_satelliteScenario) :: this
    type(type_satellite), intent(in) :: sat
    type(type_satellite), allocatable :: tempSatellites(:)
    integer :: iErr
    
    if (this%nSat == 0) then
        if (allocated(this%sat)) deallocate(this%sat)
        allocate(this%sat(1), stat=ierr)
        if (ierr /= 0) then
            print *, "ERROR: Memory allocation failed for satellites."
            return
        end if
        this%nSat = 1
        this%sat(1) = sat
    else
        allocate(tempSatellites(this%nSat + 1), stat=ierr)
        if (ierr /= 0) then
            print *, "ERROR: Memory allocation failed for satellite array expansion."
            return
        end if
        
        tempSatellites(1:this%nSat) = this%sat(1:this%nSat)
        tempSatellites(this%nSat + 1) = sat
        
        call move_alloc(tempSatellites, this%sat)
        this%nSat = this%nSat + 1
    end if
    
end subroutine addSatellite

! ============================================================================
! @brief Add a ground station to the scenario
! @param[inout] this    Satellite scenario object
! @param[in]    station Station object to add
! ============================================================================
subroutine addStation(this, station)
    implicit none
    class(type_satelliteScenario), intent(inout) :: this
    type(type_station), intent(in) :: station
    type(type_station), allocatable :: tempStations(:)
    integer :: ierr

    if (this%nStation == 0) then
        if (allocated(this%station)) deallocate(this%station)
        allocate(this%station(1), stat=ierr)
        if (ierr /= 0) then
            print *, "ERROR: Memory allocation failed for stations."
            return
        end if
        this%nStation = 1
        this%station(1) = station
    else
        allocate(tempStations(this%nStation + 1), stat=ierr)
        if (ierr /= 0) then
            print *, "ERROR: Memory allocation failed for station array expansion."
            return
        end if
        
        tempStations(1:this%nStation) = this%station(1:this%nStation)
        tempStations(this%nStation + 1) = station
        
        call move_alloc(tempStations, this%station)
        this%nStation = this%nStation + 1
    end if
    
end subroutine addStation

! ============================================================================
! @brief Add an access (link) to the scenario
! @param[inout] this   Satellite scenario object
! @param[in]    access Access object to add
! ============================================================================
subroutine addAccess(this, access)
    implicit none
    class(type_satelliteScenario), intent(inout) :: this
    type(type_access), intent(in) :: access
    type(type_access), allocatable :: tempAccesses(:)
    integer :: ierr
    
    if (this%nAccess == 0) then
        if (allocated(this%access)) deallocate(this%access)
        allocate(this%access(1), stat=ierr)
        if (ierr /= 0) then
            print *, "ERROR: Memory allocation failed for accesses."
            return
        end if
        this%nAccess = 1
        this%access(1) = access
    else
        allocate(tempAccesses(this%nAccess + 1), stat=ierr)
        if (ierr /= 0) then
            print *, "ERROR: Memory allocation failed for access array expansion."
            return
        end if
        
        tempAccesses(1:this%nAccess) = this%access(1:this%nAccess)
        tempAccesses(this%nAccess + 1) = access
        
        call move_alloc(tempAccesses, this%access)
        this%nAccess = this%nAccess + 1
    end if
    
end subroutine addAccess

! ============================================================================
! @brief Deallocate all memory in the scenario
! @param[inout] this Satellite scenario object
! ============================================================================
subroutine freeSS(this)
    implicit none
    class(type_satelliteScenario), intent(inout) :: this
    
    if (allocated(this%sat)) deallocate(this%sat)
    if (allocated(this%station)) deallocate(this%station)
    if (allocated(this%access)) deallocate(this%access)
    
    this%nSat = 0
    this%nStation = 0
    this%nAccess = 0
    
end subroutine freeSS

! ============================================================================
! @brief Print information about access links
! @param[in] this      Satellite scenario object
! @param[in] reportAll If true, report all accesses (default: false)
! @param[in] unitNumber Output unit number (default: 6 for screen)
! ============================================================================
subroutine reportAccess(this, reportAll, unitNumber)
    implicit none
    class(type_satelliteScenario) :: this
    logical, intent(in), optional :: reportAll
    integer, intent(in), optional :: unitNumber
    
    integer :: i, fromNum, toNum, maxLines, outputUnit
    character(len=20) :: fromName, toName, fromType, toType
    logical :: showAll
    
    ! Set output unit
    outputUnit = 6  ! Default to screen
    if (present(unitNumber)) outputUnit = unitNumber
    
    ! Set display limit
    showAll = .false.
    if (present(reportAll)) showAll = reportAll
    
    if (showAll .or. this%nAccess <= MAX_DISPLAY_ACCESSES) then
        maxLines = this%nAccess
    else
        maxLines = MAX_DISPLAY_ACCESSES
    end if
    
    ! Print header
    write(outputUnit, *) ""
    
    if (this%nAccess == 0) then
        write(outputUnit, *) "No access links found."
        return
    end if
    
    write(outputUnit, "('Number of access links: ', I5)") this%nAccess
    write(outputUnit, "(2x, A4, 2x, A20, 2x, A20, 2x, A15, 2x, A15)") &
          "No.", "From Name", "To Name", "From Type", "To Type"
    write(outputUnit, "(2x, 4('-'), 2x, 20('-'), 2x, 20('-'), 2x, 15('-'), 2x, 15('-'))")
    
    ! Print access information
    do i = 1, maxLines
        ! Determine from object information
        select case (this%access(i)%fromType)
        case (OBJECT_TYPE_SATELLITE)
            fromType = STR_SATELLITE
            fromNum = this%access(i)%fromNum
            if (fromNum > 0 .and. fromNum <= this%nSat) then
                fromName = this%sat(fromNum)%name
            else
                fromName = STR_UNKNOWN
            end if
            
        case (OBJECT_TYPE_STATION)
            fromType = STR_STATION
            fromNum = this%access(i)%fromNum
            if (fromNum > 0 .and. fromNum <= this%nStation) then
                fromName = this%station(fromNum)%name
            else
                fromName = STR_UNKNOWN
            end if
            
        case default
            fromType = STR_UNKNOWN
            fromName = STR_UNKNOWN
        end select
        
        ! Determine to object information
        select case (this%access(i)%toType)
        case (OBJECT_TYPE_SATELLITE)
            toType = STR_SATELLITE
            toNum = this%access(i)%toNum
            if (toNum > 0 .and. toNum <= this%nSat) then
                toName = this%sat(toNum)%name
            else
                toName = STR_UNKNOWN
            end if
            
        case (OBJECT_TYPE_STATION)
            toType = STR_STATION
            toNum = this%access(i)%toNum
            if (toNum > 0 .and. toNum <= this%nStation) then
                toName = this%station(toNum)%name
            else
                toName = STR_UNKNOWN
            end if
            
        case default
            toType = STR_UNKNOWN
            toName = STR_UNKNOWN
        end select
        
        write(outputUnit, "(I6, 2x, A20, 2x, A20, 2x, A15, 2x, A15)") &
              i, trim(fromName), trim(toName), trim(fromType), trim(toType)
    end do
    
    ! Indicate if not all accesses are shown
    if (maxLines < this%nAccess) then
        write(outputUnit, *) "... (showing ", maxLines, " of ", this%nAccess, " accesses)"
    end if
    
    write(outputUnit, *) ""
    
end subroutine reportAccess

! ============================================================================
! @brief Perform observation calculations
! @param[in] this   Satellite scenario object
! @param[in] accIdx Access index (optional, if not present do all)
! ============================================================================
subroutine observation(this, accIdx)
    implicit none
    class(type_satelliteScenario) :: this
    integer, intent(in), optional :: accIdx
    integer :: i
    
    ! Observation for a specific access
    if (present(accIdx)) then
        if (accIdx < 1 .or. accIdx > this%nAccess) then
            print *, "ERROR: Access index out of bounds."
            return
        end if
        
        call observationAccess(this, this%access(accIdx), this%ts, .false.)
        return
    end if
    
    ! Observation for all accesses
    do i = 1, this%nAccess
        call observationAccess(this, this%access(i), this%ts, .false.)
    end do
    
end subroutine observation

! ============================================================================
! @brief Calculate observation data for a specific access
! @param[in]    ss     Satellite scenario object
! @param[inout] access Access object
! @param[in]    ts     Time series
! @param[in]    force  Force recalculation flag
! ============================================================================
subroutine observationAccess(ss, access, ts, force)
    implicit none
    type(type_satelliteScenario) :: ss
    type(type_access), intent(inout) :: access
    type(type_timeSeries), intent(in) :: ts
    logical, intent(in) :: force
    
    integer :: fromType, toType, obsType, fromNum, toNum

    fromType = access%fromType
    toType = access%toType
    obsType = access%obsType
    fromNum = access%fromNum
    toNum = access%toNum

    ! Keep the names of the two devices in the access, so that the data file
    ! written by saveData identifies the link it came from
    call setAccessDeviceNames(ss, access)

    ! Route to appropriate calculation routine based on object types
    select case (fromType)
    case (OBJECT_TYPE_SATELLITE)
        select case (toType)
        case (OBJECT_TYPE_STATION)
            call calObsData1(ts, ss%sat(fromNum), ss%station(toNum), access, force)
            
        case (OBJECT_TYPE_SATELLITE)
            call calObsData2(ts, ss%sat(fromNum), ss%sat(toNum), access, force)
            
        case (UNKNOWN_OBJECT_TYPE)
            call calObsData3(ts, ss%sat(fromNum), access, force)
            
        case default
            print *, "ERROR: Unknown toType in observationAccess."
        end select
        
    case default
        print *, "ERROR: Unknown fromType in observationAccess."
    end select

end subroutine observationAccess

! ============================================================================
! @brief Copy the names of the objects an access refers to into the access
! @param[in]    ss     Satellite scenario object
! @param[inout] access Access object
! ============================================================================
subroutine setAccessDeviceNames(ss, access)
    implicit none
    type(type_satelliteScenario), intent(in) :: ss
    type(type_access), intent(inout) :: access
    integer :: n

    select case (access%fromType)
    case (OBJECT_TYPE_SATELLITE)
        n = access%fromNum
        if (n >= 1 .and. n <= ss%nSat) access%fromName = ss%sat(n)%name
    case (OBJECT_TYPE_STATION)
        n = access%fromNum
        if (n >= 1 .and. n <= ss%nStation) access%fromName = ss%station(n)%name
    end select

    select case (access%toType)
    case (OBJECT_TYPE_SATELLITE)
        n = access%toNum
        if (n >= 1 .and. n <= ss%nSat) access%toName = ss%sat(n)%name
    case (OBJECT_TYPE_STATION)
        n = access%toNum
        if (n >= 1 .and. n <= ss%nStation) access%toName = ss%station(n)%name
    end select

end subroutine setAccessDeviceNames

end module module_satelliteScenario