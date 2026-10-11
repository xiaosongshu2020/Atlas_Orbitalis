!===============================================================================
! module_cartesian.f90
!
! Purpose: Define cartesian coordinate object with coordinate transformation methods
! Author: Shawn Hu
!
! Contents:
! - type_cartesian object definition with position, velocity, reference frame and time
! - Coordinate transformation methods between ICRS, ITRS and horizon frames
! - Integration with existing module_space transformation routines
!
! Reference Frames:
! - ICRS: International Celestial Reference System (Earth-centered inertial)
! - ITRS: International Terrestrial Reference System (Earth-centered Earth-fixed)
! - HORIZON: Local horizon frame (topocentric)
!
! Time Handling:
! - Includes type_t time object for temporal reference
! - Coordinate transformations use internal time attribute
!
! Modern Fortran features used:
! - IMPLICIT NONE for type safety
! - Type-bound procedures for object-oriented design
! - Descriptive naming conventions
!===============================================================================

module module_cartesian
    use module_const
    use module_time
    use module_space
    implicit none
    
    !===========================================================================
    ! COORDINATE FRAME TYPE CONSTANTS
    !===========================================================================
    
    integer, parameter :: COORD_ICRS    = 1  ! International Celestial Reference System
    integer, parameter :: COORD_ITRS    = 2  ! International Terrestrial Reference System  
    integer, parameter :: COORD_HORIZON = 3  ! Local horizon/topocentric frame
    
    !===========================================================================
    ! TYPE_CARTESIAN OBJECT DEFINITION
    !===========================================================================
    
    type type_cartesian
        ! Position vector (3 components) in km
        double precision :: position(3) = 0d0
        
        ! Velocity vector (3 components) in km/s
        double precision :: velocity(3) = 0d0
        
        ! Reference frame identifier
        integer :: reference = COORD_ICRS
        
        ! Time attribute (type_t object)
        type(type_t) :: t
        
    contains
        ! Initialization procedure
        procedure :: init => initializeCartesian
        
        ! Coordinate transformation procedures
        procedure :: ICRS2ITRS
        procedure :: ITRS2ICRS
        procedure :: ITRS2Horizon
        procedure :: Horizon2ITRS
        
        procedure :: getReferenceName
    end type type_cartesian
    
contains
    
    !===========================================================================
    ! INITIALIZATION SUBROUTINE
    !===========================================================================
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize a cartesian object with position, velocity, reference frame and time
    !
    ! Parameters:
    !   this    : type_cartesian object to initialize
    !   pos     : position vector (3 components) in km
    !   vel     : velocity vector (3 components) in km/s
    !   ref     : reference frame identifier (COORD_ICRS, COORD_ITRS, COORD_HORIZON)
    !   time    : time object (type_t) associated with the state
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine initializeCartesian(this, pos, vel, ref, time)
        implicit none
        class(type_cartesian) :: this
        double precision, intent(in) :: pos(3)
        double precision, intent(in) :: vel(3)
        integer, intent(in) :: ref
        type(type_t), intent(in) :: time
        
        this%position = pos
        this%velocity = vel
        this%reference = ref
        this%t = time
        
    end subroutine initializeCartesian
    
    !===========================================================================
    ! COORDINATE TRANSFORMATION METHODS
    !===========================================================================
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Transform from ICRS to ITRS frame
    !
    ! Parameters:
    !   this : type_cartesian object (input in ICRS, output in ITRS)
    !
    ! Notes:
    ! - Uses IAU1980 transformation from module_space
    ! - Checks if velocity is zero for simplified calculation
    ! - Uses internal time attribute (this%t) for transformation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ICRS2ITRS(this)
        implicit none
        class(type_cartesian) :: this
        
        double precision :: utc_time
        double precision :: rITRS(3), vITRS(3)
        
        ! Check current reference frame
        if (this%reference /= COORD_ICRS) then
            write(*, *) "Warning: ICRS2ITRS called on non-ICRS object"
            return
        endif
        
        ! Convert internal time to UTC MJD for transformation
        utc_time = this%t%val
        
        ! Check if velocity is zero for simplified calculation
        if (norm2(this%velocity) < 1d-12) then
            ! Position only transformation
            call ECI2ECF_IAU1980_position(utc_time, this%position, rITRS)
            vITRS = 0d0
        else
            ! Full position and velocity transformation
            call ECI2ECF_IAU1980_positionVelocity(utc_time, this%position, this%velocity, rITRS, vITRS)
        endif
        
        ! Update object state
        this%position = rITRS
        this%velocity = vITRS
        this%reference = COORD_ITRS
        
    end subroutine ICRS2ITRS
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Transform from ITRS to ICRS frame
    !
    ! Parameters:
    !   this : type_cartesian object (input in ITRS, output in ICRS)
    !
    ! Notes:
    ! - Uses IAU1980 transformation from module_space
    ! - Checks if velocity is zero for simplified calculation
    ! - Uses internal time attribute (this%t) for transformation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ITRS2ICRS(this)
        implicit none
        class(type_cartesian) :: this
        
        double precision :: utc_time
        double precision :: rICRS(3), vICRS(3)
        
        ! Check current reference frame
        if (this%reference /= COORD_ITRS) then
            write(*, *) "Warning: ITRS2ICRS called on non-ITRS object"
            return
        endif
        
        ! Convert internal time to UTC MJD for transformation
        utc_time = this%t%val
        
        ! Check if velocity is zero for simplified calculation
        if (norm2(this%velocity) < 1d-12) then
            ! Position only transformation
            call ECF2ECI_IAU1980_position(utc_time, this%position, rICRS)
            vICRS = 0d0
        else
            ! Full position and velocity transformation
            call ECF2ECI_IAU1980_positionVelocity(utc_time, this%position, this%velocity, rICRS, vICRS)
        endif
        
        ! Update object state
        this%position = rICRS
        this%velocity = vICRS
        this%reference = COORD_ICRS
        
    end subroutine ITRS2ICRS
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Transform from ITRS to local horizon frame
    !
    ! Parameters:
    !   this       : type_cartesian object (input in ITRS, output in HORIZON)
    !   longitude  : observer longitude in radians (East positive)
    !   latitude   : observer latitude in radians (North positive)
    !
    ! Notes:
    ! - Uses ECF2horizon_position from module_space
    ! - Velocity transformation not implemented in module_space, set to zero
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ITRS2Horizon(this, longitude, latitude)
        implicit none
        class(type_cartesian) :: this
        double precision, intent(in) :: longitude, latitude
        
        double precision :: rHorizon(3)
        
        ! Check current reference frame
        if (this%reference /= COORD_ITRS) then
            write(*, *) "Warning: ITRS2Horizon called on non-ITRS object"
            return
        endif
        
        ! Position transformation
        call ECF2horizon_position(longitude, latitude, this%position, rHorizon)
        
        ! Update object state
        this%position = rHorizon
        this%velocity = 0d0  ! Velocity in horizon frame not available
        this%reference = COORD_HORIZON
        
    end subroutine ITRS2Horizon
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Transform from local horizon frame to ITRS
    !
    ! Parameters:
    !   this       : type_cartesian object (input in HORIZON, output in ITRS)
    !   longitude  : observer longitude in radians (East positive)
    !   latitude   : observer latitude in radians (North positive)
    !
    ! Notes:
    ! - Uses horizon2ECF_position from module_space
    ! - Velocity transformation not implemented in module_space, set to zero
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Horizon2ITRS(this, longitude, latitude)
        implicit none
        class(type_cartesian) :: this
        double precision, intent(in) :: longitude, latitude
        
        double precision :: rITRS(3)
        
        ! Check current reference frame
        if (this%reference /= COORD_HORIZON) then
            write(*, *) "Warning: Horizon2ITRS called on non-HORIZON object"
            return
        endif
        
        ! Position transformation
        call horizon2ECF_position(longitude, latitude, this%position, rITRS)
        
        ! Update object state
        this%position = rITRS
        this%velocity = 0d0  ! Velocity in ITRS frame not available from horizon
        this%reference = COORD_ITRS
        
    end subroutine Horizon2ITRS
    
    !===========================================================================
    ! UTILITY FUNCTIONS (Optional - for completeness)
    !===========================================================================
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get reference frame name as string
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function getReferenceName(this) result(name)
        implicit none
        class(type_cartesian) :: this
        character(len=20) :: name
        
        select case(this%reference)
        case(COORD_ICRS)
            name = "ICRS"
        case(COORD_ITRS)
            name = "ITRS"
        case(COORD_HORIZON)
            name = "HORIZON"
        case default
            name = "UNKNOWN"
        end select
        
    end function getReferenceName
    
end module module_cartesian
