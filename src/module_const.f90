!===============================================================================
! module_const.f90
!
! Purpose: Define physical constants and parameters for the Atlas project
! Author: Shawn Hu
!
! Contents:
! - Fundamental physical constants
! - Earth and Moon parameters
! - Unit conversion factors
! - Observation and orbit determination type identifiers
!
! Modern Fortran features used:
! - IMPLICIT NONE for type safety
! - Descriptive naming conventions
! - English comments for international collaboration
!===============================================================================

module module_const
    implicit none
    
    !===========================================================================
    ! UNIT DEFINITIONS
    !===========================================================================
    
    ! units
    double precision, parameter :: LENGTH_UNIT       = 6378.1363d0          ! km
    double precision, parameter :: TIME_UNIT         = 806.81099130673260d0 ! s
    double precision, parameter :: MASS_UNIT         = 5.9742272407074334d24 ! kg
    double precision, parameter :: ACC_UNIT          = 9.79828762253515d0   ! m/s^2
    double precision, parameter :: FORCE_UNIT        = 5.853719682683597d25 ! N

    ! units -> normalized units
    double precision, parameter :: u_km = 1d0/LENGTH_UNIT
    double precision, parameter :: u_meter  = 1d0/(LENGTH_UNIT*1000.0d0)
    double precision, parameter :: u_feet = 0.3048d0 * u_meter
    double precision, parameter :: u_sec = 1d0/TIME_UNIT
    double precision, parameter :: u_min = 60d0/TIME_UNIT
    double precision, parameter :: u_hour = 3600d0/TIME_UNIT
    double precision, parameter :: u_day = 86400d0/TIME_UNIT
    double precision, parameter :: u_kg = 1d0/MASS_UNIT
    double precision, parameter :: u_m_s2 = 1d0/(LENGTH_UNIT/(TIME_UNIT**2))
    double precision, parameter :: u_Newton = 1d0/FORCE_UNIT
    
    !===========================================================================
    ! EARTH PARAMETERS
    !===========================================================================
    
    ! GGM03C Earth model
    double precision, parameter :: EARTH_GM                = 3.986004415d14      ! m³/s²
    double precision, parameter :: EARTH_RADIUS_EQUATORIAL = 6378.1363d0        ! km
    double precision, parameter :: EARTH_FLATTENING        = 1.0d0/298.256415d0
    
    ! WGS84 Earth model (commented out, alternative values)
    ! double precision, parameter :: EARTH_RADIUS_WGS84    = 6378.137d0         ! km
    ! double precision, parameter :: EARTH_FLATTENING_WGS84 = 1.0d0/298.257223563d0
    
    !===========================================================================
    ! EARTH-MOON SYSTEM (Normalized Units)
    !===========================================================================
    
    double precision, parameter :: EARTH_MASS            = 5.9742272407074334d24 ! kg
    double precision, parameter :: MOON_MASS             = 7.34832238736047d22   ! kg
    double precision, parameter :: MOON_EARTH_MASS_RATIO = 0.0121505856313668d0
    double precision, parameter :: MOON_EARTH_MASS       = 6.04771046458104d24   ! kg
    double precision, parameter :: MOON_EARTH_LENGTH = 384401000.0d0        ! m
    double precision, parameter :: MOON_EARTH_TIME   = 375163.609176656d0   ! s (~4.34 days)
    
    !===========================================================================
    ! ANGLE CONVERSIONS
    !===========================================================================
    
    double precision, parameter :: PI         = 3.141592653589793d0
    double precision, parameter :: deg2RAD    = PI/180.0d0
    double precision, parameter :: RAD2deg    = 180.0d0/PI
    double precision, parameter :: ARCSEC2RAD = 4.84813681109536d-6
    double precision, parameter :: RAD2ARCSEC = 206264.806247096d0

    double precision, parameter :: RAD = 1d0
    double precision, parameter :: deg = deg2RAD
    
    !===========================================================================
    ! FUNDAMENTAL PHYSICAL CONSTANTS
    !===========================================================================
    
    ! Speed of light (CODATA 2018 recommended value)
    double precision, parameter :: SPEED_OF_LIGHT        = 299792458.0d0        ! m/s
    double precision, parameter :: SPEED_OF_LIGHT_UNIT   = 3.792265308366992d4 ! Ae/T
    
    ! Astronomical Unit
    double precision, parameter :: ASTRONOMICAL_UNIT_KM  = 1.49597870d8        ! km
    double precision, parameter :: AU_TO_EARTH_RADIUS    = ASTRONOMICAL_UNIT_KM/EARTH_RADIUS_EQUATORIAL
    
    ! Earth's dynamical form factor (J2)
    double precision, parameter :: EARTH_J2_CONSTANT     = 1.082629989d-3
    
    ! Gravitational constant (CODATA 2018)
    double precision, parameter :: GRAVITATIONAL_CONSTANT = 6.67430d-11        ! m^3/(kg·s^2)
    
    !===========================================================================
    ! ELEMENTS TYPES
    !===========================================================================
    integer, parameter :: STATE_TYPE_KEPLERIAN   = 1
    integer, parameter :: STATE_TYPE_CARTESIAN   = 2
    integer, parameter :: STATE_TYPE_NONSINGULAR_1 = 3  ! First kind of non-singular elements
    integer, parameter :: STATE_TYPE_NONSINGULAR_2 = 4  ! Second kind of non-singular elements / equinoctial elements 
    
    !===========================================================================
    ! OBJECT TYPES
    !===========================================================================
    integer, parameter :: OBJECT_TYPE_STATION    = 101   ! Ground station
    integer, parameter :: OBJECT_TYPE_SATELLITE  = 102   ! Satellite
    integer, parameter :: OBJECT_TYPE_DRONE = 103 
    integer, parameter :: OBJECT_TYPE_SHIP = 104
    integer, parameter :: OBJECT_TYPE_PLANE = 105
    integer, parameter :: OBJECT_TYPE_ROCKET = 106
    
    
end module module_const
