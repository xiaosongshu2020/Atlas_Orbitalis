!===============================================================================
! module_CRTBP.f90
!
! Purpose: Define Circular Restricted Three-Body Problem (CRTBP) type
! Author: 
!
! Contents:
! - type_CRTBP object definition with masses, distance, and derived parameters
! - Unit system: a12 as length unit, TU as time unit, (m1+m2) as mass unit
!
! CRTBP Parameters:
! - name: System name
! - m1: Primary body mass (kg)
! - m2: Secondary body mass (kg)  
! - a12: Distance between two primary bodies (km)
! - TU: Time unit (period) calculated from Kepler's third law
! - mu: m2/(m1+m2) (dimensionless mass ratio)
!
! Modern Fortran features used:
! - IMPLICIT NONE for type safety
! - Type-bound procedures for object-oriented design
! - Descriptive naming conventions
!===============================================================================

module module_CRTBP
    use module_const
    implicit none
    
    !===========================================================================
    ! LAGRANGE POINT IDENTIFIERS
    !===========================================================================
    integer, parameter :: LAGRANGE_L1 = 1
    integer, parameter :: LAGRANGE_L2 = 2
    integer, parameter :: LAGRANGE_L3 = 3
    integer, parameter :: LAGRANGE_L4 = 4
    integer, parameter :: LAGRANGE_L5 = 5
    
    !===========================================================================
    ! TYPE_CRTBP OBJECT DEFINITION
    !===========================================================================
    
    type type_CRTBP
        ! System name
        character(len=50) :: name = ""
        
        ! Primary body mass (kg)
        double precision :: m1 = 0d0
        
        ! Secondary body mass (kg)
        double precision :: m2 = 0d0
        
        ! Distance between two primary bodies (km)
        double precision :: a12 = 0d0
        
        ! Time unit (period) - seconds (calculated from Kepler's third law)
        double precision :: TU = 0d0
        
        ! Derived parameter: mu = m2/(m1+m2) (dimensionless)
        double precision :: mu = 0d0
        
        ! Lagrange points (libration points) in rotating frame coordinates
        ! L1, L2, L3, L4, L5: each is a 2-element array (x, y)
        double precision :: L1(2) = 0d0
        double precision :: L2(2) = 0d0
        double precision :: L3(2) = 0d0
        double precision :: L4(2) = 0d0
        double precision :: L5(2) = 0d0
        
    contains
        ! Initialization procedure
        procedure :: init => initializeCRTBP
        
        ! Lagrange point characteristic parameters
        procedure :: getLagrangePointCoord
        procedure :: getLagrangePointParams

        ! Characteristic distance gamma of a collinear libration point
        procedure :: getGamma
    end type type_CRTBP
    
contains
    
    !===========================================================================
    ! INITIALIZATION SUBROUTINE
    !===========================================================================
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize a CRTBP object with name, masses, and distance
    !
    ! Parameters:
    !   this    : type_CRTBP object to initialize
    !   name_in : system name
    !   m1_in   : primary body mass (kg)
    !   m2_in   : secondary body mass (kg)
    !   a12_in  : distance between two primary bodies (km)
    !
    ! Notes:
    ! - TU is calculated using Kepler's third law:
    !   TU = 2π * sqrt(a12^3 / (G * (m1+m2)))
    !   where G is gravitational constant
    ! - mu is calculated as m2/(m1+m2)
    ! - Lagrange points (L1-L5) are calculated in rotating frame coordinates
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine initializeCRTBP(this, name_in, m1_in, m2_in, a12_in)
        implicit none
        class(type_CRTBP) :: this
        character(len=*), intent(in) :: name_in
        double precision, intent(in) :: m1_in, m2_in, a12_in
        
        double precision :: total_mass, a12_meters, G
        double precision :: sqrt3_2, xL1, xL2, xL3
        logical :: converged
        
        ! Set basic parameters
        this%name = trim(name_in)
        this%m1 = m1_in
        this%m2 = m2_in
        this%a12 = a12_in
        
        ! Calculate total mass
        total_mass = m1_in + m2_in
        
        ! Calculate mass ratio mu = m2/(m1+m2)
        if (total_mass > 0d0) then
            this%mu = m2_in / total_mass
        else
            this%mu = 0d0
        endif
        
        ! Calculate TU using Kepler's third law
        ! TU = 2π * sqrt(a^3 / (G * M))
        ! where:
        !   a = a12 in meters
        !   G = gravitational constant (m^3/(kg·s^2))
        !   M = total mass (kg)
        
        if (total_mass > 0d0 .and. a12_in > 0d0) then
            ! Convert a12 from km to meters
            a12_meters = a12_in * 1000d0
            
            ! Gravitational constant
            G = GRAVITATIONAL_CONSTANT
            
            ! Calculate TU in seconds
            this%TU = 2d0 * PI * sqrt(a12_meters**3 / (G * total_mass))
        else
            this%TU = 0d0
        endif
        
        ! Calculate Lagrange points (libration points) in rotating frame
        !
        ! L4, L5 have the exact equilateral solution.
        ! L1, L2, L3 are found by Newton iteration on the collinear condition
        ! (see solveCollinearPoint). The truncated series of the classical
        ! expansion is used only as the initial guess:
        !   gamma1,2 = (mu/3)^(1/3) * [1 -/+ (1/3)(mu/3)^(1/3) - (1/9)(mu/3)^(2/3)]
        !   gamma3   ~ 1 - 7*mu/12
        ! For the Earth-Moon system the series for gamma is accurate to only
        ! ~0.3%, so the Newton solution is used instead of the series itself.

        if (this%mu > 0d0 .and. this%mu < 1d0) then
            ! --- L1, L2, L3: collinear points ---
            call solveCollinearPoint(this%mu, LAGRANGE_L1, xL1, converged)
            if (.not. converged) write(*, *) "Warning: L1 iteration did not converge"
            call solveCollinearPoint(this%mu, LAGRANGE_L2, xL2, converged)
            if (.not. converged) write(*, *) "Warning: L2 iteration did not converge"
            call solveCollinearPoint(this%mu, LAGRANGE_L3, xL3, converged)
            if (.not. converged) write(*, *) "Warning: L3 iteration did not converge"

            this%L1 = (/ xL1, 0d0 /)
            this%L2 = (/ xL2, 0d0 /)
            this%L3 = (/ xL3, 0d0 /)

            ! --- L4, L5: equilateral points ---
            sqrt3_2 = sqrt(3d0)/2d0

            this%L4(1) = 0.5d0 - this%mu
            this%L4(2) = sqrt3_2

            this%L5(1) = 0.5d0 - this%mu
            this%L5(2) = -sqrt3_2
        else
            ! If mu is invalid, set all Lagrange points to zero
            this%L1 = 0d0
            this%L2 = 0d0
            this%L3 = 0d0
            this%L4 = 0d0
            this%L5 = 0d0
        endif
        
    end subroutine initializeCRTBP
    
    !===========================================================================
    ! LAGRANGE POINT COORDINATE GETTER
    !===========================================================================
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get coordinates of a specific Lagrange point
    !
    ! Parameters:
    !   this : type_CRTBP object
    !   LID  : Lagrange point identifier (1=L1, 2=L2, 3=L3, 4=L4, 5=L5)
    !   x    : x-coordinate of Lagrange point (output)
    !   y    : y-coordinate of Lagrange point (output)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine getLagrangePointCoord(this, LID, x, y)
        implicit none
        class(type_CRTBP) :: this
        integer, intent(in) :: LID
        double precision, intent(out) :: x, y
        
        select case(LID)
        case(LAGRANGE_L1)
            x = this%L1(1)
            y = this%L1(2)
        case(LAGRANGE_L2)
            x = this%L2(1)
            y = this%L2(2)
        case(LAGRANGE_L3)
            x = this%L3(1)
            y = this%L3(2)
        case(LAGRANGE_L4)
            x = this%L4(1)
            y = this%L4(2)
        case(LAGRANGE_L5)
            x = this%L5(1)
            y = this%L5(2)
        case default
            x = 0d0
            y = 0d0
            write(*, *) "Warning: Invalid Lagrange point ID: ", LID
        end select
        
    end subroutine getLagrangePointCoord
    
    !===========================================================================
    ! LAGRANGE POINT CHARACTERISTIC PARAMETERS
    !===========================================================================
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get characteristic parameters for collinear Lagrange points (L1, L2, L3)
    !
    ! Parameters:
    !   this : type_CRTBP object
    !   LID  : Lagrange point identifier (1=L1, 2=L2, 3=L3)
    !   x0   : x-coordinate of Lagrange point (output)
    !   C0   : C0 parameter (output)
    !   d2   : in-plane frequency (output)
    !   d3   : out-of-plane frequency (output)
    !   k    : coupling coefficient k = -a2 (output)
    !   success : whether calculation was successful (output)
    !
    ! Notes:
    ! - Based on formulas from Lissajous orbit design algorithm
    ! - Only valid for collinear points (L1, L2, L3)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine getLagrangePointParams(this, LID, x0, C0, d2, d3, k, success)
        implicit none
        class(type_CRTBP) :: this
        integer, intent(in) :: LID
        double precision, intent(out) :: x0, C0, d2, d3, k
        logical, intent(out) :: success
        
        double precision :: r1, r2, sign1, sign2
        double precision :: dummy  ! dummy variable for y-coordinate
        
        ! Initialize outputs
        x0 = 0d0
        C0 = 0d0
        d2 = 0d0
        d3 = 0d0
        k = 0d0
        success = .false.
        
        ! Check if mu is valid
        if (this%mu <= 0d0 .or. this%mu >= 1d0) then
            write(*, *) "Warning: Invalid mu value for Lagrange point parameters: ", this%mu
            return
        endif
        
        ! Get Lagrange point coordinates
        call this%getLagrangePointCoord(LID, x0, dummy)
        
        ! Determine signs for r1 and r2 based on Lagrange point type
        ! r1 = distance to primary (mass m1 at x = -mu)
        ! r2 = distance to secondary (mass m2 at x = 1-mu)
        select case(LID)
        case(LAGRANGE_L1)  ! L1: between primary and secondary
            sign1 = 1d0   ! x0 + mu > 0
            sign2 = -1d0  ! x0 - 1 + mu < 0
        case(LAGRANGE_L2)  ! L2: beyond secondary
            sign1 = 1d0   ! x0 + mu > 0
            sign2 = 1d0   ! x0 - 1 + mu > 0
        case(LAGRANGE_L3)  ! L3: beyond primary
            sign1 = -1d0  ! x0 + mu < 0
            sign2 = -1d0  ! x0 - 1 + mu < 0
        case default
            write(*, *) "Warning: getLagrangePointParams only valid for L1, L2, L3"
            return
        end select
        
        ! Calculate C0 = (1-μ)/|x₀+μ|³ + μ/|x₀-1+μ|³
        r1 = abs(x0 + this%mu)
        r2 = abs(x0 - 1d0 + this%mu)
        
        if (r1 < 1d-12 .or. r2 < 1d-12) then
            write(*, *) "Warning: Distance to primary or secondary too small"
            return
        endif
        
        C0 = (1d0 - this%mu)/(r1**3) + this%mu/(r2**3)
        
        ! Calculate characteristic frequencies
        ! d₂ = √[ (√(9C₀²-8C₀) - C₀ + 2)/2 ]
        ! d₃ = √C₀
        if (C0 <= 0d0) then
            write(*, *) "Warning: Invalid C0 value: ", C0
            return
        endif
        
        d3 = sqrt(C0)
        
        ! Check if 9C₀²-8C₀ is positive
        if (9d0*C0**2 - 8d0*C0 < 0d0) then
            write(*, *) "Warning: 9C0^2 - 8C0 is negative: ", 9d0*C0**2 - 8d0*C0
            return
        endif
        
        d2 = sqrt((sqrt(9d0*C0**2 - 8d0*C0) - C0 + 2d0)/2d0)
        
        ! Calculate coupling coefficient k = -a₂
        ! a₂ = (d₂² + 2C₀ + 1)/(2d₂)
        k = -(d2**2 + 2d0*C0 + 1d0)/(2d0*d2)
        
        success = .true.

    end subroutine getLagrangePointParams

    !===========================================================================
    ! CHARACTERISTIC DISTANCE (GAMMA) OF A COLLINEAR LIBRATION POINT
    !===========================================================================

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get the characteristic distance gamma of a collinear libration point
    !
    ! gamma_i is the distance from L_i to the nearest primary and is the length
    ! unit used by the Lissajous and Halo analytical solutions
    ! (x = x_i + gamma*xi, y = gamma*eta, z = gamma*zeta):
    !
    !   L1 (between the primaries):  gamma1 = (1 - mu) - x1
    !   L2 (beyond the secondary) :  gamma2 = x2 - (1 - mu)
    !   L3 (beyond the primary)   :  gamma3 = -(x3 + mu) = |x3 + mu|
    !
    ! Parameters:
    !   this    : type_CRTBP object
    !   LID     : Lagrange point identifier (1=L1, 2=L2, 3=L3)
    !   gamma   : characteristic distance (output, dimensionless)
    !   success : whether gamma was obtained successfully (output)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine getGamma(this, LID, gamma, success)
        implicit none
        class(type_CRTBP) :: this
        integer, intent(in) :: LID
        double precision, intent(out) :: gamma
        logical, intent(out) :: success

        double precision :: x0, y0

        gamma = 0d0
        success = .false.

        ! Check if mu is valid
        if (this%mu <= 0d0 .or. this%mu >= 1d0) then
            write(*, *) "Warning: Invalid mu value for gamma: ", this%mu
            return
        endif

        ! Get Lagrange point coordinates
        call this%getLagrangePointCoord(LID, x0, y0)

        select case(LID)
        case(LAGRANGE_L1)
            gamma = (1d0 - this%mu) - x0
        case(LAGRANGE_L2)
            gamma = x0 - (1d0 - this%mu)
        case(LAGRANGE_L3)
            gamma = -(x0 + this%mu)
        case default
            write(*, *) "Warning: gamma only defined for collinear points L1, L2, L3"
            return
        end select

        if (gamma <= 0d0) then
            write(*, *) "Warning: Non-positive gamma for LID = ", LID
            success = .false.
            return
        endif

        success = .true.

    end subroutine getGamma

    !===========================================================================
    ! COLLINEAR LIBRATION POINT SOLVER (NEWTON)
    !===========================================================================

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Solve the collinear libration point equation by Newton's method
    !
    ! With the primary of mass m1 at x = -mu and the secondary of mass m2 at
    ! x = 1 - mu, the collinear point satisfies the vanishing of the gradient
    ! of the effective potential:
    !
    !   f(x)  = x - (1-mu)(x+mu)/|x+mu|^3 - mu(x-1+mu)/|x-1+mu|^3 = 0
    !   f'(x) = 1 + 2(1-mu)/|x+mu|^3 + 2mu/|x-1+mu|^3
    !
    ! Parameters:
    !   mu      : mass ratio m2/(m1+m2)
    !   LID     : 1=L1, 2=L2, 3=L3
    !   x       : x-coordinate of the libration point (output)
    !   success : whether the iteration converged (output)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine solveCollinearPoint(mu, LID, x, success)
        implicit none
        double precision, intent(in) :: mu
        integer, intent(in) :: LID
        double precision, intent(out) :: x
        logical, intent(out) :: success

        integer, parameter :: maxIter = 100
        double precision, parameter :: tol = 1d-14
        double precision :: m13, gamma3_guess, r1, r2, d1, d2, f, df, dx
        integer :: i

        ! Initial guess from the classical series expansion
        m13 = (mu/3d0)**(1d0/3d0)

        select case(LID)
        case(LAGRANGE_L1)
            x = (1d0 - mu) - m13*(1d0 - m13/3d0 - m13**2/9d0)
        case(LAGRANGE_L2)
            x = (1d0 - mu) + m13*(1d0 + m13/3d0 - m13**2/9d0)
        case(LAGRANGE_L3)
            gamma3_guess = 1d0 - 7d0*mu/12d0
            x = -mu - gamma3_guess
        case default
            x = 0d0
            success = .false.
            write(*, *) "Warning: solveCollinearPoint only valid for L1, L2, L3"
            return
        end select

        ! Newton iteration
        success = .false.
        do i = 1, maxIter
            r1 = x + mu
            r2 = x - 1d0 + mu
            d1 = abs(r1)**3
            d2 = abs(r2)**3

            if (d1 < 1d-12 .or. d2 < 1d-12) then
                write(*, *) "Warning: Collinear point iteration hit a primary body"
                return
            endif

            f  = x - (1d0 - mu)*r1/d1 - mu*r2/d2
            df = 1d0 + 2d0*(1d0 - mu)/d1 + 2d0*mu/d2
            dx = f/df
            x  = x - dx

            if (abs(dx) < tol) then
                success = .true.
                return
            endif
        enddo

        write(*, *) "Warning: Collinear point iteration did not converge in ", maxIter, " steps"

    end subroutine solveCollinearPoint

end module module_CRTBP
