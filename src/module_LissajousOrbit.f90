!===============================================================================
! module_LissajousOrbit.f90
!
! Purpose: Define Lissajous orbit type for CRTBP systems
! Author: 
!
! Contents:
! - type_LissajousOrbit object definition for Lissajous orbits around Lagrange points
! - Linearized Lissajous orbit solution for collinear Lagrange points (L1, L2, L3)
!
! Lissajous Orbit Parameters:
! - name: Orbit name
! - LID: Lagrange point identifier (1=L1, 2=L2, 3=L3)
! - A_xy: In-plane amplitude (dimensionless)
! - A_z: Out-of-plane amplitude (dimensionless)
! - phi_xy: In-plane initial phase (radians)
! - phi_z: Out-of-plane initial phase (radians)
! - t0: Initial time (type_t)
! - crtbp_system: Pointer to CRTBP system

!===============================================================================

module module_LissajousOrbit
    use module_const
    use module_time
    use module_CRTBP
    implicit none
    
    !===========================================================================
    ! TYPE_LISSAJOUSORBIT OBJECT DEFINITION
    !===========================================================================
    
    type type_LissajousOrbit
        ! Orbit name
        character(len=50) :: name = ""
        
        ! Lagrange point identifier (1=L1, 2=L2, 3=L3)
        integer :: LID = LAGRANGE_L1
        
        ! Amplitudes (dimensionless)
        double precision :: A_xy = 0d0    ! In-plane amplitude
        double precision :: A_z = 0d0     ! Out-of-plane amplitude
        
        ! Initial phases (radians)
        double precision :: phi_xy = 0d0  ! In-plane initial phase
        double precision :: phi_z = 0d0   ! Out-of-plane initial phase
        
        ! Initial time
        type(type_t) :: t0
        
        ! Reference to CRTBP system (optional, can be null)
        type(type_CRTBP), pointer :: crtbp_system => null()
        
    contains
        ! Initialization procedure
        procedure :: init => initializeLissajous
        
        ! Get position and velocity at time t
        procedure :: getPV => getLissajousPV
        
        ! Get orbit parameters
        procedure :: getOrbitParams
    end type type_LissajousOrbit
    
contains
    
    !===========================================================================
    ! LISSAJOUS ORBIT INITIALIZATION
    !===========================================================================
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize a Lissajous orbit object
    !
    ! Parameters:
    !   this        : type_LissajousOrbit object to initialize
    !   name_in     : orbit name
    !   LID_in      : Lagrange point identifier (1=L1, 2=L2, 3=L3)
    !   A_xy_in     : in-plane amplitude (dimensionless)
    !   A_z_in      : out-of-plane amplitude (dimensionless)
    !   phi_xy_in   : in-plane initial phase (radians)
    !   phi_z_in    : out-of-plane initial phase (radians)
    !   t0_in       : initial time (type_t)
    !   crtbp_ptr   : pointer to CRTBP system (optional)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine initializeLissajous(this, name_in, LID_in, A_xy_in, A_z_in, &
                                   phi_xy_in, phi_z_in, t0_in, crtbp_ptr)
        implicit none
        class(type_LissajousOrbit) :: this
        character(len=*), intent(in) :: name_in
        integer, intent(in) :: LID_in
        double precision, intent(in) :: A_xy_in, A_z_in, phi_xy_in, phi_z_in
        type(type_t), intent(in) :: t0_in
        type(type_CRTBP), intent(in), target, optional :: crtbp_ptr
        
        ! Set basic parameters
        this%name = trim(name_in)
        this%LID = LID_in
        this%A_xy = A_xy_in
        this%A_z = A_z_in
        this%phi_xy = phi_xy_in
        this%phi_z = phi_z_in
        this%t0 = t0_in
        
        ! Set CRTBP system pointer if provided
        if (present(crtbp_ptr)) then
            this%crtbp_system => crtbp_ptr
        else
            nullify(this%crtbp_system)
        endif
        
    end subroutine initializeLissajous
    
    !===========================================================================
    ! LISSAJOUS ORBIT POSITION AND VELOCITY CALCULATION
    !===========================================================================
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get position and velocity of Lissajous orbit at time t
    !
    ! Parameters:
    !   this : type_LissajousOrbit object
    !   t    : time (type_t)
    !   xi   : x-position relative to Lagrange point (output)
    !   eta  : y-position relative to Lagrange point (output)
    !   zeta : z-position relative to Lagrange point (output)
    !   xi_dot  : x-velocity relative to Lagrange point (output)
    !   eta_dot : y-velocity relative to Lagrange point (output)
    !   zeta_dot: z-velocity relative to Lagrange point (output)
    !
    ! Notes:
    ! - Returns position and velocity in rotating frame coordinates
    ! - Coordinates are relative to the Lagrange point (ξ = x - x₀)
    ! - Based on linearized equations of motion
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine getLissajousPV(this, t, xi, eta, zeta, xi_dot, eta_dot, zeta_dot)
        implicit none
        class(type_LissajousOrbit) :: this
        type(type_t), intent(in) :: t
        double precision, intent(out) :: xi, eta, zeta, xi_dot, eta_dot, zeta_dot
        
        double precision :: tau, omega_xy, omega_z, k_coeff
        double precision :: x0, C0, d2, d3, k
        logical :: success
        double precision :: cos_xy, sin_xy, cos_z, sin_z
        
        ! Check if CRTBP system is associated
        if (.not. associated(this%crtbp_system)) then
            write(*, *) "Error: CRTBP system not associated with Lissajous orbit"
            xi = 0d0; eta = 0d0; zeta = 0d0
            xi_dot = 0d0; eta_dot = 0d0; zeta_dot = 0d0
            return
        endif
        
        ! Get time difference in dimensionless units
        ! Note: We need to convert physical time to dimensionless time
        ! For simplicity, we'll assume t and t0 are in the same time system
        ! and use the difference in MJD days, scaled by TU
        tau = (t%val - this%t0%val) * 86400d0 / this%crtbp_system%TU
        
        ! Get Lagrange point parameters
        call this%crtbp_system%getLagrangePointParams(this%LID, x0, C0, d2, d3, k, success)
        
        if (.not. success) then
            write(*, *) "Error: Failed to get Lagrange point parameters"
            xi = 0d0; eta = 0d0; zeta = 0d0
            xi_dot = 0d0; eta_dot = 0d0; zeta_dot = 0d0
            return
        endif
        
        ! Set frequencies
        omega_xy = d2  ! in-plane frequency
        omega_z = d3   ! out-of-plane frequency
        k_coeff = k    ! coupling coefficient
        
        ! Calculate trigonometric functions
        cos_xy = cos(omega_xy * tau + this%phi_xy)
        sin_xy = sin(omega_xy * tau + this%phi_xy)
        cos_z = cos(omega_z * tau + this%phi_z)
        sin_z = sin(omega_z * tau + this%phi_z)
        
        ! Position (linearized solution)
        xi = this%A_xy * cos_xy
        eta = k_coeff * this%A_xy * sin_xy
        zeta = this%A_z * cos_z
        
        ! Velocity (derivatives)
        xi_dot = -omega_xy * this%A_xy * sin_xy
        eta_dot = k_coeff * omega_xy * this%A_xy * cos_xy
        zeta_dot = -omega_z * this%A_z * sin_z
        
    end subroutine getLissajousPV
    
    !===========================================================================
    ! LISSAJOUS ORBIT PARAMETERS GETTER
    !===========================================================================
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get orbit parameters
    !
    ! Parameters:
    !   this : type_LissajousOrbit object
    !   omega_xy : in-plane frequency (output)
    !   omega_z  : out-of-plane frequency (output)
    !   k_coeff  : coupling coefficient (output)
    !   success  : whether parameters were obtained successfully (output)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine getOrbitParams(this, omega_xy, omega_z, k_coeff, success)
        implicit none
        class(type_LissajousOrbit) :: this
        double precision, intent(out) :: omega_xy, omega_z, k_coeff
        logical, intent(out) :: success
        
        double precision :: x0, C0, d2, d3, k
        
        ! Check if CRTBP system is associated
        if (.not. associated(this%crtbp_system)) then
            omega_xy = 0d0
            omega_z = 0d0
            k_coeff = 0d0
            success = .false.
            return
        endif
        
        ! Get Lagrange point parameters
        call this%crtbp_system%getLagrangePointParams(this%LID, x0, C0, d2, d3, k, success)
        
        if (success) then
            omega_xy = d2
            omega_z = d3
            k_coeff = k
        else
            omega_xy = 0d0
            omega_z = 0d0
            k_coeff = 0d0
        endif
        
    end subroutine getOrbitParams
    
end module module_LissajousOrbit