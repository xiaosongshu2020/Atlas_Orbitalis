!===============================================================================
! module_halo.f90
!
! Purpose: Define Halo orbit type for CRTBP systems
!
! Contents:
! - type_halo object definition for Halo orbits around collinear libration points
! - Third-order analytical Halo orbit solution (Richardson 1980) for L1 and L2
!
! Halo orbit parameters:
! - name         : orbit name
! - LID          : libration point identifier (LAGRANGE_L1 or LAGRANGE_L2)
! - beta         : out-of-plane amplitude, in units of gamma_i
! - phi0         : initial phase (radians)
! - t0           : initial time (type_t)
! - crtbp_system : pointer to the CRTBP system
!
! Units and coordinates:
! - All internal calculations use canonical CRTBP units (a12 = 1, TU = 1).
! - Coordinates are relative to L_i and normalized by gamma_i, the distance
!   from L_i to the nearest primary:
!       x = x_i + gamma*xi,   y = gamma*eta,   z = gamma*zeta
! - beta, and the alpha derived from it, follow the same normalization; they
!   are the amplitudes in units of gamma_i (the convention of the reference
!   note "高阶分析解的构造_晕轨道的形成").
! - The Halo family has a single parameter: alpha follows from the amplitude
!   constraint   l1*alpha^2 + l2*beta^2 + (omega0^2 - c2) = 0
!
! Conventions - do not change casually:
! - k = +alpha_2 = (omega0^2 + 2*c2 + 1)/(2*omega0)  [Halo convention, (6.38)],
!   whereas type_CRTBP%getLagrangePointParams returns k = -alpha_2
!   [Lissajous convention]. c2 is recomputed here from gamma, so the two
!   conventions cannot be mixed up.
! - The sin(tau) term of eta, [b33*alpha^3 + (b34-b35)*alpha*beta^2]*sin(tau),
!   must be kept. It is missing from the printed solution (20b) of Richardson
!   (1980), for which the source document supplies the coefficients. Those three
!   coefficients are the only ones NOT checked against the paper (they do not
!   appear there); everything else matches Table I of the paper to 4e-6.
! - All coefficients follow Richardson (1980), Appendix I, including s2 and the
!   frequency, which differ from appendix 2 of the source document (see the
!   comments at Steps 9 and 13). L3 uses c_n from (6.24); the (6.27) of the
!   source document carries a spurious (-1)^n, and Richardson's (8b) plus the
!   L3 column of Table I (c3 = +1) confirm (6.24).
!===============================================================================

module module_halo
    use module_const
    use module_time
    use module_CRTBP
    implicit none

    private
    public :: type_halo

    !===========================================================================
    ! TYPE_HALO OBJECT DEFINITION
    !===========================================================================

    type type_halo
        ! --- orbit definition ---
        character(len=50) :: name = ""              ! orbit name
        integer :: LID = LAGRANGE_L2                ! 1=L1, 2=L2, 3=L3 (unverified)
        double precision :: beta = 0d0              ! out-of-plane amplitude (units of gamma)
        double precision :: phi0 = 0d0              ! initial phase (radians)
        type(type_t) :: t0                          ! initial time
        type(type_CRTBP), pointer :: crtbp_system => null()

        ! --- geometry and potential coefficients ---
        double precision :: gamma = 0d0             ! distance from L_i to nearest primary
        double precision :: r1 = 0d0, r2 = 0d0      ! distances to primary / secondary
        double precision :: c2 = 0d0, c3 = 0d0, c4 = 0d0

        ! --- frequencies and auxiliary parameters ---
        double precision :: omega0 = 0d0            ! linear in-plane frequency
        double precision :: k = 0d0                 ! shape parameter, +alpha_2
        double precision :: e1 = 0d0, e2 = 0d0      ! denominators of 2nd/3rd order terms

        ! --- second order coefficients ---
        double precision :: a21 = 0d0, a22 = 0d0, a23 = 0d0, a24 = 0d0
        double precision :: b21 = 0d0, b22 = 0d0, d21 = 0d0

        ! --- third order coefficients ---
        double precision :: a31 = 0d0, a32 = 0d0, b31 = 0d0, b32 = 0d0
        double precision :: d31 = 0d0, d32 = 0d0, b33 = 0d0, b34 = 0d0, b35 = 0d0

        ! --- frequency correction and amplitude constraint ---
        double precision :: s1 = 0d0, s2 = 0d0, a1 = 0d0, a2 = 0d0, l1 = 0d0, l2 = 0d0

        ! --- solution ---
        double precision :: alpha = 0d0             ! in-plane amplitude (units of gamma)
        double precision :: omega = 0d0             ! orbit frequency (nonlinear correction included)
        double precision :: period = 0d0            ! 2*pi/omega, canonical time units

        logical :: initialized = .false.

    contains
        ! Initialization procedure (computes all coefficients and alpha)
        procedure :: init => initializeHalo

        ! Get position and velocity at time t (type_t)
        procedure :: getPV => getHaloPV

        ! Get position and velocity at canonical time t measured from t0
        procedure :: compute_state => computeHaloState

        ! Get orbit parameters
        procedure :: getOrbitParams => getHaloParams

        ! Show orbit parameters
        procedure :: show => showHalo
    end type type_halo

contains

    !===========================================================================
    ! HALO ORBIT INITIALIZATION
    !===========================================================================

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize a Halo orbit object
    !
    ! Parameters:
    !   this      : type_halo object to initialize
    !   crtbp_ptr : pointer to the CRTBP system (required)
    !   LID_in    : libration point identifier (LAGRANGE_L1 or LAGRANGE_L2)
    !   beta_in   : out-of-plane amplitude, in units of gamma_i
    !   phi0_in   : initial phase (radians)
    !   t0_in     : initial time (type_t)
    !   name_in   : orbit name (optional)
    !   success   : whether the orbit could be constructed (optional output)
    !
    ! Notes:
    ! - All coefficients and the in-plane amplitude alpha are computed here,
    !   so an initialized object is self-contained (except for the CRTBP
    !   pointer, which is needed later to convert t to canonical time).
    ! - beta must be small enough for the amplitude constraint to have a real
    !   solution; otherwise success is false and a warning is written.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine initializeHalo(this, crtbp_ptr, LID_in, beta_in, phi0_in, t0_in, &
                              name_in, success)
        implicit none
        class(type_halo) :: this
        type(type_CRTBP), intent(in), target :: crtbp_ptr
        integer, intent(in) :: LID_in
        double precision, intent(in) :: beta_in, phi0_in
        type(type_t), intent(in) :: t0_in
        character(len=*), intent(in), optional :: name_in
        logical, intent(out), optional :: success

        double precision :: mu, term, denom

        logical :: ok

        ! Initialize status
        this%initialized = .false.
        if (present(success)) success = .false.

        ! Store inputs
        this%crtbp_system => crtbp_ptr
        this%LID = LID_in
        this%beta = beta_in
        this%phi0 = phi0_in
        this%t0 = t0_in

        if (present(name_in)) then
            this%name = trim(name_in)
        else
            write(this%name, '(A,I1)') "Halo L", LID_in
        endif

        ! Check the libration point identifier
        select case(LID_in)
        case(LAGRANGE_L1, LAGRANGE_L2)
            ! collinear points supported by the third-order solution
        case(LAGRANGE_L3)
            ! c_n from (6.24); see the module header about the (6.27) typo.
            ! Coefficients verified against Richardson (1980), Table I.
        case default
            write(*, *) "Error: Halo orbits only exist around collinear points L1, L2, L3"
            return
        end select

        mu = crtbp_ptr%mu

        ! ---------------------------------------------------------------------
        ! Step 1: distance gamma from L_i to the nearest primary, r1 and r2
        ! ---------------------------------------------------------------------
        call crtbp_ptr%getGamma(LID_in, this%gamma, ok)
        if (.not. ok) then
            write(*, *) "Error: Cannot obtain gamma for LID = ", LID_in
            return
        endif

        select case(LID_in)
        case(LAGRANGE_L1)   ! r1 = 1 - gamma, r2 = gamma
            this%r1 = 1d0 - this%gamma
            this%r2 = this%gamma
        case(LAGRANGE_L2)   ! r1 = 1 + gamma, r2 = gamma
            this%r1 = 1d0 + this%gamma
            this%r2 = this%gamma
        case(LAGRANGE_L3)   ! r1 = gamma, r2 = 1 + gamma
            this%r1 = this%gamma
            this%r2 = 1d0 + this%gamma
        end select

        ! ---------------------------------------------------------------------
        ! Step 2: potential coefficients c2, c3, c4 (Legendre expansion)
        ! ---------------------------------------------------------------------
        this%c2 = cN(mu, this%gamma, LID_in, 2)
        this%c3 = cN(mu, this%gamma, LID_in, 3)
        this%c4 = cN(mu, this%gamma, LID_in, 4)

        ! ---------------------------------------------------------------------
        ! Step 3: basic frequency omega0 and shape parameter k
        ! ---------------------------------------------------------------------
        this%omega0 = sqrt( ( sqrt(9d0*this%c2**2 - 8d0*this%c2) + 2d0 - this%c2 ) / 2d0 )
        this%k = ( this%omega0**2 + 2d0*this%c2 + 1d0 ) / ( 2d0*this%omega0 )

        ! ---------------------------------------------------------------------
        ! Step 4: auxiliary parameters e1, e2
        ! ---------------------------------------------------------------------
        this%e1 = (3d0*this%omega0**2 / this%k) * &
                  ( this%k*(6d0*this%omega0**2 - 1d0) - 2d0*this%omega0 )
        this%e2 = (8d0*this%omega0**2 / this%k) * &
                  ( this%k*(11d0*this%omega0**2 - 1d0) - 2d0*this%omega0 )

        if (abs(this%e1) < 1d-10 .or. abs(this%e2) < 1d-10) then
            write(*, *) "Warning: Singular auxiliary parameters e1, e2: ", this%e1, this%e2
            return
        endif

        ! ---------------------------------------------------------------------
        ! Step 5: second order coefficients
        ! ---------------------------------------------------------------------
        this%a21 = ( 3d0*this%c3*(this%k**2 - 2d0) ) / ( 4d0*(1d0 + 2d0*this%c2) )
        this%a22 = ( 3d0*this%c3 ) / ( 4d0*(1d0 + 2d0*this%c2) )

        this%b21 = -( 3d0*this%c3*this%omega0 ) / ( 2d0*this%e1 ) * &
                   ( 3d0*this%k*this%omega0 - 4d0 )
        this%b22 =  ( 3d0*this%c3*this%omega0 ) / this%e1

        this%a23 = -( 3d0*this%c3*this%omega0 ) / ( 4d0*this%k*this%e1 ) * &
                   ( 3d0*this%k**3*this%omega0 - 6d0*this%k*(this%k - this%omega0) + 4d0 )
        this%a24 = -( 3d0*this%c3*this%omega0 ) / ( 4d0*this%k*this%e1 ) * &
                   ( 2d0 + 3d0*this%k*this%omega0 )

        this%d21 = -this%c3 / ( 2d0*this%omega0**2 )

        ! ---------------------------------------------------------------------
        ! Step 6: third order radial coefficients a31, a32
        ! ---------------------------------------------------------------------
        this%a31 = -( 9d0*this%omega0 ) / ( 4d0*this%e2 ) * &
                   ( 4d0*this%c3*(this%k*this%a23 - this%b21) + this%k*this%c4*(4d0 + this%k**2) ) &
                   + ( 9d0*this%omega0**2 + 1d0 - this%c2 ) / ( 2d0*this%e2 ) * &
                   ( 3d0*this%c3*(2d0*this%a23 - this%k*this%b21) + this%c4*(2d0 + 3d0*this%k**2) )

        this%a32 = -1d0 / this%e2 * ( &
                   ( 9d0*this%omega0 ) / 4d0 * &
                   ( 4d0*this%c3*(this%k*this%a24 - this%b22) + this%k*this%c4 ) &
                   + ( 3d0/2d0 )*( 9d0*this%omega0**2 + 1d0 - this%c2 ) * &
                   ( this%c3*(this%k*this%b22 + this%d21 - 2d0*this%a24) - this%c4 ) )

        ! ---------------------------------------------------------------------
        ! Step 7: third order vertical coefficients b31, b32
        ! ---------------------------------------------------------------------
        this%b31 = ( 3d0 / (8d0*this%e2) ) * ( &
                   8d0*this%omega0 * &
                   ( 3d0*this%c3*(this%k*this%b21 - 2d0*this%a23) - this%c4*(2d0 + 3d0*this%k**2) ) &
                   + ( 9d0*this%omega0**2 + 1d0 + 2d0*this%c2 ) * &
                   ( 4d0*this%c3*(this%k*this%a23 - this%b21) + this%k*this%c4*(4d0 + this%k**2) ) )

        this%b32 = ( 1d0 / this%e2 ) * ( &
                   9d0*this%omega0 * &
                   ( this%c3*(this%k*this%b22 + this%d21 - 2d0*this%a24) - this%c4 ) &
                   + ( 3d0/8d0 )*( 9d0*this%omega0**2 + 1d0 + 2d0*this%c2 ) * &
                   ( 4d0*this%c3*(this%k*this%a24 - this%b22) + this%k*this%c4 ) )

        ! ---------------------------------------------------------------------
        ! Step 8: auxiliary parameters d31, d32
        ! ---------------------------------------------------------------------
        this%d31 = ( 3d0 / (64d0*this%omega0**2) ) * ( 4d0*this%c3*this%a24 + this%c4 )
        this%d32 = ( 3d0 / (64d0*this%omega0**2) ) * &
                   ( 4d0*this%c3*(this%a23 - this%d21) + this%c4*(4d0 + this%k**2) )

        ! ---------------------------------------------------------------------
        ! Step 9: frequency correction coefficients s1, s2
        ! ---------------------------------------------------------------------
        denom = 2d0*this%omega0*( this%omega0*(1d0 + this%k**2) - 2d0*this%k )

        if (abs(denom) < 1d-10) then
            write(*, *) "Warning: Singular denominator of s1, s2: ", denom
            return
        endif

        this%s1 = ( 1d0 / denom ) * ( &
                  ( 3d0/2d0 )*this%c3 * &
                  ( 2d0*this%a21*(this%k**2 - 2d0) - this%a23*(this%k**2 + 2d0) - 2d0*this%k*this%b21 ) &
                  - ( 3d0/8d0 )*this%c4*( 3d0*this%k**4 - 8d0*this%k**2 + 8d0 ) )

        ! NOTE: appendix 2 of the source document prints this bracket as
        !   ... - 2*k*b22 + 5*d21 ] + (3/8)*c4*(12 - k^4)
        ! Both details deviate from Richardson (1980), Appendix I, which is
        !   ... + 2*k*b22 + 5*d21 ] + (3/8)*c4*(12 - k^2)
        ! Richardson is authoritative here: with the printed (12 - k^4) the
        ! analytical state is a markedly poorer Halo initial guess. See the
        ! differential correction columns of test_halo:
        !   (12 - k^4): |F|1 / |F|1_linear = 1.5e-1 / 7.1e-1  (only ~5x better)
        !   (12 - k^2): |F|1 / |F|1_linear = 1.2e-2 / 7.1e-1  (~58x, i.e. rho^2)
        this%s2 = ( 1d0 / denom ) * ( &
                  ( 3d0/2d0 )*this%c3 * &
                  ( 2d0*this%a22*(this%k**2 - 2d0) + this%a24*(this%k**2 + 2d0) &
                    + 2d0*this%k*this%b22 + 5d0*this%d21 ) &
                  + ( 3d0/8d0 )*this%c4*( 12d0 - this%k**2 ) )

        ! ---------------------------------------------------------------------
        ! Step 10: third order vertical coefficients b33, b34, b35
        ! ---------------------------------------------------------------------
        this%b33 = -this%k / (16d0*this%omega0) * ( &
                   12d0*this%c3*( this%b21 - 2d0*this%k*this%a21 + this%k*this%a23 ) &
                   + 3d0*this%k*this%c4*( 3d0*this%k**2 - 4d0 ) &
                   + 16d0*this%s1*this%omega0*( this%k*this%omega0 - 1d0 ) )

        this%b34 = this%k / (8d0*this%omega0) * ( &
                   12d0*this%k*this%c3*this%a22 - 3d0*this%k*this%c4 &
                   - 8d0*this%s2*this%omega0*( this%k*this%omega0 - 1d0 ) )

        this%b35 = -this%k / (16d0*this%omega0) * &
                   ( 12d0*this%c3*( this%b22 + this%k*this%a24 ) + 3d0*this%k*this%c4 )

        ! ---------------------------------------------------------------------
        ! Step 11: auxiliary parameters a1, a2 of the amplitude constraint
        ! ---------------------------------------------------------------------
        this%a1 = -( 3d0/2d0 )*this%c3*( 2d0*this%a21 + this%a23 + 5d0*this%d21 ) &
                  - ( 3d0/8d0 )*this%c4*( 12d0 - this%k**2 )

        this%a2 = ( 3d0/2d0 )*this%c3*( this%a24 - 2d0*this%a22 ) + ( 9d0/8d0 )*this%c4

        this%l1 = this%a1 + 2d0*this%omega0**2*this%s1
        this%l2 = this%a2 + 2d0*this%omega0**2*this%s2

        ! ---------------------------------------------------------------------
        ! Step 12: in-plane amplitude alpha from the amplitude constraint
        !          l1*alpha^2 + l2*beta^2 + (omega0^2 - c2) = 0
        ! ---------------------------------------------------------------------
        if (abs(this%l1) < 1d-14) then
            write(*, *) "Warning: Singular amplitude constraint coefficient l1 = ", this%l1
            return
        endif

        term = -( this%l2*this%beta**2 + this%omega0**2 - this%c2 ) / this%l1

        if (term < 0d0) then
            write(*, *) "Warning: beta too large for a Halo orbit (no real alpha): ", this%beta
            write(*, *) "         l1*alpha^2 + l2*beta^2 + (omega0^2 - c2) = 0 has no solution"
            return
        endif

        this%alpha = sqrt(term)

        ! ---------------------------------------------------------------------
        ! Step 13: orbit frequency and period
        !
        ! The frequency is lambda*(1 + s1*alpha^2 + s2*beta^2), NOT
        ! lambda + s1*alpha^2 + s2*beta^2 as printed in (6.38) of the source
        ! document. Richardson (1980) writes omega = 1 + sum(omega_n) with
        ! omega_2 = s1*Ax^2 + s2*Az^2, and the phase is tau_1 = lambda*omega*s,
        ! so the actual frequency is lambda*(1 + omega_2). Table I of the paper
        ! settles it: for Sun-Earth L1 with A_z = 125000 km it quotes a period
        ! of 177.704 days, which lambda*(1+omega_2) reproduces exactly while
        ! lambda + omega_2 gives 176.318 days.
        ! ---------------------------------------------------------------------
        this%omega = this%omega0*( 1d0 + this%s1*this%alpha**2 + this%s2*this%beta**2 )

        if (this%omega <= 0d0) then
            write(*, *) "Warning: Non-positive orbit frequency: ", this%omega
            return
        endif

        this%period = 2d0*PI / this%omega
        this%initialized = .true.

        if (present(success)) success = .true.

    end subroutine initializeHalo

    !===========================================================================
    ! HALO ORBIT POSITION AND VELOCITY (type_t INTERFACE)
    !===========================================================================

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get position and velocity of the Halo orbit at time t
    !
    ! Parameters:
    !   this : type_halo object
    !   t    : time (type_t); the orbit phase is measured from this%t0
    !   xi, eta, zeta          : position relative to L_i, in units of gamma_i
    !   xi_dot, eta_dot, zeta_dot : velocity, canonical units (per TU)
    !
    ! Notes:
    ! - The canonical time is  tau_c = (t - t0)[days] * 86400 / TU.
    ! - The returned vectors are the rotating (synodic) frame components
    !   relative to the libration point; to obtain the full synodic position:
    !       x = x_i + gamma*xi,  y = gamma*eta,  z = gamma*zeta
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine getHaloPV(this, t, xi, eta, zeta, xi_dot, eta_dot, zeta_dot)
        implicit none
        class(type_halo) :: this
        type(type_t), intent(in) :: t
        double precision, intent(out) :: xi, eta, zeta, xi_dot, eta_dot, zeta_dot

        double precision :: tau

        ! Check that the orbit can be evaluated
        if (.not. this%initialized .or. .not. associated(this%crtbp_system)) then
            write(*, *) "Error: Halo orbit is not initialized"
            xi = 0d0; eta = 0d0; zeta = 0d0
            xi_dot = 0d0; eta_dot = 0d0; zeta_dot = 0d0
            return
        endif

        if (this%crtbp_system%TU <= 0d0) then
            write(*, *) "Error: Invalid time unit TU of the CRTBP system"
            xi = 0d0; eta = 0d0; zeta = 0d0
            xi_dot = 0d0; eta_dot = 0d0; zeta_dot = 0d0
            return
        endif

        ! Physical time difference (days) -> canonical time
        tau = (t%val - this%t0%val) * 86400d0 / this%crtbp_system%TU

        call this%compute_state(tau, xi, eta, zeta, xi_dot, eta_dot, zeta_dot)

    end subroutine getHaloPV

    !===========================================================================
    ! HALO ORBIT POSITION AND VELOCITY (CANONICAL TIME INTERFACE)
    !===========================================================================

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get position and velocity of the Halo orbit at canonical time t
    !
    ! Parameters:
    !   this : type_halo object
    !   t    : canonical time measured from t0 (per TU)
    !   xi, eta, zeta             : position relative to L_i, in units of gamma_i
    !   xi_dot, eta_dot, zeta_dot : velocity, canonical units (per TU)
    !
    ! Notes:
    ! - Third-order analytical solution, with tau = omega*t + phi0:
    !     xi   = -alpha*cos(tau) + a21*alpha^2 + a22*beta^2
    !            + (a23*alpha^2 - a24*beta^2)*cos(2tau)
    !            + (a31*alpha^3 - a32*alpha*beta^2)*cos(3tau)
    !     eta  = [k*alpha + b33*alpha^3 + (b34-b35)*alpha*beta^2]*sin(tau)
    !            + (b21*alpha^2 - b22*beta^2)*sin(2tau)
    !            + (b31*alpha^3 - b32*alpha*beta^2)*sin(3tau)
    !     zeta = beta*cos(tau) + d21*alpha*beta*(cos(2tau) - 3)
    !            + (d32*beta*alpha^2 - d31*beta^3)*cos(3tau)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine computeHaloState(this, t, xi, eta, zeta, xi_dot, eta_dot, zeta_dot)
        implicit none
        class(type_halo), intent(in) :: this
        double precision, intent(in) :: t
        double precision, intent(out) :: xi, eta, zeta, xi_dot, eta_dot, zeta_dot

        double precision :: tau, alpha, beta, omega
        double precision :: ct, st, c2t, s2t, c3t, s3t
        double precision :: coeff_sin_eta   ! combined coefficient of sin(tau) in eta

        if (.not. this%initialized) then
            write(*, *) "Warning: Halo orbit is not initialized"
            xi = 0d0; eta = 0d0; zeta = 0d0
            xi_dot = 0d0; eta_dot = 0d0; zeta_dot = 0d0
            return
        endif

        alpha = this%alpha
        beta  = this%beta
        omega = this%omega

        ! Phase angle
        tau = omega*t + this%phi0

        ct  = cos(tau)
        st  = sin(tau)
        c2t = cos(2d0*tau)
        s2t = sin(2d0*tau)
        c3t = cos(3d0*tau)
        s3t = sin(3d0*tau)

        !------------------------------------------------------------------------
        ! Position
        !------------------------------------------------------------------------
        xi = -alpha*ct &
             + this%a21*alpha**2 + this%a22*beta**2 &
             + (this%a23*alpha**2 - this%a24*beta**2)*c2t &
             + (this%a31*alpha**3 - this%a32*alpha*beta**2)*c3t

        ! The two sin(tau) terms of eta are combined
        coeff_sin_eta = this%k*alpha &
                        + this%b33*alpha**3 &
                        + (this%b34 - this%b35)*alpha*beta**2

        eta = coeff_sin_eta*st &
              + (this%b21*alpha**2 - this%b22*beta**2)*s2t &
              + (this%b31*alpha**3 - this%b32*alpha*beta**2)*s3t

        zeta = beta*ct &
               + this%d21*alpha*beta*(c2t - 3d0) &
               + (this%d32*beta*alpha**2 - this%d31*beta**3)*c3t

        !------------------------------------------------------------------------
        ! Velocity (time derivative, canonical time)
        !------------------------------------------------------------------------
        xi_dot = omega * ( &
                 alpha*st &
                 - 2d0*(this%a23*alpha**2 - this%a24*beta**2)*s2t &
                 - 3d0*(this%a31*alpha**3 - this%a32*alpha*beta**2)*s3t )

        eta_dot = omega * ( &
                  coeff_sin_eta*ct &
                  + 2d0*(this%b21*alpha**2 - this%b22*beta**2)*c2t &
                  + 3d0*(this%b31*alpha**3 - this%b32*alpha*beta**2)*c3t )

        zeta_dot = omega * ( &
                   -beta*st &
                   - 2d0*this%d21*alpha*beta*s2t &
                   - 3d0*(this%d32*beta*alpha**2 - this%d31*beta**3)*s3t )

    end subroutine computeHaloState

    !===========================================================================
    ! HALO ORBIT PARAMETERS GETTER
    !===========================================================================

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get the characteristic parameters of the Halo orbit
    !
    ! Parameters:
    !   this    : type_halo object
    !   alpha   : in-plane amplitude, units of gamma_i (output)
    !   beta    : out-of-plane amplitude, units of gamma_i (output)
    !   omega   : orbit frequency, canonical units (output)
    !   period  : orbit period, canonical time units (output)
    !   gamma   : characteristic distance of L_i (output)
    !   success : whether parameters were obtained successfully (output)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine getHaloParams(this, alpha, beta, omega, period, gamma, success)
        implicit none
        class(type_halo) :: this
        double precision, intent(out) :: alpha, beta, omega, period, gamma
        logical, intent(out) :: success

        alpha = 0d0
        beta = 0d0
        omega = 0d0
        period = 0d0
        gamma = 0d0
        success = .false.

        if (.not. this%initialized) then
            write(*, *) "Warning: Halo orbit is not initialized"
            return
        endif

        alpha = this%alpha
        beta = this%beta
        omega = this%omega
        period = this%period
        gamma = this%gamma
        success = .true.

    end subroutine getHaloParams

    !===========================================================================
    ! HALO ORBIT PARAMETERS SUMMARY
    !===========================================================================

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print a summary of the Halo orbit
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine showHalo(this)
        implicit none
        class(type_halo) :: this

        if (.not. this%initialized) then
            write(*, *) "type_halo: orbit is not initialized"
            return
        endif

        write(*, '(A)')        "=============================================="
        write(*, '(A,A)')      " Halo orbit : ", trim(this%name)
        write(*, '(A)')        "=============================================="
        write(*, '(A,I2)')     " Libration point        : L", this%LID
        write(*, '(A,F12.8)')  " Mass ratio mu          : ", this%crtbp_system%mu
        write(*, '(A,F12.8)')  " gamma                  : ", this%gamma
        write(*, '(A,3F13.8)') " c2, c3, c4             : ", this%c2, this%c3, this%c4
        write(*, '(A,F12.8)')  " omega0 (linear)        : ", this%omega0
        write(*, '(A,F12.8)')  " omega  (nonlinear)     : ", this%omega
        write(*, '(A,F12.8)')  " Period (TU)            : ", this%period
        write(*, '(A,F12.6)')  " k                      : ", this%k
        write(*, '(A,F12.8)')  " beta  (out-of-plane)   : ", this%beta
        write(*, '(A,F12.8)')  " alpha (in-plane)       : ", this%alpha
        write(*, '(A,F12.4)')  " Period (days)          : ", &
            this%period*this%crtbp_system%TU/86400d0
        write(*, '(A)')        "=============================================="

    end subroutine showHalo

    !===========================================================================
    ! LEGENDRE EXPANSION COEFFICIENT
    !===========================================================================

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Coefficient c_n (n >= 2) of the Legendre expansion of the potential
    !
    !   L1, L2 :  c_n = [(-1)^n sign_mu * mu
    !                    + (-1)^n (1-mu) (gamma/(1 -/+ gamma))^(n+1)] / gamma^3
    !   L3     :  c_n = [(1-mu) + mu (gamma/(1+gamma))^(n+1)] / gamma^3
    !
    ! with the upper sign for L1 (sign_mu = +1) and the lower for L2.
    ! c_2 is identical to the classical C0 = (1-mu)/r1^3 + mu/r2^3.
    !
    ! Parameters:
    !   mu    : mass ratio
    !   gamma : characteristic distance of the libration point
    !   LID   : 1=L1, 2=L2, 3=L3
    !   n     : order of the coefficient (n >= 2)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pure function cN(mu, gamma, LID, n) result(cval)
        implicit none
        double precision, intent(in) :: mu, gamma
        integer, intent(in) :: LID, n
        double precision :: cval

        double precision :: sgn   ! (-1)^n

        ! (-1)^n as a real
        if (mod(n, 2) == 1) then
            sgn = -1d0
        else
            sgn = 1d0
        endif

        select case(LID)
        case(LAGRANGE_L1)
            cval = ( mu + sgn*(1d0 - mu)*(gamma/(1d0 - gamma))**(n + 1) ) / gamma**3
        case(LAGRANGE_L2)
            cval = sgn * ( mu + (1d0 - mu)*(gamma/(1d0 + gamma))**(n + 1) ) / gamma**3
        case(LAGRANGE_L3)
            cval = ( (1d0 - mu) + mu*(gamma/(1d0 + gamma))**(n + 1) ) / gamma**3
        case default
            cval = 0d0
        end select

    end function cN

end module module_halo
