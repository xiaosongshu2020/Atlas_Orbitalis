    module module_elements
    use module_toolKits
    implicit none

    !============================================================================
    ! type_orbit — encapsulated orbital elements with type-aware conversions
    !
    ! Element type codes (from module_const):
    !   1 = Keplerian (a, e, i, Omega, omega, M)
    !   2 = Cartesian   (x, y, z, vx, vy, vz)
    !   3 = NSG1 — first kind of non-singular elements (a, i, Omega, xi, eta, lambda)
    !   4 = NSG2 — second kind / equinoctial elements (a, h, k, p, q, lambda)
    !
    ! Elements are stored in NORMALIZED units internally.
    ! Use toSI() / toNormalized() to convert, or init(isNormalized=.false.) for SI input.
    !============================================================================
    type, public :: type_orbit
        double precision :: elements(6) = 0.d0
        integer          :: eleType   = 1    ! STATE_TYPE_KEPLERIAN
    contains
        procedure :: init             => orbit_init
        procedure :: toKeplerian      => orbit_toKeplerian
        procedure :: toCartesian      => orbit_toCartesian
        procedure :: toNSG1           => orbit_toNSG1
        procedure :: toNSG2           => orbit_toNSG2
        procedure :: toSI             => orbit_toSI
        procedure :: toNormalized     => orbit_toNormalized
        procedure :: get_rEf          => orbit_get_rEf
        procedure :: get_rut          => orbit_get_rut
        procedure :: propagateTwoBody => orbit_propagateTwoBody
        procedure :: show             => orbit_show
    end type type_orbit

    contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine kpl_SI2normal(x0,x1)
    use module_const
    implicit none
    double precision x0(6),x1(6)
    !-----------------------------
    x1(1) = x0(1)/length_unit
    x1(2) = x0(2)
    x1(3:6) = x0(3:6)*deg2rad
    endsubroutine
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine kpl_normal2SI(x1, x0)
    use module_const
    implicit none
    double precision x0(6), x1(6)
    !--------------------------------------
    x0(1) = x1(1)*length_unit
    x0(2) = x1(2)
    x0(3:6) = x1(3:6)*rad2deg
    endsubroutine
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine cts_SI2normal(x0, x1)
    use module_const
    implicit none
    double precision x0(6), x1(6)
    !-----------------------------
    x1(1:3) = x0(1:3)/length_unit
    x1(4:6) = x0(4:6)/length_unit*time_unit
    endsubroutine
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine cts_normal2SI(x1, x0)
    use module_const
    implicit none
    double precision x0(6), x1(6)
    !--------------------------------------
    x0(1:3) = x1(1:3)*length_unit
    x0(4:6) = x1(4:6)*length_unit/time_unit
    endsubroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine kpl2nsg1(elements, roots)
    implicit none
    double precision elements(6), roots(6)
    !--------------------------------------
    roots(1) = elements(1)
    roots(2) = elements(3)
    roots(3) = elements(4)
    roots(4) = elements(2)*dcos(elements(5))
    roots(5) = elements(2)*dsin(elements(5))
    roots(6) = elements(6) + elements(5)
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine kpl2nsg2(elements, roots)
    implicit none
    double precision elements(6), roots(6), omega_low_tiao
    !--------------------------------------
    roots(1) = elements(1)
    omega_low_tiao = elements(4) + elements(5)
    roots(2) = elements(2)*dcos(omega_low_tiao)
    roots(3) = elements(2)*dsin(omega_low_tiao)
    roots(4) = dsin(elements(3)/2d0)*dcos(elements(4))
    roots(5) = dsin(elements(3)/2d0)*dsin(elements(4))
    roots(6) = elements(6) + omega_low_tiao
    endsubroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine kpl2cts(elements, cts)
    implicit none
    !CAUTION: omega is lowercase omega and C_omega is uppercase omega
    double precision elements(6),a,e_element,i,omega,C_omega,M,M_mod,E,r_norm,pi_constant
    double precision r(3),v(3)
    double precision P(3),Q(3)
    double precision cts(6)
    integer counter
    parameter(pi_constant=3.141592653589793D0)

    a = elements(1)
    e_element = elements(2)
    i = elements(3)
    C_omega = elements(4)
    omega = elements(5)
    M = elements(6)

    P(1)=dcos(C_omega)*dcos(omega)-dsin(C_omega)*dsin(omega)*dcos(i)
    P(2)=dsin(C_omega)*dcos(omega)+dcos(C_omega)*dsin(omega)*dcos(i)
    P(3)=dsin(omega)*dsin(i)

    Q(1)=-dcos(C_omega)*dsin(omega)-dsin(C_omega)*dcos(omega)*dcos(i)
    Q(2)=-dsin(C_omega)*dsin(omega)+dcos(C_omega)*dcos(omega)*dcos(i)
    Q(3)=dcos(omega)*dsin(i)

    M_mod=modulo(M,2*pi_constant)
    E=M_mod
    do counter=1,5
        E=-(E-e_element*dsin(E)-M_mod)*1.0D0/(1-e_element*dcos(E))+E
    end do

    r=a*(dcos(E)-e_element)*P+a*(1D0-e_element**2D0)**(0.5D0)*dsin(E)*Q
    r_norm = (r(1)**2 + r(2)**2 + r(3)**2)**0.5

    v=a**0.5D0/r_norm*(-dsin(E)*P+(1D0-e_element**2D0)**0.5D0*dcos(E)*Q)

    cts(1:3) = r
    cts(4:6) = v
    return
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine nsg12kpl(roots, elements)
    implicit none
    double precision roots(6), elements(6)
    double precision pi
    !--------------------------------------
    pi = 3.141592653589793d0
    elements(1) = roots(1)
    elements(2) = dsqrt(roots(4)**2 + roots(5)**2)
    elements(3) = roots(2)
    elements(4) = roots(3)
    if (elements(2) == 0d0) then !if it is a absolutely circular orbit!
        elements(5) = 0d0
        elements(6) = roots(6)
        else
            elements(5) = dacos(roots(4)/elements(2))
            if (roots(5) < 0) elements(5) = 2d0*pi - elements(5)
            elements(6) = roots(6) - elements(5)
    end if
    elements(6) =   dmod(elements(6), 3.141592653589793d0*2d0)
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !When i = 0 or e = 0, this subroutine can't work well!
    !That is reasonable as the Omega and omega are meaningless
    !mathematically when the i and the e are zero.
    subroutine nsg22kpl(roots, elements)
    implicit none
    double precision roots(6), elements(6)
    double precision pi, omega_low_tiao
    !-------------------------------------
    pi = 3.142592653589793d0
    elements(1) = roots(1)
    elements(2) = dsqrt(roots(2)**2 + roots(3)**2)
    elements(3) = 2d0*dasin(dsqrt(roots(4)**2 + roots(5)**2))
    elements(4) = datan(roots(5)/roots(4))
    if (roots(4) < 0) then
        elements(4) = elements(4) + pi
    end if
    elements(4) = dmod(elements(4), 2d0*pi)
    call get_angle(roots(3)/elements(2), roots(2)/elements(2), omega_low_tiao)
    elements(5) = omega_low_tiao - elements(4)
    elements(6) = roots(6) - omega_low_tiao
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   subroutine cts2kpl(miu,r,v,a,e_element,i,C_omega,omega,M)
    subroutine cts2kpl(rv, elements)
    implicit none
    double precision rv(6), r(3), v(3), miu
    double precision r_norm,v_norm,a,e_element,i,omega
    double precision C_omega,M,E,f,u,P(3),Q(3),R_array(3), elements(6)
    !----------------------------------------------------------

    r = rv(1:3)
    v = rv(4:6)
    miu = 1d0

    v_norm = norm2(v)
    r_norm = norm2(r)

    a=1.0D0/(2.0D0/r_norm-v_norm**2/miu)                               
    e_element=((1.0D0-r_norm/a)**2+&
                    &(dot_product(r,v)*1.0D0/(miu*a)**0.5)**2)**0.5
    

    E=datan2(dot_product(r,v)*1.0D0/(miu*a)**0.5,1.0D0-r_norm/a)
    f=datan2((1D0-e_element**2D0)**0.5D0*dsin(E),dcos(E)-e_element)
    M=E-e_element*dsin(E)

    P=dcos(E)/r_norm*r-((a/miu)**0.5)*dsin(E)*v
    Q=(1d0-e_element**2)**(-0.5)*(dsin(E)/r_norm*r+(a/miu)**0.5*(dcos(E)-e_element)*v)
    
    R_array(1)=r(2)*v(3)-r(3)*v(2)
    R_array(2)=r(3)*v(1)-r(1)*v(3)
    R_array(3)=r(1)*v(2)-r(2)*v(1)
    R_array=R_array*1.0D0/(miu*a*(1.0D0-e_element**2))**0.5
    
    omega=datan2(P(3),Q(3))
    C_omega=datan2(R_array(1),-R_array(2))
    i=dacos(R_array(3))
    u=f+omega
    elements = (/a, e_element, i, C_omega, omega, M/)
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine nsg12cts(nsg1, PV)
    implicit none
    double precision nsg1(6)
    double precision PV(6)
    double precision kpl(6)
    call nsg12kpl(nsg1, kpl)
    call kpl2cts(kpl, PV)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine cts2nsg1(cts, nsg1)
    implicit none
    double precision, intent(in) :: cts(6)
    double precision, intent(out) :: nsg1(6)
    double precision kpl(6)
    call cts2kpl(cts, kpl)
    call kpl2nsg1(kpl, nsg1)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine nsg22cts(nsg2, PV)
    implicit none
    double precision nsg2(6), PV(6), kpl(6)
    call nsg22kpl(nsg2, kpl)
    call kpl2cts(kpl, PV)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine nsg12nsg2(nsg1, nsg2)
    implicit none
    double precision, intent(in) :: nsg1(6)
    double precision, intent(out) :: nsg2(6)
    double precision kpl(6)
    call nsg12kpl(nsg1, kpl)
    call kpl2nsg2(kpl, nsg2)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine nsg22nsg1(nsg2, nsg1)
    implicit none
    double precision, intent(in) :: nsg2(6)
    double precision, intent(out) :: nsg1(6)
    double precision kpl(6)
    call nsg22kpl(nsg2, kpl)
    call kpl2nsg1(kpl, nsg1)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! solve orbit for keplerian elements
    ! inputs:
    !   kpl     double precision(6)       keplerian
    ! outputs:
    !   r       double precision          distance of satellite from Earth center
    !   E       double precision          Eccentric Anomaly
    !   f       double precision          true anomaly
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine kepler2rEf(kpl, r, E, f)
    implicit none
    double precision, intent(in) :: kpl(6)
    double precision, intent(out) :: r, E, f
    
    !--------------------------------------
    ! Solve the kapler equation, get E
    call solve_kepler_equation(kpl(2), kpl(6), E)
    ! get r
    r = kpl(1)*(1d0 - kpl(2)*dcos(E))
    ! get f
    f = E + dasin(kpl(1)/r*kpl(2)*dsin(E)*(1d0 - 1d0/(1d0 + dsqrt(1d0 - kpl(2)**2))*kpl(2)*dcos(E)))
    
    endsubroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! solve orbit for 1st kind of non-singular elements
    ! inputs:
    !   nsg1    double precision(6)       1st kind of non-singular elements
    ! outputs:
    !   r       double precision          distance of satellite from Earth center
    !   u       double precision          \omega + f
    !   ut      double precision          \omega + E
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine nsg12rut(nsg1, r, u, ut)
    implicit none
    double precision, intent(in) :: nsg1(6)
    double precision, intent(out) :: r, u, ut
    double precision a, e, xi, eta, lambda
    !--------------------------------------
    a = nsg1(1)
    xi = nsg1(4)
    eta = nsg1(5)
    e = dsqrt(xi**2 + eta**2)
    lambda = nsg1(6)
    call SolveGeneralizedKeplerEquation(xi,eta,lambda, ut)
    r = a*(1d0 - xi*dcos(ut) - eta*dsin(ut))
    u = ut + dasin(a/r*(ut - lambda)*(1d0 - 1d0/(1d0+dsqrt(1-e**2))*(xi*dcos(ut)+eta*dsin(ut))))
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! (\xi, \eta, \lambda) -> \tilde{u}
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine SolveGeneralizedKeplerEquation(xi,eta,lambda, u_tilde)
    !----------------------------------------------------------------------------------------------------------------------------------
    !   Instruction:
    !       This subroutine solves generalized kepler equation (equivalence for non-singular elements)
    !   Date:
    !       2017-04-17
    !   Update:
    !   Inputs:
    !       xi              d               xi = e*cos(w)
    !       eta             d               eta = e*sin(w)
    !       lambda          d               lambda = M+w                                    radian
    !   Outputs:
    !       u_tilde         d               u_tilde = E+w                                   radian
    !   Note:
    !----------------------------------------------------------------------------------------------------------------------------------
    implicit none
    ! input
    double precision,intent(in):: xi,eta,lambda
    ! output
    double precision, intent(out):: u_tilde
    ! Declare local variables:
    integer:: i
    double precision:: diff,u0,u1
    double precision:: convKepEq = 1d-12
    integer::  maxIterKepEq    = 50
    
    ! Execute code:
    i = 0
    u0 = lambda
    u_tilde = u0
    diff = 1d10
    do while (dabs(diff).gt.convKepEq)
        u1 = lambda+xi*dsin(u0)-eta*dcos(u0)
        diff = u1-u0
        u0 = u1
        i = i+1
        if (i.gt.maxIterKepEq) then
            print*, "@subroutine SolveGeneralizedKeplerEquation"
            print*, "EEROR: failed to calculate u_tilde."
            return
        endif
    enddo
    u_tilde = u1
    return
    endsubroutine SolveGeneralizedKeplerEquation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! solve orbit for 2nd kind of non-singular elements
    ! inputs:
    !   nsg1    double precision(6)       2nd kind of non-singular elements
    ! outputs:
    !   r       double precision          distance of satellite from Earth center
    !   u       double precision          \Omega + \omega + f
    !   ut      double precision          \Omega + \omega + E
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine nsg22rut(nsg2, r, u, ut)
    implicit none
    double precision, intent(in) :: nsg2(6)
    double precision, intent(out) :: r, u, ut
    double precision a, e, xi, eta, lambda
    !--------------------------------------
    a = nsg2(1)
    xi = nsg2(2)
    eta = nsg2(3)
    e = dsqrt(xi**2 + eta**2)
    lambda = nsg2(6)
    call SolveGeneralizedKeplerEquation(xi,eta,lambda, ut)
    r = a*(1d0 - xi*dcos(ut) - eta*dsin(ut))
    u = ut + dasin(a/r*(ut - lambda)*(1d0 - 1d0/(1d0+dsqrt(1-e**2))*(xi*dcos(ut)+eta*dsin(ut))))
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! To solve the Kepler equation M = E - e*sin(E) as we have M and e
    ! Input:
    
    ! Output:
    
    ! Created on 2020/12/15 by Shawn Hu
    ! Modified on 2021/8/25 by Shawn Hu
    subroutine solve_kepler_equation(e, M, EE)
    implicit none
    double precision e, M, EE
    double precision deviation
    integer i
    !-----------------------------------------
    deviation = 1d0
    EE = M
    i = 0
    do while (dabs(deviation) >= 1d-13)
        deviation = M + e*dsin(EE) - EE
        EE = M + e*dsin(EE)
        i = i + 1
        if (i == 20) exit ! cycle 20 times at most
    end do
    end subroutine    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine OP2Body(x0, dt, x1)
    implicit none
    double precision x0(6), dt, x1(6)
    !---------------------------------
    x1      =   x0
    x1(6)   =   x0(6) + x0(1)**(-1.5d0)*dt
    x1(6)   =   dmod(x1(6), 3.141592653589793d0*2d0)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !============================================================================
    ! type_orbit — type-bound procedure implementations
    !============================================================================

    !----------------------------------------------------------------------------
    ! Initialize orbit from elements array
    !   elements(6)   — orbital elements
    !   eleType       — 1=Keplerian, 2=Cartesian, 3=NSG1, 4=NSG2
    !   isNormalized  — (optional) .true. if elements are already in normalized
    !                    units (default). .false. if in SI units.
    !----------------------------------------------------------------------------
    subroutine orbit_init(this, elements, eleType, isNormalized)
    use module_const, only: STATE_TYPE_KEPLERIAN, STATE_TYPE_CARTESIAN, &
                             STATE_TYPE_NONSINGULAR_1, STATE_TYPE_NONSINGULAR_2
    implicit none
    class(type_orbit), intent(inout) :: this
    double precision,  intent(in)    :: elements(6)
    integer,           intent(in)    :: eleType
    logical, optional, intent(in)    :: isNormalized
    logical :: needConversion
    double precision :: tmp(6), kpl(6)

    needConversion = .true.
    if (present(isNormalized)) needConversion = .not. isNormalized

    this%eleType = eleType

    if (needConversion) then
        select case (eleType)
        case (STATE_TYPE_KEPLERIAN)
            call kpl_SI2normal(elements, this%elements)
        case (STATE_TYPE_CARTESIAN)
            call cts_SI2normal(elements, this%elements)
        case (STATE_TYPE_NONSINGULAR_1)
            ! No direct SI→normal for NSG1; go through Keplerian
            call nsg12kpl(elements, kpl)
            call kpl_SI2normal(kpl, tmp)
            call kpl2nsg1(tmp, this%elements)
        case (STATE_TYPE_NONSINGULAR_2)
            ! No direct SI→normal for NSG2; go through Keplerian
            call nsg22kpl(elements, kpl)
            call kpl_SI2normal(kpl, tmp)
            call kpl2nsg2(tmp, this%elements)
        case default
            this%elements = elements
        end select
    else
        this%elements = elements
    end if
    end subroutine orbit_init

    !----------------------------------------------------------------------------
    ! Convert to Keplerian elements in-place
    !----------------------------------------------------------------------------
    subroutine orbit_toKeplerian(this)
    use module_const, only: STATE_TYPE_KEPLERIAN, STATE_TYPE_CARTESIAN, &
                             STATE_TYPE_NONSINGULAR_1, STATE_TYPE_NONSINGULAR_2
    implicit none
    class(type_orbit), intent(inout) :: this
    double precision :: tmp(6)

    select case (this%eleType)
    case (STATE_TYPE_KEPLERIAN)
        return
    case (STATE_TYPE_CARTESIAN)
        call cts2kpl(this%elements, tmp)
    case (STATE_TYPE_NONSINGULAR_1)
        call nsg12kpl(this%elements, tmp)
    case (STATE_TYPE_NONSINGULAR_2)
        call nsg22kpl(this%elements, tmp)
    end select
    this%elements = tmp
    this%eleType  = STATE_TYPE_KEPLERIAN
    end subroutine orbit_toKeplerian

    !----------------------------------------------------------------------------
    ! Convert to Cartesian elements in-place
    !----------------------------------------------------------------------------
    subroutine orbit_toCartesian(this)
    use module_const, only: STATE_TYPE_KEPLERIAN, STATE_TYPE_CARTESIAN, &
                             STATE_TYPE_NONSINGULAR_1, STATE_TYPE_NONSINGULAR_2
    implicit none
    class(type_orbit), intent(inout) :: this
    double precision :: tmp(6)

    select case (this%eleType)
    case (STATE_TYPE_CARTESIAN)
        return
    case (STATE_TYPE_KEPLERIAN)
        call kpl2cts(this%elements, tmp)
    case (STATE_TYPE_NONSINGULAR_1)
        call nsg12cts(this%elements, tmp)
    case (STATE_TYPE_NONSINGULAR_2)
        call nsg22cts(this%elements, tmp)
    end select
    this%elements = tmp
    this%eleType  = STATE_TYPE_CARTESIAN
    end subroutine orbit_toCartesian

    !----------------------------------------------------------------------------
    ! Convert to first kind of non-singular elements (NSG1) in-place
    !----------------------------------------------------------------------------
    subroutine orbit_toNSG1(this)
    use module_const, only: STATE_TYPE_KEPLERIAN, STATE_TYPE_CARTESIAN, &
                             STATE_TYPE_NONSINGULAR_1, STATE_TYPE_NONSINGULAR_2
    implicit none
    class(type_orbit), intent(inout) :: this
    double precision :: tmp(6)

    select case (this%eleType)
    case (STATE_TYPE_NONSINGULAR_1)
        return
    case (STATE_TYPE_KEPLERIAN)
        call kpl2nsg1(this%elements, tmp)
    case (STATE_TYPE_CARTESIAN)
        call cts2nsg1(this%elements, tmp)
    case (STATE_TYPE_NONSINGULAR_2)
        call nsg22nsg1(this%elements, tmp)
    end select
    this%elements = tmp
    this%eleType  = STATE_TYPE_NONSINGULAR_1
    end subroutine orbit_toNSG1

    !----------------------------------------------------------------------------
    ! Convert to second kind of non-singular elements (NSG2 / equinoctial) in-place
    !----------------------------------------------------------------------------
    subroutine orbit_toNSG2(this)
    use module_const, only: STATE_TYPE_KEPLERIAN, STATE_TYPE_CARTESIAN, &
                             STATE_TYPE_NONSINGULAR_1, STATE_TYPE_NONSINGULAR_2
    implicit none
    class(type_orbit), intent(inout) :: this
    double precision :: tmp(6), kpl(6)

    select case (this%eleType)
    case (STATE_TYPE_NONSINGULAR_2)
        return
    case (STATE_TYPE_KEPLERIAN)
        call kpl2nsg2(this%elements, tmp)
    case (STATE_TYPE_CARTESIAN)
        call cts2kpl(this%elements, kpl)
        call kpl2nsg2(kpl, tmp)
    case (STATE_TYPE_NONSINGULAR_1)
        call nsg12nsg2(this%elements, tmp)
    end select
    this%elements = tmp
    this%eleType  = STATE_TYPE_NONSINGULAR_2
    end subroutine orbit_toNSG2

    !----------------------------------------------------------------------------
    ! Convert elements to SI units in-place (km, km/s, degrees)
    !----------------------------------------------------------------------------
    subroutine orbit_toSI(this)
    use module_const, only: STATE_TYPE_KEPLERIAN, STATE_TYPE_CARTESIAN, &
                             STATE_TYPE_NONSINGULAR_1, STATE_TYPE_NONSINGULAR_2
    implicit none
    class(type_orbit), intent(inout) :: this
    double precision :: tmp(6), kpl(6)

    select case (this%eleType)
    case (STATE_TYPE_KEPLERIAN)
        call kpl_normal2SI(this%elements, tmp)
        this%elements = tmp
    case (STATE_TYPE_CARTESIAN)
        call cts_normal2SI(this%elements, tmp)
        this%elements = tmp
    case (STATE_TYPE_NONSINGULAR_1)
        call nsg12kpl(this%elements, kpl)
        call kpl_normal2SI(kpl, tmp)
        call kpl2nsg1(tmp, this%elements)
    case (STATE_TYPE_NONSINGULAR_2)
        call nsg22kpl(this%elements, kpl)
        call kpl_normal2SI(kpl, tmp)
        call kpl2nsg2(tmp, this%elements)
    end select
    end subroutine orbit_toSI

    !----------------------------------------------------------------------------
    ! Convert elements to normalized units in-place
    !----------------------------------------------------------------------------
    subroutine orbit_toNormalized(this)
    use module_const, only: STATE_TYPE_KEPLERIAN, STATE_TYPE_CARTESIAN, &
                             STATE_TYPE_NONSINGULAR_1, STATE_TYPE_NONSINGULAR_2
    implicit none
    class(type_orbit), intent(inout) :: this
    double precision :: tmp(6), kpl(6)

    select case (this%eleType)
    case (STATE_TYPE_KEPLERIAN)
        call kpl_SI2normal(this%elements, tmp)
        this%elements = tmp
    case (STATE_TYPE_CARTESIAN)
        call cts_SI2normal(this%elements, tmp)
        this%elements = tmp
    case (STATE_TYPE_NONSINGULAR_1)
        call nsg12kpl(this%elements, kpl)
        call kpl_SI2normal(kpl, tmp)
        call kpl2nsg1(tmp, this%elements)
    case (STATE_TYPE_NONSINGULAR_2)
        call nsg22kpl(this%elements, kpl)
        call kpl_SI2normal(kpl, tmp)
        call kpl2nsg2(tmp, this%elements)
    end select
    end subroutine orbit_toNormalized

    !----------------------------------------------------------------------------
    ! Solve Kepler equation for Keplerian elements: get r, eccentric anomaly E,
    ! and true anomaly f.
    !----------------------------------------------------------------------------
    subroutine orbit_get_rEf(this, r, E, f)
    use module_const, only: STATE_TYPE_KEPLERIAN
    implicit none
    class(type_orbit), intent(inout) :: this
    double precision,  intent(out)   :: r, E, f

    if (this%eleType /= STATE_TYPE_KEPLERIAN) then
        call this%toKeplerian()
    end if
    call kepler2rEf(this%elements, r, E, f)
    end subroutine orbit_get_rEf

    !----------------------------------------------------------------------------
    ! Solve generalized Kepler equation for NSG1/NSG2 elements: get r,
    ! argument of latitude u, and eccentric argument of latitude ut.
    !----------------------------------------------------------------------------
    subroutine orbit_get_rut(this, r, u, ut)
    use module_const, only: STATE_TYPE_NONSINGULAR_1, STATE_TYPE_NONSINGULAR_2
    implicit none
    class(type_orbit), intent(inout) :: this
    double precision,  intent(out)   :: r, u, ut

    select case (this%eleType)
    case (STATE_TYPE_NONSINGULAR_1)
        call nsg12rut(this%elements, r, u, ut)
    case (STATE_TYPE_NONSINGULAR_2)
        call nsg22rut(this%elements, r, u, ut)
    case default
        call this%toNSG1()
        call nsg12rut(this%elements, r, u, ut)
    end select
    end subroutine orbit_get_rut

    !----------------------------------------------------------------------------
    ! Two-body (Keplerian) propagation: advance mean anomaly by n*dt
    ! Only meaningful for Keplerian elements.
    !----------------------------------------------------------------------------
    subroutine orbit_propagateTwoBody(this, dt)
    use module_const, only: STATE_TYPE_KEPLERIAN
    implicit none
    class(type_orbit), intent(inout) :: this
    double precision,  intent(in)    :: dt
    double precision :: tmp(6)

    if (this%eleType /= STATE_TYPE_KEPLERIAN) then
        call this%toKeplerian()
    end if
    call OP2Body(this%elements, dt, tmp)
    this%elements = tmp
    end subroutine orbit_propagateTwoBody

    !----------------------------------------------------------------------------
    ! Display orbit information
    !----------------------------------------------------------------------------
    subroutine orbit_show(this)
    use module_const, only: STATE_TYPE_KEPLERIAN, STATE_TYPE_CARTESIAN, &
                             STATE_TYPE_NONSINGULAR_1, STATE_TYPE_NONSINGULAR_2
    implicit none
    class(type_orbit), intent(in) :: this
    character(20) :: typeName

    select case (this%eleType)
    case (STATE_TYPE_KEPLERIAN)
        typeName = "Keplerian"
    case (STATE_TYPE_CARTESIAN)
        typeName = "Cartesian"
    case (STATE_TYPE_NONSINGULAR_1)
        typeName = "NSG1"
    case (STATE_TYPE_NONSINGULAR_2)
        typeName = "NSG2"
    case default
        typeName = "Unknown"
    end select

    print '(/a)',        "========================================"
    print '(a,a)',       "  Orbit Type : ", trim(typeName)
    print '(a,6es16.8)', "  Elements   : ", this%elements
    print '(a)',         "========================================"
    end subroutine orbit_show

    endmodule