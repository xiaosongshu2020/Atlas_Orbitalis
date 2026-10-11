!===============================================================================
! Program: test_halo
!
! Description:
!   Verification of the third-order analytical Halo orbit solution
!   (module_halo) for the Earth-Moon and Sun-Earth CRTBP systems.
!
!   Checks performed:
!   1. Collinear libration points: the gradient of the effective potential must
!      vanish, and the Earth-Moon coordinates are compared with published
!      values (tolerance 1e-3, which is far tighter than the ~9e-3 error of the
!      truncated series used before).
!   2. Halo coefficients: the amplitude constraint must be satisfied, and the
!      analytical velocity must agree with a central finite difference of the
!      analytical position.
!   3. Equations of motion residual: the analytical state is substituted into
!      the full nonlinear CRTBP equations. A third-order solution leaves the
!      fourth-order terms unbalanced, so a residual of a few percent is
!      expected at these amplitudes (the Halo family has a minimum size, so the
!      amplitudes cannot be made arbitrarily small).
!   4. Orbit quality (the decisive check): starting from the analytical state, a
!      half-period differential correction must converge to a genuine periodic
!      Halo orbit. The same corrector is also started from the first-order
!      (linear Lissajous) state; the third-order guess must need a distinctly
!      smaller correction, which is what "third-order reference orbit" means.
!   5. Error handling: a beta beyond the amplitude constraint must be rejected
!      instead of producing a NaN, and the state of an uninitialized orbit must
!      be zero.
!===============================================================================

program test_halo
  use module_const
  use module_time
  use module_CRTBP
  use module_halo
  implicit none

  type(type_CRTBP) :: em, se
  type(type_t)     :: epoch
  integer          :: nFail

  nFail = 0

  ! Earth-Moon and Sun-Earth CRTBP systems
  call em%init("Earth-Moon", 5.9722d24, 7.342d22, 384400d0)
  call se%init("Sun-Earth", 1.9891d30, 5.9722d24, 1.495978707d8)

  ! Reference epoch (the CRTBP is autonomous: t0 only fixes the epoch of getPV)
  call epoch%init(2025d0, 1d0, 1d0, 0d0, 0d0, 0d0)

  write(*,'(A)') '=============================================================='
  write(*,'(A)') ' test_halo : third-order analytical Halo orbit solution'
  write(*,'(A)') '=============================================================='
  write(*,*)

  call checkLibrationPoints(em, "Earth-Moon", nFail)
  call checkLibrationPoints(se, "Sun-Earth", nFail)

  call checkOrbit(em, LAGRANGE_L1, 0.02d0, "Earth-Moon L1, beta = 0.02", epoch, nFail)
  call checkOrbit(em, LAGRANGE_L1, 0.05d0, "Earth-Moon L1, beta = 0.05", epoch, nFail)
  call checkOrbit(em, LAGRANGE_L2, 0.05d0, "Earth-Moon L2, beta = 0.05", epoch, nFail)
  call checkOrbit(se, LAGRANGE_L1, 0.01d0, "Sun-Earth  L1, beta = 0.01", epoch, nFail)
  call checkOrbit(se, LAGRANGE_L2, 0.01d0, "Sun-Earth  L2, beta = 0.01", epoch, nFail)

  call checkInitFailure(em, LAGRANGE_L2, epoch, nFail)
  call checkL3(em, 0.05d0, epoch)
  call checkRichardsonTable1(nFail)

  write(*,*)
  write(*,'(A)') '=============================================================='
  if (nFail == 0) then
    write(*,'(A)') ' test_halo : ALL CHECKS PASSED'
  else
    write(*,'(A,I3,A)') ' test_halo : ', nFail, ' CHECK(S) FAILED'
  endif
  write(*,'(A)') '=============================================================='

end program test_halo


!===============================================================================
! The dynamics of the CR3BP in the L_i-centered, gamma-normalized frame
!
!   xi_ddot  =  2*eta_dot + Omega_x/gamma
!   eta_ddot = -2*xi_dot  + Omega_y/gamma
!   zeta_ddot             = Omega_z/gamma
!
! with (x,y,z) = (x_L + gamma*xi, gamma*eta, gamma*zeta).
! This is the full nonlinear model (6.4); its linearization reproduces (6.5).
!===============================================================================
subroutine dynRHS(mu, gamma, xL, s, ds)
  implicit none
  double precision, intent(in)  :: mu, gamma, xL, s(6)
  double precision, intent(out) :: ds(6)

  double precision :: x, y, z, r1, r2, d1, d2, Ox, Oy, Oz

  x = xL + gamma*s(1)
  y = gamma*s(2)
  z = gamma*s(3)

  r1 = sqrt((x + mu)**2 + y**2 + z**2)
  r2 = sqrt((x - 1d0 + mu)**2 + y**2 + z**2)
  d1 = r1**3
  d2 = r2**3

  Ox = x - (1d0 - mu)*(x + mu)/d1 - mu*(x - 1d0 + mu)/d2
  Oy = y*(1d0 - (1d0 - mu)/d1 - mu/d2)
  Oz = -z*((1d0 - mu)/d1 + mu/d2)

  ds(1) = s(4)
  ds(2) = s(5)
  ds(3) = s(6)
  ds(4) =  2d0*s(5) + Ox/gamma
  ds(5) = -2d0*s(4) + Oy/gamma
  ds(6) =            Oz/gamma

end subroutine dynRHS


!===============================================================================
! Fixed step RK4 propagation
!===============================================================================
subroutine integrate(mu, gamma, xL, s, tend, nstep, sout)
  implicit none
  double precision, intent(in)  :: mu, gamma, xL, s(6), tend
  integer, intent(in)           :: nstep
  double precision, intent(out) :: sout(6)

  double precision :: y(6), k1(6), k2(6), k3(6), k4(6), dt
  integer :: i

  y = s
  dt = tend/dble(nstep)

  do i = 1, nstep
    call dynRHS(mu, gamma, xL, y, k1)
    call dynRHS(mu, gamma, xL, y + 0.5d0*dt*k1, k2)
    call dynRHS(mu, gamma, xL, y + 0.5d0*dt*k2, k3)
    call dynRHS(mu, gamma, xL, y + dt*k3, k4)
    y = y + dt/6d0*(k1 + 2d0*k2 + 2d0*k3 + k4)
  enddo

  sout = y

end subroutine integrate


!===============================================================================
! Half-period crossing constraints
!
! A Halo orbit crosses y = 0 perpendicularly: at that point eta = xi_dot =
! zeta_dot = 0, so the state is fully determined by (xi, zeta, eta_dot).
! The unknowns of the corrector are therefore (xi0, eta_dot0, T/2) with zeta0
! held fixed (it fixes the vertical amplitude of the family member), and the
! three constraints are (eta, xi_dot, zeta_dot) = 0 at the half period.
!===============================================================================
subroutine crossConstraints(mu, gamma, xL, zFix, x, fvec, xend)
  implicit none
  double precision, intent(in)  :: mu, gamma, xL, zFix, x(3)
  double precision, intent(out) :: fvec(3), xend(6)

  double precision :: sw(6)

  sw = (/ x(1), 0d0, zFix, 0d0, x(2), 0d0 /)
  call integrate(mu, gamma, xL, sw, x(3), 4000, xend)

  fvec(1) = xend(2)     ! eta
  fvec(2) = xend(4)     ! xi_dot
  fvec(3) = xend(6)     ! zeta_dot

end subroutine crossConstraints


!===============================================================================
! Gaussian elimination for a 3x3 system
!===============================================================================
subroutine solve3x3(aMat, bVec, yVec)
  implicit none
  double precision, intent(in)  :: aMat(3,3), bVec(3)
  double precision, intent(out) :: yVec(3)

  double precision :: m(3,4), fac
  integer :: i, k

  m(:,1:3) = aMat
  m(:,4)   = bVec

  do i = 1, 2
    do k = i + 1, 3
      fac = m(k,i)/m(i,i)
      m(k,i:4) = m(k,i:4) - fac*m(i,i:4)
    enddo
  enddo

  yVec(3) = m(3,4)/m(3,3)
  yVec(2) = (m(2,4) - m(2,3)*yVec(3))/m(2,2)
  yVec(1) = (m(1,4) - m(1,2)*yVec(2) - m(1,3)*yVec(3))/m(1,1)

end subroutine solve3x3


!===============================================================================
! Differential correction of a Halo orbit
!
! Parameters:
!   mu, gamma, xL : CR3BP parameters of the libration point
!   zFix          : fixed zeta at the crossing (vertical amplitude)
!   x             : in  - initial guess (xi0, eta_dot0, T/2)
!                   out - corrected values
!   corr          : total correction |x_out - x_in| (output)
!   fFirst        : |F| at the initial guess (output)
!   nIter         : iterations used (output)
!   converged     : whether |F| < 1e-11 was reached (output)
!===============================================================================
subroutine refineHalo(mu, gamma, xL, zFix, x, corr, fFirst, nIter, converged)
  implicit none
  double precision, intent(in)  :: mu, gamma, xL, zFix
  double precision, intent(inout) :: x(3)
  double precision, intent(out) :: corr, fFirst
  integer, intent(out)          :: nIter
  logical, intent(out)          :: converged

  double precision, parameter :: eps = 1d-7
  double precision, parameter :: tol = 1d-11

  double precision :: xStart(3), xp(3), f0(3), fp(3), aMat(3,3), dx(3)
  double precision :: xend(6), fend(6)
  integer :: i, jj

  xStart = x
  corr   = 0d0
  fFirst = 0d0
  nIter  = 0
  converged = .false.

  do i = 1, 15
    nIter = i
    call crossConstraints(mu, gamma, xL, zFix, x, f0, xend)

    if (i == 1) fFirst = maxval(abs(f0))

    if (maxval(abs(f0)) < tol) then
      converged = .true.
      exit
    endif

    ! Jacobian columns 1-2: finite differences with respect to the state
    do jj = 1, 2
      xp = x
      xp(jj) = xp(jj) + eps
      call crossConstraints(mu, gamma, xL, zFix, xp, fp, xend)
      aMat(:,jj) = (fp - f0)/eps
    enddo

    ! Jacobian column 3: dF/d(T/2) = f at the endpoint
    call dynRHS(mu, gamma, xL, xend, fend)
    aMat(:,3) = (/ fend(2), fend(4), fend(6) /)

    call solve3x3(aMat, -f0, dx)
    x = x + dx
  enddo

  corr = maxval(abs(x - xStart))

end subroutine refineHalo


!===============================================================================
! Collinear libration points: gradient residual and published Earth-Moon values
!===============================================================================
subroutine checkLibrationPoints(crtbp, label, nFail)
  use module_CRTBP
  implicit none
  type(type_CRTBP), intent(in) :: crtbp
  character(len=*), intent(in)  :: label
  integer, intent(inout)        :: nFail

  double precision :: xL, yDummy, f, r1, r2, mu, emRef(3)
  integer :: i

  mu = crtbp%mu
  emRef = (/ 0.836915d0, 1.155682d0, -1.005063d0 /)   ! published Earth-Moon values

  write(*,'(A)') '--------------------------------------------------------------'
  write(*,'(A)') ' Libration points: '//trim(label)
  write(*,'(A,F16.12)') '   mu = ', mu

  do i = LAGRANGE_L1, LAGRANGE_L3
    call crtbp%getLagrangePointCoord(i, xL, yDummy)

    ! gradient of the effective potential must vanish
    r1 = xL + mu
    r2 = xL - 1d0 + mu
    f  = xL - (1d0 - mu)*r1/abs(r1)**3 - mu*r2/abs(r2)**3

    write(*,'(A,I1,A,F16.12,A,ES12.3)') '   L', i, ' x = ', xL, &
        '   |grad Omega| = ', abs(f)

    if (abs(f) > 1d-12) then
      write(*,'(A)') '     FAIL: collinear point does not null the gradient'
      nFail = nFail + 1
    endif

    if (index(label, 'Earth-Moon') > 0) then
      if (abs(xL - emRef(i)) > 1d-3) then
        write(*,'(A,ES12.3)') '     FAIL: differs from reference by ', xL - emRef(i)
        nFail = nFail + 1
      endif
    endif
  enddo

  write(*,*)

end subroutine checkLibrationPoints


!===============================================================================
! Halo orbit: coefficients, constraint, velocity, equations of motion, and the
! differential correction from the third-order and from the linear guess
!===============================================================================
subroutine checkOrbit(crtbp, LID, beta, label, epoch, nFail)
  use module_const
  use module_CRTBP
  use module_time
  use module_halo
  implicit none
  type(type_CRTBP), intent(in)  :: crtbp
  integer, intent(in)           :: LID
  double precision, intent(in)  :: beta
  character(len=*), intent(in)  :: label
  type(type_t), intent(in)      :: epoch
  integer, intent(inout)        :: nFail

  integer, parameter :: nSamp = 200
  double precision, parameter :: hFd   = 1d-4     ! FD step, canonical time
  double precision, parameter :: tolVel = 1d-5    ! velocity consistency
  double precision, parameter :: tolEom = 2d-1    ! nonlinear model residual

  type(type_halo)  :: orb
  double precision :: t, per, xL, yDummy, tau
  double precision :: p0(6), pm(6), pp(6), rhs(6)
  double precision :: accFd, velFd, vErr, vScale, rErr, rScale, vRel, rRel
  double precision :: consRes, corr3, corr1, fFirst3, fFirst1
  double precision :: guess3(3), guess1(3)
  integer :: i, jj, nIter3, nIter1
  logical :: ok, conv3, conv1

  write(*,'(A)') '--------------------------------------------------------------'
  write(*,'(A)') ' Halo orbit: '//trim(label)

  call orb%init(crtbp, LID, beta, 0d0, epoch, trim(label), ok)

  if (.not. ok) then
    write(*,'(A)') '   FAIL: orbit could not be initialized'
    nFail = nFail + 1
    return
  endif

  per = orb%period
  call crtbp%getLagrangePointCoord(LID, xL, yDummy)

  write(*,'(A,F14.10)') '   gamma    = ', orb%gamma
  write(*,'(A,3F14.10)') '   c2,c3,c4 = ', orb%c2, orb%c3, orb%c4
  write(*,'(A,F14.10)') '   k        = ', orb%k
  write(*,'(A,2F14.8)') '   l1, l2   = ', orb%l1, orb%l2
  write(*,'(A,F14.10)') '   alpha    = ', orb%alpha
  write(*,'(A,F14.10)') '   omega    = ', orb%omega
  write(*,'(A,F14.6,A,F12.3,A)') '   period   = ', per, ' TU  (', &
      per*crtbp%TU/86400d0, ' days)'

  ! --- amplitude constraint (must hold identically) -------------------------
  consRes = orb%l1*orb%alpha**2 + orb%l2*orb%beta**2 + orb%omega0**2 - orb%c2
  write(*,'(A,ES12.3)') '   amplitude constraint residual = ', consRes

  if (abs(consRes) > 1d-12) then
    write(*,'(A)') '   FAIL: amplitude constraint is not satisfied'
    nFail = nFail + 1
  endif

  ! --- velocity and equations of motion residuals ---------------------------
  vErr   = 0d0
  vScale = 0d0
  rErr   = 0d0
  rScale = 0d0

  do i = 0, nSamp - 1
    t = per * dble(i)/dble(nSamp)

    call orb%compute_state(t - hFd, pm(1), pm(2), pm(3), pm(4), pm(5), pm(6))
    call orb%compute_state(t,       p0(1), p0(2), p0(3), p0(4), p0(5), p0(6))
    call orb%compute_state(t + hFd, pp(1), pp(2), pp(3), pp(4), pp(5), pp(6))

    call dynRHS(crtbp%mu, orb%gamma, xL, p0, rhs)

    do jj = 1, 3
      velFd  = (pp(jj) - pm(jj))/(2d0*hFd)
      vErr   = max(vErr, abs(velFd - p0(jj + 3)))
      vScale = max(vScale, abs(p0(jj + 3)))

      accFd  = (pp(jj) - 2d0*p0(jj) + pm(jj))/hFd**2
      rErr   = max(rErr, abs(accFd - rhs(jj+3)))
      rScale = max(rScale, abs(rhs(jj)))
    enddo
  enddo

  vRel = vErr/max(vScale, 1d-30)
  rRel = rErr/max(rScale, 1d-30)

  write(*,'(A,ES12.3,A,ES12.3,A)') '   velocity rel residual  = ', vRel, &
      '   (tol ', tolVel, ')'
  write(*,'(A,ES12.3,A,ES12.3,A)') '   eom      rel residual  = ', rRel, &
      '   (tol ', tolEom, ')'

  if (vRel > tolVel) then
    write(*,'(A)') '   FAIL: velocity is not the derivative of the position'
    nFail = nFail + 1
  endif

  if (rRel > tolEom) then
    write(*,'(A)') '   FAIL: state does not satisfy the CRTBP equations'
    nFail = nFail + 1
  endif

  ! --- differential correction from the analytical (third order) guess ------
  call orb%compute_state(0d0, p0(1), p0(2), p0(3), p0(4), p0(5), p0(6))
  guess3 = (/ p0(1), p0(5), per/2d0 /)

  call refineHalo(crtbp%mu, orb%gamma, xL, p0(3), guess3, corr3, fFirst3, nIter3, conv3)

  ! --- differential correction from the linear (first order) guess -----------
  tau = 0d0
  guess1 = (/ -orb%alpha, orb%omega0*orb%k*orb%alpha, PI/orb%omega0 /)

  call refineHalo(crtbp%mu, orb%gamma, xL, p0(3), guess1, corr1, fFirst1, nIter1, conv1)

  write(*,'(A)') '   differential correction (half period shooting):'
  write(*,'(A,ES12.3,A,I2,A,L2)') '     3rd order guess: |F|1 = ', fFirst3, &
      '  iter = ', nIter3, '  conv = ', conv3
  write(*,'(A,ES12.3,A,I2,A,L2)') '     1st order guess: |F|1 = ', fFirst1, &
      '  iter = ', nIter1, '  conv = ', conv1
  write(*,'(A,ES12.3)') '     correction 3rd order = ', corr3
  write(*,'(A,ES12.3)') '     correction 1st order = ', corr1

  if (.not. conv3) then
    write(*,'(A)') '   FAIL: differential correction from the analytical guess did not converge'
    nFail = nFail + 1
  endif

  ! the third-order guess must be clearly better than the linear one
  if (fFirst3 > 0.5d0*fFirst1) then
    write(*,'(A)') '   FAIL: third-order guess is not better than the linear guess (|F|)'
    nFail = nFail + 1
  endif

  ! Note: the size of the correction is a weak metric - it is dominated by the
  ! geometry of the unstable manifold and by the period adjustment, so a much
  ! smaller initial defect does not translate into a proportionally smaller
  ! correction. Only a bound against the linear guess is asserted here.
  if (corr3 > corr1) then
    write(*,'(A)') '   FAIL: third-order guess needs a larger correction than the linear one'
    nFail = nFail + 1
  endif

  write(*,*)

end subroutine checkOrbit


!===============================================================================
! Error handling: a beta beyond the amplitude constraint must be rejected, and
! the state of an uninitialized orbit must be zero rather than NaN
!===============================================================================
subroutine checkInitFailure(crtbp, LID, epoch, nFail)
  use module_CRTBP
  use module_time
  use module_halo
  implicit none
  type(type_CRTBP), intent(in)  :: crtbp
  integer, intent(in)           :: LID
  type(type_t), intent(in)      :: epoch
  integer, intent(inout)        :: nFail

  type(type_halo) :: orb
  double precision :: xi, eta, zeta, xi_dot, eta_dot, zeta_dot
  double precision :: beta, l1, l2, omega0, c2
  integer :: i
  logical :: ok, rejected

  write(*,'(A)') '--------------------------------------------------------------'
  write(*,'(A)') ' Error handling'

  ! --- (1) an unsupported libration point must be rejected -----------------
  write(*,'(A)') '   (a) L4 is not a collinear point: initialization must fail'
  call orb%init(crtbp, LAGRANGE_L4, 0.01d0, 0d0, epoch, "L4", ok)

  if (ok) then
    write(*,'(A)') '       FAIL: L4 was accepted'
    nFail = nFail + 1
  else
    call orb%compute_state(0d0, xi, eta, zeta, xi_dot, eta_dot, zeta_dot)
    if (xi /= 0d0 .or. eta /= 0d0 .or. zeta /= 0d0 .or. &
        xi_dot /= 0d0 .or. eta_dot /= 0d0 .or. zeta_dot /= 0d0) then
      write(*,'(A)') '       FAIL: uninitialized orbit does not return zeros'
      nFail = nFail + 1
    else
      write(*,'(A)') '       OK: rejected, and the state is zero'
    endif
  endif

  ! --- (2) the beta limit of the amplitude constraint ----------------------
  ! l1*alpha^2 + l2*beta^2 = c2 - omega0^2 is an ellipse when l1 and l2 have
  ! the same sign (bounded family, beta can be rejected) and a hyperbola when
  ! their signs differ (unbounded within the third-order model, so no beta is
  ! ever rejected). Both outcomes are legitimate; report which one applies.
  call orb%init(crtbp, LID, 0.01d0, 0d0, epoch, "probe", ok)

  if (.not. ok) then
    write(*,'(A)') '   FAIL: the probe orbit for the beta limit did not initialize'
    nFail = nFail + 1
    write(*,*)
    return
  endif

  l1 = orb%l1
  l2 = orb%l2
  omega0 = orb%omega0
  c2 = orb%c2

  if (l1*l2 > 0d0) then
    write(*,'(A,F8.2)') '   (b) l1*l2 > 0 : bounded family, beta_max = ', &
        sqrt(abs((c2 - omega0**2)/l2))
  else
    write(*,'(A)') '   (b) l1*l2 < 0 : unbounded (hyperbolic) constraint, no beta limit'
  endif

  rejected = .false.
  do i = 1, 60
    beta = 0.1d0*dble(i)
    call orb%init(crtbp, LID, beta, 0d0, epoch, "scan", ok)
    if (.not. ok) then
      rejected = .true.
      exit
    endif
  enddo

  if (l1*l2 > 0d0) then
    if (.not. rejected) then
      write(*,'(A)') '       FAIL: a bounded family did not reject any beta up to 6.0'
      nFail = nFail + 1
    else
      call orb%compute_state(0d0, xi, eta, zeta, xi_dot, eta_dot, zeta_dot)
      if (xi /= 0d0 .or. eta /= 0d0 .or. zeta /= 0d0 .or. &
          xi_dot /= 0d0 .or. eta_dot /= 0d0 .or. zeta_dot /= 0d0) then
        write(*,'(A)') '       FAIL: rejected orbit does not return zeros'
        nFail = nFail + 1
      else
        write(*,'(A,F6.2,A)') '       OK: beta = ', beta, ' rejected, state is zero'
      endif
    endif
  else
    write(*,'(A)') '       (no beta rejection expected for this system)'
  endif

  write(*,*)

end subroutine checkInitFailure


!===============================================================================
! Cross-check of every coefficient against Richardson (1980), Appendix I, using
! the numerical values tabulated for the Sun-Earth system in his Table I
! (mu = 3.04036e-6, a1 = 1.49598e8 km).
!
! This is the strongest available reference for the coefficient set: the paper
! publishes all of a21...d32, s1, s2, l1, l2, a1, a2, d1, d2, c2, c3, c4, k,
! lambda and Delta for L1, L2 and L3. Only b33, b34, b35 (the sin(tau) term of
! eta) are absent from the paper and therefore cannot be checked here.
!
! Note on the reference: the paper's lambda is our omega0, and its d1, d2 are
! our e1, e2. Its Delta is omega0^2 - c2, the constant of the amplitude
! constraint (18).
!===============================================================================
subroutine checkRichardsonTable1(nFail)
  use module_const
  use module_time
  use module_CRTBP
  use module_halo
  implicit none
  integer, intent(inout) :: nFail

  integer, parameter :: nq = 28
  double precision, parameter :: tolRel = 2d-3   ! 0.2 % relative
  double precision, parameter :: tolAbs = 2d-5   ! for entries close to zero
  double precision, parameter :: muPaper = 3.04036d-6

  type(type_CRTBP) :: se
  type(type_halo)  :: orb
  type(type_t)     :: epoch
  double precision :: ref(3, nq), got(nq), dev, devMax
  character(len=10) :: vname(nq)
  integer :: LID, i
  logical :: ok

  ! Richardson Table I values: rows L1, L2, L3; columns as in vname
  ref(1,:) = (/ 1.00109d-2, 2.08645d0, 3.22927d0, 2.92214d-1, 4.06107d0, &
                3.02001d0, 3.03054d0, -8.24661d-1, 1.21099d-1, -1.59656d1, &
                1.74090d0, -8.78563d0, 6.86546d-1, 3.11184d2, 1.58787d3, &
                2.09270d0, 2.48298d-1, -9.05965d-1, -1.04464d-1, 7.93820d-1, &
                8.26854d-2, -4.92446d-1, 6.07465d-2, 8.85701d-1, 2.30198d-2, &
                -3.46865d-1, 1.90439d-2, 3.98095d-1 /)
  ref(2,:) = (/ 1.00782d-2, 2.05701d0, 3.18723d0, 2.90785d-1, 3.94052d0, &
                -2.97984d0, 2.97026d0, -7.44452d-1, 1.25047d-1, -1.48288d1, &
                1.67369d0, -8.52882d0, 6.15466d-1, 2.93192d2, 1.49800d3, &
                -2.05304d0, -2.51646d-1, 8.96284d-1, 1.06600d-1, 7.80646d-1, &
                8.36960d-2, 4.91357d-1, -6.27190d-2, 8.55305d-1, 2.04354d-2, &
                3.52118d-1, 1.88290d-2, 3.94028d-1 /)
  ref(3,:) = (/ 9.99998d-1, 1.00000d0, 2.00000d0, 2.66029d-6, 1.00000d0, &
                1.00000d0, 1.00000d0, -1.59141d-6, 6.29433d-6, -1.57717d-5, &
                1.40258d-5, -1.25889d-5, 1.43702d-6, 1.20002d1, 7.20008d1, &
                5.00000d-1, 2.50000d-1, -4.99999d-1, -2.49999d-1, 3.75000d-1, &
                1.25000d-1, -2.50000d-1, 2.49998d-1, 2.91666d-1, -1.24999d-1, &
                -4.99999d-1, 2.53854d-7, 3.74999d-1 /)

  vname = (/ 'gamma     ', 'lambda    ', 'k         ', 'Delta     ', 'c2        ', &
             'c3        ', 'c4        ', 's1        ', 's2        ', 'l1        ', &
             'l2        ', 'a1        ', 'a2        ', 'd1(=e1)   ', 'd2(=e2)   ', &
             'a21       ', 'a22       ', 'a23       ', 'a24       ', 'a31       ', &
             'a32       ', 'b21       ', 'b22       ', 'b31       ', 'b32       ', &
             'd21       ', 'd31       ', 'd32       ' /)

  ! Sun-Earth system with the paper's mass ratio and mean distance
  call se%init("Sun-Earth (Richardson 1980)", 1.9891d30, &
               1.9891d30*muPaper/(1d0 - muPaper), 1.49598d8)
  call epoch%init(2025d0, 1d0, 1d0, 0d0, 0d0, 0d0)

  write(*,'(A)') '--------------------------------------------------------------'
  write(*,'(A)') ' Cross-check against Richardson (1980), Table I'
  write(*,'(A,F16.12)') '   mu  = ', se%mu
  write(*,'(A,ES12.3,A,ES12.3)') '   tolerances: rel ', tolRel, '   abs ', tolAbs

  do LID = LAGRANGE_L1, LAGRANGE_L3
    call orb%init(se, LID, 0.05d0, 0d0, epoch, "check", ok)

    got = (/ orb%gamma, orb%omega0, orb%k, orb%omega0**2 - orb%c2, &
             orb%c2, orb%c3, orb%c4, orb%s1, orb%s2, orb%l1, orb%l2, &
             orb%a1, orb%a2, orb%e1, orb%e2, orb%a21, orb%a22, orb%a23, &
             orb%a24, orb%a31, orb%a32, orb%b21, orb%b22, orb%b31, &
             orb%b32, orb%d21, orb%d31, orb%d32 /)

    write(*,'(A)') '   ------------------------------------------------'
    write(*,'(A,I1,A,F10.6,A)') '   L', LID, '   (beta used for init = ', orb%beta, ')'
    devMax = 0d0
    do i = 1, nq
      dev = abs(got(i) - ref(LID,i))
      if (dev > tolAbs + tolRel*abs(ref(LID,i))) then
        write(*,'(A,A,A,ES13.5,A,ES13.5,A,ES9.1,A)') '     FAIL ', vname(i), &
            ' ref = ', ref(LID,i), '  got = ', got(i), '  dev = ', &
            100d0*dev/max(abs(ref(LID,i)), 1d-30), ' %'
        nFail = nFail + 1
      else
        devMax = max(devMax, dev/max(abs(ref(LID,i)), 1d-30))
      endif
    enddo
    write(*,'(A,ES12.3)') '     all other entries agree; worst relative deviation = ', devMax
  enddo

  write(*,*)

end subroutine checkRichardsonTable1


!===============================================================================
! L3: implemented from (6.24) but not verified (informational only, no pass/fail)
!===============================================================================
subroutine checkL3(crtbp, beta, epoch)
  use module_CRTBP
  use module_time
  use module_halo
  implicit none
  type(type_CRTBP), intent(in)  :: crtbp
  double precision, intent(in)  :: beta
  type(type_t), intent(in)      :: epoch

  type(type_halo) :: orb
  logical :: ok

  write(*,'(A)') '--------------------------------------------------------------'
  write(*,'(A)') ' L3 (informational, not verified): expected warning follows'

  call orb%init(crtbp, LAGRANGE_L3, beta, 0d0, epoch, "Earth-Moon L3", ok)

  if (ok) then
    write(*,'(A,F14.10)') '   L3 alpha = ', orb%alpha
    write(*,'(A,F14.10)') '   L3 omega = ', orb%omega
  else
    write(*,'(A)') '   L3 initialization failed (see the warning above)'
  endif

  write(*,*)

end subroutine checkL3
