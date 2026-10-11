    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    module module_space
    use module_const
    use module_toolKits
    use module_matrix
    use module_planets
    use module_time
    use EOPParam

    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Transfer station position from (longitude, latitude, altitude)
    ! to Cartesian in Earth Center Fixed (ECF).
    ! Author    :   Shawn Hu
    ! Date      :   2022/10/26
    ! Input     :
    !   longitude       double precision      rad     
    !   latitude        double precision      rad
    !   altitude        double precision      Ae
    ! Output    :
    !   position        double precision(3)   Ae      Station position in ECF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  station_lla2xyz(longitude, latitude, altitude, position)
    implicit    none
    double precision      longitude, latitude, altitude
    double precision      position(3), NN, epsi
    !-----------------------------------------------------------------
    epsi        =   1d0/298.257223563d0
    NN          =   ((dcos(latitude))**2d0 + &
    &               (1d0 - epsi)**2d0*(dsin(latitude))**2d0)**(-0.5d0)
    position(1) =   (NN + altitude)*dcos(latitude)*dcos(longitude)
    position(2) =   (NN + altitude)*dcos(latitude)*dsin(longitude)
    position(3) =   (NN*(1 - epsi)**2 + altitude)*dsin(latitude)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  ECI2ECF_IAU1980_position(utc, rECI, rECF)
    
    implicit    none
    double precision      utc, HG(3, 3), rECI(3), rECF(3)
    
    !-------------------------------------------------------------------
    call        HG_get(utc, HG)
    rECF    =   matmul(HG,  rECI)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  ECI2ECF_IAU1980_positionVelocity(utc, rECI, vECI, rECF, vECF)
    implicit    none
    double precision      utc, rECI(3), vECI(3), rECF(3), vECF(3)
    double precision      HG(3, 3), HG_dot(3, 3), HG1(3, 3), utc1
    double precision      rTmp(3), vTmp(3)
    !-------------------------------------------------------------------
    call        HG_get(utc, HG)
    rTmp    =   matmul(HG,  rECI)
    utc1    =   utc + 1d-1/86400d0
    call        HG_get(utc1, HG1)
    HG_dot  =   (HG1 - HG)/(1d-1/time_unit)
    vTmp    =   matmul(HG, vECI) + matmul(HG_dot, rECI)
    rECF    =   rTmp
    vECF    =   vTmp
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  ECF2ECI_IAU1980_position(utc, rECF, rECI)
    implicit    none
    double precision      utc, HG(3, 3), rECI(3), rECF(3)
    !-------------------------------------------------------------------
    call        HG_get(utc, HG)
    rECI    =   matmul(transpose(HG),  rECF)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  ECF2ECI_IAU2006A_position(utc, rECF, rECI)
    implicit    none
    double precision      utc, HG(3, 3), rECI(3), rECF(3)
    !-------------------------------------------------------------------
    call        HG_get_sofa_06a(utc, HG)
    rECI    =   matmul(transpose(HG),  rECF)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  ECF2ECI_IAU1980_positionVelocity(utc, rECF, vECF, rECI, vECI)
    use         module_const
    implicit    none
    double precision      utc, rECF(3), vECF(3), rECI(3), vECI(3)
    double precision      HG(3, 3), HGT(3, 3), HGT_dot(3, 3), HGT1(3, 3), utc1
    !-------------------------------------------------------------------
    call        HG_get(utc, HG)
    HGT     =   transpose(HG)
    rECI    =   matmul(HGT,  rECF)
    utc1    =   utc + 1d-1/86400d0
    call        HG_get(utc1, HG)
    HGT1    =   transpose(HG)
    HGT_dot =   (HGT1 - HGT)/(1d-1/time_unit)
    vECI    =   matmul(HGT, vECF) + matmul(HGT_dot, rECF)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   rHorizon    =   ZR*rECF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  ECF2horizon_position(longi, lati, rECF, rHorizon)
    implicit    none
    double precision      longi, lati, rECF(3), rHorizon(3)
    double precision      ZR(3, 3)
    call        getZR(longi, lati, ZR)
    rHorizon    =   matmul(ZR, rECF)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   rECF    =   ZR'*rHorizon
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  horizon2ECF_position(longi, lati, rHorizon, rECF)
    implicit    none
    double precision      longi, lati, rHorizon(3), rECF(3)
    double precision      ZR(3, 3)
    call        getZR(longi, lati, ZR)
    rECF    =   matmul(transpose(ZR), rHorizon)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ECI2EarthMoonsynodic(r1, v1, rm, vm, rc, vc)
    implicit none
    double precision  r1(3), v1(3), rm(3), vm(3), rc(3), vc(3)
    double precision  mass_lunar, ro(3), vo(3), R(3, 3), Rd(3, 3)
    double precision  vmT(3), n, am(3)
    
    mass_lunar = 1.2300037226627106D-2
    ! mass center postion and velocity in ECI
    ro = mass_lunar*rm/(1d0 + mass_lunar) 
    vo = mass_lunar*vm/(1d0 + mass_lunar)
    ! rotation matrix R
    R(1, 1:3) = rm/norm2(rm)
    call cross_3d(rm, vm, am)
    R(3, 1:3) = am/norm2(am)
    call cross_3d(R(3, 1:3), R(1, 1:3), R(2, 1:3))
    ! rc = R*(r1 - ro)
    rc = matmul(R, r1 - ro)
    
    ! n is angular velocity of the moon or Earth
    call cross_3d(vm, rm/norm2(rm), vmT)
    n = norm2(vmT)/norm2(rm)
    ! Rd = dR/dt
    Rd(1, 1:3) = n*R(2, 1:3)
    Rd(2, 1:3) = -n*R(1, 1:3)
    Rd(3, 1:3) = 0d0
    ! vc = R*(v1 - vo) + Rd*(r1 - ro)
    vc = matmul(R, v1 - vo) + matmul(Rd, r1 - ro)
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine EarthMoonSynodic2ECI(r1, v1, rm, vm, rc, vc)
    implicit none
    double precision  r1(3), v1(3), rm(3), vm(3), rc(3), vc(3)
    double precision  mass_lunar, ro(3), vo(3), R(3, 3), Rd(3, 3)
    double precision  vmT(3), n, am(3)
    double precision  Rr(3, 3)
    
    mass_lunar = 1.2300037226627106D-2
    
    ! mass center postion and velocity in ECI
    ro = mass_lunar*rm/(1d0 + mass_lunar) 
    vo = mass_lunar*vm/(1d0 + mass_lunar)
    
    ! rotation matrix R
    R(1, 1:3) = rm/norm2(rm)
    call cross_3d(rm, vm, am)
    R(3, 1:3) = am/norm2(am)
    call cross_3d(R(3, 1:3), R(1, 1:3), R(2, 1:3))
    call InvMatrix_LUsolver(3, R, Rr)
    ! call matrix_inverse_2(3, R, Rr)
    
    ! rc = ro + Rr*r2
    rc = ro + matmul(Rr, r1)
    
    ! n is angular velocity of the moon or Earth
    call cross_3d(vm, rm/norm2(rm), vmT)
    n = norm2(vmT)/norm2(rm)
    
    ! Rd = dR/dt
    Rd(1, 1:3) = n*R(2, 1:3)
    Rd(2, 1:3) = -n*R(1, 1:3)
    Rd(3, 1:3) = 0d0
    
    ! vc = vo + Rr*(v1 - Rd*(rc - ro))
    vc = vo + matmul(Rr, v1 - matmul(Rd, rc - ro))
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !This subroutine is written to get the xyz coordinate when
    !the longitude, latitude and height have been known
    !unit of length is km
    subroutine get_station_position(lambda, fi, h, RO_fixed)
    implicit none
    double precision lambda, fi, h, RO_fixed(3)
    double precision N, pi_const, ae, epsi
    !-----------------------------------------------------
    epsi = 1d0/298.257223563d0
    ae = 6378.13630d0
    !-----------------------------------------------------
    N = ae*((dcos(fi))**2d0 + (1d0 - epsi)**2d0*(dsin(fi))**2d0)**(-0.5d0)
    RO_fixed(1) = (N + h)*dcos(fi)*dcos(lambda)/ae
    RO_fixed(2) = (N + h)*dcos(fi)*dsin(lambda)/ae
    RO_fixed(3) = (N*(1 - epsi)**2 + h)*dsin(fi)/ae
    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! convert position from moon-center-fixed to moon-center-inertial
    ! t : TT in JD
    ! rf : position in moon-center-fixed
    ! ri : posiiton in moon-center-inertial
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine LCF2LCI(t, rf, ri)
    implicit none
    double precision, intent(in) :: t, rf(3)
    double precision, intent(out) :: ri(3)
    double precision :: OmegaP, is, Lambda, M1(3, 3), M2(3, 3), M3(3, 3)
    double precision :: M(3, 3), PV(6)
    integer :: iPlephStatus
    !-----------------------------------------------------
    call pleph(t, 15, 0, PV)
    ! pleph no longer aborts on a bad ephemeris: it returns zeros, so check
    ! the status -- a silent zero rotation here would be an invisible error.
    call pleph_status(iPlephStatus)
    if (iPlephStatus /= 0) then
        write(0, '(2a,i0,a)') "Error @LCF2LCI: JPL librations unavailable ", &
            &   "(pleph status ", iPlephStatus, &
            &   "); the Moon-fixed frame cannot be built."
        write(0, *) "  Run from the repository root, or restore the JPL " // &
            &   "ephemeris file.  No analytical fallback exists for librations."
        flush(0)
        error stop
    endif
    OmegaP = PV(1)
    is = PV(2)
    Lambda = PV(3)
    Lambda = dmod(Lambda, 2d0*pi)
    call RzMatrix(OmegaP, M1)
    call RxMatrix(is, M2)
    call RzMatrix(Lambda, M3)
    M = matmul(M3, matmul(M2, M1))
    ri = matmul(transpose(M), rf)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! convert position from moon-center-inertal to moon-center-fixed
    ! t : TT in JD
    ! ri : posiiton in moon-center-inertial
    ! rf : position in moon-center-fixed
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine LCI2LCF(t, ri, rf)
    implicit none
    double precision, intent(in) :: t, ri(3)
    double precision, intent(out) :: rf(3)
    double precision :: OmegaP, is, Lambda, M1(3, 3), M2(3, 3), M3(3, 3)
    double precision :: M(3, 3), PV(6)
    integer :: iPlephStatus
    !-----------------------------------------------------
    call pleph(t, 15, 0, PV)
    ! see the note in LCF2LCI: pleph returns zeros instead of aborting now
    call pleph_status(iPlephStatus)
    if (iPlephStatus /= 0) then
        write(0, '(2a,i0,a)') "Error @LCI2LCF: JPL librations unavailable ", &
            &   "(pleph status ", iPlephStatus, &
            &   "); the Moon-fixed frame cannot be built."
        write(0, *) "  Run from the repository root, or restore the JPL " // &
            &   "ephemeris file.  No analytical fallback exists for librations."
        flush(0)
        error stop
    endif
    OmegaP = PV(1)
    is = PV(2)
    Lambda = PV(3)
    Lambda = dmod(Lambda, 2d0*pi)
    call RzMatrix(OmegaP, M1)
    call RxMatrix(is, M2)
    call RzMatrix(Lambda, M3)
    M = matmul(M3, matmul(M2, M1))
    rf = matmul(M, ri)
    endsubroutine
    
    subroutine NTH2ECI(r, v, M)
    implicit none
    double precision, intent(in) :: r(3), v(3)
    double precision, intent(out) :: M(3, 3)
    double precision :: h(3), vxh(3), i(3), j(3), k(3)
    double precision :: norm_h, norm_v, norm_vxh
    
    ! Calculate angular momentum vector h = r × v
    call cross_3d(r, v, h)
    
    ! Calculate v × h
    call cross_3d(v, h, vxh)
    
    ! Calculate norms
    norm_h = norm2(h)
    norm_v = norm2(v)
    norm_vxh = norm2(vxh)
    
    ! Normal vectors for NTH coordinate system
    ! i = (v × h) / |v × h|  (normal vector)
    if (norm_vxh > 0d0) then
        i = vxh / norm_vxh
    else
        i = (/1d0, 0d0, 0d0/)
    endif
    
    ! j = v / |v|  (tangential vector)
    if (norm_v > 0d0) then
        j = v / norm_v
    else
        j = (/0d0, 1d0, 0d0/)
    endif
    
    ! k = h / |h|  (horizontal vector)
    if (norm_h > 0d0) then
        k = h / norm_h
    else
        k = (/0d0, 0d0, 1d0/)
    endif
    
    ! Construct transformation matrix M
    ! M transforms from NTH coordinates to ECI coordinates
    ! Each column of M contains the NTH basis vectors in ECI coordinates
    M(1:3, 1) = i  ! normal vector
    M(1:3, 2) = j  ! tangential vector
    M(1:3, 3) = k  ! horizontal vector

    endsubroutine

    subroutine ECI2NTH(r, v, M)
    implicit none
    double precision, intent(in) :: r(3), v(3)
    double precision, intent(out) :: M(3, 3)
    double precision :: MT(3, 3)
    call NTH2ECI(r, v, MT)
    M = transpose(MT)
    endsubroutine
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ZR        =   Ry*Rz(90)*Rx(-(90 - phi))*Rz(lambda - 90)
    ! r_gound   =   ZR*r_ECF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  getZR(longi, lati, ZR)
    implicit    none
    double precision      longi, lati, ZR(3, 3)
    double precision      R1(3, 3), R2(3, 3), R3(3, 3), R4(3, 3)
    double precision      ::  pi_const    =   3.141592653589793d0
    !-------------------------------------------------------------------
    
    call        coordinate_reverse('y', R1)
    !call        Rzmatrix(pi_const/2d0, R2)                              ! Horizonal X-axis points the South
    call        Rzmatrix(pi_const*3d0/2d0, R2)                          ! Horizonal X-axis points the North
    call        Rxmatrix(-(pi_const/2d0 - lati), R3)
    call        Rzmatrix(longi - pi_const/2d0, R4)

    ZR      =   matmul(R1, matmul(R2, matmul(R3, R4)))

    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine HG_get_sofa_06a(t, HG)
    implicit none
    double precision, intent(in) :: t
    double precision, intent(out) :: HG(3, 3)
    double precision :: UTC, UT1, TT

    ! The dummy names and the ABSENCE of INTENT must match sofa/c2t06a.for
    ! exactly.  Intel's /warn:interfaces (on by default in VS Debug builds)
    ! compares this block against the interface it generates from the SOFA
    ! source and fails with error #8000 otherwise.
    ! Argument list: 1-6 are intent(in), 7 is intent(out).
    interface
    subroutine iau_C2T06A(TTA, TTB, UTA, UTB, XP, YP, RC2T)
        double precision :: TTA, TTB, UTA, UTB, XP, YP
        double precision :: RC2T(3, 3)
    end subroutine iau_C2T06A
    end interface
    
    !Synchronize EOP coefficients
    call syncEOP(t)
    !Get UT1 and UTC
    UTC = t
    call UTC2UT1(UTC, DUT1, UT1)
    call UTC2TT(UTC, TT)
    
    call iau_C2T06A(2400000.5d0, TT, 2400000.5d0, UT1, Xp, Yp, HG)
    
    !!!!!!!!!!!!!!!!!!!!
    !write(*, *) "Xp, Yp, DUT1, HG:"
    !write(*, *) Xp, Yp, DUT1
    !write(*, "(3(2x,F18.5))") HG(1, :)
    !write(*, "(3(2x,F18.5))") HG(2, :)
    !write(*, "(3(2x,F18.5))") HG(3, :)
    
    endsubroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This subroutine is written to calculate HG Matrix
    !t: MJD(UTC)
    !HG:3*3 matrix
    subroutine HG_get(t, HG)
    implicit none
    double precision t,t1,t2
    double precision zetaA,zA,thetaA,delta_mu,delta_theta,delta_epsilon,SG,x_p,y_p,pi_const
    double precision HG(3,3),HG1(3,3),EP(3,3),EP1(3,3),EP2(3,3),ER(3,3),GR(3,3),NR(3,3),NR1(3,3)&
        &,PR(3,3),PR1(3,3)
    double precision matrix_deltaepsilon(3,3),matrix_deltatheta(3,3)
    double precision matrix_deltamu(3,3),matrix_zA(3,3),matrix_thetaA(3,3),matrix_zetaA(3,3)
    double precision UTC,TT,UT1
    parameter(pi_const=3.141592653589793D0) 
    !---------------------------------------------------------------

    !Synchronize EOP coefficients
    call syncEOP(t)

    !Get UT1 and UTC
    UTC = t
    call UTC2UT1(UTC, DUT1, UT1)
    call UTC2TT(UTC, TT)
    t1=(UT1-51544.5D0)/36525D0
    t2=(TT-51544.5D0)/36525D0

    !get EP
    call EP_get(Xp, Yp, EP)

    !get ER
    call ER_get(t1, ER)

    !get NR
    call NR_get(t2, NR)

    !get PR
    call PR_get(t2, PR)

    !HG = EP*ER*NR*PR
    GR = matmul(NR, PR)
    HG1 = matmul(ER, GR)
    HG = matmul(EP, HG1)
    return
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine EP_get(Xp, Yp, EP)
    implicit none
    double precision Xp, Yp, EP1(3, 3), EP2(3, 3), EP(3, 3)
    double precision pi_const
    pi_const = 3.141592653589793d0
    call Rxmatrix(-Yp,EP1)
    call Rymatrix(-Xp,EP2)
    EP = matmul(EP2, EP1)
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ER_get(t, ER)
    implicit none
    double precision t,SG,mean_SG,delta_mu,delta_theta,delta_epsilon
    double precision pi2,pi,S
    double precision ER(3, 3)
    parameter(pi=3.141592653589793D0)
!    t=1d4/36525D0

    pi2=pi*2d0

    call nutation_calculation(t,delta_mu,delta_theta,delta_epsilon)

!    mean_SG=18.6973746d0/12D0*pi+879000.0513367d0*t/12D0*pi+0.093104d0*t**2D0*pi/43200D0
!

    mean_SG=(18.6973746D0+879000.0513367D0*t+(0.093104D0*t**2)/3600.0D0)/24D0*pi2

    mean_SG=dmod(mean_SG,pi2)

    SG=mean_SG+delta_mu

    call Rzmatrix(SG,ER)
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine NR_get(t2, NR)
    implicit none
    double precision t2, delta_mu, delta_theta, delta_epsilon
    double precision matrix_deltaepsilon(3, 3), matrix_deltatheta(3, 3), matrix_deltamu(3, 3), NR1(3, 3), NR(3, 3)

    call nutation_calculation(t2,delta_mu,delta_theta,delta_epsilon)
    call Rxmatrix(-delta_epsilon,matrix_deltaepsilon)
    call Rymatrix(delta_theta,matrix_deltatheta)
    call Rzmatrix(-delta_mu,matrix_deltamu)
    NR1 = matmul(matrix_deltatheta,matrix_deltamu)
    NR = matmul(matrix_deltaepsilon,NR1)
    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine PR_get(t2, PR)
    implicit none
    double precision zetaA, thetaA, zA, matrix_zA(3, 3)
    double precision matrix_thetaA(3, 3), matrix_zetaA(3, 3), PR1(3, 3), PR(3, 3)
    double precision pi, t2
    parameter(pi=3.141592653589793d0)

    zetaA=(2306.2181d0*t2+0.30188d0*t2**2)/648000d0*pi
    zA=(2306.2181d0*t2+1.09468d0*t2**2)/648000d0*pi
    thetaA=(2004.3109d0*t2-0.42665d0*t2**2)/648000d0*pi


    call Rzmatrix(-zA,matrix_zA)
    call Rymatrix(thetaA,matrix_thetaA)
    call Rzmatrix(-zetaA,matrix_zetaA)
    PR1 = matmul(matrix_thetaA,matrix_zetaA)
    PR = matmul(matrix_zA,PR1)

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine nutation_calculation(t,delta_mu,delta_theta,delta_epsilon)
    implicit none
    integer I,J
    double precision t,delta_psi,delta_mu,delta_theta,delta_epsilon,e
    
    double precision,save :: A0(1:20),B0(1:20),A1(1:20),B1(1:20)
    logical,save :: firstRun = .true.
    
    double precision k(1:20,1:5)
    double precision alpha(1:5)
    double precision psi_phi1(1:20)
    double precision epsilon_phi2(1:20)
    double precision pi

    parameter(pi=3.141592653589793d0)
    data ((k(J,I),I=1,5),J=1,20)&
    &/0,0,0,0,1D0,&
    &0,0,2D0,-2D0,2D0,&
    &0,0,2D0,0,2D0,&
    &0,0,0,0,2D0,&
    &0,1D0,0,0,0,&
    &1D0,0,0,0,0,&
    &0,1D0,2D0,-2D0,2D0,&
    &0,0,2D0,0,1D0,&
    &1D0,0,2D0,0,2D0,&
    &0,-1D0,2D0,-2D0,2D0,&
    &1D0,0,0,-2D0,0,&
    &0,0,2D0,-2D0,1D0,&
    &-1D0,0,2D0,0,2D0,&
    &1D0,0,0,0,1D0,&
    &0,0,0,2D0,0,&
    &-1D0,0,2D0,2D0,2D0,&
    &-1D0,0,0,0,1D0,&
    &1D0,0,2D0,0,1D0,&
    &2D0,0,0,-2D0,0,&
    &-2D0,0,2D0,0,1D0/

    if (firstRun) then
        !Get A0, A1, B0, B1 from file
        open(302,file='data/AtlasData/Nutation.dat',status='old')
        read(302,*) A0,A1,B0,B1
        close(302)
        do I = 1,20
            A0(I)=A0(I)*pi/180d0/3600d0
            A1(I)=A1(I)*pi/180d0/3600d0
            B0(I)=B0(I)*pi/180d0/3600d0
            B1(I)=B1(I)*pi/180d0/3600d0
        enddo
        firstRun = .false.
    endif
    
    alpha(1)=(134D0+57.0d0/60D0+46.733d0/3600D0+(198D0+52.0d0/60D0+2.633d0/3600D0)*t+&
                & 31.310d0*t**2.0d0/3600D0)/180D0*pi+(1325D0*2D0)*t*pi

    alpha(2)=(357D0+31.0d0/60D0+39.804d0/3600D0+(359D0+3.0d0/60D0+1.224d0/3600D0)*t&
                &-0.577d0*t**2.0d0/3600D0)/180D0*pi+99D0*2D0*t*pi+2*pi

    alpha(3)=(93D0+16.0d0/60D0+18.877d0/3600D0+(82D0+1.0d0/60D0+3.137d0/3600D0)*t&
                &-13.257d0*t**2.0d0/3600D0)/180D0*pi+1342D0*2D0*t*pi

    alpha(4)=(297D0+51.0d0/60D0+1.307d0/3600D0+(307D0+6.0d0/60D0+41.328d0/3600D0)*t-&
                & 6.891d0*t**2.0d0/3600D0)/180D0*pi+1236D0*2D0*t*pi+2*pi

    alpha(5)=(125D0+2.0d0/60D0+40.280d0/3600D0-(134D0+8.0d0/60D0+10.539d0/3600D0)*t&
                &+7.455d0*(t**2.0d0)/3600D0)/180D0*pi-5D0*2D0*t*pi+2*pi

    alpha=dmod(alpha,2D0*pi)

    delta_psi=0d0
    delta_epsilon=0d0
    do J=1,20
        psi_phi1(J)=k(J,1)*alpha(1)+k(J,2)*alpha(2)+k(J,3)*alpha(3)+ &
        k(J,4)*alpha(4)+k(J,5)*alpha(5)
        delta_psi=delta_psi+(A0(J)+A1(J)*t)*dsin(psi_phi1(J))
    end do

    do J=1,20
        epsilon_phi2(J)=k(J,1)*alpha(1)+k(J,2)*alpha(2)+k(J,3)*alpha(3)+&
        k(J,4)*alpha(4)+k(J,5)*alpha(5)
        delta_epsilon=delta_epsilon+(B0(J)+B1(J)*t)*dcos(epsilon_phi2(J))
    end do

    e=(23.0D0+26.0D0/60D0+21.448D0/3600D0-46.815D0*t/3600D0)/180D0*pi

    delta_mu=delta_psi*dcos(e)
    delta_theta=delta_psi*dsin(e)

    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine matrix_transform(matrix_in,matrix_out)
    implicit none
    double precision matrix_in(3,3),matrix_out(3,3)
    integer I,J
    do I=1,3
        do J=1,3
            matrix_out(I,J)=matrix_in(J,I)
        end do
    end do
    return
    end subroutine matrix_transform

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Get_ZR(longi, lati, ZR)
    implicit none
    double precision longi, lati, pi_const
    double precision R1(3, 3), R2(3, 3), R3(3, 3), R4(3, 3), ZR(3, 3)
    pi_const = 3.141592653589793d0
    call coordinate_reverse('y', R1)
    call Rzmatrix(pi_const*3.0d0/2d0, R2)
    call Rxmatrix(-(pi_const/2d0 - lati), R3)
    call Rzmatrix(longi - pi_const/2d0, R4)
    ZR = matmul(R1, matmul(R2, matmul(R3, R4)))
    end subroutine    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine coordinate_reverse(axis, R)
    character axis
    double precision R(3, 3)
!    DATA R / 1d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 1d0/
    R(1, 1) = 1
    R(2, 1) = 0
    R(3, 1) = 0
    R(1, 2) = 0
    R(2, 2) = 1
    R(3, 2) = 0
    R(1, 3) = 0
    R(2, 3) = 0
    R(3, 3) = 1
    if (axis == 'x') then
        R(1, 1) = -1d0
    else if (axis == 'y') then
        R(2, 2) = -1d0
    else if (axis == 'z') then
        R(3, 3) = -1d0
    else
        write(*, *) "Input error:: subroutine reverse"

    end if
    end subroutine    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Rzmatrix(theta,a)
    implicit none
    double precision theta,a(3,3)
    a(1,1)=dcos(theta)
    a(1,2)=dsin(theta)
    a(1,3)=0
    a(2,1)=-dsin(theta)
    a(2,2)=dcos(theta)
    a(2,3)=0
    a(3,1)=0
    a(3,2)=0
    a(3,3)=1D0
    return
    end subroutine Rzmatrix

    subroutine Rymatrix(theta,a)
    implicit none
    double precision theta,a(3,3)
    a(1,1)=dcos(theta)
    a(1,2)=0
    a(1,3)=-dsin(theta)
    a(2,1)=0
    a(2,2)=1D0
    a(2,3)=0
    a(3,1)=dsin(theta)
    a(3,2)=0
    a(3,3)=dcos(theta)
    return
    end subroutine Rymatrix

    subroutine Rxmatrix(theta,a)
    implicit none
    double precision theta,a(3,3)
    a(1,1)=1D0
    a(1,2)=0
    a(1,3)=0
    a(2,1)=0
    a(2,2)=dcos(theta)
    a(2,3)=dsin(theta)
    a(3,1)=0
    a(3,2)=-dsin(theta)
    a(3,3)=dcos(theta)
    return
    end subroutine Rxmatrix

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    endmodule
