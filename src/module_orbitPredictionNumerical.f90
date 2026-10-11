
    module module_orbitPredictionNumerical
    use module_const
    use module_space
    use module_atmosphere
    use module_perturbationSRP
    use module_perturbationThirdBody
    use module_perturbationPN
    use module_thruster
    use module_Gravity
    use EarthRadiationPressure

    ! 摄动力开关, 索引顺序与 module_satellite%opMod_forceSwitchArray 一一对应
    ! (地球引力场不在这组开关里, 由下面的 degree/order 直接控制)
    integer, private :: k1 = 1 ! lunar gravity
    integer, private :: k2 = 0 ! lunar tide
    integer, private :: k3 = 1 ! solar gravity
    integer, private :: k4 = 0 ! solar tide
    integer, private :: k5 = 0 ! srp
    integer, private :: k6 = 0 ! drag
    integer, private :: k7 = 0 ! post-Newton effect

    integer, private :: degree = 2 ! Earth gravity field degree n
    integer, private :: order = 2  ! Earth gravity field order m
    
    integer, private :: orderMoon = 0        ! Moon gravity order, two-body model in default
    
    integer, private :: Mercury = 0          ! Include Mercury gravity (0=No, 1=Yes)
    integer, private :: Venus = 0            ! Include Venus gravity (0=No, 1=Yes)
    integer, private :: Mars = 0             ! Include Mars gravity (0=No, 1=Yes)
    integer, private :: Jupiter = 0          ! Include Jupiter gravity (0=No, 1=Yes)
    integer, private :: Saturn = 0           ! Include Saturn gravity (0=No, 1=Yes)
    integer, private :: Uranus = 0           ! Include Uranus gravity (0=No, 1=Yes)
    integer, private :: Pluto = 0            ! Include Pluto gravity (0=No, 1=Yes)
    
    ! about atmosphere
    double precision, private :: Cd = 2.2d0 ! 
    double precision, private :: amrDrag = 0.02d0 ! 
    double precision, private :: F107 = 130d0 !
    double precision, private :: Ap = 15d0 !
    
    ! about srp
    integer :: SRP_Model = 1 ! 1 - canon ball; 2 - ECOM (9 paras)
    double precision, private :: reflectivity = 0.0d0 !
    double precision, private :: amrSRP = 0.02d0 !
    double precision, private :: ECOM(9) = 0d0

    ! about ERP
    integer, private :: key_ERP = 0
    double precision, private :: reflectivityERP = 0.0d0
    double precision, private :: amrERP = 0.02d0
    
    ! stepLength = 0 -> adaptive step length; 
    ! stepLength > 0 -> fixed step length
    double precision, private :: stepLength = 0d0 ! sec
    
    integer, public :: satNum = 0 ! 
    logical, private :: firstRun = .true. ! 
    
    integer, private :: eleType = 3 ! 
    
    double precision, public :: mass = 500d0
    type(type_thruster), public ::  thisThruster

    type(GravityField), private :: MoonGravityField
    type(GravityField), private :: EarthGravityField
    
    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine NumOP(K_, degree_, order_, &
        &   Cd_, amrDrag_, F107_, Ap_, &
        &   reflectivity_, amrSRP_, &
        &   SRP_model_, ECOM_, &
        &   step_, &
        &   orderMoon_, Mercury_, Venus_, Mars_, Jupiter_, Saturn_, Uranus_, Pluto_, &
        &   key_ERP_, reflectivityERP_, amrERP_, &
        &   eleType_, satNum_, &
        &   x0, t0, dt, x1)
    ! use OMP_LIB
    implicit none
    integer, intent(in) :: K_(7), degree_, order_
    double precision, intent(in) :: Cd_, amrDrag_, F107_, Ap_, reflectivity_, amrSRP_
    integer, intent(in) :: SRP_model_
    double precision, intent(in) :: ECOM_(9)
    double precision, intent(in) :: step_
    integer, intent(in) :: orderMoon_, Mercury_, Venus_, Mars_, Jupiter_, Saturn_, Uranus_, Pluto_
    integer, intent(in) :: key_ERP_
    double precision, intent(in) :: reflectivityERP_, amrERP_
    integer, intent(in) :: eleType_, satNum_
    double precision, intent(in) :: x0(6), t0, dt
    double precision, intent(out) :: x1(6)
    double precision PV0(6), kep0(6), pi, t1, emphemeris_step
    double precision r(3), v(3)
    double precision nsg10(6), nsg11(6)
    integer threadNum
    !----------------------------------------------------------------
    
    ! If the height is too low
    select case (eleType_)
    case (STATE_TYPE_KEPLERIAN)
        if ((x0(1)*(1d0 - x0(2)) - 1d0)*6378.1363d0 <= 100d0) then
            x1 = 0d0
            write(*, *) "stop prediction as the orbit is lower than 100 km."
            return
        endif
    case (STATE_TYPE_CARTESIAN)
        if (norm2(x0(1:3))*6378.1363d0 <= 6478.1363d0) then
            x1 = 0d0
            write(*, *) "stop prediction as the orbit is lower than 100 km."
            return
        endif
    case default
        
    endselect
    
    ! If prediction time is 0
    if (dt  ==  0d0)    then
        x1  =   x0
        return
    endif
    
    k1 = k_(1)
    k2 = k_(2)
    k3 = k_(3)
    k4 = k_(4)
    k5 = k_(5)
    k6 = k_(6)
    k7 = k_(7)
    degree = degree_
    order = order_

    orderMoon = orderMoon_
    Mercury = Mercury_
    Venus = Venus_
    Mars = Mars_
    Jupiter = Jupiter_
    Saturn = Saturn_
    Uranus = Uranus_
    Pluto = Pluto_

    Cd = Cd_
    !threadNum = OMP_get_thread_num()
    !amrDrag(threadNum + 1) = amrDrag_
    amrDrag = amrDrag_
    F107 = F107_
    Ap = Ap_
    reflectivity = reflectivity_
    amrSRP = amrSRP_
    
    SRP_Model = SRP_model_
    ECOM = ECOM_

    key_ERP = key_ERP_
    reflectivityERP = reflectivityERP_
    amrERP = amrERP_

    stepLength = step_
    
    eleType = eleType_
    satNum = satNum_

    ! initialize Moon gravity field
    if (orderMoon > 0) then
        if (.not. MoonGravityField%initialized) then
            call MoonGravityField%init( &
                "MoonGravity", "Moon", orderMoon, orderMoon, &
                "data/AtlasData/grail.txt")
        endif
    endif

    ! initialize Earth gravity field (系数数组按 degree x order 分配, 后续查询
    ! 会各自截断到这两个上界, 不会越界)
    if (degree > 1) then
        if (.not. EarthGravityField%initialized) then
            call EarthGravityField%init( &
                "EarthGravity", "Earth", degree, order)
            EarthGravityField%GM = EARTH_GM
            EarthGravityField%referenceRadius = LENGTH_UNIT * 1000d0
            !call EarthGravityField%read_gravity_coefficients_4col( &
            !    "data/AtlasData/GGM03C.txt")
            call EarthGravityField%read_gravity_coefficients_4col( &
                "data/AtlasData/EGM2008.txt")
        else if (degree > EarthGravityField%maxdegree .or. &
            &    order  > EarthGravityField%maxOrder) then
            ! re-initialize if higher degree/order is requested
            call EarthGravityField%init( &
                "EarthGravity", "Earth", degree, order)
            EarthGravityField%GM = EARTH_GM
            EarthGravityField%referenceRadius = LENGTH_UNIT * 1000d0
            !call EarthGravityField%read_gravity_coefficients_4col( &
            !    "data/AtlasData/GGM03C.txt")
            call EarthGravityField%read_gravity_coefficients_4col( &
                "data/AtlasData/EGM2008.txt")
        endif
    endif
    
    select case (eleType_) 
    case (STATE_TYPE_KEPLERIAN)
        ! kepler elements -> 1st_nonsingular_elements
        call kpl2nsg1(x0, nsg10)
    
        ! Orbit prediction
        call prediction(0, nsg10, t0, dt, nsg11)
        
        nsg11(2) = dmod(nsg11(2), 3.141592653589793d0)
        nsg11(3) = dmod(nsg11(3), 3.141592653589793d0*2d0)
        nsg11(6) = dmod(nsg11(6), 3.141592653589793d0*2d0)
        
        ! 1st_nonsingular_elements -> kepler elements
        call nsg12kpl(nsg11, x1)
    case (STATE_TYPE_NONSINGULAR_1) 
        nsg10 = x0
        ! Orbit prediction
        call prediction(0, nsg10, t0, dt, nsg11)
        nsg11(2) = dmod(nsg11(2), 3.141592653589793d0)
        nsg11(3) = dmod(nsg11(3), 3.141592653589793d0*2d0)
        nsg11(6) = dmod(nsg11(6), 3.141592653589793d0*2d0)
        x1 = nsg11
    case (STATE_TYPE_NONSINGULAR_2)
        
    case (STATE_TYPE_CARTESIAN)
        call prediction(1, x0, t0, dt, x1)
    case default
        
    end select
        
    endsubroutine 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  prediction(usePV, elements0, t0, dt, elements1)
    use, intrinsic :: ieee_arithmetic
    implicit    none
    double precision      elements0(6), t0, dt, elements1(6), life
    logical     is_finished
    double precision      step, t_1
    double precision      r0(6), r1(6), t_object, t, error_max
    double precision      deviation
    double precision      elements(6), roots(6)
    double precision      ::  pi  =   3.141592653589793d0
    double precision      t_start
    integer     usePV
    logical :: fixedStep = .false.
    integer     nReduce                       ! 当前步已连续折半步长的次数
    integer, parameter :: MAX_REDUCE = 60     ! 连续折半上限: 2^-60 仍不收敛即判失败
    !---------------------------------------------------------------
    
    if (dabs(stepLength) > 1d-16) then ! fixed step length
        step = stepLength
        fixedStep = .true.
    else
        step = 60d0
    endif
    
    if (dt      <   0d0) step = -step 

    t_object    =   dt                                              ! t0 -------- t0 + t ----------t0 + t_object
    t           =   0d0

                                                                    !normalize
    t_start     =   t0*86400d0/time_unit
    t           =   t*86400d0/time_unit
    t_object    =   t_object*86400d0/time_unit
    step        =   step/time_unit
    error_max   =   1.5d-13

    do while (dabs(t)   <   dabs(t_object))
                                                                    !if this is the last step
        if (abs(t + step)   >=  abs(t_object)) then
            step            =   t_object - t
        end if

        nReduce = 0

        do while ( .true. )

            call rkf78_forOP(usePV, t_start + t, step, elements0, elements1, deviation)

            ! 非有限值检查: 力模型一旦返回 NaN/Inf(例如星历文件读取失败, 使三体
            ! 引力里的 r/|r|**3 成为 0/0), deviation 与 elements1 都会是 NaN。
            ! NaN 参与的比较恒为假, 下面的 "deviation <= error_max" 于是永远不
            ! 成立, 程序就在本循环里无限折半步长 —— 表现为 CPU 占满却没有任何
            ! 输出, 是最难排查的失败形式。这里显式报错并终止。
            ! 诊断一律写 stderr(unit 0), 与 module_clock 的约定一致; 不要写
            ! unit 6 —— 本过程可能在 I/O 语句中被调用, 递归 I/O 会死锁。
            if (.not. ieee_is_finite(deviation) .or. any(.not. ieee_is_finite(elements1))) then
                write(0, *) "Error @prediction: non-finite value in numerical propagation."
                write(0, *) "  deviation   = ", deviation
                write(0, *) "  epoch (MJD) = ", (t_start + t)*time_unit/86400d0
                write(0, *) "  step (TU)   = ", step
                write(0, *) "  A force-model term produced NaN/Inf. Check its inputs"
                write(0, *) "  (e.g. the ephemeris data) and the force-model switches."
                stop
            endif
            
            !if the deviation is small enough
            if (deviation <= error_max .or. fixedStep) then
                t = t + step
                
                if (error_max >= deviation*100d0 .and. (.not. fixedStep)) then
                    step = step*2.0d0
                endif
                
                elements0 = elements1
                
                if (usePV == 0) elements0(6) = dmod(elements0(6), 2d0*pi)
                
                exit
            else !if the deviation is too big
                step = step/2.0d0
                nReduce = nReduce + 1
                ! 兜底: 即使 deviation 是有限值, 也可能因为力模型不连续等原因始终
                ! 超差。连续折半 60 次(约 1e-18 倍)仍不收敛即判失败, 同样避免死循环。
                if (nReduce > MAX_REDUCE) then
                    write(0, *) "Error @prediction: step size collapsed, no convergence."
                    write(0, *) "  epoch (MJD) = ", (t_start + t)*time_unit/86400d0
                    write(0, *) "  step (TU)   = ", step
                    write(0, *) "  deviation   = ", deviation, "   tolerance = ", error_max
                    stop
                endif
                continue
            end if
               
        enddo
    enddo
    
    fixedStep = .false.
        
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    subroutine rkf78_forOP(usePV, t0, h, x0, x1, deviation)
    implicit none
    double precision t0, h, x0(6), x1(6), deviation
    double precision alpha(0:12), beta(0:12, 0:11), c(0:12), c_hat(0:12)
    double precision f(0:12, 6)
    double precision t, x(6)
    integer i, j
    double precision roots(6)
    integer usePV
    !----------------------------------------------------------------------------------------------------------------

    data alpha  /0,2.5000000000000000D-01,8.3333333333333329D-02,1.2500000000000000D-01,4.0000000000000002D-01,0.5D0,&
            & 8.5714285714285710D-01,1.4285714285714285D-01,6.6666666666666663D-01,2.8571428571428570D-01,&
            & 1.0000000000000000D+00,3.3333333333333331D-01,1.0000000000000000D+00/

    data ((beta(i,j),j=0,11),i=0,12)&
    & /0,0,0,0,0,0,0,0,0,0,0,0,&
    & 2.5000000000000000D-01,0,0,0,0,0,0,0,0,0,0,0,&
    & 6.9444444444444448D-02,1.3888888888888888D-02,0,0,0,0,0,0,0,0,0,0,&
    & 3.1250000000000000D-02,0,9.3750000000000000D-02,0,0,0,0,0,0,0,0,0,&
    & 8.4799999999999998D-01,0,-3.2639999999999998D+00,2.8159999999999998D+00,0,0,0,0,0,0,0,0,&
    & 2.0833333333333332D-02,0,0,2.4242424242424243D-01,2.3674242424242425D-01,0,0,0,0,0,0,0,&
    & -5.2603082049146188D-01,0,0,1.5120972322138504D+00,-2.4279656203854456D+00,2.2990420658059141D+00,0,0,0,0,0,0,&
    & 9.4387755102040810D-02,0,0,0,1.7272534013605442D-01,-1.3333333333333333D-01,9.0773809523809531D-03,0,0,0,0,0,&
    & 6.7317264354301387D-01,0,0,-1.8461297276112092D+00,2.6289682539682540D+00,&
    & -1.6379384675680972D+00,1.3879972565157750D-01,7.0979423868312752D-01,0,0,0,0,&
    & -1.3235596279328057D-01,0,0,-3.9288680221624829D-02,-1.0043716002170813D+00,8.4647396917950857D-01,&
    & -2.6142857142857141D-02,7.2727272727272729D-01,-8.5873310363106276D-02,0,0,0,&
    & 1.6753947368421052D+00,0,0,-3.8338755980861245D+00,5.6878613437001597D+00,&
    & -5.9788070175438595D+00,3.7195942982456143D-01,-2.4617224880382776D-01,&
    & 1.2080965909090908D+00,2.1155427631578947D+00,0,0,&
    & -1.3257397609249461D-01,0,0,1.4411487004079596D-01,-1.1304972676500453D+00,&
    & 6.9222026259063296D-01,-6.1460048010973940D-02,4.3555555555555553D-01,6.2500000000000000D-02,3.2347393689986281D-01,0,0,&
    & 1.2906682027649770D+00,0,0,-3.1959782153330543D+00,5.1372560571149277D+00,-5.1641167434715820D+00,&
    & 4.1135691593352885D-01,5.5178885630498531D-01,8.6849340175953083D-01,-1.0372983870967742D+00,0,2.1378299120234603D+00/


    data c /4.5138888888888888D-02,0,0,0,0,2.5600000000000001D-01,2.1675694444444443D-01, &
    & 1.9402020202020201D-01,1.2080965909090909D-01,1.2505208333333334D-01,3.5185185185185187D-02,0,0/


    data c_hat /4.3055555555555555D-02,0,0,0,0,2.1333333333333335D-01,&
                & 2.1220959595959596D-01,2.1220959595959596D-01,&
                & 1.3806818181818181D-01,0,0,1.3806818181818181D-01,4.3055555555555555D-02/


    x1 = x0
    do i = 0, 12
        !get t and x
        if (i == 0) then
            t = t0
            x = x0
        else
            t = t0 + alpha(i)*h
            x = x0
            do j = 0, (i - 1)
                x = x + h*beta(i, j)*f(j, :)
            end do
        end if
        !get f(i)
        if (usePV == 0) then
            call right_function2(t, x, f(i, :))
        else
            call right_function3(t, x, f(i, :))
        endif
        x1 = x1 + h*c_hat(i)*f(i, :)
    end do

    !truncation error
    deviation = norm2(-2.0833333333333333D-03*f(0, :) - 4.2666666666666665D-02*f(5, :) &
    & - 4.5473484848484848D-03*f(6, :) + 1.8189393939393939D-02*f(7, :) &
    & + 1.7258522727272727D-02*f(8, :) - 1.2505208333333334D-01*f(9, :) &
    & - 4.2222222222222223D-02*f(10, :) + 1.3806818181818181D-01*f(11, :) &
    & + 4.3055555555555555D-02*f(12, :))*dabs(h)

    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    !'subroutine right_function' uses the first kind nonsingular elements    
    !'subroutine right_function2' uses the second kind nonsingular elements
    !'subroutine right_function3' uses position and velocity
    subroutine right_function(UTC, elements, f)
    implicit none
    double precision UTC, elements(6), f(6)
    double precision force(3), r(6), elements1(6)
    double precision a, e, i, omega_up, omega_low, M
    double precision EE, ff, a_r, p_r
    double precision xi, eta, h, k, lambda
    double precision n, p, omega_low_tiao, u, u_tiao
    double precision S, T, W
    !---------------------------------------------------------------------

    UTC = UTC*time_unit/86400d0     !unit -> second -> day

    !2rd_nonsingular_elements -> elements
    call nsg22kpl(elements, elements1)

    !elements -> r
    call kpl2cts(elements1, r)

    !get perturbation force
    call get_force(UTC, r, force)

    a = elements1(1)
    e = elements1(2)
    i = elements1(3)
    omega_up = elements1(4)
    omega_low = elements1(5)
    M = elements1(6)


    !Solve the kepler equation
    !E = M + e*sin(E)
    call solve_kepler_equation(e, M, EE)

    a_r = (1d0 - e*dcos(EE))**(-1d0)
    p_r = a_r*(1d0 - e**2)
    !get ƫ����� ff
    ff = EE + dasin(a_r*e*dsin(EE)*(1d0 - 1d0/(1d0 + dsqrt(1d0 - e**2))*e*dcos(EE)))

    xi = elements(2)
    eta = elements(3)
    h = elements(4)
    k = elements(5)
    lambda = elements(6)

    n = a**(-1.5d0)
    p = a*(1d0 - e**2)
    omega_low_tiao = omega_low + omega_up
    u = ff + omega_low_tiao
    u_tiao = EE + omega_low_tiao

    !cartesian -> STW
    call get_STW(elements, u, force, S, T, W)

    f(1) = 2d0/n/dsqrt(1d0 - e**2)*(S*(xi*dsin(u) - eta*dcos(u)) + T*p_r)

    f(2) = dsqrt(1d0 - e**2)/n/a*(S*dsin(u) &
    & + T*(dcos(u_tiao) + dcos(u) - eta/(dsqrt(1d0 - e**2)*(1d0 + dsqrt(1d0 - e**2)))*(xi*dsin(u_tiao) - eta*dcos(u_tiao))) &
    & + W/p_r*(-eta*(h*dsin(u) - k*dcos(u))/dcos(i/2d0)))

    f(3) = dsqrt(1d0 - e**2)/n/a*(-S*dcos(u) &
    & + T*(dsin(u_tiao) + dsin(u) + xi/(dsqrt(1d0 - e**2)*(1d0 + dsqrt(1d0 - e**2)))*(xi*dsin(u_tiao) - eta*dcos(u_tiao))) &
    & + W/p_r*(xi*(h*dsin(u) - k*dcos(u))/dcos(i/2d0)))

    f(4) = W/(2d0*n*a*dsqrt(1 - e**2)*dcos(i/2d0))/a_r*(dcos(u) - h*(k*dsin(u) + h*dcos(u)))

    f(5) = W/(2d0*n*a*dsqrt(1 - e**2)*dcos(i/2d0))/a_r*(dsin(u) - k*(k*dsin(u) + h*dcos(u)))

    f(6) = n - dsqrt(1d0 - e**2)/n/a*(2d0*S*dsqrt(1d0 - e**2)/(a_r*(1d0 - e**2)) &
    & + 1d0/(1d0 + dsqrt(1d0 - e**2))*(S*(xi*dcos(u) + eta*dsin(u)) - T*(1d0 + 1/p_r)*(xi*dsin(u) - eta*dcos(u)))&
    & - W/p_r*(h*dsin(u) - k*dcos(u))/dcos(i/2d0))

    end subroutine
    !'subroutine right_function' uses the first kind nonsingular elements    
    !'subroutine right_function2' uses the second kind nonsingular elements
    !'subroutine right_function3' uses position and velocity
    subroutine right_function2(UTC, roots, f)
    implicit none
    double precision UTC, roots(6), f(6)
    double precision a, i, omega_up, xi, eta, lambda
    double precision e, elements(6), ff, EE
    double precision a_r, p_r
    double precision u, u_tiao, n, xx
    double precision rs(6), force(3), S, T, W
    !-------------------------------------------
    UTC = UTC*time_unit/86400d0
    a = roots(1)
    i = roots(2)
    omega_up = roots(3)
    xi = roots(4)
    eta = roots(5)
    lambda = roots(6)

    call nsg12kpl(roots, elements)
    e = elements(2)
    call solve_kepler_equation(e, elements(6), EE)

    a_r = (1d0 - e*dcos(EE))**(-1d0)
    ff = EE + dasin(a_r*e*dsin(EE)*(1d0 - 1d0/(1d0 + dsqrt(1d0 - e**2))*e*dcos(EE)))
    u = ff + elements(5)
    u_tiao = EE + elements(5)

!    a_r = (1d0 - (xi*dcos(u_tiao) + eta*dsin(u_tiao)))**(-1d0)

    p_r = a_r*(1d0 - e**2)
    n = a**(-1.5d0)
    xx = dsqrt(1d0 - e**2)
    call kpl2cts(elements, rs)
    call get_force(UTC, rs, force)
    call get_STW(elements, u, force, S, T, W)
    
    f(1) = 2d0/n/xx*(S*(xi*dsin(u) - eta*dcos(u)) + T*p_r)
    f(2) = dcos(u)/n/a/a_r/xx*W
    f(3) = dsin(u)/n/a/a_r/xx/dsin(i)*W
    f(4) = eta*dcos(i)*f(3) + xx/n/a*(S*dsin(u) + T/xx*((dcos(u_tiao) + xx*dcos(u)) &
    & - xi/(1d0 + xx)*(xi*dcos(u_tiao) + eta*dsin(u_tiao))))
    f(5) = -xi*dcos(i)*f(3) -  xx/n/a*(S*dcos(u) - T/xx*((dsin(u_tiao) + xx*dsin(u)) &
    & - eta/(1d0 + xx)*(xi*dcos(u_tiao) + eta*dsin(u_tiao))))
    f(6) = n - dcos(i)*f(3) - xx/n/a*(2d0*S*xx/p_r &
    & + 1d0/(1d0 + xx)*(S*(xi*dcos(u) + eta*dsin(u)) - T*(1d0 + 1d0/p_r)*(xi*dsin(u) - eta*dcos(u))))

    end subroutine
    !'subroutine right_function' uses the first kind nonsingular elements    
    !'subroutine right_function2' uses the second kind nonsingular elements
    !'subroutine right_function3' uses position and velocity    
    subroutine right_function3(utc, x0, f)
    implicit none
    double precision, intent(in)  :: utc, x0(6)
    double precision, intent(out)  :: f(6)
    double precision  :: t, force(3), F0(3)
    
    t = UTC*time_unit/86400d0
    call get_force(t, x0, force) ! perturbation force
    
    F0 = -1d0/(norm2(x0(1:3)))**3*x0(1:3) ! two-body gravity
    
    f(1:3) = x0(4:6)
    f(4:6) = F0 + force
    endsubroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_force(UTC, rs_in, force)
    implicit none
    double precision, intent(in) :: UTC, rs_in(6)
    double precision, intent(out) :: force(3)
    double precision TT, rs(6)
    double precision force1(3), force2(3), force3(3), force4(3), force5(3), force6(3)
    double precision force7(3)
    double precision :: massNow   ! Current satellite mass, kg (thrust acceleration)
    double precision solarPV(6)
    integer threadNum
    double precision :: LunarPV(6), TT_plus, r_lunar(3), rSat_LCI(3), rSat_LCF(3)
    double precision :: force_moonAspherical(3), force_moonAspherical_LCI(3)
    double precision :: rs_ecf(3), rs_ecf_fixed(3), force_eci(3)
    !---------------------------------------------------------------
    
    ! normalized unit -> km, km/s
    rs(1:3) = rs_in(1:3)*length_unit
    rs(4:6) = rs_in(4:6)*length_unit/time_unit
    
    ! Aspherical gravity (using GravityField)
    ! degree = 阶数, order = 次数; degree < 2 时整个非球形引力项为 0
    if (degree == 0 .OR. degree == 1 .OR. (.not. EarthGravityField%initialized)) then
        force = 0d0
    else
        ! ECI -> ECEF (position in meters)
        rs_ecf = rs(1:3) * 1000d0
        call ECI2ECF_IAU1980_position(UTC, rs_ecf, rs_ecf_fixed)
        ! Compute acceleration in ECEF (m/s^2)
        call EarthGravityField%acceleration(rs_ecf_fixed, degree, order, force_eci)
        ! ECEF -> ECI
        call ECF2ECI_IAU1980_position(UTC, force_eci, force)
    endif

    ! Solar radiative pressure
    if (k5 == 1) then
        call solar_pressure_perturbation(UTC, rs, reflectivity, amrSRP, SRP_Model, ECOM, force4)
        force = force + force4
    endif

    ! Earth Radiation Pressure (ERP)
    if (key_ERP == 1) then
        call UTC2TT(UTC, TT)
        call planetsPV(TT, nEarth, nSun, solarPV)
        call calEarthRadiationPressure(UTC, &
            solarPV(1:3) / LENGTH_UNIT, &
            rs(1:3) / LENGTH_UNIT, &
            amrERP, reflectivityERP, force2)
        force2 = force2 * acc_unit
        force = force + force2
    endif

    ! Atmospherical drag
    if (k6 == 1) then
        !threadNum = OMP_get_thread_num()
        call atmospherical_drag_perturbation(UTC, rs(1:3), rs(4:6), Cd, amrDrag, F107, Ap, force5)
        force = force + force5
    endif

    ! Third bodies' gravity and tidal effect
    if (k1 == 1 .or. k3 == 1 .or. Mercury == 1 .or. Venus == 1 .or. Mars == 1 .or. &
        Jupiter == 1 .or. Saturn == 1 .or. Uranus == 1 .or. Pluto == 1) then
        call third_body_perturbation(UTC, rs(1:3), k1, k2, k3, k4, &
                                     Mercury, Venus, Mars, Jupiter, Saturn, Uranus, Pluto, force3)
        force = force + force3
    endif

    ! Lunar aspherical gravity
    if (orderMoon > 0) then
        call UTC2TT(UTC, TT)
        call planetsPV(TT, nEarth, nMoon, LunarPV)
        r_lunar = LunarPV(1:3)
        rSat_LCI = rs(1:3) - r_lunar
        call LCI2LCF(TT+2400000.5d0, rSat_LCI, rSat_LCF)
        call MoonGravityField%acceleration(rSat_LCF * 1d3, orderMoon, orderMoon, force_moonAspherical)
        call LCF2LCI(TT+2400000.5d0, force_moonAspherical, force_moonAspherical_LCI)
        force = force + force_moonAspherical_LCI
    endif

    ! Post-Newton effect
    if (k7 == 1) then
        call post_newton_effect(rs(1:3), rs(4:6), force6)
        force = force + force6
    endif

    !  thisThruster force: a = F/m, with m the CURRENT mass, i.e. the mass
    !  published in `mass` (dry + propellant at ignition) minus everything burnt
    !  so far. `mass` and `thisThruster` are a snapshot taken when the
    !  propagation started, and the thruster ledger only moves when the caller
    !  flushes it, so m is piecewise constant over the propagation -- not a
    !  per-stage value that shrinks with each RKF78 step.
    !  A non-positive current mass means the fuel budget exceeds the satellite
    !  mass: the thrust is dropped rather than divided by zero.
    if ( thisThruster%is_active) then
        force7 =  thisThruster%get_thrust_ECI(rs(1:3), rs(4:6))
        massNow = mass - thisThruster%total_fuel_used
        if (massNow > 0d0) then
            force = force + force7/massNow
        endif
    endif

    force = force/acc_unit
    
    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    subroutine get_stw(elements, u, force, S, T, W)
    implicit none
    double precision elements(6), u, force(3), S, T, W
    double precision P(3), Q(3), rxt(3)
    double precision r_hat(3), t_hat(3), w_hat(3)
    !--------------------------------------------
    P(1) = dcos(elements(4))
    P(2) = dsin(elements(4))
    P(3) = 0d0

    Q(1) = -dsin(elements(4))*dcos(elements(3))
    Q(2) = dcos(elements(4))*dcos(elements(3))
    Q(3) = dsin(elements(3))

    rxt(1) = dsin(elements(4))*dsin(elements(3))
    rxt(2) = -dcos(elements(4))*dsin(elements(3))
    rxt(3) = dcos(elements(3))

    r_hat = dcos(u)*P + dsin(u)*Q
    t_hat = -dsin(u)*P + dcos(u)*Q
    w_hat = rxt

    call dot(3, force, r_hat, S)
    call dot(3, force, t_hat, T)
    call dot(3, force, w_hat, W)

    end subroutine
    
    subroutine get_stw_PV(x0, F, S, T, W)
    implicit none
    double precision  x0(6), F(3), S, T, W
    double precision  sUnit(3), tUnit(3), wUnit(3)
    sUnit = x0(1:3)/norm2(x0(1:3))
    call cross_3d(x0(1:3), x0(4:6), wUnit)
    wUnit = wUnit/norm2(wUnit)
    call cross_3d(wUnit, sUnit, tUnit)
    
    S = dot_product(F, sUnit)
    T = dot_product(F, tUnit)
    W = dot_product(F, wUnit)
    endsubroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    
    
    
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    endmodule