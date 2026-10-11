
    ! module_rotator，创建一个旋转器对象
    ! created by Shawn Hu in 2025/03/18
    !
    ! 旋转器定义了坐标系1在坐标系0中的姿态，使用旋转参数来定义一个
    ! 旋转器，旋转参数的类型有：1）四元素数；2）欧拉角；3）方向余弦矩阵；
    ! 旋转器功能有：1）旋转参数类型转换；2）坐标转换；
    
    module module_rotator
    implicit none
    
    public
    private :: quaternion_EulerAngles
    private :: quaternion_DCM
    private :: EulerAngles_quaternion
    private :: DCM_quaternion
    private :: EulerAngles_DCM
    private :: DCM_EulerAngles
    
    type type_rotator
        integer :: rotType = 1
        !----------------------------------------------
        ! rotType = 
        ! 1 - quaternion
        ! 2 - Euler angles
        ! 3 - direction cosine matrix, DCM
        ! 4 - Yaw Pitch Roll, YPR
        !----------------------------------------------
        double precision :: val(9) = 0d0
        double precision :: valDot(9) = 0d0
        integer :: seqEuler = 313
        !---------------------------
        ! seqEuler = 
        ! 313 - z x z
        ! 312 - z x y
        ! 321 - z y x
        !---------------------------
    contains
        procedure :: convertRotType
        procedure :: convertPos
    endtype
    
    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine convertRotType(this, toType)
    implicit none
    class(type_rotator) :: this
    integer, intent(in) :: toType
    integer :: fromType, seqEuler
    double precision :: xq(4), xe(3), xd(9)
    fromType = this%rotType
    this%rotType = toType
    ! quaternion -> Euler angles
    if (fromType == 1 .and. toType == 2) then
        xq = this%val(1:4)
        seqEuler = this%seqEuler
        call quaternion_EulerAngles(xq, seqEuler, xe)
        this%val = 0d0
        this%val(1:3) = xe(1:3)
    endif
    ! quaternion -> DCM
    if (fromType == 1 .and. toType == 3) then
        xq = this%val(1:4)
        call quaternion_DCM(xq, xd)
        this%val = 0d0
        this%val(1:9) = xd(1:9)
    endif
    ! Euler angles -> quaternion
    if (fromType == 2 .and. toType == 1) then
        xe = this%val(1:3)
        seqEuler = this%seqEuler
        call EulerAngles_quaternion(xe, seqEuler, xq)
        this%val = 0d0
        this%val(1:4) = xq(1:4)
    endif
    ! DCM -> quaternion
    if (fromType == 3 .and. toType == 1) then
        xd = this%val(1:9)
        call DCM_quaternion(xd, xq)
        this%val = 0d0
        this%val(1:4) = xq(1:4)
    endif
    ! Euler angles -> DCM
    if (fromType == 2 .and. toType == 3) then
        xe = this%val(1:3)
        seqEuler = this%seqEuler
        call EulerAngles_DCM(xe, seqEuler, xd)
        this%val = 0d0
        this%val(1:9) = xd(1:9)
    endif
    ! DCM -> Euler angles
    if (fromType == 3 .and. toType == 2) then
        xd = this%val(1:9)
        seqEuler = this%seqEuler
        call DCM_EulerAngles(xd, seqEuler, xe)
        this%val = 0d0
        this%val(1:3) = xe(1:3)
    endif
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function convertPos(this, ri) result(ro)
    implicit none
    class(type_rotator) :: this
    double precision :: ri(3)
    double precision :: ro(3)
    double precision :: DCM(3, 3), xd(9)
    if (this%rotType == 1) then ! quaternion
        call quaternion_DCM(this%val(1:4), xd)
        ! 将xd(1:9)转换为DCM(1:3, 1:3)
        DCM = transpose(reshape(xd(1:9), [3,3]))
        ro = matmul(DCM, ri)
    elseif (this%rotType == 2) then ! Euler angles
        call EulerAngles_DCM(this%val(1:3), this%seqEuler, xd)
        DCM = transpose(reshape(xd(1:9), [3,3]))
        ro = matmul(DCM, ri)
    elseif (this%rotType == 3) then ! direction cosine matrix
        xd = this%val(1:9)
        DCM = transpose(reshape(xd(1:9), [3,3]))
        ro = matmul(DCM, ri)
    elseif (this%rotType == 4) then ! yaw, pitch, roll
        
    else
        print*, "warning @convertPos"
        print*, "unknown rotType, rotType = ", this%rotType
    endif
    
    endfunction

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine quaternion_EulerAngles(xq, seqEuler, xe)
    implicit none
    double precision, intent(in) :: xq(4)
    integer, intent(in) :: seqEuler
    double precision, intent(out) :: xe(3)
    double precision :: xd(9)
    call quaternion_DCM(xq, xd)
    call DCM_EulerAngles(xd, seqEuler, xe)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 四元素数 -> 方向余弦矩阵
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine quaternion_DCM(xq, xd)
    implicit none
    double precision, intent(in) :: xq(4)
    double precision, intent(out) :: xd(9)
    double precision :: q0, q1, q2, q3, q4
    q0 = xq(1)
    q1 = xq(2)
    q2 = xq(3)
    q3 = xq(4)
    xd(1) = q0**2 + q1**2 - q2**2 - q3**2
    xd(2) = 2d0*(q1*q2 + q3*q0)
    xd(3) = 2d0*(q1*q3 - q2*q0)
    xd(4) = 2d0*(q1*q2 - q3*q0)
    xd(5) = q0**2 - q1**2 + q2**2 - q3**2
    xd(6) = 2d0*(q2*q3 + q1*q0)
    xd(7) = 2d0*(q1*q3 + q2*q0)
    xd(8) = 2d0*(q2*q3 - q1*q0)
    xd(9) = q0**2  -q1**2 - q2**2 + q3**2
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine EulerAngles_quaternion(xe, seqEuler, xq)
    implicit none
    double precision, intent(in) :: xe(3)
    integer, intent(in) :: seqEuler
    double precision, intent(out) :: xq(4)
    double precision :: xd(9)
    call EulerAngles_DCM(xe, seqEuler, xd)
    call DCM_quaternion(xd, xq)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 方向余弦矩阵 -> 四元素数
    ! Direction cosine matrix -> quaternion
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine DCM_quaternion(xd, xq)
    implicit none
    double precision, intent(in) :: xd(9)
    double precision, intent(out) :: xq(4)

    real(8) :: c(3,3)  !  3x3矩阵
    real(8) :: qs, qx, qy, qz
    real(8) :: trc1, factor
    
    c(1, 1:3) = xd(1:3)
    c(2, 1:3) = xd(4:6)
    c(3, 1:3) = xd(7:9)

    trc1 = c(1,1) + c(2,2) + c(3,3) + 1.0d0  ! C++索引转Fortran索引+1

    if (trc1 >= 0.004d0) then
        qs = sqrt(trc1) / 2.0d0
        qx = (c(2,3) - c(3,2)) / (4.0d0 * qs)  ! (0,1)->(1,2)转(2,3)
        qy = (c(3,1) - c(1,3)) / (4.0d0 * qs)
        qz = (c(1,2) - c(2,1)) / (4.0d0 * qs)
    else if (1.0d0 - c(1,1) + c(2,2) - c(3,3) >= 0.004d0) then
        ! 处理三元运算符: (c(2,0)-c(0,2)>=0 ? 1 : -1)
        factor = 1.0d0
        if (c(3,1) - c(1,3) < 0.0d0) factor = -1.0d0
        qy = sqrt(1.0d0 - c(1,1) + c(2,2) - c(3,3)) / 2.0d0 * factor
        qx = (c(2,1) + c(1,2)) / (4.0d0 * qy)
        qz = (c(3,2) + c(2,3)) / (4.0d0 * qy)
        qs = (c(3,1) - c(1,3)) / (4.0d0 * qy)
    else if (1.0d0 + c(1,1) - c(2,2) - c(3,3) >= 0.004d0) then
        factor = 1.0d0
        if (c(2,3) - c(3,2) < 0.0d0) factor = -1.0d0
        qx = sqrt(1.0d0 + c(1,1) - c(2,2) - c(3,3)) / 2.0d0 * factor
        qy = (c(2,1) + c(1,2)) / (4.0d0 * qx)
        qz = (c(1,3) + c(3,1)) / (4.0d0 * qx)
        qs = (c(2,3) - c(3,2)) / (4.0d0 * qx)
    else
        factor = 1.0d0
        if (c(1,2) - c(2,1) < 0.0d0) factor = -1.0d0
        qz = sqrt(1.0d0 - c(1,1) - c(2,2) + c(3,3)) / 2.0d0 * factor
        qx = (c(1,3) + c(3,1)) / (4.0d0 * qz)
        qy = (c(2,3) + c(3,2)) / (4.0d0 * qz)
        qs = (c(1,2) - c(2,1)) / (4.0d0 * qz)
    end if

    call Normalize(qs, qx, qy, qz)  ! 调用归一化子程序
    
    xq = (/qs, qx, qy, qz/)
contains
    ! 四元数归一化子程序
    subroutine Normalize(qs, qx, qy, qz)
        real(8), intent(inout) :: qs, qx, qy, qz
        real(8) :: norm
        norm = sqrt(qs**2 + qx**2 + qy**2 + qz**2)
        qs = qs / norm
        qx = qx / norm
        qy = qy / norm
        qz = qz / norm
    end subroutine Normalize

    endsubroutine DCM_quaternion
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine EulerAngles_DCM(xe, seqEuler, xd)
    implicit none
    double precision, intent(in) :: xe(3)
    integer, intent(in) :: seqEuler
    double precision, intent(out) :: xd(9)
    double precision :: th1, th2, th3
    double precision :: C1, C2, C3
    double precision :: S1, S2, S3
    if (seqEuler == 313) then
        th1 = xe(1)
        th2 = xe(2)
        th3 = xe(3)
        S1 = dsin(th1)
        S2 = dsin(th2)
        S3 = dsin(th3)
        C1 = dcos(th1)
        C2 = dcos(th2)
        C3 = dcos(th3)
        xd(1) = C1*C3 - S1*C2*S3
        xd(2) = S1*C3 + C1*C2*S3
        xd(3) = S2*S3
        xd(4) = -C1*S3 - S1*C2*C3
        xd(5) = -S1*S3 + C1*C2*C3
        xd(6) = S2*C3
        xd(7) = S1*S2
        xd(8) = - C1*S2
        xd(9) = C2
    else
        print*, "error@EulerAngles_DCM"
        print*, "Sequence of Euler angles can only be z-x-z."
        stop
    endif
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !subroutine DCM_EulerAngles(xd, seqEuler, xe)
    !implicit none
    !double precision, intent(in) :: xd(9)
    !integer, intent(in) :: seqEuler
    !double precision, intent(out) :: xe(3)
    !double precision :: xdp(9)
    !if(seqEuler == 313) then
    !
    !    xdp = xd
    !    xe(1) = - datan(xdp(7)/xdp(8))
    !    xe(2) = dacos(xdp(9))
    !    xe(3) = datan(xdp(3)/xdp(6))
    !else
    !    print*, "error@DCM_EulerAngles"
    !    print*, "Sequence of Euler angles can only be z-x-z."
    !    stop
    !endif
    !
    !endsubroutine
    
    subroutine DCM_EulerAngles(xd, seqEuler, xe)
    implicit none
    double precision, intent(in) :: xd(9)
    integer, intent(in) :: seqEuler
    double precision, intent(out) :: xe(3)
    double precision :: xdp(9), pi, beta_threshold
    parameter (pi = 3.141592653589793d0, beta_threshold = 1.0d-6)

    if (seqEuler == 313) then
        xdp = xd
        xe(2) = dacos(xdp(9))  ! 计算俯仰角β

        ! 奇异性检测与处理
        if (abs(xe(2) - pi/2.0d0) < beta_threshold .or. abs(xe(2) + pi/2.0d0) < beta_threshold) then
            xe(1) = 0.0d0  ! 默认航向角
            xe(3) = datan2(xdp(3), xdp(6))  ! 横滚角
        else
            xe(1) = datan2(xdp(7), - xdp(8))  ! 航向角
            xe(3) = datan2(xdp(3), xdp(6))   ! 横滚角
        endif
    else
        print*, "error@DCM_EulerAngles"
        print*, "Sequence of Euler angles can only be z-x-z."
        stop
    endif
    endsubroutine
    
    endmodule
    