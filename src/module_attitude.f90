module module_attitude
    use module_time
    use module_rotator
    implicit none
    
    public
    
    private :: coordinateTransform_J20002COCS
    
    type :: type_attitude
        integer :: satNum
        integer :: mode = 1 ! attitude control mode
        !----------------------------------------------
        ! mode = 
        ! 1 - 对地定向：X -> Nadir aligned, X -> velocity constrained;
        ! 2 - 空间静止：inertially fixed;
        ! 3 - 对日定向：X -> Sun aligned, Z -> nadir constrained;
        !----------------------------------------------
        integer :: coordinate = 1 ! coordinate frame
        !----------------------------------------------
        !coordinate = 
        ! 1 - Earth Center Inertial, J2000
        ! 2 - Earth Center Fixed, ECF
        ! 3 - Centroid Orbit Coordinate System, COCS
        ! 4 - Spacecraft Body Coordinate System, SBCS
        !----------------------------------------------
        integer :: ntrajectory = 0
        integer :: attiType = 1 ! data type of the attitude coefficients
        !----------------------------------------------
        ! attiType = 
        ! 1 - quaternion
        ! 2 - Euler angles
        ! 3 - direction cosine matrix, DCM
        ! 4 - Yaw Pitch Roll, YPR
        !----------------------------------------------
        type(type_t) :: t0
        type(type_t), allocatable :: timeTable(:)
        type(type_rotator) :: atti0
        type(type_rotator), allocatable :: atti(:)
    contains
        procedure :: initialize => init_attitude
        procedure :: free => free_attitude
        procedure :: Convert => convert_attitude_type
        procedure :: readAtti => read_attitude_file
        procedure :: saveAtti => save_attitude_file
        procedure :: predict => predict_attitude
        procedure :: calAtti
        procedure :: coorTrans => coordinate_transform
    end type type_attitude

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 初始化type_attitude对象
    subroutine init_attitude(this, mode, ntrajectory, attiType, timeTable)
    class(type_attitude), intent(inout) :: this
    integer, intent(in) :: mode, ntrajectory, attiType
    type(type_t), intent(in) :: timeTable(:)
    integer :: npe
    !--------------------------------------------------------------
    if (attiType < 1 .or. attiType > 4) then
        write(*, *) " >_< warning@ini_attitude"
        write(*, *) "unknown attiType, attiType = ", attiType
    endif
    this%mode = mode
    this%ntrajectory = ntrajectory
    this%attiType = attiType
    select case(attiType)
    case(1)
        npe = 4
    case(2)
        npe = 3
    case(3)
        npe = 9
    case(4)
        npe = 3
    endselect
    this%timeTable = timeTable
    allocate(this%atti(ntrajectory))
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 释放内存
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine free_attitude(this)
    class(type_attitude), intent(inout) :: this
    if(allocated(this%timeTable)) deallocate(this%timeTable)
    if(allocated(this%atti)) deallocate(this%atti)
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 姿态参数类型转换
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine convert_attitude_type(this, toType)
    class(type_attitude), intent(inout) :: this
    integer, intent(in) :: toType
    integer :: fromType, ntrajectory, i
    
    if (this%attiType == toType) return
    
    ntrajectory = this%ntrajectory
    do i = 1, ntrajectory
        call this%atti(i)%convertRotType(toType)
    enddo
    this%attiType = toType
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 读取文本文件
    subroutine read_attitude_file(this, filename)
    class(type_attitude), intent(inout) :: this
    character(len=*), intent(in) :: filename
    integer :: i, unit
    open(newunit=unit, file=filename, status='old')
    do i = 1, this%ntrajectory
        !read(unit, *) this%timeTable(i)%sec, this%atti(i,:)
    end do
    close(unit)
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 保存到文本文件
    subroutine save_attitude_file(this, filename)
    class(type_attitude), intent(in) :: this
    character(len=*), intent(in) :: filename
    integer :: i, unit
    open(newunit=unit, file=filename, status='replace')
    do i = 1, this%ntrajectory
        !write(unit, '(I4,1x,I3,4F12.6)') this%timeTable(i)%sec, this%atti(i,:)
    end do
    close(unit)
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 姿态预报
    subroutine predict_attitude(this)
    class(type_attitude), intent(inout) :: this
    
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calAtti(this, ntrajectory, trajectory, suntrajectory)
    implicit none
    class(type_attitude) :: this
    integer, intent(in), optional :: ntrajectory
    double precision, intent(in), optional :: trajectory(:,:)
    double precision, intent(in), optional :: suntrajectory(:,:)
    
    ! 根据this%mode的值选择调用的子程序
    select case (this%mode)
        case (1)
            ! 检查必要参数是否存在
            if (.not. present(ntrajectory) .or. .not. present(trajectory)) then
                error stop "calAtti: ntrajectory and trajectory required for mode 1"
            endif
            call calAtti_1(this, ntrajectory, trajectory)
            
        case (2)
            call calAtti_2(this)
            
        case (3)
            ! 检查必要参数是否存在
            if (.not. present(ntrajectory) .or. .not. present(trajectory) .or. .not. present(suntrajectory)) then
                error stop "calAtti: ntrajectory, trajectory, and suntrajectory required for mode 3"
            endif
            call calAtti_3(this, ntrajectory, trajectory, suntrajectory)
            
        case default
            error stop "calAtti: invalid mode"
    end select
    
    endsubroutine calAtti
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 计算航天器的姿态，姿态模式为对地定向，即X轴指向地心，Y轴受速度方向限制；
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calAtti_1(this, ntrajectory, trajectory) 
    implicit none
    class(type_attitude) :: this
    integer, intent(in) :: ntrajectory ! 星历点数
    double precision, intent(in) :: trajectory(ntrajectory, 6) ! 星历，位置速度
    integer :: i
    double precision :: r(3), v(3), r_unit(3), x_body(3), h(3), z_body(3), y_body(3)
    double precision :: r_norm, h_norm, cross_prod(3)
    double precision, parameter :: epsilon = 1.0d-12
    double precision :: DCM(9)
    
    do i = 1, ntrajectory
        ! 提取位置和速度
        r = trajectory(i, 1:3)
        v = trajectory(i, 4:6)
    
        ! 计算位置单位向量和X轴指向
        r_norm = norm2(r)
        if (r_norm < epsilon) error stop "error@calAtti_1, r_norm = 0."
        r_unit = r / r_norm
        x_body = -r_unit
    
        ! 计算轨道角动量方向(h)
        cross_prod(1) = r(2)*v(3) - r(3)*v(2)
        cross_prod(2) = r(3)*v(1) - r(1)*v(3)
        cross_prod(3) = r(1)*v(2) - r(2)*v(1)
        h = cross_prod
        h_norm = norm2(h)

        ! 计算Z轴方向（轨道面法向）
        z_body = h / h_norm
    
        ! 计算Y轴方向（X×Z保证正交）
        cross_prod(1) = x_body(2)*z_body(3) - x_body(3)*z_body(2)
        cross_prod(2) = x_body(3)*z_body(1) - x_body(1)*z_body(3)
        cross_prod(3) = x_body(1)*z_body(2) - x_body(2)*z_body(1)
        y_body = cross_prod
    
        ! 存储方向余弦矩阵
        DCM(1:3) = x_body
        DCM(4:6) = y_body
        DCM(7:9) = z_body
        
        this%atti(i)%rotType = 3
        this%atti(i)%val(1:9) = DCM
        if (this%attiType /= 3) then
            call this%atti(i)%convertRotType(this%attiType)
        endif

    enddo

    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 空间固定姿态
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calAtti_2(this)
    implicit none
    class(type_attitude) :: this
    integer :: i
    do i = 1, this%ntrajectory
        this%atti(i) = this%atti0
    enddo
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 对日定向姿态模式，X轴指向太阳，Y轴受地心方向限制
    ! inputs:
    !   ntrajectory    星历点数
    !   sattrajectory  航天器位置速度星历
    !   suntrajectory  太阳位置星历
    ! inOut:
    !   this    type_attitude对象
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calAtti_3(this, ntrajectory, sattrajectory, suntrajectory)
    implicit none
    class(type_attitude) :: this
    integer, intent(in) :: ntrajectory
    double precision, intent(in) :: sattrajectory(ntrajectory, 6), suntrajectory(ntrajectory, 3)
    ! 声明临时变量
    real(8) :: sunVec(3), earthVec(3), xAxis(3), yAxis(3), zAxis(3)
    integer :: i
    double precision :: DCM(9)

    do i = 1, ntrajectory
        ! 计算太阳矢量（航天器质心指向太阳）
        sunVec = suntrajectory(i, 1:3) - sattrajectory(i, 1:3)
        sunVec = sunVec / norm2(sunVec)  ! 归一化

        ! 计算地球矢量（航天器质心指向地心）
        earthVec = -sattrajectory(i, 1:3)   ! 地心在J2000坐标系原点时
        earthVec = earthVec / norm2(earthVec)
        
        ! X轴指向太阳（主约束）
        xAxis = sunVec

        ! Y轴约束：垂直于X轴和地球矢量（次约束）
        yAxis = cross_product(earthVec, xAxis)
        yAxis = yAxis / norm2(yAxis)

        ! Z轴通过右手法则确定
        zAxis = cross_product(xAxis, yAxis)
        
        ! 构建DCM矩阵
        DCM(1:3) = xAxis
        DCM(4:6) = yAxis
        DCM(7:9) = zAxis
        
        this%atti(i)%rotType = 3
        this%atti(i)%val(1:9) = DCM
        if (this%attiType /= 3) then
            call this%atti(i)%convertRotType(this%attiType)
        endif

    enddo
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function cross_product(a, b) result(c)
        real(8), intent(in) :: a(3), b(3)
        real(8) :: c(3)
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 坐标系转换
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine coordinate_transform(this, toCoor, trajectory)
        implicit none
        class(type_attitude), intent(inout) :: this
        integer, intent(in) :: toCoor
        double precision, intent(in), optional :: trajectory(:, :)
        integer :: fromCoor, attiType0, i
        logical :: needtrajectory

        ! 检查目标坐标系有效性
        if (toCoor < 1 .or. toCoor > 4) then
            error stop "coordinate_transform: Invalid toCoor value."
        end if

        ! 保存当前姿态类型并转换为DCM
        attiType0 = this%attiType
        if (attiType0 /= 3) call this%convert(3)

        fromCoor = this%coordinate
        needtrajectory = .false.

        ! 处理J2000到COCS的转换
        if (fromCoor == 1 .and. toCoor == 3) then
            if (.not. present(trajectory)) then
                error stop "coordinate_transform: trajectory required for J2000->COCS"
            end if
            call coordinateTransform_J20002COCS(this%ntrajectory, this%atti, trajectory)
            this%coordinate = 3

        ! 添加其他转换情况（示例）
        ! else if (其他转换条件) then
        !     ... 

        else
            error stop "coordinate_transform: Unsupported coordinate conversion"
        end if

        ! 转换回原姿态类型
        if (attiType0 /= 3) call this%convert(attiType0)

    end subroutine coordinate_transform
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! J2000到COCS坐标系转换
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine coordinateTransform_J20002COCS(ntrajectory, atti, trajectory)
        implicit none
        integer, intent(in) :: ntrajectory
        type(type_rotator), intent(inout) :: atti(ntrajectory)
        double precision, intent(in) :: trajectory(ntrajectory, 6)
        type(type_attitude) :: a1
        integer :: i
        double precision :: R_CJ(3,3), R_J_to_C(3,3), DCM_old(3,3), DCM_new(3,3)

        ! 初始化临时姿态对象
        call a1%initialize(mode=1, ntrajectory=ntrajectory, attiType=3, timeTable=[(type_t(0.0d0), i=1,ntrajectory)])
        call a1%calAtti(ntrajectory=ntrajectory, trajectory=trajectory)  ! 生成COCS到J2000的DCM

        do i = 1, ntrajectory
            ! 获取COCS到J2000的旋转矩阵
            R_CJ = reshape(a1%atti(i)%val(1:9), [3,3])
            R_J_to_C = transpose(R_CJ)  ! 转置得到J2000到COCS的旋转

            ! 获取当前姿态的DCM并转换为矩阵形式
            DCM_old = reshape(atti(i)%val(1:9), [3,3])

            ! 应用坐标系转换: DCM_new = R_J_to_C * DCM_old
            DCM_new = matmul(R_J_to_C, DCM_old)

            ! 更新姿态数据（按列优先展开）
            atti(i)%val(1:3) = DCM_new(:,1)
            atti(i)%val(4:6) = DCM_new(:,2)
            atti(i)%val(7:9) = DCM_new(:,3)
        end do

        call a1%free()  ! 清理临时对象

    end subroutine coordinateTransform_J20002COCS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end module module_attitude
