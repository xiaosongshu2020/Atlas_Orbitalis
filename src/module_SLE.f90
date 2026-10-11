
    ! a module to define the "system of linear equations (SLE)"
    
    module module_SLE
    use module_matrix
    
    type :: type_SLE
        integer :: n ! the number of variables in SLE
        double precision, allocatable :: x(:, :) ! solution vector
        double precision, allocatable :: P(:, :) ! covariance matrix, P is the inverse matrix of the BTB
        double precision, allocatable :: BTB(:, :) ! coefficient matrix
        double precision, allocatable :: BTY(:, :) ! right-hand side (RHS)
    contains
        procedure initialize ! initialize, allocate memory
        procedure freeSLE ! free memory
        procedure addData ! add data to B and Y
        procedure solve ! solve the equation : B*X = Y; get X and P
        procedure propagate ! propagate X and P with state transition matrix (STM)
        procedure dataFusion  ! update X with new data, get new X and P, in other words, it's Kalman Filter (KF)
        procedure x2zero ! set X to 0
        procedure get_BTB_BTY ! calculate BTB and BTY from x and P
    end type
    
    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine initialize(SLE, n)
    implicit none
    class(type_SLE) :: SLE
    integer, intent(in) :: n
    !-------------------------------------------
    ! ��鲢�����ڴ�
    SLE%n = n
    if (allocated(SLE%x)) deallocate(SLE%x)
    if (allocated(SLE%BTB)) deallocate(SLE%BTB)
    if (allocated(SLE%BTY)) deallocate(SLE%BTY)
    if (allocated(SLE%P)) deallocate(SLE%P)
    
    ! �����ڴ�
    allocate(SLE%x(n, 1))  ! x is a column vector
    allocate(SLE%BTB(n, n))  ! B is a square matrix
    allocate(SLE%BTY(n, 1))  ! Y is a column vector
    allocate(SLE%P(n, n))  ! P is a square matrix
    
    ! ��ʼ��Ϊ0
    SLE%x = 0.0d0
    SLE%BTB = 0.0d0
    SLE%BTY = 0.0d0
    SLE%P = 0.0d0
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine freeSLE(SLE)
    implicit none
    class(type_SLE) :: SLE
    ! �ͷ��ڴ�ǰ����Ƿ��ѷ���
    if (allocated(SLE%x)) deallocate(SLE%x)
    if (allocated(SLE%BTB)) deallocate(SLE%BTB)
    if (allocated(SLE%BTY)) deallocate(SLE%BTY)
    if (allocated(SLE%P)) deallocate(SLE%P)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine addData(SLE, m, n, B, Y, W)
    implicit none
    class(type_SLE) :: SLE
    integer, intent(in) :: n, m
    double precision, intent(in) :: B(m, n)
    double precision, intent(in) :: Y(m, 1)
    double precision, intent(in), optional :: W(m, m)
    double precision :: Vi(m, m)
    !------------------------------------------
    if (n /= SLE%n) then
        print*, "error @subroutine addData."
        print*, "The number of columns in B should equal the dimension of x."
        return
    endif
    ! W is assumed to be the identity matrix by default
    if (.not. present(W)) then
        call IMatrix(m, Vi)
    else
        Vi = W
    endif
    SLE%BTB = SLE%BTB + matmul(transpose(B), matmul(Vi, B))
    SLE%BTY = SLE%BTY + matmul(transpose(B), matmul(Vi, Y))
    
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine solve(SLE)
    implicit none
    class(type_SLE) :: SLE
    double precision, allocatable :: d(:), M(:, :), r(:, :)
    integer :: i, j, n

    n = SLE%n
    allocate(d(n), M(n, n), r(n, 1))

    ! 先做对角均衡(equilibration)再求逆。待估参数的量纲差别很大(例如轨道状态
    ! 量的位置 ~1e0 LU 与速度 ~1e-1 LU/TU, 而它们对测距的灵敏度又相差几十倍),
    ! BTB 的对角元会跨越好几个量级, 直接求逆时弱观测方向的全部有效位都被舍掉,
    ! 迭代表现为"残差每轮只降百分之十几"的几何收敛。
    !   d_i = 1/sqrt(BTB_ii);  (D·BTB·D) y = D·BTY,  x = D·y
    ! 数学上与直接求解等价, 只改数值性质。
    do i = 1, n
        if (SLE%BTB(i, i) > 0d0) then
            d(i) = 1d0/dsqrt(SLE%BTB(i, i))
        else
            d(i) = 1d0
        endif
    enddo
    do j = 1, n
        do i = 1, n
            M(i, j) = SLE%BTB(i, j)*d(i)*d(j)
        enddo
        r(j, 1) = SLE%BTY(j, 1)*d(j)
    enddo

    call InvMatrix_LUsolver(n, M, SLE%P)      ! P = inv(D·BTB·D)
    SLE%x = matmul(SLE%P, r)
    SLE%x(1:n, 1) = SLE%x(1:n, 1)*d(1:n)      ! x = D·y

    do j = 1, n                               ! P = D·inv(D·BTB·D)·D
        do i = 1, n
            SLE%P(i, j) = SLE%P(i, j)*d(i)*d(j)
        enddo
    enddo

    ! BTB and BTY will be set to zeros as the information has been saved in x and P
    SLE%BTB = 0d0
    SLE%BTY = 0d0

    deallocate(d, M, r)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine propagate(SLE, n, STM, intPNTM)
    implicit none
    class(type_SLE) :: SLE
    integer, intent(in) :: n
    double precision, intent(in) :: STM(n, n)
    double precision, intent(in), optional :: intPNTM(n, n)
    double precision :: x1(n, 1), P1(n, n)
    !------------------------------------
    x1 = matmul(STM, SLE%x)
    if (.not. present(intPNTM)) then
        P1 = matmul(STM, matmul(SLE%P, transpose(STM)))
    else
        P1 = matmul(STM, matmul(SLE%P, transpose(STM))) + intPNTM
    endif
    SLE%x = x1
    SLE%P = P1
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine x2zero(SLE)
    implicit none
    class(type_SLE) :: SLE
    SLE%x = 0d0
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dataFusion(SLE, n, m, B, Y, W)
    implicit none
    class(type_SLE) :: SLE
    integer, intent(in) :: n, m
    double precision, intent(in) :: B(n, m), Y(n, 1)
    double precision, intent(in), optional :: W(n, n)
    double precision :: x(m, 1), P(m, m)
    double precision :: IntPNTM(m, m)
    !------------------------------------
    if (all(SLE%P == 0d0)) then
        print*, "warning @subroutine dataFusion."
        print*, "it seems that x and P are not exist, use procedure 'solve' will be better."
        print*, "procedure 'solve' will be used automatically"
        call SLE%solve
    endif
    
    ! if n < m, use kalman filter
    IntPNTM = 0d0
    call KalmanFilter(m, n, Y, SLE%P, IntPNTM, B, W, x, P)    
    SLE%x = x
    SLE%P = P
    ! SLE%BTB and BTY is invalid now, if you want to use them, use subroutine SLE%get_BTB_BTY()
    SLE%BTB = 0d0
    SLE%BTY = 0d0
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine KalmanFilter(dimState,dimOBS,RES0,apCovMat,IntPNTM, OBSMat,OBSVAR, &
                                &postStateCorr,postCovMat)
    ! -------------------------------------------------------------------------------------------------------------
    !   Introduction: 
    !       This subroutine solves the state correction with the extended Kalman filter.
    !   Input:
    !       dimState    integer                         dimension of the state vector
    !       dimOBS      integer                         dimension of the observation
    !       RES0        double precision (dimOBS, 1)             observation residual
    !       apCovMat    double precision (dimState,dimState)      a prior covariance matrix
    !       OBSMat      double precision (dimOBS,dimState)        dependence of observation on the states
    !       OBSVAR      double precision (dimOBS,dimOBS)          observation error variance
    !   Output:
    !       postStateCorr   double precision (dimState)           a posteriori correction to the states
    !       postCovMat      double precision (dimState,dimState)  a posteriori covariance matrix of current epoch
    !   Note:
    !       This revision is intended for the inaccuracy during matrix multiplication, that costs various matrices their symmetries.
    !       The problem is alleviated by frequently applying Cholesky decomposition wherever possible to formally guarantee symmetry.
    !       Among various adaptions, the Kalman gain is broken into parts where Cholesky decomposition can be used.
    !       To safely test this implementation, this new revision is forked.
    ! -------------------------------------------------------------------------------------------------------------
    implicit none
    ! input
    integer,intent(in):: dimState,dimOBS
    double precision,intent(in):: RES0(dimOBS,1),apCovMat(dimState,dimState)
    double precision, intent(in) :: IntPNTM(dimState, dimState)
    double precision,intent(in):: OBSMat(dimOBS,dimState)
    double precision,intent(in):: OBSVAR(dimOBS,dimOBS)
    ! output
    double precision,intent(out):: postStateCorr(dimState),postCovMat(dimState,dimState)
    ! local
    double precision:: P(dimState, dimState)
    integer:: i,j
    double precision:: FPFT(dimState,dimState),KalmanGain(dimState,dimOBS),postStateCorr0(dimState,1)
    double precision:: PHT(dimState,dimOBS),HPHT(dimOBS,dimOBS),HPHTR(dimOBS,dimOBS),HPHTRinv(dimOBS,dimOBS)
    double precision:: PLo(dimState,dimState),FPLo(dimState,dimState)
    double precision:: HPLo(dimOBS,dimState),s,GUp(dimOBS,dimOBS),GUpHP(dimOBS,dimState)
    double precision:: IMat(dimState,dimState),KH(dimState,dimState),IKH(dimState,dimState),OBSVARinv(dimOBS,dimOBS)
    double precision:: apPinv(dimState,dimState),Rinv(dimOBS,dimOBS),poPinv(dimState,dimState)

    postStateCorr = 0d0
    postCovMat = 0d0
    s = 1d0
    
    ! apCovMat is a prior covariance matrix
    P = apCovMat + IntPNTM
    
    
    ! HPHT = B*P*BT, HPHTR = HPHT + WT
    HPHT  = matmul(OBSMat,matmul(P,transpose(OBSMat)))
    HPHTR = HPHT + OBSVAR

    ! >> do practical things
    if (dimOBS.eq.1) then
        if (HPHTR(1,1).eq.0d0) then
            !call RaiseError(1001,subName)
            return
        endif
        KalmanGain = matmul(P,transpose(OBSMat))/HPHTR(1,1)
        ! >> measurement update
        call IMatrix(dimState,IMat)
        KH = matmul(KalmanGain,OBSMat)
        IKH = IMat - KH
        postStateCorr0 = matmul(KalmanGain,RES0)
        postCovMat = matmul(IKH,P)
        postStateCorr(1:dimState) = postStateCorr0(1:dimState,1)
    else
        s = matTrace(HPHTR)
        call InvMatrix_LUsolver(dimOBS,HPHTR/s,HPHTRinv)
        !call matrix_inverse_2(dimOBS, HPHTR, HPHTRinv)
        HPHTRinv = HPHTRinv/s
        ! >> measurement update
        ! compute Kalman gain using a priori covariance matrix
        KalmanGain = matmul(matmul(P,transpose(OBSMat)),HPHTRinv)
        postStateCorr0 = matmul(KalmanGain,RES0)
        postStateCorr(1:dimState) = postStateCorr0(1:dimState,1)

        call IMatrix(dimState,IMat)
        KH = matmul(KalmanGain,OBSMat)
        IKH = IMat - KH
        postCovMat = matmul(IKH,P)

    endif
    
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_BTB_BTY(SLE)
    implicit none
    class(type_SLE) :: SLE
    call InvMatrix_LUsolver(SLE%n, SLE%P, SLE%BTB)
    SLE%BTY = matmul(SLE%BTB, SLE%x)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function matTrace(mat) result(tr)
    implicit none
    ! input
    double precision,intent(in):: mat(:,:)
    ! output
    double precision:: tr
    ! local
    integer:: n,j
    n = size(mat,dim=1)
    tr = 0d0
    do j = 1,n
        tr = tr+mat(j,j)
    enddo
    return
    end function matTrace
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine IMatrix(n,x)
    implicit none
    ! input
    integer,intent(in):: n
    ! input/output
    double precision,intent(out):: x(n,n)
    ! local
    integer:: i
    x = 0d0
    do i = 1,n
        x(i,i) = 1d0
    enddo
    return
    end subroutine IMatrix
    endmodule