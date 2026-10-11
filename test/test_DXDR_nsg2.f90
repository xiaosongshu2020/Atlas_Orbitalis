!===============================================================================
! test_DXDR_nsg2.f90
!
! 验证 calDRDX_NSG2 (dcts/dnsg2) 与 calDXDR_nsg2 (dnsg2/dcts)
! 互为逆矩阵，即 matmul(DRDX, DXDR) 应等于单位矩阵。
!
! Author: Shawn Hu
! Date: 2026-06-01
!===============================================================================

program test_DXDR_nsg2
    use module_const
    use module_elements
    use module_ODor
    implicit none

    double precision :: kpl(6), nsg2(6), cts(6)
    double precision :: DRDX(6, 6), DXDR(6, 6)
    double precision :: product1(6, 6), product2(6, 6), identity(6, 6)
    double precision :: max_err1, max_err2
    integer :: i, j, imax1, jmax1, imax2, jmax2

    !===========================================================================
    ! 步骤1: 初始化 Keplerian 轨道根数
    !        a = 7000 km, e = 0.01, i = 45°, Omega = 30°, omega = 60°, M = 90°
    !===========================================================================
    kpl(1) = 7000.0d0                           ! 半长轴 a (km)
    kpl(2) = 0.01d0                             ! 偏心率 e
    kpl(3) = 45.0d0 * deg                       ! 倾角 i (rad)
    kpl(4) = 30.0d0 * deg                       ! 升交点赤经 Omega (rad)
    kpl(5) = 60.0d0 * deg                       ! 近地点幅角 omega (rad)
    kpl(6) = 90.0d0 * deg                       ! 平近点角 M (rad)

    write(*, '(A)') ""
    write(*, '(A)') "=== Initial Keplerian Elements ==="
    write(*, '(A, F10.3)') "  a        (km) = ", kpl(1)
    write(*, '(A, F10.6)') "  e             = ", kpl(2)
    write(*, '(A, F10.4)') "  i        (deg)= ", kpl(3) / deg
    write(*, '(A, F10.4)') "  Omega    (deg)= ", kpl(4) / deg
    write(*, '(A, F10.4)') "  omega    (deg)= ", kpl(5) / deg
    write(*, '(A, F10.4)') "  M        (deg)= ", kpl(6) / deg

    !===========================================================================
    ! 步骤2: 转换为 nsg2 和 cts 形式
    !===========================================================================
    call kpl2nsg2(kpl, nsg2)
    call nsg22cts(nsg2, cts)

    write(*, '(A)') ""
    write(*, '(A)') "=== Non-Singular Elements Type 2 (nsg2) ==="
    write(*, '(A, ES20.12)') "  a             = ", nsg2(1)
    write(*, '(A, ES20.12)') "  xi            = ", nsg2(2)
    write(*, '(A, ES20.12)') "  eta           = ", nsg2(3)
    write(*, '(A, ES20.12)') "  h             = ", nsg2(4)
    write(*, '(A, ES20.12)') "  k             = ", nsg2(5)
    write(*, '(A, ES20.12)') "  lambda (rad)  = ", nsg2(6)

    write(*, '(A)') ""
    write(*, '(A)') "=== Cartesian Coordinates (cts) ==="
    write(*, '(A, 3ES20.12)') "  r  (km)  = ", cts(1:3)
    write(*, '(A, 3ES20.12)') "  v  (km/s)= ", cts(4:6)

    !===========================================================================
    ! 步骤3: 计算两个雅可比矩阵
    !===========================================================================
    call calDRDX_NSG2(nsg2, DRDX)      ! DRDX = d(cts)/d(nsg2)
    call calDXDR_nsg2(cts, DXDR)       ! DXDR = d(nsg2)/d(cts)

    write(*, '(A)') ""
    write(*, '(A)') "=== Jacobian DRDX = d(cts)/d(nsg2) ==="
    do i = 1, 6
        write(*, '(6ES16.8)') (DRDX(i, j), j = 1, 6)
    end do

    write(*, '(A)') ""
    write(*, '(A)') "=== Jacobian DXDR = d(nsg2)/d(cts) ==="
    do i = 1, 6
        write(*, '(6ES16.8)') (DXDR(i, j), j = 1, 6)
    end do

    !===========================================================================
    ! 步骤4: 计算乘积并验证是否为单位矩阵
    !===========================================================================
    product1 = matmul(DRDX, DXDR)      ! dcts/dnsg2 * dnsg2/dcts = ?
    product2 = matmul(DXDR, DRDX)      ! dnsg2/dcts * dcts/dnsg2 = ?

    ! 构造单位矩阵用于比较
    identity = 0.0d0
    do i = 1, 6
        identity(i, i) = 1.0d0
    end do

    write(*, '(A)') ""
    write(*, '(A)') "=== product1 = DRDX * DXDR (expect identity) ==="
    do i = 1, 6
        write(*, '(6ES16.8)') (product1(i, j), j = 1, 6)
    end do

    write(*, '(A)') ""
    write(*, '(A)') "=== product2 = DXDR * DRDX (expect identity) ==="
    do i = 1, 6
        write(*, '(6ES16.8)') (product2(i, j), j = 1, 6)
    end do

    !===========================================================================
    ! 步骤5: 计算与单位矩阵的最大偏差
    !===========================================================================
    max_err1 = 0.0d0
    imax1 = 1
    jmax1 = 1
    do i = 1, 6
        do j = 1, 6
            if (abs(product1(i, j) - identity(i, j)) > max_err1) then
                max_err1 = abs(product1(i, j) - identity(i, j))
                imax1 = i
                jmax1 = j
            end if
        end do
    end do

    max_err2 = 0.0d0
    imax2 = 1
    jmax2 = 1
    do i = 1, 6
        do j = 1, 6
            if (abs(product2(i, j) - identity(i, j)) > max_err2) then
                max_err2 = abs(product2(i, j) - identity(i, j))
                imax2 = i
                jmax2 = j
            end if
        end do
    end do

    write(*, '(A)') ""
    write(*, '(A)') "=============================================================="
    write(*, '(A)') "  Verification Result"
    write(*, '(A)') "=============================================================="
    write(*, '(A, ES12.4, A, I0, A, I0)') &
        "  product1 (DRDX*DXDR) max error: ", max_err1, &
        "  at (", imax1, ", ", jmax1, ")"
    write(*, '(A, ES12.4, A, I0, A, I0)') &
        "  product2 (DXDR*DRDX) max error: ", max_err2, &
        "  at (", imax2, ", ", jmax2, ")"

    if (max_err1 < 1.0d-10 .and. max_err2 < 1.0d-10) then
        write(*, '(A)') ""
        write(*, '(A)') "  PASS: The two Jacobians are inverses of each other."
    else
        write(*, '(A)') ""
        write(*, '(A)') "  FAIL: The product is not the identity matrix!"
    end if

    write(*, '(A)') ""

end program test_DXDR_nsg2