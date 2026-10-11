module module_shadow
  implicit none
  private
  public :: inShadow, SUN_RADIUS

  ! 太阳半径常数（单位必须与位置向量一致，此处默认使用 km）
  real(8), parameter :: SUN_RADIUS = 696340.0d0

contains

  function inShadow(rSun, rSat, rPlanet, radiusPlanet) result(flag)
    ! 输入：
    !   rSun        : 惯性系下太阳位置向量 (3), km
    !   rSat        : 惯性系下卫星位置向量 (3), km
    !   rPlanet     : 惯性系下行星中心位置向量 (3), km
    !   radiusPlanet: 行星半径, km
    ! 输出：
    !   flag = 1 : 在本影中
    !   flag = 2 : 在半影中
    !   flag = 0 : 不在影子中
    real(8), intent(in) :: rSun(3), rSat(3), rPlanet(3)
    real(8), intent(in) :: radiusPlanet
    integer             :: flag

    real(8) :: s(3), e(3), r(3)
    real(8) :: d, proj, x, rho_sq, rho
    real(8) :: u_shadow, p_shadow

    ! 行星指向太阳的向量
    s = rSun - rPlanet
    d = sqrt(dot_product(s, s))

    ! 太阳与行星重合，无阴影
    if (d <= 0.0d0) then
      flag = 0
      return
    end if

    ! 单位向量：行星 -> 太阳
    e = s / d

    ! 行星指向卫星的向量
    r = rSat - rPlanet

    ! 卫星在行星-太阳连线上的投影（沿 e 方向）
    proj = dot_product(r, e)

    ! x 为卫星沿背离太阳方向的投影距离
    x = -proj

    ! 卫星位于行星朝向太阳的一侧，不可能在阴影中
    if (x <= 0.0d0) then
      flag = 0
      return
    end if

    ! 卫星到行星-太阳连线的垂直距离
    rho_sq = dot_product(r, r) - proj**2
    rho = sqrt(max(rho_sq, 0.0d0))

    ! 距离行星中心 x 处的本影半径
    u_shadow = max(0.0d0, radiusPlanet - (SUN_RADIUS - radiusPlanet) / d * x)

    ! 距离行星中心 x 处的半影外边界半径
    p_shadow = radiusPlanet + (SUN_RADIUS + radiusPlanet) / d * x

    ! 本影优先
    if (rho <= u_shadow) then
      flag = 1
    else if (rho <= p_shadow) then
      flag = 2
    else
      flag = 0
    end if

  end function inShadow

end module module_shadow