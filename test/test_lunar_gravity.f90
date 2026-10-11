program test_lunar_gravity
   use module_Gravity
   implicit none
   
   type(GravityField) :: moon_gravity
   real*8 :: position(3)
   real*8 :: potential
   real*8 :: acceleration(3)
   integer :: degree, order
   integer :: MAX_degREE = 20
   integer :: MAX_ORDER = 20
   
   ! 初始化测试
   write(*,*) "=== Lunar Gravity Module Test ==="
   
   ! 初始化月球引力场对象
   call moon_gravity%init("Moon Gravity Field", "Moon", MAX_degREE, MAX_ORDER, 'data/AtlasData/grail.txt')

   ! 显示引力场信息
   call moon_gravity%show()
   
   ! 测试位置：月球表面上方100km处
   position = (/1838d3, 0.0d0, 0.0d0/)  ! [m]
   
   ! 设置计算阶数
   degree = 10
   order = 10
   
   ! 计算引力加速度
   call moon_gravity%acceleration(position, degree, order, acceleration)

   ! 输出结果
   write(*,*) "=== Calculation Results ==="
   write(*,*) "Position (m): ", position
   write(*,*) "Acceleration (m/s^2): ", acceleration
   write(*,*) "Acceleration magnitude (m/s^2): ", sqrt(dot_product(acceleration, acceleration))
   
   ! 计算标准月球引力加速度进行比较
   write(*,*) "Standard lunar gravity (m/s^2): ", moon_gravity%GM / (norm2(position)**2)
   
end program test_lunar_gravity
