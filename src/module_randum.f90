    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Module module_randum
    Implicit None
    contains 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !returns random number between 0 - 1 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function ran()    
    implicit none 
    integer , save :: flag = 0
    double precision :: ran 
    if(flag==0) then 
      call random_seed()
      flag = 1 
    endif 
    call random_number(ran)     ! built in fortran 90 random number function  
    end function ran
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !returns randum with mean and standard deviation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function normal(mean,sigma) 
    implicit none 
    integer :: flag 
    double precision, parameter :: pi = 3.141592653589793239  
    double precision :: u1, u2, y1, y2, normal, mean, sigma 
    save flag 
    data flag /0/ 
    u1 = ran(); u2 = ran() 
    if (flag.eq.0) then 
      y1 = sqrt(-2.0d0*log(u1))*cos(2.0d0*pi*u2) 
      normal = mean + sigma*y1 
      flag = 1 
    else 
      y2 = sqrt(-2.0d0*log(u1))*sin(2.0d0*pi*u2) 
      normal = mean + sigma*y2 
      flag = 0 
    endif  
    end function normal 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !returns a randum between vmin and vmax
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function ran2(vmin, vmax)
    implicit none
    double precision ran2
    double precision vmin, vmax
    ran2 = vmin + (vmax - vmin)*ran()
    end function ran2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 在n维向量var中随机选择一个元素，每个元素被选到的概率相同
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function ran3(n, var)
    implicit none
    integer, intent(in) :: n
    double precision var(n)
    double precision ran3, ran1
    integer m
    ran1 = ran() ! ran1 is a randum between 0 and 1
    m = floor(ran1/(1d0/dble(n))) + 1
    ran3 = var(m)
    endfunction ran3
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 在n维向量var中随机挑选一个元素，每个元素被选到的概率与p中对应的数值成正比；
    ! 例如：n = 3，var = [2, 5, 7], p = [3, 3, 5], 那么返回2的概率是3/11，返回
    ! 5的概率是3/11，返回7的概率是5/11；
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double precision function ran4(n, var, p)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: var(n)
    double precision, intent(in) :: p(n)
    double precision q(n), s
    integer i
    !----------------------------------------------------------------------
    
    ! n >= 2, 否则报错
    
    ! p 中均为非负值，否则报错
    
    ! 生成q向量
    do i = 1, n
        if (i == 1) then
            q(1) = p(1)/sum(p)
            cycle
        endif
        q(i) = q(i - 1) + p(i)/sum(p)
    enddo
    
    ! 生成随机结果
    s = ran()
    do i = n, 1, -1
        if (i == 1) then
            ran4 = var(1)
            exit
        endif
        if (s > q(i-1)) then
            ran4 = var(i)
            exit
        endif
    enddo
    
    endfunction ran4
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 生成一个n维向量v，v中的每一个元素都是0-1之间的随机数
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function ran5(n) result(v)
    implicit none
    integer, intent(in) :: n
    double precision :: v(n)
    integer i
    do i = 1, n
        v(i) = ran()
    enddo
    endfunction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    End Module module_randum
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine test_ran
    use module_randum
    Implicit None
    Integer , parameter :: N = 10000
    Real( Kind = 8 ) :: a( N )
    integer :: i
    Open( 12 , File = 'out.txt' )
    Do i = 1 , N
        a( i ) = normal( 5.0D0 , 2.0D0 )
        write( 12 , * ) a(i)
    End Do
    Close( 12 )
    End subroutine test_ran