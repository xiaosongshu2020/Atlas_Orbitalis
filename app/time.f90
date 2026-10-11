    ! Demonstration of module_time
    ! 2025-12-19 by Shawn Hu
    program demo_time
    use module_time
    implicit none
    
    type(type_t) :: time1, time2, time3, time4, time5
    real*8 :: week, sow, diff
    
    ! 1. All initialization methods create UTC type objects
    call time1%init(51544.5d0)                     ! MJD -> UTC
    call time2%init(51544.0d0, 43200.0d0)         ! day-second -> UTC  
    call time3%init(2000.0d0,1.0d0,1.0d0,12.0d0,0.0d0,0.0d0)  ! calendar -> UTC
    print*, "time1:", time1%day, time1%sec, time1%val, time1%ttype
    
    ! 2. Convert UTC to other time systems
    time4 = time1
    call time4%toTT()                           ! UTC -> TT
    print*, "time4 (TT):", time4%day, time4%sec, time4%val, time4%ttype
    time5 = time1
    call time5%toBDT(week, sow)                  ! TT -> BDT
    print*, "time5 (BDT):", week, sow
    
    ! 3. Time arithmetic works with any type
    diff = timeMinus(time2, time1)              ! time2 - time1 (days)
    time2 = timeAdd(time1, 1.5d0)              ! time1 + 1.5 days
    
    end program demo_time
