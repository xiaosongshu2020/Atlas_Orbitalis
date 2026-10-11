program detect_os_env
    implicit none
    character(len=64) :: os_value
    integer :: status

    ! 尝试读取环境变量 OS (Windows 系统常用)
    call get_environment_variable("OS", os_value, status=status)

    if (status == 0) then
        ! 成功读取到 OS 变量
        if (trim(os_value) == "Windows_NT") then
            print *, "Running on Windows"
            ! 这里放 Windows 专用代码
        else
            print *, "Unknown OS (OS variable not Windows_NT)"
            print *, "OS: ", trim(os_value)
        end if
    else
        ! 读取失败，大概率是 Linux / Unix 系统
        print *, "Running on Linux / Unix (OS env var not found)"
        ! 这里放 Linux 专用代码
    end if

end program detect_os_env