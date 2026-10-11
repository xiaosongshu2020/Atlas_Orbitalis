!==============================================================================
!  mod_cfgfile.f90  ---  简易配置文件读取模块
!------------------------------------------------------------------------------
!  文件格式:
!     1) 空行                                   -> 忽略
!     2) 以 '#' 开头的行(前面允许有空格)        -> 注释, 忽略
!     3) 标量项:   key = value
!     4) 块(数组)项:
!            begin_xxx
!            <一行或多行内容>
!            end_xxx
!
!  说明:
!     * 键名查找不区分大小写; 同名键出现多次时取最后一条(后面的覆盖前面的)
!     * 自动去掉每行首尾的空格/Tab, 自动兼容 Windows(CRLF) / Unix(LF) 换行
!     * 数组下标从 1 开始(Fortran 惯例)
!     * stat/errmsg 记录第一次出错的信息, 用 ok() 判断
!
!  【重要】Fortran 的字符数组所有元素长度相同, 因此 get(begin,end) 返回的数组
!     元素为定长字符串(长度 = 块内最长行), 较短的行尾部补空格。取真实内容请用
!     trim(); 参与拼接(如拼路径)时尤其要注意, 例如 open(file=trim(eph(i))).
!
!  用法示例:
!     use mod_cfgfile
!     type(type_file) :: cfg
!     character(len=:), allocatable :: ts, te
!     character(len=:), allocatable :: eph(:)
!
!     call cfg%load('OD_input.cfg')
!     if (.not. cfg%ok()) print '(a)', cfg%errmsg
!
!     ts  = cfg%get('ts')                                ! '2030-01-03T00:00:00.000'
!     te  = cfg%get('te')                                ! '2030-01-04T00:00:00.000'
!     n   = cfg%get_int('nNavSat')                       ! 15
!     eph = cfg%get('begin_navSatEph', 'end_navSatEph')  ! 字符串数组(15 行)
!==============================================================================
module module_input
   implicit none
   private

   ! 单行最大长度(超出报错, 避免静默截断)
   integer, parameter, public :: MAXLINE = 2048

   !-----------------------------------------------------------------------
   !  单个条目 : 标量项 或 块项
   !-----------------------------------------------------------------------
   type :: type_item
      character(len=:), allocatable :: key        ! 标量: '='左边的键名; 块: begin_xxx 标志行
      character(len=:), allocatable :: value      ! 标量值
      character(len=:), allocatable :: vals(:)    ! 块内容, 每行一个元素
      character(len=:), allocatable :: key_end    ! 块的 end_xxx 标志行(文件中实际出现的)
      logical :: is_block = .false.
   end type type_item

   !-----------------------------------------------------------------------
   !  文件对象
   !-----------------------------------------------------------------------
   type, public :: type_file
      character(len=:), allocatable :: filename   ! 文件名
      integer :: stat = 0                         ! 0=正常, 非0=出错
      character(len=:), allocatable :: errmsg     ! 错误信息
      logical :: loaded = .false.                 ! 是否已成功加载
      integer :: n = 0                            ! 有效条目个数
      type(type_item), allocatable :: entry(:)    ! 条目表
   contains
      procedure :: load                           ! 加载并解析文件
      procedure :: has                            ! 键(标量或块)是否存在
      procedure :: get_int                        ! 取整数
      procedure :: get_real                       ! 取实数
      procedure :: ok                             ! 是否无错
      procedure :: show                           ! 打印解析结果(调试用)
      procedure, private :: get_one               ! get('key')         -> 字符串
      procedure, private :: get_many              ! get('begin','end') -> 字符串数组
      generic   :: get => get_one, get_many
   end type type_file

contains

   !========================================================================
   !  加载文件
   !========================================================================
   subroutine load(this, filename)
      class(type_file), intent(inout) :: this
      character(len=*), intent(in)    :: filename

      integer :: u, ios, nlin, ieq
      character(len=MAXLINE) :: raw
      character(len=:), allocatable :: line, k, v
      character(len=256) :: buf
      type(type_item) :: cur
      logical :: exist, in_block

      call reset(this)
      this%filename = trim(filename)

      inquire(file=trim(filename), exist=exist)
      if (.not. exist) then
         call set_err(this, 1, 'cfg file not found: ' // trim(filename))
         return
      end if

      open(newunit=u, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call set_err(this, 2, 'cannot open cfg file: ' // trim(filename))
         return
      end if

      in_block = .false.
      nlin = 0
      do
         read(u, '(A)', iostat=ios) raw
         if (ios /= 0) exit                    ! ios<0: 读到文件尾
         nlin = nlin + 1

         line = trim_line(raw)                 ! 去 CR/LF/Tab 及首尾空白
         if (len(line) == 0) cycle             ! 空行
         if (line(1:1) == '#') cycle           ! 注释行

         ! 读满整个缓冲区(且末尾不是 CR) => 该行可能已被静默截断
         if (len_trim(raw) >= MAXLINE .and. raw(MAXLINE:MAXLINE) /= achar(13)) then
            write(buf, '(a,i0)') 'line too long (>= 2048 chars) at line ', nlin
            call set_err(this, 3, trim(buf))
            exit
         end if

         !---- 块内部 : 收内容, 直到遇到 end_xxx ----
         if (in_block) then
            if (is_end_marker(line)) then
               cur%key_end = line
               if (.not. allocated(cur%vals)) allocate(character(len=1) :: cur%vals(0))
               call push(this, cur)
               in_block = .false.
            else
               call append_line(cur%vals, line)
            end if
            cycle
         end if

         !---- 标量项 : key = value ----
         ieq = index(line, '=')
         if (ieq > 0) then
            k = trim(line(:ieq-1))
            v = trim(adjustl(line(ieq+1:)))
            if (len(k) == 0) then
               write(buf, '(a,i0)') 'missing key name before "=" at line ', nlin
               call set_err(this, 4, trim(buf))
               exit
            end if
            call clear_entry(cur)
            cur%key   = k
            cur%value = v
            call push(this, cur)

         !---- 块开始标志行 ----
         else if (index(lower(line), 'begin_') == 1) then
            call clear_entry(cur)
            cur%key      = line
            cur%is_block = .true.
            in_block     = .true.

         else if (is_end_marker(line)) then
            write(buf, '(a,i0)') 'unexpected "end_" marker (no matching "begin_") at line ', nlin
            call set_err(this, 4, trim(buf))
            exit

         else
            write(buf, '(a,i0,a)') 'bad line (not "key = value" nor block marker) at line ', &
                                   nlin, ': ' // line
            call set_err(this, 4, trim(buf))
            exit
         end if
      end do
      close(u)

      if (in_block) then
         call set_err(this, 5, 'missing "end_" marker for block: ' // cur%key)
         return
      end if
      if (this%stat == 0) this%loaded = .true.
   end subroutine load

   !========================================================================
   !  get('key') -> 字符串
   !========================================================================
   function get_one(this, key) result(val)
      class(type_file), intent(inout) :: this
      character(len=*), intent(in)    :: key
      character(len=:), allocatable   :: val
      integer :: idx

      val = ''
      if (.not. this%loaded) return
      idx = find_key(this, key)
      if (idx == 0) then
         call set_err(this, 11, 'key not found: ' // trim(key))
      else if (this%entry(idx)%is_block) then
         call set_err(this, 12, 'key is a block, use get(begin_key, end_key): ' // trim(key))
      else
         val = this%entry(idx)%value
      end if
   end function get_one

   !========================================================================
   !  get('begin_xxx','end_xxx') -> 字符串数组(元素定长, 用 trim() 取内容)
   !========================================================================
   function get_many(this, kbegin, kend) result(vals)
      class(type_file), intent(inout) :: this
      character(len=*), intent(in)    :: kbegin
      character(len=*), intent(in)    :: kend
      character(len=:), allocatable   :: vals(:)
      integer :: idx

      allocate(character(len=1) :: vals(0))     ! 默认返回空数组
      if (.not. this%loaded) return
      idx = find_key(this, kbegin)
      if (idx == 0) then
         call set_err(this, 11, 'block not found: ' // trim(kbegin))
         return
      end if
      if (.not. this%entry(idx)%is_block) then
         call set_err(this, 12, 'key is not a block: ' // trim(kbegin))
         return
      end if
      if (lower(trim(kend)) /= lower(trim(this%entry(idx)%key_end))) then
         call set_err(this, 13, 'end marker mismatch: expect "' // trim(this%entry(idx)%key_end) &
                           // '", but got "' // trim(kend) // '"')
         return
      end if
      vals = this%entry(idx)%vals               ! 整体赋值, 自动重分配
   end function get_many

   !========================================================================
   !  其它访问方法
   !========================================================================
   ! 键(标量或块)是否存在
   logical function has(this, key)
      class(type_file), intent(in) :: this
      character(len=*), intent(in) :: key
      has = (find_key(this, key) > 0)
   end function has

   ! 取整数; 键不存在或无法转换时返回 default(未给 default 则报错)
   integer function get_int(this, key, default) result(v)
      class(type_file), intent(inout) :: this
      character(len=*), intent(in)    :: key
      integer, intent(in), optional   :: default
      integer :: idx, ios
      character(len=:), allocatable :: s

      v = 0
      if (present(default)) v = default
      if (.not. this%loaded) return
      idx = find_key(this, key)
      if (idx == 0) then
         if (.not. present(default)) call set_err(this, 11, 'key not found: ' // trim(key))
         return
      else if (this%entry(idx)%is_block) then
         call set_err(this, 12, 'key is a block: ' // trim(key))
         return
      end if
      s = this%entry(idx)%value
      read(s, *, iostat=ios) v
      if (ios /= 0) then
         if (present(default)) then
            v = default
         else
            call set_err(this, 14, 'cannot convert to integer: ' // trim(key) // ' = ' // s)
         end if
      end if
   end function get_int

   ! 取实数; 键不存在或无法转换时返回 default(未给 default 则报错)
   real(kind(1.0d0)) function get_real(this, key, default) result(v)
      class(type_file), intent(inout) :: this
      character(len=*), intent(in)    :: key
      real(kind(1.0d0)), intent(in), optional :: default
      integer :: idx, ios
      character(len=:), allocatable :: s

      v = 0.0d0
      if (present(default)) v = default
      if (.not. this%loaded) return
      idx = find_key(this, key)
      if (idx == 0) then
         if (.not. present(default)) call set_err(this, 11, 'key not found: ' // trim(key))
         return
      else if (this%entry(idx)%is_block) then
         call set_err(this, 12, 'key is a block: ' // trim(key))
         return
      end if
      s = this%entry(idx)%value
      read(s, *, iostat=ios) v
      if (ios /= 0) then
         if (present(default)) then
            v = default
         else
            call set_err(this, 14, 'cannot convert to real: ' // trim(key) // ' = ' // s)
         end if
      end if
   end function get_real

   logical function ok(this)
      class(type_file), intent(in) :: this
      ok = (this%stat == 0)
   end function ok

   ! 打印解析结果(调试用), unit 缺省为标准输出
   subroutine show(this, unit)
      class(type_file), intent(in) :: this
      integer, intent(in), optional :: unit
      integer :: u, i, j, ne
      character(len=:), allocatable :: fn

      u = 6
      if (present(unit)) u = unit
      fn = '<none>'
      if (allocated(this%filename)) fn = this%filename
      write(u, '(a)') '---- cfg file: ' // fn // ' ----'
      do i = 1, this%n
         if (this%entry(i)%is_block) then
            ne = 0
            if (allocated(this%entry(i)%vals)) ne = size(this%entry(i)%vals)
            write(u, '(a,i0,a)') trim(this%entry(i)%key) // '   [block, ', ne, ' lines]'
            do j = 1, ne
               write(u, '(4x,a)') trim(this%entry(i)%vals(j))
            end do
            write(u, '(a)') trim(this%entry(i)%key_end)
         else
            write(u, '(a)') trim(this%entry(i)%key) // ' = ' // this%entry(i)%value
         end if
      end do
      write(u, '(a)') '--------------------------------------'
   end subroutine show

   !========================================================================
   !  内部工具(模块过程, 不对外)
   !========================================================================
   ! 记录第一次出错的信息
   subroutine set_err(this, code, msg)
      class(type_file), intent(inout) :: this
      integer, intent(in)          :: code
      character(len=*), intent(in) :: msg
      if (this%stat == 0) then
         this%stat   = code
         this%errmsg = trim(msg)
      end if
   end subroutine set_err

   ! 清空对象状态
   subroutine reset(this)
      class(type_file), intent(inout) :: this
      if (allocated(this%entry))    deallocate(this%entry)
      if (allocated(this%errmsg))   deallocate(this%errmsg)
      if (allocated(this%filename)) deallocate(this%filename)
      this%n      = 0
      this%stat   = 0
      this%loaded = .false.
   end subroutine reset

   subroutine clear_entry(e)
      type(type_item), intent(inout) :: e
      if (allocated(e%key))     deallocate(e%key)
      if (allocated(e%value))   deallocate(e%value)
      if (allocated(e%vals))    deallocate(e%vals)
      if (allocated(e%key_end)) deallocate(e%key_end)
      e%is_block = .false.
   end subroutine clear_entry

   ! 追加一个条目(容量不够时自动翻倍)
   subroutine push(this, e)
      class(type_file), intent(inout) :: this
      type(type_item), intent(in)     :: e
      type(type_item), allocatable    :: tmp(:)
      if (.not. allocated(this%entry)) then
         allocate(this%entry(16))
      else if (this%n >= size(this%entry)) then
         allocate(tmp(2 * size(this%entry)))
         tmp(1:this%n) = this%entry(1:this%n)
         call move_alloc(tmp, this%entry)
      end if
      this%n = this%n + 1
      this%entry(this%n) = e
   end subroutine push

   ! 给块追加一行。注意: 字符数组所有元素等长, 一旦分配就无法改变元素长度,
   ! 因此每次都按"已有最长行 + 新行"重新分配整个数组, 保证不截断。
   subroutine append_line(vals, line)
      character(len=:), allocatable, intent(inout) :: vals(:)
      character(len=*), intent(in) :: line
      character(len=:), allocatable :: tmp(:)
      integer :: n, nlen
      n = 0
      nlen = len_trim(line)
      if (allocated(vals)) then
         n = size(vals)
         if (n > 0) nlen = max(nlen, len(vals(1)))
      end if
      if (nlen < 1) nlen = 1
      allocate(character(len=nlen) :: tmp(n + 1))
      if (n > 0) tmp(1:n) = vals(1:n)
      tmp(n + 1) = trim(line)
      call move_alloc(tmp, vals)
   end subroutine append_line

   ! 查键, 返回最后一条匹配项的下标(0 = 未找到), 不区分大小写
   pure function find_key(this, key) result(idx)
      class(type_file), intent(in) :: this
      character(len=*), intent(in) :: key
      integer :: idx, i
      idx = 0
      do i = this%n, 1, -1
         if (allocated(this%entry(i)%key)) then
            if (lower(trim(this%entry(i)%key)) == lower(trim(key))) then
               idx = i
               return
            end if
         end if
      end do
   end function find_key

   ! 是否为块结束标志行(以 end_ 开头, 不区分大小写)
   pure logical function is_end_marker(line)
      character(len=*), intent(in) :: line
      is_end_marker = (index(lower(trim(line)), 'end_') == 1)
   end function is_end_marker

   ! 去 CR/LF/Tab/NUL 及首尾空白
   pure function trim_line(raw) result(line)
      character(len=*), intent(in) :: raw
      character(len=:), allocatable :: line
      character(len=len(raw)) :: buf
      integer :: i, n
      buf = raw
      do i = 1, len(buf)
         select case (buf(i:i))
         case (achar(9), achar(10), achar(13), achar(0))   ! Tab, LF, CR, NUL
            buf(i:i) = ' '
         end select
      end do
      buf = adjustl(buf)
      n = len_trim(buf)
      line = buf(:n)
   end function trim_line

   ! 转小写
   pure function lower(s) result(t)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: t
      integer :: i, c
      t = s
      do i = 1, len(s)
         c = iachar(t(i:i))
         if (c >= iachar('A') .and. c <= iachar('Z')) t(i:i) = achar(c + 32)
      end do
   end function lower

end module module_input
