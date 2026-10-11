   
! module operator_i

!     interface operator(.i.) !!// �Զ������ؾ������������
!     module procedure brinv !!// ʵ��������
!     module procedure bcinv !!// ����������
!     end interface

! contains

! function brinv(re) result(r)
!     use module_terminalOutput
!     double precision,intent(in) :: re(:,:) !!// ԭ����
!     double precision :: r(size(re,1),size(re,1)) !!// �����
!     integer*4 :: flag,n !!// flag�ж������ԣ�n�Ǿ���ά��
!     double precision :: t,d !!// �м����
!     integer*4 :: is(size(re,1)),js(size(re,1)) !!// �м����

!     n=size(re,1) !!// ��ȡ����ά��
!     r=re !!// ����ֵ
!     flag=1
!     do k=1,n
!         d=0.0
!         do i=k,n
!             do j=k,n
!                 if (abs(r(i,j)).gt.d) then
!                     d=abs(r(i,j))
!                     is(k)=i
!                     js(k)=j
!                 end if
!             end do
!         end do
!         if (d+1.0.eq.1.0) then
!             flag=0
!             write(*,*) "Error: failed to inverse a matrix as it's singular."
!             errCode = 1001
!             return !!// �������죬�˳���������
!         end if
!         do j=1,n
!             t=r(k,j)
!             r(k,j)=r(is(k),j)
!             r(is(k),j)=t
!         end do
!         do i=1,n
!             t=r(i,k)
!             r(i,k)=r(i,js(k))
!             r(i,js(k))=t
!         end do
!         r(k,k)=1/r(k,k)
!         do j=1,n
!             if (j.ne.k) r(k,j)=r(k,j)*r(k,k)
!         end do
!         do i=1,n
!             if (i.ne.k) then
!                 do j=1,n
!                     if (j.ne.k) r(i,j)=r(i,j)-r(i,k)*r(k,j)
!                 end do
!             end if
!         end do
!         do i=1,n
!             if (i.ne.k) r(i,k)=-r(i,k)*r(k,k)
!         end do
!     end do
!     do k=n,1,-1
!         do j=1,n
!             t=r(k,j)
!             r(k,j)=r(js(k),j)
!             r(js(k),j)=t
!         end do
!         do i=1,n
!             t=r(i,k)
!             r(i,k)=r(i,is(k))
!             r(i,is(k))=t
!         end do
!     end do
! end function
! !!// ������������Ĵ��� trans from ��ʿ����Fortran�����㷨����
! !!// ���ߣ�zuozhihua ʱ�䣺2020/8/10 �ص㣺����
! function bcinv(cpx)
!     complex*16,intent(in) :: cpx(:,:) !!// ԭ����
!     complex*16 :: bcinv(size(cpx,1),size(cpx,2)) !!// �����
!     integer*4 :: flag,n !!// flag�ж������ԣ�n�Ǿ���ά��
!     double precision :: ar(size(cpx,1),size(cpx,1)),ai(size(cpx,1),size(cpx,1)) !!// ʵ������ar���鲿����ai
!     double precision :: d,p,t,q,s,b !!// �м����
!     integer*4 :: is(size(cpx,1)),js(size(cpx,1)) !!// �м����

!     n=size(cpx,1)
!     forall(i=1:n,j=1:n)
!         ar(i,j) = real(cpx(i,j));ai(i,j) = imag(cpx(i,j))
!     end forall
!     flag=1
!     do k=1,n
!         d=0.0
!         do i=k,n
!             do j=k,n
!                 p=ar(i,j)*ar(i,j)+ai(i,j)*ai(i,j)
!                 if(p.gt.d) then
!                     d=p
!                     is(k)=i
!                     js(k)=j
!                 end if
!             end do
!         end do
!         if(d+1.0.eq.1.0) then
!             flag=0
!             write(*,*) "flag=0,���������죡"
!             return !!// �������죬�˳���������
!         end if
!         do j=1,n
!             t=ar(k,j)
!             ar(k,j)=ar(is(k),j)
!             ar(is(k),j)=t
!             t=ai(k,j)
!             ai(k,j)=ai(is(k),j)
!             ai(is(k),j)=t
!         end do
!         do i=1,n
!             t=ar(i,k)
!             ar(i,k)=ar(i,js(k))
!             ar(i,js(k))=t
!             t=ai(i,k)
!             ai(i,k)=ai(i,js(k))
!             ai(i,js(k))=t
!         end do
!         ar(k,k)=ar(k,k)/d
!         ai(k,k)=-ai(k,k)/d
!         do j=1,n
!             if(j.ne.k) then
!                 p=ar(k,j)*ar(k,k)
!                 q=ai(k,j)*ai(k,k)
!                 s=(ar(k,j)+ai(k,j))*(ar(k,k)+ai(k,k))
!                 ar(k,j)=p-q
!                 ai(k,j)=s-p-q
!             end if
!         end do
!         do i=1,n
!             if(i.ne.k) then
!                 do j=1,n
!                     if (j.ne.k) then
!                         p=ar(k,j)*ar(i,k)
!                         q=ai(k,j)*ai(i,k)
!                         s=(ar(k,j)+ai(k,j))*(ar(i,k)+ai(i,k))
!                         t=p-q
!                         b=s-p-q
!                         ar(i,j)=ar(i,j)-t
!                         ai(i,j)=ai(i,j)-b
!                     end if
!                 end do
!             end if
!         end do
!         do i=1,n
!             if (i.ne.k) then
!                 p=ar(i,k)*ar(k,k)
!                 q=ai(i,k)*ai(k,k)
!                 s=(ar(i,k)+ai(i,k))*(ar(k,k)+ai(k,k))
!                 ar(i,k)=q-p
!                 ai(i,k)=p+q-s
!             end if
!         end do
!     end do
!     do k=n,1,-1
!         do j=1,n
!             t=ar(k,j)
!             ar(k,j)=ar(js(k),j)
!             ar(js(k),j)=t
!             t=ai(k,j)
!             ai(k,j)=ai(js(k),j)
!             ai(js(k),j)=t
!         end do
!         do i=1,n
!             t=ar(i,k)
!             ar(i,k)=ar(i,is(k))
!             ar(i,is(k))=t
!             t=ai(i,k)
!             ai(i,k)=ai(i,is(k))
!             ai(i,is(k))=t
!         end do
!     end do
!     forall(i=1:n,j=1:n)
!         bcinv(i,j) = cmplx(ar(i,j),ai(i,j),8)
!     end forall

! end function
! end module


! ! ���������㷨
! !Input
! !   n   integer    ��������ά��
! !   m   double     n*nά�������
! !Output
! !   mr  double     n*n��m�������

! subroutine matrix_inverse_2(n, m, mr)
!     use operator_i
!     implicit none
!     integer n
!     double precision m(n, n), mr(n, n)
!     !-------------------------------------------------------------

!     mr = .i.m

! end subroutine

    module module_matrix
    implicit none
    contains
    !   This code is obtained from https://github.com/cfinch/Shocksolution_Examples/tree/master/FORTRAN/LinearEquationSolver,
    !       which, as is introduced, is implemented based on a NR algorithm.
    !   The code is orignally designed as a module to be a linear equation solver,
    !       while in GIANT (former Delta), the code is re-wrapped as a function to inverse a matrix.
    !       Nevertheless the original functionality is still fully accessible.
    !   The only revision to the original code is to change "real" and "double precision" to "double precision".

    ! >>>>>>>>>>>>>>>>>>>>>> GIANT wrapper beyond this line <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    subroutine InvMatrix_LUsolver(N,A,Ainv)
        !
        !   Introduction:
        !       This subroutine inverses square matrix A with a dimension of N, and Ainv holds the result.
        !   Date:
        !       2017/09/09
        !   Note:
        !       1. The core solver is an open-source code, obtained from github. The credit goes to the original author.
        !       2. The inverse code is implemneted by solving the linear equation per column of the identity matrix.
        !
        use MathsMatrix
        implicit none
        ! input
        integer,intent(in):: N
        double precision,intent(in):: A(N,N)
        ! output
        double precision,intent(out):: Ainv(N,N)
        ! local
        integer:: i,D,rc,INDX(N)
        double precision:: A0(N,N),B(N)
        integer clock_rate, t1, t2, t3

        A0 = A
        ! call LU decomposition routine
        call LUDCMP(A0,N,INDX,D,rc)
        if (rc.eq.1) then
            !call RaiseError(1016,subName)
            write(*, *) "Error:: failed to inverse a matrix. (subroutine InvMatrix_LUsolver)"
            write(*, *) "Input matrix is:"
            do i = 1, N
                write(*, *) A(i, :)
            enddo
            STOP
        endif

        ! call appropriate solver if previous return code is ok
        if (rc.eq.0) then
            !do i = 1,N
            !    B = 0d0
            !    B(i) = 1d0
            !    call LUBKSB(A0,N,INDX,B)
            !    Ainv(1:N,i) = B(1:N)
            !enddo
            do i = 1, N
                B = 0d0
                B(i) = 1d0
                call LUBKSB(A0, N, INDX, B)
                Ainv(1:N, i) = B(1:N)
            enddo
        endif      

        return
    end subroutine InvMatrix_LUsolver


    ! >>>>>>>>>>>>>>>>>>>>>> Original code beyond this line <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !*******************************************************
    !*    LU decomposition routines used by test_lu.f90    *
    !*                                                     *
    !*                 F90 version by J-P Moreau, Paris    *
    !* --------------------------------------------------- *
    !* Reference:                                          *
    !*                                                     *
    !* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
    !*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
    !*  University Press, 1986" [BIBLI 08].                *
    !*                                                     *
    !*******************************************************
    !
    ! NOTE: third-party code.  LUDCMP/LUBKSB below are J-P Moreau's F90
    ! translation of the LU algorithms from "Numerical Recipes"; they are
    ! kept here with their original attribution and are not relicensed by
    ! this project -- their own source terms apply.

    !  ***************************************************************
    !  * Given an N x N matrix A, this routine replaces it by the LU *
    !  * decomposition of a rowwise permutation of itself. A and N   *
    !  * are input. INDX is an output vector which records the row   *
    !  * permutation effected by the partial pivoting; D is output   *
    !  * as -1 or 1, depending on whether the number of row inter-   *
    !  * changes was even or odd, respectively. This routine is used *
    !  * in combination with LUBKSB to solve linear equations or to  *
    !  * invert a matrix. Return code is 1, if matrix is singular.   *
    !  ***************************************************************
     Subroutine LUDCMP(A,N,INDX,D,CODE)
     IMPLICIT NONE
     integer, intent(in) :: N
     double precision, intent(inout) :: A(N, N)
     integer, intent(out) :: D, CODE
     integer, intent(out) :: INDX(N)
     
     integer, parameter :: nmax = 100
     double precision, parameter :: tiny = 1.5D-16
     double precision  :: AMAX, DUM, SUMM, VV(N)
     INTEGER :: i, j, k, imax

     D=1; CODE=0
    
     ! �ɲ���ѭ��
     DO I=1,N
       AMAX=0.d0
       DO J=1,N
         IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
       END DO ! j loop
       IF(AMAX.LT.TINY) THEN
         CODE = 1
         RETURN
       END IF
       VV(I) = 1.d0 / AMAX
     END DO ! i loop

     ! ���ɲ���ѭ��
     DO J=1,N
       DO I=1,J-1
         SUMM = A(I,J)
         DO K=1,I-1
           SUMM = SUMM - A(I,K)*A(K,J)
         END DO ! k loop
         A(I,J) = SUMM
       END DO ! i loop
       AMAX = 0.d0
       DO I=J,N
         SUMM = A(I,J)
         DO K=1,J-1
           SUMM = SUMM - A(I,K)*A(K,J)
         END DO ! k loop
         A(I,J) = SUMM
         DUM = VV(I)*DABS(SUMM)
         IF(DUM.GE.AMAX) THEN
           IMAX = I
           AMAX = DUM
         END IF
       END DO ! i loop

       IF(J.NE.IMAX) THEN
         DO K=1,N
           DUM = A(IMAX,K)
           A(IMAX,K) = A(J,K)
           A(J,K) = DUM
         END DO ! k loop
         D = -D
         VV(IMAX) = VV(J)
       END IF

       INDX(J) = IMAX
       IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

       IF(J.NE.N) THEN
         DUM = 1.d0 / A(J,J)
         DO I=J+1,N
           A(I,J) = A(I,J)*DUM
         END DO ! i loop
       END IF
     END DO ! j loop

     RETURN
     END subroutine LUDCMP


    !  ******************************************************************
    !  * Solves the set of N linear equations A . X = B.  Here A is     *
    !  * input, not as the matrix A but rather as its LU decomposition, *
    !  * determined by the routine LUDCMP. INDX is input as the permuta-*
    !  * tion vector returned by LUDCMP. B is input as the right-hand   *
    !  * side vector B, and returns with the solution vector X. A, N and*
    !  * INDX are not modified by this routine and can be used for suc- *
    !  * cessive calls with different right-hand sides. This routine is *
    !  * also efficient for plain matrix inversion.                     *
    !  ******************************************************************
    Subroutine LUBKSB(A, N, INDX, B)
    implicit none
    integer, intent(in) :: N
    double precision, intent(in), dimension(N,N) :: A
    integer, intent(in), dimension(N) :: INDX
    double precision, intent(inout), dimension(N) :: B
    !f2py depend(N) A, INDX, B
    integer :: i, j, ii, ll

    double precision  SUMM

    II = 0

    DO I=1,N
    LL = INDX(I)
    SUMM = B(LL)
    B(LL) = B(I)
    IF(II.NE.0) THEN
        DO J=II,I-1
        SUMM = SUMM - A(I,J)*B(J)
        END DO ! j loop
    ELSE IF(SUMM.NE.0.d0) THEN
        II = I
    END IF
    B(I) = SUMM
    END DO ! i loop

    DO I=N,1,-1
    SUMM = B(I)
    IF(I < N) THEN
        DO J=I+1,N
        SUMM = SUMM - A(I,J)*B(J)
        END DO ! j loop
    END IF
    B(I) = SUMM / A(I,I)
    END DO ! i loop

    RETURN
    END subroutine LUBKSB

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! save a matrix into text file
    subroutine saveMatrix(n, m, mat, fn, fmt)
    implicit none
    integer, intent(in) :: n, m
    double precision, intent(in) :: mat(n, m)
    character(len=*), intent(in) :: fn
    character(len=*), intent(in) :: fmt
    integer i
    !-----------------------------------------------------------
    open(801, file = trim(fn))
    do i = 1, n
        write(801, fmt) mat(i, :)
    enddo
    close(801)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine readMatrix(fn, n, m, mat, headLines)
    implicit none
    integer, intent(in) :: n, m
    character(len=*), intent(in) :: fn
    double precision, intent(out) :: mat(n, m)
    integer, intent(in), optional :: headLines
    integer i
    character*1 tmp
    open(802, file = trim(fn))
    if (present(headLines)) then
        if (headLines > 0) then
            do i = 1, headLines
                read(802, *) tmp
            enddo
        endif
    endif
    do i = 1, n
        read(802, *) mat(i, :)
    enddo
    close(802)
    endsubroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This subroutine solves this equation: DY = B * dx
    !process is:
    !   DY = B * dx
    !=> BT * DY = BT * B * dx
    !=> dx = (BT * B)**(-1) * BT * DY
    subroutine least_squares_solution(n1, n2, B, DY, dx)
    implicit none
    integer, intent(in) :: n1, n2
    double precision, intent(in) :: B(n1, n2)
    double precision, intent(in) :: DY(n1)
    double precision, intent(out) :: dx(n2)
    double precision  BT(n2, n1), BTB(n2, n2), iBTB(n2, n2)
    integer nVar, nData
    !------------------------------------------------------
    nVar = size(dx)
    nData = size(DY)
    
    BT = transpose(B)
    BTB = matmul(BT, B)
    ! call matrix_inverse_2(nVar, BTB, iBTB)
    call InvMatrix_LUsolver(nVar, BTB, iBTB)
    dx = matmul(iBTB, matmul(BT, DY))

    endsubroutine    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end module
 

