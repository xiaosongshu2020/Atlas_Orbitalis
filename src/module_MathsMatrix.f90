module MathsMatrix
    !----------------------------------------------------------------------------------------------------------------------------------
    !   Instruction:
    !       This file collects very very fundamental matrix manipulations.
    !   Date:
    !       2014-06-15
    !   Update:
    !       2017/02/06   Rotation routines that were once in the rotation module are moved here.
    !       2017/03/16   This function collection is now compiled as a module.
    !----------------------------------------------------------------------------------------------------------------------------------

    contains

        function eye(n) result(idm)
            !----------------------------------------------------------------------------------------------------------------------------------
            !   Instruction:
            !       This function version of IM, with the same function name as in MATLAB.
            !   Date:
            !       2017-11-1
            !   Update:
            !
            !   Inputs:
            !       n               i               dimension of the matrix
            !   Outputs:
            !       idm             d(n,n)          initialized identity matrix
            !   Subroutine code:
            !       feye
            !----------------------------------------------------------------------------------------------------------------------------------
            use ErrorHandler
            implicit none
            ! input
            integer,intent(in):: n
            ! output
            double precision,allocatable:: idm(:,:)
            ! local
            character(len=4):: subName='feye'
            ! code
            call checkin(subName)
            if (allocated(idm)) then
                if ((size(idm,1).ne.n).or.(size(idm,2).ne.n)) then
                    call RaiseError(0016,subName)
                    return
                endif
            else
                allocate(idm(n,n))
            endif
            call IM(n,idm)
            call checkout(subName)
            return
        end function eye
        !
        subroutine IM(n,x)
            !----------------------------------------------------------------------------------------------------------------------------------
            !   Instruction:
            !       This subroutine initialized an n-by-n identity matrix
            !   Date:
            !       2014-06-15
            !   Update:
            !
            !   Inputs:
            !       n               i               dimension of the matrix
            !   Outputs:
            !       x               d(n,n)          initialized identity matrix
            !   Subroutine code:
            !       imat
            !----------------------------------------------------------------------------------------------------------------------------------
            use ErrorHandler
            implicit none
            ! input
            integer,intent(in):: n
            ! input/output
            double precision,intent(out):: x(n,n)
            ! local
            integer:: i
            character(len=4):: subName='imat'

            ! call checkin(subName)
            x = 0d0
            do i = 1,n
                x(i,i) = 1d0
            enddo
            ! call checkout(subName)
            return
        end subroutine IM

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        subroutine Rx(phi,M)
            !----------------------------------------------------------------------------------------------------------------------------------
            !   Instruction:
            !       This subroutine rotates the M matrix by phi (rad) along the x-axis.
            !   Date:
            !       2017/02/06
            !   Update:
            !   Inputs:
            !       phi             d               rotation angle
            !   Outputs:
            !       M               d(3,3)          rotation matrix
            !   Subroutine code:
            !       rotx
            !----------------------------------------------------------------------------------------------------------------------------------
            implicit none
            double precision phi,A(3,3),M(3,3)
            double precision cPhi,sPhi
            call im(3,A)
            cPhi = dcos(phi)
            sPhi = dsin(phi)
            A(2,2) = cPhi
            A(3,2) = -sPhi
            A(2,3) = sPhi
            A(3,3) = cPhi
            M = matmul(A,M)
            return
        end subroutine Rx

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        subroutine Ry(phi,M)
            !----------------------------------------------------------------------------------------------------------------------------------
            !   Instruction:
            !       This subroutine rotates the M matrix by phi (rad) along the y-axis.
            !   Date:
            !       2017/02/06
            !   Update:
            !   Inputs:
            !       phi             d               rotation angle
            !   Outputs:
            !       M               d(3,3)          rotation matrix
            !   Subroutine code:
            !       roty
            !----------------------------------------------------------------------------------------------------------------------------------
            implicit none
            double precision phi,A(3,3),M(3,3)
            double precision cPhi,sPhi
            call im(3,A)
            cPhi = dcos(phi)
            sPhi = dsin(phi)
            A(1,1) = cPhi
            A(3,1) = sPhi
            A(1,3) = -sPhi
            A(3,3) = cPhi
            M = matmul(A,M)
            return
        end subroutine Ry

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        subroutine Rz(phi,M)
            !----------------------------------------------------------------------------------------------------------------------------------
            !   Instruction:
            !       This subroutine rotates the M matrix by phi (rad) along the z-axis.
            !   Date:
            !       2017/02/06
            !   Update:
            !   Inputs:
            !       phi             d               rotation angle
            !   Outputs:
            !       M               d(3,3)          rotation matrix
            !   Subroutine code:
            !       rotz
            !----------------------------------------------------------------------------------------------------------------------------------
            implicit none
            double precision phi,A(3,3),M(3,3)
            double precision cPhi,sPhi
            call im(3,A)
            cPhi = dcos(phi)
            sPhi = dsin(phi)
            A(1,1) = cPhi
            A(2,1) = -sPhi
            A(1,2) = sPhi
            A(2,2) = cPhi
            M = matmul(A,M)
            return
        end subroutine Rz

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        subroutine dRx(phi,M)
            !----------------------------------------------------------------------------------------------------------------------------------
            !   Instruction:
            !       This subroutine computes the derivative of the rotation matrix along the x-axis.
            !   Date:
            !       2017/02/06
            !   Update:
            !   Inputs:
            !       phi             d               rotation angle
            !   Outputs:
            !       M               d(3,3)          derivative of rotation matrix
            !   Subroutine code:
            !       drox
            !----------------------------------------------------------------------------------------------------------------------------------
            implicit none
            double precision phi,A(3,3),M(3,3)
            double precision cPhi,sPhi
            cPhi = dcos(phi)
            sPhi = dsin(phi)
            A = 0d0
            A(2,2) = -sPhi
            A(3,2) = -cPhi
            A(2,3) = cPhi
            A(3,3) = -sPhi
            M = matmul(A,M)
            return
        end subroutine dRx

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        subroutine dRy(phi,M)
            !----------------------------------------------------------------------------------------------------------------------------------
            !   Instruction:
            !       This subroutine computes the derivative of the rotation matrix along the y-axis.
            !   Date:
            !       2017/02/06
            !   Update:
            !   Inputs:
            !       phi             d               rotation angle
            !   Outputs:
            !       M               d(3,3)          derivative of rotation matrix
            !   Subroutine code:
            !       droy
            !----------------------------------------------------------------------------------------------------------------------------------
            implicit none
            double precision phi,A(3,3),M(3,3)
            double precision cPhi,sPhi
            cPhi = dcos(phi)
            sPhi = dsin(phi)
            A = 0d0
            A(1,1) = -sPhi
            A(3,1) = cPhi
            A(1,3) = -cPhi
            A(3,3) = -sPhi
            M = matmul(A,M)
            return
        end subroutine dRy

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        subroutine dRz(phi,M)
            !----------------------------------------------------------------------------------------------------------------------------------
            !   Instruction:
            !       This subroutine computes the derivative of the rotation matrix along the z-axis.
            !   Date:
            !       2017/02/06
            !   Update:
            !   Inputs:
            !       phi             d               rotation angle
            !   Outputs:
            !       M               d(3,3)          derivative of rotation matrix
            !   Subroutine code:
            !       droz
            !----------------------------------------------------------------------------------------------------------------------------------
            implicit none
            double precision phi,A(3,3),M(3,3)
            double precision cPhi,sPhi
            cPhi = dcos(phi)
            sPhi = dsin(phi)
            A = 0d0
            A(1,1) = -sPhi
            A(2,1) = -cPhi
            A(1,2) = cPhi
            A(2,2) = -sPhi
            M = matmul(A,M)
            return
        end subroutine dRz

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        SUBROUTINE InvSPDMatrix(N,A0,A)
            use ErrorHandler
            implicit none
            ! input
            integer,intent(in):: N
            double precision,intent(in):: A0(N,N)
            ! output
            double precision,intent(out):: A(N,N)
            ! local
            integer:: K,M,I,J
            double precision:: B(N),W,G
            character(len=4):: subName='ispm'

            call checkin(subName)
            A = A0
            DO K=1,N
                M=N-K+1
                W=A(1,1)
                IF (W+1d0.EQ.1d0) THEN
                    call RaiseError(1003,subName)
                    RETURN
                END IF
                DO I=2,N
                    G=A(I,1)
                    B(I)=G/W
                    IF (I.LE.M) B(I)=-B(I)
                    DO J=2,I
                        A(I-1,J-1)=A(I,J)+G*B(J)
                    enddo
                enddo
                A(N,N)=1d0/W
                DO I=2,N
                    A(N,I-1)=B(I)
                enddo
            enddo
            DO I=1,N-1
                DO J=I+1,N
                    A(I,J)=A(J,I)
                enddo
            enddo

            call checkout(subName)
            RETURN
        END


! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        function checkSymmetry(mat) result(score)
            ! 2017/09/06
            use ErrorHandler
            implicit none
            ! input
            double precision,intent(in):: mat(:,:)
            ! output
            double precision:: score
            ! local
            integer:: n,j,k
            character(len=4):: subName='cksy'
            ! code
            call checkin(subName)
            n = size(mat,dim=1)
            if (size(mat,dim=2).ne.n) then
                score = -9999d0
                return
            endif
            score = 0d0
            do j = 1,n
                do k = j+1,n
                    if ((mat(j,k).eq.0d0).and.(mat(k,j).eq.0d0)) cycle       ! both are zero
                    if (mat(j,k)*mat(k,j).ne.0d0) then      ! neither is zero
                        score = score+((mat(j,k)-mat(k,j))/dmax1(dabs(mat(j,k)),dabs(mat(k,j))))**2
                    else                                    ! either is zero
                        ! score = score+1d0
                        score = score+dabs(mat(j,k)+mat(k,j))
                    endif
                enddo
            enddo
            score = dsqrt(score/((1+n)*n/2))
            call checkout(subName)
            return
        end function checkSymmetry

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        function Trace(mat) result(tr)
            ! 2017/09/06
            use ErrorHandler
            implicit none
            ! input
            double precision,intent(in):: mat(:,:)
            ! output
            double precision:: tr
            ! local
            integer:: n,j
            character(len=4):: subName='trac'
            ! code
            call checkin(subName)
            n = size(mat,dim=1)
            if (size(mat,dim=2).ne.n) then
                call RaiseError(1005,subName)
                return
            endif
            tr = 0d0
            do j = 1,n
                tr = tr+mat(j,j)
            enddo
            call checkout(subName)
            return
        end function Trace

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        function CholeskyDecomposition(mat,LU) result(matTri)
            ! 2017/09/07
            ! the output triangle is by default lower
            use ErrorHandler
            use String
            implicit none
            ! input
            double precision,intent(in):: mat(:,:)
            character(len=*),intent(in),optional:: LU
            ! output
            double precision,allocatable:: matTri(:,:)
            ! local
            integer:: n,j
            double precision,allocatable:: A(:,:)
            logical:: isLower
            character(len=4):: subName='chol'
            ! code
            call checkin(subName)
            n = size(mat,dim=1)
            if (n.ne.size(mat,dim=2)) then
                call RaiseError(1005,subName)
                return
            endif
            isLower = .true.
            if (present(LU)) then
                select case (lower(LU))
                    case ('u','up','upper')
                        isLower = .false.
                    case ('l','lo','low','lower')
                        isLower = .true.
                    case default
                        call RaiseWarning('Unknown LU flag.',subName)
                        isLower = .true.
                end select
            endif
            if (allocated(matTri)) then
                if ((n.ne.size(matTri,dim=1)).or.(n.ne.size(matTri,dim=2))) then
                    call RaiseWarning('Output matrix has unexpected dimensions. ',subName)
                    deallocate(matTri)
                    allocate(matTri(n,n))
                endif
            else
                allocate(matTri(n,n))
            endif
            allocate(A(n,n))
            matTri = 0d0
            A = mat
            ! begin loop
            do j = 1, n
                ! perform diagonal component
                A(j,j) = sqrt(A(j,j) - dot_product(A(j,1:j-1),A(j,1:j-1)))
                ! perform off-diagonal component
                if (j.lt.n) A(j+1:n,j) = (A(j+1:n,j) - matmul(A(j+1:n,1:j-1),A(j,1:j-1)))/A(j,j)
            enddo
            if (isLower) then
                matTri = A
            else
                matTri = transpose(A)
            endif
            if (allocated(A)) deallocate(A)
            call checkout(subName)
            return
        end function CholeskyDecomposition


! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        function Cholesky(mat,LU,Normalize) result(matTri)
            !   2017/09/09
            !   Given a positive-definite symmetric matrix a(1:n,1:n), this routine constructs its Cholesky decomposition, A = L·LT .
            !   On input, only the upper triangle of a need be given; it is not modified.
            !   The Cholesky factor L is returned in the lower triangle of a, except for its diagonal elements which are returned in p(1:n).
            use ErrorHandler
            use String
            implicit none
            ! input
            double precision,intent(in):: mat(:,:)
            character(len=*),optional,intent(in):: LU
            logical,optional,intent(in):: Normalize
            ! output
            double precision,allocatable:: matTri(:,:)
            ! local
            integer:: n,i,j,k,q
            double precision:: sum,s
            double precision,allocatable:: a(:,:),p(:)
            logical:: isLower,isNorm
            character(len=4):: subName='cho2'
            ! code
            call checkin(subName)
            isLower = .true.
            if (present(LU)) then
                select case (lower(LU))
                    case ('u','up','upper')
                        isLower = .false.
                    case ('l','lo','low','lower')
                        isLower = .true.
                    case default
                        call RaiseWarning('Unknown LU flag.',subName)
                        isLower = .true.
                end select
            endif
            isNorm = present(Normalize)
            if (present(Normalize)) isNorm = Normalize

            n = size(mat,dim=1)
            if (n.ne.size(mat,dim=2)) then
                call RaiseError(1005,subName)
                return
            endif
            if (allocated(matTri)) deallocate(matTri)
            allocate(matTri(n,n),a(n,n),p(n))

            s = 1d0
            if (isNorm) s = Trace(mat)
            a = mat/s
            matTri = 0d0
            do i = 1, n         ! column-wise
                do j = i, n     ! row-wise
                    sum = a(i,j)
                    do k = i-1,1,-1
                        sum = sum - a(i,k)*a(j,k)
                    enddo
                    if (i.eq.j) then
                        if (sum.le.0d0) then
                            ! do k = 1,n
                            !     write(999,'(<n>d25.15)') mat(k,1:n)
                            ! enddo
                            ! print '(1x,a,f)','Trace = ',Trace(a)
                            ! print '(1x,a,f)','Corrupted post-transformation diagonal element squared = ',sum
                            ! print '(1x,a,i<idint(dlog10(dble(i)))+1>)','Check column #',i
                            call RaiseWarning('Cholesky decomposition failed due to column #'//int2str(i)//'.',subName)
                            call RaiseError(1016,subName)
                            return
                        endif
                        p(i) = dsqrt(sum)
                    else
                        a(j,i) = sum/p(i)
                    endif
                enddo
            enddo
            do i = 1,n
                matTri(i,1:i-1) = a(i,1:i-1)
                matTri(i,i) = p(i)
            enddo
            matTri = matTri*dsqrt(s)
            deallocate(a,p)
            if (.not.isLower) matTri = transpose(matTri)
            call checkout(subName)
            return
        end function Cholesky

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        function MatVec(m,vin) result(vout)
            ! 2018/12/04
            ! Returns the vector of matrix x vector
            use ErrorHandler
            implicit none
            ! input
            double precision,intent(in):: m(:,:),vin(:)
            ! out
            double precision,allocatable:: vout(:)
            ! local
            integer:: i,j
            character(len=4):: subName='mxv'
            ! code
            call checkin(subName)
            if (.not.allocated(vout)) allocate(vout(size(vin)))
            if (size(m,1).ne.size(m,2).or.size(m,1).ne.size(vin).or.size(m,1).ne.size(vout)) then
                call RaiseError(0016,subName)
                return
            endif
            do i = 1,size(m,1)
                vout(i) = 0
                do j = 1,size(m,2)
                    vout(i) = vout(i)+m(i,j)*vin(j)
                enddo
            enddo
            call checkout(subName)
            return
        end function MatVec


! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>






end module MathsMatrix
