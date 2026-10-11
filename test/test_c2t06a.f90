! test_c2t06a.f90
! Test program for iau_C2T06A subroutine

program test_c2t06a
use module_space

    
    implicit none
    
    ! Declare variables
    double precision :: tta, ttb    ! TT as 2-part Julian Date
    double precision :: uta, utb    ! UT1 as 2-part Julian Date
    double precision :: taa, tab    ! TAI as 2-part Julian Date
    double precision :: xp_tmp, yp_tmp      ! Pole coordinates (radians)
    double precision :: rc2t(3,3)   ! Celestial-to-terrestrial matrix
    
    ! For date conversion
    integer :: iy, im, id           ! Year, month, day
    double precision :: djm0, djm   ! MJD parts
    double precision :: DUT1_tmp
    integer :: jstat                ! Status
    
    ! For matrix verification
    double precision :: det, trace
    integer :: i, j
    
    ! For coordinate transformation
    double precision :: celestial_vec(3)   ! Vector in celestial coordinates
    double precision :: terrestrial_vec(3) ! Vector in terrestrial coordinates
    
    real*8 :: HG(3, 3)

    
    ! Explicit interfaces for SOFA subroutines.
    ! The dummy names and the ABSENCE of INTENT must match sofa/*.for exactly:
    ! Intel's /warn:interfaces compares these blocks against the interfaces it
    ! generates from the SOFA sources and fails with error #8000 otherwise.
    ! Argument roles are 1-5 in, 6 out for the time-scale routines; the last
    ! argument of each (j) is the status.
    interface
        subroutine iau_CAL2JD(iy, im, id, djm0, djm, j)
            integer :: iy, im, id
            double precision :: djm0, djm
            integer :: j
        end subroutine iau_CAL2JD

        subroutine iau_UTCUT1(utc1, utc2, dut1, ut11, ut12, j)
            double precision :: utc1, utc2, dut1, ut11, ut12
            integer :: j
        end subroutine iau_UTCUT1

        subroutine iau_UTCTAI(utc1, utc2, tai1, tai2, j)
            double precision :: utc1, utc2, tai1, tai2
            integer :: j
        end subroutine iau_UTCTAI

        subroutine iau_TAITT(tai1, tai2, tt1, tt2, j)
            double precision :: tai1, tai2, tt1, tt2
            integer :: j
        end subroutine iau_TAITT

        subroutine iau_C2T06A(tta, ttb, uta, utb, xp, yp, rc2t)
            double precision :: tta, ttb, uta, utb, xp, yp
            double precision :: rc2t(3, 3)
        end subroutine iau_C2T06A
    end interface
    
    ! Print header
    write(*,*) '========================================'
    write(*,*) 'Test program for iau_C2T06A'
    write(*,*) '========================================'
    write(*,*)
    
    ! Set test date: 2025-01-01 00:00:00 UTC
    iy = 2025
    im = 1
    id = 1
    
    ! Convert calendar date to MJD for TT (Terrestrial Time)
    call iau_CAL2JD(iy, im, id, djm0, djm, jstat)
    if (jstat /= 0) then
        write(*,*) 'Error: Failed to convert calendar date to JD'
        stop
    endif
    
    ! Set pole coordinates (typ_tmpical values in radians)
    ! xp_tmp: -0.140 arcseconds = -0.140 * (π/648000) radians
    ! yp_tmp: 0.300 arcseconds = 0.300 * (π/648000) radians
    xp_tmp = 6.984371334205433D-007
    yp_tmp = 1.479205326167545D-006
    DUT1_tmp = 4.626730000000000D-002 ! seconds
    
    write(*, *) "xp_tmp, yp_tmp, dut1_tmp"
    write(*, *) xp_tmp, yp_tmp, dut1_tmp

    call iau_UTCUT1(djm0, djm, DUT1_tmp, uta, utb, jstat)

    call iau_UTCTAI(djm0, djm, taa, tab, jstat)
    call iau_TAITT(taa, tab, tta, ttb, jstat)

    
    ! Call the C2T06A function
    call iau_C2T06A(tta, ttb, uta, utb, xp_tmp, yp_tmp, rc2t)
    
    ! Print the resulting matrix
    write(*,*) 'Celestial-to-terrestrial matrix RC2T:'
    write(*,*) '--------------------------------------'
    do i = 1, 3
        write(*,'(A,I1,A,3F20.15)') ' Row ', i, ': ', (rc2t(i,j), j=1,3)
    end do
    write(*,*)
    
    ! Verify matrix properties
    ! 1. Check if it's approximately orthogonal (R * R^T ≈ I)
    write(*,*) 'Matrix verification:'
    write(*,*) '-------------------'
    
    ! Calculate determinant (should be close to 1 for rotation matrix)
    det = rc2t(1,1)*(rc2t(2,2)*rc2t(3,3) - rc2t(2,3)*rc2t(3,2)) &
        - rc2t(1,2)*(rc2t(2,1)*rc2t(3,3) - rc2t(2,3)*rc2t(3,1)) &
        + rc2t(1,3)*(rc2t(2,1)*rc2t(3,2) - rc2t(2,2)*rc2t(3,1))
    
    write(*,'(A,F20.15)') ' Determinant: ', det
    write(*,'(A,F20.15)') ' Deviation from 1: ', abs(det - 1.0d0)
    
    ! 2. Check trace
    trace = rc2t(1,1) + rc2t(2,2) + rc2t(3,3)
    write(*,'(A,F20.15)') ' Trace: ', trace
    
    ! 3. Check orthogonality: R * R^T should be identity matrix
    write(*,*) 'R * R^T (should be close to identity):'
    do i = 1, 3
        write(*,'(A,I1,A,3F12.8)') ' Row ', i, ': ', &
            dot_product(rc2t(i,:), rc2t(1,:)), &
            dot_product(rc2t(i,:), rc2t(2,:)), &
            dot_product(rc2t(i,:), rc2t(3,:))
    end do
    
    write(*,*)
    write(*,*) '========================================'
    write(*,*) '3D Coordinate Transformation Test'
    write(*,*) '========================================'
    write(*,*)
    
    ! Set a test vector in celestial coordinates
    ! Example: A unit vector pointing to a celestial position
    ! Right Ascension = 0 hours, Declination = 0 degrees
    ! In Cartesian coordinates: [1, 0, 0] (pointing to vernal equinox)
    celestial_vec(1) = 4d4
    celestial_vec(2) = 2d3
    celestial_vec(3) = 1d3
    
    write(*,*) 'Test vector in celestial coordinates:'
    write(*,'(A,3F15.6)') '  Vector = ', celestial_vec
    write(*,*)
    
    ! Transform the vector to terrestrial coordinates
    ! terrestrial_vec = RC2T * celestial_vec
    do i = 1, 3
        terrestrial_vec(i) = dot_product(rc2t(i,:), celestial_vec)
    end do
    
    write(*,*) 'Transformed vector in terrestrial coordinates:'
    write(*,'(A,3F15.6)') '  Vector = ', terrestrial_vec
    write(*,*)
    
    ! Verify that the length is preserved (should be 1.0)
    write(*,*) 'Verification:'
    write(*,'(A,F20.6)') '  Original length: ', &
        sqrt(dot_product(celestial_vec, celestial_vec))
    write(*,'(A,F20.6)') '  Transformed length: ', &
        sqrt(dot_product(terrestrial_vec, terrestrial_vec))
    write(*,'(A,F20.6)') '  Length difference: ', &
        abs(sqrt(dot_product(celestial_vec, celestial_vec)) - &
            sqrt(dot_product(terrestrial_vec, terrestrial_vec)))
    
    write(*,*)
    write(*,*) '========================================'
    write(*,*) 'Test completed successfully!'
    write(*,*) '========================================'
    
    call HG_get_sofa_06a((djm0+djm-2400000.5d0), HG)
    terrestrial_vec = matmul(HG, celestial_vec)
    write(*,*) 'Transformed vector in terrestrial coordinates:'
    write(*,'(A,3F15.6)') '  Vector = ', terrestrial_vec
    write(*,*)
    
end program test_c2t06a