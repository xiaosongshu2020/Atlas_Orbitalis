    module EarthRadiationPressure
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   Introduction
    !       This module's function is to calculate the perturbation
    !       of the Earth radiation pressure. Its algorithm is based
    !       on this paper -- "Earth Radiation Pressure Effects on 
    !       Satellites" P.C. Knocke, J, C. Ries, and B. D. Tapley, 
    !       in 1988.
    !   Auther  :   Shawn Hu
    !   Date    :   2022-10-19
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    
    real*8,  protected :: degPerGrid                            ! Angle per grid (deg)
    integer, protected :: nGridLong, nGridLati                  ! Grid number in longitude and latitude
    logical, protected :: isSphereGridSet = .False.             ! Whether the grid has been set
    
    real*8,  protected :: aGrid(-90:90)                         ! The albedo whose index is latitude
    real*8,  protected :: eGrid(-90:90)                         ! The emissivity whose index is latitude
    logical, protected :: isaeGridSet = .False.                 ! Whether the ae Grid has been set
    
    contains
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine setaeGrid(MJD)
        
        implicit none
        real*8, intent(in)  ::  MJD
        real*8  ::  a0 = 0.34d0, a1 = 0.00d0, a2 = 0.29d0
        real*8  ::  c0 = 0.00d0, c1 = 0.10d0, c2 = 0.00d0
        real*8  ::  e0 = 0.68d0, e1 = 0.00d0, e2 =-0.18d0
        real*8  ::  k0 = 0.00d0, k1 =-0.07d0, k2 = 0.00d0
        real*8  ::  pi      =   3.141592653589793d0
        real*8  ::  w       =   2d0*3.141592653589793d0/365.25d0
        real*8  ::  MJD0    =   44960d0                         ! 22 Dec 1981
        real*8  ::  phi, delta_t
        integer ::  i
        
        !-------------------------------------------------------
        
        do i = -90, 90
            phi     =   dble(i)*pi/180d0
            delta_t =   MJD - MJD0
            aGrid(i)=   a0 + (c0 + c1*dcos(w*delta_t) + &
                &       c2*dsin(w*delta_t))*dsin(phi) + &
                &       a2*0.5d0*(3d0*dsin(phi)**2 - 1d0)
            eGrid(i)=   e0 + (k0 + k1*dcos(w*delta_t) + &
                &       k2*dsin(w*delta_t))*dsin(phi) + &
                &       e2*0.5d0*(3d0*dsin(phi)**2 - 1d0)
        enddo
        
        isaeGridSet = .True.
        
        endsubroutine
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine setSphereGrid(height, degPerGrid1)
        
        implicit none
        real*8, intent(in)      ::  height                              ! Satellite height (km)
        real*8, intent(in), optional  ::  degPerGrid1                   ! set angle per grid (deg)
        real*8  ::  angle
        
        !-------------------------------------------------------------
        
        if (present(degPerGrid1)) then                                  ! Spliding Earth's surface using input
            if (180d0/degPerGrid1 == dble(int(180d0/degPerGrid1))) then ! If 180 can be divided by degPerGrid1
                angle   =   degPerGrid1
            else
                write(*, *) "Input error -- degPerGrid1 should be a ", &
                    & "factor of 180"
            endif
        else                                                            ! Spliding Earth's surface using the automatically selected angle
            if (height > 2000d0) then
                angle   =   10d0
            elseif (height > 1500d0) then
                angle   =   8d0
            elseif (height > 1000d0) then
                angle   =   6d0
            elseif (height > 800d0) then
                angle   =   4d0
            elseif (height > 500d0) then
                angle   =   3d0
            elseif (height > 300d0) then
                angle   =   2d0
            else
                angle   =   1d0
            endif
        endif
        
        nGridLong       =   int(360d0/angle)
        nGridLati       =   int(180d0/angle)
        isSphereGridset =   .True.
        
        endsubroutine
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calEarthRadiationPressure(t, rSun, rSat0, S_m, reflectivity, F)
        
        implicit none
        real*8, intent(in)  ::  t, rSat0(3), rSun(3), S_m, reflectivity
        real*8  ::  F(3), dF(3)
        integer ::  i, j
        real*8  ::  length_unit =   6378.1363d0                         ! km
        real*8  ::  mass_unit   =   5.9742272407074334d24               ! kg
        real*8  ::  pi          =   3.141592653589793d0
        real*8  ::  Es          =   0.3169d-17                          ! Dimension is [F]/[L]^2, normalized units
        real*8  ::  da, rSatNorm, rSphere
        real*8  ::  alphaA, betaA, rA(3), rSat(3)
        real*8  ::  areaA, cosThetaS, cosAlpha, a, e
        real*8  ::  tau
        !------------------------------------------------------------
        
        if (.NOT.isaeGridSet)       call    setaeGrid(t)                ! Check whether the ae grid has been set
        if (.NOT.isSphereGridSet)   call    setSphereGrid((norm2(rSat0)-1d0)*length_unit)
                                                                        ! Check whether the Sphere Grid has been set
        F           =   0d0
        da          =   360d0/nGridLong*pi/180d0                        ! Angel per grid
        rSphere     =   6397.38533d0/length_unit                        ! TOA(Top of atmosphere) is defined as a sphere with 
        do i = 1, nGridLong                                             ! a radius of 6397.38533d0 km (30 km above the ground).
            do j = 1, nGridlati                                         ! For every single area, cal F and add them up
                
                alphaA  =   dble(i)*da - 0.5d0*da
                betaA   =   dble(j)*da - 0.5d0*da - pi/2d0
                rA      =   (/dcos(betaA)*dcos(alphaA), dcos(betaA)*dsin(alphaA), &
                    &       dsin(betaA)/)*rSphere
                rSat    =   rSat0 - rA                                  ! rSat0 is satellite's position relative to geocenter
                                                                        ! rSat is satellite's position relative to area A
                
                if (dot_product(rA, rSat) <= 0d0) cycle                 ! If the erea A is invisible.
                
                areaA       =   rSphere**2*da*da*dcos(betaA)            ! The area of A
                cosThetaS   =   dot_product(rA, rSun)/norm2(rA)/norm2(rSun)
                cosAlpha    =   dot_product(rA, rSat)/norm2(rA)/norm2(rSat)
                a           =   aGrid(int(betaA*180/pi))                ! get a and e from aeGrid
                e           =   eGrid(int(betaA*180/pi))
                
                if (cosThetaS <= 0d0) then                              ! If the area A is darkness tau = 0 and None albedo
                    tau     =   0d0
                else
                    tau     =   1d0
                endif
                
                dF          =   (tau*a*Es*cosThetaS + &
                    &           e*Es/4d0)*cosAlpha*areaA/pi/norm2(rSat)**3*rSat
                F           =   F + dF
            enddo
        enddo
        
        F = F * S_m * (1d0 + reflectivity) * mass_unit / (length_unit * 1000d0)**2
        
        endsubroutine
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end module