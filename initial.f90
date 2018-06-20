! ***** initial **********
!       Initializes the flow field

module INITIAL
contains

        subroutine initial1

        use GLOBAL_VARIABLES
        use mkl_dfti
        use mkl_trig_transforms

        implicit none

        integer n, k, tt_type
        integer ir, ipar(128)
        integer :: i,j,jj
        real :: y,yh,yph,phi,sech2
        real :: x,xph,rer,ell,aap 
        double precision :: spar(3*jmax/2+2)
        type(dfti_descriptor), pointer :: handle

! *** ymask will be used to weight forcing toward the center latitude ***
        do j = 1,jmax+1
         y = float(j-1)*dy-0.5*width
         ymask(j) = exp(-y*y/(sigma*sigma))
        enddo

! *** Blow qymean1 and qymean2 are the part of zonal-mean PV gradient that 
! excludes the uniform component, namely beta + u0/(rd*rd) and
! beta - u0/(rd*rd).  Since initially the total mean PV gradient is
! uniform and equal to beta + u0/(rd*rd) and beta - u0/(rd*rd), qymean1 and
! qymean2 are set to zero. Departure from the initial condition will be
! introduced later through eddy PV flux and forcing/damping.
        qymean1(:,:) = 0.
        qymean2(:,:) = 0.

! *** Below eddy streamfunction psi_1 and psi_2 and eddy PV vort_1 and vort_2
! will be initialized to zero. Eddy will be introduced later through stochastic
! forcing. 
   !    psi_1 = zero
   !    vort_1 = zero
   !    psi_2 = zero
   !    vort_2 = zero

   ! test 
       do i = 2,46
        do j = 2,46
          vort_1(i,j,:) = ur*1.e-10
          vort_2(i,j,:) = ur*1.e-10
        enddo
       enddo
        
      return
      end subroutine
    
end module

