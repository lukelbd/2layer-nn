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

        real :: r(imax,jmax+1)
        double precision :: spar(3*jmax/2+2)
        type(dfti_descriptor), pointer :: handle


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
        psi_1 = zero
        vort_1 = zero
        psi_2 = zero
        vort_2 = zero


!   *** flow (uniform PV perturbation, symmetric mode only) ****
    if(init_jet) then

      do j = 1,nmax/2
         jj = 2*j-1
         ell = el * float(jj)
         aap = ((-1)**(j+1))*exp(-ell*ell*sigma*sigma)
         do i = 2,9
            vort_1(i,2*j,1) = ur*init_jet_amp*aap/float(nmax/2)
            vort_2(i,2*j,1) = ur*init_jet_amp*aap/float(nmax/2)
         enddo
      enddo

    endif
!   *** flow add random initial noise in upper layer ****
! 
    if(random_seed) then
        CALL RANDOM_NUMBER(r)
        vort_2(:,:,1)= vort_2(:,:,1) + (r-0.5) *rand_seed_amp
    endif 

    vort_1(:,:,2) = vort_1(:,:,1)
    vort_1(:,:,3) = vort_1(:,:,1)

    vort_2(:,:,2) = vort_2(:,:,1)
    vort_2(:,:,3) = vort_2(:,:,1)


  !  test 
  !do i = 6,6
  !   do j = 2,2
  !     vort_1(i,j,:) = ur*amp
  !     vort_2(i,j,:) = -ur*amp
  !   enddo
  !enddo
        
      return
      end subroutine
    
end module

