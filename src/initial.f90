!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial flow field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module initial
contains
subroutine init

  use global_variables
  use mkl_dfti
  use mkl_trig_transforms

  implicit none

  integer n, k, tt_type
  integer ir, ipar(128)
  integer :: i,j,jj
  real :: y,yh,yph,phi,sech2
  real :: x,xph,rer,ell,aap 

  real :: r(idft,jmax)
  double precision :: spar(3*(jmax-1)/2+2)
  type(dfti_descriptor), pointer :: handle

  ! ymask used to weight forcing toward center latitdue
  do j = 1,jmax
    y = float(j-1)*dy-0.5*width
    ymask(j) = exp(-y*y/(sigma_i*sigma_i))
  enddo

  !    ---- Zonal mean pv ----
  ! qymean1 and qymean2 are the part of zonal-mean PV gradient that 
  ! excludes the uniform component [i.e. beta + u0/(rd*rd) and
  ! beta - u0/(rd*rd)]. Since initially the total mean PV gradient is
  ! uniform and equal to beta + u0/(rd*rd) and beta - u0/(rd*rd), qymean1 and
  ! qymean2 are set to zero. Departure from the initial condition will be
  ! introduced later through eddy PV flux and forcing/damping.
  qymean1 = 0.
  qymean2 = 0.

  !    ---- Vorticity and streamfunc ----
  ! eddy streamfunction psi1_c and psi2_c and eddy PV q1_c and q2_c
  ! initialized to zero; eddies introduced later through stochastic forcing
  psi1_c  = zero
  q1_c = zero
  psi2_c  = zero
  q2_c = zero

  !    ---- Basic state jet ----
  if (init_jet) then
    do j = 1,jtrunc/2
      jj = 2*j-1
      ell = el * float(jj)
      aap = ((-1)**(j+1))*exp(-ell*ell*sigma_jet*sigma_jet)
      do i = 2,9
        q1_c(i,2*j,1) = one_r*amp_jet*aap/float(jtrunc/2)
        q2_c(i,2*j,1) = one_r*amp_jet*aap/float(jtrunc/2)
      enddo
    enddo
  endif

  !    ---- Random initial noise in lower layer ----
  ! necessary to get eddies going, otherwise will just stay symmetric forever
  if (ll_seed_on) then
    call random_number(r)
    q2_c(:,:,1)= q2_c(:,:,1) + (r-0.5) * ll_seed_amp
  endif 
  q1_c(:,:,2) = q1_c(:,:,1)
  q1_c(:,:,3) = q1_c(:,:,1)
  q2_c(:,:,2) = q2_c(:,:,1)
  q2_c(:,:,3) = q2_c(:,:,1)
  !  test 
  !do i = 6,6
  !   do j = 2,2
  !     q1_c(i,j,:) = one_r*amp
  !     q2_c(i,j,:) = -one_r*amp
  !   enddo
  !enddo
  return
end subroutine
end module

