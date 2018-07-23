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

  real :: r(idft,jmax+1)
  double precision :: spar(3*jmax/2+2)
  type(dfti_descriptor), pointer :: handle

  ! ymask used to weight forcing toward center latitdue
  do j = 1,jmax+1
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
  ! eddy streamfunction psi_1 and psi_2 and eddy PV vort_1 and vort_2
  ! initialized to zero; eddies introduced later through stochastic forcing
  psi_1  = zero
  vort_1 = zero
  psi_2  = zero
  vort_2 = zero

  !    ---- Basic state jet ----
  if (init_jet) then
    do j = 1,jtrunc/2
      jj = 2*j-1
      ell = el * float(jj)
      aap = ((-1)**(j+1))*exp(-ell*ell*sigma_jet*sigma_jet)
      do i = 2,9
        vort_1(i,2*j,1) = one_r*amp_jet*aap/float(jtrunc/2)
        vort_2(i,2*j,1) = one_r*amp_jet*aap/float(jtrunc/2)
      enddo
    enddo
  endif

  !    ---- Random initial noise in lower layer ----
  ! necessary to get eddies going, otherwise will just stay symmetric forever
  if (ll_seed_on) then
    call random_number(r)
    vort_2(:,:,1)= vort_2(:,:,1) + (r-0.5) * ll_seed_amp
  endif 
  vort_1(:,:,2) = vort_1(:,:,1)
  vort_1(:,:,3) = vort_1(:,:,1)
  vort_2(:,:,2) = vort_2(:,:,1)
  vort_2(:,:,3) = vort_2(:,:,1)
  !  test 
  !do i = 6,6
  !   do j = 2,2
  !     vort_1(i,j,:) = one_r*amp
  !     vort_2(i,j,:) = -one_r*amp
  !   enddo
  !enddo
  return
end subroutine
end module

