!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module calculates a bunch of diagnostics required for prognostics.f90
! including advection and PV injection forcing. Also returns cartesian x-y
! data suitable for output.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note that all q variables contain *dq/dy*, not q itself.
! Walkthrough of steps:
!  1. Invert the dq'/dy and dabar/dy equations to obtain wind, streamfunction, et cetera
!     from the fully spectral data.
!  2. Perform *backward* trigonometric (sines only or cosines only) transforms
!     on the y-dimension of the inverted data.
!  3. Perform *backward* Fourier transform on the x-dimension of the inverted
!     data.
!  4. Now we have Cartesian data; calculate PV injection, store some values
!     in arrays suitable for output.
!  5. Calculate advection terms (for q') and flux convergence terms (for qbar).
!     Requires another *forward* Fourier transform of the u' data and
!     *forward* trigonometric transforms?
!     This is done because we need to calculate ***products*** of two terms
!     in Cartesian coordinates, but calculate ***derivs*** in spectral space.
!     Why do we perform forward transform here? Didn't we alredy have that data?
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For more information on MKL utilities:
!  1. run 'module load mkl'
!  2. go to $MKLROOT/include, and check out the
!     mkl_dfti.f90 and mkl_trig_transforms.f90 files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Warning:
! Variables are often not what they seem. The 'tt' variables can be in cartesian
! y coordinates or spectral sine/cosine y coordinates.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module diagnostics
contains
subroutine diag(hx,hy)
  use global_variables
  use mkl_dfti            ! this one just includes some boilerplate stuff; is the 'interface'
  use mkl_trig_transforms ! this one contains functions we actually use
  use ftm

  implicit none

  ! Random stuff
  type(dfti_descriptor), pointer :: hx,hy,handle
  integer ret,dsadata,shape(2)
  integer tt_type ! trigonometric transformation type for mkl_trig_transforms
                  ! from file: 0 == sine, 1 == cosine
  ! Other stuff
  integer ir, ipar(128)
  double precision :: spar(3*(jmax-1)/2+2)
  real :: y,rkk,ell,fac,as,umax
  integer :: i,j,wcos,wsin
  complex :: pb, pc ! temporary scalars during PV inversion

  ! Arrays storing fully spectral data (real/imaginary coeffs of
  ! x-Fourier decomposition, and trig transform in y)
  ! Sometimes u/v will hold wind, sometimes flux
  real :: u1_r(idft,jmax),  u1_i(idft,jmax),  &
          v1_r(idft,jmax),  v1_i(idft,jmax),  &
          q1_r(idft,jmax),  q1_i(idft,jmax),  &
          qx1_r(idft,jmax), qx1_i(idft,jmax), &
          qy1_r(idft,jmax), qy1_i(idft,jmax), &
          p1_r(idft,jmax),  p1_i(idft,jmax), &
          adv1_r(idft,jmax), adv1_i(idft,jmax), &
          f1_r(idft,jmax),  f1_i(idft,jmax)
  real :: u2_r(idft,jmax),  u2_i(idft,jmax),  &
          v2_r(idft,jmax),  v2_i(idft,jmax),  &
          q2_r(idft,jmax),  q2_i(idft,jmax),  &
          qx2_r(idft,jmax), qx2_i(idft,jmax), &
          qy2_r(idft,jmax), qy2_i(idft,jmax), &
          p2_r(idft,jmax),  p2_i(idft,jmax), &
          adv2_r(idft,jmax), adv2_i(idft,jmax)

  ! Variables storing x-Fourier coefficients, for iterating over i
  ! Stores data cartesian in y *or* trig transform (sine/cosine)
  ! The 'tt' stands for trigonometric transform
  double precision :: v1_r_tt(jmax), v2_r_tt(jmax), v1_i_tt(jmax), v2_i_tt(jmax), &
    u1_r_tt(jmax),   u1_i_tt(jmax),  u2_r_tt(jmax),  u2_i_tt(jmax),  &
    q1_r_tt(jmax),   q2_r_tt(jmax),  q1_i_tt(jmax),  q2_i_tt(jmax),  &
    qx1_r_tt(jmax),  qx2_r_tt(jmax), qx1_i_tt(jmax), qx2_i_tt(jmax), &
    qy1_r_tt(jmax),  qy2_r_tt(jmax), qy1_i_tt(jmax), qy2_i_tt(jmax), &
    p1_r_tt(jmax),   p2_r_tt(jmax),  p1_i_tt(jmax),  p2_i_tt(jmax),  &
    adv1_r_tt(jmax), adv1_i_tt(jmax), adv2_r_tt(jmax), adv2_i_tt(jmax), &
    f1_r_tt(jmax),   f1_i_tt(jmax)

  ! Zonal mean/zonal flux values
  ! Stored data cartesian in y *or* trig transform (sine/cosine)
  double precision :: ub_tt(jmax), uc_tt(jmax), &
    qflux1_tt(jmax), qybar1_tt(jmax), qbar1_tt(jmax), ubar1_tt(jmax), &
    qflux2_tt(jmax), qybar2_tt(jmax), qbar2_tt(jmax), ubar2_tt(jmax)

  ! Truncated Fourier coefficients, for iterating over j
  double complex :: u1_dft(idft),  u2_dft(idft), &
    v1_dft(idft),   v2_dft(idft),   &
    q1_dft(idft),   q2_dft(idft),   &
    qx1_dft(idft),  qx2_dft(idft),  &
    qy1_dft(idft),  qy2_dft(idft),  &
    p1_dft(idft),   p2_dft(idft),   &
    adv1_dft(idft), adv2_dft(idft), &
    f1_dft(idft)

  ! Inverse transforms of the above values
  ! This data is always cartesian
  double precision :: u1_cart(imax),  u2_cart(imax), &
    v1_cart(imax),   v2_cart(imax),   &
    q1_cart(imax),   q2_cart(imax),   &
    vq1_cart(imax),  vq2_cart(imax),  &
    qx1_cart(imax),  qx2_cart(imax),  &
    qy1_cart(imax),  qy2_cart(imax),  &
    p1_cart(imax),   p2_cart(imax),   &
    adv1_cart(imax), adv2_cart(imax), &
    f1_cart(imax)

  ! For PV injection
  integer,allocatable :: iseed(:)
  integer             :: isize,idate(8)
  real                :: anglex,angley,amp_rand

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Invert d/dy zonal mean PV equation to get zonal mean u
  ! Also integrate d/dy zonal mean PV to get zonal mean PV
  ! * The umean are needed for radiation forcing in prognostics.f90, which is why
  !   we create the seemingly identical variables 'umean1/ubar1_tt' and 'umean2/ubar2_tt'
  ! * The latter will undergo in-place sine transformation in this module
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  umean1(:) = 0.
  umean2(:) = 0.
  ub_tt(:) = 0.
  uc_tt(:) = 0.
  qbar1_tt(:) = 0.
  qbar2_tt(:) = 0.
  ubar1_tt(:) = 0.
  ubar2_tt(:) = 0.
  do j = 2,jtrunc
    ell = el*float(j-1) ! angular wavenumber
    ub_tt(j) = (qymean1(j,3)+qymean2(j,3))/(ell**2)
    uc_tt(j) = (qymean1(j,3)-qymean2(j,3))/(ell**2 + (2./(rd*rd)))
    umean1(j) = (ub_tt(j)+uc_tt(j))*0.5
    umean2(j) = (ub_tt(j)-uc_tt(j))*0.5
    qbar1_tt(j) = qymean1(j,3)/ell ! d/dx(f(t)) = wavenum*F(wavenum) for sine transform (note wavenum is angular here)
    qbar2_tt(j) = qymean2(j,3)/ell
    ubar1_tt(j) = umean1(j)
    ubar2_tt(j) = umean2(j)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Invert PV equation to get streamfunction
  ! Also obtain u/v, and q derivatives in x/y (derivative in spectral coordinates)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transform eddy into physical space
  ! Used to just assign (1,:) and (:,1) zero but am fairly certain
  ! this was minor bug, and we were transforming arrays populated by 'empty'
  ! values above the trunc number (generally random, extremely small numbers)
  psi1_c(:,:) = zero
  psi2_c(:,:) = zero
  u1_c(:,:)   = zero
  u2_c(:,:)   = zero
  v1_c(:,:)   = zero
  v2_c(:,:)   = zero
  qx1_c(:,:)  = zero
  qx2_c(:,:)  = zero
  qy1_c(:,:)  = zero
  qy2_c(:,:)  = zero
  psi1_c(:,:) = zero
  psi2_c(:,:) = zero
  u1_c(:,:)   = zero
  u2_c(:,:)   = zero
  v1_c(:,:)   = zero
  v2_c(:,:)   = zero
  qx1_c(:,:)  = zero
  qx2_c(:,:)  = zero
  qy1_c(:,:)  = zero
  qy2_c(:,:)  = zero
  do j = 2,jtrunc
    ell = el*float(j-1)
    do i = 2,itrunc
      rkk = rk*float(i-1)
      pb  = -(q1_c(i,j,3)+q2_c(i,j,3))/(rkk**2 + ell**2)
      pc  = -(q1_c(i,j,3)-q2_c(i,j,3))/(rkk**2 + ell**2  + (2./(rd*rd)))
      psi1_c(i,j) = 0.5*(pb+pc)
      psi2_c(i,j) = 0.5*(pb-pc)
      u1_c(i,j) = -ell*psi1_c(i,j)
      u2_c(i,j) = -ell*psi2_c(i,j)
      v1_c(i,j) = one_i*rkk*psi1_c(i,j)
      v2_c(i,j) = one_i*rkk*psi2_c(i,j)
      qx1_c(i,j) = one_i*rkk*q1_c(i,j,3)
      qx2_c(i,j) = one_i*rkk*q2_c(i,j,3)
      qy1_c(i,j) = ell*q1_c(i,j,3)
      qy2_c(i,j) = ell*q2_c(i,j,3)
    enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transform zonal mean to physical space 
  ! * Apply inverse trigonometric transforms
  ! * Seems to perform ****in place**** modification of array values
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Inverse sine transforms
  ! Suitable if edges are always small/zero
  qybar1_tt(:) = qymean1(:,3) ! below functions need a 1D series
  qybar2_tt(:) = qymean2(:,3)
  tt_type=0
  call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
  call d_commit_trig_transform(ubar1_tt,handle,ipar,spar,ir)
  call d_backward_trig_transform(ubar1_tt,handle,ipar,spar,ir)
  call free_trig_transform(handle,ipar,ir)

  call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
  call d_commit_trig_transform(ubar2_tt,handle,ipar,spar,ir)
  call d_backward_trig_transform(ubar2_tt,handle,ipar,spar,ir)
  call free_trig_transform(handle,ipar,ir)

  call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
  call d_commit_trig_transform(qybar1_tt,handle,ipar,spar,ir)
  call d_backward_trig_transform(qybar1_tt,handle,ipar,spar,ir)
  call free_trig_transform(handle,ipar,ir)

  call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
  call d_commit_trig_transform(qybar2_tt,handle,ipar,spar,ir)
  call d_backward_trig_transform(qybar2_tt,handle,ipar,spar,ir)
  call free_trig_transform(handle,ipar,ir)

  ! Inverse cosine transforms
  ! Suitable if edges are always non-zero
  tt_type=1
  call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
  call d_commit_trig_transform(qbar1_tt,handle,ipar,spar,ir)
  call d_backward_trig_transform(qbar1_tt,handle,ipar,spar,ir)
  call free_trig_transform(handle,ipar,ir)

  call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
  call d_commit_trig_transform(qbar2_tt,handle,ipar,spar,ir)
  call d_backward_trig_transform(qbar2_tt,handle,ipar,spar,ir)
  call free_trig_transform(handle,ipar,ir)

  ! Adds in the background component by integrating the dq/dy
  ! equation from lower to upper boundary.
  do j = 1,jmax
    y = dy*float(j-1)-0.5*width 
    qbar1_tt(j) = qbar1_tt(j)+(beta+u0/(rd*rd))*y
    qbar2_tt(j) = qbar2_tt(j)+(beta-u0/(rd*rd))*y
  enddo
  if(mod(t,10000).eq.1) then
    do j = 1,jmax
      write(*,*) j,'u1:',ubar1_tt(j)+u0,' qbar1:',qbar1_tt(j)
    enddo
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Perform inverse y-sine/cosine transformations
  ! of the x-direction Fourier coefficients for 2D params
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i = 1,idft
    ! First prepare 1D vectors containing the complex x-direction
    ! Fourier coefficients. Will supply these to inverse sine/cosine transform functions.
    u1_r_tt(:)  = real(u1_c(i,:))
    u1_i_tt(:)  = aimag(u1_c(i,:))
    u2_r_tt(:)  = real(u2_c(i,:))
    u2_i_tt(:)  = aimag(u2_c(i,:))
    v1_r_tt(:)  = real(v1_c(i,:))
    v1_i_tt(:)  = aimag(v1_c(i,:))
    v2_r_tt(:)  = real(v2_c(i,:))
    v2_i_tt(:)  = aimag(v2_c(i,:))
    q1_r_tt(:)  = real(q1_c(i,:,3))
    q1_i_tt(:)  = aimag(q1_c(i,:,3))
    qx1_r_tt(:) = real(qx1_c(i,:))
    qx1_i_tt(:) = aimag(qx1_c(i,:))
    qy1_r_tt(:) = real(qy1_c(i,:))
    qy1_i_tt(:) = aimag(qy1_c(i,:))
    q2_r_tt(:)  = real(q2_c(i,:,3))
    q2_i_tt(:)  = aimag(q2_c(i,:,3))
    qx2_r_tt(:) = real(qx2_c(i,:))
    qx2_i_tt(:) = aimag(qx2_c(i,:))
    qy2_r_tt(:) = real(qy2_c(i,:))
    qy2_i_tt(:) = aimag(qy2_c(i,:))
    p1_r_tt(:)  = real(psi1_c(i,:))
    p1_i_tt(:)  = aimag(psi1_c(i,:))
    p2_r_tt(:)  = real(psi2_c(i,:))
    p2_i_tt(:)  = aimag(psi2_c(i,:))

    ! Inverse sine transforms
    ! Suitable if edges are always small/zero
    ! Note we are transforming the x-Fourier coefficients here
    tt_type=0
    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(v1_r_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(v1_r_tt,handle,ipar,spar,ir)
    v1_r(i,:) = v1_r_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(v1_i_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(v1_i_tt,handle,ipar,spar,ir)
    v1_i(i,:) = v1_i_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(q1_r_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(q1_r_tt,handle,ipar,spar,ir)
    q1_r(i,:) = q1_r_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(q1_i_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(q1_i_tt,handle,ipar,spar,ir)
    q1_i(i,:) = q1_i_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(qx1_r_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(qx1_r_tt,handle,ipar,spar,ir)
    qx1_r(i,:) = qx1_r_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(qx1_i_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(qx1_i_tt,handle,ipar,spar,ir)
    qx1_i(i,:) = qx1_i_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(p1_r_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(p1_r_tt,handle,ipar,spar,ir)
    p1_r(i,:) = p1_r_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(p1_i_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(p1_i_tt,handle,ipar,spar,ir)
    p1_i(i,:) = p1_i_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(v2_r_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(v2_r_tt,handle,ipar,spar,ir)
    v2_r(i,:) = v2_r_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(v2_i_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(v2_i_tt,handle,ipar,spar,ir)
    v2_i(i,:) = v2_i_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(q2_r_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(q2_r_tt,handle,ipar,spar,ir)
    q2_r(i,:) = q2_r_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(q2_i_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(q2_i_tt,handle,ipar,spar,ir)
    q2_i(i,:) = q2_i_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(qx2_r_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(qx2_r_tt,handle,ipar,spar,ir)
    qx2_r(i,:) = qx2_r_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(qx2_i_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(qx2_i_tt,handle,ipar,spar,ir)
    qx2_i(i,:) = qx2_i_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(p2_r_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(p2_r_tt,handle,ipar,spar,ir)
    p2_r(i,:) = p2_r_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(p2_i_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(p2_i_tt,handle,ipar,spar,ir)
    p2_i(i,:) = p2_i_tt(:)
    call free_trig_transform(handle,ipar,ir)

    ! Inverse cosine transforms
    ! Suitble if edges are always non-zero
    tt_type=1
    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(u1_r_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(u1_r_tt,handle,ipar,spar,ir)
    u1_r(i,:) = u1_r_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(u1_i_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(u1_i_tt,handle,ipar,spar,ir)
    u1_i(i,:) = u1_i_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(qy1_r_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(qy1_r_tt,handle,ipar,spar,ir)
    qy1_r(i,:) = qy1_r_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(qy1_i_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(qy1_i_tt,handle,ipar,spar,ir)
    qy1_i(i,:) = qy1_i_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(u2_r_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(u2_r_tt,handle,ipar,spar,ir)
    u2_r(i,:) = u2_r_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(u2_i_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(u2_i_tt,handle,ipar,spar,ir)
    u2_i(i,:) = u2_i_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(qy2_r_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(qy2_r_tt,handle,ipar,spar,ir)
    qy2_r(i,:) = qy2_r_tt(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(qy2_i_tt,handle,ipar,spar,ir)
    call d_backward_trig_transform(qy2_i_tt,handle,ipar,spar,ir)
    qy2_i(i,:) = qy2_i_tt(:)
    call free_trig_transform(handle,ipar,ir)
  enddo ! loop over idft

  do j = 1,jmax
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Inverse x-Fourier transforms
    ! * The y dimension now has *cartesian* units, not sine/cosine transform coeffs
    ! * Note the full Fourier transform equations take input/output arrays; they
    !   do *not* change arrays in place.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Prepare 1D vectors that will be passed to Fourier inversion functions
    ! Also apply truncation here
    v1_dft(:)  = zero
    v2_dft(:)  = zero
    u1_dft(:)  = zero
    u2_dft(:)  = zero
    q1_dft(:)  = zero
    q2_dft(:)  = zero
    qx1_dft(:) = zero
    qx2_dft(:) = zero
    qy1_dft(:) = zero
    qy2_dft(:) = zero
    p1_dft(:)  = zero
    p2_dft(:)  = zero
    v1_dft(2:itrunc)  = 0.5*(one_r*v1_r(2:itrunc,j)  + one_i*v1_i(2:itrunc,j))
    v2_dft(2:itrunc)  = 0.5*(one_r*v2_r(2:itrunc,j)  + one_i*v2_i(2:itrunc,j))
    u1_dft(2:itrunc)  = 0.5*(one_r*u1_r(2:itrunc,j)  + one_i*u1_i(2:itrunc,j))
    u2_dft(2:itrunc)  = 0.5*(one_r*u2_r(2:itrunc,j)  + one_i*u2_i(2:itrunc,j))
    q1_dft(2:itrunc)  = 0.5*(one_r*q1_r(2:itrunc,j)  + one_i*q1_i(2:itrunc,j))
    q2_dft(2:itrunc)  = 0.5*(one_r*q2_r(2:itrunc,j)  + one_i*q2_i(2:itrunc,j))
    qx1_dft(2:itrunc) = 0.5*(one_r*qx1_r(2:itrunc,j) + one_i*qx1_i(2:itrunc,j))
    qx2_dft(2:itrunc) = 0.5*(one_r*qx2_r(2:itrunc,j) + one_i*qx2_i(2:itrunc,j))
    qy1_dft(2:itrunc) = 0.5*(one_r*qy1_r(2:itrunc,j) + one_i*qy1_i(2:itrunc,j))
    qy2_dft(2:itrunc) = 0.5*(one_r*qy2_r(2:itrunc,j) + one_i*qy2_i(2:itrunc,j))
    p1_dft(2:itrunc)  = 0.5*(one_r*p1_r(2:itrunc,j)  + one_i*p1_i(2:itrunc,j))
    p2_dft(2:itrunc)  = 0.5*(one_r*p2_r(2:itrunc,j)  + one_i*p2_i(2:itrunc,j))
    v1_dft(itrunc)  = one_r*real(v1_dft(itrunc))
    v2_dft(itrunc)  = one_r*real(v2_dft(itrunc))
    u1_dft(itrunc)  = one_r*real(u1_dft(itrunc))
    u2_dft(itrunc)  = one_r*real(u2_dft(itrunc))
    q1_dft(itrunc)  = one_r*real(q1_dft(itrunc))
    q2_dft(itrunc)  = one_r*real(q2_dft(itrunc))
    qx1_dft(itrunc) = one_r*real(qx1_dft(itrunc))
    qx2_dft(itrunc) = one_r*real(qx2_dft(itrunc))
    qy1_dft(itrunc) = one_r*real(qy1_dft(itrunc))
    qy2_dft(itrunc) = one_r*real(qy2_dft(itrunc))
    p1_dft(itrunc)  = one_r*real(p1_dft(itrunc))
    p2_dft(itrunc)  = one_r*real(p2_dft(itrunc))

    ! Inverse Fourier Transform
    ! Note we ***write to global variables here***
    ! The scrfr means complex-to-real fourier transform
    as = 1.
    v1_cart(:)  = 0.
    v2_cart(:)  = 0.
    u1_cart(:)  = 0.
    u2_cart(:)  = 0.
    q1_cart(:)  = 0.
    q2_cart(:)  = 0.
    qx1_cart(:) = 0.
    qx2_cart(:) = 0.
    qy1_cart(:) = 0.
    qy2_cart(:) = 0.
    p1_cart(:)  = 0.
    p2_cart(:)  = 0.
    call scrft(v1_dft,v1_cart,imax,as,hx)
    call scrft(v2_dft,v2_cart,imax,as,hx)
    call scrft(u1_dft,u1_cart,imax,as,hx)
    call scrft(u2_dft,u2_cart,imax,as,hx)
    call scrft(q1_dft,q1_cart,imax,as,hx)
    call scrft(q2_dft,q2_cart,imax,as,hx)
    call scrft(qx1_dft,qx1_cart,imax,as,hx)
    call scrft(qx2_dft,qx2_cart,imax,as,hx)
    call scrft(qy1_dft,qy1_cart,imax,as,hx)
    call scrft(qy2_dft,qy2_cart,imax,as,hx)
    call scrft(p1_dft,p1_cart,imax,as,hx)
    call scrft(p2_dft,p2_cart,imax,as,hx)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Finally have the 2D real data in this block!!!
    ! Do various things to it
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Stochastic forcing (small-scale narrow band)
    ! Injects cosines and sines onto the raw 2D data, with fixed memory 0.5
    ! Note this injects the *forcing*
    ! * Must appear here, because we have a real ymask that limits appearence
    !   of PV perturbations into a narrow 'jet' band
    ! * Below can be used to have arbitrary memory; just replace where you see 0.5
    !   fxy1(i,j,2) = exp(-dt/tau_i)*fxy1(i,j,1)
    !   + (1.-exp(-dt/tau_i))*amp_i*ymask(j)*  &
    iseed=iseed*(idate(8)-500)
    call date_and_time(values=idate)
    call random_seed(size=isize)
    call random_seed(put=iseed)
    call random_number(anglex)
    call random_number(angley)
    do i = 1,imax
      f1_cart(i) = 0.5*fxy1(i,j,1) ! initialize with previous state
      do wcos = wmin_i,wmax_i
        rkk = rk*float(wcos-1)
        do wsin = wmin_i,wmax_i
          ell = el*float(wsin-1)
          call random_number(amp_rand)
          f1_cart(i) = f1_cart(i) + 0.5*amp_i*amp_rand*ymask(j)  &
            * sin(ell*float(j-1)/float(jmax-1)+angley)  &
            * cos(rkk*float(i-1)/float(imax)+anglex)
        enddo
      enddo
    enddo

    ! Write data to some convenient arrays that can be output in io.f90
    ! q', dq'/dx, dq'/dy fields
    qxy1(:,j) = q1_cart(:)
    qxy2(:,j) = q2_cart(:)
    qxx1(:,j) = qx1_cart(:)
    qxx2(:,j) = qx2_cart(:)
    qyy1(:,j) = qy1_cart(:)
    qyy2(:,j) = qy2_cart(:)
    ! v' field (always has zero mean)
    vxy1(:,j) = v1_cart(:)
    vxy2(:,j) = v2_cart(:)
    ! u' field
    uxy1(:,j) = u1_cart(:)
    uxy2(:,j) = u2_cart(:)
    ! u field
    ufull1(:,j) = u1_cart(:)+ubar1_tt(j)+u0 ! add basic state shear
    ufull2(:,j) = u2_cart(:)+ubar2_tt(j)
    ! p' field
    pxy1(:,j) = p1_cart(:)
    pxy2(:,j) = p2_cart(:)
    ! forcing field (update latest forcing values)
    fxy1(:,j,2) = f1_cart(:)

    ! Next get some derived ***advection*** quantities
    ! The meridional pv flux
    vq1_cart(:) = v1_cart(:)*q1_cart(:)
    vq2_cart(:) = v2_cart(:)*q2_cart(:) ! pv flux
    ! The x,y resolution advection terms
    adv1_cart(:) = qx1_cart(:)*(u1_cart(:) + ubar1_tt(j) + u0) &
         + v1_cart(:)*(qy1_cart(:) + beta + (u0/(rd*rd)) + qybar1_tt(j))
    adv2_cart(:) = qx2_cart(:)*(u2_cart(:) + ubar2_tt(j)) &
         + v2_cart(:)*(qy2_cart(:) + beta - (u0/(rd*rd)) + qybar2_tt(j))
    ! The total meridional flux here
    qflux1_tt(j) = sum(vq1_cart(:))/float(imax) ! meridional PV flux maybe
    qflux2_tt(j) = sum(vq2_cart(:))/float(imax)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Forward FT of zonal wind and forcing
    ! The srcft stands for real-to-complex Fourier transform
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    as = 1./float(imax)
    adv1_dft(:) = zero
    adv2_dft(:) = zero
    f1_dft(:) = zero
    call srcft(adv1_cart,adv1_dft,imax,as,hy)
    call srcft(adv2_cart,adv2_dft,imax,as,hy)
    call srcft(f1_cart,f1_dft,imax,as,hy)
    adv1_r(2:itrunc,j) = 2.0*real(adv1_dft(2:itrunc))
    adv1_i(2:itrunc,j) = 2.0*aimag(adv1_dft(2:itrunc))
    adv2_r(2:itrunc,j) = 2.0*real(adv2_dft(2:itrunc))
    adv2_i(2:itrunc,j) = 2.0*aimag(adv2_dft(2:itrunc))
    f1_r(2:itrunc,j) = 2.0*real(f1_dft(2:itrunc))
    f1_i(2:itrunc,j) = 2.0*aimag(f1_dft(2:itrunc))
  enddo ! loop over y in cartesian space

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Get zonal advection and forcing
  ! Will have to transform u*q flux
  ! What about the v*q flux? Does not appear here?
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i = 1,idft
    ! Preparation
    adv1_r_tt(:) = adv1_r(i,:)
    adv1_i_tt(:) = adv1_i(i,:)
    adv2_r_tt(:) = adv2_r(i,:)
    adv2_i_tt(:) = adv2_i(i,:)
    f1_r_tt(:) = f1_r(i,:)
    f1_i_tt(:) = f1_i(i,:)

    ! Forward sine transforms
    tt_type=0
    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(adv1_r_tt,handle,ipar,spar,ir)
    call d_forward_trig_transform(adv1_r_tt,handle,ipar,spar,ir)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(adv1_i_tt,handle,ipar,spar,ir)
    call d_forward_trig_transform(adv1_i_tt,handle,ipar,spar,ir)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(f1_r_tt,handle,ipar,spar,ir)
    call d_forward_trig_transform(f1_r_tt,handle,ipar,spar,ir)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(f1_i_tt,handle,ipar,spar,ir)
    call d_forward_trig_transform(f1_i_tt,handle,ipar,spar,ir)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(adv2_r_tt,handle,ipar,spar,ir)
    call d_forward_trig_transform(adv2_r_tt,handle,ipar,spar,ir)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(adv2_i_tt,handle,ipar,spar,ir)
    call d_forward_trig_transform(adv2_i_tt,handle,ipar,spar,ir)
    call free_trig_transform(handle,ipar,ir)

    ! Finally reconstruct nonlinear jacobian and forcing terms in spectral space
    ! Zonal advection of vorticity below (consider Dq/Dt equation, and we
    ! move advection terms to RHS so they are ***negative***)
    do j = 2,jtrunc
      adv1_c(i,j,3) = -(one_r*adv1_r_tt(j)+one_i*adv1_i_tt(j)) 
      adv2_c(i,j,3) = -(one_r*adv2_r_tt(j)+one_i*adv2_i_tt(j)) 
      force1_c(i,j) = one_r*f1_r_tt(j)+one_i*f1_i_tt(j) 
    enddo
  enddo ! loop over x values

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Get flux convergence
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Sine transforms
  tt_type=0
  call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
  call d_commit_trig_transform(qflux1_tt,handle,ipar,spar,ir)
  call d_forward_trig_transform(qflux1_tt,handle,ipar,spar,ir)
  call free_trig_transform(handle,ipar,ir)

  call d_init_trig_transform(jmax-1,tt_type,ipar,spar,ir)
  call d_commit_trig_transform(qflux2_tt,handle,ipar,spar,ir)
  call d_forward_trig_transform(qflux2_tt,handle,ipar,spar,ir)
  call free_trig_transform(handle,ipar,ir)

  do j = 2,jtrunc
    ell = el*float(j-1)
    qyyflux1(j,3) = -ell*ell*qflux1_tt(j) ! these are 2nd derivatives
    qyyflux2(j,3) = -ell*ell*qflux2_tt(j) ! convergence of the flux?
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write diagnostic output, useful for reference; consider testing v in future
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! write(*,*) 'max flux',maxval(qflux1_tt),minval(qflux1_tt)
  if(mod(int(t),int(dt)*1000).eq.1) then
    write(*,*) 'Printing zonal mean diagnostics.'
    do j = 1,jmax
      write(*,*) j,'ubar = ',ubar1_tt(j)+u0,' qbar = ',qbar1_tt(j)
    enddo
  endif
  umax = max(maxval(ufull1),maxval(ufull2))
  write(*,*) ' umax = ',umax,' cfl = ',umax/(dx/dt)
  return

end subroutine
end module
