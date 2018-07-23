!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Invert streamfunction from vorticity
! Note all temporary variables declared in this file are 1D
! Any 2D variable appearing in code below is global
! For more information on MKL utilities:
!  1. run 'module load mkl'
!  2. go to $MKLROOT/include, and check out the
!     mkl_dfti.f90 and mkl_trig_transforms.f90 files
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
  double precision :: spar(3*jmax/2+2)
  real :: rk2,rkn,x,y,phi,uav1,uav2,beta1
  real :: rkk,ell,fac,as
  integer :: i,j,wcos,wsin
  complex :: pb(idft,jmax+1), pc(idft,jmax+1)

  ! Temporary variables, for storing sine transform perhaps
  double precision :: v1r(jmax+1), v2r(jmax+1), v1i(jmax+1), v2i(jmax+1), &
    u1r(jmax+1),   u1i(jmax+1),  u2r(jmax+1),  u2i(jmax+1),  &
    q1r(jmax+1),   q2r(jmax+1),  q1i(jmax+1),  q2i(jmax+1),  &
    qx1r(jmax+1),  qx2r(jmax+1), qx1i(jmax+1), qx2i(jmax+1), &
    qy1r(jmax+1),  qy2r(jmax+1), qy1i(jmax+1), qy2i(jmax+1), &
    p1r(jmax+1),   p2r(jmax+1),  p1i(jmax+1),  p2i(jmax+1),  &
    f1r(jmax+1),   f1i(jmax+1), &
    vqz1(jmax+1),  vqz2(jmax+1), &
    ub(jmax+1),    uc(jmax+1)

  ! Truncated Fourier coefficients
  double complex :: u1c(idft),  u2c(idft), &
    v1c(idft),  v2c(idft), &
    q1c(idft),  q2c(idft), &
    qx1c(idft), qx2c(idft), &
    qy1c(idft), qy2c(idft), &
    p1c(idft),  p2c(idft), &
    f1c(idft)

  ! Inverse transforms of the above values
  double precision :: u1(imax),  u2(imax), &
    v1(imax),  v2(imax), &
    q1(imax),  q2(imax), &
    qx1(imax), qx2(imax), &
    qy1(imax), qy2(imax), &
    p1(imax),  p2(imax), &
    f1(imax)

  ! For PV injection
  integer,allocatable :: iseed(:)
  integer             :: isize,idate(8)
  real                :: anglex,angley,amp_rand

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Mean PV analysis
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Invert d/dy zonal mean PV equation to get zonal mean u
  ! Also integrate d/dy zonal mean PV to get zonal mean PV
  ! * The umean are needed for radiation forcing in prognostics.f90, which is why
  !   we create the seemingly identical variables 'umean1/ubar1' and 'umean2/ubar2'
  ! * The latter will undergo in-place sine transformation in this module
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ub(:) = 0.
  uc(:) = 0.
  qbar1(:) = 0.
  qbar2(:) = 0.
  ubar1(:) = 0.
  ubar2(:) = 0.
  umean1(:) = 0.
  umean2(:) = 0.
  do j = 2,jtrunc
    ell = el*float(j-1) ! angular wavenumber
    ub(j) = (qymean1(j,3)+qymean2(j,3))/(ell**2)
    uc(j) = (qymean1(j,3)-qymean2(j,3))/(ell**2 + (2./(rd*rd)))
    umean1(j) = (ub(j)+uc(j))*0.5
    umean2(j) = (ub(j)-uc(j))*0.5
    qbar1(j) = qymean1(j,3)/ell ! d/dx(f(t)) = wavenum*F(wavenum) for sine transform (note wavenum is angular here)
    qbar2(j) = qymean2(j,3)/ell
    ubar1(j) = umean1(j)
    ubar2(j) = umean2(j)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transform sponge damping to spectral space 
  ! * Copied from stochastic forcing section, needs work
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! if (t.eq.0) then
  !   do j = 2,jtrunc
  !     tau_farray(j) = tau_farray
  !   enddo
  !   do j = jtrunc+1,jmax+1
  !     tau_farray(j) = 0. !
  !   end do
  !   vqz1(j) = 0.
  !   vqz2(j) = 0.
  !   uav1 = 0.
  !   uav2 = 0.
  !   do i = 1,imax
  !     vqz1(j) = vqz1(j) + v1(i)/float(imax)
  !     vqz2(j) = vqz2(j) + v2(i)/float(imax)
  !     uav1 = uav1 + u1(i)/float(imax)
  !     uav2 = uav2 + u2(i)/float(imax)
  !   enddo
  !   tt_type=0
  !   call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
  !   call d_commit_trig_transform(vqz1,handle,ipar,spar,ir)
  !   call d_forward_trig_transform(vqz1,handle,ipar,spar,ir)
  !   call free_trig_transform(handle,ipar,ir)
  ! end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transform zonal mean to physical space 
  ! * Apply inverse trigonometric transforms
  ! * Seems to perform ****in place**** modification of array values
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Inverse sine transforms
  ! Suitable if edges are always small/zero
  qby_1(:) = qymean1(:,3) ! below functions need a 1D series
  qby_2(:) = qymean2(:,3)
  tt_type=0
  call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
  call d_commit_trig_transform(ubar1,handle,ipar,spar,ir)
  call d_backward_trig_transform(ubar1,handle,ipar,spar,ir)
  call free_trig_transform(handle,ipar,ir)

  call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
  call d_commit_trig_transform(ubar2,handle,ipar,spar,ir)
  call d_backward_trig_transform(ubar2,handle,ipar,spar,ir)
  call free_trig_transform(handle,ipar,ir)

  call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
  call d_commit_trig_transform(qby_1,handle,ipar,spar,ir)
  call d_backward_trig_transform(qby_1,handle,ipar,spar,ir)
  call free_trig_transform(handle,ipar,ir)

  call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
  call d_commit_trig_transform(qby_2,handle,ipar,spar,ir)
  call d_backward_trig_transform(qby_2,handle,ipar,spar,ir)
  call free_trig_transform(handle,ipar,ir)

  ! Inverse cosine transforms
  ! Suitble if edges are always non-zero
  tt_type=1
  call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
  call d_commit_trig_transform(qbar1,handle,ipar,spar,ir)
  call d_backward_trig_transform(qbar1,handle,ipar,spar,ir)
  call free_trig_transform(handle,ipar,ir)

  call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
  call d_commit_trig_transform(qbar2,handle,ipar,spar,ir)
  call d_backward_trig_transform(qbar2,handle,ipar,spar,ir)
  call free_trig_transform(handle,ipar,ir)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute full zonal-mean PV gradient in physical space
  ! * Adds in the background component by integrating the dq/dy
  !   equation from lower to upper boundary.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do j = 1,jmax+1
    y = dy*float(j-1)-0.5*width 
    qbar1(j) = qbar1(j)+(beta+u0/(rd*rd))*y
    qbar2(j) = qbar2(j)+(beta-u0/(rd*rd))*y
  enddo
  if(mod(t,10000).eq.1) then
    do j = 1,jmax+1
      write(*,*) j,'u1 :',ubar1(j)+u0,'   qbar1 :',qbar1(j)
    enddo
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 2D vort analysis
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Obtain streamfunction from full vorticity
  ! Also obtain u/v, and q derivatives in x/y (derivative in spectral coordinates)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transform eddy into physical space
  ! Used to just assign (1,:) and (:,1) zero but am fairly certain
  ! this was minor bug, and we were transforming arrays populated by 'empty'
  ! values above the trunc number (generally random, extremely small numbers)
  psi_1(:,:) = zero
  psi_2(:,:) = zero
  u_1(:,:)   = zero
  u_2(:,:)   = zero
  v_1(:,:)   = zero
  v_2(:,:)   = zero
  qx_1(:,:)  = zero
  qx_2(:,:)  = zero
  qy_1(:,:)  = zero
  qy_2(:,:)  = zero
  psi_1(:,:) = zero
  psi_2(:,:) = zero
  u_1(:,:)   = zero
  u_2(:,:)   = zero
  v_1(:,:)   = zero
  v_2(:,:)   = zero
  qx_1(:,:)  = zero
  qx_2(:,:)  = zero
  qy_1(:,:)  = zero
  qy_2(:,:)  = zero
  do j = 2,jtrunc
    ell = el*float(j-1)
    do i = 2,itrunc
      rkk = rk*float(i-1)
      pb(i,j) = -(vort_1(i,j,3)+vort_2(i,j,3))/(rkk**2 + ell**2)
      pc(i,j) = -(vort_1(i,j,3)-vort_2(i,j,3))/(rkk**2 + ell**2  &
        + (2./(rd*rd)))
      psi_1(i,j) = (pb(i,j)+pc(i,j))*0.5
      psi_2(i,j) = (pb(i,j)-pc(i,j))*0.5
      u_1(i,j) = -ell*psi_1(i,j)
      u_2(i,j) = -ell*psi_2(i,j)
      v_1(i,j) = one_i*rkk*psi_1(i,j)
      v_2(i,j) = one_i*rkk*psi_2(i,j)
      qx_1(i,j) = one_i*rkk*vort_1(i,j,3)
      qx_2(i,j) = one_i*rkk*vort_2(i,j,3)
      qy_1(i,j) = ell*vort_1(i,j,3)
      qy_2(i,j) = ell*vort_2(i,j,3)
    enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Perform inverse y-sine/cosine transformations
  ! of the x-direction Fourier coefficients for 2D params
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i = 1,idft
    ! First prepare 1D vectors containing the complex x-direction
    ! Fourier coefficients. Will supply these to inverse sine/cosine transform functions.
    u1r(:)  = real(u_1(i,:))
    u1i(:)  = aimag(u_1(i,:))
    u2r(:)  = real(u_2(i,:))
    u2i(:)  = aimag(u_2(i,:))
    v1r(:)  = real(v_1(i,:))
    v1i(:)  = aimag(v_1(i,:))
    v2r(:)  = real(v_2(i,:))
    v2i(:)  = aimag(v_2(i,:))
    q1r(:)  = real(vort_1(i,:,3))
    q1i(:)  = aimag(vort_1(i,:,3))
    qx1r(:) = real(qx_1(i,:))
    qx1i(:) = aimag(qx_1(i,:))
    qy1r(:) = real(qy_1(i,:))
    qy1i(:) = aimag(qy_1(i,:))
    q2r(:)  = real(vort_2(i,:,3))
    q2i(:)  = aimag(vort_2(i,:,3))
    qx2r(:) = real(qx_2(i,:))
    qx2i(:) = aimag(qx_2(i,:))
    qy2r(:) = real(qy_2(i,:))
    qy2i(:) = aimag(qy_2(i,:))
    p1r(:)  = real(psi_1(i,:))
    p1i(:)  = aimag(psi_1(i,:))
    p2r(:)  = real(psi_2(i,:))
    p2i(:)  = aimag(psi_2(i,:))

    ! Inverse sine transforms
    ! Suitable if edges are always small/zero
    ! Note we are transforming the x-Fourier coefficients here
    tt_type=0
    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(v1r,handle,ipar,spar,ir)
    call d_backward_trig_transform(v1r,handle,ipar,spar,ir)
    v_1_r(i,:) = v1r(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(v1i,handle,ipar,spar,ir)
    call d_backward_trig_transform(v1i,handle,ipar,spar,ir)
    v_1_i(i,:) = v1i(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(q1r,handle,ipar,spar,ir)
    call d_backward_trig_transform(q1r,handle,ipar,spar,ir)
    q_1_r(i,:) = q1r(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(q1i,handle,ipar,spar,ir)
    call d_backward_trig_transform(q1i,handle,ipar,spar,ir)
    q_1_i(i,:) = q1i(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(qx1r,handle,ipar,spar,ir)
    call d_backward_trig_transform(qx1r,handle,ipar,spar,ir)
    qx_1_r(i,:) = qx1r(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(qx1i,handle,ipar,spar,ir)
    call d_backward_trig_transform(qx1i,handle,ipar,spar,ir)
    qx_1_i(i,:) = qx1i(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(p1r,handle,ipar,spar,ir)
    call d_backward_trig_transform(p1r,handle,ipar,spar,ir)
    p_1_r(i,:) = p1r(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(p1i,handle,ipar,spar,ir)
    call d_backward_trig_transform(p1i,handle,ipar,spar,ir)
    p_1_i(i,:) = p1i(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(v2r,handle,ipar,spar,ir)
    call d_backward_trig_transform(v2r,handle,ipar,spar,ir)
    v_2_r(i,:) = v2r(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(v2i,handle,ipar,spar,ir)
    call d_backward_trig_transform(v2i,handle,ipar,spar,ir)
    v_2_i(i,:) = v2i(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(q2r,handle,ipar,spar,ir)
    call d_backward_trig_transform(q2r,handle,ipar,spar,ir)
    q_2_r(i,:) = q2r(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(q2i,handle,ipar,spar,ir)
    call d_backward_trig_transform(q2i,handle,ipar,spar,ir)
    q_2_i(i,:) = q2i(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(qx2r,handle,ipar,spar,ir)
    call d_backward_trig_transform(qx2r,handle,ipar,spar,ir)
    qx_2_r(i,:) = qx2r(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(qx2i,handle,ipar,spar,ir)
    call d_backward_trig_transform(qx2i,handle,ipar,spar,ir)
    qx_2_i(i,:) = qx2i(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(p2r,handle,ipar,spar,ir)
    call d_backward_trig_transform(p2r,handle,ipar,spar,ir)
    p_2_r(i,:) = p2r(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(p2i,handle,ipar,spar,ir)
    call d_backward_trig_transform(p2i,handle,ipar,spar,ir)
    p_2_i(i,:) = p2i(:)
    call free_trig_transform(handle,ipar,ir)

    ! Inverse cosine transforms
    ! Suitble if edges are always non-zero
    tt_type=1
    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(u1r,handle,ipar,spar,ir)
    call d_backward_trig_transform(u1r,handle,ipar,spar,ir)
    u_1_r(i,:) = u1r(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(u1i,handle,ipar,spar,ir)
    call d_backward_trig_transform(u1i,handle,ipar,spar,ir)
    u_1_i(i,:) = u1i(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(qy1r,handle,ipar,spar,ir)
    call d_backward_trig_transform(qy1r,handle,ipar,spar,ir)
    qy_1_r(i,:) = qy1r(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(qy1i,handle,ipar,spar,ir)
    call d_backward_trig_transform(qy1i,handle,ipar,spar,ir)
    qy_1_i(i,:) = qy1i(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(u2r,handle,ipar,spar,ir)
    call d_backward_trig_transform(u2r,handle,ipar,spar,ir)
    u_2_r(i,:) = u2r(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(u2i,handle,ipar,spar,ir)
    call d_backward_trig_transform(u2i,handle,ipar,spar,ir)
    u_2_i(i,:) = u2i(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(qy2r,handle,ipar,spar,ir)
    call d_backward_trig_transform(qy2r,handle,ipar,spar,ir)
    qy_2_r(i,:) = qy2r(:)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(qy2i,handle,ipar,spar,ir)
    call d_backward_trig_transform(qy2i,handle,ipar,spar,ir)
    qy_2_i(i,:) = qy2i(:)
    call free_trig_transform(handle,ipar,ir)
  enddo ! loop over idft

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Inverse x-Fourier transforms for 2D params
  ! * The y dimension now has *physical* units
  ! * Note the full Fourier transform equations take input/output arrays; they
  !   do *not* change arrays in place.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do j = 1,jmax+1
    ! Prepare 1D vectors that will be passed to Fourier inversion functions
    ! Also apply truncation here
    v1c(:)  = zero
    v2c(:)  = zero
    u1c(:)  = zero
    u2c(:)  = zero
    q1c(:)  = zero
    q2c(:)  = zero
    qx1c(:) = zero
    qx2c(:) = zero
    qy1c(:) = zero
    qy2c(:) = zero
    p1c(:)  = zero
    p2c(:)  = zero
    v1c(2:itrunc)  = 0.5 * (one_r*v_1_r(2:itrunc,j)  + one_i*v_1_i(2:itrunc,j))
    v2c(2:itrunc)  = 0.5 * (one_r*v_2_r(2:itrunc,j)  + one_i*v_2_i(2:itrunc,j))
    u1c(2:itrunc)  = 0.5 * (one_r*u_1_r(2:itrunc,j)  + one_i*u_1_i(2:itrunc,j))
    u2c(2:itrunc)  = 0.5 * (one_r*u_2_r(2:itrunc,j)  + one_i*u_2_i(2:itrunc,j))
    q1c(2:itrunc)  = 0.5 * (one_r*q_1_r(2:itrunc,j)  + one_i*q_1_i(2:itrunc,j))
    q2c(2:itrunc)  = 0.5 * (one_r*q_2_r(2:itrunc,j)  + one_i*q_2_i(2:itrunc,j))
    qx1c(2:itrunc) = 0.5 * (one_r*qx_1_r(2:itrunc,j) + one_i*qx_1_i(2:itrunc,j))
    qx2c(2:itrunc) = 0.5 * (one_r*qx_2_r(2:itrunc,j) + one_i*qx_2_i(2:itrunc,j))
    qy1c(2:itrunc) = 0.5 * (one_r*qy_1_r(2:itrunc,j) + one_i*qy_1_i(2:itrunc,j))
    qy2c(2:itrunc) = 0.5 * (one_r*qy_2_r(2:itrunc,j) + one_i*qy_2_i(2:itrunc,j))
    p1c(2:itrunc)  = 0.5 * (one_r*p_1_r(2:itrunc,j)  + one_i*p_1_i(2:itrunc,j))
    p2c(2:itrunc)  = 0.5 * (one_r*p_2_r(2:itrunc,j)  + one_i*p_2_i(2:itrunc,j))
    v1c(itrunc)  = one_r*real(v1c(itrunc))
    v2c(itrunc)  = one_r*real(v2c(itrunc))
    u1c(itrunc)  = one_r*real(u1c(itrunc))
    u2c(itrunc)  = one_r*real(u2c(itrunc))
    q1c(itrunc)  = one_r*real(q1c(itrunc))
    q2c(itrunc)  = one_r*real(q2c(itrunc))
    qx1c(itrunc) = one_r*real(qx1c(itrunc))
    qx2c(itrunc) = one_r*real(qx2c(itrunc))
    qy1c(itrunc) = one_r*real(qy1c(itrunc))
    qy2c(itrunc) = one_r*real(qy2c(itrunc))
    p1c(itrunc)  = one_r*real(p1c(itrunc))
    p2c(itrunc)  = one_r*real(p2c(itrunc))

    ! Inverse Fourier Transform
    ! Note we ***write to global variables here***
    ! The scrfr means complex-to-real fourier transform
    as = 1.
    v1(:)  = 0.
    v2(:)  = 0.
    u1(:)  = 0.
    u2(:)  = 0.
    q1(:)  = 0.
    q2(:)  = 0.
    qx1(:) = 0.
    qx2(:) = 0.
    qy1(:) = 0.
    qy2(:) = 0.
    p1(:)  = 0.
    p2(:)  = 0.
    call scrft(v1c,v1,imax,as,hx)
    call scrft(v2c,v2,imax,as,hx)
    call scrft(u1c,u1,imax,as,hx)
    call scrft(u2c,u2,imax,as,hx)
    call scrft(q1c,q1,imax,as,hx)
    call scrft(q2c,q2,imax,as,hx)
    call scrft(qx1c,qx1,imax,as,hx)
    call scrft(qx2c,qx2,imax,as,hx)
    call scrft(qy1c,qy1,imax,as,hx)
    call scrft(qy2c,qy2,imax,as,hx)
    call scrft(p1c,p1,imax,as,hx)
    call scrft(p2c,p2,imax,as,hx)
    ! q', dq'/dx, dq'/dy fields
    qxy1(:,j) = q1(:)
    qxy2(:,j) = q2(:)
    qxx1(:,j) = qx1(:)
    qxx2(:,j) = qx2(:)
    qyy1(:,j) = qy1(:)
    qyy2(:,j) = qy2(:)
    ! v' field
    vxy1(:,j) = v1(:)
    vxy2(:,j) = v2(:)
    ! u' field
    uxy1(:,j) = u1(:)
    uxy2(:,j) = u2(:)
    ! p' field
    pxy1(:,j) = p1(:)
    pxy2(:,j) = p2(:)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Data in real x-y coordinates found here vvvvv

    ! Stochastic forcing (small-scale narrow band)
    ! Injects cosines and sines onto the raw 2D data, with fixed memory 0.5
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
      fxy1(i,j,2) = 0.5*fxy1(i,j,1)
      do wcos = wmin_i,wmax_i
        rkk = rk*float(wcos-1)
        do wsin = wmin_i,wmax_i
          ell = el*float(wsin-1)
          call random_number(amp_rand)
          fxy1(i,j,2) = fxy1(i,j,2) &
          + 0.5*amp_i*amp_rand*ymask(j)*  &
            sin(ell*float(j-1)/float(jmax)+angley)*  &
            cos(rkk*float(i-1)/float(imax)+anglex)
        enddo
      enddo
    enddo

    ! Apply ***advection*** and ***Coriolis turning*** to
    ! the wind here. Super important step.
    do i = 1,imax
      v1(i) = vxy1(i,j)*qxy1(i,j)
      v2(i) = vxy2(i,j)*qxy2(i,j) ! pv flux
      u1(i) = (uxy1(i,j)+ubar1(j)+u0)*qxx1(i,j)
      u1(i) = u1(i) + vxy1(i,j)*(qyy1(i,j)+beta+(u0/(rd*rd))+qby_1(j))
      u2(i) = (uxy2(i,j)+ubar2(j))*qxx2(i,j)
      u2(i) = u2(i) + vxy2(i,j)*(qyy2(i,j)+beta-(u0/(rd*rd))+qby_2(j))
      uf(i,j) = uxy1(i,j)+ubar1(j)+u0
      f1(i) = fxy1(i,j,2)
    enddo

    ! Get zonal averages of real values
    vqz1(j) = sum(v1)/float(imax) ! meridional PV flux maybe
    vqz2(j) = sum(v2)/float(imax)
    uav1 = sum(u1)/float(imax)
    uav2 = sum(u2)/float(imax)

    ! Data in real x-y coordinates found here ^^^^
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Forward FT of zonal wind and forcing
    ! * The srcft stands for real-to-complex Fourier transform
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    as = 1./float(imax)
    u1c(:) = zero
    u2c(:) = zero
    f1c(:) = zero
    call srcft(u1,u1c,imax,as,hy)
    call srcft(u2,u2c,imax,as,hy)
    call srcft(f1,f1c,imax,as,hy)
    u_1_r(2:itrunc,j) = 2.0*real(u1c(2:itrunc))
    u_1_i(2:itrunc,j) = 2.0*aimag(u1c(2:itrunc))
    u_2_r(2:itrunc,j) = 2.0*real(u2c(2:itrunc))
    u_2_i(2:itrunc,j) = 2.0*aimag(u2c(2:itrunc))
    f_1_r(2:itrunc,j) = 2.0*real(f1c(2:itrunc))
    f_1_i(2:itrunc,j) = 2.0*aimag(f1c(2:itrunc))
  enddo ! loop over y values

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Forward sine/cosine transforms on y values
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i = 1,idft
    ! Preparation
    u1r(:) = u_1_r(i,:)
    u1i(:) = u_1_i(i,:)
    u2r(:) = u_2_r(i,:)
    u2i(:) = u_2_i(i,:)
    f1r(:) = f_1_r(i,:)
    f1i(:) = f_1_i(i,:)

    ! Inverse sine transforms
    tt_type=0
    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(u1r,handle,ipar,spar,ir)
    call d_forward_trig_transform(u1r,handle,ipar,spar,ir)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(u1i,handle,ipar,spar,ir)
    call d_forward_trig_transform(u1i,handle,ipar,spar,ir)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(f1r,handle,ipar,spar,ir)
    call d_forward_trig_transform(f1r,handle,ipar,spar,ir)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(f1i,handle,ipar,spar,ir)
    call d_forward_trig_transform(f1i,handle,ipar,spar,ir)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(u2r,handle,ipar,spar,ir)
    call d_forward_trig_transform(u2r,handle,ipar,spar,ir)
    call free_trig_transform(handle,ipar,ir)

    call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
    call d_commit_trig_transform(u2i,handle,ipar,spar,ir)
    call d_forward_trig_transform(u2i,handle,ipar,spar,ir)
    call free_trig_transform(handle,ipar,ir)

    ! Finally reconstruct nonlinear jacobian and forcing terms in spectral space
    ! Advection of vorticity is below
    do j = 1,jtrunc
      adv_1(i,j,3) = -(one_r*u1r(j)+one_i*u1i(j)) 
      adv_2(i,j,3) = -(one_r*u2r(j)+one_i*u2i(j)) 
      force_1(i,j) = one_r*f1r(j)+one_i*f1i(j) 
    enddo
    do j = jtrunc+1,itrunc
      adv_1(i,j,3) = zero
      adv_2(i,j,3) = zero
      force_1(i,j) = zero
    enddo
  enddo ! loop over x values

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transform zonal mean PV flux into spectral space 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Sine transforms
  tt_type=0
  call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
  call d_commit_trig_transform(vqz1,handle,ipar,spar,ir)
  call d_forward_trig_transform(vqz1,handle,ipar,spar,ir)
  call free_trig_transform(handle,ipar,ir)

  call d_init_trig_transform(jmax,tt_type,ipar,spar,ir)
  call d_commit_trig_transform(vqz2,handle,ipar,spar,ir)
  call d_forward_trig_transform(vqz2,handle,ipar,spar,ir)
  call free_trig_transform(handle,ipar,ir)

  do j = 2,jtrunc
    ell = el*float(j-1)
    vqm_1(j,3) = -ell*ell*vqz1(j)
    vqm_2(j,3) = -ell*ell*vqz2(j)
  enddo
  do j = jtrunc+1,jmax+1
    vqm_1(j,3) = 0.
    vqm_2(j,3) = 0.
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write diagnostic output, useful for reference; consider testing v in future
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! write(*,*) 'max flux',maxval(vqz1),minval(vqz1)
  if(mod(int(t),int(dt)*1000).eq.1) then
    write(*,*) 'Printing zonal mean diagnostics.'
    do j = 1,jmax+1
      write(*,*) j,'ubar = ',ubar1(j)+u0,' qbar = ',qbar1(j)
    enddo
  endif
  write(*,*) ' umax = ',maxval(uf),' cfl = ',maxval(uf)/(dx/dt)
  return

end subroutine
end module
