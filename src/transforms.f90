!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fourier (x-direction) and trigonometric (y-direction; sines/cosines only) transforms
! Uses MKL Intel libraries
! * Each of these accepts ***real*** arguments (want to store the 2D data single
!   precision for memory concerns), but performs transforms on ***double*** data.
!   This is more or less what happened in Noboru's original version -- double-precision
!   vectors were generated from single-precision 2D data, then run through transform
!   algorithms. So we aren't creating more arrays than before, except
!   in the trig algorithms where we create that middle temporary variable.
! * The ftt_rcft and btt_crft are convenient wrappers to transform 2D spectral
!   data to fully cartesian coordinates, and vice versa. Does not create more
!   arrays then before -- no copies are made in those functions.
! * For more information on MKL libraries, see:
!   https://scc.ustc.edu.cn/zlsc/sugon/intel/mkl/mkl_manual
!   Note trig transform info is under "Partial Differential Equations support"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module transforms
use mkl_dfti            ! includes some boilerplate stuff, and Fourier transforms
use mkl_trig_transforms ! includes trig transforms
implicit none
type(dfti_descriptor), pointer :: hftt, hbtt ! cannot 'initiate' data structure for trig transforms
  ! like we do for Fourier transforms, so just declare handles here instead of passing
  ! them to the function after being initialized in in main.f90

! Variables for trig and Fourier transforms
! Careful; names have to be different if 2 modules loaded
integer    :: ipar(128), stat, tt_type ! trigonometric transformation type for mkl_trig_transforms
complex(8) :: one_real=(1.,0.), one_imag=(0.,1.), zero_complex=(0.,0.)
contains

! Forward trig transform
! Usage: ftt(<input array>, <output array>, <transform type>, <array size>)
! Transform type can be 0 (sine transform) or 1 (cosine transform)
subroutine ftt(x_in, x_out, tt_type, nj)
  integer, intent(in) :: tt_type
  integer, intent(in) :: nj
  ! real(8) :: dpar(3*n/2+2)
  real(8) :: dpar(5*nj/2+2) ! documentation says we need 5*
  real, dimension(nj), intent(in)  :: x_in
  real(8), dimension(nj) :: tmp
  real, dimension(nj), intent(out) :: x_out
  tmp(:) = x_in(:)   ! implicitly converts to double
  call d_init_trig_transform(nj-1, tt_type, ipar, dpar, stat) ! 'n' passed to function should be array length minus 1
  call d_commit_trig_transform(tmp, hftt, ipar, dpar, stat)
  call d_forward_trig_transform(tmp, hftt, ipar, dpar, stat)
  call free_trig_transform(hftt, ipar, stat)
  x_out(:) = tmp(:)
end subroutine

! Backward trig transform
! Same syntax as forward version
subroutine btt(x_in, x_out, tt_type, nj)
  integer, intent(in) :: tt_type
  integer, intent(in) :: nj
  ! real :: dpar(3*n/2+2)
  real(8) :: dpar(5*nj/2+2) ! documentation says we need 5*
  real, dimension(nj), intent(in)  :: x_in
  real(8), dimension(nj) :: tmp
  real, dimension(nj), intent(out) :: x_out
  tmp(:) = x_in(:)
  call d_init_trig_transform(nj-1, tt_type, ipar, dpar, stat)
  call d_commit_trig_transform(tmp, hbtt, ipar, dpar, stat)
  call d_backward_trig_transform(tmp, hbtt, ipar, dpar, stat)
  call free_trig_transform(hbtt, ipar, stat)
  x_out(:) = tmp(:)
end subroutine

! Forward transform, 'real to complex'
! Remember for real signal, the N/2-N coefficients are complex conjugates
! of the 1-N coefficients. We can combine them.
subroutine rcft(x_r, x_c, n, as, h)
  real, intent(in)    :: as
  integer, intent(in) :: n
  real, intent(in)     :: x_r(n)
  complex, intent(out) :: x_c(1+n/2)
  real(8) :: tmp1(n+2)
  real(8) :: tmp2(n)
  integer :: i, ir, ic, l(1)
  type(dfti_descriptor), pointer, intent(in) :: h
  l(1) = n
  ! Execute cartesian-->spectral transform
  tmp2(:) = x_r(:)
  stat = DftiComputeForward(h, tmp2, tmp1)
  ! Construct complex coefficients from real array
  ! call f2trf(n,tmp2,tmp2,wfftr) ! alternative FFTRB library, from Rogue Wave software
  x_c(1) = one_real*tmp1(1) ! back to single
  do i = 2,n/2
    ir = 2*i-1
    ic = 2*i
    x_c(i) = one_real*tmp1(ir)+one_imag*tmp1(ic) ! back to single
  enddo
  x_c(1+n/2) = one_real*tmp1(n+1) ! back to single
  x_c(:) = x_c(:)*as
end subroutine

! Inverse transform, 'complex to real'
subroutine crft(x_c, x_r, n, as, h)
  real, intent(in)    :: as
  integer, intent(in) :: n
  real, intent(out)   :: x_r(n)
  complex, intent(in) :: x_c(1+n/2)
  real(8) :: tmp1(n+2)
  real(8) :: tmp2(n)
  integer :: i, k, ir, ic, l(1)
  type(dfti_descriptor), pointer, intent(in) :: h
  l(1) = n
  ! Decompose complex coefficients into real array
  tmp1(1) = real(x_c(1))
  tmp1(2) = 0.
  do i = 2,n/2
    ir = 2*i-1
    ic = 2*i
    tmp1(ir) = real(x_c(i))
    tmp1(ic) = aimag(x_c(i))
  enddo
  tmp1(n+1) = real(x_c(1+n/2))
  tmp1(n+2) = 0.
  ! Execute spectral-->cartesian transform
  ! call f2trb(n,tmp1,tmp1,wfftr) ! alternative FFTRB library, from Rogue Wave software
  stat = DftiComputeBackward(h, tmp1, tmp2)
  x_r(:) = tmp2(:)*as ! back to single
end subroutine

! Forward Fourier transform and trig transform
! Usage: ftt_rcft(<input array>, <output array>, <transform type>,
!                 <dim 1 size>, <dim 2 size>, <dim 1 trunc number>, <dft handle>)
subroutine ftt_rcft(x_in, x_out, tt_type, ni, nj, trunc, h)
  integer, intent(in) :: ni, nj, trunc, tt_type
  real :: as
  integer :: i, j
  real, dimension(ni,nj), intent(in)     :: x_in
  complex, dimension(1+ni/2,nj), intent(out) :: x_out
  complex, dimension(1+ni/2)    :: dft_out
  complex, dimension(1+ni/2,nj) :: tt_in
  real, dimension(nj)       :: tt_r, tt_i
  type(dfti_descriptor), pointer, intent(in) :: h
  ! Forward Fourier transform
  as = 1./float(ni)
  tt_in   = zero_complex
  do j = 1,nj ! real-to-complex (y-direction is in cartesian units)
    dft_out = zero_complex
    call rcft(x_in(:,j), dft_out, ni, as, h)
    tt_in(:,j) = dft_out
  enddo
  ! Inverse trig transform
  x_out = zero_complex
  do i = 2,trunc ! trig transform the Fourier coefficients we care about (un-truncated ones)
    tt_r = 0.0
    tt_i = 0.0
    call ftt( real(tt_in(i,:)), tt_r, tt_type, nj)
    call ftt(aimag(tt_in(i,:)), tt_i, tt_type, nj)
    x_out(i,:) = 2.0*(one_real*tt_r + one_imag*tt_i)
  enddo
end subroutine

! Backward Fourier transform and trig transform
subroutine btt_crft(x_in, x_out, tt_type, ni, nj, trunc, h)
  integer, intent(in) :: ni, nj, trunc, tt_type
  real :: as
  integer :: i, j
  complex, dimension(1+ni/2,nj), intent(in) :: x_in
  real, dimension(ni,nj), intent(out)   :: x_out
  complex, dimension(1+ni/2,nj) :: dft_in
  real, dimension(ni)       :: dft_out
  real, dimension(nj)       :: tt_r, tt_i
  type(dfti_descriptor), pointer, intent(in) :: h
  ! Inverse trig transform
  dft_in = zero_complex
  do i = 2,trunc ! trig transform the Fourier coefficients we care about (un-truncated ones)
    tt_r = 0.0
    tt_i = 0.0
    call btt( real(x_in(i,:)), tt_r, tt_type, nj)
    call btt(aimag(x_in(i,:)), tt_i, tt_type, nj)
    dft_in(i,:) = 0.5*(one_real*tt_r + one_imag*tt_i)
  enddo
  ! Inverse Fourier transform
  as = 1.0
  x_out = 0.0
  do j = 1,nj
    dft_out = 0.0
    call crft(dft_in(:,j), dft_out, ni, as, h)
    x_out(:,j) = dft_out
  enddo
end subroutine

end module
