!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modules for Fourier and trigonometric (sines/cosines only) transforms
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
! * The j-direction truncation number is passed, but so far does nothing. May
!   be used in future.
! * For more information on MKL libraries, see:
!   https://scc.ustc.edu.cn/zlsc/sugon/intel/mkl/mkl_manual
!   Note trig transform info is under "Partial Differential Equations support"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module spectral
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
! Transform type can be 0 (cosine transform) or 1 (sine transform)
subroutine ftt(x_in, x_out, tt_type, nj, jtrunc)
  integer :: tt_type
  integer :: nj, jtrunc
  ! real(8) :: dpar(3*n/2+2)
  real(8) :: dpar(5*nj/2+2) ! documentation says we need 5*
  real, dimension(nj)    :: x_in
  real(8), dimension(nj) :: tmp
  real, dimension(nj)    :: x_out
  tmp(:) = x_in(:)   ! implicitly converts to double
  call d_init_trig_transform(nj-1, tt_type, ipar, dpar, stat) ! 'n' passed to function should be array length minus 1
  call d_commit_trig_transform(tmp, hftt, ipar, dpar, stat)
  call d_forward_trig_transform(tmp, hftt, ipar, dpar, stat)
  call free_trig_transform(hftt, ipar, stat)
  x_out(:) = tmp(:)
end subroutine

! Backward trig transform
! Same syntax as forward version
subroutine btt(x_in, x_out, tt_type, nj, jtrunc)
  integer :: tt_type
  integer :: nj, jtrunc
  ! real :: dpar(3*n/2+2)
  real(8) :: dpar(5*nj/2+2) ! documentation says we need 5*
  real, dimension(nj)    :: x_in
  real(8), dimension(nj) :: tmp
  real, dimension(nj)    :: x_out
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
subroutine rcft(x, y, n, trunc, as, hy)
  real    :: as
  integer :: n, trunc, i, iim, iip, l(1)
  real    :: x(n)
  complex :: y(1+n/2)
  real(8) :: tmp1(n+2)
  real(8) :: tmp2(n)
  type(dfti_descriptor), pointer :: hy
  l(1) = n
  ! Execute cartesian-->spectral transform
  tmp2(:) = x(:)
  stat = dfticomputeforward(hy, tmp2, tmp1)
  ! Construct complex coefficients from real array
  ! call f2trf(n,tmp2,tmp2,wfftr) ! alternative FFTRB library, from Rogue Wave software
  y(1) = one_real*tmp1(1) ! back to single
  do 67 i = 2,n/2
  iim = 2*i-1
  iip = 2*i
  y(i) = one_real*tmp1(iim)+one_imag*tmp1(iip) ! back to single
  67 continue
  y(1+n/2) = one_real*tmp1(n+1) ! back to single
  y(:) = y(:)*as
end subroutine

! Inverse transform, 'complex to real'
subroutine crft(y, x, n, trunc, as, hx)
  real :: as
  integer :: n, trunc, i, k, iim, iip, l(1)
  real    :: x(n)
  complex :: y(1+n/2)
  real(8) :: tmp1(n+2)
  real(8) :: tmp2(n)
  type(dfti_descriptor), pointer :: hx
  l(1) = n
  ! Decompose complex coefficients into real array
  tmp1(1) = real(y(1))
  tmp1(2) = 0.
  do 67 i = 2,n/2
  iim = 2*i-1
  iip = 2*i
  tmp1(iim) = real(y(i))
  tmp1(iip) = aimag(y(i))
  67 continue
  tmp1(n+1) = real(y(1+n/2))
  tmp1(n+2) = 0.
  ! Execute spectral-->cartesian transform
  ! call f2trb(n,tmp1,tmp1,wfftr) ! alternative FFTRB library, from Rogue Wave software
  stat = dfticomputebackward(hx, tmp1, tmp2)
  x(:) = tmp2(:)*as ! back to single
end subroutine

! Forward Fourier transform and trig transform
! Usage: ftt(<input array>, <output array>, <transform type>, <dim 1 size>, <dim 2 size>, <trunc number>)
subroutine ftt_rcft(x_in, x_out, tt_type, ni, nj, itrunc, jtrunc, hy)
  real    :: as
  integer :: tt_type
  integer :: i, j, ni, nj, itrunc, jtrunc
  real, dimension(ni,nj)    :: x_in
  complex, dimension(ni,nj) :: x_out
  complex, dimension(ni)    :: dft_out
  complex, dimension(ni,nj) :: tt_in
  real, dimension(nj)       :: tt_r, tt_i
  type(dfti_descriptor), pointer :: hy
  ! Forward Fourier transform
  as = 1./float(ni)
  tt_in   = zero_complex
  dft_out = zero_complex
  do j = 1,nj
    call rcft(x_in(:,j), dft_out, ni, itrunc, as, hy)
    tt_in(:,j) = dft_out
  enddo
  ! Inverse trig transform
  x_out = zero_complex
  tt_r = 0.0
  tt_i = 0.0
  do i = 1,itrunc ! handle truncation here? why not
    call ftt( real(tt_in(i,:)), tt_r, tt_type, nj, jtrunc)
    call ftt(aimag(tt_in(i,:)), tt_i, tt_type, nj, jtrunc)
    if (i.eq.itrunc .or. i.eq.1) then
      x_out(i,:) = one_real*tt_r
    else
      x_out(i,:) = 2.0*(one_real*tt_r + one_imag*tt_i)
    end if
  enddo
end subroutine

! Backward Fourier transform and trig transform
! subroutine btt_crft(tt_type, x_in, xout_r, xout_i)
subroutine btt_crft(x_in, x_out, tt_type, ni, nj, itrunc, jtrunc, hx)
  real    :: as
  integer :: tt_type
  integer :: i, j, ni, nj, itrunc, jtrunc
  complex, dimension(ni,nj) :: x_in
  real, dimension(ni,nj)    :: x_out
  complex, dimension(ni,nj) :: dft_in
  real, dimension(ni)       :: dft_out
  real, dimension(nj)       :: tt_r, tt_i
  type(dfti_descriptor), pointer :: hx
  ! Inverse trig transform
  dft_in = zero_complex
  tt_r = 0.0
  tt_i = 0.0
  do i = 1,itrunc
    call btt( real(x_in(i,:)), tt_r, tt_type, nj, jtrunc)
    call btt(aimag(x_in(i,:)), tt_i, tt_type, nj, jtrunc)
    if (i.eq.itrunc .or. i.eq.1) then
      dft_in(i,:) = one_real*tt_r
    else
      dft_in(i,:) = 0.5*(one_real*tt_r + one_imag*tt_i)
    end if
  enddo
  ! Inverse Fourier transform
  as = 1.0
  x_out   = 0.0
  dft_out = 0.0
  do j = 1,nj
    call crft(dft_in(:,j), dft_out, ni, itrunc, as, hx)
    x_out(:,j) = dft_out
  enddo
end subroutine

end module
