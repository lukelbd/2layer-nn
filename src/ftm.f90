module ftm
contains
subroutine srcft(x,y,n,as,hx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Forward transform, 'real to complex'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use mkl_dfti
  implicit none
  type(dfti_descriptor), pointer :: hx
  integer :: status
  integer :: n, i, iim, iip, l(1)
  double precision :: x(n), x1(n+2)
  double complex :: y(1+n/2), one_r, one_i
  real :: as
  l(1) = n
  one_r = (1.,0)
  one_i = (0.,1.)

  !   status = dfticreatedescriptor( hx, dfti_single, &
  !      dfti_real, 1, l)
  !   status = dfticommitdescriptor( hx )
  status = dfticomputeforward( hx, x, x1)
  !   status = dftifreedescriptor( hx )
  !        call f2trf(n,x,x,wfftr)

  y(1) = one_r*x1(1) 
  y(1) = y(1)*as

  do 67 i = 2,n/2
  iim = 2*i-1
  iip = 2*i
  y(i) = one_r*x1(iim)+one_i*x1(iip)
  y(i) = y(i)*as
  67     continue
  y(1+n/2) = one_r*x1(n+1)*as

  return
end subroutine srcft
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inverse transform, 'complex to real'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine scrft(y,x1,n,as,hx)
  use mkl_dfti
  implicit none
  type(dfti_descriptor), pointer :: hx
  integer :: status
  integer :: n,i,k,iim,iip
  double precision :: x(n+2),x1(n)
  double complex :: y(1+n/2)
  real :: as

  x(1) = real(y(1))
  x(2) = 0.
  do 67 i = 2,n/2
  iim = 2*i-1
  iip = 2*i
  x(iim) = real(y(i))
  x(iip) = aimag(y(i))
  67     continue

  x(n+1) = real(y(1+n/2))
  x(n+2) = 0.

  !   status = dfticreatedescriptor( hx, dfti_single, &
  !      dfti_real, 1, l)
  !   status = dfticommitdescriptor( hx )
  status = dfticomputebackward(hx, x, x1)
  !   status = dftifreedescriptor( hx )
  !        call f2trb(n,x,x,wfftr)
  x1 = x1*as
  !        do 168 k = 1,n
  !        x1(k) = x(k)*as
  ! 168     continue
  return
end subroutine scrft
end module
