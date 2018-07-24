program test
  real :: a
  integer :: i=2
  integer :: j=0
  integer :: tstep
  integer :: n, m
  integer :: array(5) = 2
  complex :: comp(2) = (/ (20,-1), (0,5) /)
  integer :: grid1(3,3) = 5
  integer :: grid2(3,3) = 5
  real :: empty(10)
  ! weird(:,1:1) = (/1, 2, 3/)
  grid1(:,1:1) = reshape((/1, 2, 3/),(/3, 1/))
  grid2(:,1) = (/1, 2, 3/)
  ! minval(array)
  print *, grid1
  print *, grid2
  ! print *, real(comp(:))
  ! print *, sum(array)
  ! array(2:3) = 0
  ! print *, array
  ! print *,empty
  ! do i=1,5
  !   print *,i
  ! enddo
  ! array(:,:) = 10
  ! n = 50
  ! m = 40
  ! print *, n, m
  ! print *, 3*2.5
  ! a = 5.
  ! print *, a*i
  ! ! a = 2
  ! print *, 1.e-5
  ! ! loop
  ! do j = 0,10+2,2
  !   print *, j
  ! enddo
  ! print *, 'array', minval( (/5,minval(array)/) )
end program
