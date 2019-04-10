! Test whether successive subroutinecalls modify data
module procedures
  implicit none
  integer :: global = 100
  contains
  subroutine sub0(x)
    integer x(:)
    x(3) = 10
    print *,global
  end subroutine
  subroutine sub1(x, n)
    ! Driver
    integer :: n
    integer :: x(:)
    integer, allocatable :: z(:)
    integer :: y(n)
    ! integer :: i = size(x)
    x(2) = 5
    allocate(z(3))
    z(1:3) = x(1:3)
    ! i = size(x)
    print *, 'size: ', size(y)
    print *, 'dynamical size: ', size(z)
    call sub0(x)
  end subroutine
  subroutine sub2(x)
    ! Interface
    integer :: x(:)
    call sub1(x, 2)
  end subroutine
  subroutine sub3(x)
    real :: x(2)
    double precision :: y(2)
    y = x
    print *, 'Single:', x, 'Double:', y
  end subroutine
  subroutine sub4(x)
    real       :: x(2)
    complex(8) :: y(2)
    complex(8) :: one = (1., 0.)
    y(1) = one*x(1)
    ! y(:) = x(2)
    print *,'On the fly:', x, y
    ! x(:) = 5 ! causes segfault! array was created on the fly, has no place in memory!
  end subroutine
end module

program main
  ! Driver
  use procedures
  implicit none
  real :: single(2)
  character(len=*), parameter :: string = "hello, world!"
  integer, dimension(5) :: array = (/1, 1, 1, 1, 1/)
  integer, dimension(5) :: copy
  integer, dimension(2) :: rep = (/1, 2/)
  integer, dimension(2,2) :: res
  integer, parameter :: x = 1
  integer :: y(x) = (/x/)
  print *,string
  copy = array
  call sub2(array)
  call sub3(single)
  call sub4((/4.0, 5.0/)) ! on-the-fly generation
  ! call sub2((/2, 2, 2, 2/)) ! fails, invalid memory reference -- if we write to it, must be declared
  ! res = spread(reshape(res, (/1,2/)), 2, 1)
  ! print *, res
  ! print *, res(:,1)
  print *, 'Copy:', copy, 'Array:', array
end program
