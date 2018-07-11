! Scope experiment
! Can subroutines access parameters?
module globals
  integer, parameter :: foo = 10
end module

module routines
  ! use globals
  integer, parameter :: baz = 20
  ! common /attempt/ bar
  contains
  subroutine routine
    ! implicit none
    common /attempt/ bar
    integer :: bar = 50
    integer :: BAR = 100
    ! common /common_block/ foo
    ! print *, foo
    print *, bar
    print *, BAR
    print *, baz
    print *, foo
  end subroutine
end module

program main
  ! Load shit
  use globals
  use routines
  integer :: a = 0
  real :: b = 0.0
  common /attempt/ bar

  ! Call shit
  ! common /common_block/ foo
  call routine
  print *, "Main program."
  print *, bar
  print *, a .eq. b
end program

