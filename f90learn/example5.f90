! Module
module routines
  integer, parameter :: local = 50
  contains
  subroutine routine
    ! implicit none
    common /commonvars/ global
    print *, global
    print *, local
  end subroutine
end module

program main
  ! Load shit
  use globals
  use routines

  ! Call shit
  common /commonvars/ global
  call routine
end program

