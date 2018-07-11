! modify *existing* variables from *subroutine* call
module variables
  implicit none
  integer x
  integer y
end module

module test
  implicit none
  integer :: q = x*y
end module

subroutine init()
  use variables
  
  integer :: z = 50

  x = x + 2
  y = y + 1
endsubroutine

program main
  use variables
  use test
  implicit none

  print *, x, y

  call init()

  print *, x, y
  print *, z
  
  print *, q
  
endprogram
