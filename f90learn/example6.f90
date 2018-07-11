program main
  integer :: i = 1
  integer :: j = 20
  character (len=*), parameter :: foo='fart'
  character (len=20) :: bar
  character (len=12) :: baz
  write(bar, '(a,i3.3,i3.3)'), foo, i, j
  write(baz, '(i8.8)'), i
  print *, foo
  print *, bar
  print *, baz
end program
