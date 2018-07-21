! communicate files declared in subroutine to 
program myprog 
  integer i1, i2
  real    f1, f2 

  common /john/ i1, f1    ! must be in same order !!!          


  i1 = 1111
  f1 = 2222
  i2 = 3333
  f2 = 4444

  call mysubr

  print *, i1, f1, i2, f2
end program

subroutine mysubr
  integer i1, i2
  real    f1, f2 

  common /john/ i1, f1    ! must be in same order !!!


  i1 = 5555
  f1 = 6666
  i2 = 7777
  f2 = 8888

  return
end subroutine
