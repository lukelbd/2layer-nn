!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Data output
! Note dsadata is a ***function***, and you gotta declare functions
! just like you declare a variable. Subroutines are different, just have
! to be 'call'ed.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module io
contains
subroutine dataio
  use global_variables
  implicit none
  integer :: i,j,ret,shape1(2),shape2(2),shape3(1),dsadata
  ! character (len=8) :: tstring ! time string with leading zeros
  ! write(tstring, '(i8.8)'), t  ! flexible, allows for a few hundred days
  ! character (len=6) :: tstring   ! or just use timesteps
  ! write(tstring, '(i6.6)'), t/dt

  shape1(1) = imax
  shape1(2) = jmax

  shape2(1) = idft
  shape2(2) = jmax 

  shape3(1) = jmax

  ret = dsadata('q1.df',2,shape1,qxy1)
  ret = dsadata('q2.df',2,shape1,qxy2)

  ret = dsadata('f1.df',2,shape1,fxy1) ! forcing

  ! do i = 1,imax
  !   do j = 1,jmax
  !     qxy1(i,j) = qxy1(i,j) + qbar1(j) ! total PV
  !     qxy2(i,j) = qxy2(i,j) + qbar2(j) ! total PV
  !   enddo
  ! enddo

  ret = dsadata('v1.df',2,shape1,vxy1)
  ret = dsadata('v2.df',2,shape1,vxy2)

  ret = dsadata('u1.df',2,shape1,uxy1)
  ret = dsadata('u2.df',2,shape1,uxy2)

  ! ret = dsadata('u1_zonalmean.df',1,shape3,ubar1)
  ! ret = dsadata('u2_zonalmean.df',1,shape3,ubar2)

  ! ret = dsadata('psi1.df',2,shape2,psi1_c)
  ! ret = dsadata('psi2.df',2,shape2,psi2_c)

  ! ret = dsadata('q1_zonalmean.df',1,shape3,qbar1)
  ! ret = dsadata('q2_zonalmean.df',1,shape3,qbar2)

  ret = dsadata('eke.df',1, 1, energy2)

  return
end subroutine
end module
