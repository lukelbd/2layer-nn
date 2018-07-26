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
  integer :: i, j, ret, shape1(2), shape2(2), shape3(1), dsadata
  shape1(1) = imax
  shape1(2) = jmax
  shape2(1) = idft
  shape2(2) = jmax 
  shape3(1) = jmax

  ! Write data
  ! Total PV
  ret = dsadata('q1.df',2,shape1,qfull1_cart)
  ret = dsadata('q2.df',2,shape1,qfull2_cart)

  ! Total relative vorticity
  ret = dsadata('vor1.df',2,shape1,vorfull1_cart)
  ret = dsadata('vor2.df',2,shape1,vorfull2_cart)

  ! Total u wind
  ret = dsadata('u1.df',2,shape1,ufull1_cart)
  ret = dsadata('u2.df',2,shape1,ufull2_cart)

  ! Total v wind (always anomalous, b/c mass conservation)
  ret = dsadata('v1.df',2,shape1,v1_cart)
  ret = dsadata('v2.df',2,shape1,v2_cart)

  ! Streamfunction anomaly (TODO: get full version)
  ! ret = dsadata('psi1.df',2,shape1,psi1_cart)
  ! ret = dsadata('psi2.df',2,shape1,psi2_cart)

  ! Forcing
  ret = dsadata('f1.df',2,shape1,force1_cart(:,:,2))

  ! Energy
  ret = dsadata('eke.df',1,1,energy)

  ! ret = dsadata('u1_zonalmean.df',1,shape3,ubar1)
  ! ret = dsadata('u2_zonalmean.df',1,shape3,ubar2)

  ! ret = dsadata('q1_zonalmean.df',1,shape3,qbar1)
  ! ret = dsadata('q2_zonalmean.df',1,shape3,qbar2)


  return
end subroutine
end module
