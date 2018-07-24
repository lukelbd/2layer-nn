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
  shape1(1) = imax
  shape1(2) = jmax
  shape2(1) = idft
  shape2(2) = jmax 
  shape3(1) = jmax

  ! pv anomaly
  ret = dsadata('q1.df',2,shape1,q1_out)
  ret = dsadata('q2.df',2,shape1,q2_out)

  ! relative vorticity anomaly
  ret = dsadata('vor1.df',2,shape1,vor1_out)
  ret = dsadata('vor2.df',2,shape1,vor2_out)

  ! streamfunction
  ret = dsadata('psi1.df',2,shape1,psi1_out)
  ret = dsadata('psi2.df',2,shape1,psi2_out)

  ! forcing
  ret = dsadata('f1.df',2,shape1,f1_out)

  ! u wind
  ret = dsadata('u1.df',2,shape1,u1_out)
  ret = dsadata('u2.df',2,shape1,u2_out)

  ! v wind
  ret = dsadata('v1.df',2,shape1,v1_out)
  ret = dsadata('v2.df',2,shape1,v2_out)

  ! energy
  ret = dsadata('eke.df',1,1,energy)

  ! ret = dsadata('u1_zonalmean.df',1,shape3,ubar1)
  ! ret = dsadata('u2_zonalmean.df',1,shape3,ubar2)

  ! ret = dsadata('psi1.df',2,shape2,psi1_c)
  ! ret = dsadata('psi2.df',2,shape2,psi2_c)

  ! ret = dsadata('q1_zonalmean.df',1,shape3,qbar1)
  ! ret = dsadata('q2_zonalmean.df',1,shape3,qbar2)


  return
end subroutine
end module
