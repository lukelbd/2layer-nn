!**** Controls I/O *******

module IO
contains

  subroutine dataio

    use GLOBAL_VARIABLES
    implicit none
    integer :: i,j,ret,shape(2),shape2(2),dsadata
    ! character (len=8) :: tstring ! time string with leading zeros
    ! write(tstring, '(i8.8)'), t  ! flexible, allows for a few hundred days
    ! character (len=6) :: tstring   ! or just use timesteps
    ! write(tstring, '(i6.6)'), t/dt

    shape(1) = imx
    shape(2) = jmax+1

    shape2(1) = imax
    shape2(2) = jmax+1 

    ret = dsadata("q1.df",2,shape,qxy1)
    ret = dsadata('q2.df',2,shape,qxy2)

    do i = 1,imx
      do j = 1,jmax+1
        qxy1(i,j) = qxy1(i,j) + qbar1(j) ! total PV
        qxy2(i,j) = qxy2(i,j) + qbar2(j) ! total PV
      enddo
    enddo
    ret = dsadata('v1.df',2,shape,vxy1)
    ret = dsadata('v2.df',2,shape,vxy2)

    ret = dsadata('u1.df',2,shape,uxy1)
    ret = dsadata('u2.df',2,shape,uxy2)

    ret = dsadata('u1_total.df',2,shape,u_1_r)
    ret = dsadata('u2_total.df',2,shape,u_2_r)

    ret = dsadata('psi1.df',2,shape2,psi_1)
    ret = dsadata('psi2.df',2,shape2,psi_2)

    !ret = dsadata('q1_zonalmean.df',1,(/jmax+1/),qbar1)
    !ret = dsadata('q2_zonalmean.df',1,(/jmax+1/),qbar2)

    !ret = dsadata('dqdy1.df',1,jmax+1,qbar1)
    !ret = dsadata('dqdy2.df',1,jmax+1,qbar2)

    ret = dsadata('eke.df',1, 1, energy2)

    return
  end subroutine

end module
