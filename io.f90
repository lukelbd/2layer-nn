!**** Controls I/O *******

module IO
contains

        subroutine dataio

        use GLOBAL_VARIABLES
        implicit none
        
        integer :: i,j,ret,shape(2), dsadata  

        shape(1) = imx
        shape(2) = jmax+1
        
       do i = 1,imx
       do j = 1,jmax+1
          qxy1(i,j) = qxy1(i,j) + qbar1(j)
          qxy2(i,j) = qxy2(i,j) + qbar2(j)
       enddo
       enddo
      ret = dsadata('qd1.df',2,shape,qxy1)
      ret = dsadata('qd2.df',2,shape,qxy2)
!     ret = dsadata('v1.df',2,shape,vxy1*qxy1)
!     ret = dsadata('v2.df',2,shape,vxy2*qxy2)
!     ret = dsadata('u1.df',2,shape,uxy1*qxy1)
!     ret = dsadata('u2.df',2,shape,uxy2*qxy2)
      ret = dsadata('qdb1.df',1,jmax+1,qbar1)
      ret = dsadata('qdb2.df',1,jmax+1,qbar2)


     return
     end subroutine

end module
