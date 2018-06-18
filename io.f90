!**** Controls I/O *******


module IO
contains

     subroutine dataio

        use GLOBAL_VARIABLES
        implicit none
        
        integer :: i,j,ret,shape(2), dsadata  

        ! define output variables and places
        !opath = '/scratch/midway2/t-970c07/exp01'
        

        shape(1) = imx
        shape(2) = jmax+1
        
        do i = 1,imx
        do j = 1,jmax+1
          qxy1(i,j) = qxy1(i,j) + qbar1(j)
          qxy2(i,j) = qxy2(i,j) + qbar2(j)
        enddo
        enddo
        ret = dsadata('/scratch/midway2/t-970c07/exp01'//'/qd1.df',2,shape,qxy1)
        ret = dsadata('/scratch/midway2/t-970c07/exp01'//'/qd2.df',2,shape,qxy2)
        ret = dsadata('/scratch/midway2/t-970c07/exp01'//'/v1.df',2,shape,vxy1*qxy1)
        ret = dsadata('/scratch/midway2/t-970c07/exp01'//'/v2.df',2,shape,vxy2*qxy2)
        ret = dsadata('/scratch/midway2/t-970c07/exp01'//'/u1.df',2,shape,uxy1*qxy1)
        ret = dsadata('/scratch/midway2/t-970c07/exp01'//'/u2.df',2,shape,uxy2*qxy2)
        ret = dsadata('/scratch/midway2/t-970c07/exp01'//'/qdb1.df',1,jmax+1,qbar1)
        ret = dsadata('/scratch/midway2/t-970c07/exp01'//'/qdb2.df',1,jmax+1,qbar2)
        
        ret = dsadata('/scratch/midway2/t-970c07/exp01'//'/dqdy1.df',1,jmax+1,qbar2)
        ret = dsadata('/scratch/midway2/t-970c07/exp01'//'/dqdy2.df',1,jmax+1,qbar2)

     return
     end subroutine

end module
