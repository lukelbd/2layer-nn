!**** Controls I/O *******

module IO
contains

        subroutine dataio

        use GLOBAL_VARIABLES
        implicit none
        
        integer :: i,j,ret,shape(2), shape2(2), dsadata  
        character :: expname  

        shape(1) = imx
        shape(2) = jmax+1
       
        shape2(1) = imax
        shape2(2) = jmax+1 

       do i = 1,imx
       do j = 1,jmax+1
          qxy1(i,j) = qxy1(i,j) + qbar1(j) ! total PV
          qxy2(i,j) = qxy2(i,j) + qbar2(j) ! total PV
       enddo
       enddo
        PRINT*,exp_name
        PRINT*,"/scratch/midway2/t-970c07/"//exp_name
        !print('/scratch/midway2/t-970c07/'//expname//'/')
        ret = dsadata("/scratch/midway2/t-970c07/"//exp_name//"/q1_total.df",2,shape,qxy1)
        ret = dsadata('/scratch/midway2/t-970c07/'//exp_name//'/q2_total.df',2,shape,qxy2)
        
        ret =dsadata('/scratch/midway2/t-970c07/'//exp_name//'/v1.df',2,shape,vxy1*qxy1)
        ret =dsadata('/scratch/midway2/t-970c07/'//exp_name//'/v2.df',2,shape,vxy2*qxy2)
        
        ret =dsadata('/scratch/midway2/t-970c07/'//exp_name//'/u1.df',2,shape,uxy1*qxy1)
        ret =dsadata('/scratch/midway2/t-970c07/'//exp_name//'/u2.df',2,shape,uxy2*qxy2)

        ret=dsadata('/scratch/midway2/t-970c07/'//exp_name//'/u1_total.df',2,shape2,u_1_r)
        ret=dsadata('/scratch/midway2/t-970c07/'//exp_name//'/u2_total.df',2,shape2,u_2_r)

        ret=dsadata('/scratch/midway2/t-970c07/'//exp_name//'/psi1.df',2,shape2,psi_1)
        ret=dsadata('/scratch/midway2/t-970c07/'//exp_name//'/psi2.df',2,shape2,psi_2)
 
        ret =dsadata('/scratch/midway2/t-970c07/'//exp_name//'/q1_zonalmean.df',1,jmax+1,qbar1)
        ret =dsadata('/scratch/midway2/t-970c07/'//exp_name//'/q2_zonalmean.df',1,jmax+1,qbar2)

        !ret =dsadata('/scratch/midway2/t-970c07/'//exp_name//'/dqdy1.df',1,jmax+1,qbar1)
        !ret =dsadata('/scratch/midway2/t-970c07/'//exp_name//'/dqdy2.df',1,jmax+1,qbar2)      

        ret =dsadata('/scratch/midway2/t-970c07/'//exp_name//'/eke.df',1, 1, energy2)

     return
     end subroutine

end module
