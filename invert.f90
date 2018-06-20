! ****** invert streamfunction from vorticity ****

module INVERT
contains

        subroutine invert1(hx,hy)

        use GLOBAL_VARIABLES
        use mkl_dfti
        use mkl_trig_transforms
        use FTM

        implicit none

        integer tt_type,ret,dsadata,shape(2)
        integer ir, ipar(128)
        double precision :: spar(3*jmax/2+2)
        real :: rk2,rkn,x,y,phi,uav1,uav2,beta1
        type(dfti_descriptor), pointer :: hx,hy,handle

        integer :: i,j,k,l
        real :: rkk,ell,fac,as
        double precision :: v1r(jmax+1),v2r(jmax+1),v1i(jmax+1),v2i(jmax+1)
        double precision :: f1r(jmax+1),f1i(jmax+1)
        double precision :: u1r(jmax+1),u1i(jmax+1),u2r(jmax+1),u2i(jmax+1)
        double precision :: q1r(jmax+1),q2r(jmax+1),q1i(jmax+1),q2i(jmax+1)
        double precision :: qx1r(jmax+1),qx2r(jmax+1),qx1i(jmax+1),qx2i(jmax+1)
        double precision :: qy1r(jmax+1),qy2r(jmax+1),qy1i(jmax+1),qy2i(jmax+1)
        double precision :: p1r(jmax+1),p2r(jmax+1),p1i(jmax+1),p2i(jmax+1)
        double precision :: ub(jmax+1),uc(jmax+1)
        double precision :: vqz_1(jmax+1),vqz_2(jmax+1)
        complex :: pb(imax,jmax+1),pc(imax,jmax+1)
        double precision :: u1(imx),u2(imx)
        double precision :: v1(imx),v2(imx)
        double precision :: q1(imx),q2(imx)
        double precision :: qx1(imx),qx2(imx)
        double precision :: qy1(imx),qy2(imx)
        double precision :: p1(imx),p2(imx)
        double precision :: f1(imx)
        double complex :: um1(imax),um2(imax)
        double complex :: vm1(imax),vm2(imax)
        double complex :: qm1(imax),qm2(imax)
        double complex :: qxm1(imax),qxm2(imax)
        double complex :: qym1(imax),qym2(imax)
        double complex :: pm1(imax),pm2(imax)
        double complex :: fm1(imax)

        integer::isize,idate(8)
        integer,allocatable::iseed(:)
        real::anglex,angley,ampl

!  Obtain zonal mean u from PV gradient

         qbar1(:) = 0.
         qbar2(:) = 0.
         ubar1(:) = 0.
         ubar2(:) = 0.
         ub(:) = 0.
         uc(:) = 0.

         do j = 2,nmax
           ell = el*float(j-1)
           ub(j) = (qymean1(j,3)+qymean2(j,3))/(ell**2)
           uc(j) = (qymean1(j,3)-qymean2(j,3))/(ell**2    &
                    + (2./(rd*rd)))
   
           umean1(j) = (ub(j)+uc(j))*0.5
           umean2(j) = (ub(j)-uc(j))*0.5
           qbar1(j) = qymean1(j,3)/ell
           qbar2(j) = qymean2(j,3)/ell
           ubar1(j) = umean1(j)
           ubar2(j) = umean2(j)
         enddo
         do j = nmax+1,jmax+1
           ubar1(j) = 0.
           ubar2(j) = 0.
           qbar1(j) = 0.
           qbar2(j) = 0.
         enddo

!  Transform zonal mean to physical space 

        qby_1(:) = qymean1(:,3)
        qby_2(:) = qymean2(:,3)

        tt_type=0
        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(ubar1,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(ubar1,handle,ipar,spar,ir)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(ubar2,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(ubar2,handle,ipar,spar,ir)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(qby_1,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(qby_1,handle,ipar,spar,ir)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(qby_2,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(qby_2,handle,ipar,spar,ir)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        tt_type=1
        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(qbar1,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(qbar1,handle,ipar,spar,ir)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(qbar2,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(qbar2,handle,ipar,spar,ir)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

! *** compute full zonal-mean PV gradient in physical space ****
      do j = 1,jmax+1
         y = dy*float(j-1)-0.5*width 
         qbar1(j) = qbar1(j)+(beta+u0/(rd*rd))*y 
         qbar2(j) = qbar2(j)+(beta-u0/(rd*rd))*y 
      enddo
    
    if(mod(m,10000).eq.1) then
      do j = 1,jmax+1
         write(6,*) j,'u1 :',ubar1(j)+u0,'   qbar1 :',qbar1(j)
      enddo
    endif

!  Obtain streamfunction from PV and 
!  transform eddy into physical space
         do j = 2,nmax
            ell = el*float(j-1)
         do i = 2,mmax
            rkk = rk*float(i-1)
          pb(i,j) = -(vort_1(i,j,3)+vort_2(i,j,3))/(rkk**2 + ell**2)
          pc(i,j) = -(vort_1(i,j,3)-vort_2(i,j,3))/(rkk**2 + ell**2  &
              + (2./(rd*rd)))
          psi_1(i,j) = (pb(i,j)+pc(i,j))*0.5
          psi_2(i,j) = (pb(i,j)-pc(i,j))*0.5
          u_1(i,j) = -ell*psi_1(i,j)
          u_2(i,j) = -ell*psi_2(i,j)
          v_1(i,j) = ui*rkk*psi_1(i,j)
          v_2(i,j) = ui*rkk*psi_2(i,j)
          qx_1(i,j) = ui*rkk*vort_1(i,j,3)
          qx_2(i,j) = ui*rkk*vort_2(i,j,3)
          qy_1(i,j) = ell*vort_1(i,j,3)
          qy_2(i,j) = ell*vort_2(i,j,3)
         enddo
         enddo

         psi_1(:,1) = zero
         psi_2(:,1) = zero
         u_1(:,1) = zero
         u_2(:,1) = zero
         v_1(:,1) = zero
         v_2(:,1) = zero
         qx_1(:,1) = zero
         qx_2(:,1) = zero
         qy_1(:,1) = zero
         qy_2(:,1) = zero
         psi_1(1,:) = zero
         psi_2(1,:) = zero
         u_1(1,:) = zero
         u_2(1,:) = zero
         v_1(1,:) = zero
         v_2(1,:) = zero
         qx_1(1,:) = zero
         qx_2(1,:) = zero
         qy_1(1,:) = zero
         qy_2(1,:) = zero

!  **** Transform u,v,vort to physical space ****

       do i = 1,imax
        do j = 1,jmax+1
          u1r(j) = real(u_1(i,j))
          u1i(j) = aimag(u_1(i,j))
          u2r(j) = real(u_2(i,j))
          u2i(j) = aimag(u_2(i,j))
          v1r(j) = real(v_1(i,j))
          v1i(j) = aimag(v_1(i,j))
          v2r(j) = real(v_2(i,j))
          v2i(j) = aimag(v_2(i,j))
          q1r(j) = real(vort_1(i,j,3))
          q1i(j) = aimag(vort_1(i,j,3))
          qx1r(j) = real(qx_1(i,j))
          qx1i(j) = aimag(qx_1(i,j))
          qy1r(j) = real(qy_1(i,j))
          qy1i(j) = aimag(qy_1(i,j))
          q2r(j) = real(vort_2(i,j,3))
          q2i(j) = aimag(vort_2(i,j,3))
          qx2r(j) = real(qx_2(i,j))
          qx2i(j) = aimag(qx_2(i,j))
          qy2r(j) = real(qy_2(i,j))
          qy2i(j) = aimag(qy_2(i,j))
          p1r(j) = real(psi_1(i,j))
          p1i(j) = aimag(psi_1(i,j))
          p2r(j) = real(psi_2(i,j))
          p2i(j) = aimag(psi_2(i,j))
        enddo

        tt_type=0
        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(v1r,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(v1r,handle,ipar,spar,ir)
        v_1_r(i,:) = v1r(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(v1i,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(v1i,handle,ipar,spar,ir)
        v_1_i(i,:) = v1i(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(q1r,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(q1r,handle,ipar,spar,ir)
        q_1_r(i,:) = q1r(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(q1i,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(q1i,handle,ipar,spar,ir)
        q_1_i(i,:) = q1i(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(qx1r,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(qx1r,handle,ipar,spar,ir)
        qx_1_r(i,:) = qx1r(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(qx1i,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(qx1i,handle,ipar,spar,ir)
        qx_1_i(i,:) = qx1i(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(p1r,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(p1r,handle,ipar,spar,ir)
        p_1_r(i,:) = p1r(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(p1i,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(p1i,handle,ipar,spar,ir)
        p_1_i(i,:) = p1i(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(v2r,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(v2r,handle,ipar,spar,ir)
        v_2_r(i,:) = v2r(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(v2i,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(v2i,handle,ipar,spar,ir)
        v_2_i(i,:) = v2i(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(q2r,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(q2r,handle,ipar,spar,ir)
        q_2_r(i,:) = q2r(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(q2i,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(q2i,handle,ipar,spar,ir)
        q_2_i(i,:) = q2i(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(qx2r,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(qx2r,handle,ipar,spar,ir)
        qx_2_r(i,:) = qx2r(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(qx2i,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(qx2i,handle,ipar,spar,ir)
        qx_2_i(i,:) = qx2i(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(p2r,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(p2r,handle,ipar,spar,ir)
        p_2_r(i,:) = p2r(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(p2i,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(p2i,handle,ipar,spar,ir)
        p_2_i(i,:) = p2i(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        tt_type=1
        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(u1r,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(u1r,handle,ipar,spar,ir)
        u_1_r(i,:) = u1r(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(u1i,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(u1i,handle,ipar,spar,ir)
        u_1_i(i,:) = u1i(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(qy1r,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(qy1r,handle,ipar,spar,ir)
        qy_1_r(i,:) = qy1r(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(qy1i,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(qy1i,handle,ipar,spar,ir)
        qy_1_i(i,:) = qy1i(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(u2r,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(u2r,handle,ipar,spar,ir)
        u_2_r(i,:) = u2r(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(u2i,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(u2i,handle,ipar,spar,ir)
        u_2_i(i,:) = u2i(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(qy2r,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(qy2r,handle,ipar,spar,ir)
        qy_2_r(i,:) = qy2r(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(qy2i,handle,ipar,spar,ir)
        CALL D_BACKWARD_TRIG_TRANSFORM(qy2i,handle,ipar,spar,ir)
        qy_2_i(i,:) = qy2i(:)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

     enddo

     do j = 1,jmax+1
            vm1 = zero
            vm2 = zero
            um1 = zero
            um2 = zero
            qm1 = zero
            qm2 = zero
            qxm1 = zero
            qxm2 = zero
            qym1 = zero
            qym2 = zero
            pm1 = zero
            pm2 = zero

        do i = 1,imax
            vm1(i) = ur*v_1_r(i,j)+ui*v_1_i(i,j)
            vm2(i) = ur*v_2_r(i,j)+ui*v_2_i(i,j)
            um1(i) = ur*u_1_r(i,j)+ui*u_1_i(i,j)
            um2(i) = ur*u_2_r(i,j)+ui*u_2_i(i,j)
            qm1(i) = ur*q_1_r(i,j)+ui*q_1_i(i,j)
            qm2(i) = ur*q_2_r(i,j)+ui*q_2_i(i,j)
            qxm1(i) = ur*qx_1_r(i,j)+ui*qx_1_i(i,j)
            qxm2(i) = ur*qx_2_r(i,j)+ui*qx_2_i(i,j)
            qym1(i) = ur*qy_1_r(i,j)+ui*qy_1_i(i,j)
            qym2(i) = ur*qy_2_r(i,j)+ui*qy_2_i(i,j)
            pm1(i) = ur*p_1_r(i,j)+ui*p_1_i(i,j)
            pm2(i) = ur*p_2_r(i,j)+ui*p_2_i(i,j)

            if(i.gt.1.and.i.lt.imax) then
               fac = 0.5
            else
               fac = 1.0
            endif

             vm1(i) = vm1(i)*fac
             vm2(i) = vm2(i)*fac
             um1(i) = um1(i)*fac
             um2(i) = um2(i)*fac
             qm1(i) = qm1(i)*fac
             qm2(i) = qm2(i)*fac
             qxm1(i) = qxm1(i)*fac
             qxm2(i) = qxm2(i)*fac
             qym1(i) = qym1(i)*fac
             qym2(i) = qym2(i)*fac
             pm1(i) = pm1(i)*fac
             pm2(i) = pm2(i)*fac

           if(i.eq.mmax) then
             vm1(i) = ur*real(vm1(i))
             vm2(i) = ur*real(vm2(i))
             um1(i) = ur*real(um1(i))
             um2(i) = ur*real(um2(i))
             qm1(i) = ur*real(qm1(i))
             qm2(i) = ur*real(qm2(i))
             qxm1(i) = ur*real(qxm1(i))
             qxm2(i) = ur*real(qxm2(i))
             qym1(i) = ur*real(qym1(i))
             qym2(i) = ur*real(qym2(i))
             pm1(i) = ur*real(pm1(i))
             pm2(i) = ur*real(pm2(i))
           endif
        enddo
             vm1(1) = zero
             vm2(1) = zero
             um1(1) = zero
             um2(1) = zero
             qm1(1) = zero
             qm2(1) = zero
             qxm1(1) = zero
             qxm2(1) = zero
             qym1(1) = zero
             qym2(1) = zero
             pm1(1) = zero
             pm2(1) = zero

        do i = mmax+1,imax
             vm1(i) = zero
             vm2(i) = zero
             um1(i) = zero
             um2(i) = zero
             qm1(i) = zero
             qm2(i) = zero
             qxm1(i) = zero
             qxm2(i) = zero
             qym1(i) = zero
             qym2(i) = zero
             pm1(i) = zero
             pm2(i) = zero
        enddo

      ! *** Inverse Fourier Transform ****

             v1 = 0.
             v2 = 0.
             u1 = 0.
             u2 = 0.
             q1 = 0.
             q2 = 0.
             qx1 = 0.
             qx2 = 0.
             qy1 = 0.
             qy2 = 0.
             p1 = 0.
             p2 = 0.

             as = 1.
             call scrft(vm1,v1,imx,as,hx)
             call scrft(vm2,v2,imx,as,hx)
             call scrft(um1,u1,imx,as,hx)
             call scrft(um2,u2,imx,as,hx)
             call scrft(qm1,q1,imx,as,hx)
             call scrft(qm2,q2,imx,as,hx)
             call scrft(qxm1,qx1,imx,as,hx)
             call scrft(qxm2,qx2,imx,as,hx)
             call scrft(qym1,qy1,imx,as,hx)
             call scrft(qym2,qy2,imx,as,hx)
             call scrft(pm1,p1,imx,as,hx)
             call scrft(pm2,p2,imx,as,hx)

!   q', dq'/dx, dq'/dy fields
         qxy1(:,j) = q1(:)
         qxy2(:,j) = q2(:)
         qxx1(:,j) = qx1(:)
         qxx2(:,j) = qx2(:)
         qyy1(:,j) = qy1(:)
         qyy2(:,j) = qy2(:)
!   v' field
         vxy1(:,j) = v1(:)
         vxy2(:,j) = v2(:)
!   u' field
         uxy1(:,j) = u1(:)
         uxy2(:,j) = u2(:)
!   p' field
         pxy1(:,j) = p1(:)
         pxy2(:,j) = p2(:)

!  ****  stochastic forcing (small-scale narrow band) ****

         if(m.eq.0) then
           call random_seed(size=isize)
           allocate(iseed(isize))
         endif
         call date_and_time(values=idate)
         call random_seed(size=isize)
         iseed=iseed*(idate(8)-500)
         if (m.eq.0) then
           iseed(1)=6453
           iseed(2)=887
         endif 
         call random_seed(put=iseed)
         call random_number(anglex)  
         call random_number(angley)  

     do i = 1,imx
!        fxy1(i,j,2) = exp(-dt/tau_fc)*fxy1(i,j,1)
         fxy1(i,j,2) = 0.5*fxy1(i,j,1)
     do k = 41,46
            rkk = rk*float(k-1)
     do l = 41,46
            ell = el*float(l-1)
         call random_number(ampl)  
         fxy1(i,j,2) = fxy1(i,j,2)             &
!       + (1.-exp(-dt/tau_fc))*famp*ymask(j)*  &
        + 0.5*famp*ampl*ymask(j)*  &
          sin(float(j-1)*ell/float(jmax)+angley)*  &
          cos(rkk*(float(i-1)/float(imx)+anglex)) 
     enddo
     enddo
     enddo


     do i = 1,imx
         v1(i) = vxy1(i,j)*qxy1(i,j)
         v2(i) = vxy2(i,j)*qxy2(i,j)
         u1(i) = (uxy1(i,j)+ubar1(j)+u0)*qxx1(i,j)
         u1(i) = u1(i) + vxy1(i,j)*(qyy1(i,j)+beta+(u0/(rd*rd))+qby_1(j))
         u2(i) = (uxy2(i,j)+ubar2(j))*qxx2(i,j)
         u2(i) = u2(i) + vxy2(i,j)*(qyy2(i,j)+beta-(u0/(rd*rd))+qby_2(j))
         uf(i,j) = uxy1(i,j)+ubar1(j)+u0
         f1(i) = fxy1(i,j,2)
     enddo

          vqz_1(j) = 0.
          vqz_2(j) = 0.
          uav1 = 0. 
          uav2 = 0. 
        do i = 1,imx
          vqz_1(j) = vqz_1(j) + v1(i)/float(imx)
          vqz_2(j) = vqz_2(j) + v2(i)/float(imx)
          uav1 = uav1 + u1(i)/float(imx)
          uav2 = uav2 + u2(i)/float(imx)
        enddo

             um1 = zero
             um2 = zero
             fm1 = zero

   !   Forward FT   
            as = 1./float(imx)
            call srcft(u1,um1,imx,as,hy)
            call srcft(u2,um2,imx,as,hy)
            call srcft(f1,fm1,imx,as,hy)

         do i = 1,imax
            if(i.gt.1.and.i.lt.imax) then
               fac = 2.0
            else
               fac = 1.0
            endif
            if(i.gt.mmax) fac = 0.
            if(i.eq.1) fac = 0.
  
            u_1_r(i,j) = fac*real(um1(i))
            u_1_i(i,j) = fac*aimag(um1(i))
            u_2_r(i,j) = fac*real(um2(i))
            u_2_i(i,j) = fac*aimag(um2(i))
            f_1_r(i,j) = fac*real(fm1(i))
            f_1_i(i,j) = fac*aimag(fm1(i))

          enddo

       enddo

       do i = 1,imax
        rkk = rk*float(i-1)
        u1r(:) = u_1_r(i,:)
        u1i(:) = u_1_i(i,:)
        u2r(:) = u_2_r(i,:)
        u2i(:) = u_2_i(i,:)
        f1r(:) = f_1_r(i,:)
        f1i(:) = f_1_i(i,:)

 !  Transform back to spectral space

        tt_type=0
        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(u1r,handle,ipar,spar,ir)
        CALL D_FORWARD_TRIG_TRANSFORM(u1r,handle,ipar,spar,ir)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(u1i,handle,ipar,spar,ir)
        CALL D_FORWARD_TRIG_TRANSFORM(u1i,handle,ipar,spar,ir)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(f1r,handle,ipar,spar,ir)
        CALL D_FORWARD_TRIG_TRANSFORM(f1r,handle,ipar,spar,ir)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(f1i,handle,ipar,spar,ir)
        CALL D_FORWARD_TRIG_TRANSFORM(f1i,handle,ipar,spar,ir)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(u2r,handle,ipar,spar,ir)
        CALL D_FORWARD_TRIG_TRANSFORM(u2r,handle,ipar,spar,ir)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(u2i,handle,ipar,spar,ir)
        CALL D_FORWARD_TRIG_TRANSFORM(u2i,handle,ipar,spar,ir)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

!  finally reconstruct nonlinear jacoian and forcing terms in spectral space

    do j = 1,nmax
        adv_1(i,j,3) = -(ur*u1r(j)+ui*u1i(j)) 
        adv_2(i,j,3) = -(ur*u2r(j)+ui*u2i(j)) 
        force_1(i,j) = ur*f1r(j)+ui*f1i(j) 
    enddo

    do j = nmax+1,mmax
        adv_1(i,j,3) = (0.,0.)
        adv_2(i,j,3) = (0.,0.)
        force_1(i,j) = (0.,0.)
    enddo

    enddo

!  Transform zonal mean PV flux into spectral space 
      
 !   write(6,*) 'max flux',maxval(vqz_1),minval(vqz_1)

        tt_type=0
        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(vqz_1,handle,ipar,spar,ir)
        CALL D_FORWARD_TRIG_TRANSFORM(vqz_1,handle,ipar,spar,ir)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(vqz_2,handle,ipar,spar,ir)
        CALL D_FORWARD_TRIG_TRANSFORM(vqz_2,handle,ipar,spar,ir)
        CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)

       do j = 2,nmax
        ell = el*float(j-1)
        vqm_1(j,3) = -ell*ell*vqz_1(j)
        vqm_2(j,3) = -ell*ell*vqz_2(j)
       enddo
       do j = nmax+1,jmax+1
        vqm_1(j,3) = 0.
        vqm_2(j,3) = 0.
       enddo
       if(mod(m,1000).eq.1) then
       !if((mod(m,md).eq.0).and.(m.ge.mds)) then
        write(6,*) ' max u = ',maxval(uf)
       endif
          return 
         end subroutine
end module
