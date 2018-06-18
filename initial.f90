! ***** initial **********
!       Initializes the flow field

module INITIAL
contains

        subroutine initial1

        use GLOBAL_VARIABLES
        use mkl_dfti
        use mkl_trig_transforms

        implicit none

        integer n, k, tt_type
        integer ir, ipar(128)
        integer :: i,j,jj
        real :: y,yh,yph,phi,sech2
        real :: x,xph,rer,ell,aap 
        double precision :: spar(3*jmax/2+2)
        type(dfti_descriptor), pointer :: handle

        dx = wlength/float(imx)
        dy = width/float(jmax)

! **** hyperviscosity ****
        damp = (dx**ndeg)/(pi**ndeg)
        damp = dampcoeff*damp/dt
        write(6,897) damp
 897 format('damp =',1p1e13.5)        

!  **** Initial zonal-mean state profile ****
        do j = 1,jmax+1
         y = dy*float(j-1)-0.5*width
         phi = y/sigma
         ueq(j) = 40.*(1.-tanh(phi)*tanh(phi))
         sech2 = (1.-tanh(phi)*tanh(phi))
         
!        ueq(j) = 0.
!        sech2 = 0.

         qby_1(j) = -40.*(2./(sigma*sigma))*sech2*(2.*tanh(phi)*tanh(phi)   &
                   -sech2)+(1./(rd*rd))*ueq(j)
         qby_2(j) = -(1./(rd*rd))*ueq(j) 

!        qy_1(j) = 0.
!        qy_2(j) = 0.

        enddo 
!        u0 = ueq(1)
        write(6,*) '   '
        do j=1,jmax+1
          ueq(j) = ueq(j)-u0
           write(6,*) j,"test", ueq(j), beta+qby_1(j),beta+qby_2(j)
        end do

!  **** Transform zonal-mean to spectral space ***
        tt_type=0     ! sine transform

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(qby_1,handle,ipar,spar,ir)
        CALL D_FORWARD_TRIG_TRANSFORM(qby_1,handle,ipar,spar,ir)

        CALL D_INIT_TRIG_TRANSFORM(jmax,tt_type,ipar,spar,ir)
        CALL D_COMMIT_TRIG_TRANSFORM(qby_2,handle,ipar,spar,ir)
        CALL D_FORWARD_TRIG_TRANSFORM(qby_2,handle,ipar,spar,ir)
        
        qymean1(:,1) = 0.
        qymean2(:,1) = 0.
      do j = 2,nmax
        qymean1(j,1) = qby_1(j) 
        qymean2(j,1) = qby_2(j) 
      enddo
        qymean1(:,2) = qymean1(:,1)
        qymean1(:,3) = qymean1(:,2)
        qymean1(:,4) = qymean1(:,3)
        qymean2(:,2) = qymean2(:,1)
        qymean2(:,3) = qymean2(:,2)
        qymean2(:,4) = qymean2(:,3)

        psi_1 = zero
        vort_1 = zero
        psi_2 = zero
        vort_2 = zero
        
!   *** flow (uniform PV perturbation, symmetric mode only) ****
!  
       do j = 1,nmax/2
!      do j = 1,1
         jj = 2*j-1
         ell = el * float(jj) 
         aap = ((-1)**(j+1))*exp(-ell*ell*sigma*sigma)
       do i = 2,9
           vort_1(i,2*j,1) = ur*amp*aap/float(nmax/2)
           vort_2(i,2*j,1) = ur*amp*aap/float(nmax/2)
       enddo
       enddo

           vort_1(:,:,2) = vort_1(:,:,1)
           vort_1(:,:,3) = vort_1(:,:,1)
           vort_2(:,:,2) = vort_2(:,:,1)
           vort_2(:,:,3) = vort_2(:,:,1)

      return
      end subroutine
    
end module

