! ***** Prognosis of vorticity *****

module EXEC
contains

    subroutine prog

        use GLOBAL_VARIABLES
        use mkl_dfti
        use mkl_trig_transforms

     implicit none
     integer :: i,j,i1,j1,n1,m1,ii,jj
     real :: ell,rkk,kkk,filter
     

   ! ***** compute energy *******
   
      energy2 = 0.5*sum(vxy1*vxy1 + vxy2*vxy2 +  &
                    uxy1*uxy1 + uxy2*uxy2)/float(imx*(jmax+1))
 
    ! ***** evaluate advection terms and viscosity in the specral space ****

!    write(6,*) 'test ',maxval(abs(vqm_1(:,3))),maxval(abs(vqm_2(:,3)))
        
   
!  do j = 1,jmax+1
!    write(6,*) 'test test test ',j,vqm_1(j,3),vqm_1(j,2),vqm_1(j,1)
!  enddo

   do j = 2,nmax
            ell = el*float(j-1)
!    if(m.eq.mstart) then
!      vqm_1(j,2) = vqm_1(j,3)
!      vqm_1(j,1) = vqm_1(j,2)
!      vqm_2(j,2) = vqm_2(j,3)
!      vqm_2(j,1) = vqm_2(j,2)
!    endif

   do i = 2,mmax
            rkk = rk*float(i-1)
     visc_1(i,j) = -damp*((rkk**2+ell**2)**3)*vort_1(i,j,3)
     visc_2(i,j) = -damp*((rkk**2+ell**2)**3)*vort_2(i,j,3)
     if(m.eq.mstart) then
       adv_1(i,j,2) = adv_1(i,j,3)
       adv_1(i,j,1) = adv_1(i,j,2)
       adv_2(i,j,2) = adv_2(i,j,3)
       adv_2(i,j,1) = adv_2(i,j,2)
     endif
     vort_1(i,j,4) = vort_1(i,j,3) +     &
dt*(23.*adv_1(i,j,3)-16.*adv_1(i,j,2)+5.*adv_1(i,j,1))/12.
     vort_1(i,j,4) = vort_1(i,j,4) + dt*visc_1(i,j)
     vort_2(i,j,4) = vort_2(i,j,3) +     &
dt*(23.*adv_2(i,j,3)-16.*adv_2(i,j,2)+5.*adv_2(i,j,1))/12.
     vort_2(i,j,4) = vort_2(i,j,4) + dt*visc_2(i,j)
   enddo

   do i = mmax+1,imax
    vort_1(i,j,4) = (0.,0.)
    vort_2(i,j,4) = (0.,0.)
   enddo
   do i = 1,1
    vort_1(i,j,4) = (0.,0.)
    vort_2(i,j,4) = (0.,0.)
   enddo

    qymean1(j,4) = qymean1(j,3) -  &
dt*(23.*vqm_1(j,3)-16.*vqm_1(j,2)+5.*vqm_1(j,1))/12.
    qymean1(j,4) = qymean1(j,4) - dt*damp*(ell**6)*qymean1(j,3)
    qymean2(j,4) = qymean2(j,3) -  &
dt*(23.*vqm_2(j,3)-16.*vqm_2(j,2)+5.*vqm_2(j,1))/12.
    qymean2(j,4) = qymean2(j,4) - dt*damp*(ell**6)*qymean2(j,3)
!   qymean1(j,4) = qymean1(j,4)
!   qymean2(j,4) = qymean2(j,4)
   enddo

   do j = nmax+1,jmax+1
    vort_1(:,j,4) = (0.,0.)
    vort_2(:,j,4) = (0.,0.)
    qymean1(j,4) = 0.
    qymean2(j,4) = 0.
   enddo

    return
    end subroutine
end module
