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

    ! ------ Compute energy -------
    energy2 = 0.5*sum(vxy1*vxy1 + vxy2*vxy2 +  &
      uxy1*uxy1 + uxy2*uxy2)/float(imx*(jmax+1))

    do j = 2,nmax
    ell = el*float(j-1)
    !    if(t.eq.tstart) then
    !      vqm_1(j,2) = vqm_1(j,3)
    !      vqm_1(j,1) = vqm_1(j,2)
    !      vqm_2(j,2) = vqm_2(j,3)
    !      vqm_2(j,1) = vqm_2(j,2)
    !    endif

    do i = 2,mmax
    rkk = rk*float(i-1)
    visc_1(i,j) = -damp*((rkk**2+ell**2)**3)*vort_1(i,j,3)
    visc_2(i,j) = -damp*((rkk**2+ell**2)**3)*vort_2(i,j,3)
    rad_1(i,j) = -(psi_2(i,j)-psi_1(i,j))/((rd*rd)*tau_r)
    rad_2(i,j) = (psi_2(i,j)-psi_1(i,j))/((rd*rd)*tau_r)
    fric_2(i,j) = (rkk**2+ell**2)*psi_2(i,j)/tau_f
    !    fric_2(i,j) = fric_2(i,j)-vort_2(i,j,3)/tau_2
    if(t.eq.tstart) then
      adv_1(i,j,2) = adv_1(i,j,3)
      adv_1(i,j,1) = adv_1(i,j,2)
      adv_2(i,j,2) = adv_2(i,j,3)
      adv_2(i,j,1) = adv_2(i,j,2)
    endif
    vort_1(i,j,4) = vort_1(i,j,3) +     &
      dt*(23.*adv_1(i,j,3)-16.*adv_1(i,j,2)+5.*adv_1(i,j,1))/12.
    vort_1(i,j,4) = vort_1(i,j,4) + dt*visc_1(i,j)
    vort_1(i,j,4) = vort_1(i,j,4) + dt*rad_1(i,j)
    vort_2(i,j,4) = vort_2(i,j,3) +     &
      dt*(23.*adv_2(i,j,3)-16.*adv_2(i,j,2)+5.*adv_2(i,j,1))/12.
    vort_2(i,j,4) = vort_2(i,j,4) + dt*visc_2(i,j)
    vort_2(i,j,4) = vort_2(i,j,4) + dt*rad_2(i,j)
    vort_2(i,j,4) = vort_2(i,j,4) + dt*fric_2(i,j)
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
    qymean1(j,4) = qymean1(j,4) - dt*(umean1(j)-umean2(j))/(rd*rd*tau_r)
    qymean2(j,4) = qymean2(j,3) -  &
      dt*(23.*vqm_2(j,3)-16.*vqm_2(j,2)+5.*vqm_2(j,1))/12.
    qymean2(j,4) = qymean2(j,4) - dt*damp*(ell**6)*qymean2(j,3)
    qymean2(j,4) = qymean2(j,4) + dt*(umean1(j)-umean2(j))/(rd*rd*tau_r)
    qymean2(j,4) = qymean2(j,4) - dt*(ell*ell*umean2(j))/tau_f
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
