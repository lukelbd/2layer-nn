!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve prognostic vorticity equation in spectral coordinates
! List of parameters from diagnostics.f90 that this function reads from:
!  * adv_[12]_c:  the spectral advection terms for (anomalous?) pv
!  * q[12]_c:     the (anomalous?) pv
!  * qyyflux[12]: the total meridional pv flux convergence (eddy-->mean)
!  * umean[12]:   the zonal mean wind
! These are all ***transformed*** values; trigonometric in y, Fourier in x.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module prognostics
contains
subroutine prog
  use global_variables

  implicit none
  integer :: i,j,i1,j1,n1,m1,ii,jj
  real :: ell,rkk,kkk,filter

  ! Compute energy
  energy2 = 0.5*sum(vxy1*vxy1 + vxy2*vxy2 +  &
    uxy1*uxy1 + uxy2*uxy2)/float(imax*jmax)

  ! Populate older positions with current position
  if(t.eq.tstart) then
    adv1_c(:,:,2) = adv1_c(:,:,3)
    adv1_c(:,:,1) = adv1_c(:,:,2)
    adv2_c(:,:,2) = adv2_c(:,:,3)
    adv2_c(:,:,1) = adv2_c(:,:,2)
    qyyflux1(:,2) = qyyflux1(:,3)
    qyyflux1(:,1) = qyyflux1(:,2)
    qyyflux2(:,2) = qyyflux2(:,3)
    qyyflux2(:,1) = qyyflux2(:,2)
  endif

  ! Solve prognostic vorticity in spectral coordinates
  do j = 2,jtrunc
    ell = el*float(j-1)
    do i = 2,itrunc
      rkk = rk*float(i-1)
      ! Get *hyperdiffusion*, *radiation*, and *friction* terms
      visc1_c(i,j) = -damp*((rkk**2+ell**2)**3)*q1_c(i,j,3)
      visc2_c(i,j) = -damp*((rkk**2+ell**2)**3)*q2_c(i,j,3)
      rad1_c(i,j)  = -(psi2_c(i,j)-psi1_c(i,j))/((rd*rd)*tau_r) ! relax toward default shear
      rad2_c(i,j)  = (psi2_c(i,j)-psi1_c(i,j))/((rd*rd)*tau_r)  ! relax toward default shear
      fric_2(i,j) = (rkk**2+ell**2)*psi2_c(i,j)/tau_f ! relative vorticity; we subtract (-du/dy), but add (+dv/dx)?
      ! fric_2(i,j) = fric_2(i,j)-q2_c(i,j,3)/tau_2
      ! Apply *hyperdiffusion*, *radiation*, and *pv injection* in upper layer
      q1_c(i,j,4) = q1_c(i,j,3) &
        + dt*(23.*adv1_c(i,j,3)-16.*adv1_c(i,j,2)+5.*adv1_c(i,j,1))/12. &
        + dt*visc1_c(i,j) &
        + dt*force1_c(i,j) &
        + dt*rad1_c(i,j)
      ! Apply *hyperdiffusion*, *radiation*, and *friction* in lower layer
      q2_c(i,j,4) = q2_c(i,j,3) &
        + dt*(23.*adv2_c(i,j,3)-16.*adv2_c(i,j,2)+5.*adv2_c(i,j,1))/12. &
        + dt*visc2_c(i,j) &
        + dt*rad2_c(i,j) &
        + dt*fric_2(i,j)
    enddo
    ! Appply truncation
    ! Is this really necessary? Never read from those positions anyway.
    q1_c(1,j,4) = zero
    q2_c(1,j,4) = zero
    q1_c(itrunc+1:idft,j,4) = zero
    q2_c(itrunc+1:idft,j,4) = zero

    !    ---- Prognostic meridional derivative of zonal mean vorticity ----
    ! Iterate, apply *hyperdiffusion* and *radiation* in upper layer
    qymean1(j,4) = qymean1(j,3) &
      - dt*(23.*qyyflux1(j,3)-16.*qyyflux1(j,2)+5.*qyyflux1(j,1))/12. &
      - dt*damp*(ell**6)*qymean1(j,3) &
      - dt*(umean1(j)-umean2(j))/(rd*rd*tau_r)
    ! Iterate, apply *hyperdiffusion*, *radiation*, and *friction* in lower layer
    qymean2(j,4) = qymean2(j,3) &
      - dt*(23.*qyyflux2(j,3)-16.*qyyflux2(j,2)+5.*qyyflux2(j,1))/12. &
      - dt*damp*(ell**6)*qymean2(j,3) &
      + dt*(umean1(j)-umean2(j))/(rd*rd*tau_r) &
      - dt*(ell*ell*umean2(j))/tau_f
  enddo

  !    ---- Apply truncation ----
  ! Is this really necessary? Never read from those positions anyway.
  q1_c(:,jtrunc+1:jmax,4) = zero
  q2_c(:,jtrunc+1:jmax,4) = zero
  qymean1(jtrunc+1:jmax,4) = 0.
  qymean2(jtrunc+1:jmax,4) = 0.

  return
end subroutine
end module
