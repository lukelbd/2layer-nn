!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve prognostic vorticity equation in spectral coordinates
! List of parameters from diagnostics.f90 that this function reads from:
!  * adv_[12]_c:  the spectral advection terms for (anomalous?) pv
!  * q[12]_c:     the (anomalous?) pv
!  * qyyflux_out[12]: the total meridional pv flux convergence (eddy-->mean)
!  * ubar_out[12]:   the zonal mean wind
! These are all ***transformed*** values; trigonometric in y, Fourier in x.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module prognostics
contains
subroutine prog
  use global_variables

  implicit none
  integer :: i,j,i1,j1,n1,m1,ii,jj
  real :: ell,rkk,kkk,filter

  ! Populate older positions with current position
  if(t.eq.tstart) then
    adv1_c(:,:,2) = adv1_c(:,:,3)
    adv1_c(:,:,1) = adv1_c(:,:,2)
    adv2_c(:,:,2) = adv2_c(:,:,3)
    adv2_c(:,:,1) = adv2_c(:,:,2)
    qyyflux1_out(:,2) = qyyflux1_out(:,3)
    qyyflux1_out(:,1) = qyyflux1_out(:,2)
    qyyflux2_out(:,2) = qyyflux2_out(:,3)
    qyyflux2_out(:,1) = qyyflux2_out(:,2)
  endif

  do j = 2,jtrunc
    !    ---- Prognostic zonal pv anomalies ----
    ell = el*float(j-1)
    do i = 2,itrunc
      rkk = rk*float(i-1)
      ! Get *hyperdiffusion*, *radiation*, and *friction* terms
      ! Note vorticity is del^2(psi)
      ! Note psi are *anomalous* relative to background shear, so want
      ! radiation to relax their *difference* to zero
      visc1_c(i,j) = -damp*((rkk**2+ell**2)**3)*q1_c(i,j,3) ! hyperdiffusion
      visc2_c(i,j) = -damp*((rkk**2+ell**2)**3)*q2_c(i,j,3)
      rad1_c(i,j)  = -(psi2_c(i,j)-psi1_c(i,j))/((rd*rd)*tau_r) ! radiation
      rad2_c(i,j)  = (psi2_c(i,j)-psi1_c(i,j))/((rd*rd)*tau_r)
      fric_2(i,j)  = (rkk**2+ell**2)*psi2_c(i,j)/tau_f ! friction (equals negative of the vorticity anomaly)
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

    !    ---- Prognostic meridional derivative of zonal mean pv ----
    ! Iterate, apply *hyperdiffusion* and *radiation* in upper layer
    qybar1_out(j,4) = qybar1_out(j,3) &
      - dt*(23.*qyyflux1_out(j,3)-16.*qyyflux1_out(j,2)+5.*qyyflux1_out(j,1))/12. &
      - dt*damp*(ell**6)*qybar1_out(j,3) &
      - dt*(ubar1_out(j)-ubar2_out(j))/(rd*rd*tau_r)
    ! Iterate, apply *hyperdiffusion*, *radiation*, and *friction* in lower layer
    qybar2_out(j,4) = qybar2_out(j,3) &
      - dt*(23.*qyyflux2_out(j,3)-16.*qyyflux2_out(j,2)+5.*qyyflux2_out(j,1))/12. &
      - dt*damp*(ell**6)*qybar2_out(j,3) &
      + dt*(ubar1_out(j)-ubar2_out(j))/(rd*rd*tau_r) &
      - dt*(ell*ell*ubar2_out(j))/tau_f
  enddo

  !    ---- Apply truncation ----
  ! Is this really necessary? Never read from those positions anyway.
  q1_c(:,jtrunc+1:jmax,4) = zero
  q2_c(:,jtrunc+1:jmax,4) = zero
  qybar1_out(jtrunc+1:jmax,4) = 0.
  qybar2_out(jtrunc+1:jmax,4) = 0.

  return
end subroutine
end module
