!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rearrange arrays forward in time
! * Store the previous *two* advection timesteps for fancy schmancy
!   advection algorithm.
! * The model solves for q (zonal q anomalies) and qybar (the meridional
!   gradient in zonal mean q); just need 'previous' and 'current' timesteps there.
! * Zonal means are trangular transforms in y, and 2D resolution are 
!   Fourier transformed and triangular transformed.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module forward
contains
subroutine iterate
  use global_variables
  implicit none
  q1_sp(:,:,1)       = q1_sp(:,:,2)
  q2_sp(:,:,1)       = q2_sp(:,:,2)
  qybar1_tt(:,1)     = qybar1_tt(:,2)
  qybar2_tt(:,1)     = qybar2_tt(:,2)
  adv1_sp(:,:,1)     = adv1_sp(:,:,2)
  adv1_sp(:,:,2)     = adv1_sp(:,:,3)
  adv2_sp(:,:,1)     = adv2_sp(:,:,2)
  adv2_sp(:,:,2)     = adv2_sp(:,:,3)
  qyyflux1_tt(:,1)   = qyyflux1_tt(:,2)
  qyyflux1_tt(:,2)   = qyyflux1_tt(:,3)
  qyyflux2_tt(:,1)   = qyyflux2_tt(:,2)
  qyyflux2_tt(:,2)   = qyyflux2_tt(:,3)
  force1_cart(:,:,1) = force1_cart(:,:,2)
end subroutine
end module
