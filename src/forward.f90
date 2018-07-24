!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rearrange arrays forward in time
! * The last index (4/3) is timestep after running through
!   *prognostics.f90*. Immediately shifted by *forward.f90* to n-1 position.
! * Therefore, in diagnostics.f90, we always work with the n-1 (3/2) data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module forward
contains
subroutine iterate
  use global_variables
  implicit none
  q1_c(:,:,1) = q1_c(:,:,2)
  q1_c(:,:,2) = q1_c(:,:,3)
  q1_c(:,:,3) = q1_c(:,:,4)
  q2_c(:,:,1) = q2_c(:,:,2)
  q2_c(:,:,2) = q2_c(:,:,3)
  q2_c(:,:,3) = q2_c(:,:,4)
  qybar1_out(:,1)  = qybar1_out(:,2)
  qybar1_out(:,2)  = qybar1_out(:,3)
  qybar1_out(:,3)  = qybar1_out(:,4)
  qybar2_out(:,1)  = qybar2_out(:,2)
  qybar2_out(:,2)  = qybar2_out(:,3)
  qybar2_out(:,3)  = qybar2_out(:,4)
  adv1_c(:,:,1)  = adv1_c(:,:,2)
  adv1_c(:,:,2)  = adv1_c(:,:,3)
  adv2_c(:,:,1)  = adv2_c(:,:,2)
  adv2_c(:,:,2)  = adv2_c(:,:,3)
  qyyflux1_out(:,1)    = qyyflux1_out(:,2)
  qyyflux1_out(:,2)    = qyyflux1_out(:,3)
  qyyflux2_out(:,1)    = qyyflux2_out(:,2)
  qyyflux2_out(:,2)    = qyyflux2_out(:,3)
  return
end subroutine
end module
