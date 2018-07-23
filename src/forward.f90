!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rearrange arrays forward in time
! Newest index is in position 1, will write to it next
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module forward
contains
subroutine iterate
  use global_variables
  implicit none
  vort_1(:,:,1) = vort_1(:,:,2)
  vort_1(:,:,2) = vort_1(:,:,3)
  vort_1(:,:,3) = vort_1(:,:,4)
  vort_2(:,:,1) = vort_2(:,:,2)
  vort_2(:,:,2) = vort_2(:,:,3)
  vort_2(:,:,3) = vort_2(:,:,4)
  qymean1(:,1)  = qymean1(:,2)
  qymean1(:,2)  = qymean1(:,3)
  qymean1(:,3)  = qymean1(:,4)
  qymean2(:,1)  = qymean2(:,2)
  qymean2(:,2)  = qymean2(:,3)
  qymean2(:,3)  = qymean2(:,4)
  adv_1(:,:,1)  = adv_1(:,:,2)
  adv_1(:,:,2)  = adv_1(:,:,3)
  adv_2(:,:,1)  = adv_2(:,:,2)
  adv_2(:,:,2)  = adv_2(:,:,3)
  vqm_1(:,1)    = vqm_1(:,2)
  vqm_1(:,2)    = vqm_1(:,3)
  vqm_2(:,1)    = vqm_2(:,2)
  vqm_2(:,2)    = vqm_2(:,3)
  return
end subroutine
end module
