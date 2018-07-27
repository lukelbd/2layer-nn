!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial flow field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module initial
contains
subroutine init
  use spectral
  use global_variables
  implicit none
  integer :: i, j, jj
  real :: y, fact, offset, offset_sp, offset_wll
  real :: r(idft,jmax) ! for seeding

  !    ---- Set initial anomalous mean-pv gradient to zero ----
  qybar1_tt = 0.0
  qybar2_tt = 0.0

  !    ---- Set initial anomalous pv to zero ----
  q1_sp = zero
  q2_sp = zero

  !    ---- Y in physical coordinates ('zero' is center of channel)  ----
  do j = 1,jmax
    y_cart(j) = float(j-1)*dy - 0.5*width 
  enddo

  !    ---- Initial jet ----
  ! * If we want a perpetual jet background state, will have to create a new
  !   variable *separate* from the tracked qybar array, and add advection
  !   of vorticity due to the jet/extra wind from the jet when calculating
  !   advective terms in diagnostics.f90 -- not trivial!
  ! * Will want to define zonal jet wind on *cartesian* grid, then transform to
  !   spectral space and un-invert the dqbar/dy equation to get the jet's contribution
  !   to dqbar/dy (we already do this with zonal mean winds, so not too hard).

  !    ---- Random initial noise in lower layer ----
  ! Necessary to get eddies going, otherwise (without pv injection in
  ! upper layer) will just stay symmetric forever
  if (ll_seed_on) then
    call random_number(r)
    q2_sp(:,:,1) = q2_sp(:,:,1) + (r-0.5) * ll_seed_amp
  endif 

  !    ---- Mask for upper layer pv injection ----
  ! Will be a simple Gaussian curve
  ! The 'i' is for 'injection'
  mask_i(:) = 0.0
  do j = 1,jmax
    y = y_cart(j) ! distance from center, in physical units
    if (abs(y) >= 3*sigma_i) then
      mask_i(j) = 0.0
    else
      mask_i(j) = exp(-y*y/(sigma_i*sigma_i)) ! weighting to apply to pv injection anomalies
    endif
  enddo

  !    ---- Mask for sponge layer ----
  ! Will be a quadratic scale factor up to each wall
  mask_sp(:) = 0.0
  offset_sp  = float(jmax-1)*0.5*(1.0-y_sp) ! distance from center point in grid cell units
  offset_wll = float(jmax-1)*0.5            ! the wall location
  do j = 1,jmax
    offset     = abs(float(j-1) - 0.5*float(jmax-1)) ! distance from center point (think about it, with jmax=5)
    fact       = (offset-offset_sp) / (offset_wll-offset_sp)
    mask_sp(j) = -(1.0/tau_sp)*fact*fact ! multiplier to apply to pv anomalies, in 1/seconds
  enddo

  !    ---- Sponge layer in spectral space ----
  ! Want a cosine transform, since mask at top/bottom boundary
  ! is certainly not zero!
  tt_type = 1 ! cosine
  call ftt(mask_sp, mask_sp_tt, tt_type, jmax, jtrunc) ! jtrunc is ignored for now

  !    ---- Test ----
  ! do i = 6,6
  !    do j = 2,2
  !      q1_sp(i,j,:) = one_r*amp
  !      q2_sp(i,j,:) = -one_r*amp
  !    enddo
  ! enddo
  return
end subroutine
end module

