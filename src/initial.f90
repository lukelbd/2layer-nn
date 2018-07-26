!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial flow field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module initial
contains
subroutine init
  use global_variables
  implicit none
  integer :: n, k, tt_type
  integer :: ir, ipar(128)
  integer :: i, j, jj
  real :: y, yh,  yph, phi, sech2
  real :: x, xph, rer, ell, aap
  real :: offset, offset_sp, offset_wall, fact
  real :: r(idft,jmax)

  !    ---- Zonal mean pv ----
  ! Since initially the mean PV is spatially uniform and equal to
  ! beta +/- shear/(rd*rd), its initial meridional gradient is set to zero.
  qybar1_tt = 0.0
  qybar2_tt = 0.0

  !    ---- Anomalous pv ----
  ! The initial pv is initially set to zero
  q1_sp = zero
  q2_sp = zero

  !    ---- Y in physical coordinates ('zero' is center of channel)  ----
  do j = 1,jmax
    y_cart(j) = float(j-1)*dy - 0.5*width 
  enddo

  !    ---- Basic state jet ----
  if (init_jet) then
    do j = 1,jtrunc/2
      jj = 2*j-1
      ell = el * float(jj)
      aap = ((-1)**(j+1))*exp(-ell*ell*sigma_jet*sigma_jet)
      do i = 2,9
        q1_sp(i,2*j,1) = one_r*amp_jet*aap/float(jtrunc/2)
        q2_sp(i,2*j,1) = one_r*amp_jet*aap/float(jtrunc/2)
      enddo
    enddo
  endif

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
  offset_sp   = float(jmax-1)*0.5*(1.0-y_sp) ! distance from center point in grid cell units
  offset_wall = float(jmax-1)*0.5            ! the wall location
  do j = 1,jmax
    offset     = abs(float(j-1) - 0.5*float(jmax-1)) ! distance from center point (think about it, with jmax=5)
    fact       = (offset-offset_sp) / (offset_wall-offset_sp)
    mask_sp(j) = -(1.0/tau_sp)*fact*fact ! multiplier to apply to pv anomalies, in 1/seconds
  enddo

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

