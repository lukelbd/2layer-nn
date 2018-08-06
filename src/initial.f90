!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read namelist and transform units
! Initialize Fourier descriptor handles
! Initial flow field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module initial
contains
subroutine init
  use global_variables
  use transforms
  use mkl_dfti
  implicit none
  integer :: i, j, l(1)
  real :: x, y, fact, offset, offset_sp, offset_wll
  real :: r(idft,jmax) ! for seeding
  l(1) = imax ! array specifying length

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    ---- Declare namelist ----
  namelist /input_nml/ &
    dt, dt_io, days, days_spinup, &
    width, wlength, rd, shear, beta, &
    tau_f, tau_r, tau_sp, y_sp, &
    ndeg, visc, &
    amp_i, tau_i, y_i, wmin_i, wmax_i, &
    ll_seed_on, ll_seed_amp

  !    ---- Read namelist ----
  ! Will overwrite the above defaults **only if** the variables appear there
  open(1, file='input.nml', status='old', action='read') ! status='old' requires that file exists, wut
  read(1, nml=input_nml) ! note 'input.nml' is a filename, input_nml is a declared variable
  close(1) ! close file

  !    ---- In-place unit scaling ----
  ! Note want time steps to be integer because we iterate over them
  ! and test for equality sometimes
  t_end    = int(days*3600.*24.)        ! days to s
  t_spinup = int(days_spinup*3600.*24.) ! days to s
  tau_r    = tau_r*24.*3600.            ! days to s
  tau_f    = tau_f*24.*3600.            ! days to s
  tau_sp   = tau_sp*24.*3600.           ! days to s
  rd       = rd*1.e3                    ! km to m
  width    = width*1.e3                 ! km to m
  wlength  = wlength*1.e3               ! km to m
  y_i      = 0.5*width*y_i              ! to m

  !    ---- Calculate dependent variables ----
  dx = wlength/float(imax) ! (m) grid resolution in x
  dy = width/float(jmax-1) ! (m) grid resolution in y; non-cyclic so use -1, or something
  el = pi/width            ! (1/m) converts radians to meridional wavenums
  rk = 2.*pi/wlength       ! (1/m) converts radians to zonal wavenums
  damp = visc*(dx**ndeg)/(dt*(pi**ndeg))  ! (m^ndeg*s^-1 hyperviscocity

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    ---- Initialize handles for Fourier calls ----
  ! Arg 1 is handle, Arg 2 is precision of transform, Arg 3 is
  ! 'forward domain' (i.e. the 'forward' transform takes place on only
  ! real-valued sequences), Arg 4/5 are the number of dims/dimension lengths
  ret = DftiCreateDescriptor(hcr, dfti_double, dfti_real, 1, l)
  ret = DftiSetValue(hcr, dfti_placement, dfti_not_inplace)
  ret = DftiCommitDescriptor(hcr)
  ret = DftiCreateDescriptor(hrc, dfti_double, dfti_real, 1, l)
  ret = DftiSetValue(hrc, dfti_placement, dfti_not_inplace)
  ret = DftiCommitDescriptor(hrc)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    ---- Set initial anomalous mean-pv gradient to zero ----
  qybar1_tt = 0.0
  qybar2_tt = 0.0

  !    ---- Set initial anomalous pv to zero ----
  q1_sp = zero
  q2_sp = zero

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
    ! q1_sp(:,:,1) = q1_sp(:,:,1) + (r-0.5) * ll_seed_amp ! upper layer seed
    q2_sp(:,:,1) = q2_sp(:,:,1) + (r-0.5) * ll_seed_amp
  endif 

  !    ---- Vectors x/y in physical units ('zero' is center of channel/left boundary) ----
  do i = 1,imax
    x_cart(i) = float(i-1)*dx ! from left boundary
  enddo
  do j = 1,jmax
    y_cart(j) = float(j-1)*dy - 0.5*width ! from channel center
  enddo
  print *,y_cart

  !    ---- Mask for upper layer pv injection ----
  ! Simple Gaussian curve
  do j = 1,jmax
    y = y_cart(j) ! distance from center, in physical units
    if (abs(y) >= 3*y_i) then
      mask_i(j) = 0.0
    else
      mask_i(j) = exp(-y*y/(y_i*y_i)) ! weighting to apply to pv injection anomalies
    endif
  enddo

  !    ---- Mask for sponge layer ----
  ! Will be a quadratic scale factor up to each wall
  offset_sp  = float(jmax-1)*0.5*(1.0-y_sp) ! distance from center point in grid cell units
  offset_wll = float(jmax-1)*0.5            ! the wall location
  do j = 1,jmax
    offset     = abs(float(j-1) - 0.5*float(jmax-1)) ! distance from center point (think about it, with jmax=5)
    fact       = max(0.0, (offset-offset_sp) / (offset_wll-offset_sp)) ! i.e. zero damping in center
    mask_sp(j) = fact*fact ! multiplier to apply to pv anomalies, in 1/seconds
  enddo

  !    ---- Mask in spectral space ----
  ! Want a cosine transform, since mask at top/bottom boundary
  ! is certainly not zero!
  tt_type = 1 ! cosines
  call ftt(mask_sp, mask_sp_tt, tt_type, jmax)

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

