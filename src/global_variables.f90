module global_variables
  implicit none
  !----------------------------------------------------------------------------!
  !    ---- Initial stuff, and dimensions that must be hard-coded ----
  real, parameter :: pi=3.141592653589793
  complex :: one_r = (1.,0.), one_i = (0.,1.), zero = (0.,0.)
  integer :: t         ! time tracker
  real :: energy2, cfl ! for monitoring integration
  ! Timing
  integer, parameter :: tstart=0
  ! Grid resolution (number of cells)
  integer, parameter :: imax=512, jmax=256
  ! integer, parameter :: imax=256, jmax=256
  ! Number complex coefficients for cyclic DFT in x-direction
  ! Signal is real, so need only N/2-1 complex coeffs plus the A_0 and A_(N/2) coeffs
  ! For y sine transform, coefficients are all real; need all <jmax> of them
  ! See: https://en.wikipedia.org/wiki/Discrete_Fourier_transform#The_real-input_DFT
  integer, parameter :: idft=1+imax/2
  ! Truncation of y-sine ft and x-dft to resolve (avoid aliasing) quadratic terms
  ! See: 10.1098/rsta.2013.0289
  ! Requirement is: imax >= (3*xtrunc+1), jmax >= (3*ytrunc+1)/2
  integer, parameter :: itrunc=170, jtrunc=170
  ! integer, parameter :: itrunc=85, jtrunc=170

  !    ---- To be calculated in subroutine below ----
  real :: damp, dx, dy, el, rk

  !    ---- To be supplied by namelist ----
  logical :: ll_seed_on, init_jet         ! flags
  integer :: dt, td                       ! time steps
  real :: tend, tchange, tds              ! timing
  real :: width, wlength                  ! channel width, channel length
  real :: rd                              ! radius of deformation
  real :: tau_r, tau_f, tau_2, tau_sponge ! damping stuff
  real :: u0, beta                        ! background state
  real :: amp_jet, sigma_jet              ! initial jet
  real :: amp_i, tau_i, sigma_i, wmin_i, wmax_i ! Noboru's stoshcastic forcing
  real :: ll_seed_amp                     ! lower level forcing
  real :: visc, ndeg                      ! hyperviscocity

  !    ---- Noboru's arrays ----
  ! Arrays storing Fourier coefficients from x-direction decomposition
  complex :: vort_1(idft,jmax+1,4), psi_1(idft,jmax+1),  &
         adv_1(idft,jmax+1,3), u_1(idft,jmax+1),   v_1(idft,jmax+1), &
         visc_1(idft,jmax+1),  rad_1(idft,jmax+1), force_1(idft,jmax+1), &
         q_1(idft,jmax+1),     qx_1(idft,jmax+1),  qy_1(idft,jmax+1)
  complex :: vort_2(idft,jmax+1,4), psi_2(idft,jmax+1),  &
         adv_2(idft,jmax+1,3), u_2(idft,jmax+1),   v_2(idft,jmax+1), &
         visc_2(idft,jmax+1),  rad_2(idft,jmax+1), fric_2(idft,jmax+1), &
         q_2(idft,jmax+1),     qx_2(idft,jmax+1),  qy_2(idft,jmax+1)

  ! Special arrays
  real :: ueq(jmax+1), uf(imax,jmax+1), ymask(jmax+1)

  ! Arrays storing real and imaginary components of Fourier decomposition
  real :: u_1_r(idft,jmax+1),  u_1_i(idft,jmax+1),  &
          v_1_r(idft,jmax+1),  v_1_i(idft,jmax+1),  &
          q_1_r(idft,jmax+1),  q_1_i(idft,jmax+1),  &
          f_1_r(idft,jmax+1),  f_1_i(idft,jmax+1),  &
          qx_1_r(idft,jmax+1), qx_1_i(idft,jmax+1), &
          qy_1_r(idft,jmax+1), qy_1_i(idft,jmax+1), &
          p_1_r(idft,jmax+1),  p_1_i(idft,jmax+1)
  real :: u_2_r(idft,jmax+1),  u_2_i(idft,jmax+1),  &
          v_2_r(idft,jmax+1),  v_2_i(idft,jmax+1),  &
          q_2_r(idft,jmax+1),  q_2_i(idft,jmax+1),  &
          qx_2_r(idft,jmax+1), qx_2_i(idft,jmax+1), &
          qy_2_r(idft,jmax+1), qy_2_i(idft,jmax+1), &
          p_2_r(idft,jmax+1),  p_2_i(idft,jmax+1)

  ! Derivatives and other diagnostics in real space
  real :: vqm_1(jmax+1,3), fxy1(imax,jmax+1,2), &
          uxy1(imax,jmax+1), vxy1(imax,jmax+1), pxy1(imax,jmax+1), &
          qxy1(imax,jmax+1), qxx1(imax,jmax+1), qyy1(imax,jmax+1)
  real :: vqm_2(jmax+1,3), &
          uxy2(imax,jmax+1), vxy2(imax,jmax+1), pxy2(imax,jmax+1), &
          qxy2(imax,jmax+1), qxx2(imax,jmax+1), qyy2(imax,jmax+1)
  ! Special; stores mean and mean gradient
  real :: umean1(jmax+1), qymean1(jmax+1,4), &
          umean2(jmax+1), qymean2(jmax+1,4)
  ! Other stuff that has to be double precision for some reason
  double precision :: qby_1(jmax+1),qbar1(jmax+1),ubar1(jmax+1)
  double precision :: qby_2(jmax+1),qbar2(jmax+1),ubar2(jmax+1)

  contains

  subroutine read_namelist
    implicit none
    !    ---- Declare namelist ----
    ! real :: tau_farray(jmax)
    namelist /input_nml/ &
      ll_seed_on, init_jet, &
      width, wlength, &
      dt, td, &
      tend, tchange, tds, &
      tau_r, tau_f, tau_sponge, tau_2, &
      u0, beta, rd, &
      amp_jet, sigma_jet, &
      amp_i, tau_i, sigma_i, wmin_i, wmax_i, &
      ll_seed_amp, &
      visc, ndeg
    !    ---- Read namelist ----
    open(1, file='input.nml', status='old', action='read') ! status='old' requires that file exists, wut
    read(1, nml=input_nml)
      ! note 'input.nml' is a filename, input_nml is a declared variable
    close(1) ! close file

    !    ---- In-place unit scaling ----
    tend      = tend*3600.*24.    ! days to s
    tchange   = tchange*3600.*14. ! days to s
    tds       = tds*3600.*24.     ! days to s
    tau_r     = tau_r*24.*3600.   ! days to s
    tau_f     = tau_f*24.*3600.   ! days to s
    tau_2     = tau_2*24.*3600.   ! days to s
    rd        = rd*1.e3      ! km to m
    width     = width*1.e3   ! km to m
    wlength   = wlength*1.e3 ! km to m
    sigma_i   = rd*sigma_i   ! 'rossby radii' to m
    sigma_jet = rd*sigma_jet ! 'rossby radii' to m

    !    ---- Calculate dependent variables ----
    dx = wlength/float(imax) ! (m) grid resolution in x
    dy = width/float(jmax)   ! (m) grid resolution in y
    el = pi/width            ! (1/m) converts radians to meridional wavenums
    rk = 2.*pi/wlength       ! (1/m) converts radians to zonal wavenums
    damp = visc*(dx**ndeg)/(dt*(pi**ndeg))  ! (m**ndeg/s) hyperviscocity

  end subroutine
end module
