module GLOBAL_VARIABLES
  implicit none
  !----------------------------------------------------------------------------!
  !    ---- Initial stuff, and dimensions that must be hard-coded ----
  real, parameter :: pi=3.141592653589793
  complex :: ur = (1.,0.), ui = (0.,1.), zero = (0.,0.)
  integer :: t         ! time tracker
  real :: energy2, cfl ! for monitoring integration
  integer, parameter :: tstart=0          ! initial time
  integer, parameter :: imax=256, jmax=256 ! grid size
  integer, parameter :: mmax=85, nmax=170 ! 170+85 = 255, dunno what these are for
  integer, parameter :: idft=1+imax/2      ! signal is real, so need only N/2-1 complex coeffs
                                          ! plus the A_0 and A_(N/2) coeffs

  !    ---- To be calculated in subroutine below ----
  real :: damp, dx, dy, el, rk, bounds_l, bounds_k

  !    ---- To be supplied by namelist ----
  logical :: ll_seed_on, init_jet ! flags
  integer :: dt, td                                   ! time steps
  real :: tend, tchange, tds                          ! timing
  real :: width, wlength                              ! physical dimensions
  real :: rd, tau_r, tau_f, tau_sponge, tau_2                     ! damping stuff
  real :: bdel, beta, u0, jet_amp, jet_sigma          ! background state
  real :: y0, amp_k, sigma_y, sigma_l, sigma_k, tau_i ! stochastic forcing
  real :: famp, tau_fc                                ! Noboru's stoshcastic forcing
  real :: ll_seed_amp                                 ! lower level forcing
  real :: visc, ndeg                                  ! hyperviscocity

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
          qxy1(imax,jmax+1), qxx1(imax,jmax+1), qyy1(imax,jmax+1), &
          umean1(jmax+1), qymean1(jmax+1,4)
  real :: vqm_2(jmax+1,3), &
          uxy2(imax,jmax+1), vxy2(imax,jmax+1), pxy2(imax,jmax+1), &
          qxy2(imax,jmax+1), qxx2(imax,jmax+1), qyy2(imax,jmax+1), &
          umean2(jmax+1), qymean2(jmax+1,4)

  ! Other stuff that has to be double precision for some reason
  double precision :: qby_1(jmax+1),qbar1(jmax+1),ubar1(jmax+1)
  double precision :: qby_2(jmax+1),qbar2(jmax+1),ubar2(jmax+1)

  contains

  subroutine read_namelist
    implicit none
    !    ---- Declare namelist ----
    real :: tau_farray(jmax)
    namelist /input_nml/ &
      ll_seed_on, init_jet, &
      width, wlength, &
      dt, td, &
      tend, tchange, tds, &
      rd, tau_r, tau_f, tau_sponge, tau_2, &
      u0, jet_amp, jet_sigma, &
      famp, tau_fc, ll_seed_amp, &
      y0, amp_k, sigma_y, sigma_l, sigma_k, tau_i, &
      visc, ndeg, beta, bdel
    !    ---- Read namelist ----
    open(1, file='input.nml', status='old', action='read') ! status='old' requires that file exists, wut
    read(1, nml=input_nml)
      ! note 'input.nml' is a filename, input_nml is a declared variable
    close(1) ! close file

    !    ---- In-place unit scaling ----
    y0 = jmax/2 + y0*jmax ! make coordinates relative to center
    tend      = tend*3600.*24.
    tchange   = tchange*3600.*24.
    tds       = tds*3600.*24.
    tau_r     = tau_r*24.*3600. ! radiation
    tau_f     = tau_f*24.*3600. ! friction
    tau_i     = tau_i*24.*3600. ! injection
    tau_2     = tau_2*24.*3600.
    rd        = rd*1.e3
    width     = width*1.e3
    wlength   = wlength*1.e3
    jet_sigma = rd*jet_sigma ! (m) from rossby radii to sigma

    !    ---- Calc previously declared, empty variables ----
    dx = wlength/float(imax) ! (1000km) grid resolution in x
    dy = width/float(jmax)  ! (1000km) grid resolution in y
    el = pi/width           ! (1/m) meridional wavenum scaling
    rk = 2.*pi/wlength      ! (1/m) zonal wavenum scaling
    bounds_l = sigma_l*3.   ! bounds for sampling k, el wavenumbers; will be truncated
    bounds_k = sigma_k*3.
    damp = visc*(dx**ndeg)/(dt*(pi**ndeg))  ! (m**ndeg/s) hyperviscocity

  end subroutine
end module
