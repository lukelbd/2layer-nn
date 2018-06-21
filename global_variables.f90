module GLOBAL_VARIABLES
  implicit none
  ! public :: init_variables
  !----------------------------------------------------------------------------!
  !    ---- Initial stuff, and dimensions that must be hard-coded ----
  real, parameter :: pi=3.141592653589793                      ! pi
  integer :: t                                                 ! time tracker
  complex :: ur = (1.,0.), ui = (0.,1.), zero = (0.,0.)        ! initial flow
  integer, parameter :: tstart=0                               ! initial time
  integer, parameter :: imx=256, jmax=256, mmax=85, nmax=170  ! dimensions
  integer, parameter :: imax=1+imx/2                           ! myserious
  real :: energy2, cfl ! for monitoring integration

  !    ---- To be calculated in subroutine below ----
  real :: damp, dx, dy, el, rk, bounds_l, bounds_k

  !    ---- To be supplied by namelist ----
  logical :: ll_seed_on, init_jet ! flags
  integer :: dt, td                                   ! time steps
  real :: tend, tchange, tds                          ! timing
  real :: width, wlength                              ! physical dimensions
  real :: rd, tau_r, tau_f, tau_2                     ! damping stuff
  real :: bdel, beta, u0, jet_amp, jet_sigma          ! background state
  real :: y0, amp_k, sigma_y, sigma_l, sigma_k, tau_i ! stochastic forcing
  real :: famp, tau_fc                                ! Noboru's stoshcastic forcing
  real :: ll_seed_amp                                 ! lower level forcing
  real :: visc, ndeg                                  ! hyperviscocity

  !    ---- Noboru's arrays ----
  complex :: vort_1(imax,jmax+1,4),force_1(imax,jmax+1), &
         psi_1(imax,jmax+1), u_1(imax,jmax+1), v_1(imax,jmax+1), &
         adv_1(imax,jmax+1,3),visc_1(imax,jmax+1),q_1(imax,jmax+1), &
         qx_1(imax,jmax+1),qy_1(imax,jmax+1),rad_1(imax,jmax+1)
  complex :: vort_2(imax,jmax+1,4), &
         psi_2(imax,jmax+1), u_2(imax,jmax+1), v_2(imax,jmax+1), &
         adv_2(imax,jmax+1,3),visc_2(imax,jmax+1),q_2(imax,jmax+1), &
         qx_2(imax,jmax+1),qy_2(imax,jmax+1),rad_2(imax,jmax+1),&
         fric_2(imax,jmax+1)

  real :: ueq(jmax+1),uf(imx,jmax+1),ymask(jmax+1), &
          u_1_r(imax,jmax+1),u_1_i(imax,jmax+1),&
          v_1_r(imax,jmax+1),v_1_i(imax,jmax+1),&
          q_1_r(imax,jmax+1),q_1_i(imax,jmax+1),&
          f_1_r(imax,jmax+1),f_1_i(imax,jmax+1),&
          qx_1_r(imax,jmax+1),qx_1_i(imax,jmax+1),&
          qy_1_r(imax,jmax+1),qy_1_i(imax,jmax+1),&
          p_1_r(imax,jmax+1),p_1_i(imax,jmax+1),&
          vqm_1(jmax+1,3), fxy1(imx,jmax+1,2),&
          qxy1(imx,jmax+1),vxy1(imx,jmax+1),&
          pxy1(imx,jmax+1),uxy1(imx,jmax+1),umean1(jmax+1),&
          qxx1(imx,jmax+1),qyy1(imx,jmax+1),qymean1(jmax+1,4)
  double precision :: qby_1(jmax+1),qbar1(jmax+1),ubar1(jmax+1)

  real ::  &
          u_2_r(imax,jmax+1),u_2_i(imax,jmax+1),&
          v_2_r(imax,jmax+1),v_2_i(imax,jmax+1),&
          q_2_r(imax,jmax+1),q_2_i(imax,jmax+1),&
          qx_2_r(imax,jmax+1),qx_2_i(imax,jmax+1),&
          qy_2_r(imax,jmax+1),qy_2_i(imax,jmax+1),&
          p_2_r(imax,jmax+1),p_2_i(imax,jmax+1),&
          vqm_2(jmax+1,3),&
          qxy2(imx,jmax+1),vxy2(imx,jmax+1),&
          qxx2(imx,jmax+1),qyy2(imx,jmax+1),umean2(jmax+1),&
          pxy2(imx,jmax+1),uxy2(imx,jmax+1),qymean2(jmax+1,4)
  double precision :: qby_2(jmax+1),qbar2(jmax+1),ubar2(jmax+1)

  contains

  subroutine init_variables
    implicit none
    !    ---- Declare namelist ----
    namelist /input_nml/ &
      ll_seed_on, init_jet, &
      width, wlength, &
      dt, td, &
      tend, tchange, tds, &
      rd, tau_r, tau_f, tau_2, &
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
    dx = wlength/float(imx) ! (1000km) grid resolution in x
    dy = width/float(jmax)  ! (1000km) grid resolution in y
    el = pi/width           ! (1/m) meridional wavenum scaling
    rk = 2.*pi/wlength      ! (1/m) zonal wavenum scaling
    bounds_l = sigma_l*3.   ! bounds for sampling k, el wavenumbers; will be truncated
    bounds_k = sigma_k*3.
    damp = visc*(dx**ndeg)/(dt*(pi**ndeg))  ! (m**ndeg/s) hyperviscocity

  end subroutine
end module
