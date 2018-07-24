module global_variables
  implicit none
  !----------------------------------------------------------------------------!
  !    ---- Initial stuff, and dimensions that must be hard-coded ----
  real, parameter :: pi=3.141592653589793
  complex :: one_r = (1.,0.), one_i = (0.,1.), zero = (0.,0.)
  integer :: day, t   ! time tracker
  real :: energy, cfl ! for monitoring integration
  ! Timing
  integer, parameter :: tstart=0
  ! Grid resolution (number of cells)
  integer, parameter :: imax=512, jmax=256+1 ! we want a point exactly at y=0
  ! integer, parameter :: imax=256, jmax=256
  ! Number complex coefficients for cyclic DFT in x-direction
  ! Signal is real, so need only N/2-1 complex coeffs plus the A_0 and A_(N/2) coeffs
  ! For y sine transform, coefficients are all real; need all <jmax> of them plus the constant offset
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
  real :: shear, beta                        ! background state
  real :: amp_jet, sigma_jet              ! initial jet
  real :: amp_i, tau_i, sigma_i, wmin_i, wmax_i ! Noboru's stoshcastic forcing
  real :: ll_seed_amp                     ! lower level forcing
  real :: visc, ndeg                      ! hyperviscocity

  !    ---- Noboru's arrays ----
  ! Arrays storing Fourier coefficients from x-direction decomposition
  complex :: q1_c(idft,jmax,4), psi1_c(idft,jmax), vor1_c(idft,jmax), &
         adv1_c(idft,jmax,3), u1_c(idft,jmax),   v1_c(idft,jmax),     &
         visc1_c(idft,jmax),  rad1_c(idft,jmax), force1_c(idft,jmax), &
         qx1_c(idft,jmax),    qy1_c(idft,jmax)
  complex :: q2_c(idft,jmax,4), psi2_c(idft,jmax), vor2_c(idft,jmax), &
         adv2_c(idft,jmax,3), u2_c(idft,jmax),   v2_c(idft,jmax),     &
         visc2_c(idft,jmax),  rad2_c(idft,jmax), fric_2(idft,jmax),   &
         qx2_c(idft,jmax),    qy2_c(idft,jmax)
  ! Special array
  real :: ymask(jmax)
  ! Arrays stored in real space, suitable for output
  real :: f1_out(imax,jmax,2), utot1_out(imax,jmax), &
          u1_out(imax,jmax), v1_out(imax,jmax), q1_out(imax,jmax), &
          psi1_out(imax,jmax), vor1_out(imax,jmax)
  real :: utot2_out(imax,jmax), &
          u2_out(imax,jmax), v2_out(imax,jmax), q2_out(imax,jmax), &
          psi2_out(imax,jmax), vor2_out(imax,jmax)
  ! Special; stores mean and mean gradient
  real :: ubar1_out(jmax), qyyflux1_out(jmax,3), qybar1_out(jmax,4), &
          ubar2_out(jmax), qyyflux2_out(jmax,3), qybar2_out(jmax,4)

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
      shear, beta, rd, &
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
    dy = width/float(jmax-1) ! (m) grid resolution in y; non-cyclic so use -1, or something
    el = pi/width            ! (1/m) converts radians to meridional wavenums
    rk = 2.*pi/wlength       ! (1/m) converts radians to zonal wavenums
    damp = visc*(dx**ndeg)/(dt*(pi**ndeg))  ! (m**ndeg/s) hyperviscocity

  end subroutine
end module
