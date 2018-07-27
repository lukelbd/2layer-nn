module global_variables
  implicit none
  !----------------------------------------------------------------------------!
  !    ---- Initial stuff/stuff that has to be hard-coded ----
  real, parameter :: pi=3.141592653589793
  complex :: one_r = (1.,0.), one_i = (0.,1.), zero = (0.,0.)
  integer :: day, t   ! time tracker
  real :: energy, cfl ! for monitoring integration
  ! Timing
  integer, parameter :: tstart=0
  ! Grid resolution (number of cells)
  ! Probably want approx 1:1 aspect ratio, then add one cell to y because
  ! top and bottom boundary aren't same point (whereas left/right boundary are)
  ! Without 1:1 aspect ratio PV injection is more complicated
  integer, parameter :: imax=256, jmax=256+1 ! we want a point exactly at y=0
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
  logical :: ll_seed_on                         ! flags
  integer :: dt, td                             ! time steps
  real :: tend, tchange, tds                    ! timing
  real :: width, wlength                        ! channel width, channel length
  real :: rd                                    ! radius of deformation
  real :: tau_r, tau_f, tau_2, tau_sp           ! damping stuff
  real :: y_sp                                  ! percentage of top half/bottom half that we want covered by sponge layer
  real :: shear, beta                           ! background state
  real :: amp_i, tau_i, sigma_i, wmin_i, wmax_i ! Noboru's stoshcastic forcing
  real :: ll_seed_amp                           ! lower level forcing
  real :: visc, ndeg                            ! hyperviscocity

  !    ---- Noboru's arrays ----
  ! Arrays storing Fourier coefficients from x-direction decomposition
  complex :: q1_sp(idft,jmax,2), adv1_sp(idft,jmax,3), &
         tmp1_sp(idft,jmax),  psi1_sp(idft,jmax), vor1_sp(idft,jmax), &
         visc1_sp(idft,jmax), rad1_sp(idft,jmax), sponge1_sp(idft,jmax), &
         force1_sp(idft,jmax), &
         u1_sp(idft,jmax),    v1_sp(idft,jmax),   &
         qx1_sp(idft,jmax),   qy1_sp(idft,jmax)
  complex :: q2_sp(idft,jmax,2), adv2_sp(idft,jmax,3), &
         tmp2_sp(idft,jmax),  psi2_sp(idft,jmax), vor2_sp(idft,jmax), &
         visc2_sp(idft,jmax), rad2_sp(idft,jmax), sponge2_sp(idft,jmax), &
         fric2_sp(idft,jmax), &
         u2_sp(idft,jmax),    v2_sp(idft,jmax),   &
         qx2_sp(idft,jmax),   qy2_sp(idft,jmax)
  ! Arrays storing values in 2D cartesian space
  real :: force1_cart(imax,jmax,2), adv1_cart(imax,jmax), &
          ufull1_cart(imax,jmax), vorfull1_cart(imax,jmax), qfull1_cart(imax,jmax), &
          u1_cart(imax,jmax),   v1_cart(imax,jmax),   q1_cart(imax,jmax), &
          psi1_cart(imax,jmax), vor1_cart(imax,jmax), &
          qx1_cart(imax,jmax),  qy1_cart(imax,jmax)
  real :: adv2_cart(imax,jmax), &
          ufull2_cart(imax,jmax), qfull2_cart(imax,jmax), vorfull2_cart(imax,jmax), &
          u2_cart(imax,jmax),   v2_cart(imax,jmax),   q2_cart(imax,jmax), &
          psi2_cart(imax,jmax), vor2_cart(imax,jmax), &
          qx2_cart(imax,jmax),  qy2_cart(imax,jmax)
  ! Means in tt space
  real :: tmp1_tt(jmax), tmp2_tt(jmax), &
          qyyflux1_tt(jmax,3), qyyflux2_tt(jmax,3), &
          qybar1_tt(jmax,2), qybar2_tt(jmax,2), &
          vorbar1_tt(jmax), qbar1_tt(jmax), ubar1_tt(jmax), &
          vorbar2_tt(jmax), qbar2_tt(jmax), ubar2_tt(jmax)
  ! Means in cartesian space
  real :: qybar1_cart(jmax),   qybar2_cart(jmax), &
          qflux1_cart(jmax),   qflux2_cart(jmax), &
          vorbar1_cart(jmax),  qbar1_cart(jmax),  ubar1_cart(jmax), &
          vorbar2_cart(jmax),  qbar2_cart(jmax),  ubar2_cart(jmax)
  ! Physical y coordinate, and masks for sponge and pv injection
  real :: y_cart(jmax), mask_i(jmax), mask_sp(jmax), mask_sp_tt(jmax)

  contains

  subroutine read_namelist
    implicit none
    !    ---- Declare namelist ----
    ! real :: tau_farray(jmax)
    namelist /input_nml/ &
      ll_seed_on, &
      width, wlength, &
      dt, td, &
      tend, tchange, tds, &
      tau_r, tau_f, tau_sp, tau_2, &
      shear, beta, rd, &
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
    tau_sp    = tau_sp*24.*3600.   ! days to s
    rd        = rd*1.e3      ! km to m
    width     = width*1.e3   ! km to m
    wlength   = wlength*1.e3 ! km to m
    sigma_i   = rd*sigma_i   ! 'rossby radii' to m

    !    ---- Calculate dependent variables ----
    dx = wlength/float(imax) ! (m) grid resolution in x
    dy = width/float(jmax-1) ! (m) grid resolution in y; non-cyclic so use -1, or something
    el = pi/width            ! (1/m) converts radians to meridional wavenums
    rk = 2.*pi/wlength       ! (1/m) converts radians to zonal wavenums
    damp = visc*(dx**ndeg)/(dt*(pi**ndeg))  ! (m**ndeg/s) hyperviscocity

  end subroutine
end module
