!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Global variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module global_variables
  use mkl_dfti
  implicit none
  complex :: one_r=(1.,0.), one_i=(0.,1.), zero=(0.,0.) ! helpful constants
  real, parameter :: pi=3.141592653589793 ! pi
  !----------------------------------------------------------------------------!
  !    ---- Important spacing params that have to be hard-coded ----
  ! Grid resolution (number of cells)
  ! Probably want approx 1:1 aspect ratio, then add one cell to y because
  ! top and bottom boundary aren't same point (whereas left/right boundary are)
  ! Without 1:1 aspect ratio PV injection is more complicated
  integer, parameter :: imax=512, jmax=256+1 ! we want a point exactly at y=0
  ! Number complex coefficients for cyclic DFT in x-direction
  ! Signal is real, so need only N/2-1 complex coeffs plus the A_0 and A_(N/2) coeffs
  ! For y sine transform, coefficients are all real; need all <jmax> of them plus the constant offset
  ! See: https://en.wikipedia.org/wiki/Discrete_Fourier_transform#The_real-input_DFT
  integer, parameter :: idft=1+imax/2
  ! Truncation of y-sine ft and x-dft to resolve (avoid aliasing) quadratic terms
  ! See: 10.1098/rsta.2013.0289
  ! Requirement is: imax >= (3*xtrunc+1), jmax >= (3*ytrunc+1)/2
  integer, parameter :: itrunc=170, jtrunc=170 ! for 512 by 256
  ! integer, parameter :: itrunc=85, jtrunc=170 ! for 256 by 256

  !    ---- Namelist accessible params ----
  ! Temporal
  integer :: dt=300, dt_io=21600 ! (s) integration time step and data-save timestep steps
  real :: days=100.0      ! (days) total integration time
  real :: days_spinup=0.0 ! (days) 'spinup' time before which no data is saved; useful to save disk space/ignore useless data
  ! Spatial
  real :: width=72.e3, wlength=36.e3 ! (km) channel width, channel length
  real :: rd=800.0                   ! (km) radius of deformation
  ! Background
  real :: shear=3.0    ! (m/s) basic state shear
  real :: beta=1.6e-11 ! (1/(m*s)) background state; default is value at 45deg latitude
  ! Basic damping
  real :: tau_r=30.0, tau_f=6.0, tau_sp=1.0 ! (days) radiation, friction, sponge (at the edge) damping timescales
  real :: y_sp=0.3                          ! (unitless) percentage of top half/bottom half of channel that we want covered by sponge layer
  integer :: ndeg=6                         ! (unitless) degree of hyperdifussion (must be *even*)
  real :: visc = 0.04                       ! (m^ndeg/s) viscocity coefficient (probably best not to touch this one)
  ! PV injection
  real :: amp_i=3.0e-8            ! (1/s^2) amplitude of dq/dt injections
  real :: tau_i=600.0             ! (s) forcing correlation timescale
  real :: y_i=0.2                 ! (unitless) e-folding width of injection band, as percentage of top/bottom half of channel
  integer :: wmin_i=41, wmax_i=46 ! (unitless) min and max injection wavenumbers
  ! Initial low-level forcing
  logical :: ll_seed_on=.true. ! (logical) flags
  real :: ll_seed_amp=1.0e-7   ! (s^-1) amplitude of initial injections
  ! Stochastic forcing parameters Sam, Momme, and Luke came up with
  ! Consider implementing in future; some of the above
  ! real :: amp_i = 1.e-5, ! (1/s) amplitude of pv injections
  ! real :: tau_i = 1.,    ! (days) injection timecale, CDF of binomial distribution characterizing the <yes/no> trials for whether we inject
  ! integer :: center_l = 50, ! (unitless) central injection wavenumber
  ! integer :: center_k = 50,
  ! real :: sigma_l = 5.,  ! (unitless) amplitude decay away from central wavenumber
  ! real :: sigma_k = 5.,

  !    ---- Other params ----
  ! Some derived from above namelist parameters
  integer :: t, t_start=0, t_end, t_io=1, t_spinup
  real :: day, damp, dx, dy, el, rk

  !    ---- Scalars ----
  real :: energy(1), umax(1), cfl(1)

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
  real :: x_cart(imax), y_cart(jmax), mask_i(jmax), mask_sp(jmax), mask_sp_tt(jmax)
  real :: rkks(idft), ells(jmax) ! will populate these with rkk and ell

  !    ---- Misc ----
  ! Note netcdf handles ***have*** to be here or they will get overwritten every time!
  type(dfti_descriptor), pointer :: hcr, hrc ! handles for real-to-complex and complex-to-real fourier transforms
  integer :: file_id ! netcdf file descriptor
  integer :: ret ! dummy variable for storing return value
  integer :: s4d(4), s3d(3), s1d(1)
  integer :: x_id, y_id, z_id, t_id ! dimension variable descriptors
  integer :: q_id, u_id, v_id, f_id, vor_id, psi_id
  integer :: qx_id, qy_id, adv_id, eke_id, cfl_id ! data variable descriptors
  integer :: xdim_id, ydim_id, zdim_id, tdim_id ! dimension descriptors

end module
