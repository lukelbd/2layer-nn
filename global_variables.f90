module GLOBAL_VARIABLES

! ******* global variables for the baroclinic 2 layer model *****

  integer :: imx = 256, jmax = 256, mmax = 85, nmax = 170
  integer :: imax = 1+imx/2
  integer :: ndeg = 6  ! degree of hyperdiffusion

! ******* constants *******

  real, parameter :: pi=3.141592653589793

!    ---- run time parameters ---

  real :: dt=15.             !integration time increment(s)
  integer :: md = 1600        !dataio timestep
  integer :: mds = 0 ! start saving file
  integer :: mstart=1        !starting time step
  integer :: mend=160001        !terminating time step

! modify damping rate at certain timestep
  integer :: tchange = dt

!    ---- model dimensions ----
! (change model size ix, iz above)
  real :: width = 72000000.                  !(m) channel width in y
  real :: wlength = 16000000.                !(m) channel length in x
  real :: dx = wlength/float(imx) !(m) grid resolution in x
  real :: dy = width/float(jmax)  !(m) grid resolution in y

!    **** hyperviscosity ****
  real :: damp = 0.04*(dx**ndeg)/(dt*(pi**ndeg))  !(m**ndeg/s)
!    ---- basic flow parameters ----

  real :: beta=1.6e-11  !background vorticity gradient (1/(s*m))
  real :: delta = 4.4e-18  !linear decrease in beta (1/(s*m*m))
  real :: pi=3.141592653589793
  complex :: ur = (1.,0.), ui = (0.,1.), zero = (0.,0.)

!    ---- energy parameters ----
  real :: damp_coeff = 0.04
  real :: energy2,u0
  integer :: m

!    ---- initial amplitude ----
  real :: amp = 1.e-5  !(1/s)  initial amplitude in vorticity
! real :: amp = 1.e-3  !(1/s)  initial amplitude in vorticity
  real :: el = pi/width   !(1/m) initial meridional wavenumber
  real :: rk = 2.*pi/wlength  !(1/m) initial zonal wavenumber
  real :: rd = 800000.  !(m) internal rossby radius
  real :: tau_r = 30.*24.*3600.  !(s) radiative damping timescale
  real :: tau_f = 6.*24.*3600.  !(s) fricional damping timescale
  real :: tau_2 = 1000.*24.*3600.  !(s) additional layer2 damping timescale
  real :: u0 = 5. !(m/s) background (uniform) shear
  real :: sigma = rd*2.  !(m) width of the jet
  real :: del = 1.

!    ---- read namelist, potentially overwrite parameters ---
namelist /input_nml/ tchange, dt, mstart, mend
read(unit=2, nml=input_nml) ! simpler way
! open(unit=2, file='input.nml', form='formatted')

!    ---- arrays ----
!    wheeeee!
  complex :: vort_1(imax,jmax+1,4), &
         psi_1(imax,jmax+1), u_1(imax,jmax+1), v_1(imax,jmax+1), &
         adv_1(imax,jmax+1,3),visc_1(imax,jmax+1),q_1(imax,jmax+1), &
         qx_1(imax,jmax+1),qy_1(imax,jmax+1),rad_1(imax,jmax+1)
  complex :: vort_2(imax,jmax+1,4), &
         psi_2(imax,jmax+1), u_2(imax,jmax+1), v_2(imax,jmax+1), &
         adv_2(imax,jmax+1,3),visc_2(imax,jmax+1),q_2(imax,jmax+1), &
         qx_2(imax,jmax+1),qy_2(imax,jmax+1),rad_2(imax,jmax+1),&
         fric_2(imax,jmax+1)

  real :: ueq(jmax+1),uf(imx,jmax+1), &
          u_1_r(imax,jmax+1),u_1_i(imax,jmax+1),&
          v_1_r(imax,jmax+1),v_1_i(imax,jmax+1),&
          q_1_r(imax,jmax+1),q_1_i(imax,jmax+1),&
          qx_1_r(imax,jmax+1),qx_1_i(imax,jmax+1),&
          qy_1_r(imax,jmax+1),qy_1_i(imax,jmax+1),&
          p_1_r(imax,jmax+1),p_1_i(imax,jmax+1),&
          vqm_1(jmax+1,3), &
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
end module
