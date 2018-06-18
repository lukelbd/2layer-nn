module GLOBAL_VARIABLES

! ******* global variables for the baroclinic 2 layer model *****

  integer, parameter :: imx = 256, jmax = 256, mmax = 85, nmax = 170
  integer, parameter :: imax = 1+imx/2
  integer, parameter :: ndeg = 6  ! degree of hyperdiffusion

! ******* constants *******

!    ---- run time parameters ---

  real, parameter :: dt=15.             !integration time inclement(s)
  integer, parameter :: md = 1600        !dataio timestep
  integer, parameter :: mstart=1        !starting time step
  integer, parameter :: mend=160001        !terminating time step
 
!    ---- model dimensions ----
!        (change model size ix, iz above)

  real, parameter :: width = 72000000.          !(m)y
  real, parameter :: wlength = 16000000.          !(m)x
!    ---- basic flow parameters ----

  real, parameter :: beta=1.6e-11  !background vorticity gradient (1/(s*m))
  real, parameter :: delta = 4.4e-18  !linear devrease in beta (1/(s*m*m))
  real, parameter :: pi=3.141592653589793
  complex, parameter :: ur = (1.,0.), ui = (0.,1.), zero = (0.,0.)

!    ---- initial amplitude ----

  real, parameter :: amp = 1.e-5  !(1/s)  initial amplitude in vorticity
! real, parameter :: amp = 1.e-3  !(1/s)  initial amplitude in vorticity
  real, parameter :: el = pi/width   !(1/m) initial meridional wavenumber
  real, parameter :: rk = 2.*pi/wlength  !(1/m) initial zonal wavenumber
  real, parameter :: rd = 800000.  !(m) internal rossby radius
  real, parameter :: sigma = rd*2.  !(m) width of the jet
  real, parameter :: del = 1.

!    ---- arrays ----

  complex :: vort_1(imax,jmax+1,4), &
         psi_1(imax,jmax+1), u_1(imax,jmax+1), v_1(imax,jmax+1), &
         adv_1(imax,jmax+1,3),visc_1(imax,jmax+1),q_1(imax,jmax+1), &
         qx_1(imax,jmax+1),qy_1(imax,jmax+1) 
  complex :: vort_2(imax,jmax+1,4), &
         psi_2(imax,jmax+1), u_2(imax,jmax+1), v_2(imax,jmax+1), &
         adv_2(imax,jmax+1,3),visc_2(imax,jmax+1),q_2(imax,jmax+1), &
         qx_2(imax,jmax+1),qy_2(imax,jmax+1) 
            
  real :: ueq(jmax+1),uf(imx,jmax+1), & 
          u_1_r(imax,jmax+1),u_1_i(imax,jmax+1),&
          v_1_r(imax,jmax+1),v_1_i(imax,jmax+1),&
          q_1_r(imax,jmax+1),q_1_i(imax,jmax+1),&
          qx_1_r(imax,jmax+1),qx_1_i(imax,jmax+1),&
          qy_1_r(imax,jmax+1),qy_1_i(imax,jmax+1),&
          p_1_r(imax,jmax+1),p_1_i(imax,jmax+1),&
          vqm_1(jmax+1,3), &
          qxy1(imx,jmax+1),vxy1(imx,jmax+1),&
          pxy1(imx,jmax+1),uxy1(imx,jmax+1),&
          qxx1(imx,jmax+1),qyy1(imx,jmax+1),qymean1(jmax+1,4)
  double precision :: qby_1(jmax+1),qbar1(jmax+1)

  real ::  & 
          u_2_r(imax,jmax+1),u_2_i(imax,jmax+1),&
          v_2_r(imax,jmax+1),v_2_i(imax,jmax+1),&
          q_2_r(imax,jmax+1),q_2_i(imax,jmax+1),&
          qx_2_r(imax,jmax+1),qx_2_i(imax,jmax+1),&
          qy_2_r(imax,jmax+1),qy_2_i(imax,jmax+1),&
          p_2_r(imax,jmax+1),p_2_i(imax,jmax+1),&
          vqm_2(jmax+1,3),&
          qxy2(imx,jmax+1),vxy2(imx,jmax+1),&
          qxx2(imx,jmax+1),qyy2(imx,jmax+1),&
          pxy2(imx,jmax+1),uxy2(imx,jmax+1),qymean2(jmax+1,4)
  double precision :: qby_2(jmax+1),qbar2(jmax+1)

!    ---- energy parameters ----

  real :: dx,dy,energy2,damp,u0
  integer :: m

end module
