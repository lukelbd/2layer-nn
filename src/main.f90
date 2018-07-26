!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Driver program for model, calls other submodules
! Nonlinear 2 layer baroclinic instability on the beta-plane.
! Periodic in x and rigid walls at north and south.
! Uses MKL Intel libraries
! For more information, see: https://scc.ustc.edu.cn/zlsc/sugon/intel/mkl/mkl_manual
! For some reason trig transform info is under "Partial Differential Equations support"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
  ! implicit none ! must be commented out, due to global variables error
  use global_variables
  use initial
  use forward
  use prognostics
  use diagnostics
  use mkl_dfti            ! includes some boilerplate stuff, and Fourier transforms
  use mkl_trig_transforms ! includes trig transforms
  use io
  type(dfti_descriptor), pointer :: hcr, hrc ! handles for real-to-complex and complex-to-real fourier transforms
  integer :: ret, shp(1)
  shp(1) = imax ! shape of arrays sent to Fourier methods

  !    ---- Initialize fields ----
  call read_namelist ! (global_variables.f90)
  call init          ! (initial.f90)

  !    ---- Initialize handles for Fourier calls ----
  ! Arg 1 is handle, Arg 2 is precision of transform, Arg 3 is
  ! 'forward domain' (i.e. the 'forward' transform takes place on only
  ! real-valued sequences), Arg 4/5 are the number of dims/dimension lengths
  ret = dfticreatedescriptor(hcr, dfti_double, dfti_real, 1, shp)
  ret = dftisetvalue(hcr, dfti_placement, dfti_not_inplace)
  ret = dfticommitdescriptor(hcr)
  ret = dfticreatedescriptor(hrc, dfti_double, dfti_real, 1, shp)
  ret = dftisetvalue(hrc, dfti_placement, dfti_not_inplace)
  ret = dfticommitdescriptor(hrc)

  !    ---- Integration ----
  write(*,*) "Starting integration."
  do t = tstart, int(tend), int(dt) ! note that loop is normally end-inclusive
    !    ---- Invert spectral data, get diagnostic params ----
    call diag(hcr,hrc) ! (diagnostics.f90) pass existing handles 
    !    ---- Save data ----
    if ((mod(t,td).eq.0).and.(t.ge.int(tds))) then 
      write(*,*) "Writing data."
      call dataio ! (io.f90)
    end if
    !    ---- Move timestep forward and evaluate prognostic data ----
    ! call energy1  ! (energy.f90) eddy energy calculation
    call prog    ! (prognostics.f90) eddy prognostic equation
    call iterate ! (forward.f90) time forwarding
    ! Testing things
    ! print *, y0, tend, tchange, tds, tau_r, tau_f, tau_2, rd, width, wlength, sigma_jet
    stop
  enddo
  close(11)
  ret = dftifreedescriptor(hcr)
  ret = dftifreedescriptor(hrc)
  stop
end program
