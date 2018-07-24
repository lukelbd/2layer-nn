!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Driver program for model, calls other submodules
! Nonlinear 2 layer baroclinic instability on the beta-plane.
! Periodic in x and rigid walls at north and south.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
  ! implicit none ! must be commented out, due to global variables error
  use mkl_dfti
  use global_variables
  use initial
  use forward
  use prognostics
  use diagnostics
  use io
  type(dfti_descriptor), pointer :: hx,hy
  integer :: l(1)
  l(1) = imax

  !    ---- Initialize fields ----
  call read_namelist ! (global_variables.f90)
  call init          ! (initial.f90)

  status = dfticreatedescriptor( hx, dfti_double, &
    dfti_real, 1, l)
  status = dftisetvalue(hx, dfti_placement, dfti_not_inplace)
  status = dfticommitdescriptor( hx )
  status = dfticreatedescriptor( hy, dfti_double, &
    dfti_real, 1, l)
  status = dftisetvalue(hy, dfti_placement, dfti_not_inplace)
  status = dfticommitdescriptor( hy )

  !    ---- Integration ----
  write(*,*) "Starting integration."
  do t = tstart, int(tend), int(dt) ! note that loop is normally end-inclusive
    !    ---- Invert spectral data, get diagnostic parameters, print diagnostics ----
    call diag(hx,hy) ! (diagnostics.f90) invert streamfunction from vorticity
    write(*,677) float(t)/(3600.*24.), energy2, cfl
    677 format("days = ", 1f8.3, " eke = ", 1p1e13.5, 1f3.3) ! 1p ensures non-zero digit to left of decimal
    !    ---- Save data ----
    if ((mod(t,td).eq.0).and.(t.ge.int(tds))) then 
      write(*,*) "Writing data."
      call dataio ! (io.f90) save data
    end if
    !    ---- Move timestep forward and evaluate prognostic data ----
    ! call energy1  ! (energy.f90) eddy energy calculation
    call prog    ! (prognostics.f90) eddy prognostic equation
    call iterate ! (forward.f90) time forwarding
    ! Testing things
    ! print *, y0, tend, tchange, tds, tau_r, tau_f, tau_2, rd, width, wlength, sigma_jet
    stop
    ! Test
  enddo
  write(*,*) 'Final Day = ', float(t)/(3600.*24.)
  close(11)

  status = dftifreedescriptor( hx )
  status = dftifreedescriptor( hy )

  stop
end program
