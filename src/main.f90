!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Driver program for model, calls other submodules and prints information
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
  use cleanup
  use netcdf
  use forward
  use prognostics
  use diagnostics
  use io

  !    ---- Initialize fields and other stuff ----
  call init ! (initial.f90)
  ! Printout; should definitely add this
  ! print *, y0, tend, tds, tau_r, tau_f, rd, width, wlength, sigma_jet

  !    ---- Integration ----
  write(*,*) "Starting integration."
  do t = tstart, int(tend), int(dt) ! note that loop is normally end-inclusive
    !    ---- Day ----
    day = float(t)/(3600.*24.)

    !    ---- Invert spectral data, get diagnostic params ----
    call diag ! (diagnostics.f90) get diagnostic stuff

    !    ---- Move timestep forward and evaluate prognostic data ----
    call prog    ! (prognostics.f90) eddy prognostic equation
    call iterate ! (forward.f90) time forwarding

    !    ---- Save data ----
    if (t.ge.int(tds) .and. mod(t,td).eq.0) then 
      write(*,*) "Writing data."
      call ncwrite ! (io.f90) save data; if this is first time, will initialize handles
    endif

    !    ---- Print diagnostic output ----
    if(mod(int(t),int(dt)*1000).eq.1) then
      write(*,*) 'Printing zonal mean diagnostics.'
      do j = 1,jmax
        write(*,*) j, 'ubar1 = ', ubar1_cart(j)+shear, 'qbar1 = ', qbar1_cart(j)
      enddo
    endif
    write(*,677) day, energy(1), umax(1), cfl(1)
    677 format("days = ", 1f8.3, " eke = ", 1p1e13.5, " umax = ", 1f3.3, " cfl = ", 1f3.3)
    ! 1p ensures non-zero digit to left of decimal
  enddo

  !    ---- Cleanup ----
  call clean ! (cleanup.f90)
end program
