!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Driver program for model, calls other submodules and prints information
! Nonlinear 2 layer baroclinic instability on the beta-plane.
! Periodic in x and rigid walls at north and south.
! Uses MKL Intel libraries
! For more information, see: https://scc.ustc.edu.cn/zlsc/sugon/intel/mkl/mkl_manual
! For some reason trig transform info is under "Partial Differential Equations support"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
  use global_variables
  use initial
  use cleanup
  use netcdf
  use forward
  use prognostics
  use diagnostics
  use io
  ! implicit none ! cannot use this, due to global variables error

  !    ---- Initialize fields and other stuff ----
  call init ! (initial.f90)

  !    ---- Printout ----
  ! write(*,400,advance='no') mask_sp
  write(*,*) "sponge mask = "
  do j=1,jmax
    if (mod(j,10).eq.0) then
      write(*,'(f6.3)',advance='yes') mask_sp(j)
    else
      write(*,'(f6.3)',advance='no') mask_sp(j)
    end if
  enddo
  write(*,*)
  write(*,*) "inject mask = "
  do j=1,jmax
    if (mod(j,10).eq.0) then
      write(*,'(f6.3)',advance='yes') mask_i(j)
    else
      write(*,'(f6.3)',advance='no') mask_i(j)
    end if
  enddo
  write(*,*)
  write(*,500) shear, shear_phillips
  500 format("background shear = ", f4.1, " unstable shear = ", f4.1)
  write(*,200) tau_f, tau_r, tau_sp, tau_i
  200 format("friction = ", f12.1, " radiation = ", f12.1, " sponge = ", f12.1, " forcing = ", f12.1)
  write(*,100) dt, dt_io, t_end
  100 format("time step = ", i8, " io step = ", i8, " final time = ", i12)

  !    ---- Integration ----
  write(*,*) "Starting integration."
  do t = t_start, t_end, dt ! note that loop is normally end-inclusive
    !    ---- Day ----
    day = float(t)/(3600.*24.)

    !    ---- Invert spectral data, get diagnostic params ----
    call diag ! (diagnostics.f90) get diagnostic stuff

    !    ---- Move timestep forward and evaluate prognostic data ----
    call prog    ! (prognostics.f90) eddy prognostic equation
    call iterate ! (forward.f90) time forwarding

    !    ---- Save data ----
    ! If we are past spinup 'days_spinup', and we are on the data save interval, save
    if (t.ge.t_spinup .and. mod(t,dt_io).eq.0) then 
      write(*,400) t_io
      400 format(" Writing data on index:", i4)
      call ncwrite ! (io.f90) save data; if this is first time, will initialize handles
    endif

    !    ---- Print diagnostic output ----
    if (mod(t,dt*1000).eq.1) then
      write(*,*) 'Printing upper level zonal means.'
      do j = 1,jmax
        write(*,*) j, 'ubar1 = ', ubar1_cart(j)+shear, 'qbar1 = ', qbar1_cart(j)
      enddo
    endif
    write(*,300) day, energy(1), umax(1), cfl(1)
    300 format("days = ", f6.3, " eke = ", e12.5, " umax = ", f6.3, " cfl = ", f5.3)
    ! 1p ensures non-zero digit to left of decimal
  enddo

  !    ---- Cleanup ----
  call clean ! (cleanup.f90)
end program
