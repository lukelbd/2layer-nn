!****** Nonlinear 2 layer baroclinic instability on the beta-plane.  Periodic in x and rigid walls at north and south *****   

program main
  ! implicit none
  use MKL_DFTI
  use GLOBAL_VARIABLES
  use INITIAL
  use FORWARD
  use EXEC
  use IO
  use INVERT
  type(DFTI_DESCRIPTOR), POINTER :: HX,HY
  integer :: L(1)
  L(1) = imx

  ! ****** Initialize fields ****  
  call read_namelist ! (global_variables.f90)
  call init_data     ! (initial.f90)

  Status = DftiCreateDescriptor( HX, DFTI_DOUBLE, &
    DFTI_REAL, 1, L)
  Status = DftiSetValue(HX, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
  Status = DftiCommitDescriptor( HX )
  Status = DftiCreateDescriptor( HY, DFTI_DOUBLE, &
    DFTI_REAL, 1, L)
  Status = DftiSetValue(HY, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
  Status = DftiCommitDescriptor( HY )

  ! ****** Integration *****
  write(6,*) "Starting integration."
  do t = tstart, int(tend), int(dt) ! note that loop is normally end-inclusive
  ! Invert data from previous timestep, and print diagnostics
  call invert1(HX,HY)                             ! (invert.f90) invert streamfunction from vorticity
  write(6,677) float(t)/(3600.*24.), energy2, cfl
  677 format("days = ", 1f8.3, " eke = ", 1p1e13.5, 1f3.3) ! 1p ensures non-zero digit to left of decimal
  ! Input data
  if ((mod(t,td).eq.0).and.(t.ge.int(tds))) then 
    write(6,*) "Writing data."
    call dataio ! (io.f90) save data
  end if
  ! Iterate and evaluate prognostic data
  ! call energy1                              ! (energy.f90) eddy energy calculation
  call prog                                       ! (exec.f90) eddy prognostic equation
  call forward1                                   ! (forward.f90) time forwarding
  ! Test
  print *, y0, tend, tchange, tds, tau_r, tau_f, tau_i, tau_2, rd, width, wlength, jet_sigma
  stop
  ! Test
  enddo
  write(6,*) 'Final Day = ', float(t)/(3600.*24.)
  close(11)

  Status = DftiFreeDescriptor( HX )
  Status = DftiFreeDescriptor( HY )

  stop
end program
