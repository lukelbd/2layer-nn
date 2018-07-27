!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! So far not much here, feel free to add to it
! Just deallocates memory locations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module cleanup
use netcdf
use global_variables
contains
subroutine clean
  implicit none
  !    ---- Clear Fourier descriptor ----
  ret = DftiFreeDescriptor(hcr)
  ret = DftiFreeDescriptor(hrc)

  !    ---- Close netcdf file ----
  ret = nf90_close(file_id)
end subroutine
end module
