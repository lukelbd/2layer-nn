!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! File for 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module io
use netcdf
use global_variables
contains
subroutine ncwrite
  implicit none
  integer :: s4d(4), s3d(3), s1d(1)
  integer :: x_id, y_id, z_id, t_id
  integer :: q_id, u_id, v_id, f_id, vor_id, psi_id, eke_id
  integer :: xdim_id, ydim_id, zdim_id, tdim_id
  !    ---- Create netcdf file and initialze variables ----
  ! Probably makes sense to put this here instead of in initialize.f90, because
  ! perhaps we want to run the model without saving any data -- just want to read
  ! some print statements/play with it.
  if (nctime.eq.0) then
    ! Create file
    ret = nf90_create(path='data.nc', cmode=nf90_clobber, ncid=file_id)
    ! Define dimensions
    ret = nf90_def_dim(file_id, 'x', imax, xdim_id)
    ret = nf90_def_dim(file_id, 'y', jmax, ydim_id)
    ret = nf90_def_dim(file_id, 'z', jmax, zdim_id)
    ret = nf90_def_dim(file_id, 'time', nf90_unlimited, tdim_id)
    ! Dimension variables
    ret = nf90_def_var(file_id, 'x',  nf90_float, (/xdim_id/), x_id)
    ret = nf90_def_var(file_id, 'y',  nf90_float, (/ydim_id/), y_id)
    ret = nf90_def_var(file_id, 'z',  nf90_int,  (/zdim_id/), z_id)
    ret = nf90_def_var(file_id, 't',  nf90_float, (/tdim_id/), t_id)
    ret = nf90_put_var(file_id, x_id, x_cart, start=(/1/), count=(/imax/))
    ret = nf90_put_var(file_id, y_id, y_cart, start=(/1/), count=(/jmax/))
    ret = nf90_put_var(file_id, z_id, (/1,2/), start=(/1/), count=(/2/))
    ! Data variables
    s4d = (/xdim_id, ydim_id, zdim_id, tdim_id/) ! 2 levels
    ret = nf90_def_var(file_id, 'q',   nf90_float, s4d, q_id)
    ret = nf90_def_var(file_id, 'u',   nf90_float, s4d, u_id)
    ret = nf90_def_var(file_id, 'v',   nf90_float, s4d, v_id)
    ret = nf90_def_var(file_id, 'vor', nf90_float, s4d, vor_id)
    s3d = (/xdim_id, ydim_id, tdim_id/) ! 1 level
    ret = nf90_def_var(file_id, 'f',   nf90_float, s3d, f_id)   ! forcing perturbations
    s1d = (/tdim_id/)
    ret = nf90_def_var(file_id, 'eke', nf90_float, s1d, eke_id)
    ! Global attributes
    ret = nf90_put_att(file_id, nf90_global, 'title', '2-layer QG model results')
    ! Dimension attributes
    ret = nf90_put_att(file_id, x_id, 'units', 'm')
    ret = nf90_put_att(file_id, y_id, 'units', 'm')
    ret = nf90_put_att(file_id, t_id, 'units', 'days')
    ret = nf90_put_att(file_id, z_id, 'units', 'z(1)=top z(2)=bottom')
    ! Variable attributes
    ret = nf90_put_att(file_id, q_id,   'units',     '1/s')
    ret = nf90_put_att(file_id, u_id,   'units',     'm/s')
    ret = nf90_put_att(file_id, v_id,   'units',     'm/s')
    ret = nf90_put_att(file_id, vor_id, 'units',     '1/s')
    ret = nf90_put_att(file_id, f_id,   'units',     '1/s^2')
    ret = nf90_put_att(file_id, eke_id, 'units',     'J/m^2')
    ret = nf90_put_att(file_id, q_id,   'long_name', 'total potential vorticity')
    ret = nf90_put_att(file_id, u_id,   'long_name', 'total zonal wind')
    ret = nf90_put_att(file_id, v_id,   'long_name', 'total meridional wind')
    ret = nf90_put_att(file_id, vor_id, 'long_name', 'total relative vorticity')
    ret = nf90_put_att(file_id, f_id,   'long_name', 'total potential vorticity forcing in upper layer')
    ret = nf90_put_att(file_id, eke_id, 'long_name', 'total eddy kinetic energy')
    ! Finished defining
    ret = nf90_enddef(file_id)
  endif

  !    ---- Save data ----
  ! Add current day to time variable
  ret = nf90_put_var(file_id, t_id, (/day/), start=(/nctime/), count=(/1/))
  ! If we are past spinup 'tds', and we are on the data save interval, save
  ! Note if you mess up the dimensionality of input arrays, error message will be mysterious:
  ! "There is no matching specific function for this generic function reference."
  ret = nf90_put_var(file_id, q_id, qfull1_cart, start=(/1,1,1,nctime/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, q_id, qfull2_cart, start=(/1,1,2,nctime/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, u_id, ufull1_cart, start=(/1,1,1,nctime/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, u_id, ufull2_cart, start=(/1,1,2,nctime/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, v_id, v1_cart, start=(/1,1,1,nctime/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, v_id, v2_cart, start=(/1,1,2,nctime/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, vor_id, vorfull1_cart, start=(/1,1,1,nctime/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, vor_id, vorfull2_cart, start=(/1,1,2,nctime/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, f_id, force1_cart(:,:,1), start=(/1,1,nctime/), count=(/imax,jmax,1/))
  ret = nf90_put_var(file_id, eke_id, energy, start=(/nctime/), count=(/1/))

  !    ---- Increment save time ----
  nctime = nctime+1

end subroutine
end module
