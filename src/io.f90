!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! File for 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module io
use netcdf
use global_variables
contains
subroutine ncwrite
  implicit none
  !    ---- Create netcdf file and initialze variables ----
  ! Probably makes sense to put this here instead of in initialize.f90, because
  ! perhaps we want to run the model without saving any data -- just want to read
  ! some print statements/play with it.
  if (t_io.eq.1) then
    ! Create file (also 'opens' it and turns on 'defmode'; no need to call those functions)
    ! The clobber mode says to overwrite any existing dataset (this is the default)
    write(*,*) 'Creating dataset.'
    ret = nf90_create(path='data.nc', cmode=nf90_clobber, ncid=file_id)
    ! Define dimensions
    ret = nf90_def_dim(file_id, 'x', imax, xdim_id)
    ret = nf90_def_dim(file_id, 'y', jmax, ydim_id)
    ret = nf90_def_dim(file_id, 'z', 2, zdim_id)
    ret = nf90_def_dim(file_id, 'time', nf90_unlimited, tdim_id)
    ! Dimension variables
    ret = nf90_def_var(file_id, 'x',  nf90_float, (/xdim_id/), x_id)
    ret = nf90_def_var(file_id, 'y',  nf90_float, (/ydim_id/), y_id)
    ret = nf90_def_var(file_id, 'z',  nf90_int,  (/zdim_id/), z_id)
    ret = nf90_def_var(file_id, 'time',  nf90_float, (/tdim_id/), t_id)
    ! Data variables
    s4d = (/xdim_id, ydim_id, zdim_id, tdim_id/) ! 2 levels
    ret = nf90_def_var(file_id, 'q',   nf90_float, s4d, q_id)
    ret = nf90_def_var(file_id, 'u',   nf90_float, s4d, u_id)
    ret = nf90_def_var(file_id, 'v',   nf90_float, s4d, v_id)
    ret = nf90_def_var(file_id, 'vor', nf90_float, s4d, vor_id)
    ret = nf90_def_var(file_id, 'adv', nf90_float, s4d, adv_id)
    ret = nf90_def_var(file_id, 'qx', nf90_float, s4d, qx_id)
    ret = nf90_def_var(file_id, 'qy', nf90_float, s4d, qy_id)
    s3d = (/xdim_id, ydim_id, tdim_id/) ! 1 level
    ret = nf90_def_var(file_id, 'f',   nf90_float, s3d, f_id)   ! forcing perturbations
    s1d = (/tdim_id/)
    ret = nf90_def_var(file_id, 'eke', nf90_float, s1d, eke_id)
    ret = nf90_def_var(file_id, 'cfl', nf90_float, s1d, cfl_id)
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
    ret = nf90_put_att(file_id, adv_id, 'units',     '1/s^2')
    ret = nf90_put_att(file_id, qx_id, 'units',      '1/s*m')
    ret = nf90_put_att(file_id, qy_id, 'units',      '1/s*m')
    ret = nf90_put_att(file_id, f_id,   'units',     '1/s^2')
    ret = nf90_put_att(file_id, eke_id, 'units',     'J/m^2')
    ret = nf90_put_att(file_id, q_id,   'long_name', 'total potential vorticity')
    ret = nf90_put_att(file_id, u_id,   'long_name', 'total zonal wind')
    ret = nf90_put_att(file_id, v_id,   'long_name', 'total meridional wind')
    ret = nf90_put_att(file_id, vor_id, 'long_name', 'total relative vorticity')
    ret = nf90_put_att(file_id, f_id,   'long_name', 'potential vorticity forcing in upper layer')
    ret = nf90_put_att(file_id, adv_id, 'long_name', 'eddy advection')
    ret = nf90_put_att(file_id, qx_id, 'long_name',      'eddy zonal pv gradient')
    ret = nf90_put_att(file_id, qy_id, 'long_name',      'eddy meridional pv gradient')
    ret = nf90_put_att(file_id, eke_id, 'long_name', 'eddy kinetic energy')
    ret = nf90_put_att(file_id, cfl_id, 'long_name', 'cfl number')
    ! Finished defining
    ret = nf90_enddef(file_id)
    ! Write dimension variables (must come after define mode or they will be empty)
    ret = nf90_put_var(file_id, x_id, x_cart, start=(/1/), count=(/imax/))
    ret = nf90_put_var(file_id, y_id, y_cart, start=(/1/), count=(/jmax/))
    ret = nf90_put_var(file_id, z_id, (/1,2/), start=(/1/), count=(/2/))
  endif

  !    ---- Save data ----
  ! Note if you mess up the dimensionality of input arrays, error message will be mysterious:
  ! "There is no matching specific function for this generic function reference."
  ! Add current day to time variable
  ret = nf90_put_var(file_id, t_id, (/day/), start=(/t_io/), count=(/1/))
  ! Add model data
  ret = nf90_put_var(file_id, q_id, qfull1_cart, start=(/1,1,1,t_io/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, q_id, qfull2_cart, start=(/1,1,2,t_io/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, u_id, ufull1_cart, start=(/1,1,1,t_io/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, u_id, ufull2_cart, start=(/1,1,2,t_io/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, v_id, v1_cart, start=(/1,1,1,t_io/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, v_id, v2_cart, start=(/1,1,2,t_io/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, vor_id, vorfull1_cart, start=(/1,1,1,t_io/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, vor_id, vorfull2_cart, start=(/1,1,2,t_io/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, adv_id, adv1_cart, start=(/1,1,1,t_io/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, adv_id, adv2_cart, start=(/1,1,2,t_io/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, qx_id, qx1_cart, start=(/1,1,1,t_io/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, qx_id, qx2_cart, start=(/1,1,2,t_io/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, qy_id, qy1_cart, start=(/1,1,1,t_io/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, qy_id, qy2_cart, start=(/1,1,2,t_io/), count=(/imax,jmax,1,1/))
  ret = nf90_put_var(file_id, f_id, force1_cart(:,:,1), start=(/1,1,t_io/), count=(/imax,jmax,1/))
  ret = nf90_put_var(file_id, eke_id, energy, start=(/t_io/), count=(/1/))
  ret = nf90_put_var(file_id, cfl_id, cfl, start=(/t_io/), count=(/1/))

  !    ---- Increment save time ----
  t_io = t_io+1

end subroutine
end module
