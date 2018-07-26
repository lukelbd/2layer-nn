!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module calculates a bunch of diagnostics required for prognostics.f90
! including advection and PV injection forcing. Also returns cartesian x-y
! data suitable for output.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note that all q variables contain *dq/dy*, not q itself.
! Walkthrough of steps:
!  1. Invert the dq'/dy and dabar/dy equations to obtain wind, streamfunction, et cetera
!     from the fully spectral data.
!  2. Perform *backward* trigonometric (sines only or cosines only) transforms
!     on the y-dimension of the inverted data.
!  3. Perform *backward* Fourier transform on the x-dimension of the inverted
!     data.
!  4. Now we have Cartesian data; calculate PV injection, store some values
!     in arrays suitable for output.
!  5. Calculate advection terms (for q') and flux convergence terms (for qbar).
!     Requires another *forward* Fourier transform of the u' data and
!     *forward* trigonometric transforms?
!     This is done because we need to calculate ***products*** of two terms
!     in Cartesian coordinates, but calculate ***derivs*** in spectral space.
!     Why do we perform forward transform here? Didn't we alredy have that data?
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For more information on MKL utilities:
!  1. run 'module load mkl'
!  2. go to $MKLROOT/include, and check out the
!     mkl_dfti.f90 and mkl_trig_transforms.f90 files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Warning: Several tricks to be aware of
! We always calculate u ***relative to a background shear***. Need to add that
! to wind in the upper layer to get its true value.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module diagnostics
use global_variables
use spectral ! includes utilities: srcft, scrft, ftt, btt, ftt_rcft, btt_crft
use mkl_dfti            ! includes some boilerplate stuff, and Fourier transforms
use mkl_trig_transforms ! includes trig transforms
contains

! Get diagnostics
subroutine diag(hcr, hrc)

  ! Scalars
  implicit none
  type(dfti_descriptor), pointer :: hcr, hrc ! handles for real-to-complex and complex-to-real fourier transforms
  real :: y, rkk, ell, fac, s, umax
  integer :: i, j, wcos, wsin
  complex :: pb, pc ! temporary scalars during PV inversion

  ! For PV injection
  integer,allocatable :: iseed(:)
  integer             :: isize,idate(8)
  real                :: anglex,angley,amp_rand

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Reset previous inverted data
  ! Is this necessary? Don't think so, consider deleting.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ub_tt      = 0.
  uc_tt      = 0.
  ubar1_cart = 0.
  ubar2_cart = 0.
  qbar1_tt   = 0.
  qbar2_tt   = 0.
  ubar1_tt   = 0.
  ubar2_tt   = 0.
  psi1_sp    = zero
  psi2_sp    = zero
  vor1_sp    = zero
  vor2_sp    = zero
  u1_sp      = zero
  u2_sp      = zero
  v1_sp      = zero
  v2_sp      = zero
  qx1_sp     = zero
  qx2_sp     = zero
  qy1_sp     = zero
  qy2_sp     = zero

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Invert d/dy zonal mean PV equation to get zonal mean u
  ! Also integrate d/dy zonal mean PV to get zonal mean PV
  ! Note we truncate in j here
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do j = 2,jtrunc
    ell = el*float(j-1) ! angular wavenumber
    ub_tt(j) = (qybar1_tt(j,3) + qybar2_tt(j,3))/(ell*ell)
    uc_tt(j) = (qybar1_tt(j,3) - qybar2_tt(j,3))/(ell*ell + (2./(rd*rd)))
    qbar1_tt(j) = qybar1_tt(j,3)/ell ! d/dx(f(t)) = wavenum*F(wavenum) for sine transform (note wavenum is angular here)
    qbar2_tt(j) = qybar2_tt(j,3)/ell
    ubar1_tt(j) = 0.5*(ub_tt(j) + uc_tt(j))
    ubar2_tt(j) = 0.5*(ub_tt(j) - uc_tt(j))
    vorbar1_tt(j) = (ell*ell)*ubar1_tt(j) ! shear vorticity, e.g. due to jet
    vorbar2_tt(j) = (ell*ell)*ubar2_tt(j)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Invert PV equation to get streamfunction
  ! Also obtain u/v, and q derivatives in x/y (derivative in spectral coordinates)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transform eddy into cartesian space
  ! Used to just assign (1,:) and (:,1) zero but am fairly certain
  ! this was minor bug, and we were transforming arrays populated by 'empty'
  ! values above the trunc number (generally random, extremely small numbers)
  do j = 2,jtrunc
    ell = el*float(j-1)
    do i = 2,itrunc
      rkk = rk*float(i-1)
      pb  = -(q1_sp(i,j,3)+q2_sp(i,j,3))/(rkk**2 + ell**2)
      pc  = -(q1_sp(i,j,3)-q2_sp(i,j,3))/(rkk**2 + ell**2  + (2./(rd*rd)))
      psi1_sp(i,j)  = 0.5*(pb+pc)
      psi2_sp(i,j)  = 0.5*(pb-pc)
      vor1_sp(i,j) = -(rkk**2+ell**2)*psi1_sp(i,j)
      vor2_sp(i,j) = -(rkk**2+ell**2)*psi2_sp(i,j)
      u1_sp(i,j)    = -ell*psi1_sp(i,j)
      u2_sp(i,j)    = -ell*psi2_sp(i,j)
      v1_sp(i,j)    = one_i*rkk*psi1_sp(i,j)
      v2_sp(i,j)    = one_i*rkk*psi2_sp(i,j)
      qx1_sp(i,j)   = one_i*rkk*q1_sp(i,j,3)
      qx2_sp(i,j)   = one_i*rkk*q2_sp(i,j,3)
      qy1_sp(i,j)   = ell*q1_sp(i,j,3)
      qy2_sp(i,j)   = ell*q2_sp(i,j,3)
    enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transform zonal mean to cartesian space 
  ! Also, after transforming mean qbar anomalies, add the basic-state
  ! qbar from background. Can be used to output data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Inverse sine transforms
  ! Suitable if edges are always small/zero
  tt_type=0
  call btt(ubar1_tt, ubar1_cart, tt_type, jmax, jtrunc)
  call btt(ubar2_tt, ubar2_cart, tt_type, jmax, jtrunc)
  call btt(vorbar1_tt, vorbar1_cart, tt_type, jmax, jtrunc)
  call btt(vorbar2_tt, vorbar2_cart, tt_type, jmax, jtrunc)
  call btt(qybar1_tt(:,3), qybar1_cart, tt_type, jmax, jtrunc) ! needed for total advection
  call btt(qybar2_tt(:,3), qybar2_cart, tt_type, jmax, jtrunc)

  ! Inverse cosine transforms
  ! Suitable if edges are always non-zero
  tt_type=1
  call btt(qbar1_tt, qbar1_cart, tt_type, jmax, jtrunc)
  call btt(qbar2_tt, qbar2_cart, tt_type, jmax, jtrunc)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transform 2D data to cartesian space
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Inverse sine transforms
  ! Suitable if edges are always small/zero
  ! Note we are transforming the x-Fourier coefficients here
  tt_type=0
  call btt_crft(v1_sp, v1_cart, tt_type, imax, jmax, itrunc, jtrunc, hcr)
  call btt_crft(q1_sp, q1_cart, tt_type, imax, jmax, itrunc, jtrunc, hcr)
  call btt_crft(qx1_sp, qx1_cart, tt_type, imax, jmax, itrunc, jtrunc, hcr)
  call btt_crft(psi1_sp, psi1_cart, tt_type, imax, jmax, itrunc, jtrunc, hcr)
  call btt_crft(vor1_sp, vor1_cart, tt_type, imax, jmax, itrunc, jtrunc, hcr)
  call btt_crft(v2_sp, v2_cart, tt_type, imax, jmax, itrunc, jtrunc, hcr)
  call btt_crft(q2_sp, q2_cart, tt_type, imax, jmax, itrunc, jtrunc, hcr)
  call btt_crft(qx2_sp, qx2_cart, tt_type, imax, jmax, itrunc, jtrunc, hcr)
  call btt_crft(psi2_sp, psi2_cart, tt_type, imax, jmax, itrunc, jtrunc, hcr)
  call btt_crft(vor2_sp, vor2_cart, tt_type, imax, jmax, itrunc, jtrunc, hcr)
  ! Inverse cosine transforms
  ! Suitble if edges are always non-zero
  tt_type=1
  call btt_crft(u1_sp, u1_cart, tt_type, imax, jmax, itrunc, jtrunc, hcr)
  call btt_crft(qy1_sp, qy1_cart, tt_type, imax, jmax, itrunc, jtrunc, hcr)
  call btt_crft(u2_sp, u2_cart, tt_type, imax, jmax, itrunc, jtrunc, hcr)
  call btt_crft(qy2_sp, qy2_cart, tt_type, imax, jmax, itrunc, jtrunc, hcr)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Construct 'full' values from anomalous values
  ! Not used to derive stuff, but probably will want to output these
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do j = 1,jmax
    ! Total wind
    ufull1_cart(:,j) = u1_cart(:,j) + ubar1_cart(j) + shear ! add basic state shear
    ufull2_cart(:,j) = u2_cart(:,j) + ubar2_cart(j)
    ! Total PV
    qfull1_cart(:,j) = q1_cart(:,j) + qbar1_cart(j) + (beta + shear/(rd*rd))*y_cart(j)
    qfull2_cart(:,j) = q2_cart(:,j) + qbar2_cart(j) + (beta - shear/(rd*rd))*y_cart(j)
    ! Total relative vorticity
    vorfull1_cart(:,j) = vor1_cart(:,j) + vorbar1_cart(j)
    vorfull2_cart(:,j) = vor2_cart(:,j) + vorbar2_cart(j)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Stochastic forcing (small-scale narrow band)
  ! Inject cosines and sines onto the raw 2D data, with fixed memory 0.5
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do j = 1,jmax
    ! Note this injects the *forcing*, not q anomalies themselves
    ! * Must appear here, because we have a real mask that limits appearence
    !   of PV perturbations into a narrow 'jet' band
    ! * Below can be used to have arbitrary memory; just replace where you see 0.5
    !   force1_cart(i) = exp(-dt/tau_i)*force1_cart(i,j,1)
    !   + (1.-exp(-dt/tau_i))*amp_i*mask_i(j)*  &
    iseed = iseed*(idate(8)-500)
    call date_and_time(values=idate)
    call random_seed(size=isize)
    call random_seed(put=iseed)
    call random_number(anglex)
    call random_number(angley)
    do i = 1,imax
      force1_cart(i,j,2) = 0.5*force1_cart(i,j,1) ! initialize with previous state
      do wcos = wmin_i,wmax_i
        rkk = rk*float(wcos-1)
        do wsin = wmin_i,wmax_i
          ell = el*float(wsin-1)
          call random_number(amp_rand)
          force1_cart(i,j,2) = force1_cart(i,j,2) + 0.5*amp_i*amp_rand*mask_i(j)  &
            * sin(ell*float(j-1)/float(jmax-1)+angley)  &
            * cos(rkk*float(i-1)/float(imax)+anglex)
        enddo
      enddo
    enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Get zonal advection and forcing
  ! Also get some other stuff
  ! Will have to transform u*q flux back to spectral space
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  tt_type=0
  do j = 1,jmax
    adv1_cart(:,j) = qx1_cart(:,j)*(u1_cart(:,j) + ubar1_cart(j) + shear) &
                 + v1_cart(:,j)*(qy1_cart(:,j) + beta + (shear/(rd*rd)) + qybar1_cart(j))
    adv2_cart(:,j) = qx2_cart(:,j)*(u2_cart(:,j) + ubar2_cart(j)) &
                 + v2_cart(:,j)*(qy2_cart(:,j) + beta - (shear/(rd*rd)) + qybar2_cart(j))
  enddo
  call ftt_rcft(force1_cart(:,:,2), force1_sp, tt_type, imax, jmax, itrunc, jtrunc, hrc)
  call ftt_rcft(adv1_cart, tmp1_sp, tt_type, imax, jmax, itrunc, jtrunc, hrc)
  call ftt_rcft(adv2_cart, tmp2_sp, tt_type, imax, jmax, itrunc, jtrunc, hrc)
  adv1_sp(:,:,3) = tmp1_sp(:,:)
  adv2_sp(:,:,3) = tmp2_sp(:,:)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Get meridional q flux convergence
  ! Have to transform back to trig space
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  tt_type=0
  call ftt(sum(v1_cart*q1_cart,1)/float(imax), tmp1_tt, tt_type, jmax, jtrunc)
  call ftt(sum(v2_cart*q2_cart,1)/float(imax), tmp2_tt, tt_type, jmax, jtrunc)
  do j = 2,jtrunc
    ell = el*float(j-1)
    qyyflux1_tt(j,3) = -ell*ell*tmp1_tt(j) ! these are 2nd derivatives
    qyyflux2_tt(j,3) = -ell*ell*tmp2_tt(j) ! convergence of the flux?
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Print diagnostic output, useful for monitoring the situation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! write(*,*) 'max flux',maxval(qyyflux1_tt),minval(qyyflux1_tt)
  if(mod(int(t),int(dt)*1000).eq.1) then
    write(*,*) 'Printing zonal mean diagnostics.'
    do j = 1,jmax
      write(*,*) j, 'ubar1 = ', ubar1_tt(j)+shear, 'qbar1 = ', qbar1_tt(j)
    enddo
  endif
  day    = float(t)/(3600.*24.)
  cfl    = umax*dt/dx
  umax   = max(maxval(ufull1_cart),maxval(ufull2_cart))
  energy = 0.5*sum(v1_cart**2 + v2_cart**2 + u1_cart**2 + u2_cart**2)/float(imax*jmax)
  write(*,677) day, energy, umax, cfl
  677 format("days = ", 1f8.3, " eke = ", 1p1e13.5, " umax = ", 1f3.3, " cfl = ", 1f3.3)
  ! 1p ensures non-zero digit to left of decimal
  return
end subroutine

end module
