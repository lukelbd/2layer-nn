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
use transforms ! includes utilities: srcft, scrft, ftt, btt, ftt_rcft, btt_crft
use mkl_dfti ! includes some boilerplate stuff, and Fourier transforms
contains

! Get diagnostics
subroutine diag
  ! Scalars
  implicit none
  real :: y, rkk, ell, fac, s, qplus, qminus
  integer :: i, j, wx_i, wy_i
  complex :: pb_sp, pc_sp ! temporary scalars during PV inversion
  real    :: pb_tt, pc_tt
  ! For PV injection
  integer,allocatable :: iseed(:)
  integer             :: isize, idate(8)
  real                :: fact, ang_x, ang_y, amp_rand

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Reset previous inverted data
  ! Is this necessary? Don't think so, consider deleting.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  qbar1_tt   = 0.
  qbar2_tt   = 0.
  ubar1_tt   = 0.
  ubar2_tt   = 0.
  vorbar1_tt = 0.
  vorbar2_tt = 0.
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
  ! Invert d/dy zonal mean PV equation to get zonal mean u, vor, and pv
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do j = 2,jtrunc
    ell = el*float(j-1) ! angular wavenumber
    pb_tt = (qybar1_tt(j,2) + qybar2_tt(j,2))/(ell*ell)
    pc_tt = (qybar1_tt(j,2) - qybar2_tt(j,2))/(ell*ell + (2./(rd*rd)))
    ubar1_tt(j) = 0.5*(pb_tt + pc_tt)
    ubar2_tt(j) = 0.5*(pb_tt - pc_tt)
    qbar1_tt(j) = qybar1_tt(j,2)/ell ! d/dx(f(t)) = wavenum*F(wavenum) for sine transform (note wavenum is angular here)
    qbar2_tt(j) = qybar2_tt(j,2)/ell
    vorbar1_tt(j) = (ell*ell)*ubar1_tt(j) ! shear vorticity, e.g. due to jet
    vorbar2_tt(j) = (ell*ell)*ubar2_tt(j)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Invert PV equation to get u, v, streamfunction, vorticity, and q gradient
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transform eddy into cartesian space
  ! Used to just assign (1,:) and (:,1) zero but am fairly certain
  ! this was minor bug, and we were transforming arrays populated by 'empty'
  ! values above the trunc number (generally random, extremely small numbers)
  do j = 2,jtrunc
    ell = el*float(j-1)
    do i = 2,itrunc
      rkk = rk*float(i-1)
      qx1_sp(i,j)  = one_i*rkk*q1_sp(i,j,2)
      qx2_sp(i,j)  = one_i*rkk*q2_sp(i,j,2)
      qy1_sp(i,j)  = ell*q1_sp(i,j,2)
      qy2_sp(i,j)  = ell*q2_sp(i,j,2)
      pb_sp = -(q1_sp(i,j,2) + q2_sp(i,j,2))/(rkk**2 + ell**2)
      pc_sp = -(q1_sp(i,j,2) - q2_sp(i,j,2))/(rkk**2 + ell**2  + (2./(rd*rd)))
      psi1_sp(i,j) = 0.5*(pb_sp + pc_sp)
      psi2_sp(i,j) = 0.5*(pb_sp - pc_sp)
      vor1_sp(i,j) = -(rkk**2 + ell**2)*psi1_sp(i,j)
      vor2_sp(i,j) = -(rkk**2 + ell**2)*psi2_sp(i,j)
      u1_sp(i,j)   = -ell*psi1_sp(i,j)
      u2_sp(i,j)   = -ell*psi2_sp(i,j)
      v1_sp(i,j)   = one_i*rkk*psi1_sp(i,j)
      v2_sp(i,j)   = one_i*rkk*psi2_sp(i,j)
    enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transform zonal means to cartesian space 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Inverse sine transforms (suitable for (d/dy)^n of qybar, n *even*)
  tt_type=0
  call btt(qybar1_tt(:,2), qybar1_cart, tt_type, jmax) ! needed for total advection
  call btt(qybar2_tt(:,2), qybar2_cart, tt_type, jmax)
  call btt(vorbar1_tt, vorbar1_cart, tt_type, jmax)
  call btt(vorbar2_tt, vorbar2_cart, tt_type, jmax)
  call btt(ubar1_tt, ubar1_cart, tt_type, jmax)
  call btt(ubar2_tt, ubar2_cart, tt_type, jmax)
  ! Inverse cosine transforms (suitable for (d/dy)^n of dqbardy, n *odd*)
  tt_type=1
  call btt(qbar1_tt, qbar1_cart, tt_type, jmax)
  call btt(qbar2_tt, qbar2_cart, tt_type, jmax)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transform 2D data to cartesian space
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Inverse sine transforms (suitable for (d/dy)^n of q, n *even*)
  ! Note we are transforming the x-Fourier coefficients here
  tt_type=0
  call btt_crft(q1_sp(:,:,2), q1_cart, tt_type, imax, jmax, itrunc, hcr)
  call btt_crft(q2_sp(:,:,2), q2_cart, tt_type, imax, jmax, itrunc, hcr)
  call btt_crft(qx1_sp, qx1_cart, tt_type, imax, jmax, itrunc, hcr)
  call btt_crft(qx2_sp, qx2_cart, tt_type, imax, jmax, itrunc, hcr)
  call btt_crft(vor1_sp, vor1_cart, tt_type, imax, jmax, itrunc, hcr)
  call btt_crft(vor2_sp, vor2_cart, tt_type, imax, jmax, itrunc, hcr)
  call btt_crft(psi1_sp, psi1_cart, tt_type, imax, jmax, itrunc, hcr)
  call btt_crft(psi2_sp, psi2_cart, tt_type, imax, jmax, itrunc, hcr)
  call btt_crft(v1_sp, v1_cart, tt_type, imax, jmax, itrunc, hcr)
  call btt_crft(v2_sp, v2_cart, tt_type, imax, jmax, itrunc, hcr)
  ! Inverse cosine transforms (suitable for (d/dy)^n of q, n *odd*)
  tt_type=1
  call btt_crft(u1_sp, u1_cart, tt_type, imax, jmax, itrunc, hcr)
  call btt_crft(u2_sp, u2_cart, tt_type, imax, jmax, itrunc, hcr)
  call btt_crft(qy1_sp, qy1_cart, tt_type, imax, jmax, itrunc, hcr)
  call btt_crft(qy2_sp, qy2_cart, tt_type, imax, jmax, itrunc, hcr)

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
  ! Inject cosines and sines onto the raw 2D data, with some autocorrelation timescale.
  ! * Note this injects the *forcing*, not q anomalies themselves
  ! * Why do this in cartesian space instead of spectral? Slightly easier
  !   since the mask is also cartesian.
  ! * Below version used forcing with memory at every timestep; was way too
  !   often and generated totally weird zonal winds.
  ! force1_cart(:,:,2) = exp(-dt/tau_i)*force1_cart(:,:,1) ! with autocorrelation
  ! + (1.0-exp(-dt/tau_i))*amp_i*amp_rand*mask_i(j) & ! with autocorrelation
  ! call random_number(trial)
  ! if (t.eq.0 .or. trial.le.p) then
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  force1_sp = zero ! initialize spectral forcing array
  call date_and_time(values=idate)
  iseed = iseed*(idate(8)-500) ! initialize with pseudo-random 'seed'
  call random_seed(size=isize)
  call random_seed(put=iseed)
  if (contin_i .or. mod(t,dt_i).eq.0) then
    ! Initialize; the interval-inject version has no memory
    if (contin_i) then
      force1_cart(:,:,2) = exp(-dt/tau_i)*force1_cart(:,:,1) ! with autocorrelation
      fact = (1.0-exp(-dt/tau_i))
    else
      print *, 'Injecting pv.'
      force1_cart(:,:,2) = 0.0
      fact = 1.0
    endif
    do wx_i = wmin_i,wmax_i
      do wy_i = wmin_i,wmax_i
        ! Assign random phase/amplitude to this particular zonal/meridional wavenumber
        call random_number(ang_x)
        call random_number(ang_y)
        call random_number(amp_rand)
        ! Apply waves to Cartesian grid
        ! print *, 'wavenums', wx_i, wy_i
        ! do j = 1,jmax
        !   print *, 'y', sin(wy_i*2*pi*(float(j-1)/float(jmax-1)+ang_y))
        ! enddo
        ! do i = 1,imax
        !   print *, 'x', sin(wx_i*2*pi*(float(i-1)/float(imax)+ang_x))
        ! enddo
        do j = 1,jmax
          do i = 1,imax
            force1_cart(i,j,2) = force1_cart(i,j,2) &
              + fact*amp_i*amp_rand*mask_i(j) &
              * sin(wy_i*2*pi*(float(j-1)/float(jmax-1)+ang_y))  &
              * sin(wx_i*2*pi*(float(i-1)/float(imax)+ang_x))
          enddo
        enddo
      enddo
    enddo
    call ftt_rcft(force1_cart(:,:,2), force1_sp, tt_type, imax, jmax, itrunc, hrc)
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Get advection terms in spectral space
  ! Note advection is positive here, i.e. is on 'LHS' of dq/dt equation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  tt_type=0
  do j = 1,jmax
    adv1_cart(:,j) = qx1_cart(:,j)*(u1_cart(:,j) + ubar1_cart(j) + shear) &
       + v1_cart(:,j)*(qy1_cart(:,j) + beta + (shear/(rd*rd)) + qybar1_cart(j))
    adv2_cart(:,j) = qx2_cart(:,j)*(u2_cart(:,j) + ubar2_cart(j)) &
       + v2_cart(:,j)*(qy2_cart(:,j) + beta - (shear/(rd*rd)) + qybar2_cart(j))
  enddo
  call ftt_rcft(adv1_cart, tmp1_sp, tt_type, imax, jmax, itrunc, hrc)
  call ftt_rcft(adv2_cart, tmp2_sp, tt_type, imax, jmax, itrunc, hrc)
  do i=2,itrunc
    do j=2,jtrunc
      adv1_sp(i,j,3) = tmp1_sp(i,j)
      adv2_sp(i,j,3) = tmp2_sp(i,j)
    enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Get meridional q flux convergence for each zonal band
  ! Note this is the only way the q' array communicates with the dqbar/dy array
  ! Have to transform back to trig space
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  tt_type=0
  qflux1_cart = sum(v1_cart*q1_cart,1)/float(imax)
  qflux2_cart = sum(v2_cart*q2_cart,1)/float(imax)
  call ftt(qflux1_cart, tmp1_tt, tt_type, jmax)
  call ftt(qflux2_cart, tmp2_tt, tt_type, jmax)
  do j = 2,jtrunc
    ell = el*float(j-1)
    qyyflux1_tt(j,3) = -ell*ell*tmp1_tt(j) ! these are 2nd derivatives
    qyyflux2_tt(j,3) = -ell*ell*tmp2_tt(j) ! convergence of the flux?
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Scalars, including eke
  ! Store in arrays so they can be saved easily
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  energy(1) = 0.5*sum(v1_cart**2 + v2_cart**2 + u1_cart**2 + u2_cart**2)/float(imax*jmax)
  umax(1)   = max(maxval(ufull1_cart),maxval(ufull2_cart))
  cfl(1)    = umax(1)*dt/dx

end subroutine
end module
