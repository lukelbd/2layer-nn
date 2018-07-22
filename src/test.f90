program test 

  use MKL_DFTI
  use FTM

  integer,parameter :: idft = 129, imax = 256, mmax = 85
  complex :: ur,ui
  double precision :: p1(imax),p2(imax)
  double complex :: pm1(idft),pm2(idft)
  type(DFTI_DESCRIPTOR), POINTER :: HX
  integer :: L(1)

  L(1) = imax

  Status = DftiCreateDescriptor( HX, DFTI_DOUBLE, &
    DFTI_REAL, 1, L)
  Status = DftiSetValue(HX, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
  Status = DftiCommitDescriptor( HX )


  ur = (1.,0.)
  ui = (0.,1.)

  !  **** Transform u,v,vort to physical space ****

  pm1 = zero
  pm1(2) = ur
  pm2 = zero
  pm2(2) = ui

  do i = 1,idft
  if(i.gt.1.and.i.lt.idft) then
    fac = 0.5
  else
    fac = 1.0
  endif

  pm1(i) = pm1(i)*fac
  pm2(i) = pm2(i)*fac

  if(i.eq.mmax) then
    pm1(i) = ur*real(pm1(i))
    pm2(i) = ur*real(pm2(i))
  endif
  enddo

  do i = mmax+1,idft
  pm1(i) = zero
  pm2(i) = zero
  enddo

  ! *** Inverse Fourier Transform ****

  as = 1.
  call scrft(pm1,p1,imax,as,hx)
  call scrft(pm2,p2,imax,as,hx)

  do i = 1,imax
  write(6,*) 'itest ',i,p1(i),p2(i)
  enddo
  Status = DftiFreeDescriptor( HX )

  stop
  end

