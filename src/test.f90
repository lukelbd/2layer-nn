program test 

  use MKL_DFTI
  use FTM

  integer,parameter :: imax = 129, imx = 256, mmax = 85
  complex :: ur,ui
  double precision :: p1(imx),p2(imx)
  double complex :: pm1(imax),pm2(imax)
  type(DFTI_DESCRIPTOR), POINTER :: HX
  integer :: L(1)

  L(1) = imx

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

  do i = 1,imax
  if(i.gt.1.and.i.lt.imax) then
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

  do i = mmax+1,imax
  pm1(i) = zero
  pm2(i) = zero
  enddo

  ! *** Inverse Fourier Transform ****

  as = 1.
  call scrft(pm1,p1,imx,as,hx)
  call scrft(pm2,p2,imx,as,hx)

  do i = 1,imx
  write(6,*) 'itest ',i,p1(i),p2(i)
  enddo
  Status = DftiFreeDescriptor( HX )

  stop
  end

