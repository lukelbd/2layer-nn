program test 

  use mkl_dfti
  use ftm

  integer,parameter :: idft = 129, imax = 256, itrunc = 85
  complex :: one_r,one_i
  double precision :: p1(imax),p2(imax)
  double complex :: pm1(idft),pm2(idft)
  type(dfti_descriptor), pointer :: hx
  integer :: L(1)

  L(1) = imax

  Status = DftiCreateDescriptor( HX, DFTI_DOUBLE, &
    DFTI_REAL, 1, L)
  Status = DftiSetValue(HX, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
  Status = DftiCommitDescriptor( HX )


  one_r = (1.,0.)
  one_i = (0.,1.)

  !  **** Transform u,v,vort to physical space ****

  pm1 = zero
  pm1(2) = one_r
  pm2 = zero
  pm2(2) = one_i

  do i = 1,idft
  if(i.gt.1.and.i.lt.idft) then
    fac = 0.5
  else
    fac = 1.0
  endif

  pm1(i) = pm1(i)*fac
  pm2(i) = pm2(i)*fac

  if(i.eq.itrunc) then
    pm1(i) = one_r*real(pm1(i))
    pm2(i) = one_r*real(pm2(i))
  endif
  enddo

  do i = itrunc+1,idft
  pm1(i) = zero
  pm2(i) = zero
  enddo

  ! *** Inverse Fourier Transform ****

  as = 1.
  call scrft(pm1,p1,imax,as,hx)
  call scrft(pm2,p2,imax,as,hx)

  do i = 1,imax
  write(*,*) 'itest ',i,p1(i),p2(i)
  enddo
  Status = DftiFreeDescriptor( HX )

  stop
  end

