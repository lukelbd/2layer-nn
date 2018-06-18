module FTM
contains

         subroutine srcft(x,y,n,as,hx)

        Use MKL_DFTI

        implicit none
        
        type(DFTI_DESCRIPTOR), POINTER :: HX
        Integer :: Status
         integer :: n,i,iim,iip,L(1)
         double precision :: x(n),x1(n+2)
         double complex :: y(1+n/2),ur,ui
         real :: as

         L(1) = n
         ur = (1.,0)
         ui = (0.,1.)

!   Status = DftiCreateDescriptor( HX, DFTI_SINGLE, &
!      DFTI_REAL, 1, L)
!   Status = DftiCommitDescriptor( HX )
    Status = DftiComputeForward( HX, X, X1)
!   Status = DftiFreeDescriptor( HX )
!        call f2trf(n,x,x,wfftr)

         y(1) = ur*x1(1) 
         y(1) = y(1)*as

         do 67 i = 2,n/2
          iim = 2*i-1
          iip = 2*i
         y(i) = ur*x1(iim)+ui*x1(iip)
         y(i) = y(i)*as
  67     continue
         y(1+n/2) = ur*x1(n+1)*as

         return
         end subroutine srcft
!***********************************************************
         subroutine scrft(y,x1,n,as,hx)

        Use MKL_DFTI
        implicit none

        type(DFTI_DESCRIPTOR), POINTER :: HX
        Integer :: Status
         integer :: n,i,k,iim,iip
         double precision :: x(n+2),x1(n)
         double complex :: y(1+n/2)
         real :: as

         x(1) = real(y(1))
         x(2) = 0.
         do 67 i = 2,n/2
          iim = 2*i-1
          iip = 2*i
         x(iim) = real(y(i))
         x(iip) = aimag(y(i))
  67     continue
       
         x(n+1) = real(y(1+n/2))
         x(n+2) = 0.

!   Status = DftiCreateDescriptor( HX, DFTI_SINGLE, &
!      DFTI_REAL, 1, L)
!   Status = DftiCommitDescriptor( HX )
    Status = DftiComputeBackward(HX, X, X1)
!   Status = DftiFreeDescriptor( HX )
!        call f2trb(n,x,x,wfftr)

      x1 = x1*as
!        do 168 k = 1,n
!        x1(k) = x(k)*as
! 168     continue

         return
         end subroutine scrft
         
end module
