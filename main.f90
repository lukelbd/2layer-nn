!****** Nonlinear 2 layer baroclinic instability on the beta-plane.  Periodic in x and rigid walls at north and south *****   

        program main

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

!       implicit none

! ****** Initialize fields ****  
  
         call initial1

        Status = DftiCreateDescriptor( HX, DFTI_DOUBLE, &
          DFTI_REAL, 1, L)
        Status = DftiSetValue(HX, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
        Status = DftiCommitDescriptor( HX )
        Status = DftiCreateDescriptor( HY, DFTI_DOUBLE, &
          DFTI_REAL, 1, L)
        Status = DftiSetValue(HY, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
        Status = DftiCommitDescriptor( HY )

! ****** Integration starts *****


        do m = mstart, mend 
!       do m = 1,20
         call invert1(HX,HY)              !invert streamfunction from vorticity
          if((mod(m,md).eq.0).and.(m.ge.mds)) call dataio   !save data
    !     call energy1                     !eddy energy calculation
          call prog                        !eddy prognostic equation
          call forward1                     !time forwarding

         write(6,677) m, energy2
  677   format("normal end m = ",i8,"   eke =",1p1e13.5)
        enddo
        
        write(6,*) 'm =',m,'  Yeh'
        close(11)
 
         Status = DftiFreeDescriptor( HX )
         Status = DftiFreeDescriptor( HY )

      stop
    end program
