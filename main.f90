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
  
         call initial1      ! (initial.f90)

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
        
         call invert1(HX,HY)      ! (invert.f90) invert streamfunction from vorticity
          if((mod(m,md).eq.0).and.(m.ge.mds)) call dataio   ! (io.f90) save data
    !     call energy1            ! (energy.f90) eddy energy calculation
          call prog               ! (exec.f90) eddy prognostic equation
          call forward1           ! (forward.f90) time forwarding

        if((mod(m,md).eq.0).and.(m.ge.mds)) then
                  write(6,677) m, energy2
 677   format("normal end m = ",i8,"   eke =",1p1e13.5)
                write(6,*) ' max u = ',maxval(uf)
        endif
        enddo
        
        write(6,*) 'm =',m,'  Yeh'
        close(11)
 
         Status = DftiFreeDescriptor( HX )
         Status = DftiFreeDescriptor( HY )

      stop
    end program
