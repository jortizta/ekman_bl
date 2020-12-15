! C******************************************************************************|
! C duct.f, the duct-flow solvers for diablo.                        VERSION 0.9
! C These solvers were written by ? and ? (spring 2001).
! C******************************************************************************|

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine INIT_DUCT
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! C Initialize any constants here
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
  USE mg_vari, 		ONLY : INIT_FLAG
  USE mpi_var

    IMPLICIT NONE 

    INTEGER N
    
    PI=4.D0*ATAN(1.D0)

! ! At YMIN Location
! c        write(*,*) 'U_BC_YMIN: ',U_BC_YMIN
    IF (W_BC_YMIN.EQ.0) THEN
      !IF (IC_TYPE == 7 ) THEN
      JSTART=1
      !ELSE
      !JSTART=2
      !ENDIF
    ELSE IF (W_BC_YMIN.EQ.1) THEN
      JSTART=1
    ELSE IF (W_BC_YMIN.EQ.6) THEN
      JSTART=1  
    ELSE IF (W_BC_YMIN.EQ.7) THEN
      JSTART=1
    ELSE IF (W_BC_YMIN.EQ.3) THEN
      JSTART=1
    ELSE
      JSTART=2
    ENDIF
    
! ! At ZMIN Location  
! c        write(*,*) 'U_BC_ZMIN: ',U_BC_ZMIN
      IF (U_BC_ZMIN.EQ.0) THEN
	IF (IC_TYPE == 7 ) THEN
	ZSTART=2
	ELSE
	ZSTART=2
	ENDIF   
      ELSE IF (U_BC_ZMIN.EQ.1) THEN	
	ZSTART=1
      ELSE IF (U_BC_ZMIN.EQ.6) THEN	
	ZSTART=1  
      ELSE
	ZSTART=2
      ENDIF

! Now, set the indexing for the scalar equations
! At YMIN Location
      DO N=1,N_TH
	IF (TH_BC_YMIN(N).EQ.0) THEN
	  JSTART_TH(N)=1
	ELSE IF (TH_BC_YMIN(N).EQ.1) THEN
	  JSTART_TH(N)=1
        ELSE IF (TH_BC_YMIN(N).EQ.5) THEN
          JSTART_TH(N)=1
	ELSE IF (TH_BC_YMIN(N).EQ.2) THEN
	  JSTART_TH(N)=1
	ELSE IF (TH_BC_YMIN(N).EQ.6) THEN
	  JSTART_TH(N)=1
	ELSE IF (TH_BC_YMIN(N).EQ.7) THEN
	  JSTART_TH(N)=1
	ELSE IF (TH_BC_YMIN(N).EQ.8) THEN
	  JSTART_TH(N)=1
	ELSE
	  JSTART_TH(N)=2
	ENDIF
      ENDDO

! At ZMIN Location    
      DO N=1,N_TH
	IF (TH_BC_ZMIN(N).EQ.0) THEN
	  IF (IC_TYPE == 7 ) THEN
	  ZSTART_TH(N)=1
	  ELSE
	  ZSTART_TH(N)=2
	  ENDIF
	ELSE IF (TH_BC_ZMIN(N).EQ.1) THEN
	  ZSTART_TH(N)=1
	ELSE IF (TH_BC_ZMIN(N).EQ.6) THEN
	  ZSTART_TH(N)=1    
	ELSE
	  ZSTART_TH(N)=2
	ENDIF
      ENDDO      

! At YMAX Location

      IF (W_BC_YMAX.EQ.0) THEN
	IF (IC_TYPE==7)THEN
	JEND=NY
	ELSE
	JEND=NY
	ENDIF 
      ELSE IF (W_BC_YMAX.EQ.1) THEN
	JEND=NY
      ELSE IF (W_BC_YMAX.EQ.6) THEN
	JEND=NY  
      ELSE
	JEND=NY-1
      ENDIF

! At ZMAX Location
      IF (U_BC_ZMAX.EQ.0) THEN
	IF (IC_TYPE == 7)THEN 
	ZEND=NZ
	ELSE
	ZEND=NZ-1
	ENDIF 
      ELSE IF (U_BC_ZMAX.EQ.1) THEN
	ZEND=NZ
      ELSE IF (U_BC_ZMAX.EQ.6) THEN
	ZEND=NZ  
      ELSE
	ZEND=NZ-1
      ENDIF

! Set the upper and lower limits of timestepping of the scalar equations
! At YMAX Location 
        DO N=1,N_TH
	  IF (TH_BC_YMAX(N).EQ.0) THEN
	  IF (IC_TYPE == 7)THEN
	    JEND_TH(N)=NY
	  ELSE 
	    JEND_TH(N)=NY-1
	  ENDIF
	  ELSE IF (TH_BC_YMAX(N).EQ.1) THEN
	    JEND_TH(N)=NY
	  ELSE IF (TH_BC_YMAX(N).EQ.6) THEN
	    JEND_TH(N)=NY
	  ELSE IF (TH_BC_YMAX(N).EQ.5) THEN
	    JEND_TH(N)=NY
	  ELSE
	    JEND_TH(N)=NY-1
	  ENDIF
        ENDDO
        
! At ZMAX Location
        DO N=1,N_TH
	  IF (TH_BC_ZMAX(N).EQ.0) THEN
	  IF (IC_TYPE == 7 ) THEN
	    ZEND_TH(N)=NZ
	  ELSE
	    ZEND_TH(N)=NZ-1
	  ENDIF
	  ELSE IF (TH_BC_ZMAX(N).EQ.1) THEN
	    ZEND_TH(N)=NZ
	  ELSE IF (TH_BC_ZMAX(N).EQ.5) THEN
	    ZEND_TH(N)=NZ
	  ELSE IF (TH_BC_ZMAX(N).EQ.6) THEN
	    ZEND_TH(N)=NZ  
	  ELSE
	    ZEND_TH(N)=NZ-1
	  ENDIF
        ENDDO

       IF (RANK .eq. 0 ) THEN  
	  WRITE(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
		      
	  WRITE(6,*)'Boundary condition for U'
	  WRITE(6,*)'U_BC_YMIN',U_BC_YMIN,'U_BC_YMAX',U_BC_YMAX
	  WRITE(6,*)'U_BC_ZMIN',U_BC_ZMIN,'U_BC_ZMAX',U_BC_ZMAX
	  WRITE(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
		      
		      
	  WRITE(6,*)'U_BC_YMIN_C1',U_BC_YMIN_C1,'U_BC_YMAX_C1',U_BC_YMAX_C1
	  WRITE(6,*)'U_BC_ZMIN_C1',U_BC_ZMIN_C1,'U_BC_ZMAX_C1',U_BC_ZMAX_C1

	  WRITE(6,*)'Boundary condition for W'
	  WRITE(6,*)'W_BC_YMIN',W_BC_YMIN,'W_BC_YMAX', W_BC_YMAX
	  WRITE(6,*)'W_BC_ZMIN',W_BC_ZMIN,'W_BC_ZMAX', W_BC_ZMAX
	  WRITE(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
		      
		      
	  WRITE(6,*)'W_BC_YMIN_C1',W_BC_YMIN_C1,'W_BC_YMAX_C1',W_BC_YMAX_C1
	  WRITE(6,*)'W_BC_ZMIN_C1',W_BC_ZMIN_C1,'W_BC_ZMAX_C1',W_BC_ZMAX_C1


	  WRITE(6,*)'Boundary condition for V'
	  WRITE(6,*)'V_BC_YMIN',V_BC_YMIN,'V_BC_YMAX',V_BC_YMAX
	  WRITE(6,*)'V_BC_ZMIN',V_BC_ZMIN,'V_BC_ZMAX',V_BC_ZMAX
	  WRITE(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
		      
		      
	  WRITE(6,*)'V_BC_YMIN_C1',V_BC_YMIN_C1,'V_BC_YMAX_C1',V_BC_YMAX_C1
	  WRITE(6,*)'V_BC_ZMIN_C1',V_BC_ZMIN_C1,'V_BC_ZMAX_C1',V_BC_ZMAX_C1


	  WRITE(6,*)'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
		      
		    
	  WRITE(6,*) 'JSATRT', JSTART, 'JEND', JEND
	  WRITE(6,*) 'ZSATRT', ZSTART, 'ZEND', ZEND 
	  WRITE(6,*)'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
		      
		    
	  DO N=1,N_TH
	    WRITE(6,*)'Boundary condition for TH'
	    WRITE(6,*)'TH_BC_YMIN',TH_BC_YMIN(N),'TH_BC_YMAX', TH_BC_YMAX(N)
	    WRITE(6,*)'TH_BC_ZMIN',TH_BC_ZMIN(N),'TH_BC_ZMAX', TH_BC_ZMAX(N)
	    WRITE(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
			
			
	    WRITE(6,*)'TH_BC_YMIN_C1',TH_BC_YMIN_C1(N),   &
		      'TH_BC_YMAX_C1',TH_BC_YMAX_C1(N)
	    WRITE(6,*)'TH_BC_ZMIN_C1',TH_BC_ZMIN_C1(N),   &
		      'TH_BC_ZMAX_C1',TH_BC_ZMAX_C1(N)
	    WRITE(6,*)'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
		      
		      
	    WRITE(6,*) 'JSATRT_TH', JSTART_TH(N), 'JEND_TH', JEND_TH(N)
	    WRITE(6,*) 'ZSATRT_TH', ZSTART_TH(N), 'ZEND_TH', ZEND_TH(N) 
	    WRITE(6,*)'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
	  ENDDO
       ENDIF
 

  RETURN
END
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|


! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE RK_DUCT
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! C Alternative time-stepping algorithm for the duct-flow case with ADI.
! C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
! C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
  USE ADI_var
  USE mg_vari, 		ONLY : INIT_FLAG, MUDPACK_BC
  ! USE omp_lib      
  USE les_chan_var
  USE mpi_var

    IMPLICIT NONE
    
    INTEGER I,J,K,N      
    REAL*8 TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, D
    REAL(r8) CTHX_MEAN(1:NZ,1:NY)
 
!       REAL*8 gt1,gt2,deltaU
!       
!       PI = ATAN(1.d0)*4.d0
!       deltaU = 0.565d0 - 0.125d0
!       gt1=43.942d0
!       gt2=gt1+PI/2.d0
!       !Ramping up the barotropic velocity       
!       IF( TIME.LE.gt2) THEN
!       	IF( TIME.GE.gt1 ) THEN
! 	  		AMP_OMEGA0 = - ( 0.125d0 + deltaU*dtanh( 2.5d0*(time-gt1)/(gt2-gt1) )  )	!0.565-0.125=0.44
! 	  		if( RANK.EQ.0 ) write(6,*) 'AMP_OMEGA0=', AMP_OMEGA0
! 	  	ENDIF
! 	  ELSE
! 	  		AMP_OMEGA0 = - ( 0.125d0 + deltaU )
! 	  ENDIF
	  
! C Communicate the information between ghost cells 
! 
! c      CALL GHOST_CHAN_MPI
! 
! C Define the constants that are used in the time-stepping
! C For reference, see Numerical Renaissance
    TEMP1 = NU * H_BAR(RK_STEP) / 2.0d0
    TEMP2 = H_BAR(RK_STEP) / 2.0d0
    TEMP3 = ZETA_BAR(RK_STEP) * H_BAR(RK_STEP)
    TEMP4 = H_BAR(RK_STEP)
    TEMP5 = BETA_BAR(RK_STEP) * H_BAR(RK_STEP)

    
! C 	this is required convective boundary condition 
    dtc = TEMP4      
!       wtime =  omp_get_wtime ( )
!       
! C 	First, we will compute the explicit RHS terms and store in Ri
! C 	Note, Momentum equation and hence the RHS is evaluated at the
! C 	corresponding velocity points.

    
    CR1X(:,:,:) = (0.d0,0.d0)
    CR2X(:,:,:) = (0.d0,0.d0)
    CR3X(:,:,:) = (0.d0,0.d0)    
    

    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NX2P
	  CR1X(I,K,J)=INT_JACOB(K,J)*CU1X(I,K,J)
	ENDDO
      ENDDO
    ENDDO

    DO J= JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NX2P
	  CR2X(I,K,J)=INT_JACOB(K,J)*CU2X(I,K,J)
	ENDDO
      ENDDO
    ENDDO

    DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CR3X(I,K,J)=INT_JACOB(K,J)*CU3X(I,K,J)
          ENDDO
        ENDDO
    ENDDO


! C Add the R-K term from the rk-1 step

    IF (RK_STEP .GT. 1) THEN
        DO J=JSTART,JEND
          DO K=ZSTART,ZEND
            DO I=0,NX2P
              CR1X(I,K,J)=CR1X(I,K,J)+TEMP3*CF1X(I,K,J)
            ENDDO
          ENDDO
        ENDDO
        DO J=JSTART,JEND
          DO K=ZSTART,ZEND
            DO I=0,NX2P
              CR2X(I,K,J)=CR2X(I,K,J)+TEMP3*CF2X(I,K,J)
            ENDDO
          ENDDO
        ENDDO
        
        DO J=JSTART,JEND
          DO K=ZSTART,ZEND
            DO I=0,NX2P
              CR3X(I,K,J)=CR3X(I,K,J)+TEMP3*CF3X(I,K,J)
            ENDDO
          ENDDO
        ENDDO
    ENDIF

      
! c       DO J=JSTART,JEND
! c         DO K=ZSTART,ZEND 
! c          CR1X(0,K,J)=CR1X(0,K,J)-TEMP4*INT_JACOB(K,J)*PX0
! c	  CR2X(0,K,J)=CR2X(0,K,J)-TEMP4*INT_JACOB(K,J)*PX0/2.0
! C	  CR3X(0,K,J)=CR3X(0,K,J)-TEMP4*INT_JACOB(K,J)*PX0
! c         ENDDO
! c        ENDDO   
! c      IF (F_TYPE == 1) THEN
! C Add the pressure gradient to the RHS as explicit Euler	
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NX2P  
	  ! C  dp/dz => d[J-1*\zeta_x*p]/d\zeta + d[J-1*\eta_x*p]/d\eta	 
	    CS1X(I,K,J) = 0.5d0 * CJOB_11(K+1,J,2) * ( CPX(I,K,J) + CPX(I,K+1,J) ) & 
			- 0.5d0 * CJOB_11(K,J,2) * ( CPX(I,K,J) + CPX(I,K-1,J) )  &
			+ 0.5d0 * CJOB_21(K,J+1,1) * ( CPX(I,K,J) + CPX(I,K,J+1) )  &
			- 0.5d0 * CJOB_21(K,J,1) * ( CPX(I,K,J) + CPX(I,K,J-1) )  
	  
	  ! C  dp/dy => d[J-1*\zeta_y*p]/d\zeta + d[J-1*\eta_y*p]/d\eta          
	    
	    CS2X(I,K,J) = 0.5d0 * CJOB_12(K+1,J,2) * ( CPX(I,K,J) + CPX(I,K+1,J) ) &
			  - 0.5d0 * CJOB_12(K,J,2) * ( CPX(I,K,J) + CPX(I,K-1,J) ) &
			  + 0.5d0 * CJOB_22(K,J+1,1) * ( CPX(I,K,J) + CPX(I,K,J+1) ) &
			  - 0.5d0 * CJOB_22(K,J,1) * ( CPX(I,K,J) + CPX(I,K,J-1) )
	ENDDO
      ENDDO
    ENDDO	 
      
    ! C Add the pressure gradient to the RHS as explicit Euler
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NX2P
      ! C  dp/dx => iKX(I)*J-1*{P}	  
	  CR1X(I,K,J) = CR1X(I,K,J) - TEMP4 * CIKXP(I) * INT_JACOB(K,J) * CPX(I,K,J)
	ENDDO
      ENDDO
    ENDDO

    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NX2P
	  CR2X(I,K,J) = CR2X(I,K,J) - TEMP4 * CS2X(I,K,J)
	ENDDO
      ENDDO
    ENDDO

    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NX2P
	  CR3X(I,K,J) = CR3X(I,K,J) - TEMP4 * CS1X(I,K,J)
	ENDDO
      ENDDO
    ENDDO      
			
    IF (RANK .EQ. 0) THEN
    
    
     IF (F_TYPE.EQ.0) THEN
       
       Qn1=0.0d0
!       DO J=JSTART,JEND-1
!         Qn1=Qn1+(dble(CU3X(0,3,J+1))+dble(CU3X(0,3,J)))*(ypoint(1,J+1)-ypoint(1,J))/2
!       ENDDO
      
       DO J=JSTART,JEND-1
        DO K=ZSTART+1, ZEND-1
          Qn1=Qn1+(dble(CU3X(0,K,J+1))+dble(CU3X(0,K,J)))*(ypoint(K,J+1)-ypoint(K,J))/2
        ENDDO
       ENDDO
       
      ! Qn1=Qn1/(ypoint(1,JEND)-ypoint(1,JSTART))  
       Qn1=Qn1/(ypoint(1,JEND)-ypoint(1,JSTART))/float(ZEND-ZSTART-1)  
       PXV=PX0-(weight/TEMP4)*(Qn+Q0-2*Qn1)
       PX0=PXV
       Qn=Qn1

       DO J=JSTART,JEND
         DO K=ZSTART,ZEND
             CR3X(0,K,J)=CR3X(0,K,J)-TEMP4*PX0*INT_JACOB(K,J)
         END DO
       END DO 

       ELSEIF ( F_TYPE .EQ. 1) THEN
	DO J=JSTART,JEND
	  DO K=ZSTART,ZEND
	      CR3X(0,K,J) = CR3X(0,K,J) - TEMP4 * PX0 * INT_JACOB(K,J)
	  ENDDO
	ENDDO
      
      ELSE IF (F_TYPE.EQ.2) THEN
      ! C If oscillatory pressure gradient
	DO J=JSTART,JEND
	  DO K=ZSTART,ZEND
	    CR3X(0,K,J) = CR3X(0,K,J) - TEMP4 * (PX0 + AMP_OMEGA0  &
			 * cos(OMEGA0*TIME)) * INT_JACOB(K,J)
	  ENDDO
	ENDDO 
      ENDIF
    ENDIF 

    ! C Now compute the term R-K term Ai
    ! C Compile terms of Ai in CFi which will be saved for next time step
    ! C First, store the horizontal viscous terms in CFi

    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NX2P
	
	  CF1X(I,K,J) = -NU * KX2P(I) * INT_JACOB(K,J) * CU1X(I,K,J)
	  CF2X(I,K,J) = -NU * KX2P(I) * INT_JACOB(K,J) * CU2X(I,K,J) 
	  CF3X(I,K,J) = -NU * KX2P(I) * INT_JACOB(K,J) * CU3X(I,K,J)
	  
	  !  Adding rotation when velocities are in fourier space
	  
	  CF3X(I,K,J)=CF3X(I,K,J) + F_0 * (CU1X(I,K,J) - U_BC_YMAX_C1)*INT_JACOB(K,J)
!         CF1X(I,K,J)=CF1X(I,K,J) - F_0 * (CU3X(I,K,J) - W_BC_YMAX_C1)*INT_JACOB(K,J)
          CF1X(I,K,J)=CF1X(I,K,J) - F_0 * (CU3X(I,K,J) - WBULK0)*INT_JACOB(K,J)

	  
	ENDDO 
      ENDDO
    ENDDO

  ! Do for each scalar
    DO N=1,N_TH
  ! If a scalar contributes to the denisty, RI_TAU is not equal to zero and
  ! add the buoyancy term as explicit R-K.  Don't add the 0,0 mode which 
  ! corresponds to a plane average.  The plane averaged density balances
  ! the hydrostratic pressure component.
      IF ((F_TYPE .EQ. 4).OR.(F_TYPE.EQ.5)) THEN        
	    DO J=2,NY 
	    DO K=ZSTART,ZEND
	      DO I=0,NX2P
    ! Use second order interpolation
		  CF2X(I,K,J)=CF2X(I,K,J) - RI_TAU(N)*cos(ANG_BETA)* &
		(CTHX(I,K,J,N)*DYF(J-1)+CTHX(I,K,J-1,N)*DYF(J))/(2.d0*DY(J))
	      ENDDO
	    ENDDO
	    ENDDO        
	  
	    DO J=JSTART,JEND
	    DO K=ZSTART,ZEND
	      DO I=0,NX2P
		CF1X(I,K,J)=CF1X(I,K,J)-RI_TAU(N)*sin(ANG_BETA)*CTHX(I,K,J,N)
	      ENDDO
	    ENDDO
	    ENDDO

      ELSE


!        IF ( TH_BC_ZMAX(N) .EQ. 6 ) THEN
!	  DO J=JSTART,JEND
!	   DO K=ZSTART-1,ZEND+1
!	     DO I=1,NX2P
!             CF2X(I,K,J) = CF2X(I,K,J) - INT_JACOB(K,J) * RI_TAU(N) *CTHX(I,K,J,N)
!	     ENDDO
!	   ENDDO
!	  ENDDO
!        ELSE
          DO J=JSTART,JEND
           DO K=ZSTART,ZEND
             DO I=0,NX2P
               CF2X(I,K,J) = CF2X(I,K,J) - INT_JACOB(K,J) * RI_TAU(N)*CTHX(I,K,J,N) 
             ENDDO
           ENDDO
          ENDDO
!        ENDIF

      ENDIF
	
      
! Now, compute the RHS vector for the scalar equations
! Since TH is defined at horizontal velocity points, the
! scalar update equation will be very similar to the horizontal
! velocity update.

! We will store the RHS scalar terms in CRTH, RTHX
! The k-1 term for the R-K stepping is saved in FTHX, CFTHX

      CRTHX(:,:,:,N) = (0.d0,0.d0)

! First, build the RHS vector, use CRTH
      DO J=JSTART_TH(N),JEND_TH(N)
	DO K=ZSTART_TH(N),ZEND_TH(N)
	  DO I=0,NX2P
	    CRTHX(I,K,J,N) = INT_JACOB(K,J) * CTHX(I,K,J,N)
	  ENDDO
	ENDDO
      ENDDO
      
! Add term from k-2 step to free up CFTHX variable
      IF (RK_STEP .GT. 1) THEN
	DO J=JSTART_TH(N),JEND_TH(N)
	  DO K=ZSTART_TH(N),ZEND_TH(N)
	    DO I=0,NX2P
		CRTHX(I,K,J,N) = CRTHX(I,K,J,N) + TEMP3 * CFTHX(I,K,J,N)
	    ENDDO
	  ENDDO
	ENDDO
      ENDIF
   
! c       DO J=JSTART_TH(N),JEND_TH(N)
! c         DO K=ZSTART_TH(N),ZEND_TH(N) 
! c          CRTH(0,K,J,N)=CRTH(0,K,J,N)-TEMP4*INT_JACOB(K,J)*PX0
! c         ENDDO
! c        ENDDO
       
       
! Now compute the explicit R-K term Ai
! Compile terms of Ai in CFi which will be saved for next time step
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NX2P
           CFTHX(I,K,J,N) = -(NU/PR(N)) * INT_JACOB(K,J) * KX2P(I) * CTHX(I,K,J,N)
          ENDDO
        ENDDO
      ENDDO
      
! ! c      IF (F_TYPE.EQ.4) THEN
!        DO J=JSTART_TH(N),JEND_TH(N)
!         DO K=ZSTART_TH(N),ZEND_TH(N)
!           DO I=0,NX2P
!             CFTHX(I,K,J,N)=CFTHX(I,K,J,N) + INT_JACOB(K,J)*  &
!                     CU2X(I,K,J)!/(1.d0 + SPONGE_SIGMA_OUT(K))
!           ENDDO
!         ENDDO
!       ENDDO
! ! End of  loop for passive scalars (N_TH)

    IF ( CONT_STRAT ) THEN
	DO J=JSTART_TH(N),JEND_TH(N)
	  DO K=ZSTART_TH(N),ZEND_TH(N)
	    DO I=0,NX2P
	      CFTHX(I,K,J,N) = CFTHX(I,K,J,N) + INT_JACOB(K,J) * CU2X(I,K,J)!/(1.d0 + SPONGE_SIGMA_OUT(K))
	    ENDDO
	  ENDDO
	ENDDO
    ELSE
      DO J=JSTART_TH(N),JEND_TH(N)
	DO K=ZSTART_TH(N),ZEND_TH(N)
	  DO I=0,NX2P
	    CFTHX(I,K,J,N) = CFTHX(I,K,J,N) - CU2X(I,K,J) 	& 
		  * (  &
		      0.5 * CJOB_12(K+1,J,2) * ( TH_BAR(K,J) + TH_BAR(K+1,J) ) &
		    - 0.5 * CJOB_12(K,J,2) * ( TH_BAR(K,J) + TH_BAR(K-1,J) )   &
		    + 0.5 * CJOB_22(K,J+1,1) * ( TH_BAR(K,J) + TH_BAR(K,J+1) ) &
		    - 0.5 * CJOB_22(K,J,1) * ( TH_BAR(K,J) + TH_BAR(K,J-1) )   & 
		    )

    !!!!   Addition of kappa*d^2(rho_bar(x,z))/dz^2 + kappa*d^2(rho_bar(x,z))/dx^2
    !!!!   at right hand side of density equation

              IF (LES) THEN
!in case of LES, Kappa_t*d^2(th)/dz^2 (eddy diffusivity term of background) should be included 
               CFTHX(I,K,J,N) = CFTHX(I,K,J,N) +  0.25d0 * ( NU/PR(N) +KAPPA_T(I,K,J,N)) * (                                           &
                + GMAT_12(K+1,J,2) * (TH_BAR(K+1,J+1) + TH_BAR(K,J+1)-TH_BAR(K,J-1) - TH_BAR(K+1,J-1) )                               &
                - GMAT_12(K,J,2) * (TH_BAR(K-1,J+1) + TH_BAR(K,J+1)-TH_BAR(K,J-1) - TH_BAR(K-1,J-1))                                  &
                + GMAT_12(K,J+1,1) * ( TH_BAR(K+1,J+1) + TH_BAR(K+1,J)-TH_BAR(K-1,J+1) - TH_BAR(K-1,J) )                              &
                - GMAT_12(K,J,1) * ( TH_BAR(K+1,J) + TH_BAR(K+1,J-1)-TH_BAR(K-1,J) - TH_BAR(K-1,J-1) )                                &
                + 4.0d0 * ( GMAT_11(K+1,J,2) * ( TH_BAR(K+1,J) - TH_BAR(K,J) )-GMAT_11(K,J,2) * ( TH_BAR(K,J) - TH_BAR(K-1,J)) )      &
                + 4.0d0 * ( GMAT_22(K,J+1,1)*(TH_BAR(K,J+1) - TH_BAR(K,J))-GMAT_22(K,J,1) * ( TH_BAR(K,J) - TH_BAR(K,J-1)) )          &
                  )
              ELSE
	       CFTHX(I,K,J,N) = CFTHX(I,K,J,N) +  0.25d0 * ( NU/PR(N) ) * (								&
		+ GMAT_12(K+1,J,2) * (TH_BAR(K+1,J+1) + TH_BAR(K,J+1) - TH_BAR(K,J-1) - TH_BAR(K+1,J-1) )				&
		- GMAT_12(K,J,2) * (TH_BAR(K-1,J+1) + TH_BAR(K,J+1) - TH_BAR(K,J-1) - TH_BAR(K-1,J-1))				& 
		+ GMAT_12(K,J+1,1) * ( TH_BAR(K+1,J+1) + TH_BAR(K+1,J) - TH_BAR(K-1,J+1) - TH_BAR(K-1,J) )				&
		- GMAT_12(K,J,1) * ( TH_BAR(K+1,J) + TH_BAR(K+1,J-1) - TH_BAR(K-1,J) - TH_BAR(K-1,J-1) )				&
		+ 4.0d0 * ( GMAT_11(K+1,J,2) * ( TH_BAR(K+1,J) - TH_BAR(K,J) ) - GMAT_11(K,J,2) * ( TH_BAR(K,J) - TH_BAR(K-1,J)) ) 	&
		+ 4.0d0 * ( GMAT_22(K,J+1,1)*(TH_BAR(K,J+1) - TH_BAR(K,J)) - GMAT_22(K,J,1) * ( TH_BAR(K,J) - TH_BAR(K,J-1)) )	& 
		  )
              ENDIF
	  ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF (SPONGE_TYPE .NE. 0) CALL sponge_th(N)   
  ENDDO
            
! c      CALL sponge_vel
      
      
  IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
  
     CALL les_chan
     
     !CALL REAL_FOURIER_TRANS_U1 (.false.)
     !CALL REAL_FOURIER_TRANS_U2 (.false.)
     !CALL REAL_FOURIER_TRANS_U3 (.false.)
     
     
! Add the subgrid scale scalar flux to the scalar equations
      DO N=1,N_TH
	CALL les_chan_th(N)
!         CALL FFT_X_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N),0,NY+1,0,NZ+1)
      ENDDO

    ELSE 
! C If the subgrid model hasn't been called, then it is necessary to 
! C convert to physical space.
      CALL REAL_FOURIER_TRANS_U1 (.false.)
      CALL REAL_FOURIER_TRANS_U2 (.false.)
      CALL REAL_FOURIER_TRANS_U3 (.false.)

!C        CALL FFT_X_TO_PHYSICAL(CU1,U1,0,NY+1,0,NZ+1)
!C        CALL FFT_X_TO_PHYSICAL(CU2,U2,0,NY+1,0,NZ+1)
!C        CALL FFT_X_TO_PHYSICAL(CU3,U3,0,NY+1,0,NZ+1)

! Transform THETA to physical space for computation of nonlinear terms
! Here pass the first location in memory of the array for scalar n
      DO N=1,N_TH
       CALL REAL_FOURIER_TRANS_TH (.false.)
      ENDDO
	
    ENDIF
      
      
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CCCCC        Non-linear terms calculation   CCCCC
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C actual nonlinear terms
! c      IF ( F_TYPE == 4 ) THEN
! C     For x-momentum equation
! C d(U1*u1)/dx
    S1X(:,:,:) =0.d0

    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NXP
	  S1X(I,K,J) = INT_JACOB(K,J) * U1X(I,K,J) * U1X(I,K,J)
	ENDDO
      ENDDO
    ENDDO

    CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
    varp(:,:,:) = 0.d0

    DO I=0,NXM
      varp(I,:,:)=S1Z(I,:,:)
    ENDDO

    CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)

    DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
    ENDDO

    CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
    CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)

!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)
      
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NX2P
	  CF1X(I,K,J) = CF1X(I,K,J) - CIKXP(I) * CS1X(I,K,J) 
	ENDDO
      ENDDO
    ENDDO

! C Calculation of Contravariant velocity

    DO J=JSTART,JEND+1
      DO K=ZSTART,ZEND+1
	DO I=0,NXP
  ! C      U2 required to define at bottom cell face 	  
	  U2bX(I,K,J) = 0.5d0 * CJOB_22(K,J,1) * ( U2X(I,K,J)+U2X(I,K,J-1) )  &
		      + 0.5d0 * CJOB_21(K,J,1) * ( U3X(I,K,J)+U3X(I,K,J-1) )
  ! C      U3 required to define at side cell face     	    
	  U3bX(I,K,J) = 0.5d0 * CJOB_11(K,J,2) * ( U3X(I,K,J)+U3X(I,K-1,J) ) &
		      + 0.5d0 * CJOB_12(K,J,2) * ( U2X(I,K,J)+U2X(I,K-1,J) )
	ENDDO
      ENDDO
    ENDDO

  ! C d(U2*u1)/d\eta and d(U3*u1)/d\zeta together       
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NXP
  ! C d(U2*u1)/d\eta	  
	  S1X(I,K,J) = 0.5d0 * U2bX(I,K,J+1) * ( U1X(I,K,J)+U1X(I,K,J+1) ) &
		     - 0.5d0 * U2bX(I,K,J) * ( U1X(I,K,J)+U1X(I,K,J-1) )
  ! C d(U3*u1)/d\zeta   
	  S2X(I,K,J) = 0.5d0 * U3bX(I,K+1,J) * ( U1X(I,K,J)+U1X(I,K+1,J) ) &
		     - 0.5d0 * U3bX(I,K,J) * ( U1X(I,K,J)+U1X(I,K-1,J) )              
    

! ! No-linear terms for U1 comes from the diffusion terms
! 
! ! Kappa*d/d\zeta(GMAT_12(:,:,2)dU1/d\eta

	S1X(I,K,J) = -S1X(I,K,J) + 0.25d0 * NU * &  
		     ( GMAT_12(K+1,J,2)* ( U1X(I,K+1,J+1) + U1X(I,K,J+1) - U1X(I,K,J-1) - U1X(I,K+1,J-1) )  &
		     - GMAT_12(K,J,2)*(U1X(I,K-1,J+1) + U1X(I,K,J+1) - U1X(I,K,J-1) - U1X(I,K-1,J-1))  & 
		     )
! Kappa*d/d\eta(GMAT_12(:,:,1))dU1/d\zeta
 	S2X(I,K,J) = -S2X(I,K,J) + 0.25d0 * NU * & 
		    ( GMAT_12(K,J+1,1) * ( U1X(I,K+1,J+1) + U1X(I,K+1,J) - U1X(I,K-1,J+1) - U1X(I,K-1,J) ) &
		    - GMAT_12(K,J,1) * ( U1X(I,K+1,J) + U1X(I,K+1,J-1) - U1X(I,K-1,J) - U1X(I,K-1,J-1) ) &
		    )
        ENDDO
      ENDDO
    ENDDO

    S1X(:,:,:)=S1X(:,:,:)+S2X(:,:,:)
  
    CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
    varp(:,:,:) = 0.d0
  
    DO I=0,NXM
      varp(I,:,:) = S1Z(I,:,:)
    ENDDO
    
    CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
    
    DO I=0,NKX
      CS1Z(I,:,:) = cvarp(I,:,:)
    ENDDO
    
    CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
    CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)


     DO J=JSTART,JEND
       DO K=ZSTART,ZEND
         DO I=0,NX2P
           CF1X(I,K,J) = CF1X(I,K,J) + CS1X(I,K,J)
         ENDDO
       ENDDO
     ENDDO

     
  ! C     For y-momentum equation
  ! C d(U1*u2)/dx

    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NXP
	  S1X(I,K,J) = INT_JACOB(K,J) * U1X(I,K,J) * U2X(I,K,J)
	ENDDO
      ENDDO
    ENDDO

    CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
    varp(:,:,:) = 0.d0
    
    DO I=0,NXM
      varp(I,:,:)=S1Z(I,:,:)
    ENDDO
    
    CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
    
    DO I=0,NKX
      CS1Z(I,:,:) = cvarp(I,:,:)
    ENDDO
    
    CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
    CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)

    
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NX2P
	  CF2X(I,K,J) = CF2X(I,K,J) - CIKXP(I) * CS1X(I,K,J) 
	ENDDO
      ENDDO
    ENDDO

  ! C d(U2*u2)/d\eta and d(U3*u2)/d\zeta together 
  ! 
  ! C Now at this time U2 at cell bottom and U3 at cell
  ! C    side face have been calculated previous step.  
    
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NXP
! C d(U2*u2)/d\eta	  
	  S1X(I,K,J) = 0.5d0 * U2bX(I,K,J+1) * ( U2X(I,K,J)+U2X(I,K,J+1) ) &
		     - 0.5d0 * U2bX(I,K,J) * ( U2X(I,K,J)+U2X(I,K,J-1) )
! C d(U3*u2)/d\zeta   
	  S2X(I,K,J) = 0.5 * U3bX(I,K+1,J) * ( U2X(I,K,J)+U2X(I,K+1,J) ) &
		     - 0.5 * U3bX(I,K,J) * ( U2X(I,K,J)+U2X(I,K-1,J) )   

! No-linear terms for U2 comes from the diffusion terms

! Kappa*d/d\zeta(GMAT_12(:,:,2)dU2/d\eta

	S1X(I,K,J) = -S1X(I,K,J) + 0.25d0 * NU * & 
		    ( GMAT_12(K+1,J,2) * ( U2X(I,K+1,J+1) + U2X(I,K,J+1) - U2X(I,K,J-1) - U2X(I,K+1,J-1) )   &
		    - GMAT_12(K,J,2) * ( U2X(I,K-1,J+1) + U2X(I,K,J+1) - U2X(I,K,J-1) - U2X(I,K-1,J-1)) & 
		    )
! Kappa*d/d\eta(GMAT_12(:,:,1))dU2/d\zeta

	S2X(I,K,J) = -S2X(I,K,J) + 0.25d0 * NU * & 
		    ( GMAT_12(K,J+1,1) * ( U2X(I,K+1,J+1) + U2X(I,K+1,J) - U2X(I,K-1,J+1) - U2X(I,K-1,J) ) &
		    - GMAT_12(K,J,1) * ( U2X(I,K+1,J) + U2X(I,K+1,J-1)   - U2X(I,K-1,J) - U2X(I,K-1,J-1)) & 
		    )
           
	ENDDO
      ENDDO
    ENDDO

    S1X(:,:,:)=S1X(:,:,:)+S2X(:,:,:)
    CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
    
    varp(:,:,:) = 0.d0
    
    DO I=0,NXM
      varp(I,:,:)=S1Z(I,:,:)
    ENDDO
    
    CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
    
    DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
    ENDDO
    
    CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
    
    CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)

	
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NX2P
	  CF2X(I,K,J) = CF2X(I,K,J) + CS1X(I,K,J) 
	ENDDO
      ENDDO
    ENDDO      

  ! C     For z-momentum equation
  ! C d(U1*u3)/dx

    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NXP
	  S1X(I,K,J) = INT_JACOB(K,J) * U1X(I,K,J) * U3X(I,K,J)
	ENDDO
      ENDDO
    ENDDO
    
    CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
    varp(:,:,:) = 0.d0
    
    DO I=0,NXM
      varp(I,:,:)=S1Z(I,:,:)
    ENDDO
    
    CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
    
    DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
    ENDDO
    
    CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
    CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)

  !C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)
    
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NX2P
	  CF3X(I,K,J) = CF3X(I,K,J) - CIKXP(I) * CS1X(I,K,J) 
	ENDDO
      ENDDO
    ENDDO

  ! C d(U2*u3)/d\eta and d(U3*u3)/d\zeta together 
  ! 
  ! C Now at this time U2 at cell bottom and U3 at cell
  ! C    side face have been calculated previous step.  
    
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NXP
  ! C d(U2*u3)/d\eta	  
	  S1X(I,K,J) = 0.5d0 * U2bX(I,K,J+1) * ( U3X(I,K,J) + U3X(I,K,J+1) ) &
		     - 0.5d0 * U2bX(I,K,J) * ( U3X(I,K,J) + U3X(I,K,J-1) )
  ! C d(U3*u3)/d\zeta   
	  S2X(I,K,J) = 0.5d0 * U3bX(I,K+1,J) * ( U3X(I,K,J)+U3X(I,K+1,J) ) &
		     - 0.5d0 * U3bX(I,K,J) * ( U3X(I,K,J)+U3X(I,K-1,J) )      

! No-linear terms for U3 comes from the diffusion terms
! Kappa*d/d\zeta(GMAT_12(:,:,2)dU3/d\eta

	S1X(I,K,J)= -S1X(I,K,J) + 0.25d0 * NU * & 
		    ( GMAT_12(K+1,J,2) * ( U3X(I,K+1,J+1) + U3X(I,K,J+1) - U3X(I,K,J-1) - U3X(I,K+1,J-1) )  &
		    - GMAT_12(K,J,2) * ( U3X(I,K-1,J+1) + U3X(I,K,J+1) - U3X(I,K,J-1) - U3X(I,K-1,J-1)) & 
		    )
      
! Kappa*d/d\eta(GMAT_12(:,:,1))dU3/d\zeta

	S2X(I,K,J)= -S2X(I,K,J) + 0.25d0 * NU * & 
		    ( GMAT_12(K,J+1,1) * ( U3X(I,K+1,J+1) + U3X(I,K+1,J) - U3X(I,K-1,J+1) - U3X(I,K-1,J) )  &
		    - GMAT_12(K,J,1) * ( U3X(I,K+1,J) + U3X(I,K+1,J-1) - U3X(I,K-1,J) - U3X(I,K-1,J-1)) & 
		    )
        
	ENDDO
      ENDDO
    ENDDO

    S1X(:,:,:)=S1X(:,:,:)+S2X(:,:,:)
    CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
    
    varp(:,:,:) = 0.d0
    
    DO I=0,NXM
      varp(I,:,:) = S1Z(I,:,:)
    ENDDO
    
    CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
    
    DO I=0,NKX
      CS1Z(I,:,:) = cvarp(I,:,:)
    ENDDO
    
    CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
    
    CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)

	      
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NX2P
	  CF3X(I,K,J)=CF3X(I,K,J) + CS1X(I,K,J)
	ENDDO
      ENDDO
    ENDDO      
  
! c      ENDIF      

  IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
	DO J=JSTART,JEND
	  DO K=ZSTART,ZEND
	    DO I=0,NXP
	    S1X(I,K,J)=  0.25d0 * ( &  
				    GMAT_12(K+1,J,2) * 0.5d0 * ( NU_T(I,K+1,J)+NU_T(I,K,J) )			&
				    * ( U1X(I,K+1,J+1) + U1X(I,K,J+1) - U1X(I,K,J-1) - U1X(I,K+1,J-1) )	&
				    - GMAT_12(K,J,2) * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K-1,J) )			&
				    * ( U1X(I,K-1,J+1) + U1X(I,K,J+1) - U1X(I,K,J-1) - U1X(I,K-1,J-1))	& 
				    
				    + GMAT_12(K,J+1,1) * 0.5d0 * (NU_T(I,K,J)+NU_T(I,K,J+1))			&
				    * ( U1X(I,K+1,J+1) + U1X(I,K+1,J) - U1X(I,K-1,J+1) - U1X(I,K-1,J))	&
				    - GMAT_12(K,J,1) * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K,J-1) )			&
				    * ( U1X(I,K+1,J) + U1X(I,K+1,J-1) - U1X(I,K-1,J) - U1X(I,K-1,J-1) )	&
				    )
	    ENDDO  
	   ENDDO
	  ENDDO

	CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
	varp(:,:,:) = 0.d0
	
	DO I=0,NXM
	varp(I,:,:)=S1Z(I,:,:)
	ENDDO
	
	CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
	
	DO I=0,NKX
	CS1Z(I,:,:)=cvarp(I,:,:)
	ENDDO
	
	CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
	CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)

  !C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)

	DO J=JSTART,JEND
	  DO K=ZSTART,ZEND
	    DO I=0,NX2P
	      CF1X(I,K,J)=CF1X(I,K,J) + CS1X(I,K,J)
	    ENDDO
	  ENDDO
	ENDDO

	DO J=JSTART,JEND
	  DO K=ZSTART,ZEND
	    DO I=0,NXP
	      S1X(I,K,J)=  0.25d0 * ( & 
				    GMAT_12_y(K+1,J,2) * 0.5d0 * ( NU_T(I,K+1,J) + NU_T(I,K,J) )		&
				    * ( U2X(I,K+1,J+1) + U2X(I,K,J+1) - U2X(I,K,J-1) - U2X(I,K+1,J-1) )	&
				    - GMAT_12_y(K,J,2) * 0.5d0 * ( NU_T(I,K,J) + NU_T(I,K-1,J) )		&
				    * ( U2X(I,K-1,J+1) + U2X(I,K,J+1) - U2X(I,K,J-1) - U2X(I,K-1,J-1) ) 	& 
				    
				    + GMAT_12_y(K,J+1,1) * 0.5d0 * ( NU_T(I,K,J) + NU_T(I,K,J+1) )		&
				    * ( U2X(I,K+1,J+1) + U2X(I,K+1,J) - U2X(I,K-1,J+1) - U2X(I,K-1,J) )	&
				    - GMAT_12_y(K,J,1) * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K,J-1) )		&
				    * ( U2X(I,K+1,J) + U2X(I,K+1,J-1) - U2X(I,K-1,J) - U2X(I,K-1,J-1) )	& 
				    )

	    ENDDO
	  ENDDO
	ENDDO

	CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
	varp(:,:,:) = 0.d0
	DO I=0,NXM
	varp(I,:,:)=S1Z(I,:,:)
	ENDDO
	CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
	DO I=0,NKX
	CS1Z(I,:,:)=cvarp(I,:,:)
	ENDDO
	CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
	CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)

  !C       CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)

	DO J=JSTART,JEND
	  DO K=ZSTART,ZEND
	    DO I=0,NX2P
	      CF2X(I,K,J)=CF2X(I,K,J) + CS1X(I,K,J)
	    ENDDO
	  ENDDO
	ENDDO
		  
	DO J=JSTART,JEND
	  DO K=ZSTART,ZEND
	    DO I=0,NXP
	      S1X(I,K,J)=  0.25d0 * ( & 
			    GMAT_12_z(K+1,J,2) * 0.5d0 * ( NU_T(I,K+1,J) + NU_T(I,K,J) )		&
			    * ( U3X(I,K+1,J+1) + U3X(I,K,J+1) - U3X(I,K,J-1) - U3X(I,K+1,J-1) )	&
			    - GMAT_12_z(K,J,2) * 0.5d0 * ( NU_T(I,K,J) + NU_T(I,K-1,J) )		&
			    * ( U3X(I,K-1,J+1) + U3X(I,K,J+1) - U3X(I,K,J-1) - U3X(I,K-1,J-1) ) 	& 
			    
			    + GMAT_12_z(K,J+1,1) * 0.5d0 * ( NU_T(I,K,J) + NU_T(I,K,J+1) )		&
			    * ( U3X(I,K+1,J+1) + U3X(I,K+1,J) - U3X(I,K-1,J+1) - U3X(I,K-1,J) )	&
			    - GMAT_12_z(K,J,1) * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K,J-1) )		&
			    * ( U3X(I,K+1,J) + U3X(I,K+1,J-1) - U3X(I,K-1,J) - U3X(I,K-1,J-1) )	& 
			    )
	      ENDDO
	    ENDDO
	  ENDDO

	CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
	varp(:,:,:) = 0.d0
	
	DO I=0,NXM
	varp(I,:,:)=S1Z(I,:,:)
	ENDDO
	
	CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
	
	DO I=0,NKX
	CS1Z(I,:,:)=cvarp(I,:,:)
	ENDDO
	
	CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
	CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)

  !C       CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)

	DO J=JSTART,JEND
	  DO K=ZSTART,ZEND
	    DO I=0,NX2P
	      CF3X(I,K,J)=CF3X(I,K,J) + CS1X(I,K,J)
	    ENDDO
	  ENDDO
	ENDDO

    ENDIF
      
! C -- At this point, we are done computing the nonlinear terms --
! 
!       
! C Finally, Add CFi to CRi

  DO J=JSTART,JEND
    DO K=ZSTART,ZEND
      DO I=0,NX2P
	CR1X(I,K,J) = CR1X(I,K,J) + TEMP5 * CF1X(I,K,J)
      END DO
    END DO
  END DO
  
  DO J=JSTART,JEND
    DO K=ZSTART,ZEND
      DO I=0,NX2P
	CR2X(I,K,J) = CR2X(I,K,J) + TEMP5 * CF2X(I,K,J)
      ENDDO
    ENDDO
  ENDDO

  DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NX2P
	  CR3X(I,K,J) = CR3X(I,K,J) + TEMP5 * CF3X(I,K,J)
	ENDDO
      ENDDO
    ENDDO

! C Convert RHS terms to physical space
! C made a change in FFT_X_TO_PHYSICAL(CR2X,R2X,0,NY+1,0,NZ+1)
! C
! 
! c      IF (MOD(TIME_STEP,5).EQ.0) THEN 
! c       CALL sponge_vel
! c      endif

  CALL REAL_FOURIER_TRANS_R1 (.false.)
  CALL REAL_FOURIER_TRANS_R2 (.false.)
  CALL REAL_FOURIER_TRANS_R3 (.false.)

!      CALL FFT_X_TO_PHYSICAL(CR1X,R1X,0,NY+1,0,NZ+1)                 
!      CALL FFT_X_TO_PHYSICAL(CR2X,R2X,0,NY+1,0,NZ+1)                 
!      CALL FFT_X_TO_PHYSICAL(CR3X,R3X,0,NY+1,0,NZ+1)      
      
      
          
  DO N=1,N_TH
! ! Do for each scalar:
! 
! c      IF ( F_TYPE .EQ. 4 ) THEN
! ! Compute the nonlinear terms that are present in the explicit term A
! U1*TH
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            S1X(I,K,J)=INT_JACOB(K,J)*THX(I,K,J,N)*U1X(I,K,J)
          ENDDO
        ENDDO
      ENDDO

      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      
      DO I=0,NXM
       varp(I,:,:)=S1Z(I,:,:)
      ENDDO
      
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      
      DO I=0,NKX
       CS1Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)

!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)

      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NX2P
            CFTHX(I,K,J,N)=CFTHX(I,K,J,N) - CIKXP(I) * CS1X(I,K,J)
          ENDDO
        ENDDO
      ENDDO
! 
! ! U3*TH and U2*TH together
! C d(U2*TH)/d\eta and d(U3*TH)/d\zeta together       
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
! C d(U2*TH)/d\eta	  
            S1X(I,K,J) = 0.5d0 * U2bX(I,K,J+1) * ( THX(I,K,J,N) + THX(I,K,J+1,N) ) &
		        - 0.5d0 * U2bX(I,K,J) * ( THX(I,K,J,N) + THX(I,K,J-1,N) )
! C d(U3*TH)/d\zeta   
	    S2X(I,K,J)= 0.5d0 * U3bX(I,K+1,J) * ( THX(I,K,J,N)+THX(I,K+1,J,N) ) &
		      - 0.5d0 * U3bX(I,K,J) * ( THX(I,K,J,N)+THX(I,K-1,J,N) )     

! No-linear terms comes from the diffusion terms

! Kappa*d/d\zeta(GMAT_12(:,:,2)dTH/d\eta
	    S1X(I,K,J)= -S1X(I,K,J) + 0.25d0 * ( NU/PR(N) ) *  & 
			( GMAT_12(K+1,J,2) * (THX(I,K+1,J+1,N) + THX(I,K,J+1,N)- THX(I,K,J-1,N) - THX(I,K+1,J-1,N))  &
			- GMAT_12(K,J,2) * (THX(I,K-1,J+1,N) + THX(I,K,J+1,N) - THX(I,K,J-1,N) - THX(I,K-1,J-1,N)) 	& 
			)
         
! Kappa*d/d\eta(GMAT_12(:,:,1))dTH/d\zeta

           S2X(I,K,J)= -S2X(I,K,J) + 0.25d0 * ( NU/PR(N) ) * & 
			( GMAT_12(K,J+1,1) * (THX(I,K+1,J+1,N) + THX(I,K+1,J,N) - THX(I,K-1,J+1,N) - THX(I,K-1,J,N)) &
			- GMAT_12(K,J,1) * (THX(I,K+1,J,N) + THX(I,K+1,J-1,N) - THX(I,K-1,J,N) - THX(I,K-1,J-1,N)) 	& 
			)         
          ENDDO
        ENDDO
      ENDDO
  
      S1X(:,:,:)=S1X(:,:,:)+S2X(:,:,:)
      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      
      DO I=0,NXM
       varp(I,:,:) = S1Z(I,:,:)
      ENDDO
      
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      
      DO I=0,NKX
       CS1Z(I,:,:) = cvarp(I,:,:)
      ENDDO
      
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)

    
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NX2P
            CFTHX(I,K,J,N) = CFTHX(I,K,J,N) + CS1X(I,K,J) 
          ENDDO
        ENDDO
      ENDDO
! Done with non-linear terms for theta
  
    IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
! adding two reminding LES terms from les_chan_th sobrutine
	DO J=JSTART_TH(N),JEND_TH(N)
	  DO K=ZSTART_TH(N),ZEND_TH(N)
	    DO I=0,NXP
	    S1X(I,K,J)= 0.25d0 * ( 		& 
				  + GMAT_12(K+1,J,2) * 0.5d0 * ( KAPPA_T(I,K+1,J,N)+KAPPA_T(I,K,J,N) )		&
				  * ( THX(I,K+1,J+1,N) + THX(I,K,J+1,N)- THX(I,K,J-1,N) - THX(I,K+1,J-1,N) )		&
				  - GMAT_12(K,J,2) * 0.5d0 * ( KAPPA_T(I,K,J,N)+KAPPA_T(I,K-1,J,N) )			&
				  * ( THX(I,K-1,J+1,N) + THX(I,K,J+1,N) - THX(I,K,J-1,N) - THX(I,K-1,J-1,N))		& 
				  
				  + GMAT_12(K,J+1,1) * 0.5d0 * ( KAPPA_T(I,K,J,N)+KAPPA_T(I,K,J+1,N))			&
				  * ( THX(I,K+1,J+1,N) + THX(I,K+1,J,N) - THX(I,K-1,J+1,N) - THX(I,K-1,J,N))		&
				  - GMAT_12(K,J,1) * 0.5d0 * ( KAPPA_T(I,K,J,N)+KAPPA_T(I,K,J-1,N) )			&
				  * ( THX(I,K+1,J,N) + THX(I,K+1,J-1,N) - THX(I,K-1,J,N) - THX(I,K-1,J-1,N))		& 
				  )
	    ENDDO
	   ENDDO
	  ENDDO

	CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
	varp(:,:,:) = 0.d0
	
	DO I=0,NXM
	 varp(I,:,:) = S1Z(I,:,:)
	ENDDO
	
	CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
	
	DO I=0,NKX
	 CS1Z(I,:,:) = cvarp(I,:,:)
	ENDDO
	
	CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
	CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)
  

	DO J=JSTART_TH(N),JEND_TH(N)
	  DO K=ZSTART_TH(N),ZEND_TH(N)
	    DO I=0,NX2P
	      CFTHX(I,K,J,N) = CFTHX(I,K,J,N) + CS1X(I,K,J)
	    ENDDO
	  ENDDO
	ENDDO

    ENDIF


! Add CFTHX to the RHS vector CRTH
    DO J=JSTART_TH(N),JEND_TH(N)
      DO K=ZSTART_TH(N),ZEND_TH(N)
	DO I=0,NX2P
	  CRTHX(I,K,J,N) = CRTHX(I,K,J,N) + TEMP5 * CFTHX(I,K,J,N)
	ENDDO
      ENDDO
    ENDDO

! Done with computation of RHS, explicit terms for the THETA equation
! Transform back to physical space

    CALL REAL_FOURIER_TRANS_Rth (.false.)

! 
! C Before going to  ADI1 we have to store RHS after subtracting previous TH_n 
! C at previous time step to bulid new RHS for  ADI2 based on TH_n+1/2 at intermidiate 



    DO J=JSTART_TH(N),JEND_TH(N)
      DO K=ZSTART_TH(N),ZEND_TH(N)
	DO I=0,NXP
	  S1X(I,K,J)= RTHX(I,K,J,N)-INT_JACOB(K,J)*THX(I,K,J,N)
	ENDDO
      ENDDO
    ENDDO


! C Compute the vertical viscous term in physical space and add to RHS
! C  at  ADI1 step .... Important in this step we have to add TH_n
! C at previous time step so that after multiplying a factor 1/2 into the 
! C RHS side during ADI1 operation RHS will be taken care TH_n not 1/2*TH_n


    IF (WAVE_ABS) THEN
	DO J=JSTART_TH(N),JEND_TH(N)
	  DO K=ZSTART_TH(N),ZEND_TH(N)
	    DO I=0,NXP
	      RTHX(I,K,J,N) = 0.5d0 * S1X(I,K,J) + INT_JACOB(K,J) * THX(I,K,J,N) &
			    + ( TEMP1/PR(N) ) * (1.d0 + ART_VISCOSITY * SPONGE_SIGMA_OUT(K)) &
			    * ( GMAT_11(K+1,J,2) * ( THX(I,K+1,J,N) - THX(I,K,J,N) )    &
			      - GMAT_11(K,J,2) * (THX(I,K,J,N) - THX(I,K-1,J,N)) )  
	    ENDDO
	  ENDDO
	ENDDO
    ELSE 
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
	    RTHX(I,K,J,N) = 0.5d0 * S1X(I,K,J) + INT_JACOB(K,J) * THX(I,K,J,N) &
		  + ( TEMP1/PR(N) ) &
		  * ( GMAT_11(K+1,J,2) * ( THX(I,K+1,J,N) - THX(I,K,J,N) )    &
		    - GMAT_11(K,J,2) * (THX(I,K,J,N) - THX(I,K-1,J,N)) )  
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
       DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            RTHX(I,K,J,N) = RTHX(I,K,J,N) + TEMP2 * ( 		& 
			    0.5d0 * ( KAPPA_T(I,K+1,J,N)+KAPPA_T(I,K,J,N) )			&
			    * ( GMAT_11(K+1,J,2) * ( THX(I,K+1,J,N)-THX(I,K,J,N)) )            &
			    -0.5d0 * ( KAPPA_T(I,K,J,N) + KAPPA_T(I,K-1,J,N) )                &
			    * ( GMAT_11(K,J,2) * (THX(I,K,J,N) - THX(I,K-1,J,N)) ) 		& 
			    )
         ENDDO
        ENDDO
       ENDDO
    ENDIF


! c      IF ( TH_BC_ZMAX(N)  .EQ. 1 ) THEN
! c       DO J=JSTART_TH(N),JEND_TH(N)
! c	   DO I=0,NXM
! c             RTHX(I,NZ+1,J,N)=RTHX(I,NZ,J,N)
! c       	   END DO
! c       END DO  
! c      ENDIF
!  
! c      IF ( TH_BC_ZMIN(N)  .EQ. 1 ) THEN
! c       DO J=JSTART_TH(N),JEND_TH(N)
! c	   DO I=0,NXM
! c             RTHX(I,0,J,N)=RTHX(I,1,J,N)
! c       	   END DO
! c       END DO  
! c      ENDIF

      
! C Solve for TH
! C Note, here the matrix will be indexed from 1...NY+1 corresponding to U1(0:NY)

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXP
          MATLY(I,J) = 0.d0
          MATDY(I,J) = 1.d0
          MATUY(I,J) = 0.d0
          VECY(I,J) = 0.d0
        ENDDO
      ENDDO 

! C Build the implicit system of equations for U1 
!$OMP PARALLEL DO PRIVATE(K,J,I)
   DO K=ZSTART_TH(N),ZEND_TH(N)

      DO J=JSTART_TH(N),JEND_TH(N)
	DO I=0,NXP
	  MATLY(I,J) = -( TEMP1/PR(N) ) * GMAT_22(K,J,1)
	  
	  MATDY(I,J) = INT_JACOB(K,J) + ( TEMP1/PR(N) ) * ( GMAT_22(K,J+1,1) + GMAT_22(K,J,1) )   
	  
	  MATUY(I,J) = - ( TEMP1/PR(N) ) * GMAT_22(K,J+1,1)
	  
	  VECY(I,J) = RTHX(I,K,J,N)
	ENDDO
      ENDDO

! IF using a subgrid model (LES) then add the eddy viscosity part implicitly

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN

	DO J=JSTART_TH(N),JEND_TH(N)
	  DO I=0,NXP
		MATLY(I,J) = MATLY(I,J) - TEMP2 * 0.5d0 * ( KAPPA_T(I,K,J,N) + KAPPA_T(I,K,J-1,N) ) * GMAT_22(K,J,1)

		MATDY(I,J) = MATDY(I,J) + TEMP2 * 0.5d0 * ( KAPPA_T(I,K,J,N) + KAPPA_T(I,K,J+1,N) ) * GMAT_22(K,J+1,1)  &
					 + TEMP2 * 0.5d0 * ( KAPPA_T(I,K,J,N) + KAPPA_T(I,K,J-1,N) ) * GMAT_22(K,J,1)
		
		MATUY(I,J) = MATUY(I,J) - TEMP2 * 0.5d0 * ( KAPPA_T(I,K,J,N) + KAPPA_T(I,K,J+1,N) ) * GMAT_22(K,J+1,1)
          ENDDO
        ENDDO

      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C Set the boundary conditions for TH
      CALL APPLY_BC_TH_LOWER(N,K)
      CALL APPLY_BC_TH_UPPER(N,K)
! C Now, solve the tridiagonal system for TH(:,k,:)
      CALL THOMAS_REAL(MATLY,MATDY,MATUY,VECY,NY+1,NXP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO J=0,NY+1
	DO I=0,NXP
	  THX(I,K,J,N)=VECY(I,J)
	ENDDO
      ENDDO
! End do k
   ENDDO
!$OMP END PARALLEL DO

      
      
! C Compute the horizontal viscous term in physical space and add to new RHS
! C already store in S1 
! C Important we have to multiply 1/2 with S1 before adding to RHS due
! C half time splitting during AD1

    DO J=JSTART_TH(N),JEND_TH(N)
      DO K=ZSTART_TH(N),ZEND_TH(N)
	DO I=0,NXP
	  RTHX(I,K,J,N) = 0.5d0 * S1X(I,K,J) + INT_JACOB(K,J)*THX(I,K,J,N) + ( TEMP1/PR(N) ) * ( & 
		    + GMAT_22(K,J+1,1) * (THX(I,K,J+1,N) - THX(I,K,J,N))	&
		    - GMAT_22(K,J,1)*(THX(I,K,J,N)   - THX(I,K,J-1,N))  	& 
		    )

	ENDDO
      ENDDO
    ENDDO


      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
	DO J=JSTART_TH(N),JEND_TH(N)
	  DO K=ZSTART_TH(N),ZEND_TH(N)
	    DO I=0,NXP
	      RTHX(I,K,J,N) = RTHX(I,K,J,N) + TEMP2 * (  & 
			0.5d0 * ( KAPPA_T(I,K,J+1,N) + KAPPA_T(I,K,J,N) ) 		&
		      * ( GMAT_22(K,J+1,1) * (THX(I,K,J+1,N) - THX(I,K,J,N)) )        &
		      - 0.5d0 * ( KAPPA_T(I,K,J-1,N)+KAPPA_T(I,K,J,N) )               &
		      * ( GMAT_22(K,J,1) * (THX(I,K,J,N) - THX(I,K,J-1,N)) )		& 
		      )
	    ENDDO
	  ENDDO
	ENDDO
      ENDIF


! c      IF ( IC_TYPE == 7 ) THEN
! c       DO K=ZSTART_TH(N),ZEND_TH(N)
! c	   DO I=0,NXM
! c             RTHX(I,K,NY+1,N)=RTHX(I,K,NY,N)
! c       	   END DO
! c       END DO  
! c      ENDIF
! Initialize the matrix used to store implicit coefficients
      DO K=0,NZ+1
        DO I=0,NXP
          MATLX(I,K)=0.d0
          MATDX(I,K)=1.d0
          MATUX(I,K)=0.d0
          VECX(I,K)=0.d0
        ENDDO
      ENDDO 

! C Build the implicit system of equations for TH 
!$OMP PARALLEL DO PRIVATE(K,J,I)
   DO J=JSTART_TH(N),JEND_TH(N)
	IF(WAVE_ABS) THEN
	  DO K=ZSTART_TH(N),ZEND_TH(N)
	    DO I=0,NXP
		  MATLX(I,K) = -( TEMP1/PR(N) ) * GMAT_11(K,J,2)  &
			       * ( 1.d0 + ART_VISCOSITY*SPONGE_SIGMA_OUT(K) )           
		  
		  MATDX(I,K) = INT_JACOB(K,J) + (TEMP1/PR(N))     &
			      * (1.d0 + ART_VISCOSITY*SPONGE_SIGMA_OUT(K))   &
			      * (GMAT_11(K+1,J,2) + GMAT_11(K,J,2))
		  
		  MATUX(I,K)= - ( TEMP1/PR(N) ) * GMAT_11(K+1,J,2)   &
			      * (1.d0 + ART_VISCOSITY*SPONGE_SIGMA_OUT(K))
		  
		  VECX(I,K)=RTHX(I,K,J,N)
	    ENDDO
	  ENDDO
	ELSE
	  DO K=ZSTART_TH(N),ZEND_TH(N)
	    DO I=0,NXP
		MATLX(I,K) = -( TEMP1/PR(N) ) * GMAT_11(K,J,2)
      
		
		MATDX(I,K) = INT_JACOB(K,J) + (TEMP1/PR(N))     &
			    * (GMAT_11(K+1,J,2) + GMAT_11(K,J,2))
		
		MATUX(I,K)= - ( TEMP1/PR(N) ) * GMAT_11(K+1,J,2)   
		
		VECX(I,K) = RTHX(I,K,J,N)
	    ENDDO
	  ENDDO
	ENDIF 

! IF using a subgrid model (LES) then add the eddy viscosity part implicitly

        IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
          DO K=ZSTART_TH(N),ZEND_TH(N)
            DO I=0,NXP
		MATLX(I,K) = MATLX(I,K)- TEMP2 * 0.5d0 * ( KAPPA_T(I,K,J,N) + KAPPA_T(I,K-1,J,N) ) * GMAT_11(K,J,2)
		
		MATDX(I,K) = MATDX(I,K) + TEMP2 * 0.5d0 * ( KAPPA_T(I,K+1,J,N) + KAPPA_T(I,K,J,N) ) * GMAT_11(K+1,J,2)   &
					 + TEMP2 * 0.5d0 * ( KAPPA_T(I,K,J,N) + KAPPA_T(I,K-1,J,N) ) * GMAT_11(K,J,2)
		
		MATUX(I,K) = MATUX(I,K) - TEMP2 * 0.5d0 * ( KAPPA_T(I,K+1,J,N) + KAPPA_T(I,K,J,N) ) * GMAT_11(K+1,J,2)
            ENDDO
          ENDDO
        ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C Set the boundary conditions for TH
	CALL APPLY_BC_TH_LEFT(N,J)
	CALL APPLY_BC_TH_RIGHT(N,J)
! C Now, solve the tridiagonal system for TH(:,k,:)
        IF ( TH_BC_ZMAX(N)  .EQ. 5 ) THEN
         D=1.0
         CALL THOMAS_REAL_SP(MATLX,MATDX,MATUX,VECX,D,NZ,NZ+1,NXP)
        ELSE   
         CALL THOMAS_REAL(MATLX,MATDX,MATUX,VECX,NZ+1,NXP)
        ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

        DO K=0,NZ+1
          DO I=0,NXP
            THX(I,K,J,N) = VECX(I,K)
          ENDDO
        ENDDO
! End do J
   ENDDO
!$OMP END PARALLEL DO

!      IF ( TH_BC_ZMAX(N)  .EQ. 6 ) THEN
!	Do J=0,NY+1
!	  DO I=0,NXP
!	      THX(I,NZ+1,J,N)=THX(I,1,J,N)
!	  ENDDO
!	ENDDO

!	Do J=0,NY+1
!	  DO I=0,NXP
!	      THX(I,0,J,N)=THX(I,NZ,J,N)
!	  ENDDO
!	ENDDO
!      ENDIF

! C NEAD TO UPDATE THE CORNERS
      IF (TH_BC_ZMAX(N).EQ.1) THEN
	 DO I=0,NXP
          THX(I,NZ+1,NY+1,N)=THX(I,NZ,NY+1,N)
          THX(I,NZ+1,0,N)=THX(I,NZ,0,N)
          THX(I,0,NY+1,N) = THX(I,1,NY+1,N)
          THX(I,0,0,N) = THX(I,1,0,N)
       	 ENDDO
      ENDIF
 
      IF (TH_BC_ZMIN(N).EQ.6) THEN
	 DO I=0,NXP
          THX(I,NZ+1,NY+1,N)=THX(I,1,NY+1,N)
          THX(I,NZ+1,0,N)=THX(I,1,0,N)
          THX(I,0,NY+1,N) = THX(I,NZ,NY+1,N)
          THX(I,0,0,N) = THX(I,NZ,0,N)
       	 ENDDO  
      ENDIF
      
  ENDDO

! c      read(6,*)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CCCCCCCCC                   Now same for velocity          CCCCCCCCCCCCCCCCCC
! 
! 
! C ADI STEP FOR  U1 
! C Before going to  ADI1 we have to store RHS after subtracting previous U1_n 
! C at previous time step to bulid new RHS for  ADI2 based on U1_n+1/2 at intermidiate 



  DO J=JSTART,JEND
    DO K=ZSTART,ZEND
      DO I=0,NXP
	S1X(I,K,J)= R1X(I,K,J) - INT_JACOB(K,J)*U1X(I,K,J)
      ENDDO
    ENDDO
  ENDDO


! C Compute the horizontal viscous term in physical space and add to RHS
! C  at  ADI1 step .... Important in this step we have to add U1_n
! C at previous time step so that after multiplying a factor 1/2 into the 
! C RHS side during ADI1 operation RHS will be taken care U1_n not 1/2*U1_n
 

  DO J=JSTART,JEND
    DO K=ZSTART,ZEND
      DO I=0,NXP
	R1X(I,K,J) = 0.5d0 * S1X(I,K,J) + INT_JACOB(K,J) * U1X(I,K,J)   &
		    + TEMP1 * (  & 
		      GMAT_11(K+1,J,2) * (U1X(I,K+1,J) - U1X(I,K,J) ) &
		    - GMAT_11(K,J,2) * ( U1X(I,K,J)- U1X(I,K-1,J) ) & 
		    )
      ENDDO
    ENDDO
  ENDDO

  IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
    DO J=JSTART,JEND
     DO K=ZSTART,ZEND
      DO I=0,NXP
	R1X(I,K,J) = R1X(I,K,J) + TEMP2 * (  		& 
					    0.5d0 * ( NU_T(I,K+1,J)+NU_T(I,K,J) )			&
					    * ( GMAT_11(K+1,J,2) * ( U1X(I,K+1,J)- U1X(I,K,J) ) )	&
					    - 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K-1,J) )			&
					    * ( GMAT_11(K,J,2) * ( U1X(I,K,J)-U1X(I,K-1,J)) )		& 
					    )
      ENDDO
     ENDDO
    ENDDO
  ENDIF

      
! C Solve for U1
! C Note, here the matrix will be indexed from 1...NY+1 corresponding to U1(0:NY)

! applying the wall model
! it has to be before APPL_BC calls, to make sure that BC values are computed
! correctly
  IF ( W_BC_YMIN.EQ.3 ) THEN
    CALL LES_WALL_MODEL
  ENDIF

! Initialize the matrix used to store implicit coefficients
  DO J=0,NY+1
    DO I=0,NXP
      MATLY(I,J)=0.d0
      MATDY(I,J)=1.d0
      MATUY(I,J)=0.d0
      VECY(I,J)=0.d0
    ENDDO
  ENDDO 
! 
! C Build the implicit system of equations for U1 
!$OMP PARALLEL DO PRIVATE(K,J,I)
  DO K=ZSTART,ZEND
      DO J=JSTART,JEND
	DO I=0,NXP
	  MATLY(I,J) = -TEMP1 * GMAT_22(K,J,1)
	  
	  MATDY(I,J)= INT_JACOB(K,J) + TEMP1 * (GMAT_22(K,J+1,1) + GMAT_22(K,J,1)) 
	  
	  MATUY(I,J)= -TEMP1 * GMAT_22(K,J+1,1)
	  
	  VECY(I,J)=R1X(I,K,J)
	END DO
      END DO

! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
       
    IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
      DO J=JSTART,JEND
	DO I=0,NXP
	  MATLY(I,J) = MATLY(I,J) - TEMP2 * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K,J-1) ) * GMAT_22(K,J,1)
	  
	  MATDY(I,J) = MATDY(I,J)                &
		      + TEMP2 * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K,J+1) ) * GMAT_22(K,J+1,1) &
		      + TEMP2 * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K,J-1) ) * GMAT_22(K,J,1)
	  
	  MATUY(I,J) = MATUY(I,J) - TEMP2 * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K,J+1) ) * GMAT_22(K,J+1,1)
	ENDDO
      ENDDO
    ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C Else, we are running in serial mode
! C Set the boundary conditions for U1
    CALL APPLY_BC_1_LOWER(K)
    CALL APPLY_BC_1_UPPER(K)
! C Now, solve the tridiagonal system for U1(:,k,:)
    CALL THOMAS_REAL(MATLY,MATDY,MATUY,VECY,NY+1,NXP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO J=0,NY+1
      DO I=0,NXP
	U1X(I,K,J)=VECY(I,J)
      ENDDO
    ENDDO
! End do k
  ENDDO
!$OMP END PARALLEL DO

     
! C Compute the horizontal viscous term in physical space and add to new RHS
! C already store in S1 
! C Important we have to multiply 1/2 with S1 before adding to RHS due
! C half time splitting during AD1

  DO J=JSTART,JEND
    DO K=ZSTART,ZEND
      DO I=0,NXP
	  R1X(I,K,J) = 0.5d0 * S1X(I,K,J) + INT_JACOB(K,J) * U1X(I,K,J)   &
	      + TEMP1 * (  		& 
		GMAT_22(K,J+1,1)*( U1X(I,K,J+1) - U1X(I,K,J) )		&
	      - GMAT_22(K,J,1) * ( U1X(I,K,J)   - U1X(I,K,J-1) )& 
	      )
      ENDDO
    ENDDO
  ENDDO

  
  IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
    DO J=JSTART,JEND
     DO K=ZSTART,ZEND
      DO I=0,NXP
	R1X(I,K,J) = R1X(I,K,J) + TEMP2 * ( & 
					    0.5d0 * (NU_T(I,K,J+1)+NU_T(I,K,J))			&
					    * ( GMAT_22(K,J+1,1) * ( U1X(I,K,J+1)-U1X(I,K,J) ) )	&
					    - 0.5d0 * (NU_T(I,K,J-1)+NU_T(I,K,J))			&
					    * ( GMAT_22(K,J,1) * ( U1X(I,K,J)- U1X(I,K,J-1) ) )	& 
					    )
      ENDDO
     ENDDO
    ENDDO
  ENDIF


  IF ( IC_TYPE == 7 ) THEN
    DO K=ZSTART,ZEND
	DO I=0,NXP
	  R1X(I,K,NY+1) = R1X(I,K,NY)
	ENDDO
    ENDDO  
  ENDIF
      
! Initialize the matrix used to store implicit coefficients
  DO K=0,NZ+1
    DO I=0,NXP
      MATLX(I,K)=0.d0
      MATDX(I,K)=1.d0
      MATUX(I,K)=0.d0
      VECX(I,K)=0.d0
    END DO
  END DO 

! C Build the implicit system of equations for TH 
!$OMP PARALLEL DO PRIVATE(K,J,I)
  DO J=JSTART,JEND
    DO K=ZSTART,ZEND
      DO I=0,NXP
	  MATLX(I,K) = -TEMP1 * GMAT_11(K,J,2)
	  
	  MATDX(I,K) = INT_JACOB(K,J) + TEMP1 * ( GMAT_11(K+1,J,2)+GMAT_11(K,J,2) ) 
	  
	  MATUX(I,K) = -TEMP1 * GMAT_11(K+1,J,2)
	  
	  VECX(I,K) = R1X(I,K,J)
      ENDDO
    ENDDO
! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
         
    IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
      DO K=ZSTART,ZEND
	DO I=0,NXP
	      MATLX(I,K) = MATLX(I,K)- TEMP2 * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K-1,J) ) * GMAT_11(K,J,2)
	      
	      MATDX(I,K) = MATDX(I,K)    &
			   + TEMP2 * 0.5d0 * ( NU_T(I,K+1,J)+NU_T(I,K,J) ) * GMAT_11(K+1,J,2) &
			   + TEMP2 * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K-1,J) ) * GMAT_11(K,J,2)
	      
	      MATUX(I,K) = MATUX(I,K) - TEMP2 * 0.5d0 * ( NU_T(I,K+1,J)+NU_T(I,K,J) ) * GMAT_11(K+1,J,2)
	ENDDO
      ENDDO
    ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C Set the boundary conditions for U1
    CALL APPLY_BC_1_LEFT(J)
    CALL APPLY_BC_1_RIGHT(J)
! C Now, solve the tridiagonal system for TH(:,k,:)
    IF ( U_BC_ZMAX  .EQ. 5 ) THEN
      D=GMAT_11(NZ,J,2)
      CALL THOMAS_REAL_SP(MATLX,MATDX,MATUX,VECX,D,NZ,NZ+1,NXP)
    ELSE   
      CALL THOMAS_REAL(MATLX,MATDX,MATUX,VECX,NZ+1,NXP)
    ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO K=0,NZ+1
      DO I=0,NXP
	U1X(I,K,J) = VECX(I,K)
      ENDDO
    ENDDO
! End do J
ENDDO
!$OMP END PARALLEL DO
   
!   IF ( U_BC_ZMIN .EQ. 6 ) THEN
!     Do J=0,NY+1
!       DO I=0,NXP
!          U1X(I,NZ+1,J)=U1X(I,1,J)
!       ENDDO
!     ENDDO

!     Do J=0,NY+1
!       DO I=0,NXP
!          U1X(I,0,J)=U1X(I,NZ,J)
!       ENDDO
!     ENDDO
!   ENDIF


! C ADI STEP FOR  U2 
! C Before going to  ADI1 we have to store RHS after subtracting previous U2_n 
! C at previous time step to bulid new RHS for  ADI2 based on U2_n+1/2 at intermidiate 


  DO J=JSTART,JEND
    DO K=ZSTART,ZEND
      DO I=0,NXP
	S1X(I,K,J) = R2X(I,K,J) - INT_JACOB(K,J)*U2X(I,K,J)
      ENDDO
    ENDDO
  ENDDO

! C Compute the horizontal viscous term in physical space and add to RHS
! C  at  ADI1 step .... Important in this step we have to add U2_n
! C at previous time step so that after multiplying a factor 1/2 into the 
! C RHS side during ADI1 operation RHS will be taken care U2_n not 1/2*U2_n
 
  IF (WAVE_ABS) THEN
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NXP
  ! C    wave absorb layer
	    R2X(I,K,J) = 0.5d0 * S1X(I,K,J) + INT_JACOB(K,J) * U2X(I,K,J) 		& 
			+ TEMP1 * (1.d0 + ART_VISCOSITY*SPONGE_SIGMA_OUT(K) )		&
				* ( GMAT_11(K+1,J,2) * ( U2X(I,K+1,J)-U2X(I,K,J) )	&
				    - GMAT_11(K,J,2) * ( U2X(I,K,J)-U2X(I,K-1,J) )	& 
				  )
	ENDDO
      ENDDO
    ENDDO
  ELSE
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NXP
	  R2X(I,K,J) = 0.5d0 * S1X(I,K,J) + INT_JACOB(K,J) * U2X(I,K,J) 	& 
			+ TEMP1 * ( GMAT_11(K+1,J,2) * ( U2X(I,K+1,J)-U2X(I,K,J) )	&
				    - GMAT_11(K,J,2) * ( U2X(I,K,J)-U2X(I,K-1,J) )	& 
				  )
	 ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NXP
	    R2X(I,K,J) = R2X(I,K,J) 			& 
			 + TEMP2 * ( 			& 
				    0.5d0 * ( NU_T(I,K+1,J)+NU_T(I,K,J) )			&
				    * ( GMAT_11_y(K+1,J,2) * ( U2X(I,K+1,J)-U2X(I,K,J) ) )	&
				    - 0.5d0 *( NU_T(I,K,J)+NU_T(I,K-1,J) )			&
				    * ( GMAT_11_y(K,J,2) * (U2X(I,K,J)-U2X(I,K-1,J) ) )	& 
				    )
	ENDDO
      ENDDO
    ENDDO
  ENDIF

! C Solve for U2
! C Note, here the matrix will be indexed from 1...NY+1 corresponding to U1(0:NY)

! Initialize the matrix used to store implicit coefficients
  DO J=0,NY+1
    DO I=0,NXP
      MATLY(I,J)=0.d0
      MATDY(I,J)=1.d0
      MATUY(I,J)=0.d0
      VECY(I,J)=0.d0
    END DO
  END DO 

! C Build the implicit system of equations for U1 
!$OMP PARALLEL DO PRIVATE(K,J,I)
  DO K=ZSTART,ZEND

      DO J=JSTART,JEND
	DO I=0,NXP
	      MATLY(I,J) = -TEMP1 * GMAT_22(K,J,1)  
	      
	      MATDY(I,J) = INT_JACOB(K,J) + TEMP1 * (GMAT_22(K,J+1,1)+GMAT_22(K,J,1) ) 
	      
	      MATUY(I,J) = -TEMP1 * GMAT_22(K,J+1,1)
	      
	      VECY(I,J) = R2X(I,K,J)
	ENDDO
      ENDDO

  ! IF using a subgrid model (LES) then add the eddy viscosity part implicitly

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
	DO J=JSTART,JEND
	  DO I=0,NXP
	      MATLY(I,J) = MATLY(I,J) - TEMP2 * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K,J-1) ) * GMAT_22_y(K,J,1)
	      
	      MATDY(I,J) = MATDY(I,J)   &
			    + TEMP2 * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K,J+1) ) * GMAT_22_y(K,J+1,1) &
			    + TEMP2 * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K,J-1) ) * GMAT_22_y(K,J,1)
	      
	      MATUY(I,J) = MATUY(I,J)- TEMP2 * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K,J+1) ) * GMAT_22_y(K,J+1,1)
	  ENDDO
	ENDDO
      ENDIF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! C Set the boundary conditions for U1
      CALL APPLY_BC_2_LOWER(K)
      CALL APPLY_BC_2_UPPER(K)
  ! C Now, solve the tridiagonal system for U1(:,k,:)
      CALL THOMAS_REAL(MATLY,MATDY,MATUY,VECY,NY+1,NXP)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      DO J=0,NY+1
	DO I=0,NXP
	  U2X(I,K,J)=VECY(I,J)
	ENDDO
      ENDDO
  ! End do k
  ENDDO
!$OMP END PARALLEL DO
     
! C Compute the horizontal viscous term in physical space and add to new RHS
! C already store in S1 
! C Important we have to multiply 1/2 with S1 before adding to RHS due
! C half time splitting during AD1

  DO J=JSTART,JEND
    DO K=ZSTART,ZEND
	DO I=0,NXP
	R2X(I,K,J) = 0.5d0 * S1X(I,K,J) + INT_JACOB(K,J) * U2X(I,K,J)  &
		    + TEMP1 * ( 	&	 
			        GMAT_22(K,J+1,1) * ( U2X(I,K,J+1)-U2X(I,K,J) )		&
			      - GMAT_22(K,J,1) * ( U2X(I,K,J)-U2X(I,K,J-1) )		& 
			      )

      ENDDO
    ENDDO
  ENDDO

  IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NXP
	    R2X(I,K,J) = R2X(I,K,J)  	& 
			+ TEMP2 * (   	& 
				    0.5d0 * ( NU_T(I,K,J+1)+NU_T(I,K,J) )			&
				  * ( GMAT_22_y(K,J+1,1) * ( U2X(I,K,J+1)-U2X(I,K,J) ) )	&
				  - 0.5d0*(NU_T(I,K,J-1)+NU_T(I,K,J))                   	&
				  * ( GMAT_22_y(K,J,1) * ( U2X(I,K,J)- U2X(I,K,J-1) ) ) 	& 
				  )
	ENDDO
      ENDDO
    ENDDO
  ENDIF


  IF ( IC_TYPE == 7 ) THEN
    DO K=ZSTART,ZEND
	DO I=0,NXP
	  R2X(I,K,NY+1) = R2X(I,K,NY)
	ENDDO
    ENDDO  
  ENDIF 

! Initialize the matrix used to store implicit coefficients
  DO K=0,NZ+1
  DO I=0,NXP
    MATLX(I,K)=0.d0
    MATDX(I,K)=1.d0
    MATUX(I,K)=0.d0
    VECX(I,K)=0.d0
  ENDDO
ENDDO 

! C Build the implicit system of equations for U2 
!$OMP PARALLEL DO PRIVATE(K,J,I)
  DO J=JSTART,JEND
      IF (WAVE_ABS) THEN
	DO K=ZSTART,ZEND
	  DO I=0,NXP
  ! C wave absorver layer
	    MATLX(I,K) = -TEMP1 * (1.d0 + ART_VISCOSITY * SPONGE_SIGMA_OUT(K)) * GMAT_11(K,J,2)
	    
	    MATDX(I,K) = INT_JACOB(K,J) - TEMP1 * ( 1.d0 + ART_VISCOSITY*SPONGE_SIGMA_OUT(K) )* &
			(-GMAT_11(K+1,J,2) - GMAT_11(K,J,2))
	    
	    MATUX(I,K) = -TEMP1 * ( 1.d0 + ART_VISCOSITY*SPONGE_SIGMA_OUT(K) ) * GMAT_11(K+1,J,2)
	    
	    VECX(I,K) = R2X(I,K,J)

	  ENDDO
	ENDDO
      ELSE 
	DO K=ZSTART,ZEND
	  DO I=0,NXP
	    MATLX(I,K) = -TEMP1 * GMAT_11(K,J,2)
	    
	    MATDX(I,K) = INT_JACOB(K,J)  + TEMP1 * ( GMAT_11(K+1,J,2)+GMAT_11(K,J,2) ) 
	    
	    MATUX(I,K) = -TEMP1 * GMAT_11(K+1,J,2)
	    
	    VECX(I,K) = R2X(I,K,J)
	  ENDDO
	ENDDO
      ENDIF
  ! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
      
      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
	DO K=ZSTART,ZEND
	  DO I=0,NXP
		MATLX(I,K) = MATLX(I,K) - TEMP2 * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K-1,J) ) * GMAT_11_y(K,J,2)
		
		MATDX(I,K) = MATDX(I,K) &
			      + TEMP2 * 0.5d0 * ( NU_T(I,K+1,J)+NU_T(I,K,J) ) * GMAT_11_y(K+1,J,2) &
			      + TEMP2 * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K-1,J) ) * GMAT_11_y(K,J,2)
		
		MATUX(I,K) = MATUX(I,K) - TEMP2 * 0.5d0 * ( NU_T(I,K+1,J)+NU_T(I,K,J) ) * GMAT_11_y(K+1,J,2)
	  ENDDO
	ENDDO
      ENDIF 


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! C Set the boundary conditions for U1
      CALL APPLY_BC_2_LEFT(J)
      CALL APPLY_BC_2_RIGHT(J)
  ! C Now, solve the tridiagonal system for TH(:,k,:)
      IF ( V_BC_ZMAX  .EQ. 5 ) THEN
	D=1.0d0
	CALL THOMAS_REAL_SP(MATLX,MATDX,MATUX,VECX,D,NZ,NZ+1,NXP)
      ELSE   
	CALL THOMAS_REAL(MATLX,MATDX,MATUX,VECX,NZ+1,NXP)
      ENDIF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO K=0,NZ+1
	DO I=0,NXP
	  U2X(I,K,J) = VECX(I,K)
	END DO
      END DO
  ! End do J
  ENDDO
!$OMP END PARALLEL DO
     
!  IF ( V_BC_ZMIN .EQ. 6 ) THEN
!    Do J=0,NY+1
!      DO I=0,NXP
!         U2X(I,NZ+1,J)=U2X(I,1,J)
!      ENDDO
!    ENDDO

!    Do J=0,NY+1
!     DO I=0,NXP
!         U2X(I,0,J)=U2X(I,NZ,J)
!      ENDDO
!    ENDDO
!  ENDIF
      
! C ADI STEP FOR  U3 
! C Before going to  ADI1 we have to store RHS after subtracting previous U3_n 
! C at previous time step to bulid new RHS for  ADI2 based on U3_n+1/2 at intermidiate 


  DO J=JSTART,JEND
    DO K=ZSTART,ZEND
      DO I=0,NXP
	S1X(I,K,J) = R3X(I,K,J) - INT_JACOB(K,J) * U3X(I,K,J)
      END DO
    END DO
  END DO

! C Compute the horizontal viscous term in physical space and add to RHS
! C  at  ADI1 step .... Important in this step we have to add U3_n
! C at previous time step so that after multiplying a factor 1/2 into the 
! C RHS side during ADI1 operation RHS will be taken care U3_n not 1/2*U3_n
 
  IF (WAVE_ABS) THEN
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NXP
	R3X(I,K,J) = 0.5d0 * S1X(I,K,J) + INT_JACOB(K,J) * U3X(I,K,J)		&
		    + TEMP1 * (1.d0 + ART_VISCOSITY*SPONGE_SIGMA_OUT(K) )   		&
		    * ( GMAT_11(K+1,J,2) * ( U3X(I,K+1,J)-U3X(I,K,J) )		&
		    - GMAT_11(K,J,2) * ( U3X(I,K,J)-U3X(I,K-1,J) )		& 
		      )
	ENDDO
      ENDDO
    ENDDO	
  ELSE
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NXP
	  R3X(I,K,J) = 0.5 * S1X(I,K,J) + INT_JACOB(K,J) * U3X(I,K,J)    &
		      + TEMP1 * ( 	& 
				 GMAT_11(K+1,J,2) * ( U3X(I,K+1,J)-U3X(I,K,J) )	&
				-GMAT_11(K,J,2) * ( U3X(I,K,J)-U3X(I,K-1,J) )	& 
				)
	ENDDO
      ENDDO
    ENDDO
  ENDIF

      
  IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NXP
	  R3X(I,K,J) = R3X(I,K,J) & 
		      + TEMP2 * (  & 
				  0.5d0 * ( NU_T(I,K+1,J)+NU_T(I,K,J) ) &
				  * ( GMAT_11_z(K+1,J,2) * ( U3X(I,K+1,J)-U3X(I,K,J)) )        &
				  - 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K-1,J) )                  &
				  *( GMAT_11_z(K,J,2) * ( U3X(I,K,J)-U3X(I,K-1,J)) )   & 
				 )
	ENDDO
      ENDDO
    ENDDO
  ENDIF
 
! C Solve for U3
! C Note, here the matrix will be indexed from 1...NY+1 corresponding to U3(0:NY)

! Initialize the matrix used to store implicit coefficients
  DO J=0,NY+1
    DO I=0,NXP
      MATLY(I,J)=0.d0
      MATDY(I,J)=1.d0
      MATUY(I,J)=0.d0
      VECY(I,J)=0.d0
    END DO
  END DO 

! C Build the implicit system of equations for U3 
!$OMP PARALLEL DO PRIVATE(K,J,I)
  DO K=ZSTART,ZEND

    DO J=JSTART,JEND
      DO I=0,NXP
	MATLY(I,J) = -TEMP1 * GMAT_22(K,J,1)
	
	MATDY(I,J) = INT_JACOB(K,J) + TEMP1 * ( GMAT_22(K,J+1,1)+GMAT_22(K,J,1) ) 
	
	MATUY(I,J) = -TEMP1 * GMAT_22(K,J+1,1)
	
	VECY(I,J) = R3X(I,K,J)
      ENDDO
    ENDDO

  ! IF using a subgrid model (LES) then add the eddy viscosity part implicitly

    IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
      DO J=JSTART,JEND
	DO I=0,NXP
	  MATLY(I,J) = MATLY(I,J) - TEMP2 * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K,J-1) ) * GMAT_22_z(K,J,1)
	  
	  MATDY(I,J) = MATDY(I,J) &
		      + TEMP2 * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K,J+1) ) * GMAT_22_z(K,J+1,1) &
		      + TEMP2 * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K,J-1) ) * GMAT_22_z(K,J,1)
	  
	  MATUY(I,J) = MATUY(I,J) - TEMP2 * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K,J+1) ) * GMAT_22_z(K,J+1,1)
       ENDDO
      ENDDO
    ENDIF        


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! C Set the boundary conditions for U1
    CALL APPLY_BC_3_LOWER(K)
    CALL APPLY_BC_3_UPPER(K)
  ! C Now, solve the tridiagonal system for U1(:,k,:)
    CALL THOMAS_REAL(MATLY,MATDY,MATUY,VECY,NY+1,NXP)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO J=0,NY+1
      DO I=0,NXP
	U3X(I,K,J)=VECY(I,J)
      ENDDO
    ENDDO
  ! End do k
  ENDDO
!$OMP END PARALLEL DO

      
! C     requred for convective bc
  C_avg = 0.0
  DO J=JSTART,JEND
    DO I=0,NXP
      C_avg = C_avg + U3X(I,NZ+1,J)
    ENDDO
  ENDDO	
		    
  C_avg = C_avg/real(NX*(JEND-JSTART+1))

  DO J=JSTART,JEND
    C_int(j) = 0.d0
    C_int_le(j) = 0.d0
    DO I=0,NXP
      C_int(J) = C_int(J) + U3X(I,NZ,J)
      C_int_le(J) = C_int_le(J) + U3X(I,1,J)
    ENDDO
    C_int(j) = C_int(j)/real(NX)
    C_int_le(j) = C_int_le(j)/real(NX)
  ENDDO
   
! c      DO J=JSTART,JEND  	        
! c       C_int(j) = C_int(j)/real(NX)
! c      ENDDO
!            
! C Compute the horizontal viscous term in physical space and add to new RHS
! C already store in S1 
! C Important we have to multiply 1/2 with S1 before adding to RHS due
! C half time splitting during AD1

  DO J=JSTART,JEND
    DO K=ZSTART,ZEND
	DO I=0,NXP
	R3X(I,K,J) = 0.5d0 * S1X(I,K,J) + INT_JACOB(K,J) * U3X(I,K,J)  &
		  + TEMP1 * (  & 
			      GMAT_22(K,J+1,1) * ( U3X(I,K,J+1)-U3X(I,K,J) ) &
			    - GMAT_22(K,J,1) * ( U3X(I,K,J)-U3X(I,K,J-1) )  & 
			    )

      END DO
    END DO
  END DO

  IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
    DO J=JSTART,JEND
     DO K=ZSTART,ZEND
      DO I=0,NXP
	R3X(I,K,J) = R3X(I,K,J) 			& 
		     + TEMP2 * (			& 
				+ 0.5d0 * ( NU_T(I,K,J+1)+NU_T(I,K,J) )			&
				* ( GMAT_22_z(K,J+1,1) * ( U3X(I,K,J+1)-U3X(I,K,J)) )		&
				- 0.5d0 * ( NU_T(I,K,J-1)+NU_T(I,K,J) )			&
				* ( GMAT_22_z(K,J,1)*( U3X(I,K,J)-U3X(I,K,J-1) ) )		& 
				)
      END DO
     END DO
    END DO
  ENDIF 


! Initialize the matrix used to store implicit coefficients
  DO K=0,NZ+1
    DO I=0,NXP
      MATLX(I,K)=0.d0
      MATDX(I,K)=1.d0
      MATUX(I,K)=0.d0
      VECX(I,K)=0.d0
    ENDDO
  ENDDO 

! C Build the implicit system of equations for U3
!$OMP PARALLEL DO PRIVATE(K,J,I)
  DO J=JSTART,JEND
    IF (WAVE_ABS) THEN
      DO K=ZSTART,ZEND
	DO I=0,NXP
	  MATLX(I,K) = -TEMP1*GMAT_11(K,J,2) * ( 1.d0 + ART_VISCOSITY*SPONGE_SIGMA_OUT(K) )
	  
	  MATDX(I,K) = INT_JACOB(K,J) + TEMP1 * ( 1.d0 + ART_VISCOSITY*SPONGE_SIGMA_OUT(K) )*  &
		      ( GMAT_11(K+1,J,2)+GMAT_11(K,J,2) )
	  	  
	  MATUX(I,K) = -TEMP1*GMAT_11(K+1,J,2) * ( 1.d0 + ART_VISCOSITY*SPONGE_SIGMA_OUT(K) )
	  
	  VECX(I,K)=R3X(I,K,J)
	ENDDO
      ENDDO
    ELSE
      DO K=ZSTART,ZEND
	DO I=0,NXP
	  MATLX(I,K) = -TEMP1 * GMAT_11(K,J,2) 
	  
	  MATDX(I,K) = INT_JACOB(K,J) + TEMP1 * ( GMAT_11(K+1,J,2)+GMAT_11(K,J,2) ) 
	  
	  MATUX(I,K) = -TEMP1 * GMAT_11(K+1,J,2)
	  
	  VECX(I,K)=R3X(I,K,J)
	END DO
    END DO
    ENDIF
  ! ! IF using a subgrid model (LES) then add the eddy viscosity part implicitly

    IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
      DO K=ZSTART,ZEND
	DO I=0,NXP
	  MATLX(I,K) = MATLX(I,K) - TEMP2 * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K-1,J) ) * GMAT_11_z(K,J,2)
	  
	  MATDX(I,K) = MATDX(I,K) &
			+ TEMP2 * 0.5d0 * ( NU_T(I,K+1,J)+NU_T(I,K,J) ) * GMAT_11_z(K+1,J,2) &
			+ TEMP2 * 0.5d0 * ( NU_T(I,K,J)+NU_T(I,K-1,J) ) * GMAT_11_z(K,J,2)
	  
	  MATUX(I,K) = MATUX(I,K) - TEMP2 * 0.5d0 * ( NU_T(I,K+1,J)+NU_T(I,K,J) ) * GMAT_11_z(K+1,J,2)
	END DO
      END DO
    ENDIF


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! C Set the boundary conditions for U1
    CALL APPLY_BC_3_LEFT(J)
    CALL APPLY_BC_3_RIGHT(J)
  ! C Now, solve the tridiagonal system for TH(:,k,:)
    IF ( W_BC_ZMAX  .EQ. 5 ) THEN
      D=1.
      CALL THOMAS_REAL_SP(MATLX,MATDX,MATUX,VECX,D,NZ,NZ+1,NXP)
    ELSE   
      CALL THOMAS_REAL(MATLX,MATDX,MATUX,VECX,NZ+1,NXP)
    ENDIF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO K=0,NZ+1
      DO I=0,NXP
	U3X(I,K,J) = VECX(I,K)
      ENDDO
    ENDDO
  ! End do J
  ENDDO
!$OMP END PARALLEL DO     

!   IF ( W_BC_ZMIN .EQ. 6 ) THEN
!     Do J=0,NY+1
!       DO I=0,NXP
!          U3X(I,NZ+1,J)=U3X(I,1,J)
!       ENDDO
!     ENDDO

!     Do J=0,NY+1
!       DO I=0,NXP
!          U3X(I,0,J)=U3X(I,NZ,J)
!       ENDDO
!     ENDDO
!   ENDIF


! updating lower corners 
! U3X(:,0,0)=U3X(:,1,0); U3X(:,NZ+1,0)=U3X(:,NZ,0);

! C If Variable timestepping and done with one full R-K step, update
! C DELTA_T based on the specified CFL number
! C This is not parallelized and should be used only in the serial
! C version to ensure that each process uses the same timestep
  IF ((VARIABLE_DT).and.(RK_STEP.eq.3).and.(MOD(TIME_STEP,UPDATE_DT).EQ.0)) THEN
    CALL COURANT_CURV_MPI
  ENDIF 

  CALL REAL_FOURIER_TRANS_U1 (.TRUE.)
  CALL REAL_FOURIER_TRANS_U2 (.TRUE.) 
  CALL REAL_FOURIER_TRANS_U3 (.TRUE.) 

  CALL allocation_R1 (.FALSE.)
  CALL allocation_R2 (.FALSE.)
  CALL allocation_R3 (.FALSE.)  

!saving Temp variables for wall model type 1
  IF ( (W_BC_YMIN.EQ.3) .AND. (WALL_MODEL_TYPE.EQ.2 .OR. WALL_MODEL_TYPE.EQ.4) ) THEN
        CU1X_WALL_MODEL(:,:,:)=CU1X(:,:,:)
        CALL MPI_BCAST_COMPLEX(CU1X_WALL_MODEL(0,:,:),NZV,NY+2)
        CU2X_WALL_MODEL(:,:,:)=CU2X(:,:,:)
        CALL MPI_BCAST_COMPLEX(CU2X_WALL_MODEL(0,:,:),NZV,NY+2)
        CU3X_WALL_MODEL(:,:,:)=CU3X(:,:,:)
        CALL MPI_BCAST_COMPLEX(CU3X_WALL_MODEL(0,:,:),NZV,NY+2)
  ENDIF
!C      CALL FFT_X_TO_FOURIER(U1,CU1,0,NY+1,0,NZ+1)
!C      CALL FFT_X_TO_FOURIER(U2,CU2,0,NY+1,0,NZ+1)
!C      CALL FFT_X_TO_FOURIER(U3,CU3,0,NY+1,0,NZ+1)
  DO N=1,N_TH
    CALL REAL_FOURIER_TRANS_TH (.TRUE.)
    CALL allocation_Rth (.FALSE.) 
    
    !saving Temp variables for wall model type 1
    IF ( (W_BC_YMIN.EQ.3) .AND. (WALL_MODEL_TYPE.EQ.2) ) THEN
        CTHX_WALL_MODEL(:,:,:,N)=CTHX(:,:,:,N)
        CALL MPI_BCAST_COMPLEX(CTHX_WALL_MODEL(0,:,:,N),NZV,NY+2)
    ENDIF
!C        CALL FFT_X_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N),0,NY+1,0,NZ+1)
  ENDDO

! C Begin second step of the Fractional Step algorithm, making u divergence free
! C The following subroutine projects Uhat onto divergence free space
  CALL REM_DIV_CURV
  
 
! C Now, phi is stored in CR1X, use this to update the pressure field
! C Note, here we divide by H_BAR since it was absorbed into PHI in REM_DIV
   DO J=JSTART-1,JEND+1
     DO K=ZSTART-1,ZEND+1
       DO I=0,NX2P
 	 CPX(I,K,J) = CPX(I,K,J) + CR1X(I,K,J) / TEMP4
       ENDDO
     ENDDO
   ENDDO

  RETURN
END

! C--.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE REM_DIV_CURV
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
  USE mg_vari, 		ONLY : INIT_FLAG, MUDPACK_BC
  USE pr_rem
  USE mpi_var      
  USE, INTRINSIC :: IEEE_ARITHMETIC

    IMPLICIT NONE 

  !      REAL*8 P_TMP(0:NZ+1,0:NY+1), P_TMP2(0:NZ+1,0:NY+1), &
  !       RHS_TMP(0:NZ+1,0:NY+1)
  !      REAL*8 RTMP(0:NXP,0:NZV-1,0:NY+1)
  !      COMPLEX*16 CRTMP(0:NX2P,0:NZV-1,0:NY+1)
  !      EQUIVALENCE (RTMP,CRTMP)

    REAl*8 DIV  
    INTEGER I,J,K

! c     Set BCs for phi (has been read from input file)
! c     MUDPACK_BC(1) : left, MUDPACK_BC(2) : right, MUDPACK_BC(3) : bottom, MUDPACK_BC(4) : top
! c     For periodic and dirchlet set RHS to corresponding values
! c     For neumann, set RHS to the derivative times the grid spacing
! c     normal to the wall


    CALL allocation_ub(.FALSE.)

    CR1X = (0.0d0,0.0d0)
!      CRTMP(:,:,:) = (0.0d0,0.0d0)
    CS1X = (0.0d0,0.0d0)
    CU2bX = (0.0d0,0.0d0)
    CU3bX = (0.0d0,0.0d0)
      
    DO J=JSTART,JEND+1
      DO K=ZSTART,ZEND+1
	DO I=0,NX2P
! C      U2 required to define at bottom cell face 	  
	  CU2bX(I,K,J) = 0.5d0 * CJOB_22(K,J,1) * ( CU2X(I,K,J)+CU2X(I,K,J-1) )  &
		       + 0.5d0 * CJOB_21(K,J,1) * ( CU3X(I,K,J)+CU3X(I,K,J-1) )
! C      U3 required to define at side cell face     	    
	  CU3bX(I,K,J) = 0.5d0 * CJOB_11(K,J,2) * ( CU3X(I,K,J)+CU3X(I,K-1,J) ) &
		       + 0.5d0 * CJOB_12(K,J,2) * ( CU2X(I,K,J)+CU2X(I,K-1,J) )
	ENDDO
      ENDDO
    ENDDO
      
     
      
! C Now, create the RHS vector
    DO J=JSTART,JEND
      DO K=ZSTART-1,ZEND
       DO I=0,NX2P
	  CR1X(I,K,J) = CIKXP(I) * INT_JACOB(K,J) * CU1X(I,K,J) &
			+ ( CU2bX(I,K,J+1)-CU2bX(I,K,J) ) + ( CU3bX(I,K+1,J)-CU3bX(I,K,J) )
	  CS1X(I,K,J) = CR1X(I,K,J)
       ENDDO
      ENDDO
    ENDDO
!  		print *, 'check 3'

 
!       CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS1X,CS1Z)
!       cvarp(:,:,:)=(0.d0,0.d0)
!       DO I=0,NKX
!       cvarp(I,:,:)=CS1Z(I,:,:)
!       ENDDO
!       CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
!       DO I=0,NXM
!       S1Z(I,:,:)=varp(I,:,:)
!       ENDDO
!       S1Z(NXM+1:NXV-1,:,:)=0.0
!       CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)
! 
! !C      CALL FFT_X_TO_PHYSICAL(CRTMP,RTMP,0,NY+1,0,NZ+1)
! 
!       DIV = 0.0D0
!       DO J = 1,NY
!        DO K = 1,NZ
!         DO I = 1,NXP
!          DIV = DIV + DABS(S1X(I,K,J))
!         ENDDO
!        ENDDO
!       ENDDO
!       DIV = DIV/(DBLE(NXP*NP)*DBLE(NY)*DBLE(NZ))
!       CALL MPI_COMBINE_STATS(DIV,1,1)
! 
!       IF (rank .eq. 0) THEN
! 
!       write(6,*) "Divergence before projection step: ", DIV
!      
!       open(99,file='out_screen.txt',form='formatted', &
!       status='unknown', position='append')
! 
!       write(99,*) "Divergence before projection step: ", DIV
!       ENDIF

    IF ( INIT_FLAG .EQ. .TRUE. ) THEN
      IF (MUDPACK_USE) THEN
	
	DO I = 0,NX2P
	  kx2_c = KX2P(I)
	  RHS_TMP(:,:) = 0.0d0
	  
	  DO K = ZSTART,ZEND
	    DO J = JSTART,JEND
	      RHS_TMP(K,J) = dble(CR1X(I,K,J))
	    ENDDO
	  ENDDO

!       IF((MUDPACK_BC(1) .EQ. 0) .AND. (MUDPACK_BC(2) .EQ. 0)) THEN
!          RHS_TMP(NZ+1,:) = RHS_TMP(1,:)
!          RHS_TMP(0,:) = RHS_TMP(NZ,:)
!       ELSE
!          RHS_TMP(NZ+1,:) = RHS_TMP(NZ,:)
!          RHS_TMP(0,:) = RHS_TMP(1,:)      
!       ENDIF

	  P_TMP(:,:) = 0.0d0
	  CALL MULTIGRID_MUDPACK(P_TMP)
	ENDDO
	
      ELSE
	
	DO I = 0,NX2P
	P_TMP(:,:) = 0.0d0; RHS_TMP(:,:) = 0.0d0
	DO K = ZSTART,ZEND
	  DO J = JSTART,JEND
	  RHS_TMP(K,J) = dble(CR1X(I,K,J))
  ! C           RHS_TMP(K,J) =1.0
	  ENDDO
	ENDDO
	CALL MULTIGRID_DMGD9V(P_TMP,RHS_TMP,I)
	ENDDO
	
      ENDIF
    ENDIF
    INIT_FLAG = .FALSE.



! c  Depends on which Divergence_removal package you have choosen in INPUT file

    IF (MUDPACK_USE) THEN
     
     DO I = 0,NX2P
      kx2_c = KX2P(I)
! c     Solve for complex part of phi
       RHS_TMP(:,:) = 0.0d0
       DO K = ZSTART,ZEND
        DO J = JSTART,JEND
         RHS_TMP(K,J) = -dimag(CR1X(I,K,J))
        ENDDO
       ENDDO

!       IF((MUDPACK_BC(1) .EQ. 0) .AND. (MUDPACK_BC(2) .EQ. 0)) THEN
!          RHS_TMP(NZ+1,:) = RHS_TMP(1,:)
!          RHS_TMP(0,:) = RHS_TMP(NZ,:)
!       ELSE
!          RHS_TMP(NZ+1,:) = RHS_TMP(NZ,:)
!          RHS_TMP(0,:) = RHS_TMP(1,:) 
!       ENDIF

       P_TMP(:,:) = 0.0d0
       CALL MULTIGRID_MUDPACK(P_TMP)

! C     Solve for real part of phi
       RHS_TMP(:,:) = 0.0d0
       DO K = ZSTART,ZEND
        DO J = JSTART,JEND
         RHS_TMP(K,J) = -dble(CR1X(I,K,J))
        ENDDO
       ENDDO

!       IF((MUDPACK_BC(1) .EQ. 0) .AND. (MUDPACK_BC(2) .EQ. 0)) THEN
!          RHS_TMP(NZ+1,:) = RHS_TMP(1,:)
!          RHS_TMP(0,:) = RHS_TMP(NZ,:)
!       ELSE
!          RHS_TMP(NZ+1,:) = RHS_TMP(NZ,:)
!          RHS_TMP(0,:) = RHS_TMP(1,:)
!       ENDIF

       P_TMP2(:,:) = 0.0d0
       CALL MULTIGRID_MUDPACK(P_TMP2)

       DO J = 0,NY+1
        DO K = 0,NZ+1
         CR1X(I,K,J) = CMPLX(P_TMP2(K,J),P_TMP(K,J))
        ENDDO
       ENDDO 

     ENDDO
    
    
    ELSE
	DO I = 0,NX2P
    ! c     Solve for complex part of phi
	  RHS_TMP(:,:) = 0.0d0
	  DO K = ZSTART,ZEND
	    DO J = JSTART,JEND
	      RHS_TMP(K,J) =-dimag(CR1X(I,K,J))
	    ENDDO
	  ENDDO
	  P_TMP(:,:) = 0.0d0
	
	  CALL MULTIGRID_DMGD9V(P_TMP,RHS_TMP,I)


    ! C     Solve for real part of phi
	  RHS_TMP(:,:) = 0.0d0
	  DO K = ZSTART,ZEND
	    DO J = JSTART,JEND
	      RHS_TMP(K,J) = -dble(CR1X(I,K,J))
	    ENDDO
	  ENDDO

	  P_TMP2(:,:) = 0.0d0

	  CALL MULTIGRID_DMGD9V(P_TMP2,RHS_TMP,I)


	  DO J = 0,NY+1
	    DO K = 0,NZ+1
	      CR1X(I,K,J) = CMPLX(P_TMP2(K,J),P_TMP(K,J))
	    ENDDO
	  ENDDO 

	ENDDO 
    ENDIF

    
!END OF I LOOP 


! C Now, Solve for CUi at cell face the divergenceless velocity field

    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NX2P
! c	    CU1(I,K,J) = CU1(I,K,J)- CIKXP(I)*CR1X(I,K,J)
	  CU2bX(I,K,J) = CU2bX(I,K,J) - GMAT_22(K,J,1) * ( CR1X(I,K,J)-CR1X(I,K,J-1) ) &
			- 0.25d0 * GMAT_12(K,J,1) * ( CR1X(I,K+1,J) + CR1X(I,K+1,J-1) - CR1X(I,K-1,J) - CR1X(I,K-1,J-1) )
	  
	  CU3bX(I,K,J) = CU3bX(I,K,J) - GMAT_11(K,J,2) * ( CR1X(I,K,J)-CR1X(I,K-1,J) ) &
			-0.25d0 * GMAT_12(K,J,2) * (CR1X(I,K,J+1) + CR1X(I,K-1,J+1) - CR1X(I,K,J-1) - CR1X(I,K-1,J-1) )
	ENDDO
      ENDDO
    ENDDO

! C  Calculating velocity at the cell center from the 
! C  velocity at the cell face. 

    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NX2P
	  CU1X(I,K,J) = CU1X(I,K,J) - CIKXP(I)*CR1X(I,K,J)
	  CU2X(I,K,J) = CU2X(I,K,J)  &
			- 0.5d0 * (  & 
				    CJOB_22(K,J+1,1) * ( CR1X(I,K,J)+CR1X(I,K,J+1) ) &
				  - CJOB_22(K,J,1) * ( CR1X(I,K,J)+CR1X(I,K,J-1) ) 	& 
				  + CJOB_12(K+1,J,2) * ( CR1X(I,K+1,J)+CR1X(I,K,J) ) &
				  - CJOB_12(K,J,2) * ( CR1X(I,K,J)+CR1X(I,K-1,J) )	& 
				  ) / INT_JACOB(K,J)
    
	  CU3X(I,K,J) = CU3X(I,K,J)   &
			- 0.5d0 * (   & 
				    CJOB_11(K+1,J,2)*(CR1X(I,K,J) + CR1X(I,K+1,J)) &
				  - CJOB_11(K,J,2)*(CR1X(I,K,J) + CR1X(I,K-1,J))   &
				  + CJOB_21(K,J+1,1)*(CR1X(I,K,J+1) + CR1X(I,K,J)) &
				  - CJOB_21(K,J,1)*(CR1X(I,K,J) + CR1X(I,K,J-1))   & 
				  ) / INT_JACOB(K,J)
	ENDDO
      ENDDO
    ENDDO
      
    
!     DO I=0,NX2P
!       DO K=ZSTART-1,ZEND+1
!       CU1X(I,K,NY+1) = CU1X(I,K,NY) !- CU1X(I,NZ-1,J)
! c        CU2X(I,K,NY+1) = CU2X(I,K,NY) !- CU2X(I,NZ-1,J)
!       CU3X(I,K,NY+1) = CU3X(I,K,NY) !- CU3X(I,NZ-1,J)
!       ENDDO
!     ENDDO

! C  Now, create the RHS vector
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NX2P
	    CS1X(I,K,J)= CIKXP(I) * INT_JACOB(K,J) * CU1X(I,K,J) &
			+ ( CU2bX(I,K,J+1)-CU2bX(I,K,J) )   &
			+ ( CU3bX(I,K+1,J)-CU3bX(I,K,J) )
	ENDDO
      ENDDO
    ENDDO

    CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS1X,CS1Z)
    cvarp(:,:,:)=(0.d0,0.d0)
    
    DO I=0,NKX
      cvarp(I,:,:)=CS1Z(I,:,:)
    ENDDO
    
    CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
    
    DO I=0,NXM
      S1Z(I,:,:)=varp(I,:,:)
    ENDDO
    
    S1Z(NXM+1:NXV-1,:,:)=0.0
    CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)

!      CALL FFT_X_TO_PHYSICAL(CRTMP,RTMP,0,NY+1,0,NZ+1)

    DIV = 0.0D0
    
    DO J = 2,NYM
      DO K = 2,NZM
	DO I = 1,NXP
	  DIV = DIV + DABS(S1X(I,K,J))
	ENDDO
      ENDDO
    ENDDO
    
    DIV = DIV / ( DBLE(NXP*NP) * DBLE(NY) * DBLE(NZ) )

    CALL MPI_COMBINE_STATS(DIV,1,1)

    IF (rank .eq. 0) THEN
      WRITE(6,*) "Divergence after the projection step:", DIV
    ENDIF
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!     deallocating cubs
!ccccccccccccccccccccccccccccccccccccccccccccccccc

     IF (IEEE_IS_NAN(DIV)) STOP 'DIV is a NaN' 
     IF (DELTA_T .LT. 1.0E-4) STOP 'DELTA_T is zero'

    CALL allocation_ub(.TRUE.)

		
  RETURN
END      
     
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE CREATE_TH_DUCT
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! C Initialize the scalar fields
! C In this subroutine, you should initialize each scalar field for the
! C particular problem of interest

  USE ntypes
  USE Domain
  USE Grid
  !USE Fft_var, only : NKX
  USE TIME_STEP_VAR
  USE run_variable
  USE mg_vari, 		ONLY : INIT_FLAG
  USE mpi_var, 		ONLy: rank
      
    IMPLICIT NONE 

    INTEGER I,J,K,N

    REAL*8  RNUM1, DAMP_FACT 

    DO N=1,N_TH
	IF (CREATE_NEW_TH(N)) THEN
	  IF(rank.EQ.0) write(*,*) 'A new thfield has been created'
	  IF ( IC_TYPE.eq.0 ) THEN
	  ELSE IF ( (IC_TYPE.eq.5).OR.(IC_TYPE.eq.6) )then
	    DO J=0,NY+1
	      DO K=0,NZ+1
		DO I=0,NXP
		  THX(I,K,J,N)=U3_BAR(1,J) 
		END DO
	      END DO
	    ENDDO
	  ELSE IF ( IC_TYPE.eq.7 ) THEN
	    
	    DO K=0,NZ+1
	      DO I=0,NXP
		DO J=0,NY+1
		  THX(I,K,J,N)=0.d0 !TH_BC_YMAX_C1(N)*(ypoint(K,J)-ypoint(1,0))
! & 			              +TH_BC_YMIN_C1(N) 
!             			TH_BAR(K,J)= (-ypoint(k,j)+ypoint(0,NY+1))
		ENDDO
	      ENDDO
	    ENDDO
	  ENDIF

	  CALL REAL_FOURIER_TRANS_TH (.true.)      
	  CALL REAL_FOURIER_TRANS_Rth (.true.)

	  CALL allocation_Fth(.false.)

!! initial value for random numbers
!          RNUM1 = 0.d0
!     DO I=1,min(NX2P,NX2P_L)
!      DO J=1,NY
!        DO K=1,NZ
!      ! C Now, give the velocity field a random perturbation
!
!          CALL RANDOM_NUMBER(RNUM1)
!          CTHX(I,K,J,N)=CTHX(I,K,J,N)+(RNUM1-0.5d0) * KICK
!          IF (J .EQ. 0) THEN
!            DAMP_FACT = 1.0
!          ELSE
!            DAMP_FACT = exp(-2 * dble((J)**2)/dble((NY)**2))
!          ENDIF
!        
!          CTHX(I,K,J,N)=CTHX(I,K,J,N) * DAMP_FACT
!        ENDDO
!      ENDDO
!    ENDDO

      ENDIF
    ENDDO      

  RETURN
END 


! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE CREATE_FLOW_DUCT
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var, 		ONLY : pi
  USE TIME_STEP_VAR
  USE run_variable
  USE mg_vari, 		ONLY : INIT_FLAG

    IMPLICIT NONE 

    INTEGER I,J,K,MM,NN,M,N
    
    LOGICAL dir_exists, INPUT_V_INI
    
    REAL*8  RNUM1,RNUM2,RNUM3,SUM1, &
	    FACT1,FACT2,FACT3,FACT5,FACT, DAMP_FACT
    REAL*8  NGY(NY), NGZ(NZ)
    
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    REAL(r8),ALLOCATABLE,DIMENSION(:,:)   :: tmpU,tmpW
    ALLOCATE (tmpU(0:NZ+1,0:NY+1),tmpW(0:NZ+1,0:NY+1))

! C Initialize the random number generator
    CALL RANDOM_SEED(SIZE = K)
    ALLOCATE (seed(1:K))
    seed(1:K)=10
    CALL RANDOM_SEED(PUT = seed)

! C WBULK0 and KICK should be set in input.dat

! C Set the laminar velocity profile in physical space
    WRITE(*,*) 'UBULK: ',WBULK0

    IF (IC_TYPE.eq.0) THEN
    ! C For closed duct flow
      mm = 5 
      nn = 5 
      fact  = LZ/LY

      DO K =1,NZ
	  NGZ(K) = 2.0*GZF(K)/LZ ;
      ENDDO 
      
      DO J =1,NY
	  NGY(J) = 2.0*GYF(J)/LY ;
      ENDDO

      DO J=1,NY
	DO K=1,NZ
	  SUM1 = 0.0 ;

	  DO m=0,mm
	    DO n=0,nn
	      fact1 = (2.0*m+1) ;
	      fact2 = (2.0*n+1) ;
	      fact3 = fact1*fact2**3.0 + fact**2.0*fact2*fact1**3.0 ;
	      fact5 = ((-1.0)**(m+n))*fact/fact3 ;
	      sum1 = sum1 + fact5*cos(fact1*pi*NGZ(K)/2)*  &
		    cos(fact2*pi*NGY(J)/2) ;
	    ENDDO
	  ENDDO

	  DO I=0,NXP
	    U1X(I,K,J)=WBULK0*sum1 
	    U2X(I,K,J)=0.d0
	    U3X(I,K,J)=0.d0
	  ENDDO
	  
	ENDDO
      ENDDO
    ELSE IF ((IC_TYPE.eq.1) .OR. (IC_TYPE.eq.4) ) THEN
    ! C For open channel flow :
      DO K=0,NZM
	DO I=0,NXP
	  DO J=1,NY
      !            U1(I,K,J)=-(3./2.)*WBULK0*GYF(J)**2.+3.*WBULK0*GYF(J)
      !            U1(I,K,J)=(-GYF(J)**2.d0+(NU+LY**2.d0)*GYF(J)/LY)/NU
	    U1X(I,K,J)= (-GYF(J)**2.d0+ 2.0*LY*GYF(J))/LY**2.0
      !             U1(I,K,J)=0.d0
	    U2X(I,K,J)=0.d0
	    U3X(I,K,J)=0.d0
	  ENDDO
	  U1X(I,K,0)=0.d0
	  U3X(I,K,0)=0.d0
	  U1X(I,K,NY+1)=0.d0
	  U3X(I,K,NY+1)=0.d0
	ENDDO
      ENDDO
      
    ELSE IF (IC_TYPE.eq.2) THEN
    ! C For Couette flow:
      DO J=0,NY
	DO K=0,NZM
	  DO I=0,NXP
	    U1X(I,K,J)=gyf(j)
	    U2X(I,K,J)=0.d0
	    U3X(I,K,J)=0.d0
	  ENDDO
	ENDDO
      ENDDO

    ELSE IF (IC_TYPE.eq.3) THEN
    ! Shear layer
      DO J=0,NY+1
	DO K=0,NZ+1
	  DO I=0,NXP
	    U1X(I,K,J)=TANH(GYF(J)*20.d0)
	    U2X(I,K,J)=0.d0
	    U3X(I,K,J)=0.d0
	  END DO
	END DO
      END DO
      
    ELSE IF ( (IC_TYPE.eq.5).OR.(IC_TYPE.eq.6) ) THEN
      DO J=0,NY+1
	DO K=0,NZ+1
	  DO I=0,NXP
	      U1X(I,K,J)=U1_BAR(J)
	      U2X(I,K,J)=0.d0 !U2_BAR(K,J)
	      U3X(I,K,J)=1.d0 !U3_BAR(1,J) 
	  ENDDO
	ENDDO	 
      ENDDO
      
    ELSE IF (IC_TYPE.eq.7) THEN

      INPUT_V_INI= .FALSE.  
      
      IF (INPUT_V_INI) THEN 
      
	WRITE(6,*) 'VELOCITY feild has been read from input file'
	
	INQUIRE(FILE='velocity_data_input.dat', EXIST=dir_exists) 
	
	IF ( .NOT. dir_exists) THEN
	  WRITE(6,*) ' "velocity_data_input.dat" does not exists'
	  WRITE(*,*) ' "velocity_data_input.dat" does not exists'
	  STOP
	ENDIF  
	
	
	OPEN(1231,file='velocity_data_input.dat',form='formatted', status='old')    
	      
	DO J=0,NY+1
	  DO K=0,NZ+1
	    READ(1231,*) tmpW(K,J),tmpU(K,J)              
	  ENDDO
	ENDDO 
           
	CLOSE(1231)
	
	
	DO J=0,NY+1
          DO K=0,NZ+1
           DO I=0,NXP
	      U3X(I,K,J)= tmpW(K,J)
	      U1X(I,K,J)= 0.0d0 +tmpU(K,J)
	      U2X(I,K,J)= 0.0d0 
           ENDDO
         ENDDO	 
        ENDDO
	
      ELSE
       
      
	DO J=0,NY+1
	  DO K=0,NZ+1
	    DO I=0,NXP
	      !U3X(I,K,J)= 1.0d0 * WBULK0 * (1.d0-ypoint(I,J)**2.d0)   !U3_BAR(k,J)
	      U3X(I,K,J)= 1.0d0 * WBULK0
	      U1X(I,K,J)= 0.0d0 
	      U2X(I,K,J)= 0.0d0 
	    ENDDO
	  ENDDO	 
	ENDDO
      
      ENDIF
      
      WBULK0=U3X(0,1,NY)
    ! c	DO J=0,NY+1
    ! c          DO K=0,NZ+1
    ! c           DO I=0,NXM
    ! c	    U3(I,K,J)=U3(I,K,J)/30.0
    ! c	   END DO
    ! c         END DO	 
    ! c        ENDDO   
    ENDIF

    CALL REAL_FOURIER_TRANS_U1 (.true.)
    CALL REAL_FOURIER_TRANS_U2 (.true.)      
    CALL REAL_FOURIER_TRANS_U3 (.true.)
    CALL REAL_FOURIER_TRANS_P (.true.)

    CALL REAL_FOURIER_TRANS_R1 (.true.)
    CALL REAL_FOURIER_TRANS_R2 (.true.)
    CALL REAL_FOURIER_TRANS_R3 (.true.)

    CALL allocation_F1 (.false.)
    CALL allocation_F2 (.false.)
    CALL allocation_F3 (.false.)

!C       CALL FFT_X_TO_FOURIER(U1,CU1,0,NY+1,0,NZ+1)
!C       CALL FFT_X_TO_FOURIER(U2,CU2,0,NY+1,0,NZ+1)
!C       CALL FFT_X_TO_FOURIER(U3,CU3,0,NY+1,0,NZ+1)

! initial value for random numbers
    RNUM1 = 0.d0
    RNUM2 = 0.d0
    RNUM3 = 0.d0

    DO I=1,min(NX2P,NX2P_L)
      DO J=JSTART,JEND
	DO K=ZSTART,ZEND
      ! C Now, give the velocity field a random perturbation

 	  CALL RANDOM_NUMBER(RNUM1)
	  CALL RANDOM_NUMBER(RNUM2)
	  CALL RANDOM_NUMBER(RNUM3)
	  IF (IC_TYPE.eq.3) THEN
      ! C If we are initializing with a shear layer 
	    CU1X(I,K,J)=CU1X(I,K,J) &
			+(RNUM1-0.5d0)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)
	    CU2X(I,K,J)=CU2X(I,K,J) &
			+(RNUM1-0.5d0)*KICK*EXP(-(GY(J)*20.d0)**2.d0)
	    CU3X(I,K,J)=CU3X(I,K,J) &
			+(RNUM1-0.5d0)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)

	  ELSE
! 	      Adding a random number to initial contdition based on this energy spectrum: E(k) = (K/K0)^4 * exp (-2*(K/K0)^2)
! 	      CU1X(I,K,J)=CU1X(I,K,J)+(RNUM1-0.5d0) * KICK * (5.d0*I/NX)**4.d0 * EXP(-2.d0 * (5.d0*I/NX)**2.d0 )
! 	      CU2X(I,K,J)=CU2X(I,K,J)+(RNUM2-0.5d0) * KICK * (5.d0*I/NX)**4.d0 * EXP(-2.d0 * (5.d0*I/NX)**2.d0 )
! 	      CU3X(I,K,J)=CU3X(I,K,J)+(RNUM3-0.5d0) * KICK * (5.d0*I/NX)**4.d0 * EXP(-2.d0 * (5.d0*I/NX)**2.d0 )
 
	      CU1X(I,K,J)=CU1X(I,K,J)+(RNUM1-0.5d0) * KICK 
	      CU2X(I,K,J)=CU2X(I,K,J)+(RNUM2-0.5d0) * KICK 
	      CU3X(I,K,J)=CU3X(I,K,J)+(RNUM3-0.5d0) * KICK 
	      
	  ENDIF
	    
	    IF (J .EQ. 0) THEN
		DAMP_FACT = 1.0
	    ELSE
  ! 	      DAMP_FACT = 1
  ! 	      DAMP_FACT = exp(-20.d0*dble(J)/dble(NY))
		DAMP_FACT = exp(-6 * dble((J)**2)/dble((NY)**2))
	    ENDIF
	    
	    CU1X(I,K,J)=CU1X(I,K,J) * DAMP_FACT
	    CU2X(I,K,J)=CU2X(I,K,J) * DAMP_FACT
	    CU3X(I,K,J)=CU3X(I,K,J) * DAMP_FACT
	    
	ENDDO
      ENDDO
    ENDDO

!      DO J=0,0.5*NY
!        DO K=1,NZ
!          CALL RANDOM_NUMBER(RNUM1)
!          CALL RANDOM_NUMBER(RNUM2)
!          KICK=0.10
!          CU1X(0,K,J)=CU1X(0,K,J)*(1.0+(RNUM1-0.5d0) * KICK * exp(-10*dble((J)**2)/dble((NY)**2)))
!          CU3X(0,K,J)=CU3X(0,K,J)*(1.0+(RNUM2-0.5d0) * KICK * exp(-10*dble((J)**2)/dble((NY)**2)))
!        ENDDO
!      ENDDO




    WRITE(*,*) 'KICK is : ',KICK, 'Bellow j', RNUM1

!      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CU1X,CU1Z)
!      CU1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
!      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CU1Z,CU1X)

!      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CU2X,CU2Z)
!      CU2Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
!      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CU2Z,CU2X)

!      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CU3X,CU3Z)
!      CU3Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
!      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CU3Z,CU3X)

! C Remove the divergence of the velocity field
    CALL SAVE_STATS_CURVI(.FALSE.)

    INIT_FLAG = .TRUE.
    
    CALL REM_DIV_CURV
    CALL SAVE_STATS_CURVI(.FALSE.)

  RETURN
END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE sponge_th(N)
! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state for the temperature field
! The intention is to allow an open boundary
  USE ntypes
  USE Domain
  USE Grid
  !USE Fft_var, only : NKX
  USE TIME_STEP_VAR
  USE run_variable
  USE mg_vari, ONLY : INIT_FLAG
  USE mpi_var  

    IMPLICIT NONE 

    ! The following variables will store the background state
    INTEGER i,j,k,N, NY_sponge
    REAL*8 Sponge_Amp, Dynamic_sponge_threshold
    REAL*8 CU1X_sponge(0:NZ), CU2X_sponge, CU3X_sponge
    REAL*8 CU1X_target, CU2X_target, CU3X_target    
    LOGICAL DYNAMIC_SPONE
    COMPLEX(r8) tmp
    
! Dynamic sponge, applies the spone whenever the vertically averaged zonal velocity is higher than a specified value 
    DYNAMIC_SPONE = .TRUE.
    
!   Sponge_Amp = 1.0d0
    Sponge_Amp = 0.1d0/DELTA_T
    Dynamic_sponge_threshold = 0.001d0    !2.5% tollerance for dynamic sponge
    
! Damp fluctuation flowfield
	!TOP AND BOTTOM 
    DO j=jstart,jend
      DO k=zstart,zend
	DO i=1,NX2P
	    CF1X(i,k,j)=CF1X(i,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU1X(i,k,j)-0.d0)* INT_JACOB(K,J)  
	    CF2X(i,k,j)=CF2X(i,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU2X(i,k,j)-0.d0)* INT_JACOB(K,J)
	    CF3X(i,k,j)=CF3X(i,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU3X(i,k,j)-0.d0)* INT_JACOB(K,J)
	    CFTHX(i,k,j,N)=CFTHX(i,k,j,N) - Sponge_Amp * SPONGE_SIGMA(j) * (CTHX(i,k,j,N)-0.d0)*INT_JACOB(K,J)
	ENDDO
      ENDDO
    ENDDO


!INFLOW OUTFLOW
    DO j=jstart,jend
      DO k=zstart,zend
	DO i=1,NX2P
	    CF1X(i,k,j)=CF1X(i,k,j)- Sponge_Amp * SPONGE_SIGMA_OUT(k)*(CU1X(i,k,j)-0.d0) * INT_JACOB(K,J)
	    CF2X(i,k,j)=CF2X(i,k,j)- Sponge_Amp * SPONGE_SIGMA_OUT(k)*(CU2X(i,k,j)-0.d0) * INT_JACOB(K,J)
	    CF3X(i,k,j)=CF3X(i,k,j)- Sponge_Amp * SPONGE_SIGMA_OUT(k)*(CU3X(i,k,j)-0.d0) * INT_JACOB(K,J)
	    CFTHX(i,k,j,N)=CFTHX(i,k,j,N) - Sponge_Amp * SPONGE_SIGMA_OUT(k) * (CTHX(i,k,j,N)-0.d0)*INT_JACOB(K,J)
	ENDDO
      ENDDO
    ENDDO

! Appling equavalent pressure gradient in span-wise direction at top of atmospher! the span-wise velocity has been foreced to be zero in top of atmospher!     
    IF (NY_S_T .GT. NY) THEN
    ! deafult value
      NY_sponge = INT (2 * NY/3 )
    ELSE
    ! specified value
      NY_sponge = NY_S_T
    ENDIF

    CU1X_sponge(:) = 0.d0
!    CU1X_sponge = 0.d0
    CU2X_sponge = 0.d0
    CU3X_sponge = 0.d0

    DO J = NY_sponge, JEND
      DO K = ZSTART, ZEND
	CU1X_sponge(K) = CU1X_sponge(K) + CU1X(0,K,J)
!	CU1X_sponge = CU1X_sponge + CU1X(0,K,J)
	CU2X_sponge = CU2X_sponge + CU2X(0,K,J)
	CU3X_sponge = CU3X_sponge + CU3X(0,K,J)
      ENDDO
    ENDDO

!    CU1X_sponge = CU1X_sponge / (NZ * (NY - NY_sponge +1))
    CU1X_sponge(:) = CU1X_sponge(:) / (JEND - NY_sponge +1)
    CU2X_sponge = CU2X_sponge / ((ZEND -ZSTART +1) * (JEND - NY_sponge +1))
    CU3X_sponge = CU3X_sponge / ((ZEND -ZSTART +1) * (JEND - NY_sponge +1))

    CU1X_target = V_BC_YMAX_C1 !CU1X_sponge(:)
    CU2X_target = U_BC_YMAX_C1 
    CU3X_target = WBULK0 !W_BC_YMAX_C1
    
! Damp mean flow
	!TOP AND BOTTOM 
    IF (RANK .EQ. 0) THEN
    ! in case of dynamic sponge, top of atmospher would be sponged when zonal velocity is higher than the specified threshold * target velocity
      IF (DYNAMIC_SPONE) THEN
       IF (ABS(CU3X_sponge-CU3X_target).GT.(CU3X_target*Dynamic_sponge_threshold)) THEN
	  IF (RK_STEP.EQ.1)  WRITE(6,*) "**WARNING: Dynamic Sponeg is used. U_sponge=", CU3X_sponge, "U_target=", CU3X_target, "threshold=", Dynamic_sponge_threshold*CU3X_target
	  DO j=jstart,jend
	    DO k=zstart,zend
		CF1X(0,k,j)=CF1X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU1X(0,k,j)-CU1X_target) * INT_JACOB(K,J) 
!		CF1X(0,k,j)=CF1X(0,k,j) - CU1X_target * (CU3X(0,k,j)/CU3X_target) * INT_JACOB(K,J) 
		CF2X(0,k,j)=CF2X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU2X(0,k,j)-CU2X_target) * INT_JACOB(K,J) 
		CF3X(0,k,j)=CF3X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU3X(0,k,j)-CU3X_target) * INT_JACOB(K,J)  
		CFTHX(0,k,j,N)=CFTHX(0,k,j,N) - Sponge_Amp * SPONGE_SIGMA(j) * (CTHX(0,k,j,N)-0.d0)*INT_JACOB(K,J)
	    ENDDO
	  ENDDO
	ELSE
	  ! Forcing at top of atmospher to make the Meridonal componet to be zero
	  DO j=jstart,jend
	    DO k=zstart,zend
!		CF1X(0,k,j)=CF1X(0,k,j) - CU1X_target *(CU3X(0,k,j)/CU3X_target) *INT_JACOB(K,J) 
                CF2X(0,k,j)=CF2X(0,k,j) - Sponge_Amp *SPONGE_SIGMA(j)*(CU2X(0,k,j)-CU2X_target) * INT_JACOB(K,J)
                CF1X(0,k,j)=CF1X(0,k,j) - Sponge_Amp *SPONGE_SIGMA(j)*(CU1X(0,k,j)-CU1X_target) * INT_JACOB(K,J)
	    ENDDO
	  ENDDO
	ENDIF  
      ELSEIF (.NOT. DYNAMIC_SPONE) THEN
	
	DO j=jstart,jend
	  DO k=zstart,zend
	    IF (F_TYPE == 1) THEN
	      CF1X(0,k,j)=CF1X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU1X(0,k,j)-CU1X_target) * INT_JACOB(K,J) 
	     ! CF1X(0,k,j)=CF1X(0,k,j) - CU1X_target(k) * INT_JACOB(K,J) 
	      CF2X(0,k,j)=CF2X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU2X(0,k,j)-CU2X_target) * INT_JACOB(K,J) 
	      CF3X(0,k,j)=CF3X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU3X(0,k,j)-CU3X_target) * INT_JACOB(K,J)  
	      CFTHX(0,k,j,N)=CFTHX(0,k,j,N) - Sponge_Amp * SPONGE_SIGMA(j) * (CTHX(0,k,j,N)-0.d0)*INT_JACOB(K,J)
	    ELSEIF (F_TYPE == 0) THEN
	      CF1X(0,k,j)=CF1X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU1X(0,k,j)-CU1X_target) * INT_JACOB(K,J) 
	     ! CF1X(0,k,j)=CF1X(0,k,j) - CU1X_target(k) * INT_JACOB(K,J) 
	      CF2X(0,k,j)=CF2X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU2X(0,k,j)-CU2X_target) * INT_JACOB(K,J) 
	      CF3X(0,k,j)=CF3X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU3X(0,k,j)-CU3X_target) * INT_JACOB(K,J)  
	      CFTHX(0,k,j,N)=CFTHX(0,k,j,N) - Sponge_Amp * SPONGE_SIGMA(j) * (CTHX(0,k,j,N)-0.d0)*INT_JACOB(K,J)
	    ENDIF
	  ENDDO
	ENDDO
	
      ENDIF
! In the case of doubly periodic case, mean of vertical velocity is zero, therefore 
! In the nuetral case it works fine, so I just did it for stratified cases
!      IF ( V_BC_ZMAX.EQ.6 .AND. U_BC_ZMAX.EQ.6 .AND. W_BC_ZMAX.EQ.6 .AND. RI_TAU(1).NE.0) THEN
!       CF2X(0,:,:)=CF2X(0,:,:)-Sponge_Amp*(CU2X(0,:,:)-CU2X_target)*INT_JACOB(K,J)
!       IF (MOD(2*TIME_STEP,SAVE_STATS_INT).EQ.0) CU2X(0,:,:)=0.0d0
!      ENDIF

      IF ( V_BC_ZMAX.EQ.6 .AND. U_BC_ZMAX.EQ.6 .AND. W_BC_ZMAX.EQ.6 .AND. RI_TAU(1).NE.0) THEN
       CF2X(0,:,:)=CF2X(0,:,:)-Sponge_Amp*(CU2X(0,:,:)-CU2X_target)*INT_JACOB(K,J)
       IF (MOD(2*TIME_STEP,SAVE_STATS_INT).EQ.0) THEN
        CU2X(0,:,:)=0.0d0
!        DO J=0,NY+1
!         tmp=SUM(CPX(0,1:NZ,J))/dble(NZ);
!         CPX(0,:,J)=tmp;
!        ENDDO
       ENDIF
      ENDIF  

	!INFLOW OUTFLOW
	DO j=jstart,jend
	  DO k=zstart,zend
	    CF1X(0,k,j)=CF1X(0,k,j)- 2*Sponge_Amp * SPONGE_SIGMA_OUT(k)*(CU1X(0,k,j)-0.d0) * INT_JACOB(K,J)
	    CF2X(0,k,j)=CF2X(0,k,j)- 2*Sponge_Amp * SPONGE_SIGMA_OUT(k)*(CU2X(0,k,j)-0.d0) * INT_JACOB(K,J)
            CF3X(0,k,j)=CF3X(0,k,j)- 2*Sponge_Amp * SPONGE_SIGMA_OUT(k)*(CU3X(0,k,j)-CU3X_target) * INT_JACOB(K,J)
	    CFTHX(0,k,j,N)=CFTHX(0,k,j,N) - 2*Sponge_Amp * SPONGE_SIGMA_OUT(k) * (CTHX(0,k,j,N)-0.d0)*INT_JACOB(K,J)
	  ENDDO
	ENDDO
      
    ENDIF

  RETURN
END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE sponge_th_2(N)
! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state for the temperature field
! The intention is to allow an open boundary
  USE ntypes
  USE Domain
  USE Grid
  !USE Fft_var, only : NKX
  USE TIME_STEP_VAR
  USE run_variable
  USE mg_vari, ONLY : INIT_FLAG
  USE mpi_var  

    IMPLICIT NONE 

    ! The following variables will store the background state
    INTEGER i,j,k,N, NY_sponge
    REAL*8 Sponge_Amp, Dynamic_sponge_threshold
    REAL*8 CU1X_sponge, CU2X_sponge, CU3X_sponge
    REAL*8 CU1X_target, CU2X_target, CU3X_target    
    LOGICAL DYNAMIC_SPONE
    
! Dynamic sponge, applies the spone whenever the vertically averaged zonal velocity is higher than a specified value 
    DYNAMIC_SPONE = .TRUE. 
    
!    Sponge_Amp = 1.0d0
    Sponge_Amp = 0.1d0/DELTA_T
    Dynamic_sponge_threshold = 0.0125d0     !2.5% tollerance for dynamic sponge
    
! Damp fluctuation flowfield
	!TOP AND BOTTOM 
    DO j=jstart,jend
      DO k=zstart,zend
	DO i=1,NX2P
	    CF1X(i,k,j)=CF1X(i,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU1X(i,k,j)-0.d0)* INT_JACOB(K,J)  
	    CF2X(i,k,j)=CF2X(i,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU2X(i,k,j)-0.d0)* INT_JACOB(K,J)
	    CF3X(i,k,j)=CF3X(i,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU3X(i,k,j)-0.d0)* INT_JACOB(K,J)
	    CFTHX(i,k,j,N)=CFTHX(i,k,j,N) - Sponge_Amp * SPONGE_SIGMA(j) * (CTHX(i,k,j,N)-0.d0)*INT_JACOB(K,J)
	ENDDO
      ENDDO
    ENDDO


!INFLOW OUTFLOW
    DO j=jstart,jend
      DO k=zstart,zend
	DO i=1,NX2P
	    CF1X(i,k,j)=CF1X(i,k,j)- Sponge_Amp * SPONGE_SIGMA_OUT(k)*(CU1X(i,k,j)-0.d0) * INT_JACOB(K,J)
	    CF2X(i,k,j)=CF2X(i,k,j)- Sponge_Amp * SPONGE_SIGMA_OUT(k)*(CU2X(i,k,j)-0.d0) * INT_JACOB(K,J)
	    CF3X(i,k,j)=CF3X(i,k,j)- Sponge_Amp * SPONGE_SIGMA_OUT(k)*(CU3X(i,k,j)-0.d0) * INT_JACOB(K,J)
	    CFTHX(i,k,j,N)=CFTHX(i,k,j,N) - Sponge_Amp * SPONGE_SIGMA_OUT(k) * (CTHX(i,k,j,N)-0.d0)*INT_JACOB(K,J)
	ENDDO
      ENDDO
    ENDDO

! Appling equavalent pressure gradient in span-wise direction at top of atmospher! the span-wise velocity has been foreced to be zero in top of atmospher!     
    IF (NY_S_T .GT. NY) THEN
    ! deafult value
      NY_sponge = INT (2 * NY/3 )
    ELSE
    ! specified value
      NY_sponge = NY_S_T
    ENDIF

    CU1X_sponge = 0.d0
    CU2X_sponge = 0.d0
    CU3X_sponge = 0.d0

    DO J = NY_sponge, JEND
      DO K = ZSTART, ZEND
	CU1X_sponge = CU1X_sponge + CU1X(0,K,J)
	CU2X_sponge = CU2X_sponge + CU2X(0,K,J)
	CU3X_sponge = CU3X_sponge + CU3X(0,K,J)
      ENDDO
    ENDDO

    CU1X_sponge = CU1X_sponge / (NZ * (NY - NY_sponge +1))
    CU2X_sponge = CU2X_sponge / ((ZEND -ZSTART +1) * (JEND - NY_sponge +1))
    CU3X_sponge = CU3X_sponge / ((ZEND -ZSTART +1) * (JEND - NY_sponge +1))

    CU1X_target = U_BC_YMAX_C1
    CU2X_target = CU2X_sponge
    CU3X_target = W_BC_YMAX_C1
    
! Damp mean flow
	!TOP AND BOTTOM 
    IF (RANK .EQ. 0) THEN
    ! in case of dynamic sponge, top of atmospher would be sponged when zonal velocity is higher than the specified threshold * target velocity
      IF (DYNAMIC_SPONE) THEN
       IF (ABS(CU3X_sponge-CU3X_target).GT.(CU3X_target*Dynamic_sponge_threshold)) THEN
	  WRITE(6,*) "**WARNING: Dynamic Sponeg is used. U_sponge=", CU3X_sponge, "U_target=", CU3X_target, "threshold=", Dynamic_sponge_threshold
	  DO j=jstart,jend
	    DO k=zstart,zend
		CF1X(0,k,j)=CF1X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU1X(0,k,j)-CU1X_target) * INT_JACOB(K,J) 
		CF2X(0,k,j)=CF2X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU2X(0,k,j)-CU2X_target) * INT_JACOB(K,J) 
		CF3X(0,k,j)=CF3X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU3X(0,k,j)-CU3X_target) * INT_JACOB(K,J)  
		!CFTHX(0,k,j,N)=CFTHX(0,k,j,N) - Sponge_Amp * SPONGE_SIGMA(j) * (CTHX(0,k,j,N)-0.d0)*INT_JACOB(K,J)
	    ENDDO
	  ENDDO
	ELSEIF (ABS(CU1X_sponge-CU1X_target).GT.(CU1X_target*Dynamic_sponge_threshold)) THEN
	WRITE(6,*) "**WARNING: Dynamic Sponeg is used. V_sponge=", CU1X_sponge, "V_target=", CU1X_target, "threshold=", Dynamic_sponge_threshold
	  DO j=jstart,jend
	    DO k=zstart,zend
		CF1X(0,k,j)=CF1X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU1X(0,k,j)-CU1X_target) * INT_JACOB(K,J) 
		CF2X(0,k,j)=CF2X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU2X(0,k,j)-CU2X_target) * INT_JACOB(K,J) 
		CF3X(0,k,j)=CF3X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU3X(0,k,j)-CU3X_target) * INT_JACOB(K,J)   
	    ENDDO
	  ENDDO
	ENDIF  
      ELSEIF (.NOT. DYNAMIC_SPONE) THEN
	
	DO j=jstart,jend
	  DO k=zstart,zend
	    IF (F_TYPE == 1) THEN
	      CF1X(0,k,j)=CF1X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU1X(0,k,j)-CU1X_target) * INT_JACOB(K,J) 
	      CF2X(0,k,j)=CF2X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU2X(0,k,j)-CU2X_target) * INT_JACOB(K,J) 
	      CF3X(0,k,j)=CF3X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU3X(0,k,j)-CU3X_target) * INT_JACOB(K,J)  
	      !CFTHX(0,k,j,N)=CFTHX(0,k,j,N) - Sponge_Amp * SPONGE_SIGMA(j) * (CTHX(0,k,j,N)-0.d0)*INT_JACOB(K,J)
	    ELSEIF (F_TYPE == 0) THEN
	      CF1X(0,k,j)=CF1X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU1X(0,k,j)-CU1X_target) * INT_JACOB(K,J) 
	      !CF2X(0,k,j)=CF2X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU2X(0,k,j)-CU2X_target) * INT_JACOB(K,J) 
	      CF3X(0,k,j)=CF3X(0,k,j) - Sponge_Amp * SPONGE_SIGMA(j)*(CU3X(0,k,j)-CU3X_target) * INT_JACOB(K,J)  
	      !CFTHX(0,k,j,N)=CFTHX(0,k,j,N) - Sponge_Amp * SPONGE_SIGMA(j) * (CTHX(0,k,j,N)-0.d0)*INT_JACOB(K,J)
	    ENDIF
	  ENDDO
	ENDDO
	
      ENDIF
           
      !INFLOW OUTFLOW
      DO j=jstart,jend
	DO k=zstart,zend
	  CF1X(0,k,j)=CF1X(0,k,j)- Sponge_Amp * SPONGE_SIGMA_OUT(k)*(CU1X(0,k,j)-0.d0) * INT_JACOB(K,J)
	  CF2X(0,k,j)=CF2X(0,k,j)- Sponge_Amp * SPONGE_SIGMA_OUT(k)*(CU2X(0,k,j)-0.d0) * INT_JACOB(K,J)
	  CFTHX(0,k,j,N)=CFTHX(0,k,j,N) - Sponge_Amp * SPONGE_SIGMA_OUT(k) * (CTHX(0,k,j,N)-0.d0)*INT_JACOB(K,J)
	ENDDO
      ENDDO
      
    ENDIF

  RETURN
END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE sponge_vel
! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state
! The intention is to allow an open boundary
  USE ntypes
  USE Domain
  USE Grid
  !USE Fft_var, only : NKX
  USE TIME_STEP_VAR
  USE run_variable
  USE mg_vari, 		ONLY : INIT_FLAG
  USE mpi_var      
  
    IMPLICIT NONE 

! The following variables will store the background state
    REAL*8 Sponge_Amp
    INTEGER j,k

    Sponge_Amp =	0.067d0/DELTA_T

    IF (rank .eq. 0) THEN
      DO j=jstart,jend
	DO k=zstart,zend
	  CRTHX(0,K,J,1) = CRTHX(0,K,J,1)/(1.d0 + Sponge_Amp * SPONGE_SIGMA_OUT(K))
	ENDDO
      ENDDO
    ENDIF


  RETURN
END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE courant_curv
! This subroutine sets the timestep based on the specified CFL number
! The subroutine should be called with the velocity in physical space

  USE Domain
  USE Grid
  !USE Fft_var, only : NKX
  USE TIME_STEP_VAR
  USE run_variable
  USE mg_vari, 		ONLY : INIT_FLAG

    IMPLICIT NONE 

    REAL*8 dt
    REAL*8 dt_x,dt_y,dt_z,min_x,dt_dif
    INTEGER i,j,k

! Set the initial dt to some arbitrary large number
    dt=999.d0

!    dt based on Diffusion no 
    min_x = LX/dble(NX)       
    dt_dif    = DFN*min_x**2/NU

    DO j=1,NY
      DO k=1,NZ
	DO i=0,NXP
      ! C      U2 required to define at bottom cell face 	  
	  U2bX(I,K,J) = 0.5d0 * CJOB_22(K,J,1) * ( U2X(I,K,J)+U2X(I,K,J-1) ) &
		      + 0.5d0 * CJOB_21(K,J,1) * ( U3X(I,K,J)+U3X(I,K,J-1) )
      ! C      U3 required to define at side cell face     	    
	  U3bX(I,K,J) = 0.5d0 * CJOB_11(K,J,2) * ( U3X(I,K,J)+U3X(I,K-1,J) ) &
		      + 0.5d0 * CJOB_12(K,J,2) * ( U2X(I,K,J)+U2X(I,K-1,J) )

	  dt_x = cfl * dx(i) / abs( U1X(i,k,j) )
	  dt_y = cfl * INT_JACOB(K,J) / abs( U2bX(i,k,j) )
	  dt_z = cfl * INT_JACOB(K,J) / abs( U3bX(i,k,j) )
	  dt = min(dt,dt_x,dt_y,dt_z)
	ENDDO
      ENDDO
    ENDDO


    IF (dt.le.0) THEN
      WRITE(*,*) 'Error: dt<=0 in courant, dt=', dt
! 		Set DELTA_T to some small default value
!    	DELTA_T=0.001d0
    ELSE IF (dt.ge.DELTA_T) THEN
      WRITE(*,*) '**WARNING: dt=',dt,' capped at ', DELTA_T, ' **'
    ELSE
      
      DELTA_T = dt
      H_BAR(1) = DELTA_T * (8.d0/15.d0)
      H_BAR(2) = DELTA_T * (2.d0/15.d0)
      H_BAR(3) = DELTA_T * (5.d0/15.d0)
      
    ENDIF

  RETURN
END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE courant_curv_mpi
! This subroutine sets the timestep based on the specified CFL number
! The subroutine should be called with the velocity in physical space

  USE Domain
  USE Grid
  !USE Fft_var, only : NKX
  USE TIME_STEP_VAR
  USE run_variable
  USE mg_vari, 		ONLY : INIT_FLAG
  USE mpi_var, 		ONLY : rank      

    IMPLICIT NONE 

    REAL*8 dt
    REAL*8 dt_x, dt_y, dt_z, min_x, dt_dif
    INTEGER i, j, k

! Set the initial dt to some arbitrary large number
    dt = 999.d0

!    dt based on Diffusion no 
    min_x = LX/dble(NX)       
    dt_dif = DFN * min_x**2/NU
    
    Do J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO i=0,NXP
! C      U2 required to define at bottom cell face 	  
	  U2bX(I,K,J) = 0.5d0 * CJOB_22(K,J,1) * ( U2X(I,K,J)+U2X(I,K,J-1) ) &
		      + 0.5d0 * CJOB_21(K,J,1) * ( U3X(I,K,J)+U3X(I,K,J-1) )
! C      U3 required to define at side cell face     	    
	  U3bX(I,K,J) = 0.5d0 * CJOB_11(K,J,2) * ( U3X(I,K,J)+U3X(I,K-1,J) ) &
		      + 0.5d0 * CJOB_12(K,J,2) * ( U2X(I,K,J)+U2X(I,K-1,J) )
	  
	  dt_x = cfl * dx(i) / abs( U1X(i,k,j) )
	  dt_y = cfl * INT_JACOB(K,J) / abs( U2bX(i,k,j) )
	  dt_z = cfl * INT_JACOB(K,J) / abs( U3bX(i,k,j) )

!          IF (dt_x .LT. dt)   write(6,*) "dtx= ", dt_x, i,k,j, rank
!          IF (dt_y .LT. dt)   write(6,*) "dty= ", dt_y, i,k,j, rank
!          IF (dt_z .LT. dt)   write(6,*) "dtz= ", dt_z, i,k,j, rank

  	 dt = min(dt,dt_x,dt_y,dt_z)
	  
	ENDDO
      ENDDO
    ENDDO

    CALL MPI_COURANT(dt)

     IF (dt.le.0) THEN
       WRITE(*,*) 'Error: dt<=0 in courant, dt=', dt
     ELSE IF (  dt .ge. DELTA_T_SET  )  THEN
      IF (rank .eq. 0) THEN  
	WRITE(*,*) 'WARNING: dt=',dt,' capped at ', DELTA_T_SET
      ENDIF
       DELTA_T = DELTA_T_SET
       H_BAR(1) = DELTA_T * (8.d0/15.d0)
       H_BAR(2) = DELTA_T * (2.d0/15.d0)
       H_BAR(3) = DELTA_T * (5.d0/15.d0)
     ELSE
	DELTA_T = dt
	H_BAR(1) = DELTA_T * (8.d0/15.d0)
	H_BAR(2) = DELTA_T * (2.d0/15.d0)
	H_BAR(3) = DELTA_T * (5.d0/15.d0)
     ENDIF

  RETURN
END


! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE VEL_PROFILE
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
  USE mg_vari, 		ONLY : INIT_FLAG
!  USE, INTRINSIC :: IEEE_ARITHMETIC
	
    IMPLICIT NONE 

    INTEGER I,J,K 
    
    REAL*8 phi_int(0:NZ+1,0:NY+1)       
    
    REAL*8  div,Q,H1, H2,sigma,omega_0,z1,y1 

111   format(2f16.12)
      
    DO J=0,NY+1
      DO K=0,NZ+1
	  U3_BAR(K,J) = 0.d0
	  U2_BAR(K,J) = 0.d0
	  phi_int(k,j)= 0.d0
	  u2bx(:,k,j)  = 0.d0
	  u3bx(:,k,j)  = 0.d0
      ENDDO
    ENDDO
     
    H1 = (ypoint(1,NY) - ypoint(1,1))
    Q =  0.d0 !-(PX0/(12.d0*NU))*H1**3.0
    H2 = (ypoint(NZ,NY) - ypoint(NZ,1))
            

    DO J=1,NY
      DO K=0,1
	U3_BAR(K,J) = 6.d0 * (Q / H1**3.0) * ( (ypoint(1,NY)+ypoint(1,1) &
		  - ypoint(1,J) )*ypoint(1,J) - ypoint(1,NY)*ypoint(1,1) )
      ENDDO
      
      DO K=NZ-1,NZ+1
	U3_BAR(K,J)=6.d0*(Q/H2**3.0)*( (ypoint(NZ,NY)+ypoint(NZ,1) &
		  - ypoint(NZ,J) )*ypoint(NZ,J) - ypoint(NZ,NY)*ypoint(NZ,1) )
      ENDDO
    ENDDO      
       
    DO J=1,NY
      DO K=2,NZ-2
	H1 = (ypoint(k,NY) - ypoint(k,1))
	U3_BAR(K,J) = 6.d0 * (Q/H1**3.0) * ( (ypoint(k,NY)+ypoint(k,1) &
		    - ypoint(k,J) )*ypoint(k,J) - ypoint(k,NY)*ypoint(k,1) )  
      ENDDO
    ENDDO     
 
      z1 = 3.d0 
      y1 = 1.5d0 
      omega_0 = 0.5d0 
      sigma = 0.3d0	 

! c       DO J=2,NY-1
! c       DO K=2,NZ+1
! c        din1 = (ypoint(k,j)-y1)**2.0 + (xpoint(k,j)-z1)**2.0 ;
! c        U2_BAR(K,J)=U2_BAR(K,J)-(sigma**2.0*omega_0*(xpoint(k,j)-z1))
! c     &   *( 1.0 - exp(-(din1/(2.0*sigma**2.0))) )/din1
! c        U3_BAR(K,J)=U3_BAR(K,J)+(sigma**2.0*omega_0*(ypoint(k,j)-y1))
! c     &   *( 1.0 - exp(-(din1/(2.0*sigma**2.0))) )/din1
! 
! c       ENDDO
! c       ENDDO

    DO j=1,NY
      DO k=1,NZ
	DO i=0,NXP
! C      U2 required to define at bottom cell face        
	  U2bx(I,K,J) = 0.5d0 * CJOB_22(K,J,1) * ( U2_bar(K,J)+U2_bar(K,J-1) ) &
		      + 0.5d0 * CJOB_21(K,J,1) * ( U3_bar(K,J)+U3_bar(K,J-1) )
! C      U3 required to define at side cell face          
	  U3bx(I,K,J) = 0.5d0 * CJOB_11(K,J,2) * ( U3_bar(K,J)+U3_bar(K-1,J) ) &
		      + 0.5d0 * CJOB_12(K,J,2) * ( U2_bar(K,J)+U2_bar(K-1,J) )

	ENDDO
      ENDDO
    ENDDO
       
       
    DO J=1,NY  
      DO i=0,NXP
	U3bx(I,0,J) = U3bx(I,1,J)
      ENDDO
    ENDDO
      
      
!       
! c      call SOR (phi_int,rhs_w)
!       
! c      DO J=2,NY-1
! c        DO K=2,NZ
! c         U2_bar(K,J)=0.5*(CJOB_22(K,J+1,1)*(phi_int(K,J)+phi_int(K,J+1))
! c     &             - CJOB_22(K,J,1)*(phi_int(K,J)+phi_int(K,J-1)))
! c     &              /INT_JACOB(K,J)
! c     &             - 0.5*(CJOB_12(K+1,J,2)*(phi_int(K+1,J)+phi_int(K,J))
! c     &             - CJOB_12(K,J,2)*(phi_int(K,J)+phi_int(K-1,J)))
! c     &              /INT_JACOB(K,J)
! 
! c         U3_bar(K,J)=0.5*(CJOB_11(K+1,J,2)*(phi_int(K,J)+phi_int(K+1,J))
! c     &             - CJOB_11(K,J,2)*(phi_int(K,J)+phi_int(K-1,J)))
! c     &               /INT_JACOB(K,J)
! c     &             - 0.5*(CJOB_21(K,J+1,1)*(phi_int(K,J+1)+phi_int(K,J))
! c     &             - CJOB_21(K,J,1)*(phi_int(K,J)+phi_int(K,J-1)))
! c     &              /INT_JACOB(K,J)
! c        END DO
! c      END DO

    DO j=1,NY
      DO k=1,NZ
	DO i=0,NXP
! C      U2 required to define at bottom cell face        
	  U2bx(I,K,J) = 0.5d0 * CJOB_22(K,J,1) * ( U2_bar(K,J)+U2_bar(K,J-1) ) &
		      + 0.5d0 * CJOB_21(K,J,1) * ( U3_bar(K,J)+U3_bar(K,J-1) )
! C      U3 required to define at side cell face           
	  U3bx(I,K,J) = 0.5d0 * CJOB_11(K,J,2) * ( U3_bar(K,J)+U3_bar(K-1,J) ) &
		      + 0.5d0 * CJOB_12(K,J,2) * ( U2_bar(K,J)+U2_bar(K-1,J) )

	ENDDO
      ENDDO
    ENDDO

! c       DO J=1,NY
! c        DO K=1,NZ
! c          DO I=0,NXM
! c            U2b(I,K,J)= GMAT_22(K,J,1)*(phi_int(K,J)-phi_int(K,J-1))
! c     &          + 0.25*GMAT_12(K,J,1)*(phi_int(K+1,J)+phi_int(K+1,J-1)
! c     &             - phi_int(K-1,J)-phi_int(K-1,J-1))
! c            U3b(I,K,J) = GMAT_11(K,J,2)*(phi_int(K,J)-phi_int(K-1,J))
! c     &          + 0.25*GMAT_12(K,J,2)*(phi_int(K,J+1)+phi_int(K-1,J+1)
! c     &             - phi_int(K,J-1)-phi_int(K-1,J-1))
! c          END DO
! c        END DO
! c      END DO


    S1X(:,:,:)=0.d0

    DO J=2,NY-1
      DO K=2,NZ-1
	Do i=0,NXP
	  S1X(I,K,J)=  (U2bx(I,K,J+1) - U2bx(I,K,J)) &
		      + (U3bx(I,K+1,J) - U3bx(I,K,J))
	ENDDO
      ENDDO
    ENDDO

    DIV = 0.0D0
    
    DO J = 1,NY
      DO K = 1,NZ
	DO i=0,NXP
	    S1X(I,K,J)=  (U2bx(I,K,J+1) - U2bx(I,K,J)) &
			 + (U3bx(I,K+1,J) - U3bx(I,K,J))
	ENDDO
      ENDDO
    ENDDO

      DIV = 0.0D0
      
    DO J = 1,NY
      DO K = 1,NZ
	DO I = 1,NXP
	  DIV = DIV + DABS(S1X(I,K,J))
	ENDDO
      ENDDO
    ENDDO
    
    DIV = DIV / ( DBLE(NX) * DBLE(NY) * DBLE(NZ) )

      
    WRITE(6,*) "THE DIVERGENCE IS ", DIV, "at INITIAL FLOW field"

!    IF (IEEE_IS_NAN(DIV)) STOP 'DIV is a NaN'
!    IF (DELTA_T .EQ. 0.d0) STOP 'DELTA_T is zero'
 
  RETURN
END  

! C     input periodic boundary at any grid point (x,y) in the solution region
! C     (xa.le.x.le.xb,yc.le.y.le.yd) for mud2cr (MULTIGRID_MUDPACK)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE phi_bc(rhs_mg, nxa,nxb,nya,nyb)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  use ntypes
  use Domain
  use Grid
  use run_variable
  use mg_vari
  use pr_rem

      IMPLICIT NONE
      
      REAL*8  rhs_mg(0:nz+1,0:ny+1)
      INTEGER k,j
      INTEGER nxa, nxb, nya, nyb

      nxa = MUDPACK_BC(1)
      nxb = MUDPACK_BC(2)
      nya = MUDPACK_BC(3) ! =nyc
      nyb = MUDPACK_BC(4) !=nyd


      DO k=0,nz+1
	DO j=0,ny+1
	  rhs_mg(k,j) = - RHS_TMP(k,j) !cxx*pxx+cxy*pxy+cyy*pyy+cx*px+cy*py+ce*pe
	ENDDO
      ENDDO
      
   RETURN	
END

! C     input pde coefficients at any grid point (x,y) in the solution region
! C     (xa.le.x.le.xb,yc.le.y.le.yd) for mud2cr (MULTIGRID_MUDPACK)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE cofcr(x,y,cxx,cxy,cyy,cx,cy,ce)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

  USE ntypes
  USE Domain
  USE Grid
  USE TIME_STEP_VAR
  USE run_variable
  USE fft_var, 		ONLY: kx2_c
  USE les_chan_var
  USE forcing

      IMPLICIT NONE
      
      REAL*8 x,y,cxx,cxy,cyy,cx,cy,ce
      INTEGER k,j
      
      k = int(x) 
      j = int(y)
      
      
      cxx = 0.25d0 * ( GMAT_11(k,j,2) + GMAT_11(k+1,j,2) + GMAT_11(k,j,1) + GMAT_11(k,j+1,1))
      
      cxy = 0.5d0 * ( GMAT_12(k,j,2) + GMAT_12(k+1,j,2) + GMAT_12(k,j,1) + GMAT_12(k,j+1,1))
      
      cyy = 0.25d0 * ( GMAT_22(k,j,2) + GMAT_22(k+1,j,2) + GMAT_22(k,j,1) + GMAT_22(k,j+1,1))
      
      cx = GMAT_11(k+1,j,2)-GMAT_11(k,j,2)+GMAT_12(k,j+1,1)-GMAT_12(k,j,1)
      
      cy = GMAT_22(k,j+1,1)-GMAT_22(k,j,1)+GMAT_12(k+1,j,2)-GMAT_12(k,j,2)
      
      ce = -kx2_c * INT_JACOB(k,j)
      
  RETURN	
END
