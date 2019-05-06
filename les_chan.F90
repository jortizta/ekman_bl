SUBROUTINE les_chan
! This subroutine models the terms owing to the subgrid scale stress
! if the computation is to be treated as an LES not a DNS
!C This subroutine should be called when the velocity is in fourier space 
!C in the periodic directions, on output, the velocity will be 
!C in physical space.
!C It is assumed that the test filter and the LES filter are performed
!C by the same operation
!C On output S1 should contain |S| which may be used again in les_chan_th
!C if for the subgrid scalar dissipation

  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
  USE les_chan_var
  USE mpi_var, only : rank
    
    IMPLICIT NONE

    INTEGER i,j,k,ij
    CHARACTER*28 file_name, file_name_plt

! Variables for Dynamic Smagorinsky model:
    REAL*8 C_SMAG
    PARAMETER (C_SMAG=0.085d0)
      
    REAL*8 alpha,beta 
! Array to store the velocity index for each component of the strain rate tensor
    INTEGER U_index1(6)
    INTEGER U_index2(6)
    LOGICAL LES_IMPLICIT, dir_exists

! Here, alpha is the test/LES filter width ratio
    PARAMETER (alpha=2.44950)
! beta is the LES/grid filter width ratio
    PARAMETER (beta=1.d0)
    
    LES_IMPLICIT =.FALSE.		! has to be false for LES_MODEL_TYPE=1 and LES_MODEL_TYPE=2 


! Set the velocity index for each component of the stress tensor
    U_index1(1)=1
    U_index2(1)=1

    U_index1(2)=2
    U_index2(2)=2 

    U_index1(3)=3
    U_index2(3)=3

    U_index1(4)=1
    U_index2(4)=2

    U_index1(5)=1
    U_index2(5)=3

    U_index1(6)=2
    U_index2(6)=3


! When using a Near-wall model, DOn't use LES at the wall
!      IF ((U_BC_YMIN.EQ.3).or.(W_BC_YMIN.EQ.3)) then
!C         J1=2
!C      ELSE
!C         J1=JSTART
!C      END IF

!C      IF ((U_BC_YMAX.EQ.3).or.(W_BC_YMAX.EQ.3)) then
!C         J2=NY-1
!C      ELSE
!C         J2=JEND
!C      END IF


! First, for all models, apply boundary conditions to the velocity field
! (fill ghost cells) to ensure accurate calculation of gradients
!C Apply Boundary conditions to velocity field
    J1  = JSTART
    J2  = JEND
    J1i= 1
    J2e=NY+1


    IF (LES_MODEL_TYPE.EQ.1) THEN
!     Constant Smagorinsky model
! First, compute the rate of strain tensor S_ij
      IF (LES_IMPLICIT) STOP 'for LES_MODEL_TYPE 1, LES_IMPLICIT should be false, please correct it in les_ch.F90'
      CALL compute_strain_chan

! Compute |S| at GYF points, store in S2X
! Interpolation to GYF points is easy since by definition
! GYF points are exactly midway between neighboring GY points
      DO J=JSTART,JEND
	DO K=ZSTART,ZEND
	  DO I=0,NXP
	    S2X(I,K,J)=SQRT(                               &
			    2.d0 * Sij(I,K,J,1)**2.d0        & 
			    + 4.d0 * Sij(I,K,J,4)**2.d0       &
			    + 4.d0 * Sij(I,K,J,5)**2.d0       &
			    + 2.d0 * Sij(I,K,J,2)**2.d0       & 
			    + 4.d0 * Sij(I,K,J,6)**2.d0       &
			    + 2.d0 * Sij(I,K,J,3)**2.d0       &
			    )
	  ENDDO
	ENDDO
      ENDDO
! ExtEND |S| to ghost cells
      DO K=0,NZ+1
	DO I=0,NXP
	  S2X(I,K,0) = S2X(I,K,1)
	  S2X(I,K,NY+1) = S2X(I,K,NY)
	ENDDO
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      Storing for tke les
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO ij=1,6
          DO J=0,NY+1
            DO K=0,NZ+1
              DO I=0,NXP
                St_rij(I,K,J,IJ) = Sij(I,K,J,IJ)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

! Need to calculate Mean Strain Rate

        DO ij=1,6
          DO j=0,NY+1
            DO k=0,NZ+1
              Sij_mean(k,j,ij) = SUM(Sij(0:min(NXP,NXP_L),k,j,ij))/dble(NX)
            ENDDO
          ENDDO
          CALL MPI_COMBINE_STATS(Sij_mean(0,0,ij),NZ+2,NY+2)
        ENDDO



! Convert the velocities to physical space
      CALL REAL_FOURIER_TRANS_U1 (.FALSE.)
      CALL REAL_FOURIER_TRANS_U2 (.FALSE.)
      CALL REAL_FOURIER_TRANS_U3 (.FALSE.)

! Compute the filter lengthscale
      DO J=1,NY+1
	DO K=1,NZ+1
!	DELTA_Y(K,J) = (C_SMAG *(DX(1) * INT_JACOB(K,J))**(1.d0/3.d0))**(2.d0)
        DELTA_Y(K,J) = (beta*DX(1) * INT_JACOB(K,J))**(1.d0/3.d0)
	ENDDO
      ENDDO

! Get the eddy viscosity at GY points
! NU_T = (C_S^2 * DELTA^2)*|S|
      DO J=0,NY+1
        DO K=0,NZ+1
          DO I=0,NXP
            NU_T(I,K,J) = C_SMAG**2.d0 * DELTA_Y(K,J)**2.d0 * S2X(I,K,J)
          ENDDO
        ENDDO
      ENDDO

 !Now, compute Tau_ij=-2.0 NU_T*Sij      
      DO ij=1,6
	  IF ( ij.eq.1 ) THEN
    ! Here, Sij is defined at GYF points 
	    DO j=1,NY+1
	      DO k=1,NZ+1
		DO i=0,NXP
		  Sij(i,k,j,ij)= -2.d0*NU_T(I,K,J)*Sij(I,K,J,ij)
		ENDDO 
	      ENDDO
	    ENDDO
	  
	  ELSEIF (ij.EQ.2) THEN
    ! Sij(:,:,:,2) = du2/dy, but this term will be added implicitly through
    ! an eddy viscosity, so set it equal to zero here
	      DO j=1,NY+1
		DO k=1,NZ+1
		  DO i=0,NXP
		    Sij(i,k,j,ij) = 0.d0
		  ENDDO
		ENDDO
	      ENDDO
	      
	  ELSEIF (ij.EQ.3) THEN
      ! Sij(:,:,:,3) = du3/dz, but this term will be added implicitly through
      ! an eddy viscosity, so set it equal to zero here
	      DO j=0,NY+1
		DO k=1,NZ+1
		  DO i=0,NXP
		    Sij(i,k,j,ij) = 0.d0
		  ENDDO
		ENDDO
	      ENDDO
	    
	  ELSEIF (ij.EQ.4) THEN 
    ! Sij(:,:,:,4)=0.5*(dU1/dy + dU2/dx)
    ! But dU1/dy will be accounted for as an implicit eddy viscosity term,
    ! So, subtract if off from Sij here
	  
	    DO j=JSTART,JEND
	      DO k=ZSTART,ZEND
		DO i=0,NXP
		  Sij(i,k,j,ij)=-2.d0 *NU_T(I,K,J)*(Sij(i,k,j,ij)-              &
			      0.5* ( 0.5*CJOB_12(K+1,J,2)*(U1X(I,K,J)   + U1X(I,K+1,J))     &
			      - 0.5*CJOB_12(K,J,2)*(U1X(I,K,J)   + U1X(I,K-1,J))            &
			      + 0.5*CJOB_22(K,J+1,1)*(U1X(I,K,J) + U1X(I,K,J+1))            &
			      - 0.5*CJOB_22(K,J,1)*(U1X(I,K,J)   + U1X(I,K,J-1)))/          &
			      INT_JACOB(K,J)  )
		ENDDO
	      ENDDO
	    ENDDO
	  
    
	  ELSEIF (ij.EQ.5) THEN
    ! Sij(:,:,:,5)=0.5*(dU1/dz + dU3/dx)
    ! But dU1/dz will be accounted for as an implicit eddy viscosity term,
    ! So, subtract if off from Sij here
	    DO j=JSTART,JEND
	      DO k=ZSTART,ZEND
		DO i=0,NXP
		  Sij(i,k,j,ij)=-2.d0 * NU_T(I,K,J)*(Sij(i,k,j,ij)-0.5d0*      &
			      ( 0.5*CJOB_11(K+1,J,2)*(U1X(I,K,J)+ U1X(I,K+1,J))          &
			      - 0.5*CJOB_11(K,J,2)*(U1X(I,K,J)   + U1X(I,K-1,J))         &
			      + 0.5*CJOB_21(K,J+1,1)*(U1X(I,K,J) + U1X(I,K,J+1))         &
			      - 0.5*CJOB_21(K,J,1)*(U1X(I,K,J)   + U1X(I,K,J-1)) )/      & 
			    INT_JACOB(K,J)   )                         
		ENDDO
	      ENDDO
	    ENDDO
	    
	  ELSEIF (ij.EQ.6) THEN     
    ! Sij(:,:,:,6)=0.5*(dU3/dy + dU2/dz)
    ! But dU3/dy as well as  dU2/dz will be accounted for as an implicit eddy viscosity term,
    ! So, subtract if off from Sij here
	    DO j=JSTART,JEND
	      DO k=ZSTART,ZEND
		DO i=0,NXP
		  Sij(i,k,j,ij)=-2.d0 *NU_T(I,K,J)*(Sij(i,k,j,ij)-           &
				  0.5*( 0.5*CJOB_11(K+1,J,2)*(U2X(I,K,J)+ U2X(I,K+1,J))       &
				    - 0.5*CJOB_11(K,J,2)*(U2X(I,K,J)   + U2X(I,K-1,J))         &
				    + 0.5*CJOB_21(K,J+1,1)*(U2X(I,K,J) + U2X(I,K,J+1))         &
				    - 0.5*CJOB_21(K,J,1)*(U2X(I,K,J)   + U2X(I,K,J-1))) / INT_JACOB(K,J)  )
		ENDDO
	      ENDDO
	    ENDDO 
	  ENDIF
	      
      ENDDO   
      
! Now, compute |S|*S_ij, storing in Sij
! First compute at GYF points 
! 	DO J=JSTART,JEND
! 	  DO K=ZSTART,ZEND
! 	    DO I=0,NXP
! 	      Sij(I,K,J,1) = S2X(I,K,J) * Sij(I,K,J,1)
! 	      Sij(I,K,J,5) = S2X(I,K,J) * Sij(I,K,J,5)
!   ! CSij(:,:,:,2) is added through an implicit eddy viscosity
! 	      Sij(I,K,J,2) = 0.d0
! 	      Sij(I,K,J,3) = S2X(I,K,J) * Sij(I,K,J,3)
! 	    ENDDO
! 	  ENDDO
! 	ENDDO
      
! ! Now, compute at GY points, interpolating |S|
! 	DO J=1,NY
! 	  DO K=1,NZ
! 	    DO I=0,NXP
!   ! |S| interpolated to GY point 
! 	      !TEMP(I,K,J) = (S2X(I,K,J)*DYF(j-1)+S2X(I,K,J-1)*DYF(j)) /(2.d0*DY(j))
! 	      TEMP(I,K,J) = S2X(I,K,J)
!   ! The terms dU1/dy and dU3/dy in CSij(:,:,:,4) and CSij(:,:,:,6) respectively
!   ! are subtracted out from Sij here since they are treated implicitly
!   ! in eddy viscosity terms
! 	      !Sij(I,K,J,4) = TEMP(I,K,J) * ( Sij(I,K,J,4)-0.5*(CU1X(I,K,J)-CU1X(I,K,J-1))/DY(j) )
! 	      !Sij(I,K,J,6) = TEMP(I,K,J) * ( Sij(I,K,J,6)-0.5*(CU3X(I,K,J)-CU3X(I,K,J-1))/DY(j) )
! 	    ENDDO
! 	  ENDDO
! 	ENDDO


! We now have Tau_ij stored in Sij in Physical space
! Convert Tau_ij to Fourier space
	DO ij=1,6
	  S1X=Sij(:,:,:,ij)    
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
	  CSij(:,:,:,ij)=CS1X  
	ENDDO

! ! ! ! Now, compute TAU, store in the corresponging Sij
! ! ! 	DO J=1,NY
! ! ! 	  DO K=1,NZ
! ! ! 	    DO I=0,NX2P
! ! ! 	      CSij(I,K,J,1) = DELTA_Y(k,J) * CSij(I,K,J,1)
! ! ! 	      CSij(I,K,J,5) = DELTA_Y(k,J) * CSij(I,K,J,5)
! ! !   ! CSij(:,:,:,2) is added through an implicit eddy viscosity
! ! !   !            CSij(I,K,J,2)=DELTA_Y(K,J)*CSij(I,K,J,2)
! ! ! 	      CSij(I,K,J,3) = DELTA_Y(K,J) * CSij(I,K,J,3)
! ! ! 	    ENDDO
! ! ! 	  ENDDO
! ! ! 	ENDDO
! ! ! 	
! ! ! 	DO J=1,NY+1 
! ! ! 	  DO K=1,NZ+1
! ! ! 	    DO I=0,NX2P
! ! ! 	      CSij(I,K,J,4) = DELTA_Y(K,J) * CSij(I,K,J,4)
! ! ! 	      CSij(I,K,J,6) = DELTA_Y(K,J) * CSij(I,K,J,6)
! ! ! 	    ENDDO
! ! ! 	  ENDDO
! ! ! 	ENDDO

! tau_ij is now contained in CSij in Fourier space

    ELSEIF ((LES_MODEL_TYPE.EQ.2).or.(LES_MODEL_TYPE.eq.3)) THEN
! Here, use a dynamic smagorinsky model with or without scale similar part

! Compute the filter width
! 	DO J=0,NY
! ! At GYF points:
! 	  DELTA_YF(J)=(beta*DX(1)*DYF(J)*beta*DZ(1))**(1.d0/3.d0)
! !        DELTA_YF(J)=sqrt((beta*DX(1))**2.d0+(DYF(J)*2.d0)**2.d0
! !     &      +(beta*DZ(1))**2.d0)
! 	ENDDO

	DO J=0,NY+1
	  DO K=0,NZ+1
! At cell centered  points:
	    DELTA_Y(K,J) = (beta*DX(1) * INT_JACOB(K,J))**(1.d0/3.d0)
	  ENDDO
	ENDDO

! We need to calculate the components of C, the dynamic coefficient
	
! Compute the rate of strain tensor, store in Sij
	CALL compute_strain_chan
	
  
! Compute |S| , store in S2X, bcz S2X will be used for ffts and filter variables
      DO J=JSTART,JEND
	  DO K=ZSTART,ZEND
	    DO I=0,NXP
	      S2X(I,K,J)=SQRT(            	                                   &
			      2.d0*Sij(I,K,J,1)**2.d0                              &
			      + 4.d0*Sij(I,K,J,4)**2.d0                            &
			      + 4.d0*Sij(I,K,J,5)**2.d0                            &
			      + 2.d0*Sij(I,K,J,2)**2.d0                            &
			      + 4.d0*Sij(I,K,J,6)**2.d0                            &
			      + 2.d0*Sij(I,K,J,3)**2.d0                            &
			      )
	    ENDDO
	  ENDDO
	ENDDO
! ExtEND |S| to ghost cells
	DO K=0,NZ+1
	  DO I=0,NXP
	    S2X(I,K,0)    = S2X(I,K,1)
	    S2X(I,K,NY+1) = S2X(I,K,NY)
	  ENDDO
	ENDDO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      Storing for tke les
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DO ij=1,6
	  DO J=0,NY+1
	    DO K=0,NZ+1
	      DO I=0,NXP
		St_rij(I,K,J,ij) = Sij(I,K,J,ij)
	      ENDDO
	    ENDDO
	  ENDDO
	ENDDO      


! Need to calculate Mean Strain Rate

	DO ij=1,6 
	  DO j=0,NY+1
	    DO k=0,NZ+1
	      Sij_mean(k,j,ij) = SUM(Sij(0:min(NXP,NXP_L),k,j,ij))/dble(NX)
	    ENDDO
	  ENDDO
	  CALL MPI_COMBINE_STATS(Sij_mean(0,0,ij),NZ+2,NY+2)  
	ENDDO



	DO K=0,NZ+1
	  DO J=0,NY+1
	    U1_BAR_les(k,j) = CU1X(0,k,j)
	    U2_BAR_les(k,j) = CU2X(0,k,j)
	    U3_BAR_les(k,j) = CU3X(0,k,j)
	  ENDDO
	ENDDO

	CALL MPI_BCAST_REAL(U1_BAR_les,NZ+2,NY+2)
	CALL MPI_BCAST_REAL(U2_BAR_les,NZ+2,NY+2)
	CALL MPI_BCAST_REAL(U3_BAR_les,NZ+2,NY+2)     

! Convert Ui to physical space


	CALL REAL_FOURIER_TRANS_U1 (.FALSE.)
	CALL REAL_FOURIER_TRANS_U2 (.FALSE.)
	CALL REAL_FOURIER_TRANS_U3 (.FALSE.)

!      CALL FFT_X_TO_PHYSICAL(CU1,U1,0,NY+1,0,NZ+1)
!      CALL FFT_X_TO_PHYSICAL(CU2,U2,0,NY+1,0,NZ+1)
!      CALL FFT_X_TO_PHYSICAL(CU3,U3,0,NY+1,0,NZ+1)
	
! Apply the filter to the LES velocity and save
	DO j=0,NY+1
	  DO k=0,NZ+1
	    DO i=0,NXP
	      U_BAR_TIL(i,k,j,1) = U1X(i,k,j)
	      U_BAR_TIL(i,k,j,2) = U2X(i,k,j)
	      U_BAR_TIL(i,k,j,3) = U3X(i,k,j)
	    ENDDO
	  ENDDO
	ENDDO

	IF (les_model_type.eq.3) THEN
! storing in the U_2BAR
	  DO ij =1,3
	    DO j=0,NY+1
	      DO k=0,NZ+1
		DO i=0,NXP
		  U_2BAR(i,k,j,ij) = U_BAR_TIL(i,k,j,ij)
		ENDDO
	      ENDDO
	    ENDDO
	  ENDDO
	
	  DO i=1,3
! application of grid filter i.e. filter_type = 1
	    S1X = U_2BAR(:,:,:,i)         
	    CALL FILTER_VAR(1)
	    U_2BAR(:,:,:,i)=S1X
!        call les_filter_chan(U_2BAR(0,0,0,i),0,NY+1,1)
	  ENDDO
	ENDIF
    
	
  ! Now, filter the velocity
	DO i=1,3
  ! application of test filter i.e. filter_type = 2
	  S1X = U_BAR_TIL(:,:,:,i)         
	  CALL FILTER_VAR(2)
	  U_BAR_TIL(:,:,:,i)=S1X

  !        call les_filter_chan(U_BAR_TIL(0,0,0,i),0,NY+1,2)
	ENDDO
	

  ! Compute C_DYN only every x # of timesteps
	IF (((MOD(TIME_STEP,5).EQ.0).AND.(RK_STEP.eq.1)) .OR.FIRST_TIME) THEN

	  IF (rank .EQ. 0) THEN
	    WRITE(6,*) '##############################'
	    WRITE(6,*) 'C_dyn is recalculating' 
	    WRITE(6,*) '##############################'
	  ENDIF

	  CALL allocate_les_tmp
	  
    ! Filter |S| and store in S_2BAR
	  DO J=0,NY+1 
	    DO K=0,NZ+1
	      DO I=0,NXP
		S_2BAR(I,K,J) = S2X(I,K,J)
	      ENDDO
	    ENDDO
	  ENDDO
	  
    ! Test filtering operation filter type = 2
	  S1X = S_2BAR
	  CALL FILTER_VAR(2)
	  S_2BAR=S1X

    !      call les_filter_chan(S_2BAR,0,NY+1,2)

    ! Save a copy of the velocity which will be filtered twice more 

	  IF (les_model_type.EQ.3) THEN 
    ! Do only if a scale similar part is needed
	    DO ij=1,3  
	      DO j=0,NY+1
		DO k=0,NZ+1
		  DO i=0,NXP
		    U_4BAR(i,k,j,ij) = U_BAR_TIL(i,k,j,ij)
		  ENDDO
		ENDDO
	      ENDDO
	    ENDDO

	    DO i=1,3
    ! application of test and bar filter togather i.e. two filtering
    ! operation first filter type = 1, second filter type = 2
	      S1X = U_4BAR(:,:,:,i)
	      CALL FILTER_VAR(1)
	      CALL FILTER_VAR(2)
	      U_4BAR(:,:,:,i) = S1X
    !          call les_filter_chan(U_4BAR(0,0,0,i),0,NY+1,1)
    !          call les_filter_chan(U_4BAR(0,0,0,i),0,NY+1,2)
	    ENDDO
	  ENDIF


    ! Zero C_DYN
	  DO j=0,NY+1
	    DO k=0,NZ+1
	      C_DYN(k,j) = 0.d0
	    ENDDO
	  ENDDO

	  

    ! The prep. work is now DOne, we are ready to start the algorithm to compute
    ! the dynamic model coefficient

	  DO j =0,NY+1
	    DO k =0,NZ+1
	      denominator_sum(k,j) = 0.d0
	      numerator_sum(k,j)   = 0.d0
	    ENDDO
	  ENDDO
    
    ! Do over all non-repeating components of the stress tensor
	  DO ij=1,6

      ! Here        ij=1 -> l=1,m=1
      !             ij=2 -> l=2,m=2
      !             ij=3 -> l=3,m=3
      !             ij=4 -> l=1,m=2
      !             ij=5 -> l=1,m=3
      !             ij=6 -> l=2,m=3       
	
      ! Zero the numerator and denominator:
	      DO j=0,NY+1
		DO k=0,NZ+1
		  DO i=0,NXP
		    numerator(i,k,j)=0.d0
		    denominator(i,k,j)=0.d0
		  ENDDO
		ENDDO
	      ENDDO

      ! cross is used to multiply by two to include the contribution from the
      ! other symmetric term if we are dealing with a cross-term
	      IF (ij.le.3) THEN 
      ! We are computing a diagonal term
		cross=1.d0
	      ELSE
      ! We are computing a cross term
		cross=2.d0 
	      ENDIF 

      ! First, compute Mij

	      DO j=0,NY+1 
		DO k=0,NZ+1
		  DO i=0,NXP
	! Sij is defined at GYF points, no interpolation needed
		    temp(i,k,j) = Sij(i,k,j,ij)
		  ENDDO
		ENDDO
	      ENDDO
	    
      ! Filter temp test filter operation filter type = 2
	      S1X = temp
	      CALL FILTER_VAR(2)
	      temp = S1X
	      
      ! call les_filter_chan(temp,0,NY+1,2)
      ! Multiply by |S| filtered        
	      DO j=0,NY+1
		DO k=0,NZ+1
		  DO i=0,NXP
		    temp(i,k,j) = temp(i,k,j) * (alpha*DELTA_Y(k,j))**2.d0 * S_2BAR(i,k,j) 
		  ENDDO
		ENDDO
	      ENDDO
  
      ! Sij is used for Mij 
	      DO j=0,NY+1
		DO k=0,NZ+1 
		  DO i=0,NXP
		    Mij(i,k,j) = DELTA_Y(k,j)**2.d0*S2X(i,k,j) * Sij(i,k,j,ij)
		  ENDDO
		ENDDO
	      ENDDO

      ! Filter Mij test filter operation filter type = 2
	      S1X = Mij
	      CALL FILTER_VAR(2)
	      Mij=S1X 
      !       call les_filter_chan(Mij,0,NY+1,2)
      
      ! Add the second term of Mij stored in temp
	      DO j=0,NY+1
		DO k=0,NZ+1
		  DO i=0,NXP
		    Mij(i,k,j) = temp(i,k,j) - Mij(i,k,j)    
		  ENDDO
		ENDDO
	      ENDDO
      ! Now, compute Lij and add Lij*Mij to the numerator
      ! temp=Ui*Uj:
	    SELECT CASE (ij)
	    
	      CASE(1)
		DO j=0,NY+1
		  DO k=0,NZ+1
		    DO i=0,NXP
		      temp(i,k,j) = U1X(i,k,j)*U1X(i,k,j)
		    ENDDO
		  ENDDO 
		ENDDO
		
	      CASE(2)
		DO j=0,NY+1
		  DO k=0,NZ+1
		    DO i=0,NXP
		      temp(i,k,j) = U2X(i,k,j)*U2X(i,k,j)
		    ENDDO
		  ENDDO
		ENDDO
		
	      CASE(3)
		DO j=0,NY+1
		  DO k=0,NZ+1
		    DO i=0,NXP
		      temp(i,k,j) = U3X(i,k,j)*U3X(i,k,j)
		    ENDDO
		  ENDDO 
		ENDDO
		
	      CASE(4) 
		DO j=0,NY+1
		  DO k=0,NZ+1
		    DO i=0,NXP
		      temp(i,k,j) = U1X(i,k,j)*U2X(i,k,j)
		    ENDDO
		  ENDDO
		ENDDO
		
	      CASE(5)
		DO j=0,NY+1
		  DO k=0,NZ+1
		    DO i=0,NXP
		      temp(i,k,j) = U1X(i,k,j)*U3X(i,k,j)
		    ENDDO
		  ENDDO 
		ENDDO
		
	      CASE(6) 
		DO j=0,NY+1
		  DO k=0,NZ+1
		    DO i=0,NXP
		      temp(i,k,j) = U2X(i,k,j)*U3X(i,k,j)
		    ENDDO
		  ENDDO
		ENDDO
      !END select
	    ENDSELECT

	    DO j=0,NY+1
	      DO k=0,NZ+1
		  DO i=0,NXP
		      temp_1(i,k,j) = temp(i,k,j)
		  ENDDO
		ENDDO
	    ENDDO
      ! Filter temp: test filter operation i.e. filter type = 2
	      S1X = temp
	      CALL FILTER_VAR(2)
	      temp=S1X
      !       call les_filter_chan(temp,0,NY+1,2)
      ! Add Lij*Mij to numerator
	      DO j=0,NY+1
		DO k=0,NZ+1
		  DO i=0,NXP 
		    numerator(i,k,j) = Mij(i,k,j) * (temp(i,k,j) - U_BAR_TIL(i,k,j,U_index1(ij)) * U_BAR_TIL(i,k,j,U_index2(ij)))
		  ENDDO
		ENDDO
	      ENDDO

	    IF (LES_MODEL_TYPE.eq.3) THEN
      ! If mixed model, include the scale similar part  
      ! Add Nij*Mij to the numerator piece-by-piece
      ! Term3 

      ! similarly grid filter operation is DOne on temp i.e. filter type = 1
	      S1X = temp_1
	      CALL FILTER_VAR(1)
      !         call les_filter_chan(temp_1,0,NY+1,1)
      ! Filter temp_1 test filter operation filter type = 2
	      CALL FILTER_VAR(2)
	      temp_1=S1X
      !         call les_filter_chan(temp_1,0,NY+1,2)
	      DO j=0,NY+1
		DO k=0,NZ+1
		  DO i=0,NXP
		    numerator(i,k,j) = numerator(i,k,j)+temp_1(i,k,j)*Mij(i,k,j)
		  ENDDO
		ENDDO
	      ENDDO
      ! Term 4
	      DO j=0,NY+1
		DO k=0,NZ+1
		  DO i=0,NXP
		    temp(i,k,j)=U_2BAR(i,k,j,U_index1(ij)) * U_2BAR(i,k,j,U_index2(ij))      
		    temp_1(i,k,j) = U_BAR_TIL(i,k,j,U_index1(ij)) * U_BAR_TIL(i,k,j,U_index2(ij))
		  ENDDO
		ENDDO
	      ENDDO
      ! Filter temp test filter operation filter type = 2
	      S1X = temp
	      CALL FILTER_VAR(2)
	      temp=S1X 
      !         call les_filter_chan(temp,0,NY+1,2)
	      DO j=0,NY+1
		DO k=0,NZ+1
		  DO i=0,NXP
		    numerator(i,k,j) = numerator(i,k,j)-temp(i,k,j)*Mij(i,k,j)
		  ENDDO
		ENDDO
	      ENDDO 
      ! Terms 1 and 2
      ! Filter temp_1 grid filter operation filter type = 1
	      S1X = temp_1
	      call FILTER_VAR(1)
	      call FILTER_VAR(2)
	      temp_1=S1X
      !         call les_filter_chan(temp_1,0,NY+1,1)
      ! Filter temp test filter operation filter type = 2
      !         call les_filter_chan(temp_1,0,NY+1,2)
		DO j=0,NY+1
		  DO k=0,NZ+1
		    DO i=0,NXP
		      numerator(i,k,j) = numerator(i,k,j)-(temp_1(i,k,j) - U_4BAR(i,k,j,U_index1(ij)) * U_4BAR(i,k,j,U_index2(ij))) * Mij(i,k,j)
		    ENDDO
		  ENDDO
		ENDDO
	    ENDIF
      ! Now, the denominator for this ij
	    DO j=0,NY+1
	      DO k=0,NZ+1
		DO i=0,NXP
		  denominator(i,k,j) = Mij(i,k,j)*Mij(i,k,j)
		ENDDO
	      ENDDO
	    ENDDO

	      DO j=0,NY+1
		DO k=0,NZ+1
		  denominator_sum(k,j) = denominator_sum(k,j) + cross * SUM(denominator(0:min(NXP,NXP_L),k,j))
		  numerator_sum(k,j)   = numerator_sum(k,j)  + cross * SUM(numerator(0:min(NXP,NXP_L),k,j))
		ENDDO
	      ENDDO 

    ! End to ij
	  ENDDO


	  CALL MPI_COMBINE_STATS(denominator_sum,NZ+2,NY+2)
	  CALL MPI_COMBINE_STATS(numerator_sum,NZ+2,NY+2)

    ! Get plane average of numerator and denominator, add to C
    ! Note, since both the numerator and denominator are averaged, there
    ! is not need to divide the sum by NX*NZ

	  DO j=jstart,JEND
	    DO k=zstart,ZEND
	      IF (denominator_sum(k,j).ne.0.) THEN
		C_DYN(k,j) = -0.5d0 * numerator_sum(k,j)/denominator_sum(k,j)
	      ELSE 
		C_DYN(k,j) = 0.d0  
	      ENDIF
	    ENDDO
	  ENDDO

    ! We are now done with the dynamic procedure to calculate C

    ! If C_DYN < 0 at any level, set C_DYN=0 for numerical stability
	  DO j=0,NY+1
	    DO k=0,NZ+1
	      IF (C_DYN(k,j).LT.0) C_DYN(k,j)=0.d0
	    ENDDO
	  ENDDO


    ! At this point we have C_DYN and Sij
	  
	  CALL MPI_BCAST_REAL(C_DYN, NZ+2, NY+2)
    
	  
    ! End if compute C_DYN
	  
	  CALL deallocate_les_tmp

	ENDIF


  ! Get the eddy viscosity at cell centered points
  ! NU_T = C_DYN * (DELTA^2)*|S|

  !      DO J=J1,J2
	DO J=JSTART,JEND
	  DO K=ZSTART,ZEND
	    DO I=0,NXP
	      NU_T(I,K,J) = C_DYN(k,j) * DELTA_Y(k,j)**2.d0*S2X(i,k,j)
	    ENDDO
	  ENDDO
	ENDDO
	

  ! Calculate TAUij in physical space, stored in Sij

	DO ij=1,6
	  IF (LES_MODEL_TYPE.eq.2) THEN
    ! Dynamic Smagorinsky model, no scale similar part
	  IF ( ij.eq.1 ) THEN
    ! Here, Sij is defined at GYF points 
	    DO j=1,NY+1
	      DO k=1,NZ+1
		DO i=0,NXP
		  Sij(i,k,j,ij)= -2.d0 * C_DYN(k,j) * DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)
		ENDDO 
	      ENDDO
	    ENDDO
	  
	  ELSEIF (ij.EQ.2) THEN
    ! Sij(:,:,:,2) = du2/dy, but this term will be added implicitly through
    ! an eddy viscosity, so set it equal to zero here
	      DO j=1,NY+1
		DO k=1,NZ+1
		  DO i=0,NXP
		    Sij(i,k,j,ij) = 0.d0
		  ENDDO
		ENDDO
	      ENDDO
	      
	  ELSEIF (ij.EQ.3) THEN
      ! Sij(:,:,:,3) = du3/dz, but this term will be added implicitly through
      ! an eddy viscosity, so set it equal to zero here
	      DO j=0,NY+1
		DO k=1,NZ+1
		  DO i=0,NXP
		    Sij(i,k,j,ij) = 0.d0
		  ENDDO
		ENDDO
	      ENDDO
	    
	  ELSEIF (ij.EQ.4) THEN 
    ! Here, Sij is defined at GY points, interpolate C_DYN, etc 
    ! Use exact second order interpolation
    ! Sij(:,:,:,4)=0.5*(dU1/dy + dU2/dx)
    ! But dU1/dy will be accounted for as an implicit eddy viscosity term,
    ! So, subtract if off from Sij here
	  
	    DO j=JSTART,JEND
	      DO k=ZSTART,ZEND
		DO i=0,NXP
		  temp_1(i,k,j) = -2.d0 * C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)

		  Sij(i,k,j,ij)=-2.d0 * C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*(Sij(i,k,j,ij)-              &
			      0.5* ( 0.5*CJOB_12(K+1,J,2)*(U1X(I,K,J)   + U1X(I,K+1,J))     &
			      - 0.5*CJOB_12(K,J,2)*(U1X(I,K,J)   + U1X(I,K-1,J))            &
			      + 0.5*CJOB_22(K,J+1,1)*(U1X(I,K,J) + U1X(I,K,J+1))            &
			      - 0.5*CJOB_22(K,J,1)*(U1X(I,K,J)   + U1X(I,K,J-1)))/          &
			      INT_JACOB(K,J)  )
		ENDDO
	      ENDDO
	    ENDDO
	  
    
	  ELSEIF (ij.EQ.5) THEN
    ! Here, Sij is defined at GY points, interpolate C_DYN, etc 
    ! Use exact second order interpolation
    ! Sij(:,:,:,5)=0.5*(dU1/dz + dU3/dx)
    ! But dU1/dz will be accounted for as an implicit eddy viscosity term,
    ! So, subtract if off from Sij here
	    DO j=JSTART,JEND
	      DO k=ZSTART,ZEND
		DO i=0,NXP
		  temp_2(i,k,j)=-2.d0 * C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)

		  Sij(i,k,j,ij)=-2.d0 * C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*(Sij(i,k,j,ij)-0.5d0*      &
			      ( 0.5*CJOB_11(K+1,J,2)*(U1X(I,K,J)+ U1X(I,K+1,J))          &
			      - 0.5*CJOB_11(K,J,2)*(U1X(I,K,J)   + U1X(I,K-1,J))         &
			      + 0.5*CJOB_21(K,J+1,1)*(U1X(I,K,J) + U1X(I,K,J+1))         &
			      - 0.5*CJOB_21(K,J,1)*(U1X(I,K,J)   + U1X(I,K,J-1)) )/      & 
			    INT_JACOB(K,J)   )                         
		ENDDO
	      ENDDO
	    ENDDO
	    
	  ELSEIF (ij.EQ.6) THEN     
    ! Here, Sij is defined at GY points, interpolate C_DYN, etc 
    ! Use exact second order interpolation
    ! Sij(:,:,:,6)=0.5*(dU3/dy + dU2/dz)
    ! But dU3/dy as well as  dU2/dz will be accounted for as an implicit eddy viscosity term,
    ! So, subtract if off from Sij here
	    DO j=JSTART,JEND
	      DO k=ZSTART,ZEND
		DO i=0,NXP
		  temp_3(i,k,j)=-2.d0 * C_DYN(k,j) *DELTA_Y(k,j)**2.d0*S2X(i,k,j)*(Sij(i,k,j,ij)-           &
				    0.5* ( 0.5*CJOB_12(K+1,J,2)*(U3X(I,K,J)   + U3X(I,K+1,J))  &
				      - 0.5*CJOB_12(K,J,2)*(U3X(I,K,J)   + U3X(I,K-1,J))        &
				      + 0.5*CJOB_22(K,J+1,1)*(U3X(I,K,J) + U3X(I,K,J+1))        &
				      - 0.5*CJOB_22(K,J,1)*(U3X(I,K,J)   + U3X(I,K,J-1))) / INT_JACOB(K,J) )

		  Sij(i,k,j,ij)=-2.d0 * C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*(Sij(i,k,j,ij)-           &
				  0.5*( 0.5*CJOB_11(K+1,J,2)*(U2X(I,K,J)+ U2X(I,K+1,J))       &
				    - 0.5*CJOB_11(K,J,2)*(U2X(I,K,J)   + U2X(I,K-1,J))         &
				    + 0.5*CJOB_21(K,J+1,1)*(U2X(I,K,J) + U2X(I,K,J+1))         &
				    - 0.5*CJOB_21(K,J,1)*(U2X(I,K,J)   + U2X(I,K,J-1))) / INT_JACOB(K,J)  )
		ENDDO
	      ENDDO
	    ENDDO 
	  
	    DO k=0,NZ+1
	      DO i=0,NXP
		temp_3(i,k,NY+1) = temp_3(i,k,NY)
		temp_3(i,k,0)    = temp_3(i,k,1)
	      ENDDO
	    ENDDO

	    DO j=0,NY+1
	      DO i=0,NXP
	      temp_3(i,NZ+1,j) = temp_3(i,NZ,j)
	      temp_3(i,0,j)    = temp_3(i,1,j)
	      ENDDO
	    ENDDO

    ! End if ij
	  ENDIF

    !    EXTRAPOLATING VALUES 

	  DO k=0,NZ+1
	    DO i=0,NXP
	      Sij(i,k,NY+1,ij) = Sij(i,k,NY,ij)
	      Sij(i,k,0,ij) = Sij(i,k,1,ij)
	    ENDDO
	  ENDDO

	  DO j=0,NY+1
	    DO i=0,NXP
	    Sij(i,NZ+1,j,ij) = Sij(i,NZ,j,ij)
	    Sij(i,0,j,ij) = Sij(i,1,j,ij)
	    ENDDO
	  ENDDO


	ELSEIF (LES_MODEL_TYPE.eq.3) THEN
    ! Model type = 3, dynamic mixed model with scale similar part
    ! Always define temp at GYF points to match U_2BAR
    ! temp=Ui*Uj:
	    SELECT CASE (ij)
	      
	      CASE(1)
		DO j=0,NY+1
		  DO k=0,NZ+1
		    DO i=0,NXP 
		      temp(i,k,j) = U1X(i,k,j)*U1X(i,k,j)
		    ENDDO
		  ENDDO
		ENDDO
		
	      CASE(2)
		DO j=0,NY+1
		  DO k=0,NZ+1
		    DO i=0,NXP 
		      temp(i,k,j) = U2X(i,k,j)*U2X(i,k,j)
		    ENDDO
		  ENDDO
		ENDDO
	      
	      CASE(3)
		DO j=0,NY+1
		  DO k=0,NZ+1
		    DO i=0,NXP 
		      temp(i,k,j) = U3X(i,k,j)*U3X(i,k,j)
		    ENDDO
		  ENDDO
		ENDDO
	      
	      CASE(4) 
		DO j=0,NY+1
		  DO k=0,NZ+1
		    DO i=0,NXP 
		      temp(i,k,j) = U2X(i,k,j)*U1X(i,k,j)
		    ENDDO
		  ENDDO
		ENDDO
	      
	      CASE(5)
		DO j=0,NY+1
		  DO k=0,NZ+1
		    DO i=0,NXP 
		      temp(i,k,j) = U1X(i,k,j)*U3X(i,k,j)
		    ENDDO
		  ENDDO
		ENDDO
	      
	      CASE(6) 
		DO j=0,NY+1
		  DO k=0,NZ+1
		    DO i=0,NXP 
		      temp(i,k,j) = U2X(i,k,j)*U3X(i,k,j)
		    ENDDO
		  ENDDO
		ENDDO
		
	    ENDSELECT

    ! Filter temp grid filter operation filter type = 1
	    S1X = temp
	    CALL FILTER_VAR(1)
	    temp=S1X

    !       call les_filter_chan(temp,0,NY+1,1)


	    IF(LES_IMPLICIT) THEN


		IF ( ij.eq.1 ) THEN
	  ! Here, Sij is defined at GYF points 

		  DO j=0,NY+1
		    DO k=1,NZ+1
		      DO i=0,NXP
			Sij(i,k,j,ij) = temp(i,k,j) - U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij)) &
					-2.d0*C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)
		      ENDDO
		    ENDDO 
		  ENDDO
		  
		ELSEIF ( (ij.eq.2) .or. (ij.eq.3) ) THEN
	  ! Sij(:,:,:,2) = du2/dy, 
	  ! Sij(:,:,:,3) = du3/dz, but this term will be added implicitly through
	  ! an eddy viscosity, so set the Smagorinsky pert equal to zero here
		  DO j=0,NY+1
		    DO k=1,NZ+1
		      DO i=0,NXP
			Sij(i,k,j,ij) = temp(i,k,j) - U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij))
		      ENDDO
		    ENDDO
		  ENDDO
565	format(4f8.5)		  
		ELSEIF (ij.eq.4) THEN
	  ! Sij(:,:,:,4)=0.5*(dU1/dy + dU2/dx)
	  ! But the dU1/dy term will be accounted for as an implicit eddy viscosity
	  ! so subtract this term from Sij in the Smagorinsky part
		  DO j=1,NY
		    DO k=1,NZ
		      DO i=0,NXP
			temp_1(i,k,j)=  temp(i,k,j) -U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij)) &
					-2.d0*C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)

			Sij(i,k,j,ij)=  temp(i,k,j) -U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij))  &
					-2.d0*C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*(Sij(i,k,j,ij)-0.5d0*   &
						      ( 0.5*CJOB_12(K+1,J,2)*(U1X(I,K,J)   + U1X(I,K+1,J))     &
						      - 0.5*CJOB_12(K,J,2)*(U1X(I,K,J)   + U1X(I,K-1,J))        &
						      + 0.5*CJOB_22(K,J+1,1)*(U1X(I,K,J) + U1X(I,K,J+1))        &
						      - 0.5*CJOB_22(K,J,1)*(U1X(I,K,J)   + U1X(I,K,J-1)))/      &
						      INT_JACOB(K,J)   )
		      ENDDO
		    ENDDO
		  ENDDO 
		
		ELSEIF (ij.eq.5) THEN
	  ! Here, Sij is defined at GY points, interpolate C_DYN, etc 
	  ! Use exact second order interpolation
	  ! Sij(:,:,:,5)=0.5*(dU1/dz+dU3/dx) 
	  ! But the dU1/dz term will be accounted for as an implicit eddy viscosity
	  ! so subtract this term from Sij in the Smagorinsky part
		  DO j=1,NY
		    DO k=1,NZ
		      DO i=0,NXP
			temp_2(i,k,j)= temp(i,k,j)-U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij)) &
					-2.d0*C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)
	    
			Sij(i,k,j,ij)= temp(i,k,j)-U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij))  &
					-2.d0*C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*(Sij(i,k,j,ij)-0.5d0*  &
						    ( 0.5*CJOB_11(K+1,J,2)*(U1X(I,K,J)+ U1X(I,K+1,J))          &
						    - 0.5*CJOB_11(K,J,2)*(U1X(I,K,J)   + U1X(I,K-1,J))         &
						    + 0.5*CJOB_21(K,J+1,1)*(U1X(I,K,J) + U1X(I,K,J+1))         &
						    - 0.5*CJOB_21(K,J,1)*(U1X(I,K,J)   + U1X(I,K,J-1)) )/      &
						  INT_JACOB(K,J) )
		      ENDDO
		    ENDDO
		  ENDDO
		
		ELSEIF (ij.eq.6) THEN
	  ! Here, Sij is defined at GY points, interpolate C_DYN, etc 
	  ! Use exact second order interpolation
	  ! Sij(:,:,:,6)=0.5*(dU3/dy + dU2/dz)
	  ! But both  dU3/dy and dU2/dz  term will be accounted for as an implicit eddy viscosity
	  ! so subtract this term from Sij in the Smagorinsky part
		  DO j=1,NY
		    DO k=1,NZ
		      DO i=0,NXP
			temp_3(i,k,j)= temp(i,k,j) -U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij))  &
					-2.d0 * C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*(Sij(i,k,j,ij)-          &
						0.5* ( 0.5*CJOB_12(K+1,J,2)*(U3X(I,K,J)   + U3X(I,K+1,J))   &
						  - 0.5*CJOB_12(K,J,2)*(U3X(I,K,J)   + U3X(I,K-1,J))         &
						  + 0.5*CJOB_22(K,J+1,1)*(U3X(I,K,J) + U3X(I,K,J+1))         &
						  - 0.5*CJOB_22(K,J,1)*(U3X(I,K,J)   + U3X(I,K,J-1)))/       &
						  INT_JACOB(K,J) )

			Sij(i,k,j,ij)= temp(i,k,j)  -U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij))  &
					-2.d0*C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*(Sij(i,k,j,ij)-          &
						    0.5*( 0.5*CJOB_11(K+1,J,2)*(U2X(I,K,J)+ U2X(I,K+1,J))      &
						    - 0.5*CJOB_11(K,J,2)*(U2X(I,K,J)   + U2X(I,K-1,J))         &
						    + 0.5*CJOB_21(K,J+1,1)*(U2X(I,K,J) + U2X(I,K,J+1))         &
						    - 0.5*CJOB_21(K,J,1)*(U2X(I,K,J)   + U2X(I,K,J-1)))/       &
						  INT_JACOB(K,J)  ) 
			ENDDO
		      ENDDO
		    ENDDO 

	  ! End if ij
		ENDIF

	      ELSE
		DO j=JSTART,JEND
		  DO k=ZSTART,ZEND
		    DO i=0,NXP
		      Sij(i,k,j,ij) = temp(i,k,j) - U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij)) &
				      -2.d0*C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)
		    ENDDO
		  ENDDO
		ENDDO

		DO k=0,NZ+1
		  DO i=0,NXP
		    Sij(i,k,NY+1,ij) = Sij(i,k,NY,ij)
		    Sij(i,k,0,ij) = Sij(i,k,1,ij)
		  ENDDO 
		ENDDO

		DO j=0,NY+1
		  DO i=0,NXP
		    Sij(i,NZ+1,j,ij) = Sij(i,NZ,j,ij)
		    Sij(i,0,j,ij) = Sij(i,1,j,ij)
		  ENDDO
		ENDDO

	      ENDIF 

      ! End if Mixed Model
	    ENDIF

      ! Convert TAUij, now stored in Sij to Fourier space
	    S1X=Sij(:,:,:,ij)
	    
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
	    CSij(:,:,:,ij)=CS1X

      !       call FFT_X_TO_FOURIER(Sij(0,0,0,ij),CSij(0,0,0,ij),0,NY+1,0,NZ+1)
    ! End DO ij
	  ENDDO 
  ! Convert TAUij, now stored in Sij to Fourier space (for only u momentum equation)

	  S1X=temp_1
	  
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
	  Ctemp_1=CS1X

	  S1X=temp_2
	  
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
	  
	  Ctemp_2=CS1X
	  S1X=temp_3
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
	  Ctemp_3=CS1X

!       call FFT_X_TO_FOURIER(temp_1,ctemp_1,0,NY+1,0,NZ+1)
!       call FFT_X_TO_FOURIER(temp_2,ctemp_2,0,NY+1,0,NZ+1) 
!       call FFT_X_TO_FOURIER(s2,cs2,0,NY+1,0,NZ+1)
! End if LES_MODEL_TYPE dynamic Smagorinsky or Dynamic mixed model
	ELSE
	  PAUSE 'Error, unsupported LES_MODEL_TYPE is chosen'
	ENDIF

! Now, add the subgrid scale forcing to CFi
! (This includes the subgrid scale stress as an explicit R-K term

        IF (LES_IMPLICIT) THEN
         DO J=JSTART,JEND
          DO K=ZSTART,ZEND
            DO I=0,NX2P
             CF1X(I,K,J)=CF1X(I,K,J)                                              &
                     -INT_JACOB(K,J)*CIKXP(I)*CSij(I,K,J,1)                       &
                      -( 0.5*CJOB_12(K+1,J,2)*(CSij(I,K,J,4)  + CSij(I,K+1,J,4))  &
                        - 0.5*CJOB_12(K,J,2)*(CSij(I,K,J,4)   + CSij(I,K-1,J,4))  &
                        + 0.5*CJOB_22(K,J+1,1)*(CSij(I,K,J,4) + CSij(I,K,J+1,4))  &
                        - 0.5*CJOB_22(K,J,1)*(CSij(I,K,J,4)   + CSij(I,K,J-1,4)) )& 

                      -(  0.5*CJOB_11(K+1,J,2)*(CSij(I,K,J,5) + CSij(I,K+1,J,5))  &
                        - 0.5*CJOB_11(K,J,2)*(CSij(I,K,J,5)   + CSij(I,K-1,J,5))  &
                        + 0.5*CJOB_21(K,J+1,1)*(CSij(I,K,J,5) + CSij(I,K,J+1,5))  &
                        - 0.5*CJOB_21(K,J,1)*(CSij(I,K,J,5)   + CSij(I,K,J-1,5)) )             
   

             CF2X(I,K,J)=CF2X(I,K,J)                                              &
                     -INT_JACOB(K,J)*CIKXP(I)*Ctemp_1(i,k,j)                     &
                      -( 0.5*CJOB_12(K+1,J,2)*(CSij(I,K,J,2)  + CSij(I,K+1,J,2))  &
                        - 0.5*CJOB_12(K,J,2)*(CSij(I,K,J,2)   + CSij(I,K-1,J,2))   &
                        + 0.5*CJOB_22(K,J+1,1)*(CSij(I,K,J,2) + CSij(I,K,J+1,2))  &
                        - 0.5*CJOB_22(K,J,1)*(CSij(I,K,J,2)   + CSij(I,K,J-1,2)) )&

                       -(  0.5*CJOB_11(K+1,J,2)*(CSij(I,K,J,6) + CSij(I,K+1,J,6)) &
                        - 0.5*CJOB_11(K,J,2)*(CSij(I,K,J,6)   + CSij(I,K-1,J,6))  &
                        + 0.5*CJOB_21(K,J+1,1)*(CSij(I,K,J,6) + CSij(I,K,J+1,6))  &
                        - 0.5*CJOB_21(K,J,1)*(CSij(I,K,J,6)   + CSij(I,K,J-1,6)) ) 

           
           CF3X(I,K,J)=CF3X(I,K,J)                                                        &
                     -INT_JACOB(K,J)*CIKXP(I)*Ctemp_2(i,k,j)                              &

                      -( 0.5*CJOB_12(K+1,J,2)*(Ctemp_3(I,K,J)  + Ctemp_3(I,K+1,J))        &
                        - 0.5*CJOB_12(K,J,2)*(Ctemp_3(I,K,J)   + Ctemp_3(I,K-1,J))        &
                        + 0.5*CJOB_22(K,J+1,1)*(Ctemp_3(I,K,J) + Ctemp_3(I,K,J+1))        &
                        - 0.5*CJOB_22(K,J,1)*(Ctemp_3(I,K,J)   + Ctemp_3(I,K,J-1)) )      &

                      -(  0.5*CJOB_11(K+1,J,2)*(CSij(I,K,J,3) + CSij(I,K+1,J,3))          &
                        - 0.5*CJOB_11(K,J,2)*(CSij(I,K,J,3)   + CSij(I,K-1,J,3))          &
                        + 0.5*CJOB_21(K,J+1,1)*(CSij(I,K,J,3) + CSij(I,K,J+1,3))          &
                        - 0.5*CJOB_21(K,J,1)*(CSij(I,K,J,3)   + CSij(I,K,J-1,3)) ) 
                           
            ENDDO
          ENDDO
        ENDDO
       ELSE
        DO J=JSTART,JEND
          DO K=ZSTART,ZEND
            DO I=0,NX2P
             CF1X(I,K,J)=CF1X(I,K,J)                                             &
                     -INT_JACOB(K,J)*CIKXP(I)*CSij(I,K,J,1)                       &
                      -( 0.5*CJOB_12(K+1,J,2)*(CSij(I,K,J,4)  + CSij(I,K+1,J,4)) &
                        - 0.5*CJOB_12(K,J,2)*(CSij(I,K,J,4)   + CSij(I,K-1,J,4)) &
                        + 0.5*CJOB_22(K,J+1,1)*(CSij(I,K,J,4) + CSij(I,K,J+1,4)) &
                        - 0.5*CJOB_22(K,J,1)*(CSij(I,K,J,4)   + CSij(I,K,J-1,4)) ) &

                      -(  0.5*CJOB_11(K+1,J,2)*(CSij(I,K,J,5) + CSij(I,K+1,J,5)) &
                        - 0.5*CJOB_11(K,J,2)*(CSij(I,K,J,5)   + CSij(I,K-1,J,5)) &
                        + 0.5*CJOB_21(K,J+1,1)*(CSij(I,K,J,5) + CSij(I,K,J+1,5)) &
                        - 0.5*CJOB_21(K,J,1)*(CSij(I,K,J,5)   + CSij(I,K,J-1,5)) )


             CF2X(I,K,J)=CF2X(I,K,J)                                                &
                      -INT_JACOB(K,J)*CIKXP(I)*CSij(I,K,J,4)                      &
                      -( 0.5*CJOB_12(K+1,J,2)*(CSij(I,K,J,2)  + CSij(I,K+1,J,2)) &
                        - 0.5*CJOB_12(K,J,2)*(CSij(I,K,J,2)   + CSij(I,K-1,J,2)) &
                        + 0.5*CJOB_22(K,J+1,1)*(CSij(I,K,J,2) + CSij(I,K,J+1,2)) &
                        - 0.5*CJOB_22(K,J,1)*(CSij(I,K,J,2)   + CSij(I,K,J-1,2)) ) &

                       -(  0.5*CJOB_11(K+1,J,2)*(CSij(I,K,J,6) + CSij(I,K+1,J,6))&
                        - 0.5*CJOB_11(K,J,2)*(CSij(I,K,J,6)   + CSij(I,K-1,J,6)) &
                        + 0.5*CJOB_21(K,J+1,1)*(CSij(I,K,J,6) + CSij(I,K,J+1,6)) &
                        - 0.5*CJOB_21(K,J,1)*(CSij(I,K,J,6)   + CSij(I,K,J-1,6)) )
      
             CF3X(I,K,J)=CF3X(I,K,J)                                             &
                     -INT_JACOB(K,J)*CIKXP(I)*CSij(I,K,J,5)                       &
                     -( 0.5*CJOB_12(K+1,J,2)*(CSij(I,K,J,6)  + CSij(I,K+1,J,6))  &
                        - 0.5*CJOB_12(K,J,2)*(CSij(I,K,J,6)   + CSij(I,K-1,J,6)) &
                        + 0.5*CJOB_22(K,J+1,1)*(CSij(I,K,J,6) + CSij(I,K,J+1,6)) &
                        - 0.5*CJOB_22(K,J,1)*(CSij(I,K,J,6)   + CSij(I,K,J-1,6)) ) &

                      -(  0.5*CJOB_11(K+1,J,2)*(CSij(I,K,J,3) + CSij(I,K+1,J,3)) &
                        - 0.5*CJOB_11(K,J,2)*(CSij(I,K,J,3)   + CSij(I,K-1,J,3)) &
                        + 0.5*CJOB_21(K,J+1,1)*(CSij(I,K,J,3) + CSij(I,K,J+1,3)) &
                        - 0.5*CJOB_21(K,J,1)*(CSij(I,K,J,3)   + CSij(I,K,J-1,3)) )

            ENDDO
          ENDDO
        ENDDO

      ENDIF


! Periodically, output mean quantities
      IF ((MOD(TIME_STEP,SAVE_STATS_INT).EQ.0).AND.(RK_STEP.EQ.1)) THEN
! Get plane averages

       IF (LES_MODEL_TYPE .EQ. 1) C_DYN(:,:)=C_SMAG**2
       
       DO j=0,NY+1
        DO k=0,NZ+1
         NU_T_MEAN(k,j) = SUM(NU_T(0:min(NXP,NXP_L),k,j))/dble(NX)
        ENDDO
       ENDDO

       CALL MPI_COMBINE_STATS(NU_T_MEAN,NZ+2,NY+2)

!     Output SGS contribution to the TKE equation
!     LES term in the tke equation: -<u_i'dtau_ij'/dx_j>

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     sgs contribution to the transport: -dtau_ij'u_i'/dx_j 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

     
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     sgs contribution to the production: -<tau_ij><S_ij>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      

	DO ij = 1,6	!S_ij is already stored in st_rij(i,k,j,ij)
	  DO j=1,NY
	    DO k=1,NZ
	      DO i=0,NXP
		Sij(i,k,j,ij) = -2.d0*NU_T(I,K,J)*St_rij(i,k,j,ij)
	      ENDDO
	    ENDDO
	  ENDDO
	ENDDO

!     To get S_ij for strain rate tensor we need to store tau_ij 
  
!       DO ij = 1,6
!       
!         DO j=1,NY
!          DO k=1,NZ
!           DO i=0,NX2P
!             ctemp(i,k,j) = CSij(i,k,j,ij)   
!           ENDDO
!          ENDDO
!         ENDDO
! 
!       if ( (ij.eq. 1).or.(ij.eq.3).or.(ij.eq.5)) then
! 
!       CS1X=Ctemp
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
!       temp=S1X
! 
!        
! !        call fft_x_to_physical(ctemp,temp,0,NY+1,0,NZ+1)
!         DO j=1,NY
!          DO k=1,NZ
!           DO i=0,NXP
!             Sij(i,k,j,ij) = temp(i,k,j)
!           ENDDO
!          ENDDO
!         ENDDO
!       else if (ij.eq.2) then
! ! we need to take into account Smagorinsky part du2/dy
!       CS1X=Ctemp
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
!       temp=S1X
! 
! !        call fft_x_to_physical(ctemp,temp,0,NY+1,0,NZ+1)
!         DO j=1,NY
!          DO k=1,NZ
!           DO i=0,NXP
!             Sij(i,k,j,ij) = temp(i,k,j)-2.d0*C_DYN(k,j)*DELTA_Y(k,j)**2.d0 &
!                            *S2X(i,k,j)*st_rij(i,k,j,ij)
!           ENDDO
!          ENDDO
!         ENDDO
!       else if (ij.eq.4) then
! 
!       CS1X=Ctemp_1
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
!       temp_1=S1X
! 
! !        call fft_x_to_physical(ctemp_1,temp_1,0,NY+1,0,NZ+1)
!         DO j=1,NY
!          DO k=1,NZ
!           DO i=0,NXP
!             Sij(i,k,j,ij) = temp_1(i,k,j)
!           ENDDO
!          ENDDO
!         ENDDO
!       else if (ij.eq.6) then
! 
!       CS1X=Ctemp_2
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
!       temp_2=S1X
! 
! 
! !        call fft_x_to_physical(ctemp_2,temp_2,0,NY+1,0,NZ+1)
!         DO j=1,NY
!          DO k=1,NZ
!           DO i=0,NXP
!             Sij(i,k,j,ij) = temp_2(i,k,j)
!           ENDDO
!          ENDDO
!         ENDDO
!       ENDif
!        
!       ENDDO

! Need to calculate Strain Rate

      DO ij=1,6
       DO j=0,NY+1
        DO k=0,NZ+1 
         TAU_mean(k,j,ij) = SUM(Sij(0:min(NXP,NXP_L),k,j,ij))/dble(NX)
        ENDDO
       ENDDO
        CALL MPI_COMBINE_STATS(TAU_mean(0,0,ij),NZ+2,NY+2)
      ENDDO
         
!c      DO j=1,NY
!c          check_tau(j) = 0.d0
!c          check_tau(j)=SUM(temp_1(0:NXM,0:NZM,j))/dble(NX*NZ)
!c      ENDDO 



      DO j=1,NY
       DO k=1,NZ 
         tke_sgs_p(k,j) = 0.d0
         DO ij=1,3 
	    tke_sgs_p(k,j)=tke_sgs_p(k,j)-Sij_mean(k,j,ij)*TAU_mean(k,j,ij) 
         ENDDO 
	ENDDO
      ENDDO     

      
      DO j=1,NY
	DO k=1,NZ
	  tke_sgs_p(k,j)= tke_sgs_p(k,j)- 2.0*Sij_mean(k,j,4)*TAU_mean(k,j,4)  &
			  - 2.0*Sij_mean(k,j,6)*TAU_mean(k,j,6)       &
			  - 2.0*Sij_mean(k,j,5)*TAU_mean(k,j,5)  
	ENDDO
      ENDDO
        
         


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     sgs contribution to the dissipation: -<tau_ij*S_ij>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DO j=1,NY
       DO k=1,NZ
         
         tke_sgs_diss(k,j) = 0.d0
         
         DO ij=1,6
           IF (ij.LE.3) THEN
! We are computing a diagonal term
             cross=1.d0
           ELSE
! We are computing a cross term
            cross=2.d0
           ENDIF
           
           DO i=0,NXP
             temp(i,k,j) =   cross*Sij(i,k,j,ij)*St_rij(i,k,j,ij)
           ENDDO
	   
	   tke_sgs_diss(k,j)=tke_sgs_diss(k,j)+SUM(temp(0:min(NXP,NXP_L),k,j)) /dble(NX) 
	   
         ENDDO
	ENDDO
       ENDDO  
       
      CALL MPI_COMBINE_STATS(tke_sgs_diss,NZ+2,NY+2)


       IF (rank .eq. 0) THEN
       
	  !check for directory existence
	  INQUIRE(DIRECTORY='./plane_les/.', EXIST=dir_exists) 

	  IF ( .NOT. dir_exists) THEN
	    WRITE(6,*) 'plane_les directory DOes not exists'
	    CALL system('mkdir plane_les')
	    WRITE(6,*) 'plane_les is created!'
	  ENDIF 
	
	
	k = time_step/SAVE_STATS_INT

	file_name = 'plane_les/data_les_'  &
	      //CHAR(MOD(k,100000)/10000+48) &
	      //CHAR(MOD(k,10000)/1000+48)   &
	      //CHAR(MOD(k,1000)/100+48)     &
	      //CHAR(MOD(k,100)/10+48)       &
	      //CHAR(MOD(k,10)+48) //        &
	      '.pln'
	file_name_plt = 'plane_les/data_les_'  &
	      //CHAR(MOD(k,100000)/10000+48) &
	      //CHAR(MOD(k,10000)/1000+48)   &
	      //CHAR(MOD(k,1000)/100+48)     &
	      //CHAR(MOD(k,100)/10+48)       &
	      //CHAR(MOD(k,10)+48) //        &
	      '.plt'

	CALL plot_les_tecplot(file_name, file_name_plt)
       ENDIF 
      

!      deallocate (tke_sgs_diss)
!      deallocate (tke_sgs_p)
!      deallocate (tke_sgs_t)
!      deallocate (Sij_mean)
!      deallocate (TAU_mean)
!      deallocate (NU_T_mean)
!       write(6,*) 'Deallocating TKE_LES_vars'
 
      ENDIF

  RETURN
END



SUBROUTINE compute_strain_chan

!C This subroutine computes S_ij for the filtered velocity field
!C The input velocity field should be in fourier space in the periodic
!C directions.
!C For use in the LES model in channel flow (2 periodic directions)
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
  USE les_chan_var

    IMPLICIT NONE

    INTEGER I,J,K,ij

       
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NX2P
	  CSij(I,K,J,1)=CIKXP(I)*CU1X(I,K,J)

! du2/dy=> d[J-1*\zeta_y*u2]/d\zeta + d[J-1*\eta_y*u2]/d\eta

	  CSij(I,K,J,2)= ( 0.5*CJOB_12(K+1,J,2)*(CU2X(I,K,J)  + CU2X(I,K+1,J))     &
		      - 0.5*CJOB_12(K,J,2)*(CU2X(I,K,J)       + CU2X(I,K-1,J))     &
		      + 0.5*CJOB_22(K,J+1,1)*(CU2X(I,K,J)     + CU2X(I,K,J+1))     &
		      - 0.5*CJOB_22(K,J,1)*(CU2X(I,K,J)       + CU2X(I,K,J-1)) )/  &
			  INT_JACOB(K,J) 

!  du3/dz => d[J-1*\zeta_x*p]/d\zeta + d[J-1*\eta_x*p]/d\eta
	  CSij(I,K,J,3) = ( 0.5*CJOB_11(K+1,J,2)*(CU3X(I,K,J) + CU3X(I,K+1,J))     &
		      - 0.5*CJOB_11(K,J,2)*(CU3X(I,K,J)       + CU3X(I,K-1,J))     &
		      + 0.5*CJOB_21(K,J+1,1)*(CU3X(I,K,J)     + CU3X(I,K,J+1))     &
		      - 0.5*CJOB_21(K,J,1)*(CU3X(I,K,J)       + CU3X(I,K,J-1)) ) / &
			  INT_JACOB(K,J)

!  1/2(du1/dy+du2/dx) =>
	  CSij(I,K,J,4)=0.5d0*( ( 0.5*CJOB_12(K+1,J,2)*(CU1X(I,K,J) + CU1X(I,K+1,J))  &
		      - 0.5*CJOB_12(K,J,2)*(CU1X(I,K,J)   + CU1X(I,K-1,J))            &
		      + 0.5*CJOB_22(K,J+1,1)*(CU1X(I,K,J) + CU1X(I,K,J+1))            &
		      - 0.5*CJOB_22(K,J,1)*(CU1X(I,K,J)   + CU1X(I,K,J-1)) ) /        &
			    INT_JACOB(K,J)                                             &
		      + CIKXP(I)*CU2X(I,K,J) )

!  1/2(du1/dz+du3/dx) =>
	  CSij(I,K,J,5)=0.5d0*( ( 0.5*CJOB_11(K+1,J,2)*(CU1X(I,K,J) + CU1X(I,K+1,J)) &
		      - 0.5*CJOB_11(K,J,2)*(CU1X(I,K,J)   + CU1X(I,K-1,J))           &
		      + 0.5*CJOB_21(K,J+1,1)*(CU1X(I,K,J) + CU1X(I,K,J+1))           &
		      - 0.5*CJOB_21(K,J,1)*(CU1X(I,K,J)   + CU1X(I,K,J-1))  )/       &
			    INT_JACOB(K,J)                                            &
		      +CIKXP(I)*CU3X(I,K,J) )
  
!  1/2(du2/dz+du3/dy) =>
	  CSij(I,K,J,6)=0.5d0*( 0.5*CJOB_11(K+1,J,2)*(CU2X(I,K,J) + CU2X(I,K+1,J)) &
		      - 0.5*CJOB_11(K,J,2)*(CU2X(I,K,J)   + CU2X(I,K-1,J))         &
		      + 0.5*CJOB_21(K,J+1,1)*(CU2X(I,K,J) + CU2X(I,K,J+1))         &
		      - 0.5*CJOB_21(K,J,1)*(CU2X(I,K,J)   + CU2X(I,K,J-1))         &
			      + 0.5*CJOB_12(K+1,J,2)*(CU3X(I,K,J) + CU3X(I,K+1,J)) &
		      - 0.5*CJOB_12(K,J,2)*(CU3X(I,K,J)   + CU3X(I,K-1,J))         &
		      + 0.5*CJOB_22(K,J+1,1)*(CU3X(I,K,J) + CU3X(I,K,J+1))         &
		      - 0.5*CJOB_22(K,J,1)*(CU3X(I,K,J)   + CU3X(I,K,J-1)) )/      &
			    INT_JACOB(K,J)

	ENDDO
      ENDDO
    ENDDO


! Convert rate of strain tensor to physical space
    DO ij=1,6
      CS1X=CSij(:,:,:,ij)
      
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
      Sij(:,:,:,ij)=S1X

!        call FFT_X_TO_PHYSICAL(CSij(0,0,0,ij),Sij(0,0,0,ij),0,NY+1,0,NZ+1)
    ENDDO 
! We now have S_ij in Physical space
 RETURN
END

SUBROUTINE les_filter_chan(A,zstart,zEND,jstart,jEND, filter_type)
! This subroutine applies the les filter to the input field
! The indices to the start and END of the array in the y-direction
! The array that is passed should be in physical space
  USE ntypes
  USE Domain

    IMPLICIT NONE

    INTEGER i,k,j,zstart,zEND,jstart,jEND,filter_type


    REAL*8 A(0:NXV-1,0:NZP,0:NY+1)
    REAL*8 B(0:NXV-1,0:NZP,0:NY+1)

    INTEGER im2(0:NX-1),im1(0:NX-1),ip1(0:NX+1),ip2(0:NX+2)

! These are the weights for the filtering operation USEd
    REAL*8 W0,W1,W2,Wm1,Wm2

! filter type = 1: grid filter operation
! filter type = 2: test filter operation
! The following is for the 3-point trapezoidal rule, alpha*beta=sqrt(6)
      
    IF ( filter_type .EQ. 1 ) THEN
      Wm2=0.d0
      Wm1=1.d0/8.d0
      W0=6.d0/8.d0
      W1=1.d0/8.d0
      W2=0.d0
    ELSEIF ( filter_type .EQ. 2 ) THEN
!     else 
!      Wm1_j=1.d0/4.d0  
!      W0_j=1.d0/2.d0
!      W1_j=1.d0/4.d0
! The following is for the 5-point trapezoidal rule, alpha*beta=9
!       Wm2=1.d0/8.d0
!       Wm1=1.d0/4.d0
!       W0=1.d0/4.d0
!       W1=1.d0/4.d0
!       W2=1.d0/8.d0

       Wm2=0.d0
       Wm1=1.d0/4.d0
       W0=1.d0/2.d0
       W1=1.d0/4.d0
       W2=0.d0 
      
    ELSE
      PAUSE 'Error, unsupported LES_FILTER_TYPE chosen'
    ENDIF

!      NXM=NX-1
!      NZM=NZ-1

!      DO j=0,NY+1
!        DO k=0,NZM
!          DO i=0,NXM
!            B(i,k,j)=A(i,k,j)
!          ENDDO
!        ENDDO
!      ENDDO

! Filter in the periodic directions using cshift
! Note, cshift if not used since it appears to be slower
! Apply filter in the x-direction
!      B=Wm2*CSHIFT(B,-2,1)+Wm1*CSHIFT(B,-1,1)+W0*B+W1*CSHIFT(B,1,1)
!     &       +W2*CSHIFT(B,2,1)

! Filter using more efficient F77 syntax:
! Set up array to loop around periodic directions
    DO i=2,NXM
      im2(i)=i-2
    ENDDO
    
    im2(1)=NXM
    im2(0)=NX-2
    
    DO i=1,NXM
      im1(i)=i-1
    ENDDO
    
    im1(0)=NXM
    
    DO i=0,NX-2
      ip1(i)=i+1
    ENDDO
    
    ip1(NXM)=0
    
    DO i=0,NX-3
      ip2(i)=i+2    
    ENDDO
    
    ip2(NX-2)=0
    ip2(NXM)=1

    DO j=jstart,jEND
      DO k=0,NZP
	DO i=0,NXM
	  B(i,k,j)=Wm2*A(im2(i),k,j)+Wm1*A(im1(i),k,j)+W0*A(i,k,j) &
	    +W1*A(ip1(i),k,j)+W2*A(ip2(i),k,j)
	ENDDO
      ENDDO  
    ENDDO

    DO j=jstart,JEND
      DO k=0,NZP
	DO i=0,NXM
	  A(i,k,j) = B(i,k,j)
	ENDDO
      ENDDO
    ENDDO

 RETURN
END


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C *** Subroutine to write a .plt file (uses TecPlot tecio.a library)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE plot_les_tecplot(file_name, file_name_plt)

  USE ntypes
  USE Domain
  USE Grid
  USE run_variable
  USE les_chan_var

    IMPLICIT NONE

    INTEGER  k,ncycles
    INTEGER  j,imax,jmax,kmax
    INTEGER  debug,ier,itot
    INTEGER  tecini,tecdat,teczne,tecEND
    INTEGER  visDOuble,disDOuble
    REAL*8   phase,PI
    CHARACTER*1 nulchar
    CHARACTER(len=33) :: title
    CHARACTER*28 file_name,file_name_plt

    PI  = ATAN(1.0)*4.0
    nulchar = char(0)
    debug   = 0
    visDOuble = 0
    disDOuble = 1
    imax = NZ+2
    jmax = NY+2
    kmax = 1

    WRITE(6,*) file_name_plt 

    OPEN(22,file=file_name,status='unknown',form='unformatted')
    WRITE(22)xpoint,ypoint,U3_bar_les(0:NZ+1,0:NY+1),U2_bar_les(0:NZ+1,0:NY+1),&
	    U1_bar_les(0:NZ+1,0:NY+1),C_DYN(0:NZ+1,0:NY+1),NU_T_mean(0:NZ+1,0:NY+1),&
	    tke_sgs_p(0:NZ+1,0:NY+1),tke_sgs_diss(0:NZ+1,0:NY+1)
    CLOSE(22)  

    ncycles = time/(2.d0*pi)
    phase = (time-ncycles*2.d0*pi)*180.d0/pi
      
    WRITE(title,'(a5,f8.4,a7,f5.1)') 'time=',time,' phase=', phase
      
     ! Write the zone header information.
!      ier = tecini(trim(title)//nulchar,'x,z,ume,wme,vme,C_dyn,nu_t_mean,& 
!              &tke_sgs_p,tke_sgs_diss'&
     ier = tecini(trim(title)//nulchar,'x,z,C_dyn,nu_t_mean,tke_sgs_p,tke_sgs_diss'&
                  &//nulchar,&
                 &file_name_plt//nulchar,&
                 &'.'//nulchar,&
                 &debug,visDOuble)

     ier = teczne(trim(title)//nulchar,  &
                 imax,jmax,kmax,         &
                 'BLOCK'//nulchar,nulchar)

      
     ! Write out the field data.
     itot = imax*jmax*kmax
     ier = tecdat(itot,xpoint,disDOuble)
     ier = tecdat(itot,ypoint,disDOuble)
!      ier = tecdat(itot,U3_bar_les(0:NZ+1,0:NY+1),disDOuble)
!      ier = tecdat(itot,U2_bar_les(0:NZ+1,0:NY+1),disDOuble)
!      ier = tecdat(itot,U1_bar_les(0:NZ+1,0:NY+1),disDOuble)
     ier = tecdat(itot,SQRT(C_DYN(0:NZ+1,0:NY+1)),disDOuble)
     ier = tecdat(itot,NU_T_mean(0:NZ+1,0:NY+1),disDOuble)
     ier = tecdat(itot,tke_sgs_p(0:NZ+1,0:NY+1),disDOuble)
     ier = tecdat(itot,tke_sgs_diss(0:NZ+1,0:NY+1),disDOuble)

     ! close file
     ier = tecEND()
 
  RETURN
END

SUBROUTINE LES_WALL_MODEL
! This is wall model sunroutine which applies Monin-Obokhov smilarity solution on the first grid point above ground
! for more info please visit: http://glossary.ametsoc.org/wiki/Monin-obukhov_similarity_theory
! we (Iman and masoud) decided to get alinear relation between 3th cell and 1st cell as viscous sub-layer
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
  USE les_chan_var
  USE mpi_var, only : rank
  
    INTEGER I, K, J, height_index, npoints
    REAL*8 KAPA, Z_0, U_STAR, V_STAR, TAO_WALL, DELTA_Z, OBOKHOV_BETA, Z, OBOKHOV_L, U_1, W_1, NU_tmp, M
    REAL*8 avg_buoyant_flux, avg_th, avg_moment_flux, avg_moment_flux_1, TAO_WALL_1
    REAL*8 tmp_th(0:NZ),tmp_th_1(0:NZ),tmp_u(0:NZ),tmp_w(0:NZ)
    
    KAPA=0.41d0				! von karman constant
    OBOKHOV_BETA=4.7d0			! Monin-Obokhov constant from Businger et al (1971)
    Z_0= 0.1d0  			! surface roughness
    DELTA_Z = DY(3) + DY(2)		! what is the distance between 3th cell height and ground
    Z=DYY(2,2)               		! the height where in velocity will be updated
    height_index = FLOOR (0.09* NY)	! what is the vertical index you want to take vertical average
    alpha=20.6				! Spalart surface stress angle (1989)
    
    avg_buoyant_flux = 0.d0;		avg_moment_flux=0.d0;
    avg_moment_flux_1 = 0.d0;		avg_th = 0.d0;
    TAO_WALL = 0.d0;			TAO_WALL_1 = 0.d0;
    NU_tmp = NU			! define the pre-desribed NU (just in wall model function) or, set NU_tmp=NU
    OBOKHOV_L= 1.0E+8		! initial value for MO length scale

    IF (WALL_MODEL_TYPE.EQ.2) THEN
    ! This is Monin Obokhov similarity theory
    ! Computing horizontal average of buoyancy flux for Monin-Obokhov Lenght scale (the averaged has been computed over 5 vertical cells to make sure it's near surface)
       DO J=JSTART_TH(1),height_index
	DO K=ZSTART_TH(1),ZEND_TH(1)
	  DO I=0,MIN(NXP,NXP_L)
	   ! avg_buoyant_flux = avg_buoyant_flux + ( THX(i,k,j,1)-dble(CTHX_WALL_MODEL(0,k,j,1)) ) * ( U2X(i,k,j)-dble(CU2X_WALL_MODEL(0,k,j)) )
	    avg_moment_flux = avg_moment_flux+ ( U3X(i,k,j)-dble(CU3X_WALL_MODEL(0,k,j)) ) * ( U2X(i,k,j)-dble(CU2X_WALL_MODEL(0,k,j)) )
	    avg_moment_flux_1 = avg_moment_flux_1+ ( U1X(i,k,j)-dble(CU1X_WALL_MODEL(0,k,j)) ) * ( U2X(i,k,j)-dble(CU2X_WALL_MODEL(0,k,j)) )
	    avg_th = avg_th + THX(i,k,j,1)
	    !TAO_WALL = TAO_WALL + NU_tmp * ABS(U3X(i,k,JSTART)-U3X(i,k,JSTART-1))
            ENDDO
	ENDDO
      ENDDO
      
      npoints = (height_index-JSTART_TH(1)+1) * (ZEND_TH(1)-ZSTART_TH(1)+1) * MIN(NXP,NXP_L)

!      avg_buoyant_flux = avg_buoyant_flux / float (npoints)
      avg_moment_flux = avg_moment_flux / float (npoints)
      avg_moment_flux_1 = avg_moment_flux_1 / float (npoints)
      avg_th = avg_th / float (npoints)
      
      DO K=ZSTART,ZEND
       tmp_th(K) = SUM(THX(0:MIN(NXP,NXP_L),K,JSTART-1,1))/(MIN(NXP,NXP_L))+TH_BAR(K,0)
       tmp_th_1(K) =SUM(THX(0:MIN(NXP,NXP_L),K,JSTART,1))/(MIN(NXP,NXP_L))+TH_BAR(K,1)
!       WRITE (*,*) tmp_th(K), tmp_th_1(K)
      ENDDO

      U_STAR = 0.28d0 !SQRT(ABS(avg_moment_flux))
      V_STAR = SQRT(ABS(avg_moment_flux_1))      
      TAO_WALL = U_STAR**2.d0 
      avg_buoyant_flux = U_STAR * KAPA *(SUM(tmp_th(1:NZ))/NZ-SUM(tmp_th_1(1:NZ)/NZ))/(LOG(Z/Z_0) + 7.8d0 *Z/OBOKHOV_L )      
      
      CALL MPI_AVG (avg_buoyant_flux)
      CALL MPI_AVG (avg_moment_flux)
      CALL MPI_AVG (avg_th)
      CALL MPI_AVG (U_STAR)
      CALL MPI_AVG (V_STAR) 

      IF (RI_FINAL .NE. 0.d0) THEN 
        OBOKHOV_L =  U_STAR**3.0d0 / (KAPA * RI_FINAL * avg_buoyant_flux)
      ELSE 
        OBOKHOV_L= 1.0E+8
      ENDIF
      
      IF (OBOKHOV_L .LE. (5.d0*Z)) OBOKHOV_L=5.d0*Z

      IF ( (RK_STEP .EQ. 3) .AND. (rank .EQ. 0) .AND. (MOD(TIME_STEP,SAVE_STATS_INT).EQ.0) ) THEN
        
        OPEN(676, file='plane_data/L_OBOKHOV_USTAR.dat', form='formatted', status='unknown', position='append' )
        WRITE(676,5675) TIME, OBOKHOV_L, U_STAR,avg_buoyant_flux, V_STAR
        CLOSE(676)
5675     format(5f22.7)
      ENDIF

    ENDIF

    IF (WALL_MODEL_TYPE.EQ.3) THEN
! The Wall model proposed by John tylor in his IJHH 2008
! john taylor model from his phd thesis, pp 171	
      DO K=ZSTART,ZEND
       tmp_w(K) = SUM(U3X(0:MIN(NXP,NXP_L),K,JSTART))/(MIN(NXP,NXP_L))
       tmp_u(K) = SUM(U1X(0:MIN(NXP,NXP_L),K,JSTART))/(MIN(NXP,NXP_L))
      ENDDO
      
      U_1 = SUM(tmp_u(1:NZ))/NZ; W_1 = SUM(tmp_w(1:NZ))/NZ; 
      M = SQRT(U_1**2.d0+W_1**2.d0)
      
      ! Computing U_STAR
       U_STAR = M
       DO I=1,30
         U_STAR=M/(1/KAPA*LOG(Z/NU_tmp*U_STAR)+5.2)
       ENDDO
     
!      WRITE(*,*) U_STAR,(U_STAR**2*COS(alpha*3.1415/180)/NU_tmp)/(0.5*(CJOB_22(K,1,2)+CJOB_22(K+1,1,1)))*INT_JACOB(K,1),M,SUM(U3X(0:MIN(NXP,NXP_L),1:NZ,JSTART-1))/(NZ*MIN(NXP,NXP_L))
 
      CALL MPI_AVG (U_STAR)
      
      TAO_WALL = U_STAR**2*COS(alpha*3.1415/180)
      TAO_WALL_1 = U_STAR**2*SIN(alpha*3.1415/180)

      IF ( (RK_STEP .EQ. 3) .AND. (rank .EQ. 0) .AND. (MOD(TIME_STEP,SAVE_STATS_INT).EQ.0)) THEN
        OPEN(676, file='plane_data/L_OBOKHOV_USTAR.dat', form='formatted', status='unknown', position='append' )
        WRITE(676,567) TIME,U_STAR
        CLOSE(676)
567     format(2f22.7)
      ENDIF

    ENDIF
    
    
    IF (WALL_MODEL_TYPE.EQ.4) THEN
! This is Monin Obokhov similarity theory, second version

      DO K=ZSTART,ZEND
       tmp_th(K) = SUM(THX(0:MIN(NXP,NXP_L),K,JSTART-1,1))/(MIN(NXP,NXP_L))+TH_BAR(K,JSTART-1)
       tmp_th_1(K) =SUM(THX(0:MIN(NXP,NXP_L),K,JSTART,1))/(MIN(NXP,NXP_L))+TH_BAR(K,JSTART)

       tmp_w(K) =SUM(SQRT(U3X(0:MIN(NXP,NXP_L),K,JSTART-1)**2.0d0+U1X(0:MIN(NXP,NXP_L),K,JSTART-1)**2.d0))/(MIN(NXP,NXP_L))
!       tmp_u(K) = SUM(U1X(0:MIN(NXP,NXP_L),K,JSTART-1))/(MIN(NXP,NXP_L))
      ENDDO

!     U_STAR = SQRT( (SUM(tmp_u(1:NZ))/NZ)**2.d0+(SUM(tmp_w(1:NZ))/NZ)**2.d0) * KAPA / ( LOG(Z/Z_0) + OBOKHOV_BETA * Z/OBOKHOV_L )
!     U_STAR = (SUM(tmp_w(1:NZ))/NZ) * KAPA / ( LOG(Z/Z_0) + OBOKHOV_BETA * Z/OBOKHOV_L )
      U_STAR = 0.28d0
      avg_buoyant_flux = U_STAR * KAPA * (SUM(tmp_th(1:NZ))/NZ-SUM(tmp_th_1(1:NZ))/NZ ) /(LOG(Z/Z_0) + 7.8d0 *Z/OBOKHOV_L ) 
!     WRITE(*,*) avg_buoyant_flux, U_STAR,SUM(tmp_u)/NZ, SUM(tmp_w)/NZ

      CALL MPI_AVG (avg_buoyant_flux)
      CALL MPI_AVG (U_STAR)

      IF (RI_FINAL .NE. 0.d0) THEN
        OBOKHOV_L =  U_STAR**3.0d0 / (KAPA * RI_FINAL * avg_buoyant_flux)
      ELSE
        OBOKHOV_L= 1.0E+8
      ENDIF
      
      IF (OBOKHOV_L .LE. (5.d0*Z)) OBOKHOV_L=5.d0*Z

      IF ( (RK_STEP .EQ. 3) .AND. (rank .EQ. 0) .AND. (MOD(TIME_STEP,SAVE_STATS_INT).EQ.0)) THEN
        OPEN(676, file='plane_data/L_OBOKHOV_USTAR.dat', form='formatted', status='unknown', position='append' )
        WRITE(676,5676) TIME, OBOKHOV_L, U_STAR,avg_buoyant_flux
        CLOSE(676)
5676     format(4f22.7)
      ENDIF

    ENDIF

    
    DO K=ZSTART,ZEND
      DO I=0,MIN(NXP,NXP_L)
      !updating the height where in velocity will be updated
	Z = DYY(K,2)
        
!	IF (WALL_MODEL_TYPE.EQ.2) THEN 
	  !avg_moment_flux = ( U3X(i,k,JSTART)-dble(CU3X_WALL_MODEL(0,k,JSTART)) ) * ( U2X(i,k,JSTART)-dble(CU2X_WALL_MODEL(0,k,JSTART)) )
	  !avg_buoyant_flux = ( THX(i,k,JSTART_TH(1),1)-dble(CTHX_WALL_MODEL(0,k,JSTART_TH(1),1)) ) * ( U2X(i,k,JSTART)-dble(CU2X_WALL_MODEL(0,k,JSTART)) )
	  !U_STAR = sqrt(ABS(avg_moment_flux))   
	  !TAO_WALL = U_STAR**2.d0
	  
!	ELSEIF ( (WALL_MODEL_TYPE.EQ.1) .OR. (WALL_MODEL_TYPE.EQ.3) ) THEN
!	  TAO_WALL = NU_tmp * ABS(U3X(i,k,2)-U3X(i,k,1))
!	  U_STAR = sqrt(TAO_WALL)
!	ENDIF
        
        IF ((WALL_MODEL_TYPE.EQ.2) .AND. (TIME_STEP.GT.1)) THEN 
	  IF (RI_FINAL .NE. 0.d0) THEN
            OBOKHOV_L =  U_STAR**3.0d0 / (KAPA * RI_FINAL * avg_buoyant_flux)
          ELSE
            OBOKHOV_L= 1.0E+8
          ENDIF

	  IF (OBOKHOV_L .LE. (5.d0*Z)) OBOKHOV_L=5.d0*Z

	  W_BC_LOWER_WALLMODEL(I,K) = U_STAR / KAPA * ( LOG(Z/Z_0 + 1.0d0) + OBOKHOV_BETA * Z/OBOKHOV_L )
          
          IF (U_BC_YMIN.EQ.3) THEN
           U_BC_LOWER_WALLMODEL(I,K) = U1X(I,K,JSTART) * U3X(I,K,JSTART-1)/U3X(I,K,JSTART)
          ENDIF 
        ELSEIF ((WALL_MODEL_TYPE.EQ.1) .AND. (TIME_STEP.GT.1)) THEN
	  W_BC_LOWER_WALLMODEL(I,K) = U_STAR / KAPA * ( LOG(Z/Z_0 + 1.0d0) )
        
        ELSEIF ((WALL_MODEL_TYPE.EQ.3) .AND. (TIME_STEP.GT.1)) THEN
          
	  W_BC_LOWER_WALLMODEL(I,K) = U3X(i,k,JSTART)/W_1*(TAO_WALL/NU_tmp)/(0.5*(CJOB_22(K,1,2)+CJOB_22(K+1,1,1))) * INT_JACOB(K,1)
	  U_BC_LOWER_WALLMODEL(I,K) = MAX(U1X(i,k,JSTART)/W_1*(TAO_WALL/NU_tmp), U1X(i,k,JSTART)/U_1*(TAO_WALL_1/NU_tmp)) / (0.5*(CJOB_22(K,1,2)+CJOB_22(K+1,1,1))) * INT_JACOB(K,1)
	  
	ELSEIF ((WALL_MODEL_TYPE.EQ.4) .AND. (TIME_STEP.GT.1)) THEN 
!	  M = SQRT(U3X(I,K,JSTART-1)**2+U1X(I,K,JSTART-1)**2) 
!	  U_STAR = M * KAPA / ( LOG(Z/Z_0) + OBOKHOV_BETA * Z/OBOKHOV_L )
!         M = SUM(tmp_w(1:NZ))/NZ
!	  W_BC_LOWER_WALLMODEL(I,K) = U3X(I,K,JSTART) - U_STAR**2.d0 * ( U3X(I,K,JSTART-1)/M ) / (0.5*(CJOB_22(K,1,2)+CJOB_22(K+1,1,1))) * INT_JACOB(K,1)
!         W_BC_LOWER_WALLMODEL(I,K) = U_STAR**2.d0 * (U3X(I,K,JSTART-1)/M ) / (0.5*(CJOB_22(K,1,2)+CJOB_22(K+1,1,1))) * INT_JACOB(K,1)          

          W_BC_LOWER_WALLMODEL(I,K) = U_STAR / KAPA * ( 1.d0/Z + OBOKHOV_BETA/OBOKHOV_L )/ (0.5*(CJOB_22(K,1,2)+CJOB_22(K+1,1,1))) * INT_JACOB(K,1)

          IF (U_BC_YMIN.EQ.3) THEN 
!          U_BC_LOWER_WALLMODEL(I,K) = U1X(I,K,JSTART) - U_STAR**2.d0 * ( U1X(I,K,JSTART-1)/M ) / (0.5*(CJOB_22(K,1,2)+CJOB_22(K+1,1,1))) * INT_JACOB(K,1)
!          U_BC_LOWER_WALLMODEL(I,K) = U_STAR**2.d0 * (U1X(I,K,JSTART-1)/M ) / (0.5*(CJOB_22(K,1,2)+CJOB_22(K+1,1,1))) * INT_JACOB(K,1)
           U_BC_LOWER_WALLMODEL(I,K) = U1X(I,K,JSTART) * U3X(I,K,JSTART-1)/U3X(I,K,JSTART)
          ENDIF

        ELSE 
           IF (TIME_STEP.EQ.1) THEN
	    W_BC_LOWER_WALLMODEL(I,K) = U3X(i,k,2) 
           ELSE
	    PAUSE 'Error, unsupported WALL_MODEL_TYPE chosen'
           ENDIF
        ENDIF
        
!  Filter the stress condition at surface to make sure that, the wall model stress doesn't pass the terminal values of stress at surface which is no-slip one. 
!        IF ( W_BC_LOWER_WALLMODEL(I,K) .GT. U3X(i,k,JSTART) )  W_BC_LOWER_WALLMODEL(I,K) = U3X(i,k,JSTART)
!        IF ( ABS(U_BC_LOWER_WALLMODEL(I,K)) .GT. ABS(U1X(i,k,JSTART)) )  U_BC_LOWER_WALLMODEL(I,K) = U1X(i,k,JSTART)

      ENDDO
    ENDDO
  
  RETURN
END

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
 SUBROUTINE CMP_CS
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This is a subroutine which extracts Coherent Structures in flow and save them
! This subrotine has been called when variables are in fourier space 
  USE ntypes
  USE Domain
  USE Grid
  USE TIME_STEP_VAR
  USE FFT_var
  USE run_variable
  USE les_chan_var
  USE variable_stat  
  USE mpi_var, only : rank

    IMPLICIT NONE

    INTEGER I,J,K,ij, JEND_tmp
    LOGICAL dir_exists
        
    REAL*8  Oij(0:NXP,0:NZV-1,0:NY+1,1:3), EIG(0:NXP,0:NZV-1,0:NY+1),QQ(0:NXP,0:NZV-1,0:NY+1),P_prime(0:NXP,0:NZV-1,0:NY+1)
    COMPLEX*16 COij(0:NX2P,0:NZV-1,0:NY+1,1:3)

    REAL*8 S11,S12,S13,S22,S23,S33,O12,O13,O23,A11,A12,A13,A22,A23,A33,B,C,D,Q,R,THETA,EIG1,EIG2,EIG3
    CHARACTER*35 file_name
    CHARACTER*3 PID

! Computing the Strain rate tenstor
     CALL compute_strain_chan

! Converting pressure back to physical space (it's ignored in flow_statistics.F90)
      CALL REAL_FOURIER_TRANS_P (.false.)

!Computing the rotation tensor

JEND_tmp=80

      DO J=JSTART,JEND_tmp
	DO K=ZSTART,ZEND
	  DO I=0,NX2P
  !  1/2(du2/dz-du3/dy) =>
	    COij(I,K,J,1)=0.5d0*( (0.5*CJOB_11(K+1,J,2)*(CU2X(I,K,J) + CU2X(I,K+1,J)) &
			- 0.5*CJOB_11(K,J,2)*(CU2X(I,K,J)   + CU2X(I,K-1,J))         &
			+ 0.5*CJOB_21(K,J+1,1)*(CU2X(I,K,J) + CU2X(I,K,J+1))         &
			- 0.5*CJOB_21(K,J,1)*(CU2X(I,K,J)   + CU2X(I,K,J-1)) )       &
			- (0.5*CJOB_12(K+1,J,2)*(CU3X(I,K,J) + CU3X(I,K+1,J))         &
			- 0.5*CJOB_12(K,J,2)*(CU3X(I,K,J)   + CU3X(I,K-1,J))         &
			+ 0.5*CJOB_22(K,J+1,1)*(CU3X(I,K,J) + CU3X(I,K,J+1))         &
			- 0.5*CJOB_22(K,J,1)*(CU3X(I,K,J)   + CU3X(I,K,J-1)) ) )/      &
			      INT_JACOB(K,J)

  !  1/2(du1/dz-du3/dx) =>
	    COij(I,K,J,2)=0.5d0*( ( 0.5*CJOB_11(K+1,J,2)*(CU1X(I,K,J) + CU1X(I,K+1,J)) &
			- 0.5*CJOB_11(K,J,2)*(CU1X(I,K,J)   + CU1X(I,K-1,J))           &
			+ 0.5*CJOB_21(K,J+1,1)*(CU1X(I,K,J) + CU1X(I,K,J+1))           &
			- 0.5*CJOB_21(K,J,1)*(CU1X(I,K,J)   + CU1X(I,K,J-1))  )/       &
			      INT_JACOB(K,J)                                            &
			-CIKXP(I)*CU3X(I,K,J) )

  !  1/2(du1/dy-du2/dx) =>
	    COij(I,K,J,3)=0.5d0*( ( 0.5*CJOB_12(K+1,J,2)*(CU1X(I,K,J) + CU1X(I,K+1,J))  &
			- 0.5*CJOB_12(K,J,2)*(CU1X(I,K,J)   + CU1X(I,K-1,J))            &
			+ 0.5*CJOB_22(K,J+1,1)*(CU1X(I,K,J) + CU1X(I,K,J+1))            &
			- 0.5*CJOB_22(K,J,1)*(CU1X(I,K,J)   + CU1X(I,K,J-1)) ) /        &
			      INT_JACOB(K,J)                                             &
			- CIKXP(I)*CU2X(I,K,J) )

	  ENDDO
	ENDDO
      ENDDO
 
  ! Convert rotation tensor to physical space
      DO ij=1,3
	CS1X=COij(:,:,:,ij)

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
	Oij(:,:,:,ij)=S1X
      ENDDO

! updateing the terminal points
    Oij(:,0,:,:)=Oij(:,1,:,:)
    Oij(:,NZ+1,:,:)=Oij(:,NZ,:,:)

    Oij(:,:,0,:)=Oij(:,:,1,:)
    Oij(:,:,NY+1,:)=Oij(:,:,NY,:)

! computing CS in flow based on Lamda2 critria 
! for more info please read: "On the identidfication of a vortex" by Jinhee Jeong JFM 1995
      DO J=JSTART,JEND_tmp
	DO K=ZSTART,ZEND
	  DO I=0,MIN(NXP,NXP_L)
	    
            S11 = Sij(I,K,J,1); S12 = Sij(I,K,J,4); S13 = Sij(I,K,J,5)
            S22 = Sij(I,K,J,2); S23 = Sij(I,K,J,6); S33 = Sij(I,K,J,3)
            O12 = Oij(I,K,J,3); O13 = Oij(I,K,J,2); O23 = Oij(I,K,J,1)

	    !     S AND OMEGA ARE COMPUTED, NOW FIND S*S + OMEGA*OMEGA
	    A11 = S11*S11 + S12*S12 + S13*S13 - O12*O12 - O13*O13
	    A22 = S12*S12 + S22*S22 + S23*S23 - O12*O12 - O23*O23
	    A33 = S13*S13 + S23*S23 + S33*S33 - O13*O13 - O23*O23
	    A12 = S11*S12 + S12*S22 + S13*S23 - O13*O23
	    A13 = S11*S13 + S12*S23 + S13*S33 + O12*O23
	    A23 = S12*S13 + S22*S23 + S23*S33 - O12*O13

	    PI = 4D0*DATAN(1D0)

	    !      DETERMINE MIDDLE EIGENVALUE
	    !      COEFFICIENTS OF POLYNOMIAL X^3+B*X^2+C*X+D=0
	    B = -(A11+A22+A33);
	    C = -(A12*A12+A13*A13+A23*A23-A11*A22-A11*A33-A22*A33);
	    D = -(2D0*A12*A13*A23-A11*A23*A23-A22*A13*A13-A33*A12*A12+A11*A22*A33)

	    !      INTERMEDIATE COEFFICIENTS
	    Q = (3D0*C - B*B)/9D0; R = (9D0*C*B - 27D0*D - 2D0*B*B*B)/54D0; THETA = DACOS(R/DSQRT(-Q*Q*Q))

	    !      EIGENVALUES      
	    EIG1 = 2D0*DSQRT(-Q)*DCOS(THETA/3D0)-B/3D0; EIG2 = 2D0*DSQRT(-Q)*DCOS((THETA+2D0*PI)/3D0)-B/3D0; EIG3 = 2D0*DSQRT(-Q)*DCOS((THETA+4D0*PI)/3D0)-B/3D0
            
            ! COMPUTING Seconf invariant of velocity gradient tensor as Q=-0.5*(trace(S**2+Omega**))
            QQ(I,K,J) = -0.5 * ( EIG1 + EIG2 + EIG3) 

	    !      IF U2 IS THE MIDDLE EIGENVALUE
	    IF  ((EIG2.LT.EIG3).AND.(EIG2.GT.EIG1)) THEN
	     EIG1 = EIG2
	    ELSE IF ( (EIG2.LT.EIG1).AND.(EIG2.GT.EIG3) ) THEN
	     EIG1 = EIG2
	    !      IF U3 IS THE MIDDLE EIGENVALUE
	    ELSE IF ( (EIG3.LT.EIG1).AND.(EIG3.GT.EIG2) ) THEN
	     EIG1 = EIG3
	    ELSE IF ( (EIG3.LT.EIG2).AND.(EIG3.GT.EIG1) ) THEN
	     EIG1 = EIG3
	    ENDIF
		
	    !     STORE LAMBDA2 VALUES IN F1
	    EIG(I,K,J) = EIG1

            ! COMPUTING PRESSURE PRETURBATION
            P_PRIME(I,K,J) = PX(I,K,J)-dble(P_MEAN(K,J))
	  ENDDO
	ENDDO
      ENDDO


! saving the CS on different files
      !check for directory existence
      INQUIRE(DIRECTORY='./Chr_Strct/.', EXIST=dir_exists) 

      IF ( .NOT. dir_exists) THEN
	WRITE(6,*) 'Chr_Strct directory DOes not exists'
	CALL system('mkdir Chr_Strct')
	WRITE(6,*) 'Chr_Strct is created!'
      ENDIF 
    
      k = time_step/SAVE_STATS_INT

      PID = CHAR(MOD(Rank+1,1000)/100+48)     &
            //CHAR(MOD(Rank+1,100)/10+48)    &
            //CHAR(MOD(Rank+1,10)+48)

      file_name = 'Chr_Strct/CS_'  &
          //CHAR(MOD(k,100000)/10000+48) &
          //CHAR(MOD(k,10000)/1000+48)   &
          //CHAR(MOD(k,1000)/100+48)     &
          //CHAR(MOD(k,100)/10+48)       &
          //CHAR(MOD(k,10)+48) //        &
          '.pln_'                        &
          //PID


      IF (RANK.EQ.0) WRITE(6,*) 'start writing in pln format ',file_name

      OPEN(22,file=file_name,form='unformatted',status='unknown')
      WRITE(22) xpoint(0:NZ+1,0:JEND_tmp+1),ypoint(0:NZ+1,0:JEND_tmp+1),dx(1),EIG(0:MIN(NXP,NXP_L),0:NZ+1,:)
      WRITE(22) QQ(0:MIN(NXP,NXP_L),0:NZ+1,:)
      WRITE(22) Oij(0:MIN(NXP,NXP_L),0:NZ+1,:,1)
      WRITE(22) Oij(0:MIN(NXP,NXP_L),0:NZ+1,:,2)
      WRITE(22) Oij(0:MIN(NXP,NXP_L),0:NZ+1,:,3)
      WRITE(22) P_PRIME(0:MIN(NXP,NXP_L),0:NZ+1,:)
      CLOSE(22)
  RETURN
END

!!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! SUBROUTINE CMP_CS_tmp
!!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!! This is a subroutine which extracts Coherent Structures in flow and save them
!! This subrotine has been called when variables are in fourier space 
!  USE ntypes
!  USE Domain
!  USE Grid
!  USE TIME_STEP_VAR
!  USE FFT_var
!  USE run_variable
!  USE les_chan_var
!  USE variable_stat  
!  USE mpi_var, only : rank

!    IMPLICIT NONE

!    INTEGER I,J,K,ij, JEND_tmp, ZEND_tmp
!    LOGICAL dir_exists
!        
!    REAL*8  Oij(0:NXP,0:NZV-1,0:NY+1,1:3), EIG(0:NXP,0:NZV-1,0:NY+1),QQ(0:NXP,0:NZV-1,0:NY+1),P_prime(0:NXP,0:NZV-1,0:NY+1)
!    COMPLEX*16 COij(0:NX2P,0:NZV-1,0:NY+1,1:3)

!    REAL*8 S11,S12,S13,S22,S23,S33,O12,O13,O23,A11,A12,A13,A22,A23,A33,B,C,D,Q,R,THETA,EIG1,EIG2,EIG3
!    CHARACTER*35 file_name
!    CHARACTER*3 PID

!! Computing the Strain rate tenstor
!     CALL compute_strain_chan

!! Converting pressure back to physical space (it's ignored in flow_statistics.F90)
!      CALL REAL_FOURIER_TRANS_P (.false.)

!!Computing the rotation tensor

!JEND_tmp=100
!ZEND_tmp=160

!      DO J=JSTART,JEND_tmp
!	DO K=ZSTART,ZEND_tmp
!	  DO I=0,NX2P
!  !  1/2(du2/dz-du3/dy) =>
!	    COij(I,K,J,1)=0.5d0*( (0.5*CJOB_11(K+1,J,2)*(CU2X(I,K,J) + CU2X(I,K+1,J)) &
!			- 0.5*CJOB_11(K,J,2)*(CU2X(I,K,J)   + CU2X(I,K-1,J))         &
!			+ 0.5*CJOB_21(K,J+1,1)*(CU2X(I,K,J) + CU2X(I,K,J+1))         &
!			- 0.5*CJOB_21(K,J,1)*(CU2X(I,K,J)   + CU2X(I,K,J-1)) )       &
!			- (0.5*CJOB_12(K+1,J,2)*(CU3X(I,K,J) + CU3X(I,K+1,J))         &
!			- 0.5*CJOB_12(K,J,2)*(CU3X(I,K,J)   + CU3X(I,K-1,J))         &
!			+ 0.5*CJOB_22(K,J+1,1)*(CU3X(I,K,J) + CU3X(I,K,J+1))         &
!			- 0.5*CJOB_22(K,J,1)*(CU3X(I,K,J)   + CU3X(I,K,J-1)) ) )/      &
!			      INT_JACOB(K,J)

!  !  1/2(du1/dz-du3/dx) =>
!	    COij(I,K,J,2)=0.5d0*( ( 0.5*CJOB_11(K+1,J,2)*(CU1X(I,K,J) + CU1X(I,K+1,J)) &
!			- 0.5*CJOB_11(K,J,2)*(CU1X(I,K,J)   + CU1X(I,K-1,J))           &
!			+ 0.5*CJOB_21(K,J+1,1)*(CU1X(I,K,J) + CU1X(I,K,J+1))           &
!			- 0.5*CJOB_21(K,J,1)*(CU1X(I,K,J)   + CU1X(I,K,J-1))  )/       &
!			      INT_JACOB(K,J)                                            &
!			-CIKXP(I)*CU3X(I,K,J) )

!  !  1/2(du1/dy-du2/dx) =>
!	    COij(I,K,J,3)=0.5d0*( ( 0.5*CJOB_12(K+1,J,2)*(CU1X(I,K,J) + CU1X(I,K+1,J))  &
!			- 0.5*CJOB_12(K,J,2)*(CU1X(I,K,J)   + CU1X(I,K-1,J))            &
!			+ 0.5*CJOB_22(K,J+1,1)*(CU1X(I,K,J) + CU1X(I,K,J+1))            &
!			- 0.5*CJOB_22(K,J,1)*(CU1X(I,K,J)   + CU1X(I,K,J-1)) ) /        &
!			      INT_JACOB(K,J)                                             &
!			- CIKXP(I)*CU2X(I,K,J) )

!	  ENDDO
!	ENDDO
!      ENDDO
! 
!  ! Convert rotation tensor to physical space
!      DO ij=1,3
!	CS1X=COij(:,:,:,ij)

!	CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS1X,CS1Z)
!	cvarp(:,:,:)=(0.d0,0.d0)
!	DO I=0,NKX
!	  cvarp(I,:,:)=CS1Z(I,:,:)
!	ENDDO

!	CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
!	DO I=0,NXM
!	  S1Z(I,:,:)=varp(I,:,:)
!	ENDDO
!	S1Z(NXM+1:NXV-1,:,:)=0.0

!	CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)
!	Oij(:,:,:,ij)=S1X
!      ENDDO

!! updateing the terminal points
!    Oij(:,0,:,:)=Oij(:,1,:,:)
!    Oij(:,NZ+1,:,:)=Oij(:,NZ,:,:)

!    Oij(:,:,0,:)=Oij(:,:,1,:)
!    Oij(:,:,NY+1,:)=Oij(:,:,NY,:)

!! computing CS in flow based on Lamda2 critria 
!! for more info please read: "On the identidfication of a vortex" by Jinhee Jeong JFM 1995
!      DO J=JSTART,JEND_tmp
!	DO K=ZSTART,ZEND_tmp
!	  DO I=0,MIN(NXP,NXP_L)
!	    
!            S11 = Sij(I,K,J,1); S12 = Sij(I,K,J,4); S13 = Sij(I,K,J,5)
!            S22 = Sij(I,K,J,2); S23 = Sij(I,K,J,6); S33 = Sij(I,K,J,3)
!            O12 = Oij(I,K,J,3); O13 = Oij(I,K,J,2); O23 = Oij(I,K,J,1)

!	    !     S AND OMEGA ARE COMPUTED, NOW FIND S*S + OMEGA*OMEGA
!	    A11 = S11*S11 + S12*S12 + S13*S13 - O12*O12 - O13*O13
!	    A22 = S12*S12 + S22*S22 + S23*S23 - O12*O12 - O23*O23
!	    A33 = S13*S13 + S23*S23 + S33*S33 - O13*O13 - O23*O23
!	    A12 = S11*S12 + S12*S22 + S13*S23 - O13*O23
!	    A13 = S11*S13 + S12*S23 + S13*S33 + O12*O23
!	    A23 = S12*S13 + S22*S23 + S23*S33 - O12*O13

!	    PI = 4D0*DATAN(1D0)

!	    !      DETERMINE MIDDLE EIGENVALUE
!	    !      COEFFICIENTS OF POLYNOMIAL X^3+B*X^2+C*X+D=0
!	    B = -(A11+A22+A33);
!	    C = -(A12*A12+A13*A13+A23*A23-A11*A22-A11*A33-A22*A33);
!	    D = -(2D0*A12*A13*A23-A11*A23*A23-A22*A13*A13-A33*A12*A12+A11*A22*A33)

!	    !      INTERMEDIATE COEFFICIENTS
!	    Q = (3D0*C - B*B)/9D0; R = (9D0*C*B - 27D0*D - 2D0*B*B*B)/54D0; THETA = DACOS(R/DSQRT(-Q*Q*Q))

!	    !      EIGENVALUES      
!	    EIG1 = 2D0*DSQRT(-Q)*DCOS(THETA/3D0)-B/3D0; EIG2 = 2D0*DSQRT(-Q)*DCOS((THETA+2D0*PI)/3D0)-B/3D0; EIG3 = 2D0*DSQRT(-Q)*DCOS((THETA+4D0*PI)/3D0)-B/3D0
!            
!            ! COMPUTING Seconf invariant of velocity gradient tensor as Q=-0.5*(trace(S**2+Omega**))
!            QQ(I,K,J) = -0.5 * ( EIG1 + EIG2 + EIG3) 

!	    !      IF U2 IS THE MIDDLE EIGENVALUE
!	    IF  ((EIG2.LT.EIG3).AND.(EIG2.GT.EIG1)) THEN
!	     EIG1 = EIG2
!	    ELSE IF ( (EIG2.LT.EIG1).AND.(EIG2.GT.EIG3) ) THEN
!	     EIG1 = EIG2
!	    !      IF U3 IS THE MIDDLE EIGENVALUE
!	    ELSE IF ( (EIG3.LT.EIG1).AND.(EIG3.GT.EIG2) ) THEN
!	     EIG1 = EIG3
!	    ELSE IF ( (EIG3.LT.EIG2).AND.(EIG3.GT.EIG1) ) THEN
!	     EIG1 = EIG3
!	    ENDIF
!		
!	    !     STORE LAMBDA2 VALUES IN F1
!	    EIG(I,K,J) = EIG1

!            ! COMPUTING PRESSURE PRETURBATION
!            P_PRIME(I,K,J) = PX(I,K,J)-dble(P_MEAN(K,J))
!	  ENDDO
!	ENDDO
!      ENDDO


!! saving the CS on different files
!      !check for directory existence
!      INQUIRE(DIRECTORY='./Chr_Strct/.', EXIST=dir_exists) 

!      IF ( .NOT. dir_exists) THEN
!	WRITE(6,*) 'Chr_Strct directory DOes not exists'
!	CALL system('mkdir Chr_Strct')
!	WRITE(6,*) 'Chr_Strct is created!'
!      ENDIF 
!    
!      k = time_step/SAVE_STATS_INT

!      PID = CHAR(MOD(Rank+1,1000)/100+48)     &
!            //CHAR(MOD(Rank+1,100)/10+48)    &
!            //CHAR(MOD(Rank+1,10)+48)

!      file_name = 'Chr_Strct/CS_'  &
!          //CHAR(MOD(k,100000)/10000+48) &
!          //CHAR(MOD(k,10000)/1000+48)   &
!          //CHAR(MOD(k,1000)/100+48)     &
!          //CHAR(MOD(k,100)/10+48)       &
!          //CHAR(MOD(k,10)+48) //        &
!          '.pln_'                        &
!          //PID


!      IF (RANK.EQ.0) WRITE(6,*) 'start writing in pln format ',file_name

!      OPEN(22,file=file_name,form='unformatted',status='unknown')
!      WRITE(22) xpoint(0:ZEND_tmp+1,0:JEND_tmp+1),ypoint(0:ZEND_tmp+1,0:JEND_tmp+1),dx(1),EIG(0:MIN(NXP,NXP_L),0:ZEND_tmp+1,:)
!      WRITE(22) QQ(0:MIN(NXP,NXP_L),0:ZEND_tmp+1,:)
!      WRITE(22) Oij(0:MIN(NXP,NXP_L),0:ZEND_tmp+1,:,1)
!      WRITE(22) Oij(0:MIN(NXP,NXP_L),0:ZEND_tmp+1,:,2)
!      WRITE(22) Oij(0:MIN(NXP,NXP_L),0:ZEND_tmp+1,:,3)
!      WRITE(22) P_PRIME(0:MIN(NXP,NXP_L),0:ZEND_tmp+1,:)
!      CLOSE(22)

!  RETURN
!END
