! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE SAVE_STATS_CURVI(FINAL)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var,		 ONLY :  pi
  USE TIME_STEP_VAR
  USE run_variable
  USE variable_stat
  USE les_chan_var
  USE mpi_var
  USE mg_vari, 		ONLY : INIT_FLAG
      
    implicit none

    LOGICAL FINAL, dir_exists
    integer i,j,k,n
    REAL(r8) ubulk, vbulk, area
    REAL(r8) dummy1(0:NZ+1,0:NY+1),dummy2(0:NZ+1,0:NY+1)
    REAL(r8) dummy3(0:NZ+1,0:NY+1),dummy4(0:NZ+1,0:NY+1), dummy5(1:NZ+1,1:NY+1)
    CHARACTER*29 file_name
    CHARACTER*3 PID
    CHARACTER*5 file_num
    LOGICAL     TKE_BUDGET, SAVE_3D, SAVE_SPAN_XY, SAVE_CS, SAVE_SPAN_XZ, HIGH_MOMENT, SAVE_SPAN_ZY

    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

! C Initialize the random number generator
    CALL RANDOM_SEED(SIZE = K)
    Allocate (seed(1:K))
    seed(1:K)=10
    CALL RANDOM_SEED(PUT = seed)

    
    TKE_BUDGET = .False.
    SAVE_3D    = .True.

    SAVE_SPAN_XY    = .TRUE.
    SAVE_SPAN_XZ    = .TRUE.
    SAVE_SPAN_ZY   = .TRUE.

    SAVE_CS    = .FALSE.   
    HIGH_MOMENT    = .FALSE.   
    
    IF ( rank .eq. 0) THEN
      OPEN(99,file = 'out_screen.txt',form='formatted', status='unknown', position='append')
      WRITE(99,*) 'Saving flow statistics.'


      WRITE(6,*) 'Saving flow statistics.'
      WRITE(6,*) 'Allocate all the tmp arrays'
    ENDIF



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! need to allocate temp variables
    CALL allocate_temps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 
! ! Compute and write out the centerline velocity
!     IF (int(float(NY)/2.) .eq. float(NY)/2.) THEN
! ! IF NY is even
!       IF (int(float(NZ)/2.) .eq. float(NZ)/2.) THEN
! !          uc=dble(CU3X(0,NZ/2,int(float(NY)/2.)))
!       ELSE
! !          uc=0.5*(dble(CU3X(0,int(float(NZ)/2.)-1,NY/2.)   &
! !                +  CU3X(0,int(float(NZ)/2.),NY/2) ))
!       ENDIF  
!     ELSE
! !        uc=0.5*(dble(CU3X(0,NZ/2,int(float(NY)/2.)-1)) &
! !              +dble(CU3X(0,NZ/2,int(float(NY)/2.))))        
!     ENDIF
     
    IF (RANK .eq. 0) THEN
  !      write(*,*) 'Centerline velocity = ', uc 
  ! Compute and write out bulk velocity
  ! Integrat the instantaneous mean profile numerically at GY points
      
      UBULK=0.0d0
      VBULK=0.0d0
      DO J=1,NY
	DO K=1,NZ
	UBULK= UBULK + 0.25d0 * (dble(CU3X(0,K,J))+dble(CU3X(0,K-1,J)) + &
		dble(CU3X(0,K-1,J-1)) + dble(CU3X(0,K,J-1))) 
	
	VBULK= VBULK + 0.25d0 * (dble(CU1X(0,K,J))+dble(CU1X(0,K-1,J)) + &
		dble(CU1X(0,K-1,J-1)) + dble(CU1X(0,K,J-1)))
	ENDDO
      ENDDO
      UBULK = UBULK/real(NY*NZ)  
      VBULK = VBULK/real(NY*NZ)    
  ! Write out UBULK
      WRITE(*,*) 'UBULK: ',UBULK
    ENDIF


! Save CUi
      DO k=0,NZ+1
        DO i=0,NX2P
          DO j=0,NY+1
            CR1X(i,k,j) = CU1X(i,k,j)
            CR2X(i,k,j) = CU2X(i,k,j)
            CR3X(i,k,j) = CU3X(i,k,j)
!   THIS STEPs ARE REQURIED WHEN DERIVATIVES w.r.t X IS REQURIED LATER
            CF1X(i,k,j) = CU1X(i,k,j)
            CF2X(i,k,j) = CU2X(i,k,j)
            CF3X(i,k,j) = CU3X(i,k,j)
          ENDDO
        ENDDO
      ENDDO
!    TRANSFERRING MEAN TO ALL NODES 

      DO k=1,NZ
          DO j=1,NY
            p_mean(k,j) = CPX(0,k,j)
          ENDDO
      ENDDO
      CALL MPI_BCAST_REAL(p_mean ,NZ+2, NY+2)

      DO k=0,NZ+1
          DO j=0,NY+1
            dummy1(k,j)=CR1X(0,k,j)
            dummy2(k,j)=CR2X(0,k,j)
            dummy3(k,j)=CR3X(0,k,j)
          ENDDO
      ENDDO

      CALL MPI_BCAST_REAL(dummy1,NZ+2,NY+2)         
      CALL MPI_BCAST_REAL(dummy2,NZ+2,NY+2)
      CALL MPI_BCAST_REAL(dummy3,NZ+2,NY+2)


      DO k=0,NZ+1
       DO j=0,NY+1
!           IF (W_BC_ZMAX .NE. 6) THEN
            CR1X(0,k,j)=dummy1(k,j)
            CR2X(0,k,j)=dummy2(k,j)
            CR3X(0,k,j)=dummy3(k,j)
!           ELSE
!            CR1X(0,k,j)=SUM(dummy1(1:NZ,j))/dble(NZ)
!            CR2X(0,k,j)=SUM(dummy2(1:NZ,j))/dble(NZ)
!            CR3X(0,k,j)=SUM(dummy3(1:NZ,j))/dble(NZ)
!           ENDIF
       ENDDO
      ENDDO

! Computing CS in flow
      IF (SAVE_CS) CALL CMP_CS 

! Convert to physical space
      
      CALL REAL_FOURIER_TRANS_U1 (.false.)
      CALL REAL_FOURIER_TRANS_U2 (.false.)
      CALL REAL_FOURIER_TRANS_U3 (.false.)
! pressure has been already convereted back in COMP_CS subrotine 
      IF (.NOT. SAVE_CS) CALL REAL_FOURIER_TRANS_P (.false.)


! ! Get the turbulent kinetic energy at each level
     
      DO k=0,NZ+1
        DO j=0,NY+1
         IF (AMP_OMEGA0 .EQ. 0.d0) THEN
	  TKE(k,j)= ( dble(CR1X(0,k,j)) ** 2.0 + dble(CR2X(0,k,j)) ** 2.0 &
                   + dble(CR3X(0,k,j)) ** 2.0)
	 ELSE
	  TKE(k,j)= ( dble(CR1X(0,k,j)) ** 2.0 + dble(CR2X(0,k,j)) ** 2.0 &
                   + dble(CR3X(0,k,j)) ** 2.0) * dble(AMP_OMEGA0/OMEGA0)**(-2.0)
	 ENDIF
        ENDDO
      ENDDO

    
      DO K=0,NZ+1 
        DO J=0,NY+1
          
          urms(k,j)=0.d0
          vrms(k,j)=0.d0
          wrms(k,j)=0.d0
	  
	  DO i=0,min(NXP,NXP_L) 
	    urms(k,j) = urms(k,j) + ( U1X(i,k,j) - dble(CR1X(0,k,j)) ) ** 2.0d0
	    vrms(k,j) = vrms(k,j) + ( U2X(i,k,j) - dble(CR2X(0,k,j)) ) ** 2.0d0
	    wrms(k,j) = wrms(k,j) + ( U3X(i,k,j) - dble(CR3X(0,k,j)) ) ** 2.0d0
	  ENDDO
	  
!        urms(k,j)=dsqrt(urms(k,j)/(dble(NX)))
!        vrms(k,j)=dsqrt(vrms(k,j)/(dble(NX)))
!        wrms(k,j)=dsqrt(wrms(k,j)/(dble(NX)))
	ENDDO 
      ENDDO
      

! Compute the Reynolds stress and mean velocity gradient
      DO k=1,NZ
	DO j=0,NY+1
	  uv(k,j)=0.d0 
	  uw(k,j)=0.d0
	  wv(k,j)=0.d0
	  pv(k,j)=0.d0
	  
	  DO i=0,min(NXP,NXP_L)
	      uv(k,j) = uv(k,j) + ( U1X(i,k,j) - dble( CR1X(0,k,j)) ) * (U2X(i,k,j) - dble( CR2X(0,k,j) ))
	      wv(k,j) = wv(k,j) + ( U3X(i,k,j) - dble( CR3X(0,k,j)) ) * (U2X(i,k,j) - dble( CR2X(0,k,j) )) 
	      uw(k,j) = uw(k,j) + ( U1X(i,k,j) - dble( CR1X(0,k,j)) ) * (U3X(i,k,j) - dble( CR3X(0,k,j) )) 
	      pu(k,j) = pu(k,j) + ( U3X(i,k,j) - dble( CR3X(0,k,j)) ) * (PX(i,k,j) - dble( p_mean(k,j) ))
	      pv(k,j) = pv(k,j) + ( U2X(i,k,j) - dble( CR2X(0,k,j)) ) * (PX(i,k,j) - dble( p_mean(k,j) ))
	  ENDDO
	  
	  uv(k,j) = uv(k,j) / float(NX)
	  uw(k,j) = uw(k,j) / float(NX)
	  wv(k,j) = wv(k,j) / float(NX)
	  pu(k,j) = pu(k,j) / float(NX)
	  pv(k,j) = pv(k,j) / float(NX)
	ENDDO
      ENDDO        

      
! Get the y-derivative of the mean velocity at GYF points
!      do j=1,NY
!        dudy(j)=dble(CR1X(0,0,j+1)-CR1X(0,0,j-1))/(2.*DYF(j))
!        dwdy(j)=dble(CR3X(0,0,j+1)-CR3X(0,0,j-1))/(2.*DYF(j))
!      ENDDO
! Get the y-derivative of the mean velocity at GY points
      DO K=ZSTART,ZEND
	DO J=JSTART,JEND
	  
	  dudy(k,j) = ( 0.125 * (CJOB_12(K,J,1) + CJOB_12(K,J+1,1) + CJOB_12(K,J,2) + CJOB_12(K+1,J,2) )   &
		      * dble ( CR1X(0,k+1,j) - CR1X(0,k-1,j) )            &
		      + 0.125 * (CJOB_22(K,J,1) + CJOB_22(K,J+1,1) + CJOB_22(K,J,2) + CJOB_22(K+1,J,2) )   &
		      * dble ( CR1X(0,k,j+1) - CR1X(0,k,j-1)) )/ INT_JACOB(K,J)
		      
	  dwdy(k,j) = ( 0.125 * (CJOB_12(K,J,1) + CJOB_12(K,J+1,1) + CJOB_12(K,J,2) + CJOB_12(K+1,J,2) )  &
		      * dble ( CR3X(0,k+1,j) - CR3X(0,k-1,j) )           &
		      + 0.125 * (CJOB_22(K,J,1) + CJOB_22(K,J+1,1) + CJOB_22(K,J,2) + CJOB_22(K+1,J,2) )  &
		      * dble ( CR3X(0,k,j+1) - CR3X(0,k,j-1)) ) / INT_JACOB(K,J)
	ENDDO
      ENDDO

      DO k=1,NZ
       dudy(k,0)=dudy(k,1)
       dwdy(k,0)=dwdy(k,1) 
       dudy(k,NY+1)=dudy(k,NY)
       dwdy(k,NY+1)=dwdy(k,NY)       
      ENDDO
      

      DO k=1,NZ              
	DO j=1,NY
	  
	  dudz(k,j) = ( 0.125 * (CJOB_21(K,J,1) + CJOB_21(K,J+1,1) + CJOB_21(K,J,2) + CJOB_21(K+1,J,2) )  &
		      * dble ( CR1X(0,k,j+1) - CR1X(0,k,j-1) )           &
		      + 0.125 * (CJOB_11(K,J,1) + CJOB_11(K,J+1,1) + CJOB_11(K,J,2) + CJOB_11(K+1,J,2) )  &
		      * dble ( CR1X(0,k+1,j) - CR1X(0,k-1,j)) ) / INT_JACOB(K,J)
	  
	  dwdz(k,j) = ( 0.125 * (CJOB_21(K,J,1) + CJOB_21(K,J+1,1) + CJOB_21(K,J,2) + CJOB_21(K+1,J,2) )  &
		      * dble ( CR3X(0,k,j+1) - CR3X(0,k,j-1) )           &
		      + 0.125 * (CJOB_11(K,J,1) + CJOB_11(K,J+1,1) + CJOB_11(K,J,2) + CJOB_11(K+1,J,2) )  &
		      * dble ( CR3X(0,k+1,j) - CR3X(0,k-1,j)) )/ INT_JACOB(K,J)
	ENDDO
      ENDDO
      
      DO k=1,NZ
       dudz(k,0)=dudz(k,1)
       dwdz(k,0)=dwdz(k,1)
       dudz(k,NY+1)=dudz(k,NY)
       dwdz(k,NY+1)=dwdz(k,NY)
      ENDDO

      DO j=1,NY
       dwdz(1,j)=dwdz(2,j)
       dudz(1,j)=dudz(2,j)
      ENDDO
        
! Calculate the mean square shear
!       DO k=1,NZ
! 	DO j=1,NY
! 	  
! 	  shear(k,j)=0.d0
!           
!           DO i=0,min(NXP,NXP_L)
!             shear(k,j)=shear(k,j)         &
! 			+( (U1X(i,k,j+1) - U1X(i,k,j-1) ) * (2.d0 * DYF(j)))** (-2.d0) &
! 			+( (U3X(i,k,j+1) - U3X(i,k,j-1) ) * (2.d0 * DYF(j)))** (-2.d0) &
! 			+( (U1X(i,k+1,j) - U1X(i,k-1,j) ) * (2.d0 * DZF(k)))** (-2.d0) &
! 			+( (U3X(i,k+1,j) - U3X(i,k-1,j) ) * (2.d0 * DZF(k)))** (-2.d0)
!           ENDDO
!         ENDDO
!         
!         shear(k,j)=shear(k,j) / dble(NX)
!       ENDDO


      DO k=1,NZ              
	DO j=1,NY
	
!	 omega_x(k,j) = dwdy(k,j) -  dudz(k,j)
	
 	  omega_x(k,j) = ( 0.125 * (CJOB_12(K,J,1) + CJOB_12(K,J+1,1) + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) ) &
 			  * dble ( CR3X(0,k+1,j) - CR3X(0,k-1,j) )           &
 			  + 0.125 * (CJOB_22(K,J,1) + CJOB_22(K,J+1,1) + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) ) &
 			  * dble ( CR3X(0,k,j+1) - CR3X(0,k,j-1)) ) / INT_JACOB(K,J)                       &
			-( 0.125 * (CJOB_21(K,J,1) + CJOB_21(K,J+1,1) + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) ) &
 			  * dble ( CR2X(0,k,j+1) - CR2X(0,k,j-1) )           &
			  + 0.125 * (CJOB_11(K,J,1) + CJOB_11(K,J+1,1) + CJOB_11(K,J,2) + CJOB_11(K+1,J,2) ) &
			  * dble ( CR2X(0,k+1,j) - CR2X(0,k-1,j)) ) / INT_JACOB(K,J)
!            omega_x(k,j)=( (0.5*CJOB_11(K+1,J,2)*(CR2X(0,K,J) +CR2X(0,K+1,J)) &
!                        - 0.5*CJOB_11(K,J,2)*(CR2X(0,K,J) + CR2X(0,K-1,J))&
!                        + 0.5*CJOB_21(K,J+1,1)*(CR2X(0,K,J) + CR2X(0,K,J+1))&
!                        - 0.5*CJOB_21(K,J,1)*(CR2X(0,K,J)   + CR2X(0,K,J-1)) )&
!                        - (0.5*CJOB_12(K+1,J,2)*(CR3X(0,K,J) + CR3X(0,K+1,J))&
!                        - 0.5*CJOB_12(K,J,2)*(CR3X(0,K,J)   + CR3X(0,K-1,J))&
!                        + 0.5*CJOB_22(K,J+1,1)*(CR3X(0,K,J) + CR3X(0,K,J+1))&
!                        - 0.5*CJOB_22(K,J,1)*(CR3X(0,K,J)   + CR3X(0,K,J-1)) ))/      &
!                              INT_JACOB(K,J)


	ENDDO
      ENDDO
    
      DO k=1,NZ
       omega_x(k,0)   =omega_x(k,1)
       omega_x(k,NY+1)=omega_x(k,NY)
      ENDDO
      
! Convert to physical space
!      call fft_xz_to_physical(CS1,S1X,0,NY+1)
! Get the rms value
!      DO j=1,NY
!      omega_x(j)=0.d0
!      DO k=1,NZM
!      DO i=1,NXP
!        omega_x(j)=omega_x(j)+S1X(i,k,j)**2.d0
!      ENDDO
!      ENDDO
!      omega_x(j)=sqrt(omega_x(j)/(dble(NX-1)*dble(NZ-1)))
!      ENDDO

! Now, get the y-component in fourier space
!      DO j=1,NY
!      DO k=0,TNKZ
!      DO i=0,NX2P
!        CS1(i,k,j)=CIKZ(k)*CR1X(i,k,j)-CIKXP(i)*CR3X(i,k,j)
!      ENDDO
!      ENDDO
!      ENDDO
! Convert to physical space
!      call fft_xz_to_physical(CS1,S1X,0,NY+1)
! Get the rms value
!      DO j=1,NY
!      omega_y(j)=0.d0
!      DO k=0,NZM
!      DO i=0,NXP
!        omega_y(j)=omega_y(j)+S1X(i,k,j)**2.d0
!      ENDDO
!      ENDDO
!      omega_y(j)=sqrt(omega_y(j)/(dble(NX)*dble(NZ)))
!      ENDDO

! Now, get the y-component in fourier space
!      DO j=1,NY
!      DO k=0,TNKZ
!      DO i=0,NX2P
!        CS1(i,k,j)=CIKXP(i)*0.5d0*(CR2X(i,k,j+1)+CR2X(i,k,j))
!     &             -(CR1X(i,k,j+1)-CR1X(i,k,j-1))/(2.d0*DYF(j))
!      ENDDO
!      ENDDO
!        CS1(0,0,j)=CS1(0,0,j)-dudy(j)
!      ENDDO
! Convert to physical space
!      call fft_xz_to_physical(CS1,S1X,0,NY+1)
! Get the rms value
!      DO j=1,NY
!      omega_z(j)=0.d0
!      DO k=0,NZM
!      DO i=0,NXP
!        omega_z(j)=omega_z(j)+S1X(i,k,j)**2.d0
!      ENDDO
!      ENDDO
!      omega_z(j)=sqrt(omega_z(j)/(dble(NX)*dble(NZ)))
!      ENDDO



! Write out the mean statistics at each time as a binary file
      IF ( rank .eq. 0) THEN
      ! check if the folder exists 
	INQUIRE(DIRECTORY='./plane_data/.', EXIST=dir_exists) 
	
	IF ( .NOT. dir_exists) THEN
	  WRITE(6,*) 'Plane_data directory does not exists'
	  CALL system('mkdir plane_data')
	  WRITE(6,*) 'plane_data is created!'
	ENDIF  
	
	
	OPEN(66, file='plane_data/time_bulk.txt', form='formatted', status='unknown', position='append' )
	WRITE(6,*) TIME_STEP,TIME,DELTA_T
	WRITE(66,565) TIME,DELTA_T,UBULK,VBULK,U3X(NXP,2,NY),U1X(NXP,2,NY)
	CLOSE(66)
565	format(6f15.7)
	

	WRITE(99,*) TIME_STEP,TIME,DELTA_T
      ENDIF
 
      count_data = count_data + 1 
      
!       OPEN(60,file='count.txt',form='formatted',status='unknown')
!       WRITE(60,*) count_data
!       CLOSE(60)

501   format(5(F25.9,' '))
 
! c      call tkebudget_chan_1   

! Do over the number of passive scalars
  DO N=1,N_TH

    ! Save CTHX
	  DO k=0,NZ+1
	    DO i=0,NX2P
	      DO j=0,NY+1
		CRTHX(i,k,j,n)=CTHX(i,k,j,n)
	      ENDDO
	    ENDDO
	  ENDDO

	  DO k=0,NZ+1
	    DO j=0,NY+1
	      dummy4(k,j)=CRTHX(0,k,j,N)
	    ENDDO
	  ENDDO

          CALL MPI_BCAST_REAL(dummy4,NZ+2,NY+2)

          DO k=0,NZ+1
           DO j=0,NY+1
!            IF (TH_BC_ZMIN(N) .NE. 6) THEN
             CRTHX(0,k,j,N)=dummy4(k,j)           
!            ELSE
!             CRTHX(0,k,j,N)=SUM(dummy4(1:NZ,j))/dble(NZ)
!            ENDIF
           ENDDO
          ENDDO
    ! Convert to physical space

	  CALL REAL_FOURIER_TRANS_TH (.false.)

	  DO j=0,NY+1
	    DO k=0,NZ+1
	      thrms(k,j,n)=0.
	      DO i=0,min(NXP,NXP_L)
		thrms(k,j,n) = thrms(k,j,n) + ( abs ( THX(i,k,j,n) -dble(CRTHX(0,k,j,n)) ) )**2.0
	      ENDDO
    !        thrms(k,j,n)=sqrt(thrms(k,j,n)/float(NX))
	    ENDDO
	  ENDDO

    ! Get the y-derivative of the mean scalar 
    ! C     dt/dy=zeta_y*(T_zeta) + eta_y*(T_eta)
          DO k=1,NZ
            DO j=1,NY
                dthdy(k,j,n) = ( 0.125 * (CJOB_12(K,J,1) + CJOB_12(K,J+1,1) +CJOB_12(K,J,2) + CJOB_12(K+1,J,2) ) &
                                * dble ( CRTHX(0,k+1,j,n) - CRTHX(0,k-1,j,n))  &
                                + 0.125 * (CJOB_22(K,J,1) + CJOB_22(K,J+1,1) +CJOB_22(K,J,2) + CJOB_22(K+1,J,2) ) &
                                * dble ( CRTHX(0,k,j+1,n) - CRTHX(0,k,j-1,n)) )/ INT_JACOB(K,J)
! adding the background dth/dy
                IF (.NOT. CONT_STRAT) THEN
                  dthdy(k,j,n) =  dthdy(k,j,n) +                          &
                                  ( 0.125 * (CJOB_12(K,J,1) + CJOB_12(K,J+1,1) +CJOB_12(K,J,2) + CJOB_12(K+1,J,2) )   &
                                  * ( TH_BAR(k+1,j)-TH_BAR(k-1,j) )&
                                  + 0.125 * (CJOB_22(K,J,1) + CJOB_22(K,J+1,1) +CJOB_22(K,J,2) + CJOB_22(K+1,J,2) )   &
                                  * ( TH_BAR(k,j+1)-TH_BAR(k,j-1) ) ) / INT_JACOB(K,J)
                ENDIF
            ENDDO
          ENDDO

          DO j=1,NY
            dthdy(0,j,n)    = dthdy(1,j,n)
            dthdy(NZ+1,j,n) = dthdy(NZ,j,n)
          ENDDO

          DO k=0,NZ+1
            dthdy(k,0,n)    = 2.0 * dthdy(k,1,n)-dthdy(k,2,n)
            dthdy(k,NY+1,n) = 2.0 * dthdy(k,NY,n)-dthdy(k,NY-1,n)
          ENDDO
  
    ! Compute the Reynolds stress and mean velocity gradient
	  DO j=0,NY+1
	    DO k=0,NZ+1
		thv(k,j,n)=0.0
		thw(k,j,n)=0.0
		DO i=0,min(NXP,NXP_L)
		  thv(k,j,n) = thv(k,j,n) + ( THX(i,k,j,n) - dble(CRTHX(0,k,j,n))) * ( U2X(i,k,j) - dble(CR2X(0,k,j)) )

		  thw(k,j,n) = thw(k,j,n) + ( THX(i,k,j,n) - dble(CRTHX(0,k,j,n))) * ( U3X(i,k,j) - dble(CR3X(0,k,j)) )
		ENDDO
	      thv(k,j,n) = thv(k,j,n) / float(NX)
	      thw(k,j,n) = thw(k,j,n) / float(NX)
	    ENDDO
	  ENDDO
	  
    ! C     dt/dz=zeta_x*(T_zeta) + eta_x*(T_eta)
	  DO k=1,NZ
	    DO j=1,NY
	      dthdz(k,j,n) =  ( 0.125 * (CJOB_21(K,J,1) + CJOB_21(K,J+1,1) + CJOB_21(K,J,2) + CJOB_21(K+1,J,2) ) &
				* dble ( CRTHX(0,k,j+1,n) - CRTHX(0,k,j-1,n))  &
				+ 0.125 * (CJOB_11(K,J,1) + CJOB_11(K,J+1,1) + CJOB_11(K,J,2) + CJOB_11(K+1,J,2) ) &
				* dble ( CRTHX(0,k+1,j,n) - CRTHX(0,k-1,j,n)) ) / INT_JACOB(K,J)
	    ENDDO
	  ENDDO     

	  DO j=1,NY
	    dthdz(0,j,n)	= dthdz(1,j,n)
	    dthdz(NZ+1,j,n)	= dthdz(NZ,j,n)
	  ENDDO

	  DO k=0,NZ+1
	    dthdz(k,0,n)	= 2.0 * dthdz(k,1,n) - dthdz(k,2,n)
	    dthdz(k,NY+1,n)	= 2.0 * dthdz(k,NY,n) - dthdz(k,NY-1,n)
	  ENDDO
	  
	  UBULK=0.
	  DO J=1,NY
	    DO K=1,NZ
	      UBULK = UBULK + 0.25 * (dble(CRTHX(0,K,J,n))+dble(CRTHX(0,K-1,J,n)) + &
		      dble( CRTHX(0,K-1,J-1,n)) + dble(CRTHX(0,K,J-1,n)) )
	  
	      IF (CONT_STRAT) THEN
		IF( dwdy(k,j) .eq. 0.d0 ) THEN
		  Rig(k,j) = -RI_TAU(N)*(dthdy(k,j,n)-1.d0)*(10.0d8)**2.0
		ELSE 
		  Rig(k,j) = -RI_TAU(N)*(dthdy(k,j,n)-1.d0)/      &
			    (dwdy(k,j)**2.0 + dudy(k,j)**2.0)
		ENDIF
	      ELSE
		IF( dwdy(k,j) .eq. 0.d0 ) THEN
		  Rig(k,j) = -RI_TAU(N)*(dthdy(k,j,n))*(10.0d8)**2.0
		ELSE
		  Rig(k,j) = -RI_TAU(N)*(dthdy(k,j,n))/         & 
			    (dwdy(k,j)**2.0+ dudy(k,j)**2.0)
		ENDIF
	      ENDIF   
		  
	    ENDDO
	  ENDDO
	  
	  UBULK = UBULK / real(NY*NZ) 

	  IF (rank .eq. 0) THEN
	    WRITE(*,*) 'THBULK: ',UBULK, 'Ri', RI_TAU(1)
	    WRITE(99,*) 'THBULK: ',UBULK, 'Ri', RI_TAU(1)
	  ENDIF
    ! ! Compute the potential energy dissipation, grad(THX) \cdot grad(THX)
    ! c      DO j=1,NY
    ! c        pe_diss(j,n)=0.d0
    ! c        DO k=0,NZM
    ! c          DO i=0,NXP
    ! c            pe_diss(j,n)=pe_diss(j,n)
    ! c     &          +R1X(i,k,j)**2.d0+R2X(i,k,j)**2.d0+R3X(i,k,j)**2.d0
    ! c          ENDDO
    ! c        ENDDO
    ! c        pe_diss(j,n)=pe_diss(j,n)/dble(NX*NZ)
    ! c      ENDDO

    !!!!!!!!!!!!SAVING FULL 3D DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (SAVE_3D) THEN
	! check for directory existence
	  INQUIRE(DIRECTORY='./plane_3D/.', EXIST=dir_exists) 

	  IF ( .NOT. dir_exists) THEN
	    WRITE(6,*) 'Plane_3D directory does not exists'
	    CALL system('mkdir plane_3D')
	    WRITE(6,*) 'plane_3D is created!'
	  ENDIF  
	  
	  k = time_step/SAVE_STATS_INT

          PID = CHAR(MOD(Rank+1,1000)/100+48)     &
            //CHAR(MOD(Rank+1,100)/10+48)    &
            //CHAR(MOD(Rank+1,10)+48)

	  file_name = 'plane_3D/data_3d_'  &
		//CHAR(MOD(k,10000)/1000+48)   &
		//CHAR(MOD(k,1000)/100+48)     &
		//CHAR(MOD(k,100)/10+48)       &
		//CHAR(MOD(k,10)+48) //        &
                '.pln_'                        &
                //PID

         IF (RANK.EQ.0) WRITE(6,*) 'start writing in pln format ',file_name
	
        CALL plot_3D(file_name)     
      ENDIF
    !!!!!!!!!!END oF SAVING 3D DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! ! !!!!!!!!!!!!SAVING SPAN-WISE PLANE DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( SAVE_SPAN_XY ) THEN
      
	    I=1+RANK

	    PID = CHAR(MOD(I,1000)/100+48)     &
		  //CHAR(MOD(I,100)/10+48)    &
		  //CHAR(MOD(I,10)+48)
	  
	    k = time_step/SAVE_STATS_INT
	    file_num = CHAR(MOD(k,10000)/1000+48) &
		  //CHAR(MOD(k,1000)/100+48)     &
		  //CHAR(MOD(k,100)/10+48)       &
		  //CHAR(MOD(k,10)+48)//'_'
	  
	    !check for directory existence
	    INQUIRE(DIRECTORY='./xy_plane/.', EXIST=dir_exists) 

	    IF ( .NOT. dir_exists) THEN
	      WRITE(6,*) 'xy_plane directory does not exists'
	      CALL system('mkdir xy_plane')
	      WRITE(6,*) 'xy_plane is created!'
	    ENDIF 
	  
	    k = INT((NZ+1)/4) !you can change this deafult 
	    file_name = 'xy_plane/span1_'//file_num//PID//'.pln'
	    CALL plane_XY_binary(file_name,k)     
	    
	    k = INT((NZ+1)/2) !you can change this deafult 
	    file_name = 'xy_plane/span2_'//file_num//PID//'.pln'
	    CALL plane_XY_binary(file_name,k)     
	    
	    k = INT(3*(NZ+1)/4)   !you can change this deafult 
	    file_name = 'xy_plane/span3_'//file_num//PID//'.pln'
	    CALL plane_XY_binary(file_name,k)     
	    
      ENDIF

      IF ( SAVE_SPAN_ZY ) THEN

            k = time_step/SAVE_STATS_INT
            file_num = CHAR(MOD(k,10000)/1000+48) &
                  //CHAR(MOD(k,1000)/100+48)     &
                  //CHAR(MOD(k,100)/10+48)       &
                  //CHAR(MOD(k,10)+48)//'_'
        
            !check for directory existence
            INQUIRE(DIRECTORY='./zy_plane/.', EXIST=dir_exists)

            IF ( .NOT. dir_exists) THEN
              WRITE(6,*) 'zy_plane directory does not exists'
              CALL system('mkdir zy_plane')
              WRITE(6,*) 'zy_plane is created!'
            ENDIF
       
            !for 1/n'th plane please set rank = INT(NP/(NXP+1)/n) and I
            !=MOD(NX/n,NXP+1)

            IF (RANK.EQ. INT(NX/(NXP+1)/4) ) THEN
             ! the defult for NX/4'th ZY plane 
             k = MOD(NX/4,NXP+1)
             file_name = 'zy_plane/span1_'//file_num//'.pln'
             CALL plane_ZY_binary(file_name,k)
            ENDIF

            IF (RANK.EQ. INT(NX/(NXP+1)/2) ) THEN
             ! the defult for NX/3'th ZY plane 
             k = MOD(NX/2,NXP+1)
             file_name = 'zy_plane/span2_'//file_num//'.pln'
             CALL plane_ZY_binary(file_name,k)
            ENDIF

            IF (RANK.EQ. INT(NX/(NXP+1)/(4/3)) ) THEN
             ! the defult for NX/3'th ZY plane 
             k = MOD(4*NX/3,NXP+1)
             file_name = 'zy_plane/span3_'//file_num//'.pln'
             CALL plane_ZY_binary(file_name,k)
            ENDIF

!            IF (RANK.EQ.1) THEN 
!             !Impilicitly telling processor RANK to write the I'th ZY plane 
!             ! The overal index = I + rank * (NXP+1) 
!             k = 23
!             file_name = 'zy_plane/span1_'//file_num//'.pln'
!             CALL plane_ZY_binary(file_name,k)
!            ENDIF
!
!            IF (RANK.EQ.3) THEN
!             !Impilicitly telling processor RANK to write the I'th ZY plane 
!             ! The overal index = I + rank * (NXP+1) 
!             k = 21
!             file_name = 'zy_plane/span2_'//file_num//'.pln'
!             CALL plane_ZY_binary(file_name,k)
!            ENDIF
!
!            IF (RANK.EQ.5) THEN
!             !Impilicitly telling processor RANK to write the I'th ZY plane 
!             ! The overal index = I + rank * (NXP+1) 
!             k = 19
!             file_name = 'zy_plane/span3_'//file_num//'.pln'
!             CALL plane_ZY_binary(file_name,k)
!            ENDIF
      ENDIF

      IF ( SAVE_SPAN_XZ ) THEN

            I=1+RANK

            PID = CHAR(MOD(I,1000)/100+48)     &
                  //CHAR(MOD(I,100)/10+48)    &
                  //CHAR(MOD(I,10)+48)
        
            k = time_step/SAVE_STATS_INT
            file_num = CHAR(MOD(k,10000)/1000+48) &
                  //CHAR(MOD(k,1000)/100+48)     &
                  //CHAR(MOD(k,100)/10+48)       &
                  //CHAR(MOD(k,10)+48)//'_'
        
            !check for directory existence
            INQUIRE(DIRECTORY='./xz_plane/.', EXIST=dir_exists) 

            IF ( .NOT. dir_exists) THEN
              WRITE(6,*) 'xz_plane directory does not exists'
              CALL system('mkdir xz_plane')
              WRITE(6,*) 'xz_plane is created!'
            ENDIF
        
            J = 13 !you can change this deafult 
            file_name = 'xz_plane/span1_'//file_num//PID//'.pln'
            CALL plane_XZ_binary(file_name,J)
        
            J = 44 !you can change this deafult 
            file_name = 'xz_plane/span2_'//file_num//PID//'.pln'
            CALL plane_XZ_binary(file_name,J)
        
            J = 23   !you can change this deafult 
            file_name = 'xz_plane/span3_'//file_num//PID//'.pln'
            CALL plane_XZ_binary(file_name,J)
        
            J = 51  !you can change this deafult 
            file_name = 'xz_plane/span4_'//file_num//PID//'.pln'
            CALL plane_XZ_binary(file_name,J)

            J = 8  !you can change this deafult 
            file_name = 'xz_plane/span5_'//file_num//PID//'.pln'
            CALL plane_XZ_binary(file_name,J)

            J = 33  !you can change this deafult 
            file_name = 'xz_plane/span6_'//file_num//PID//'.pln'
            CALL plane_XZ_binary(file_name,J)
      ENDIF

  ! ! !!!!!!!!!!END oF SAVING SPAN-WISE PLANE DATA!!!!!!!!!!!!!!!!!!!!!

  
  ! Convert back to Fourier space
      CALL REAL_FOURIER_TRANS_TH (.true.)
    
    !C      call FFT_X_TO_FOURIER(THX(0,0,0,n),CTHX(0,0,0,n),0,NY+1,0,NZ+1)

    ! End do over number of passive scalars, n
  ENDDO


      CALL MPI_COMBINE_STATS(urms,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(vrms,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(wrms,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(uv,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(uw,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(wv,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(pu,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(pv,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(dudz,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(dudy,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(dwdz,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(dwdy,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(omega_x,NZ+2,NY+2)

! Get the bulk rms value

      IF (rank .eq. 0) THEN

	DO k=0,NZ+1
	  DO j=0,NY+1
	    urms(k,j) = dsqrt(urms(k,j) / dble(NX))
	    vrms(k,j) = dsqrt(vrms(k,j) / dble(NX))
	    wrms(k,j) = dsqrt(wrms(k,j) / dble(NX))
	  ENDDO
	ENDDO      
    
	urms_b = 0.d0
	vrms_b = 0.d0
	wrms_b = 0.d0
	area   = 0.d0
	
	DO k=1,NZ
	  DO j=1,NY
	    area  = area + INT_JACOB(K,J)
	    
	    urms_b= urms_b + 0.25d0 * ( urms(k,j) + urms(k,j-1) + urms(k-1,j) + urms(k-1,j-1) ) *INT_JACOB(K,J)
	    
	    vrms_b= vrms_b + 0.25d0 * ( vrms(k,j) + vrms(k,j-1) + vrms(k-1,j) + vrms(k-1,j-1) ) *INT_JACOB(K,J)
	    
	    wrms_b= wrms_b + 0.25d0 * ( wrms(k,j) + wrms(k,j-1) + wrms(k-1,j) + wrms(k-1,j-1) ) *INT_JACOB(K,J)
	    
	  ENDDO
	ENDDO
	
	urms_b = urms_b / area
	vrms_b = vrms_b / area
	wrms_b = wrms_b / area
  

  ! Write out the bulk rms velocity
	WRITE(*,*) '<U_rms_avg>: ',urms_b
	WRITE(*,*) '<V_rms_avg>: ',vrms_b
	WRITE(*,*) '<W_rms_avg>: ',wrms_b

	WRITE(99,*) '<U_rms_avg>: ',urms_b
	WRITE(99,*) '<V_rms_avg>: ',vrms_b
	WRITE(99,*) '<W_rms_avg>: ',wrms_b

      ENDIF 


      Do n =1,N_th
        dummy5(:,:) = thrms(:,:,n)
	CALL MPI_COMBINE_STATS(dummy5,NZ+2,NY+2)
        thrms(:,:,n) = dummy5(:,:)

        dummy5(:,:) = thw(:,:,n)
        CALL MPI_COMBINE_STATS(dummy5,NZ+2,NY+2)
        thw(:,:,n) = dummy5(:,:)
        dummy5(:,:) = thv(:,:,n)
        CALL MPI_COMBINE_STATS(dummy5,NZ+2,NY+2)
        thv(:,:,n) = dummy5(:,:)

!	CALL MPI_COMBINE_STATS(thv(:,:,n),NZ+2,NY+2)
!	CALL MPI_COMBINE_STATS(thw(:,:,n),NZ+2,NY+2)
!       CALL MPI_COMBINE_STATS(dthdz(0,0,n),NZ+2,NY+2)
!       CALL MPI_COMBINE_STATS(dthdy(0,0,n),NZ+2,NY+2)
!       CALL MPI_COMBINE_STATS(Rig,NZ+2,NY+2)

	IF (rank .eq. 0) THEN
	  DO j=0,NY+1
	    DO k=0,NZ+1
	      thrms(k,j,n) = sqrt( thrms(k,j,n) / float(NX) )
	    ENDDO
	  ENDDO
	ENDIF

      ENDDO

!!!!!!!!!!!!!TKE BUDGET!!!!!!!!!!!!!!!!!!!!!
     IF (TKE_BUDGET) THEN
     
      ! check for directory existence 
      INQUIRE(DIRECTORY='./plane_tke/.', EXIST=dir_exists) 

      IF ( .NOT. dir_exists) THEN
	WRITE(6,*) 'Plane_tke directory does not exists'
	CALL system('mkdir plane_tke')
	WRITE(6,*) 'plane_tke is created!'
      ENDIF  
     
      CALL tkebudget_curvi_duct
     ENDIF
!!!!!!!!!!!!END OF TKE BUDGET!!!!!!!!!!!!!!!

!!!!!!!!!!!!!HIGH MOMEMNTUMS!!!!!!!!!!!!!!!!!!!!!
     IF (HIGH_MOMENT) THEN
     
      ! check for directory existence 
      INQUIRE(DIRECTORY='./plane_mmnt/.', EXIST=dir_exists) 

      IF ( .NOT. dir_exists) THEN
	WRITE(6,*) 'Plane_mmnt directory does not exists'
	CALL system('mkdir plane_mmnt')
	WRITE(6,*) 'plane_mmnt is created!'
      ENDIF  
     
      CALL high_oder_momentums_curvi
     ENDIF
!!!!!!!!!!!!HIGH MOMEMNTUMS!!!!!!!!!!!!!!!




! C Convert velocity back to Fourier space
      CALL REAL_FOURIER_TRANS_U1 (.true.) 
      CALL REAL_FOURIER_TRANS_U2 (.true.)
      CALL REAL_FOURIER_TRANS_U3 (.true.)
      CALL REAL_FOURIER_TRANS_P (.true.)

      CF1X = (0.d0,0.d0)
      CF2X = (0.d0,0.d0)
      CF3X = (0.d0,0.d0)     
      
! 
!C       call fft_x_to_fourier(U1X,CU1X,0,NY+1,0,NZ+1)
!C       call fft_x_to_fourier(U2X,CU2X,0,NY+1,0,NZ+1)
!C       call fft_x_to_fourier(U3X,CU3X,0,NY+1,0,NZ+1)
!C       call fft_x_to_fourier(PX,CPX,0,NY+1,0,NZ+1)

!    combine all statistics and send to root==rank-0
      

	  !VELOCITY FIELD STATISTICS
      IF (rank .eq. 0 ) THEN
	k = time_step/SAVE_STATS_INT
	file_name = 'plane_data/data_tec_'    &
	      //CHAR(MOD(k,100000)/10000+48) &
	      //CHAR(MOD(k,10000)/1000+48)   &
	      //CHAR(MOD(k,1000)/100+48)     &
	      //CHAR(MOD(k,100)/10+48)       &
	      //CHAR(MOD(k,10)+48) //        &
	      '.pln'

	    ! writing the plane averaged statistics in binary format
 	CALL  plot_binary(file_name)  
	
	IF ( WRITE_VEL_TECPLOT ) THEN
	  file_name = 'plane_data/data_tec_'    &
			  //CHAR(MOD(k,100000)/10000+48) &
			  //CHAR(MOD(k,10000)/1000+48)   &
			  //CHAR(MOD(k,1000)/100+48)     &
			  //CHAR(MOD(k,100)/10+48)       &
			  //CHAR(MOD(k,10)+48) //        &
			  '.plt'
      
      ! writing the plane averaged statistics in tecplot format
 	  CALL  plot_tecplot(file_name)
	ENDIF  

	WRITE(99,*) 'done save_stats curvi'
	WRITE(*,*) 'done save_stats curvi' 

	CLOSE(99)
      ENDIF

!     need to deallocate all the arrays
      CALL deallocate_temps
!     CCCCCCCCCCCCCCCCCCCCCCCCCC

      IF (rank.EQ.0) THEN
	WRITE(*,*) 'Dealloc tmp: SAVE_STATS_CURVI'
	WRITE(*,*)
	WRITE(*,*) 'Saving Stats is Completed!'
      ENDIF
     
      IF (FINAL) THEN
	WRITE(*,*) 'LAST Saving is done! GOODBYE.'
      ENDIF
      

  RETURN
END

!-------------------------------------C---
SUBROUTINE high_oder_momentums_curvi
!-------------------------------------C--
! NOte, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps
! This subroutine should be called in SAVE_STATS_CHAN after computing
! plane averaged statistics, with the velocity in physical space, and 
! CRi containing the velocity in Fourier space
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var, 		ONLY :  CIKXP
  USE TIME_STEP_VAR
  USE run_variable
  USE mg_vari, 		ONLY : INIT_FLAG
  USE variable_stat
  USE mpi_var

    IMPLICIT NONE

    INTEGER I,J,K,N
    CHARACTER*33 file_name
    REAL*8 u_prime_tmp(0:NZ+1,0:NY+1), tmp(0:NZ+1,0:NY+1), tmp_1
     
     DO J=0,NY+1      
      DO K=0,NZ+1
	UU(K,J)=0.d0;		VV(K,J)=0.d0;
	WW(K,J)=0.d0;		UUU(K,J)=0.d0;
	VVV(K,J)=0.d0;		WWW(K,J)=0.d0;
	WWV(K,J)=0.d0;		UUV(K,J)=0.d0;
	UUUU(K,J)=0.d0;		VVVV(K,J)=0.d0;
	WWWW(K,J)=0.d0;		UUUV(K,J)=0.d0;
	WWWV(K,J)=0.d0;		
		
	DO I=0,MIN(NXP,NXP_L)
	  UU(K,J) = UU(K,J)+(U1X(I,K,J)-dble(CR1X(0,k,j)))**2
	  VV(K,J) = VV(K,J)+(U2X(I,K,J)-dble(CR2X(0,k,j)))**2
	  WW(K,J) = WW(K,J)+(U3X(I,K,J)-dble(CR3X(0,k,j)))**2
	  
	  UUU(K,J) = UUU(K,J)+(U1X(I,K,J)-dble(CR1X(0,k,j)))**3
	  VVV(K,J) = VVV(K,J)+(U2X(I,K,J)-dble(CR2X(0,k,j)))**3
	  WWW(K,J) = WWW(K,J)+(U3X(I,K,J)-dble(CR3X(0,k,j)))**3
	  
	  UUUU(K,J) = UUUU(K,J)+(U1X(I,K,J)-dble(CR1X(0,k,j)))**4
	  VVVV(K,J) = VVVV(K,J)+(U2X(I,K,J)-dble(CR2X(0,k,j)))**4
	  WWWW(K,J) = WWWW(K,J)+(U3X(I,K,J)-dble(CR3X(0,k,j)))**4
	  
	  
	  UUV(K,J) = UUV(K,J)+(U1X(I,K,J)-dble(CR1X(0,k,j)))**2*(U2X(I,K,J)-dble(CR2X(0,k,j)))
	  WWV(K,J) = WWV(K,J)+(U3X(I,K,J)-dble(CR3X(0,k,j)))**2*(U2X(I,K,J)-dble(CR2X(0,k,j)))
	  
	  
	  UUUV(K,J) = UUUV(K,J)+(U1X(I,K,J)-dble(CR1X(0,k,j)))**3*(U2X(I,K,J)-dble(CR2X(0,k,j)))
	  WWWV(K,J) = WWWV(K,J)+(U3X(I,K,J)-dble(CR3X(0,k,j)))**3*(U2X(I,K,J)-dble(CR2X(0,k,j)))
	  
	  
	ENDDO

	UU(K,J)=UU(K,J)/dble(NX); WW(K,J)=WW(K,J)/dble(NX); VV(K,J)=VV(K,J)/dble(NX);
	UUU(K,J)=UUU(K,J)/dble(NX); WWW(K,J)=WWW(K,J)/dble(NX); VVV(K,J)=VVV(K,J)/dble(NX);
	UUUU(K,J)=UUUU(K,J)/dble(NX); WWWW(K,J)=WWWW(K,J)/dble(NX); VVVV(K,J)=VVVV(K,J)/dble(NX);
	
	UUV(K,J)=UUV(K,J)/dble(NX); WWV(K,J)=WWV(K,J)/dble(NX); UUUV(K,J)=UUUV(K,J)/dble(NX); WWWV(K,J)=WWWV(K,J)/dble(NX); 
      ENDDO
    ENDDO

    CALL MPI_COMBINE_STATS(UU,NZ+2,NY+2)
    CALL MPI_COMBINE_STATS(VV,NZ+2,NY+2)
    CALL MPI_COMBINE_STATS(WW,NZ+2,NY+2)
    CALL MPI_COMBINE_STATS(UUU,NZ+2,NY+2)
    CALL MPI_COMBINE_STATS(VVV,NZ+2,NY+2)
    CALL MPI_COMBINE_STATS(WWW,NZ+2,NY+2)
    CALL MPI_COMBINE_STATS(UUUU,NZ+2,NY+2)
    CALL MPI_COMBINE_STATS(VVVV,NZ+2,NY+2)
    CALL MPI_COMBINE_STATS(WWWW,NZ+2,NY+2)
    CALL MPI_COMBINE_STATS(UUV,NZ+2,NY+2)
    CALL MPI_COMBINE_STATS(WWV,NZ+2,NY+2)
    CALL MPI_COMBINE_STATS(UUUV,NZ+2,NY+2)
    CALL MPI_COMBINE_STATS(WWWV,NZ+2,NY+2)

!!!! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   
!!!! Computing stream-wise derivative of high momentums   (dwdz, (dwdz)^2, (dwdz)^3, and so on)

    DO I=0,MIN(NXP,NXP_L)
      tmp(:,:)=U3X(I,:,:)-dble(CR3X(0,:,:))
      DO J=1,NY      
	DO K=1,NZ

	  tmp_1=0.d0
	  tmp_1 = ( 0.125 * (CJOB_11(K,J,1) + CJOB_11(K,J+1,1) + CJOB_11(K,J,2) + CJOB_11(K+1,J,2) ) &
			    * ( tmp(k+1,j) - tmp(k-1,j) )      &
			    +  0.125 * (CJOB_21(K,J,1) + CJOB_21(K,J+1,1) + CJOB_21(K,J,2) + CJOB_21(K+1,J,2) ) &
			    *  ( tmp(k,j+1) - tmp(k,j-1) ) ) / INT_JACOB(K,J)
			    
	  dWdz1(k,j) = dWdz1(k,j) + tmp_1		  
	  dWdz2(k,j) = dWdz2(k,j) + tmp_1**2
	  dWdz3(k,j) = dWdz3(k,j) + tmp_1**3
	  dWdz4(k,j) = dWdz4(k,j) + tmp_1**4
	  dWdz5(k,j) = dWdz5(k,j) + tmp_1**5
	  
	ENDDO 
      ENDDO
      dWdz1(:,:)=dWdz1(:,:)/dble(NX);    dWdz2(:,:)=dWdz2(:,:)/dble(NX);
      dWdz4(:,:)=dWdz4(:,:)/dble(NX);    dWdz3(:,:)=dWdz3(:,:)/dble(NX);
      dWdz5(:,:)=dWdz5(:,:)/dble(NX);
    ENDDO

    CALL MPI_COMBINE_STATS(dWdz1,NZ+2,NY+2)
    CALL MPI_COMBINE_STATS(dWdz2,NZ+2,NY+2)
    CALL MPI_COMBINE_STATS(dWdz3,NZ+2,NY+2)
    CALL MPI_COMBINE_STATS(dWdz4,NZ+2,NY+2)
    CALL MPI_COMBINE_STATS(dWdz5,NZ+2,NY+2)

!!!! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   
    

    IF (RANK .EQ. 0 ) THEN
      k = time_step/SAVE_STATS_INT
      file_name = 'plane_mmnt/data_mmnt_' &
	    //CHAR(MOD(k,100000)/10000+48)  &
	    //CHAR(MOD(k,10000)/1000+48)    &
	    //CHAR(MOD(k,1000)/100+48)      &
	    //CHAR(MOD(k,100)/10+48)        &
	    //CHAR(MOD(k,10)+48) //         &
	    '.plt'

      CALL plot_tec_mmnt(file_name)
    ENDIF

  RETURN
END


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE tkebudget_curvi_duct
! NOte, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps
! This subroutine should be called in SAVE_STATS_CHAN after computing
! plane averaged statistics, with the velocity in physical space, and 
! CRi containing the velocity in Fourier space

  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var, 		ONLY : CIKXP
  USE TIME_STEP_VAR
  USE run_variable
  USE mg_vari, 		ONLY : INIT_FLAG
  USE variable_stat
  USE mpi_var  

    IMPLICIT NONE

    INTEGER i,j,k,n
    CHARACTER*28 file_tke,file_tke_plt



! Define working arrays
!       real*8 tke_1(0:NZ+1,0:NY+1)
!       real*8 tke_2(0:NZ+1,0:NY+1)
!       real*8 tke_2_1(0:NZ+1,0:NY+1)
!       real*8 tke_2_2(0:NZ+1,0:NY+1)
!       real*8 tke_3(0:NZ+1,0:NY+1)
!       real*8 tke_3_1(0:NZ+1,0:NY+1)
!       real*8 tke_3_2(0:NZ+1,0:NY+1)
!       real*8 tke_3_3(0:NZ+1,0:NY+1)
!       real*8 tke_4(0:NZ+1,0:NY+1)
!       real*8 tke_5(0:NZ+1,0:NY+1,1:N_TH)
!       real*8 tke_6_1(0:NZ+1,0:NY+1)
!       real*8 tke_6_1_1(0:NZ+1,0:NY+1)
!       real*8 tke_6_1_2(0:NZ+1,0:NY+1)
!       real*8 tke_6_1_3(0:NZ+1,0:NY+1)
!       real*8 tke_6_2(0:NZ+1,0:NY+1)
!       real*8 tke_6_2_1(0:NZ+1,0:NY+1)
!       real*8 tke_6_2_2(0:NZ+1,0:NY+1)
!       real*8 tke_6_2_3(0:NZ+1,0:NY+1)
!       real*8 tke_7(0:NZ+1,0:NY+1)
!       real*8 S1_mean(0:NZ+1,0:NY+1)
!       real*8 p_mean(0:NZ+1,0:NY+1)
!       real*8 transport(0:NZ+1,0:NY+1)

      

  ! tke_mean defined at GY points
    DO j=0,NY+1
      DO k=0,NZ+1
	tke_mean_old(k,j) = tke_mean(k,j)
	tke_mean(k,j) = 0.5d0 * ( urms(k,j)**2.d0 + vrms(k,j)**2.d0 + wrms(k,j)**2.d0 ) 
	tke_1(k,j) = ( tke_mean(k,j) - tke_mean_old(k,j) ) / ( TIME-TIME_old )
      ENDDO
    ENDDO

    time_old=TIME
!updating ghost cells 
    tke_mean(0,:)=tke_mean(1,:)
    tke_mean_old(0,:)=tke_mean_old(1,:)
    tke_mean(NZ+1,:)=tke_mean(NZ,:)
    tke_mean_old(NZ+1,:)=tke_mean_old(NZ,:)
    
    !     -U2*dtke/dy
    
    DO j=1,NY
      DO k=1,NZ
      tke_2_1(k,j) = -dble(CR2X(0,K,J)) *  &
		( 0.125 * (CJOB_12(K,J,1) + CJOB_12(K,J+1,1) + CJOB_12(K,J,2) + CJOB_12(K+1,J,2) ) * ( tke_mean(k+1,j) - tke_mean(k-1,j) )         &
		+ 0.125 * (CJOB_22(K,J,1) + CJOB_22(K,J+1,1) + CJOB_22(K,J,2) + CJOB_22(K+1,J,2) ) * ( tke_mean(k,j+1) - tke_mean(k,j-1) ) )       &
		/ INT_JACOB(K,J)
      ENDDO
    ENDDO
    
   
!     -U3*dtke/dz

      DO j=1,NY
	DO k=1,NZ
	  tke_2_2(k,j) = -dble(CR3X(0,K,J)) *  &
			  ( 0.125*(CJOB_21(K,J,1) + CJOB_21(K,J+1,1) + CJOB_21(K,J,2) + CJOB_21(K+1,J,2) ) * ( tke_mean(k,j+1) - tke_mean(k,j-1) )        &
			  + 0.125*(CJOB_11(K,J,1) + CJOB_11(K,J+1,1) + CJOB_11(K,J,2) + CJOB_11(K+1,J,2) ) * ( tke_mean(k+1,j) - tke_mean(k-1,j) ) )      &
			  / INT_JACOB(K,J)
	  
	  tke_2(k,j) = tke_2_1(k,j) + tke_2_2(k,j)
	ENDDO
      ENDDO
       


! Get the production at GY points
    DO j=1,NY
      DO k=1,NZ
  !        tke_3(k,j)=-uv(k,j)*dUdy(k,j)-wv(k,j)*dWdy(k,j)-vv(k,j)*dVdy(k,j)  
  !                 -vw(k,j)*dVdz(k,j)-uw(k,j)*dUdz(k,j)-ww(k,j)*dWdz(k,j)
      
      tke_3(k,j) = -(															&
		    + uv(k,j) * (0.125 * (CJOB_12(K,J,1) + CJOB_12(K,J+1,1) + CJOB_12(K,J,2) + CJOB_12(K+1,J,2) )			&
		    * dble ( CR1X(0,k+1,j)-CR1X(0,k-1,j) )                   							&
		    +  0.125 * (CJOB_22(K,J,1)+CJOB_22(K,J+1,1) + CJOB_22(K,J,2) + CJOB_22(K+1,J,2) )                 		&
		    * dble ( CR1X(0,k,j+1)-CR1X(0,k,j-1) ) )										&
		    
		    + uw(k,j) * (0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1) + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )                 	&
		    * dble ( CR1X(0,k,j+1)-CR1X(0,k,j-1) )                   							&
		    + 0.125 * (CJOB_11(K,J,1)+CJOB_11(K,J+1,1) + CJOB_11(K,J,2) + CJOB_11(K+1,J,2) )                		&
		    * dble ( CR1X(0,k+1,j) - CR1X(0,k-1,j) ) ) 									&
		    
		    + vrms(k,j)**2.0 * (0.125*(CJOB_12(K,J,1) + CJOB_12(K,J+1,1) + CJOB_12(K,J,2) + CJOB_12(K+1,J,2) )             &
		    * dble ( CR2X(0,k+1,j) - CR2X(0,k-1,j) )                 							&
		    + 0.125 * (CJOB_22(K,J,1) + CJOB_22(K,J+1,1) + CJOB_22(K,J,2) + CJOB_22(K+1,J,2) )            			&
		    * dble ( CR2X(0,k,j+1) - CR2X(0,k,j-1) ) )									&
		    
		    + wv(k,j) * ( 0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1) + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) ) 			&
		    * dble ( CR2X(0,k,j+1) - CR2X(0,k,j-1) )                 							&
		    + 0.125 * (CJOB_11(K,J,1)+CJOB_11(K,J+1,1) + CJOB_11(K,J,2) + CJOB_11(K+1,J,2) )   				&
		    * dble ( CR2X(0,k+1,j) - CR2X(0,k-1,j) ) ) 									&
		    
		    + wv(k,j) * (0.125 * (CJOB_12(K,J,1) + CJOB_12(K,J+1,1) + CJOB_12(K,J,2) + CJOB_12(K+1,J,2) )			&
		    * dble ( CR3X(0,k+1,j)-CR3X(0,k-1,j) )										&
		    + 0.125 * (CJOB_22(K,J,1)+CJOB_22(K,J+1,1) + CJOB_22(K,J,2) + CJOB_22(K+1,J,2) ) 				&
		    * dble ( CR3X(0,k,j+1) - CR3X(0,k,j-1) ) ) 									&
		    
		    + wrms(k,j)**2.0 * ( 0.125 * (CJOB_21(K,J,1)+CJOB_21(K,J+1,1) + CJOB_21(K,J,2) + CJOB_21(K+1,J,2) )		&
		    * dble ( CR3X(0,k,j+1) - CR3X(0,k,j-1) )                								&
		    + 0.125 * ( CJOB_11(K,J,1)+CJOB_11(K,J+1,1) + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )              			&
		    * dble ( CR3X(0,k+1,j)-CR3X(0,k-1,j) ) ) 										&
		    ) / INT_JACOB(K,J)  
      

      ENDDO
    ENDDO

! Get the components of the production
! Use S1X as a working variable
! <u1*u1*dU1/dz>==0.0

! u1*u2*dU1/dy
! 
! c      do j=1,NY
! c        do k=1,NZ
! c          do i=0,NXP
! c            S1X(i,k,j)=(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)
! c     &          +  CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*
! c     &             dble(CR1X(0,k+1,j)-CR1X(0,k-1,j))
! c     &          +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)
! c     &          +  CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*
! c     &             dble(CR1X(0,k,j+1)-CR1X(0,k,j-1)))/
! c     &             INT_JACOB(K,J)
! c
! c            S1X(i,k,j)=S1X(i,k,j)*(U1(i,k,j)-dble(CR1X(0,k,j)))
! c     &                *(U2(i,k,j)-dble(CR2X(0,k,j)))
! c          ENDDO
! c        ENDDO
! c      ENDDO
! c      do j=1,NY
! c       do k=1,NZ
! c        tke_3_1(k,j)=-SUM(S1X(0:NXP,k,j))/dble(NX)
! c       ENDDO
! c      ENDDO
! c      do k=0,NZ+1 
! c       tke_3_1(k,0)=0.d0
! c       tke_3_1(k,NY+1)=tke_3_1(k,NY)
! c      ENDDO
! 
! ! u1*u3*dU1/dz
! c      do j=1,NY
! c       do k=1,NZ
! c        do i=0,NXP
! c         S1X(i,k,j)=(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1)
! c     &          + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*
! c     &            dble(CR1X(0,k,j+1)-CR1X(0,k,j-1))
! c     &          +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)
! c     &          + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*
! c     &            dble(CR1X(0,k+1,j)-CR1X(0,k-1,j)))/
! c     &             INT_JACOB(K,J)
! c         S1X(i,k,j)=S1X(i,k,j)*(U1(i,k,j)-dble(CR1X(0,k,j)))
! c     &              *(U3(i,k,j)-dble(CR3X(0,k,j)))
! c        ENDDO
! c       ENDDO
! c      ENDDO
! 
! C      do j=1,NY
! C       do k=1,NZ
! C        tke_3_1(k,j)=tke_3_1(k,j)-SUM(S1X(0:NXP,k,j))/dble(NX)
! C       ENDDO
! C      ENDDO
     
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  -ww*dWdz - vvdVdy
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DO j=1,NY
	DO k=1,NZ
	  tke_3_1(k,j)= -( & 
			  +wrms(k,j)**2.0 * ( 0.125 * (CJOB_21(K,J,1) + CJOB_21(K,J+1,1) + CJOB_21(K,J,2) 		&
					    + CJOB_21(K+1,J,2) ) * dble(CR3X(0,k,j+1)-CR3X(0,k,j-1))     		&
					    + 0.125 * (CJOB_11(K,J,1) + CJOB_11(K,J+1,1) + CJOB_11(K,J,2)		& 
					    + CJOB_11(K+1,J,2) ) * dble(CR3X(0,k+1,j)-CR3X(0,k-1,j)) )		&
					    
			  +vrms(k,j)**2.0 * ( 0.125 * (CJOB_12(K,J,1) + CJOB_12(K,J+1,1) + CJOB_12(K,J,2)		&
					      + CJOB_12(K+1,J,2) ) * dble(CR2X(0,k+1,j)-CR2X(0,k-1,j))			&
					      + 0.125 * (CJOB_22(K,J,1) + CJOB_22(K,J+1,1) + CJOB_22(K,J,2)		&
					      +CJOB_22(K+1,J,2) ) * dble(CR2X(0,k,j+1)-CR2X(0,k,j-1)) ) 		&
			  ) / INT_JACOB(K,J)
	ENDDO
      ENDDO
      
! <u2*u1*dU1dx> == 0 

! <u2*u2*dU2/dy>
! c      do j=1,NY
! c        do k=1,NZ
! c          do i=0,NXP
! c            S1X(i,k,j)=(U2(i,k,j)-dble(CR2X(0,k,j)))**2.0
! c     &              *(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)
! c     &          +  CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*
! c     &             dble(CR2X(0,k+1,j)-CR2X(0,k-1,j))
! c     &          +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)
! c     &          +  CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*
! c     &             dble(CR2X(0,k,j+1)-CR2X(0,k,j-1)))/
! c     &             INT_JACOB(K,J)
! c          ENDDO
! c        ENDDO
! c      ENDDO
! 
! c      do j=1,NY
! c       do k=1,NZ
! c        tke_3_2(k,j)=-SUM(S1X(0:NXP,k,j))/dble(NX)
! c       ENDDO
! c      ENDDO
! c      do k=1,NZ+1 
! c       tke_3_2(k,0)=0.d0
! c       tke_3_2(k,NY+1)=tke_3_2(k,NY)
! c      ENDDO
! 
! ! <u2*u3*dU2dz
! c      do j=1,NY
! c       do k=1,NZ
! c        do i=0,NXP
! c         S1X(i,k,j)=(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1)
! c     &                + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*
! c     &            dble(CR2X(0,k,j+1)-CR2X(0,k,j-1))
! c     &          +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)
! c     &                + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*
! c     &             dble(CR2X(0,k+1,j)-CR2X(0,k-1,j)))/
! c     &             INT_JACOB(K,J)
! c         S1X(i,k,j)=S1X(i,k,j)*(U2(i,k,j)-dble(CR2X(0,k,j)))
! c     &                *(U3(i,k,j)-dble(CR3X(0,k,j)))
! c         ENDDO
! c        ENDDO
! c      ENDDO
! 
! c      do j=1,NY
! c       do k=1,NZ
! c        tke_3_2(k,j)=tke_3_2(k,j)-SUM(S1X(0:NXP,k,j))/dble(NX)
! c       ENDDO
! c      ENDDO

      DO j=1,NY
	DO k=1,NZ
	  tke_3_2(k,j)=tke_3(k,j)-tke_3_1(k,j)
	ENDDO
      ENDDO


! <u3*u1*dU3/dx>==0.0

      
! <u3*u2*dU3/dy>

! c      do j=1,NY
! c        do k=1,NZ
! c          do i=0,NXP
! c            S1X(i,k,j)=(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)
! c     &          +  CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*
! c     &             dble(CR3X(0,k+1,j)-CR3X(0,k-1,j))
! c     &          +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)
! c     &          +  CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*
! c     &             dble(CR3X(0,k,j+1)-CR3X(0,k,j-1)))/
! c     &             INT_JACOB(K,J)
! c
! c            S1X(i,k,j)=S1X(i,k,j)*(U3(i,k,j)-dble(CR3X(0,k,j)))
! c     &                *(U2(i,k,j)-dble(CR2X(0,k,j)))
! c          ENDDO
! c        ENDDO
! c      ENDDO
! 
! c      do j=1,NY
! c       do k=1,NZ
! c        tke_3_3(k,j)=-SUM(S1X(0:NXP,k,j))/dble(NX)
! c       ENDDO
! c      ENDDO
! 
! c      do k=0,NZ+1 
! c       tke_3_3(k,0)=0.d0
! c       tke_3_3(k,NY+1)=tke_3_3(k,NY)
! c      ENDDO
! 
! ! <u3*u3*dU3dz
! 
! c      do j=1,NY
! c       do k=1,NZ
! c        do i=0,NXP
! c         S1X(i,k,j)=(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1)
! c     &                + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*
! c     &            dble(CR3X(0,k,j+1)-CR3X(0,k,j-1))
! c     &          +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)
! c     &                + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*
! c     &             dble(CR3X(0,k+1,j)-CR3X(0,k-1,j)))/
! c     &             INT_JACOB(K,J)
! c         S1X(i,k,j)=S1X(i,k,j)*(U3(i,k,j)-dble(CR3X(0,k,j)))
! c     &                *(U3(i,k,j)-dble(CR3X(0,k,j)))
! c         ENDDO
! c        ENDDO
! c      ENDDO
! 
! c      do j=1,NY
! c       do k=1,NZ
! c        tke_3_3(k,j)=tke_3_3(k,j)-SUM(S1X(0:NXP,k,j))/dble(NX)
! c       ENDDO
! c      ENDDO

    DO n=1,N_TH
      DO j=0,NY+1
	DO k=0,NZ+1
	  tke_3_3(k,j)=dble(CRTHX(0,k,j,n))+TH_BAR(K,J) 
	ENDDO
      ENDDO     
    ENDDO

!    viscous diffusion NU*d2tke/dxidxi
    DO j=1,NY
      DO k=1,NZ
	tke_4(k,j)= NU * ( 	&
			  GMAT_22(K,J+1,1) * ( tke_mean(K,J+1) - tke_mean(K,J) )   &
			- GMAT_22(K,J,1)   * ( tke_mean(K,J)  - tke_mean(K,J-1) )    &
			+ GMAT_11(K+1,J,2) * ( tke_mean(K+1,J) - tke_mean(K,J) )   &
			- GMAT_11(K,J,2)   * ( tke_mean(K,J) - tke_mean(K-1,J) )  &
			+ 0.25 * ( GMAT_12(K+1,J,2) * ( tke_mean(K+1,J+1) + tke_mean(K,J+1) - tke_mean(K,J-1) - tke_mean(K+1,J-1) )    &
				  - GMAT_12(K,J,2) * ( tke_mean(K-1,J+1) + tke_mean(K,J+1) - tke_mean(K,J-1) - tke_mean(K-1,J-1)) )  &
			+ 0.25 * ( GMAT_12(K,J+1,1) * (tke_mean(K+1,J+1) + tke_mean(K+1,J) - tke_mean(K-1,J+1) - tke_mean(K-1,J))   &
				  - GMAT_12(K,J,1) * ( tke_mean(K+1,J) + tke_mean(K+1,J-1) - tke_mean(K-1,J) - tke_mean(K-1,J-1))) &
			 )/ INT_JACOB(K,J)            
      ENDDO 
    ENDDO

    DO n=1,N_TH
      DO j=1,NY
	DO k=1,NZ
	  tke_5(k,j,n)=-RI_TAU(n)*thv(k,j,n)
	ENDDO
      ENDDO
    ENDDO
  
! Construct the resolved turbulent transport terms 
! Convert the pressure to physical space for computation of the
! pressure transport
! Get the mean of the pressure
!       do j=0,NY+1
!        do k=0,NZ+1
!         p_mean(k,j)=dble(CF1X(0,k,j))
!        ENDDO
!       ENDDO

    DO j=0,NY+1      
      DO k=0,NZ+1
	transport(k,j)=0.d0
	DO i=0,min(NXP,NXP_L)
	  transport(k,j) = transport(k,j)     &
! Vertical Pressure Transport term:
	      +( U2X(I,K,J)-dble(CR2X(0,k,j)) ) * (PX(I,K,J) - p_mean(k,J)) 
	  ENDDO        
	transport(k,j) = transport(k,j) / dble(NX)
      ENDDO
    ENDDO

    CALL MPI_COMBINE_STATS(transport,NZ+2,NY+2)

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
    
    DO j=1,NY
      DO k=1,NZ
      	tke_6_1_1(k,j) = - (0.125 * (CJOB_12(K,J,1) + CJOB_12(K,J+1,1) + CJOB_12(K,J,2) + CJOB_12(K+1,J,2) ) &
			  * ( transport(k+1,j) - transport(k-1,j) )      &
			  +  0.125 * (CJOB_22(K,J,1) + CJOB_22(K,J+1,1) + CJOB_22(K,J,2) + CJOB_22(K+1,J,2) ) &
			  *  ( transport(k,j+1) - transport(k,j-1) ) ) / INT_JACOB(K,J)
      ENDDO 
    ENDDO
    
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    DO j=0,NY+1      
      DO k=0,NZ+1
	transport(k,j)=0.d0
	DO i=0,min(NXP,NXP_L)
	  transport(k,j) = transport(k,j)  &
! Horizontal Pressure Transport term:
	    +( U3X(I,K,J)-dble(CR3X(0,k,j)) ) * (PX(I,K,J) - p_mean(k,J))
	ENDDO        
	transport(k,j) = transport(k,j) / dble(NX)
      ENDDO
    ENDDO

    CALL MPI_COMBINE_STATS(transport,NZ+2,NY+2)

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
    DO j=1,NY
      DO k=1,NZ
      	tke_6_1_2(k,j) = - (0.125 * (CJOB_11(K,J,1) + CJOB_11(K,J+1,1) + CJOB_11(K,J,2) + CJOB_11(K+1,J,2) ) &
		  * ( transport(k+1,j) - transport(k-1,j) )      &
		  +  0.125 * (CJOB_21(K,J,1) + CJOB_21(K,J+1,1) + CJOB_21(K,J,2) + CJOB_21(K+1,J,2) ) &
		  * ( transport(k,j+1) - transport(k,j-1) ) ) / INT_JACOB(K,J)

	tke_6_1(k,j)=tke_6_1_1(k,j) + tke_6_1_2(k,j)
      ENDDO 
    ENDDO

! Turbulent transport terms:
!   d(0.5*u_i^2.0*v)dy


    DO j=0,NY+1
      DO k=0,NZ+1
	transport(k,j)=0.d0
	  do i=0,min(NXP,NXP_L)
	    transport(k,j) = transport(k,j)             &
			      ! u1^2*u2
			     + 0.5d0 * (U1X(I,K,J)-dble(CR1X(0,k,J)))**2.d0 * (U2X(I,K,J)- dble(CR2X(0,k,J)) )       &
			      ! u2^3
			     + 0.5d0 * (U2X(I,K,J)- dble(CR2X(0,k,J)))**3.d0  		&
			      ! u3^2*u2
			     + 0.5d0 * (U3X(I,K,J)-dble(CR3X(0,k,J)))**2.d0 * (U2X(I,K,J)- dble(CR2X(0,k,J)))
	ENDDO
	
	transport(k,j) = transport(k,j) / dble(NX)
      ENDDO
    ENDDO

    CALL MPI_COMBINE_STATS(transport,NZ+2,NY+2)
      
      
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
! Now, the vertical derivative of the transport term:
    DO j=1,NY
      DO k=1,NZ
	  tke_6_2_1(k,j) = - (0.125 * (CJOB_12(K,J,1) + CJOB_12(K,J+1,1) + CJOB_12(K,J,2) + CJOB_12(K+1,J,2) ) &
		    * ( transport(k+1,j) - transport(k-1,j) )      &
		    +  0.125 * (CJOB_22(K,J,1) + CJOB_22(K,J+1,1) + CJOB_22(K,J,2) + CJOB_22(K+1,J,2) ) &
		    *  ( transport(k,j+1) - transport(k,j-1) ) ) / INT_JACOB(K,J)
      ENDDO
    ENDDO


!   d(0.5*u_i^2.0*w)dz

    DO j=0,NY+1
      DO k=0,NZ+1
	transport(k,j)=0.d0
	DO i=0,min(NXP,NXP_L)
	    transport(k,j) = transport(k,j)   		&
				! u1^2*u3
			      + 0.5d0 * (U1X(I,K,J)-dble(CR1X(0,K,J)))**2.d0 * (U3X(I,K,J)- dble(CR3X(0,K,J)) ) &
				! u3^3
			      + 0.5d0 * (U3X(I,K,J)- dble(CR3X(0,K,J)))**3.d0 		& 
				! u2^2.0*u3
			      + 0.5d0 * (U2X(I,K,J)-dble(CR2X(0,K,J)))**2.d0 * (U3X(I,K,J)- dble(CR3X(0,K,J)))
	ENDDO
	
	transport(k,j) = transport(k,j) / dble(NX)
	
      ENDDO
    ENDDO

    CALL MPI_COMBINE_STATS(transport,NZ+2,NY+2)
   
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
! Now, the horizontal derivative of the transport term:

    DO j=1,NY
      DO k=1,NZ
	tke_6_2_2(k,j) = - (0.125 * (CJOB_11(K,J,1) + CJOB_11(K,J+1,1) + CJOB_11(K,J,2) + CJOB_11(K+1,J,2) ) &
	  * ( transport(k+1,j) - transport(k-1,j) )      &
	  +  0.125 * (CJOB_21(K,J,1) + CJOB_21(K,J+1,1) + CJOB_21(K,J,2) + CJOB_21(K+1,J,2) ) &
	  *  ( transport(k,j+1) - transport(k,j-1) ) ) / INT_JACOB(K,J)

	tke_6_2(k,j) = tke_6_2_1(k,j) + tke_6_2_2(k,j)
      ENDDO
    ENDDO
       
 
! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
    DO j=0,NY+1
      DO k=0,NZ+1
      epsilon(k,j)=0.d0
      ENDDO
    ENDDO

! Store du/dx in CS1
    DO j=1,NY
      DO k=1,NZ
	DO i=0,NX2P
	  !WRONG STEP ... CF2X ... dv/dx
	  CS1X(i,k,j)=CIKXP(i)*CF1X(i,k,j) ! WE HAVE ALSO STORED CU1 IN CF1X
  !        CS1(i,k,j)=CIKXP(i)*CR1X(i,k,j) 
	ENDDO
      ENDDO
    ENDDO
! Convert to physical space

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

!C      call fft_x_to_physical(CS1,S1X,0,NY+1,0,NZ+1)


      DO j=1,NY
	DO k=1,NZ
	  DO i=0,min(NXP,NXP_L)
    !        epsilon(k,j)=epsilon(k,j)+(S1X(i,k,j)**2.0)
	    epsilon(k,j) = epsilon(k,j) + S1X(i,k,j)**2.0
	  ENDDO
	ENDDO
      ENDDO
      
      
! Store dv/dx in CS1
      Do j=0,NY
	DO k=0,NZ
	  DO i=0,NX2P
	    CS1X(i,k,j)=CIKXP(i)*CF2X(i,k,j)   ! WE HAVE ALSO STORED CU2 IN CF2X
    !        CS1(i,k,j)=CIKXP(i)*CR2X(i,k,j)
	  ENDDO
	ENDDO
      ENDDO

! Convert to physical space
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

!C       call fft_x_to_physical(CS1,S1X,0,NY+1,0,NZ+1)
    DO j=0,NY
      DO k=0,NZ
	Do i=0,min(NXP,NXP_L)
	  epsilon(k,j) = epsilon(k,j) + (S1X(i,k,j)**2.0)
	ENDDO
      ENDDO
    ENDDO
      
      
      
! Compute du/dy note remove mean
    DO j=0,NY+1
      DO k=0,NZ+1
	DO i=0,min(NXP,NXP_L)
	   S2X(i,k,j)= U1X(i,k,j)-dble(CR1X(0,k,j))
	ENDDO
      ENDDO
    ENDDO

    DO j=1,NY
      DO k=1,NZ
	DO i=0,min(NXP,NXP_L)
	  S1X(i,k,j) =  ( 0.125 * (CJOB_12(K,J,1) + CJOB_12(K,J+1,1) + CJOB_12(K,J,2) + CJOB_12(K+1,J,2) )      & 
			* ( S2X(i,k+1,j)-S2X(i,k-1,j) )               &
			+ 0.125 * (CJOB_22(K,J,1) + CJOB_22(K,J+1,1) + CJOB_22(K,J,2) + CJOB_22(K+1,J,2) )       &
			* ( S2X(i,k,j+1)-S2X(i,k,j-1)) )/ INT_JACOB(K,J)
	ENDDO
      ENDDO
    ENDDO


    DO j=1,NY
      DO k=1,NZ
	DO i=0,min(NXP,NXP_L)
	  epsilon(k,j)=epsilon(k,j)+(S1X(i,k,j)**2.0)
	ENDDO
      ENDDO
    ENDDO
      
           
! Store dw/dx in CS1
    DO j=1,NY
      DO k=1,NZ
	DO i=0,NX2P
	  CS1X(i,k,j)=CIKXP(i)*CF3X(i,k,j)  ! WE HAVE ALSO STORED CU3 IN CF3X
	ENDDO
      ENDDO
    ENDDO

! Convert to physical space
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

!C      call fft_x_to_physical(CS1,S1X,0,NY+1,0,NZ+1)
    DO j=0,NY+1
      DO k=0,NZ+1
	DO i=0,min(NXP,NXP_L)
	  epsilon(k,j)=epsilon(k,j)+ (S1X(i,k,j)**2.0)
	ENDDO
      ENDDO
    ENDDO

      
! Compute du/dz 
    DO j=1,NY
      DO k=1,NZ
	DO i=0,min(NXP,NXP_L)
	    S2X(i,k,j)= U1X(i,k,j)-dble(CR1X(0,k,j))
	ENDDO
      ENDDO
    ENDDO

    DO j=1,NY
      DO k=1,NZ
	DO i=0,min(NXP,NXP_L)
	  S1X(i,k,j) = ( 0.125 * (CJOB_21(K,J,1) + CJOB_21(K,J+1,1) + CJOB_21(K,J,2) + CJOB_21(K+1,J,2) )	&
			* ( S2X(i,k,j+1)-S2X(i,k,j-1) )                  &
		       + 0.125 * (CJOB_11(K,J,1) + CJOB_11(K,J+1,1) + CJOB_11(K,J,2) + CJOB_11(K+1,J,2) )	&
			* ( S2X(i,k+1,j)-S2X(i,k-1,j)) ) / INT_JACOB(K,J)
	ENDDO
      ENDDO
    ENDDO    

    DO j=1,NY
      DO k=1,NZ
	DO i=0,min(NXP,NXP_L)
	    epsilon(k,j)=epsilon(k,j)+ (S1X(i,k,j)**2.0)
	ENDDO
      ENDDO
    ENDDO

! Compute dv/dy  note remove mean

    DO j=0,NY+1
      DO k=0,NZ+1
	DO i=0,min(NXP,NXP_L)
	  S2X(i,k,j)= U2X(i,k,j)-dble(CR2X(0,k,j))
	ENDDO
      ENDDO
    ENDDO

    DO j=1,NY
      DO k=1,NZ
	DO i=0,min(NXP,NXP_L)
	  S1X(i,k,j) = ( 0.125 * (CJOB_12(K,J,1) + CJOB_12(K,J+1,1) + CJOB_12(K,J,2) + CJOB_12(K+1,J,2) )	&
			* ( S2X(i,k+1,j)-S2X(i,k-1,j) )               &
			+ 0.125 * (CJOB_22(K,J,1) + CJOB_22(K,J+1,1) + CJOB_22(K,J,2) + CJOB_22(K+1,J,2) )	&
			* ( S2X(i,k,j+1)-S2X(i,k,j-1)) ) / INT_JACOB(K,J)
	ENDDO
      ENDDO
    ENDDO

    DO j=1,NY
      DO k=1,NZ
	DO i=0,min(NXP,NXP_L)
	  epsilon(k,j)=epsilon(k,j)+(S1X(i,k,j)**2.0)
	ENDDO
      ENDDO
    ENDDO

! Compute dw/dy, note remove mean
    DO j=0,NY+1
      DO k=0,NZ+1
	DO i=0,min(NXP,NXP_L)
	  S2X(i,k,j)= U3X(i,k,j)-dble(CR3X(0,k,j))
	ENDDO
      ENDDO
    ENDDO

    DO j=1,NY
      DO k=1,NZ
	DO i=0,min(NXP,NXP_L)
	  S1X(i,k,j) = ( 0.125 * (CJOB_12(K,J,1) + CJOB_12(K,J+1,1) + CJOB_12(K,J,2) + CJOB_12(K+1,J,2) )	&
			* ( S2X(i,k+1,j)-S2X(i,k-1,j) )               & 	
			+ 0.125 * (CJOB_22(K,J,1) + CJOB_22(K,J+1,1) + CJOB_22(K,J,2) + CJOB_22(K+1,J,2) ) 	&
			* (S2X(i,k,j+1)-S2X(i,k,j-1)) ) / INT_JACOB(K,J)
	ENDDO
      ENDDO
    ENDDO
     

    DO j=1,NY
      DO k=1,NZ
	DO i=0,min(NXP,NXP_L)
	  epsilon(k,j) = epsilon(k,j) + (S1X(i,k,j)**2.0)
	ENDDO
      ENDDO
    ENDDO

! Store dv/dz in CF1X
    DO j=0,NY+1
      DO k=0,NZ+1
	DO i=0,min(NXP,NXP_L)
	    S2X(i,k,j)= U2X(i,k,j)-dble(CR2X(0,k,j))
	ENDDO
      ENDDO
    ENDDO

    DO j=1,NY
      DO k=1,NZ
	DO i=0,min(NXP,NXP_L)
	  S1X(i,k,j) = ( 0.125 * (CJOB_21(K,J,1) + CJOB_21(K,J+1,1) + CJOB_21(K,J,2) + CJOB_21(K+1,J,2) )	&
			* ( S2X(i,k,j+1)-S2X(i,k,j-1) )                & 
			+ 0.125 * (CJOB_11(K,J,1) + CJOB_11(K,J+1,1) + CJOB_11(K,J,2) + CJOB_11(K+1,J,2) ) &
			* ( S2X(i,k+1,j)-S2X(i,k-1,j) ) ) / INT_JACOB(K,J)

	ENDDO
      ENDDO
    ENDDO 

    DO j=1,NY
      DO k=1,NZ
	DO i=0,min(NXP,NXP_L)
	  epsilon(k,j) = epsilon(k,j) + (S1X(i,k,j)**2.0)
	ENDDO
      ENDDO
    ENDDO

! Store dw/dz in CS1
    DO j=0,NY+1
      DO k=0,NZ+1
	DO i=0,min(NXP,NXP_L)
	    S2X(i,k,j) = U3X(i,k,j)-dble(CR3X(0,k,j))
	ENDDO
      ENDDO
    ENDDO

    DO j=1,NY
      DO k=1,NZ
	DO i=0,min(NXP,NXP_L)
	  S1X(i,k,j) = ( 0.125 * (CJOB_21(K,J,1) + CJOB_21(K,J+1,1) + CJOB_21(K,J,2) + CJOB_21(K+1,J,2) )	&
			* ( S2X(i,k,j+1)-S2X(i,k,j-1) )                 &
			+ 0.125 * (CJOB_11(K,J,1) + CJOB_11(K,J+1,1) + CJOB_11(K,J,2) + CJOB_11(K+1,J,2) )	&
			* ( S2X(i,k+1,j)-S2X(i,k-1,j) ) ) / INT_JACOB(K,J)

	ENDDO
      ENDDO
    ENDDO

    DO j=1,NY
      DO k=1,NZ
	DO i=0,min(NXP,NXP_L)
	    epsilon(k,j) = epsilon(k,j) + (S1X(i,k,j)**2.0)
	ENDDO
      ENDDO
    ENDDO

!!!!!! NEED ADD FOR ALL PROCESSORS!!!!!!!!!!!!
    CALL MPI_COMBINE_STATS(epsilon,NZ+2,NY+2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO j=1,NY
      DO k=1,NZ
	tke_7(k,j) = -NU * epsilon(k,j) / float(NX)
      ENDDO
    ENDDO

    
    DO k=1,NZ
      tke_7(k,0)=2*tke_7(k,1)-2*tke_7(k,0)
      !tke_7(k,0)=tke_7(k,1)
    ENDDO

    

    DO j=0,NY+1
      DO k=0,NZ+1
	DO i=0,NX2P
	  CF1X(I,K,J)=0.
	ENDDO
      ENDDO
    ENDDO
      

    IF (RANK .eq. 0 ) THEN
      k = time_step/SAVE_STATS_INT
      
      file_tke = 'plane_tke/data_tke_'        &
	      //CHAR(MOD(k,100000)/10000+48)  &
	      //CHAR(MOD(k,10000)/1000+48)    &
	      //CHAR(MOD(k,1000)/100+48)      &
	      //CHAR(MOD(k,100)/10+48)        &
	      //CHAR(MOD(k,10)+48) //         &
	      '.pln'
      
      file_tke_plt = 'plane_tke/data_tke_'        &
	      //CHAR(MOD(k,100000)/10000+48)  &
	      //CHAR(MOD(k,10000)/1000+48)    &
	      //CHAR(MOD(k,1000)/100+48)      &
	      //CHAR(MOD(k,100)/10+48)        &
	      //CHAR(MOD(k,10)+48) //         &
	      '.plt'
      CALL plot_tec_tke(file_tke, file_tke_plt)        
    ENDIF


  RETURN 
END
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

SUBROUTINE plot_3D(file_name)

  USE ntypes
  USE Domain
  USE Grid
  USE run_variable, 	ONLY : U1X, U2X, U3X, THX, PX, TH_BAR
  USE mg_vari, 		ONLY : INIT_FLAG
      
    implicit none

    INTEGER  I,J,K
    CHARACTER*29 file_name

    REAL(r8),ALLOCATABLE,DIMENSION(:,:,:,:)   :: THX_TOT

    ALLOCATE (THX_TOT(0:NXP,0:NZV-1,0:NY+1,1:N_TH))

    DO J=0,NY+1
     DO K=0,NZ+1      
       DO I=0,NXP
        THX_TOT(I,K,J,:)=THX(I,K,J,:)+TH_BAR(K,J)
       ENDDO
     ENDDO
    ENDDO

    WRITE(6,*) file_name

    WRITE(6,*) 'Saving full 3D data'

    OPEN(22,file=file_name,form='unformatted',status='unknown')
    WRITE(22) xpoint,ypoint,dx(1),U3X(0:MIN(NXP,NXP_L),0:NZ+1,:)
    WRITE(22) U2X(0:MIN(NXP,NXP_L),0:NZ+1,:)
    WRITE(22) U1X(0:MIN(NXP,NXP_L),0:NZ+1,:)
    WRITE(22) THX_TOT(0:MIN(NXP,NXP_L),0:NZ+1,:,:)
    WRITE(22) PX(0:MIN(NXP,NXP_L),0:NZ+1,:)
    CLOSE(22)

  ! Open the file and write the tecplot datafile header information.
  !

    WRITE(6,*) 'Done with saving full 3D data'

    DEALLOCATE(THX_TOT)

  RETURN
END


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C *** Subroutine to write a .PLN file of the instantaneous data
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

SUBROUTINE plane_XY_binary(file_name,ind)

  USE ntypes
  USE Domain
  USE Grid
  USE run_variable, 	ONLY : U1X,U2X,U3X,THX,PX,TIME,TIME_STEP
  USE mg_vari, 		ONLY : INIT_FLAG
  USE TIME_STEP_VAR, 	ONLY: DELTA_T 
  USE mpi_var, ONLY: RANK

    IMPLICIT NONE

    INTEGER  ind
    INTEGER  i,j,k, NI
    CHARACTER*29 file_name

    k=ind
    NI = min(NXP,NXP_L)
	  
    IF (RANK.EQ.0) WRITE(6,*) 'Saving XZ-plane data on ', file_name
    OPEN(22,file=file_name,form='unformatted',status='unknown')
    WRITE(22) TIME,dx(1),ypoint(k,0:NY+1),U3X(0:NI,k,0:NY+1), &
              U2X(0:NI,k,0:NY+1), U1X(0:NI,k,0:NY+1), &
	      PX(0:NI,k,0:NY+1),THX(0:NI,k,0:NY+1,1)
      
    CLOSE(22)

  RETURN
END subroutine plane_XY_binary


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C *** Subroutine to write a .PLN file of the instantaneous data
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

SUBROUTINE plane_ZY_binary(file_name,ind)

  USE ntypes
  USE Domain
  USE Grid
  USE run_variable,     ONLY : U1X,U2X,U3X,THX,PX,TIME,TIME_STEP
  USE mg_vari,          ONLY : INIT_FLAG
  USE TIME_STEP_VAR,    ONLY: DELTA_T
  USE mpi_var, ONLY: RANK
    
    IMPLICIT NONE

    INTEGER  ind
    INTEGER  i,j,k, NI
    CHARACTER*29 file_name

    I=ind
       
    WRITE(6,*) 'Saving ZY-plane data on ', file_name 
    OPEN(22,file=file_name,form='unformatted',status='unknown')
    WRITE(22) TIME,xpoint(0:NZ+1,0:NY+1),ypoint(0:NZ+1,0:NY+1),U3X(I,0:NZ+1,0:NY+1), &
              U2X(I,0:NZ+1,0:NY+1), U1X(I,0:NZ+1,0:NY+1), &
              PX(I,0:NZ+1,0:NY+1),THX(I,0:NZ+1,0:NY+1,1)

    CLOSE(22)

  RETURN
END

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C *** Subroutine to write a .PLN file of the instantaneous data
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

SUBROUTINE plane_XZ_binary(file_name,ind)

  USE ntypes
  USE Domain
  USE Grid
  USE run_variable,     ONLY : U1X,U2X,U3X,THX,PX,TIME,TIME_STEP
  USE mg_vari,          ONLY : INIT_FLAG
  USE TIME_STEP_VAR,    ONLY: DELTA_T
  USE mpi_var, ONLY: RANK

    IMPLICIT NONE

    INTEGER  ind
    INTEGER  i,j,k, NI
    CHARACTER*29 file_name

    J=ind
    NI = min(NXP,NXP_L)

    IF (RANK.EQ.0) WRITE(6,*) 'Saving XZ-plane data on ', file_name
    OPEN(22,file=file_name,form='unformatted',status='unknown')
    WRITE(22) TIME,dx(1),xpoint(0:NZ+1,J),U3X(0:NI,0:NZ+1,J), &
              U2X(0:NI,0:NZ+1,J), U1X(0:NI,0:NZ+1,J), &
              PX(0:NI,0:NZ+1,J),THX(0:NI,0:NZ+1,J,1)

    CLOSE(22)
RETURN
END


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C *** Subroutine to write a .PLN file 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

SUBROUTINE plot_binary(file_name)

  USE ntypes
  USE Domain
  USE Grid,   ONLY: xpoint, ypoint
  USE run_variable, ONLY : CU1X, CU2X, CU3X, CTHX, CPX , TH_BAR
  USE variable_stat

    implicit none

    INTEGER k

    INTEGER  i,j,imax,jmax,kmax
    INTEGER  debug
    INTEGER  visdouble,disdouble
    CHARACTER*1 nulchar      
    CHARACTER*29 file_name

    REAL(r8),ALLOCATABLE,DIMENSION(:,:)   :: th_wh

    ALLOCATE (th_wh(0:NZ+1,0:NY+1))

    nulchar = char(0)
    debug   = 0
    visdouble = 0
    disdouble = 1
    imax = NZ+2
    jmax = NY+2
    kmax = 1

    DO j=0,NY+1
      DO i=0,NZ+1
!       th_wh(i,j)=dble(CTHX(0,i,j,1)) + (-ypoint(i,j)+ypoint(1,NY+1))
	th_wh(i,j)=dble(CTHX(0,i,j,1))+ TH_BAR(i,j)
      ENDDO
    ENDDO
      
    WRITE(6,*) 'start writing in pln format ',file_name

    OPEN(22,file=file_name,status='unknown',form='unformatted')
    WRITE(22)xpoint,ypoint,dble(CU3X(0,0:NZ+1,0:NY+1)),dble(CU2X(0,0:NZ+1,0:NY+1)),  &
	    dble(CU1X(0,0:NZ+1,0:NY+1)),th_wh,dble(CTHX(0,0:NZ+1,0:NY+1,1)),wrms, &
	    vrms,urms,wv,uv,pv(:,:),thrms(:,:,1),thv(:,:,1),&
	    thw(:,:,1),Rig,dthdy(:,:,1),dthdz(:,:,1),TKE,omega_x,dwdy, &
	    dwdz,dudy, dble(CPX(0,0:NZ+1,0:NY+1)) 

    CLOSE(22)

  RETURN
END subroutine plot_binary

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C *** Subroutine to write a .PLT file (uses TecPlot tecio.a library)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

SUBROUTINE plot_tecplot(file_name)

 USE ntypes
 USE Domain
 USE Grid, ONLY: xpoint, ypoint
 USE run_variable, ONLY : CU1X, CU2X, CU3X, CTHX, CPX, TIME, TH_BAR,W_BC_ZMAX
 USE variable_stat

      implicit none

      INTEGER k,ncycles
      INTEGER  i,j,imax,jmax,kmax
      INTEGER  debug,ier,itot
      INTEGER  tecini,tecdat,teczne,tecend
      INTEGER  visdouble,disdouble
      REAL*8   phase,PI
      CHARACTER*1 nulchar      
      CHARACTER*29 file_name
      CHARACTER(len=33) :: title

      REAL(r8),ALLOCATABLE,DIMENSION(:,:)   :: th_wh, CPX_tmp, CU3X_tmp, CU2X_tmp, CU1X_tmp 
      REAL(r8),ALLOCATABLE,DIMENSION(:,:)   :: THRMS_tmp, THV_tmp, THW_tmp, CTHX_tmp, DTHDY_tmp, DTHDZ_tmp

      ALLOCATE (th_wh(0:NZ+1,0:NY+1)) 
      ALLOCATE (CPX_tmp(0:NZ+1,0:NY+1), CU3X_tmp(0:NZ+1,0:NY+1), CU2X_tmp(0:NZ+1,0:NY+1), CU1X_tmp(0:NZ+1,0:NY+1))
      ALLOCATE (THRMS_tmp(0:NZ+1,0:NY+1), THV_tmp(0:NZ+1,0:NY+1), THW_tmp(0:NZ+1,0:NY+1), CTHX_tmp(0:NZ+1,0:NY+1))
      ALLOCATE (DTHDY_tmp(0:NZ+1,0:NY+1), DTHDZ_tmp(0:NZ+1,0:NY+1) ) 

      PI  = ATAN(1.0)*4.0
      nulchar = char(0)
      debug   = 0
      visdouble = 0
      disdouble = 1
      imax = NZ+2
      jmax = NY+2
      kmax = 1

! update the tmp variables     
      DO j=0,NY+1
	DO i=0,NZ+1
	    th_wh(i,j)=dble(CTHX(0,i,j,1))+ TH_BAR(i,j)
	    CPX_tmp(i,j) = dble((CPX(0,i,j)))
	    CU3X_tmp(i,j) = dble((CU3X(0,i,j)))
	    CU2X_tmp(i,j) = dble((CU2X(0,i,j)))
	    CU1X_tmp(i,j) = dble((CU1X(0,i,j)))
	    THRMS_tmp(i,j) = thrms(i,j,1)
	    THV_tmp(i,j) = thv(i,j,1)
	    THW_tmp(i,j) = thw(i,j,1)
	    CTHX_tmp (i,j) = dble(CTHX(0,i,j,1))
	    DTHDZ_tmp (i,j) = dthdz(i,j,1)
	    DTHDY_tmp (i,j) = dthdy(i,j,1)
	ENDDO
      ENDDO

!  This is update for "Letf and Right" Ghost Cells (you can comment it if you don't need it)
    
    IF (W_BC_ZMAX.EQ.6) THEN
      DO J=0,NY+1
	  urms(0,J) = urms(NZ,J)
	  vrms(0,J) = vrms(NZ,J)
	  wrms(0,J) = wrms(NZ,J)
	  THRMS_tmp(0,J) = THRMS_tmp(NZ,J)
	  wv(0,J) = wv(NZ,J)
	  uw(0,J) = uw(NZ,J)
          uv(0,J) = uv(NZ,J)
	  pv(0,J) = pv(NZ,J)
	  THV_tmp(0,J) = THV_tmp(NZ,J)
	  THW_tmp(0,J) = THW_tmp(NZ,J)
	  Rig(0,J) = Rig(NZ,J)
          dudy(0,J) = dudy(NZ,J)
          dwdy(0,J) = dwdy(NZ,J)
	  
	  vrms(NZ+1,J) = vrms(1,J)
	  urms(NZ+1,J) = urms(1,J)
	  wrms(NZ+1,J) = wrms(1,J)
	  THRMS_tmp(NZ+1,J) = THRMS_tmp(1,J)
	  wv(NZ+1,J) = wv(1,J)
	  uw(NZ+1,J) = uw(1,J)
          uv(NZ+1,J) = uv(1,J)
	  pv(NZ+1,J) = pv(1,J)
	  THV_tmp(NZ+1,J) = THV_tmp(1,J)
	  THW_tmp(NZ+1,J) = THW_tmp(1,J)
	  Rig(NZ+1,J) = Rig(1,J)
          dudy(NZ+1,J) = dudy(1,J)
          dwdy(NZ+1,J) = dwdy(1,J)
	  
      ENDDO
    ENDIF
  
  !  This is update for "TOP" Ghost Cells (you can comment it if you don't need it)    
      DO I=0,NZ+1
	  urms(I,NY+1) = urms(I,NY)
	  vrms(I,NY+1) = vrms(I,NY)
	  wrms(I,NY+1) = wrms(I,NY)
	  THRMS_tmp(I,NY+1) = THRMS_tmp(I,NY)
	  wv(I,NY+1) = wv(I,NY)
	  uw(I,NY+1) = uw(I,NY)
	  pv(I,NY+1) = pv(I,NY)
	  THV_tmp(I,NY+1) = THV_tmp(I,NY)
	  THW_tmp(I,NY+1) = THW_tmp(I,NY)
	  Rig(I,NY+1) = Rig(I,NY)
      ENDDO
 
      WRITE(6,*) 'start writing in plt format',file_name
!
! Open the file and write the tecplot datafile header information.
!
      ncycles = time/(2.d0*pi)
      phase = (time-ncycles*2.d0*pi)*180.d0/pi
      
      WRITE(*,'(a5,f15.5,a7,f5.1)') 'time=',time,' phase=', phase
      WRITE(title,'(a5,f15.5,a7,f5.1)') 'time=',time,' phase=', phase
!      
!       WRITE(title,'(a6,f13.6)') 'time =',TIME
!
       ier = tecini(trim(title)//nulchar,&
                  &'x,z,ume,wme,vme,thme,thme_d,urms,wrms,vrms,uw,&
                   &vw,pw,thrms,thw,thu,Ri_g,dthdz,dthdx,MKE,omega_y,&
                   &dudz,dudx,dvdz,pme'&
                   &//nulchar,&
                  &file_name//nulchar,&
                  &'.'//nulchar,&
                  &debug,visdouble)
                 

!
! Write the zone header information.
!
      ier = teczne(trim(title)//nulchar, &
                 imax,jmax,kmax,        &
                 'BLOCK'//nulchar,nulchar)

!
! Write out the field data.
!
     itot = imax*jmax*kmax
     ier = tecdat(itot,xpoint,disdouble)
     ier = tecdat(itot,ypoint,disdouble)
     ier = tecdat(itot,CU3X_tmp,disdouble)
     ier = tecdat(itot,CU2X_tmp,disdouble)
     ier = tecdat(itot,CU1X_tmp,disdouble)
     ier = tecdat(itot,th_wh,disdouble)
     ier = tecdat(itot,CTHX_tmp,disdouble)
     ier = tecdat(itot,wrms,disdouble)
     ier = tecdat(itot,vrms,disdouble)
     ier = tecdat(itot,urms,disdouble)
     ier = tecdat(itot,wv,disdouble)
     ier = tecdat(itot,uv,disdouble)
     ier = tecdat(itot,pv(:,:),disdouble)
     ier = tecdat(itot,THRMS_tmp,disdouble)
     ier = tecdat(itot,THV_tmp,disdouble)
     ier = tecdat(itot,THW_tmp,disdouble)
     ier = tecdat(itot,Rig,disdouble)
     ier = tecdat(itot,DTHDY_tmp,disdouble)
     ier = tecdat(itot,DTHDZ_tmp,disdouble)
     ier = tecdat(itot,TKE,disdouble)
     ier = tecdat(itot,omega_x,disdouble)
     ier = tecdat(itot,dwdy,disdouble)
     ier = tecdat(itot,dwdz,disdouble)
     ier = tecdat(itot,dudy,disdouble)
     ier = tecdat(itot,CPX_tmp,disdouble)

!
! close file
!
     ier = tecend()

     DEALLOCATE (th_wh, CPX_tmp, CU3X_tmp, CU2X_tmp, CU1X_tmp)
     DEALLOCATE (THRMS_tmp, THV_tmp, THW_tmp, CTHX_tmp, DTHDY_tmp, DTHDZ_tmp)

  RETURN
END


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C *** Subroutine to write a .plt file (uses TecPlot tecio.a library)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE plot_tec_mmnt(file_name)
  USE ntypes
  USE Domain
  USE Grid, ONLY: xpoint, ypoint
  USE run_variable, 	ONLY : CR2X, CR3X, CTHX, TIME, TH_BAR
  USE variable_stat

      IMPLICIT NONE

      INTEGER  k,ncycles
      INTEGER  i,j,imax,jmax,kmax
      INTEGER  debug,ier,itot
      INTEGER  tecini,tecdat,teczne,tecend
      INTEGER  visdouble,disdouble
      REAL*8   phase,PI
      CHARACTER*1 nulchar      
      CHARACTER(len=66) :: title   
      CHARACTER*30 file_name

      REAL(r8),ALLOCATABLE,DIMENSION(:,:)   :: th_wh, CU3X_tmp, CU2X_tmp
      
      ALLOCATE (CU3X_tmp(0:NZ+1,0:NY+1), CU2X_tmp(0:NZ+1,0:NY+1))
      ALLOCATE (th_wh(0:NZ+1,0:NY+1))

      PI  = ATAN(1.0)*4.0
      nulchar = char(0)
      debug   = 0
      visdouble = 0
      disdouble = 1
      imax = NZ+2
      jmax = NY+2
      kmax = 1

      DO j=0,NY+1
	DO i=0,NZ+1
	  th_wh(i,j)=dble(CTHX(0,i,j,1))+ TH_BAR(i,j)
	  CU3X_tmp(i,j) = dble((CR3X(0,i,j)))
	  CU2X_tmp(i,j) = dble((CR2X(0,i,j)))
	ENDDO
      ENDDO
!      
! Open the file and write the tecplot datafile header information.
!
 
      WRITE(6,*) 'writing high order momentums: ', file_name

      ncycles = time/(2.d0*pi)
      phase = (time-ncycles*2.d0*pi)*180.d0/pi
      
      WRITE(*,'(a5,f15.8,a7,f5.1)') 'time=',time,' phase=', phase
      WRITE(title,'(a5,f15.8,a7,f5.1)') 'time=',time,' phase=', phase
      !
      ! Write the zone header information.
      !
     ier = tecini(trim(title)//nulchar,'x,z,ume,wme,thme,uu,vv,ww,uuu,vvv,www, & 
             &uuuu,vvvv,wwww,uuw,vvw,uuuw,vvvw,dudx,(dudx)^2,(dudx)^3,(dudx)^4,(dudx)^5'&
                  &//nulchar,&
                 &file_name//nulchar,&
                 &'.'//nulchar,&
                 &debug,visdouble)

     ier = teczne(trim(title)//nulchar,imax,jmax,kmax,'BLOCK'//nulchar,nulchar)

      
     ! Write out the field data.
     itot = imax*jmax*kmax
     ier = tecdat(itot,xpoint,disdouble)
     ier = tecdat(itot,ypoint,disdouble)
     ier = tecdat(itot,CU3X_tmp,disdouble)
     ier = tecdat(itot,CU2X_tmp,disdouble)
     ier = tecdat(itot,th_wh,disdouble)
     ier = tecdat(itot,WW,disdouble)
     ier = tecdat(itot,VV,disdouble)
     ier = tecdat(itot,UU,disdouble)
     ier = tecdat(itot,WWW,disdouble)
     ier = tecdat(itot,VVV,disdouble)
     ier = tecdat(itot,UUU,disdouble)
     ier = tecdat(itot,WWWW,disdouble)
     ier = tecdat(itot,VVVV,disdouble)
     ier = tecdat(itot,UUUU,disdouble)
     ier = tecdat(itot,WWV,disdouble)
     ier = tecdat(itot,UUV,disdouble)
     ier = tecdat(itot,WWWV,disdouble)
     ier = tecdat(itot,UUUV,disdouble)
     ier = tecdat(itot,dWdz1,disdouble)
     ier = tecdat(itot,dWdz2,disdouble)
     ier = tecdat(itot,dWdz3,disdouble)
     ier = tecdat(itot,dWdz4,disdouble)
     ier = tecdat(itot,dWdz5,disdouble)
     
     ! close file
     ier = tecend()
      
     DEALLOCATE (th_wh, CU3X_tmp, CU2X_tmp)

  RETURN
END


  
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C *** Subroutine to write a .plt file (uses TecPlot tecio.a library)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE plot_tec_tke(file_tke, file_tke_plt)
  USE ntypes
  USE Domain
  USE Grid, ONLY: xpoint, ypoint
  USE run_variable, 	ONLY : CR2X, CR3X, CTHX, TIME, TH_BAR
  USE variable_stat

      implicit none

      INTEGER  k,ncycles
      INTEGER  i,j,imax,jmax,kmax
      INTEGER  debug,ier,itot
      INTEGER  tecini,tecdat,teczne,tecend
      INTEGER  visdouble,disdouble
      REAL*8   phase,PI
      CHARACTER*1 nulchar      
      CHARACTER(len=33) :: title
      CHARACTER*28 file_tke, file_tke_plt
      
      REAL(r8),ALLOCATABLE,DIMENSION(:,:)   :: th_wh, CU3X_tmp, CU2X_tmp
      
      ALLOCATE (CU3X_tmp(0:NZ+1,0:NY+1), CU2X_tmp(0:NZ+1,0:NY+1))
      ALLOCATE (th_wh(0:NZ+1,0:NY+1))

      PI  = ATAN(1.0)*4.0
      nulchar = char(0)
      debug   = 0
      visdouble = 0
      disdouble = 1
      imax = NZ+2
      jmax = NY+2
      kmax = 1

      DO j=0,NY+1
	DO i=0,NZ+1
	  th_wh(i,j)=dble(CTHX(0,i,j,1))+ TH_BAR(i,j)
	  CU3X_tmp(i,j) = dble((CR3X(0,i,j)))
	  CU2X_tmp(i,j) = dble((CR2X(0,i,j)))
	ENDDO
      ENDDO

!  This is update for "Letf and Right" Ghost Cells (you can comment it if you don't need it)
    
      DO J=0,NY+1
	  CU3X_tmp(0,J) = CU3X_tmp(1,J)
	  CU2X_tmp(0,J) = CU2X_tmp(1,J)
	  th_wh(0,J) = th_wh(1,J)
	  tke_mean(0,J) = tke_mean(1,J)
	  tke_1(0,J) = tke_1(1,J)
	  tke_2(0,J) = tke_2(1,J)
	  tke_2_1(0,J) = tke_2_1(1,J)
	  tke_2_2(0,J) = tke_2_2(1,J)
	  tke_3(0,J) = tke_3(1,J)
	  tke_3_1(0,J) = tke_3_1(1,J)
	  tke_3_2(0,J) = tke_3_2(1,J)
	  tke_3_3(0,J) = tke_3_3(1,J)
	  tke_4(0,J) = tke_4(1,J)
	  tke_5(0,J, 1) = tke_5(1,J, 1)
	  tke_6_1(0,J) = tke_6_1(1,J)
	  tke_6_1_1(0,J) = tke_6_1_1(1,J)
	  tke_6_1_2(0,J) = tke_6_1_2(1,J)
	  tke_6_2(0,J) = tke_6_2(1,J)
	  tke_6_2_1(0,J) = tke_6_2_1(1,J)
	  tke_6_2_2(0,J) = tke_6_2_2(1,J)
 	  tke_7(0,J) = tke_7(1,J) 

	  CU3X_tmp(NZ+1,J) = CU3X_tmp(NZ,J)
	  CU2X_tmp(NZ+1,J) = CU2X_tmp(NZ,J)
	  th_wh(NZ+1,J) = th_wh(NZ,J)
	  tke_mean(NZ+1,J) = tke_mean(NZ,J)
	  tke_1(NZ+1,J) = tke_1(NZ,J)
	  tke_2(NZ+1,J) = tke_2(NZ,J)
	  tke_2_1(NZ+1,J) = tke_2_1(NZ,J)
	  tke_2_2(NZ+1,J) = tke_2_2(NZ,J)
	  tke_3(NZ+1,J) = tke_3(NZ,J)
	  tke_3_1(NZ+1,J) = tke_3_1(NZ,J)
	  tke_3_2(NZ+1,J) = tke_3_2(NZ,J)
	  tke_3_3(NZ+1,J) = tke_3_3(NZ,J)
	  tke_4(NZ+1,J) = tke_4(NZ,J)
	  tke_5(NZ+1,J, 1) = tke_5(NZ,J, 1)
	  tke_6_1(NZ+1,J) = tke_6_1(NZ,J)
	  tke_6_1_1(NZ+1,J) = tke_6_1_1(NZ,J)
	  tke_6_1_2(NZ+1,J) = tke_6_1_2(NZ,J)
	  tke_6_2(NZ+1,J) = tke_6_2(NZ,J)
	  tke_6_2_1(NZ+1,J) = tke_6_2_1(NZ,J)
	  tke_6_2_2(NZ+1,J) = tke_6_2_2(NZ,J)
 	  tke_7(NZ+1,J) = tke_7(NZ,J) 
      ENDDO
      
!  This is update for "TOP" Ghost Cells (you can comment it if you don't need it)
      DO I=0,NZ+1
	  CU3X_tmp(I,NY+1) = CU3X_tmp(I,NY)
	  CU2X_tmp(I,NY+1) = CU2X_tmp(I,NY)
	  th_wh(I,NY+1) = th_wh(I,NY)
	  tke_mean(I,NY+1) = tke_mean(I,NY)
	  tke_1(I,NY+1) = tke_1(I,NY)
	  tke_2(I,NY+1) = tke_2(I,NY)
	  tke_2_1(I,NY+1) = tke_2_1(I,NY)
	  tke_2_2(I,NY+1) = tke_2_2(I,NY)
	  tke_3(I,NY+1) = tke_3(I,NY)
	  tke_3_1(I,NY+1) = tke_3_1(I,NY)
	  tke_3_2(I,NY+1) = tke_3_2(I,NY)
	  tke_3_3(I,NY+1) = tke_3_3(I,NY)
	  tke_4(I,NY+1) = tke_4(I,NY)
	  tke_5(I,NY+1, 1) = tke_5(I,NY, 1)
	  tke_6_1(I,NY+1) = tke_6_1(I,NY)
	  tke_6_1_1(I,NY+1) = tke_6_1_1(I,NY)
	  tke_6_1_2(I,NY+1) = tke_6_1_2(I,NY)
	  tke_6_2(I,NY+1) = tke_6_2(I,NY)
	  tke_6_2_1(I,NY+1) = tke_6_2_1(I,NY)
	  tke_6_2_2(I,NY+1) = tke_6_2_2(I,NY)
 	  tke_7(I,NY+1) = tke_7(I,NY)  
      ENDDO

!  This is update for "LOWER" Ghost Cells (you can comment it if you don't need
!  it)
      DO I=0,NZ+1
          tke_2(I,0) = 0.d0
          tke_2_1(I,0) = 0.d0
          tke_2_2(I,0)= 0.d0
          tke_3(I,0) = 0.d0
          tke_3_1(I,0) = 0.d0
          tke_3_2(I,0) = 0.d0
          tke_4(I,0) = 2*tke_4(I,1)-tke_4(I,2)
          tke_6_1(I,0) = 0.d0
          tke_6_1_1(I,0) = 0.d0
          tke_6_1_2(I,0) = 0.d0
          tke_6_2(I,0) = 0.d0
          tke_6_2_1(I,0) = 0.d0
          tke_6_2_2(I,0) = 0.d0
          tke_7(I,0) = 2*tke_7(I,1)-tke_7(I,2)
      ENDDO


!
! Open the file and write the tecplot datafile header information.
!
 
      WRITE(6,*) 'writing TKE budget: ', file_tke_plt!,file_tke
      OPEN(23,file=file_tke,status='unknown',form='unformatted')
!        write(23)xpoint,ypoint,dble(CR3X(0,0:NZ+1,0:NY+1)),dble(CR2X(0,0:NZ+1,0:NY+1)), &
      WRITE(23)xpoint,ypoint,CU3X_tmp, CU2X_tmp,th_wh, &
                tke_mean,tke_1,tke_2,tke_2_1,tke_2_2,tke_3,tke_3_1, &
                tke_3_2,tke_3_3,tke_4,tke_5(:,:,1),tke_6_1,tke_6_1_1, &
                tke_6_1_2,tke_6_2,tke_6_2_1,tke_6_2_2,tke_7
      
      CLOSE(23)

      ncycles = time/(2.d0*pi)
      phase = (time-ncycles*2.d0*pi)*180.d0/pi
      
      WRITE(*,'(a5,f15.5,a7,f5.1)') 'time=',time,' phase=', phase
      WRITE(title,'(a5,f15.5,a7,f5.1)') 'time=',time,' phase=', phase
      !
      ! Write the zone header information.
      !
     ier = tecini(trim(title)//nulchar,'x,z,ume,wme,thme,tke_m,dkdt,& 
             &advec,advec_y,advec_z,Prod,Prod_u,&
             &Prod_v, Prod_w,vis_diff,&
             &buoy_f,Pr_trns,Pr_trns_y,&
             &Pr_trns_z,tur_trns,&
             &tur_trns_y,tur_trns_z,dissip'&
                  &//nulchar,&
                 &file_tke_plt//nulchar,&
                 &'.'//nulchar,&
                 &debug,visdouble)

     ier = teczne(trim(title)//nulchar,  &
                 imax,jmax,kmax,         &
                 'BLOCK'//nulchar,nulchar)

      
     ! Write out the field data.
     itot = imax*jmax*kmax
     ier = tecdat(itot,xpoint,disdouble)
     ier = tecdat(itot,ypoint,disdouble)
     ier = tecdat(itot,CU3X_tmp,disdouble)
     ier = tecdat(itot,CU2X_tmp,disdouble)
     ier = tecdat(itot,th_wh,disdouble)
     ier = tecdat(itot,tke_mean,disdouble)
     ier = tecdat(itot,tke_1,disdouble)
     ier = tecdat(itot,tke_2,disdouble)
     ier = tecdat(itot,tke_2_1,disdouble)
     ier = tecdat(itot,tke_2_2,disdouble)
     ier = tecdat(itot,tke_3,disdouble)
     ier = tecdat(itot,tke_3_1,disdouble)
     ier = tecdat(itot,tke_3_2,disdouble)
     ier = tecdat(itot,tke_3_3,disdouble)
     ier = tecdat(itot,tke_4,disdouble)
     ier = tecdat(itot,tke_5(:,:,1),disdouble)
     ier = tecdat(itot,tke_6_1,disdouble)
     ier = tecdat(itot,tke_6_1_1,disdouble)
     ier = tecdat(itot,tke_6_1_2,disdouble)
     ier = tecdat(itot,tke_6_2,disdouble)
     ier = tecdat(itot,tke_6_2_1,disdouble)
     ier = tecdat(itot,tke_6_2_2,disdouble)
     ier = tecdat(itot,tke_7,disdouble)

     ! close file
     ier = tecend()
      
     DEALLOCATE (th_wh, CU3X_tmp, CU2X_tmp)
      
  RETURN
END
