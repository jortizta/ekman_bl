SUBROUTINE les_chan_th(n)
!C This subroutine models the subgridscale terms 
!c in the scalar advection equation for scalar number n
!C if the computation is to be treated as an LES not a DNS
!C This subroutine should be CALLed when the velocity is in fourier space 
!C   in the periodic directions
!C S1 should contain |S| which was calculated in les_chan

  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
  USE les_chan_var
  USE mpi_var, ONLY : rank
  
    IMPLICIT NONE

    INTEGER i,j,k,ij,n

! Variables for Dynamic Smagorinsky model:
    REAL*8 C_SMAG
    PARAMETER (C_SMAG=0.13d0)
    LOGICAL dir_exists
    REAL*8 alpha,beta 
    

! Here, alpha is the test/LES filter width ratio
    PARAMETER (alpha=2.44950)
! beta is the LES/grid filter width ratio
    PARAMETER (beta=1.d0)
    CHARACTER*30   file_name

    I = 1
    N=1
! First, for all models, apply boundary conditions to the velocity field
! (fill ghost cells) to ensure accurate calculation of gradients
    J1=JSTART
    J2=JEND 


    IF (LES_MODEL_TYPE_TH .EQ. 1) THEN
!     Constant Smagorinsky model

! First, compute the rate of strain tensor S_ij

      CALL compute_scalar_grad(n)

! Get the eddy diffusivity at GY points
! Kappa_T = NU_T/PR
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
!            KAPPA_T(I,K,J,N)=-1.d0*DELTA_Y(K,J)*S1(I,K,J)
             KAPPA_T(I,K,J,N) = NU_T(I,K,J)/PR(1)
          ENDDO
        ENDDO
      ENDDO

! Now, compute density flux vector, store in the corresponging Sij
     DO J=0,NY+1
      DO K=0,NZ+1
       DO I=0,NXP
             Sij(I,K,J,1) = -NU_T(I,K,J)/PR(1)*Sij(I,K,J,1)
          ENDDO
        ENDDO
      ENDDO

      DO K=0,NZ+1
        DO I=0,NXP
           Sij(i,k,NY+1,1) = Sij(i,k,NY,1)
           Sij(i,k,0,1) = Sij(i,k,1,1)
        ENDDO
      ENDDO

      DO J=0,NY+1
       DO I=0,NXP
          Sij(i,NZ+1,j,1) = Sij(i,NZ,j,1)
          Sij(i,0,j,1) = Sij(i,1,j,1)
       ENDDO
      ENDDO
        
! Convert S_ij to Fourier space
      ij=1 

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

! CSij(:,:,:,2 and 3) are added through an implicit eddy viscosity
      DO J=1,NY+1
        DO K=0,NZ+1
          DO I=0,NX2P
             CSij(I,K,J,2) = 0.d0
             CSij(I,K,J,3) = 0.d0
          ENDDO
        ENDDO
      ENDDO

! lamda_i is now contained in CSij in Fourier space

! Convert the scalar to physical space
      CALL REAL_FOURIER_TRANS_TH (.false.)

    ELSEIF ((LES_MODEL_TYPE_TH.EQ.2).or.(LES_MODEL_TYPE_TH.eq.3)) THEN
! based on this paper "An investigation of stably stratiﬁed turbulent channel ﬂow using large-eddy simulation"
! by Armenio and Sarkar, 2002

! Here, use a dynamic smagorinsky model
! Note, there is no scale similar model for the scalar,
! so model type choice 2 and 3 are identical for the scalar equation

! Compute the filter width
      
      DO J=1,NY+1
       DO K=0,NZ+1
! At GY points:
       DELTA_Y(K,J)=beta*(beta*DX(1)*INT_JACOB(K,J))**(1.d0/3.d0)
!       DELTA_Y(J)=sqrt((beta*DX(1))**2.d0+(DY(J)*2.d0)**2.d0
!     &          +(beta*DZ(1))**2.d0)
       ENDDO
      ENDDO

! We need to calculate the components of C, the dynamic coefficient
! C_DYN_TH will be defined at GYF points

! Compute the scalar gradient, store in Sij(:,:,:,1..3)
      CALL compute_scalar_grad(n)

! Convert the scalar to physical space
      CALL REAL_FOURIER_TRANS_TH (.false.)
      
! Compute C_DYN_TH only every x # of timesteps
    IF ( ((MOD(TIME_STEP,10).eq.0).AND.(RK_STEP.eq.1)).OR.FIRST_TIME) THEN

      IF (rank .eq. 0) THEN
	WRITE(6,*) '##############################'
	WRITE(6,*) 'C_dyn_th is recalculating'
	WRITE(6,*) '##############################'
      ENDIF

      CALL allocate_les_tmp

! Store TH in Sij(:,:,:,4) and apply the test filter
      DO j=0,NY+1
        DO k=0,NZ+1
          DO i=0,NXP
            Sij(i,k,j,4)=THX(i,k,j,n)
          ENDDO
        ENDDO
      ENDDO
      S1X = Sij(:,:,:,4)
      CALL FILTER_VAR(2)
      Sij(:,:,:,4)=S1X

!      CALL les_filter_chan(Sij(0,0,0,4),0,NY+1,2)

! Zero C_DYN_TH
      DO j=0,NY+1
       DO k=0,NZ+1
        C_DYN_TH(k,j,N)=0.d0
        denominator_sum(k,j)=0.d0 
        numerator_sum(k,j)  =0.d0
       ENDDO
      ENDDO

! Do over all non-repeating components of the scalar gradient
    DO ij=1,3

! Zero the numerator and denominator:
        DO j=0,NY+1
          DO k=0,NZ+1
            DO i=0,NXP
              numerator(i,k,j)=0.d0
              denominator(i,k,j)=0.d0
            ENDDO
          ENDDO
        ENDDO
        
! First, compute Mij
	DO j=0,NY+1
	  DO k=0,NZ+1
	    DO i=0,NXP
		temp(i,k,j)=Sij(i,k,j,ij)
	    ENDDO
	  ENDDO
	ENDDO
! Filter temp
	S1X = temp
	CALL FILTER_VAR(2)
	temp=S1X

! Multiply by |S| filtered        
	DO j=0,NY+1
	  DO k=0,NZ+1
	    DO i=0,NXP
	      temp(i,k,j) = temp(i,k,j)*(alpha*DELTA_Y(k,j))**2.d0 * S_2BAR(i,k,j) 
	    ENDDO
	  ENDDO
	ENDDO
! Get second term of Mij
	DO i=0,NXP
	  DO k=0,NZ+1 
	    DO j=0,NY+1
	      Mij(i,k,j) = DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)
	    ENDDO
	  ENDDO
	ENDDO
! Filter Mij
	S1X = Mij	
	CALL FILTER_VAR(2)
	Mij=S1X

!       CALL les_filter_chan(Mij,0,NY+1,2)
 
! Add the second term of Mij stored in temp
	Mij=temp-Mij    

  ! Now, compute Lij and add Lij*Mij to the numerator
  ! temp=Ui*Uj:
	SELECT CASE (ij)
	  
	  CASE(1)
	    DO j=0,NY+1
	      DO k=0,NZ+1
		DO i=0,NXP
		  temp(i,k,j) = U1X(i,k,j)*THX(i,k,j,n)
		ENDDO
	      ENDDO 
	    ENDDO
	  
	  CASE(2)
	    DO j=0,NY+1
	      DO k=0,NZ+1
		DO i=0,NXP
		  temp(i,k,j) = THX(i,k,j,n)*U2X(i,k,j)
		ENDDO
	      ENDDO
	    ENDDO
	
	  CASE(3)
	    DO j=0,NY+1
	      DO k=0,NZ+1
		DO i=0,NXP
		  temp(i,k,j) = U3X(i,k,j)*THX(i,k,j,n)
		ENDDO
	      ENDDO 
	    ENDDO
	    
	ENDSELECT
! Filter temp
	S1X = temp
	CALL FILTER_VAR(2)
	temp=S1X

  ! Add Lij*Mij to numerator
  ! ReCALL that Sij(:,:,:,4) holds TH_2BAR
	DO j=0,NY+1
	  DO k=0,NZ+1
	    DO i=0,NXP 
	      numerator(i,k,j) = Mij(i,k,j) * (temp(i,k,j)-Sij(i,k,j,4)*U_BAR_TIL(i,k,j,ij))
	    ENDDO
	  ENDDO
	ENDDO

! Now, the denominator for this ij  
	DO j=0,NY+1
	  DO k=0,NZ+1
	    DO i=0,NXP
	      denominator(i,k,j) = Mij(i,k,j)*Mij(i,k,j)
	    ENDDO
	  ENDDO
	ENDDO

! Get plane average of numerator and denominator, add to C
! Note, since both the numerator and denominator are averaged, there
! is not need to divide the sum by NX*NZ

	DO j=0,NY+1
	 DO k=0,NZ+1
	  denominator_sum(k,j) = denominator_sum(k,j) + SUM(denominator(0:min(NXP,NXP_L),k,j))
	  numerator_sum(k,j) = numerator_sum(k,j)+SUM(numerator(0:min(NXP,NXP_L),k,j))
	 ENDDO
	ENDDO

! End to ij
   ENDDO
      
   CALL MPI_COMBINE_STATS(denominator_sum,NZ+2,NY+2)
   CALL MPI_COMBINE_STATS(numerator_sum,NZ+2,NY+2)

    DO j=jstart,jEND
      DO k=zstart,zEND
	  IF (denominator_sum(k,j) .NE. 0.) THEN
	    C_DYN_TH(k,j,n) = -0.5d0 * numerator_sum(k,j)/denominator_sum(k,j)
	  ELSE
	    C_DYN_TH(k,j,n) = 0.d0
	  ENDIF
      ENDDO
    ENDDO

! We are now Done with the dynamic procedure to calculate C

! If C_DYN_TH < 0 at any level, set C_DYN_TH=0 for numerical stability
      DO j=0,NY+1
       DO k=0,NZ+1
        IF (C_DYN_TH(k,j,n).lt.0) C_DYN_TH(k,j,n)=0.d0
       ENDDO
      ENDDO

      CALL MPI_BCAST_REAL(C_DYN_TH(0,0,N), NZ+2, NY+2)

! End if compute C_DYN_TH
       CALL deallocate_les_tmp

    ENDIF

! Get the eddy diffusivity at GY points
! KAPPA_T = C_DYN_TH * (DELTA^2)*|S|
! Use exact second order interpolation to get C_DYN_TH and S1 interpolated to GY points

      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            KAPPA_T(I,K,J,N) = C_DYN_TH(k,j,n)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)
          ENDDO
        ENDDO
      ENDDO

! At this point we have C_DYN_TH and dTH/dx_i (stored in Sij(:,:,:,1...3)
! Calculate lambda_i in physical space, stored in Sij(:,:,:,1..3)

    DO ij=1,3
! Dynamic Smagorinsky model, no scale similar part
      DO j=0,NY+1
	DO k=0,NZ+1
	  DO i=0,NXP
	    Sij(i,k,j,ij) = - C_DYN_TH(k,j,n)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)
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
	

! End DO ij
      ENDDO  

! End if LES_MODEL_TYPE dynamic Smagorinsky or Dynamic model
    ELSE
     PAUSE 'Error, unsupported LES_MODEL_TYPE chosen'
    ENDIF

! Now, add the subgrid scale forcing to CFi
! (This includes the subgrid scale stress as an explicit R-K term
! Include only CSij terms 1 and 3 since term 2 is accounted for
! as an implicit eddy diffusivity through KAPPA_T

!      DO J=J2+1,NY+1
!        DO K=0,TNKZ
!          DO I=0,NKX
!            DO ij =1,6
!             CSij(I,K,J,ij) = 0.d0
!            ENDDO
!           ENDDO
!         ENDDO
!       ENDDO

      DO J=JSTART_TH(N),JEND_TH(N) 
       DO K=ZSTART_TH(N),ZEND_TH(N)
        DO I=0,NX2P
!The other two derivatives of density flux terms will be added implicitly through an eddy viscosity (duct.F90 line:~1600)
          CFTHX(I,K,J,N) = CFTHX(I,K,J,N) - INT_JACOB(K,J)*CIKXP(I)*CSij(I,K,J,1) 
        ENDDO
       ENDDO
      ENDDO

! PeriodiCALLy, output mean quantities
     IF ((MOD(TIME_STEP,SAVE_STATS_INT).EQ.0).AND.(RK_STEP.EQ.1)) THEN
	
       DO j=0,NY+1
        DO k=0,NZ+1
         KAPPA_T_MEAN(k,j) = SUM(KAPPA_T(0:min(NXP,NXP_L),k,j,n))/dble(NX)
        ENDDO
       ENDDO

       CALL MPI_COMBINE_STATS(KAPPA_T_MEAN,NZ+2,NY+2)
		
	IF (rank .eq. 0) THEN
          !check for directory existence
          INQUIRE(DIRECTORY='./plane_les_th/.', EXIST=dir_exists)

          IF ( .NOT. dir_exists) THEN
            WRITE(6,*) 'plane_les directory DOes not exists'
            CALL system('mkdir plane_les_th')
            WRITE(6,*) 'plane_les_th is created!'
          ENDIF

	  k = time_step/SAVE_STATS_INT

	  file_name = 'plane_les_th/les_tec_'  &
		//CHAR(MOD(k,100000)/10000+48) &
		//CHAR(MOD(k,10000)/1000+48)   &
		//CHAR(MOD(k,1000)/100+48)     &
		//CHAR(MOD(k,100)/10+48)       &
		//CHAR(MOD(k,10)+48) //        &
		'.plt'

	  CALL plot_les_th_tecplot(file_name)     
	ENDIF
      ENDIF
     

  RETURN
END
!***********************************
!***********************************
SUBROUTINE compute_scalar_grad(n)
!C This subroutine computes dTH/dx_i for the filtered scalar field
!C The input velocity field should be in fourier space in the periodic
!C directions.
!C For use in the LES model in channel flow (2 periodic directions)
!C Store in Sij(:,:,:,1..3)
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
  USE les_chan_var

    IMPLICIT NONE

    INTEGER I,J,K,ij,n

    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
	DO I=0,NX2P
	  CSij(I,K,J,1)=CIKXP(I)*CTHX(I,K,J,n)
	  CSij(I,K,J,2)=( 0.5*CJOB_12(K+1,J,2)*(CTHX(I,K,J,n) + CTHX(I,K+1,J,n))   &
		      - 0.5*CJOB_12(K,J,2)*(CTHX(I,K,J,n)   + CTHX(I,K-1,J,n))         &
		      + 0.5*CJOB_22(K,J+1,1)*(CTHX(I,K,J,n) + CTHX(I,K,J+1,n))         &
		      - 0.5*CJOB_22(K,J,1)*(CTHX(I,K,J,n)   + CTHX(I,K,J-1,n)) ) /     &
			    INT_JACOB(K,J)
	  CSij(I,K,J,3)= ( 0.5*CJOB_11(K+1,J,2)*(CTHX(I,K,J,n) + CTHX(I,K+1,J,n)) &
		      - 0.5*CJOB_11(K,J,2)*(CTHX(I,K,J,n)   + CTHX(I,K-1,J,n))    &
		      + 0.5*CJOB_21(K,J+1,1)*(CTHX(I,K,J,n) + CTHX(I,K,J+1,n))    &
		      - 0.5*CJOB_21(K,J,1)*(CTHX(I,K,J,n)   + CTHX(I,K,J-1,n))  )/&
			    INT_JACOB(K,J)
	ENDDO
      ENDDO
    ENDDO

! Convert the scalar gradients to physical space
    DO ij=1,3
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

    ENDDO

! adding the background-th gradient
    DO J=JSTART,JEND
      DO K=ZSTART,ZEND
        DO I=0,NXP
          Sij(I,K,J,2)= Sij(I,K,J,2) + (0.5 * CJOB_12(K+1,J,2) * (TH_BAR(K,J) + TH_BAR(K+1,J) ) &
                    - 0.5 * CJOB_12(K,J,2) * ( TH_BAR(K,J) + TH_BAR(K-1,J) )   &
                    + 0.5 * CJOB_22(K,J+1,1) * ( TH_BAR(K,J) + TH_BAR(K,J+1) ) &
                    - 0.5 * CJOB_22(K,J,1) * ( TH_BAR(K,J) + TH_BAR(K,J-1) )   &
                        )/INT_JACOB(K,J)
        ENDDO
      ENDDO
    ENDDO

! We now have dTH/dx_i in Physical space
  RETURN
END


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C *** Subroutine to write a .plt file (uses TecPlot tecio.a library)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE plot_les_th_tecplot(file_name_plt)

  USE ntypes
  USE Domain
  USE Grid
  USE run_variable
  USE les_chan_var

    IMPLICIT NONE

    INTEGER  I,J,K,N
    INTEGER  ncycles
    INTEGER  imax,jmax,kmax
    INTEGER  debug,ier,itot
    INTEGER  tecini,tecdat,teczne,tecEND
    INTEGER  visDOuble,disDOuble
    REAL*8   phase,PI,ncycle
    CHARACTER*1 nulchar
    CHARACTER(len=33) :: title
    CHARACTER*30 file_name_plt

    REAL(r8) THW_SGS(0:NZ+1,0:NY+1,1:N_TH), DTHDY
    
    PI  = ATAN(1.0)*4.0
    nulchar = char(0)
    debug   = 0
    visDOuble = 0
    disDOuble = 1
    imax = NZ+2
    jmax = NY+2
    kmax = 1

! computing SGS density flux as: thw_sgs = - Kappa * dthdy
        DO N=1,N_TH
          DO K=1,NZ
            DO J=1,NY
                dthdy = ( 0.125 * (CJOB_12(K,J,1) + CJOB_12(K,J+1,1)+CJOB_12(K,J,2) + CJOB_12(K+1,J,2) ) &
                                * dble ( CRTHX(0,k+1,j,n) - CRTHX(0,k-1,j,n))  &
                                + 0.125 * (CJOB_22(K,J,1) + CJOB_22(K,J+1,1)+CJOB_22(K,J,2) + CJOB_22(K+1,J,2) ) &
                                * dble ( CRTHX(0,k,j+1,n) - CRTHX(0,k,j-1,n)) )/INT_JACOB(K,J)
! adding the background dth/dy
                IF (.NOT. CONT_STRAT) THEN
                  dthdy =  dthdy +                          &
                                  ( 0.125 * (CJOB_12(K,J,1) + CJOB_12(K,J+1,1)+CJOB_12(K,J,2) + CJOB_12(K+1,J,2) )   &
                                  * ( TH_BAR(k+1,j)-TH_BAR(k-1,j) )&
                                  + 0.125 * (CJOB_22(K,J,1) + CJOB_22(K,J+1,1)+CJOB_22(K,J,2) + CJOB_22(K+1,J,2) )   &
                                  * ( TH_BAR(k,j+1)-TH_BAR(k,j-1) ) ) /INT_JACOB(K,J)

                  THW_SGS(K,J,N) = -(KAPPA_T_MEAN(K,J)+NU/PR(N)) * dthdy
                ENDIF
            ENDDO
          ENDDO
        ENDDO

    WRITE(6,*) file_name_plt 

    ncycles = time/(2.d0*pi)
    phase = (time-ncycles*2.d0*pi)*180.d0/pi
      
    WRITE(title,'(a5,f8.4,a7,f5.1)') 'time=',time,' phase=', phase
      
     ! Write the zone header information.
!      ier = tecini(trim(title)//nulchar,'x,z,ume,wme,vme,C_dyn,nu_t_mean,& 
!              &tke_sgs_p,tke_sgs_diss'&
     ier =tecini(trim(title)//nulchar,'x,z,nu_t_mean,Kappa_T,nu_t_mean_PR,thw_sgs '&
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
!      ier = tecdat(itot,C_DYN(0:NZ+1,0:NY+1),disDOuble)
     ier = tecdat(itot,NU_T_mean(0:NZ+1,0:NY+1),disDOuble)
     ier = tecdat(itot,KAPPA_T_MEAN(0:NZ+1,0:NY+1),disDOuble)
     ier = tecdat(itot,1/PR(1)*NU_T_mean(0:NZ+1,0:NY+1),disDOuble)
     ier = tecdat(itot,THW_SGS(0:NZ+1,0:NY+1,1),disDOuble)

     ! close file
     ier = tecEND()
 
  RETURN
END
