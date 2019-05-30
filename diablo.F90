! C***********************************************************************************************|
! C  diablo.f90 --> Stably Stratified Duct Flow DNS/LES    Last Update: Oct 2013  VERSION 1.2
! C 
! C  This Fortran 90 code computes Stably Stratified flow in a Duct. Code is 
! C  in curvilinear coordinate and has MPI parallelization. 
! C
! C  Primative variables (u,v,w,p) are used, and continuity is enforced with a
! C  fractional step algorithm.
! C 
! C  There is governing equation for boussinesq approximation for stably stratified flow.
! C
! C  SPATIAL DERIVATIVES:
! C    1 primary directions is taken to be periodic and handled spectrally 
! C    (this case is commonly referred to as "duct")
! C    The remaining directions are taken to be bounded by walls and handled with
! C    momentum- and energy-conserving second-order central finite differences.
! C
! C  Coordinates:
! C    Code is written in curvilinear colocated coordinate, but Stream-wise Direction is referred as "Z",
! C    Vertical Direction as "Y", and Span-wise as "X" (Span-wise Direction is handled spectrally) 
! C
! C  TIME ADVANCEMENT
! C    The main approach is implemented as:
! C     RKW3 on non-linear terms and CN on viscous terms over each RK sub-step.
! C
! C  A few simple high-performance programming constructs are used:
! C   -> The inner 2 loops are broken out in such a way as to enable out-of-order
! C      execution of these loops as much as possible, thereby leveraging
! C      vector and superscalar CPU architectures.
! C   -> The outer loops are fairly long (including as many operations as
! C      possible inside on a single J plane of data) in order to make effective
! C      use of cache.
! C   -> There are some careful thougths on saving flops in few loops regarding 
! C      divisons/multipications
! C
! C  Multiple time advancement algorithms are implemented in order to compare 
! C  its efficiency for various configuration on different computational architectures.
! C
! C  This is an in-house code modified and advanced in CFD Lab in UC San Digeo. 
! C  Primary contributions follow:
! C   Thomas Bewley		(2001-2004)
! C   John Taylor               (2004-2007)
! C   Bishakhdatta Gayen 	(2007-2012)
! C   Masoud Jalalli 		(2011-Present)		mjalalib@ucsd.edu
! C   S. M. Iman Gohari		(2013-Present)		sgohari@ucsd.edu
!
! C  ACKNOWLEDGEMENT: The code was originally developed as a joint project in MAE 223 (CFD),
! C  taught by Thomas Bewley, at UC San Diego (spring 2001, spring 2005).
! C***********************************************************************************************|

! C----|--.---------.---------.---------.---------.---------.---------.-|-------|
! 
! C----|--.---------.---------.---------.---------.---------.---------.-|-------|
PROGRAM DIABLO
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
! c USE omp_lib
  USE mpi_var

    IMPLICIT NONE

    INTEGER N, K
    REAL*8  TRANS_TIME, FINAL_TIME
    

    IF(rank.EQ.0) THEN
      WRITE(6,*) 
      WRITE(6,*) '             ************ WELCOME TO DIABLO ************'
      WRITE(6,*)
      WRITE(6,*) 'Initializing...'
    ENDIF

    CALL INITIALIZE

! c Initialize START_TIME for run timing

    IF(rank.EQ.0) WRITE(6,*) 'Initialization Done!'
		

    CALL DATE_AND_TIME (VALUES = TIME_ARRAY)
    START_TIME = TIME_ARRAY(5) * 3600 + TIME_ARRAY(6) * 60 + TIME_ARRAY(7) + TIME_ARRAY(8)*0.001

! c A flag to determine if we are considering the first time-step
    FIRST_TIME=.TRUE. 
    
!    INT_TREAT = .TRUE. 
    INT_TIME = TIME
!    TRANS_TIME = 900.d0
    TRANS_TIME = 0.50d0/F_0 
    FINAL_TIME = 6.00d0/F_0
   
 ! This is specific calculations just for time dependent BC for TH_LOW   
    INT_THX = 0.d0
 
    IF (Rank.eq.0) THEN
      DO N=1,N_TH 
	DO K=ZSTART,ZEND
            INT_THX = INT_THX + CTHX(0,K,JSTART_TH(N)-1,N)/ dble(ZEND-ZSTART+1)
	ENDDO
      ENDDO
    ENDIF
    CALL MPI_BCAST(INT_THX,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, IERROR)
 ! TEND the time dependent BC for TH_LOW  
    
    DO TIME_STEP = TIME_STEP+1, TIME_STEP+N_TIME_STEPS

	IF (rank .eq. 0) THEN
	  WRITE(6,*) ' '
	  WRITE(6,'(i6,a9,f15.8,a7,f15.8)') TIME_STEP,"    time=",TIME,"    dt=",DELTA_T

	  open(99,file='out_screen.txt',form='formatted', &
			status='unknown', position='append')
	  WRITE(99,'(i6,a9,f15.8,a7,f15.8)') TIME_STEP,"    time=",TIME,"    dt=",DELTA_T
	  CLOSE(99)
	ENDIF

! c	Time marching loop begins here

! Iman: I have added BC=5, for finite time cooling flux imposed
        DO N=1,N_TH
         IF ( TH_BC_YMIN(N).EQ.5 .AND. TIME.GT.FINAL_TIME) THEN
           IF ( RANK.EQ.0) WRITE (6,*) '**WARNING: Finite time cooling flux is imposed!**'
           CALL SAVE_STATS_CURVI(.FALSE.)
           CALL SAVE_FLOW(.FALSE.)
           CALL MPI_FINALIZE(IERROR)
         ENDIF
        ENDDO

!Jose
!           probably a minor bug in the computation of the transient value of the Rig when restarting. 
!           Have to reset time when computing stratified case.\\
	DO N=1,N_TH 
	  IF ( INT_TREAT ) THEN
	    IF ( TIME .LT. TRANS_TIME) THEN
!	    IF ( (TIME-INT_TIME) .LT. TRANS_TIME) THEN
              IF ( RANK.EQ.0) WRITE (6,*) '**WARNING: initial TREATMENT for stratification has been used**'
!	      RI_TAU(N) = RI_FINAL * ( TIME - INT_TIME ) / TRANS_TIME
	      RI_TAU(N) = RI_FINAL * ( TIME ) / TRANS_TIME
	    ELSE
	      RI_TAU(N) = RI_FINAL
	    ENDIF    
	  ELSE
	    RI_TAU(N) = RI_FINAL
	  ENDIF   
	ENDDO   
        
   
	DO RK_STEP=1,3
	    CALL RK_DUCT
	ENDDO
        
	TIME = TIME + DELTA_T
	FIRST_TIME = .FALSE.
  

! c Save statistics to an output file
	IF (MOD(TIME_STEP,SAVE_STATS_INT).EQ.0) CALL SAVE_STATS_CURVI(.FALSE.)

! c Save the flow to a restart file
	IF (MOD(TIME_STEP,SAVE_FLOW_INT).EQ.0) CALL SAVE_FLOW(.FALSE.)

	IF ( (MOD(TIME_STEP,BACKUP_FLOW_INT).EQ.0) .AND. (TIME_STEP .GT. 1) ) CALL BACKUP_FLOW
        
    ENDDO

! Calculate and display the runtime for the simulation
    CALL DATE_AND_TIME (VALUES=TIME_ARRAY)
    END_TIME = TIME_ARRAY(5)*3600+TIME_ARRAY(6)*60 &
		+TIME_ARRAY(7)+TIME_ARRAY(8)*0.001
    WRITE(*,*) 'Elapsed Time (sec): ',end_time-start_time
    WRITE(*,*) 'Seconds per Iteration: ' &
		,(end_time-start_time)/N_TIME_STEPS

    TIME_STEP=TIME_STEP-1
    
    CALL SAVE_FLOW(.TRUE.)
    CALL SAVE_STATS_CURVI(.TRUE.)
    
    WRITE(6,*)
    WRITE(6,*) '        ****** FINISHED! Have a nice day! ******'
    WRITE(6,*)
   
    CALL MPI_FINALIZE(IERROR)

END


! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE INITIALIZE
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
  USE mg_vari, 		ONLY : INIT_FLAG, MUDPACK_BC
  USE les_chan_var , 	ONLY : C_DYN     
  USE mpi_var 

      IMPLICIT NONE

      REAL    VERSION, CURRENT_VERSION
      LOGICAL RESET_TIME
      INTEGER I, J, K, N
     
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

! C Initialize the random number generator
      
      CALL RANDOM_SEED(SIZE = K)
      ALLOCATE (seed(1:K))
      seed(1:K)=10
      CALL RANDOM_SEED(PUT = seed)

      IF (rank .eq. 0) THEN
      
	OPEN (99, file='out_screen.txt', form='formatted', &
		      status='unknown', position='append') 
  
	OPEN (11, file='input.dat', form='formatted', status='old')      

	WRITE(6,*)  'Note that this code is advanced in CFD LAB at UC San Diego'
	WRITE(6,*)
  
	WRITE(99,*) 'Note that this code is advanced in CFD LAB at UC San Diego'
	WRITE(99,*)

      ENDIF 

! C ALLOCATE VARIABLES

      CALL   allocate_var
     
! C    Read input file.
! C   (Note - if you change the following section of code, update the
! C    CURRENT_VERSION number to make obsolete previous input files!)

      CURRENT_VERSION = 1.21
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*) VERSION
      IF (VERSION .NE. CURRENT_VERSION) STOP 'Wrong input file Version! Please check Diablo.F90 and input file '
      
      READ(11,*)
      READ(11,*) NU, LX, LY, LZ
      READ(11,*)
      READ(11,*) N_TIME_STEPS, DELTA_T_set, RESET_TIME, VARIABLE_DT, CFL, DFN, UPDATE_DT
      READ(11,*)
      READ(11,*) LES, LES_MODEL_TYPE, LES_MODEL_TYPE_TH, WALL_MODEL_TYPE

      IF (LES .AND. LES_MODEL_TYPE.EQ.1 .AND. LES_MODEL_TYPE_TH.EQ.2)  THEN
        STOP '***WORNG LES CONFIG*** (LES_MODEL_TYPE=1 & LES_MODEL_TYPE_TH=2) would not work correctly!'
      ENDIF

      READ(11,*)
      READ(11,*) MUDPACK_USE, MUDPACK_BC(1), MUDPACK_BC(2), MUDPACK_BC(3), MUDPACK_BC(4)
      
      IF (MUDPACK_USE)  THEN
	WRITE(6,*), 'MUDPACK Will be Used as Remove_Divergence Package. '
	
	WRITE(*,801) MUDPACK_BC(1), MUDPACK_BC(2), MUDPACK_BC(3), MUDPACK_BC(4)
	801   FORMAT(/' Bc  diablo nxa = ',i2,' nxb = ',i2,' nyc = ',i2,' nyd = ',i2)
      ELSE
	WRITE(6,*), 'DMGD9V Will be Used as Remove_Divergence Package. '
      ENDIF

      READ(11,*)
      READ(11,*) VERBOSITY, SAVE_FLOW_INT, SAVE_STATS_INT, BACKUP_FLOW_INT, MOVIE
      READ(11,*)
      READ(11,*) CREATE_NEW_FLOW, IC_TYPE, KICK, INT_TREAT, WRITE_VEL_TECPLOT
      READ(11,*)
      READ(11,*) F_TYPE, WBULK0, PX0, OMEGA0, AMP_OMEGA0, ANG_BETA, F_0, SPONGE_TYPE
 
!  C Reading the Bottom and TOP BCs  (Vertical)
      READ(11,*)
      READ(11,*) U_BC_YMIN, U_BC_YMIN_C1, U_BC_YMIN_C2, U_BC_YMIN_C3
      READ(11,*) 
      READ(11,*) V_BC_YMIN, V_BC_YMIN_C1, V_BC_YMIN_C2, V_BC_YMIN_C3
      READ(11,*)
      READ(11,*) W_BC_YMIN, W_BC_YMIN_C1, W_BC_YMIN_C2, W_BC_YMIN_C3
      READ(11,*)
      READ(11,*) U_BC_YMAX, U_BC_YMAX_C1, U_BC_YMAX_C2, U_BC_YMAX_C3
      READ(11,*)
      READ(11,*) V_BC_YMAX, V_BC_YMAX_C1, V_BC_YMAX_C2, V_BC_YMAX_C3
      READ(11,*)
      READ(11,*) W_BC_YMAX, W_BC_YMAX_C1, W_BC_YMAX_C2, W_BC_YMAX_C3

 !  C Reading the LEFT and RIGHT BCs  (Stream-Wise) 
      READ(11,*)
      READ(11,*) U_BC_ZMIN, U_BC_ZMIN_C1, U_BC_ZMIN_C2, U_BC_ZMIN_C3
      READ(11,*) 
      READ(11,*) V_BC_ZMIN, V_BC_ZMIN_C1, V_BC_ZMIN_C2, V_BC_ZMIN_C3
      READ(11,*)
      READ(11,*) W_BC_ZMIN, W_BC_ZMIN_C1, W_BC_ZMIN_C2, W_BC_ZMIN_C3
      READ(11,*)
      READ(11,*) U_BC_ZMAX, U_BC_ZMAX_C1, U_BC_ZMAX_C2, U_BC_ZMAX_C3
      READ(11,*)
      READ(11,*) V_BC_ZMAX, V_BC_ZMAX_C1, V_BC_ZMAX_C2, V_BC_ZMAX_C3
      READ(11,*)
      READ(11,*) W_BC_ZMAX, W_BC_ZMAX_C1, W_BC_ZMAX_C2, W_BC_ZMAX_C3

      
! Read Stochastic forcing parameters
      READ(11,*)
      READ(11,*) STOCHASTIC_FORCING
      
! Filter for the velocity field
      READ(11,*)
      READ(11,*) FILTER_VEL, FILTER_VEL_INT
      
      READ(11,*)

! Read the input parameters for the N_TH scalars
      DO N=1,N_TH
        READ(11,*)
        READ(11,*) CREATE_NEW_TH(N)
        READ(11,*)
        READ(11,*) FILTER_TH(N), FILTER_INT(N)
        READ(11,*)
        READ(11,*) RI_FINAL, PR(N), BACKGROUND_GRAD(N) , CONT_STRAT
        READ(11,*)
        READ(11,*) TH_BC_YMIN(N),TH_BC_YMIN_C1(N),TH_BC_YMIN_C2(N),TH_BC_YMIN_C3(N)
        READ(11,*)
        READ(11,*) TH_BC_YMAX(N),TH_BC_YMAX_C1(N),TH_BC_YMAX_C2(N),TH_BC_YMAX_C3(N)
        

	READ(11,*)
	READ(11,*) TH_BC_ZMIN(N),TH_BC_ZMIN_C1(N),TH_BC_ZMIN_C2(N),TH_BC_ZMIN_C3(N)
	READ(11,*)
	READ(11,*) TH_BC_ZMAX(N),TH_BC_ZMAX_C1(N),TH_BC_ZMAX_C2(N),TH_BC_ZMAX_C3(N)
      ENDDO
	
      DELTA_T=DELTA_T_set
      

! C If we are using MPI, then Initialize the MPI Variables

      CALL INT_MPI

      
      AMP_OMEGA0=AMP_OMEGA0 * OMEGA0

      DO N=1,N_TH
       RI_TAU(N) = RI_FINAL
      ENDDO
  
      ANG_BETA = ANG_BETA * 3.14159265d0 / 180.d0 
      Q_H0  = 1.d0
      H0    = 10.0*LY     
      In_H0 = 1.d0/H0 
      count_data = 0
      
      Q0=WBULK0
      Qn=Q0
      weight=0.1


 ! C Initialize grid
      IF (rank .eq. 0) THEN

	WRITE(6,*) 'Grid size: NX =',NX,', NY =',NY,', NZ =',NZ,'.'
	WRITE(6,*) 'Domain size: LX =',LX,', LY =',LY,', LZ =',LZ,'.' 
	WRITE(6,*) 'NU',NU,'DELTA_T',DELTA_T,'Delta_s',sqrt(2.0*NU/OMEGA0) 
	WRITE(6,*) 'CFL',CFL, 'Diffu No' , DFN
	WRITE(6,*) 'Variable delta_t applited',  VARIABLE_DT,U_BC_ZMAX_C1

	WRITE(99,*) 'Grid size: NX =',NX,', NY =',NY,', NZ =',NZ,'.'
	WRITE(99,*) 'Actual Grid size: NXV =',NXV,', NY =',NY,', NZV =',NZV,'.'
	WRITE(99,*) 'Domain size: LX =',LX,', LY =',LY,', LZ =',LZ,'.'
	WRITE(99,*) 'NU', NU,'DELTA_T',DELTA_T,'Delta_s', sqrt(2.0*NU/OMEGA0)
	WRITE(99,*) 'Amp of tide', AMP_OMEGA0,'Frequency of tide', OMEGA0, 'Frequency of Rot', F_0
	WRITE(99,*) 'CFL', CFL, 'Diffu No', DFN
	WRITE(99,*) 'Variable delta_t applited',  VARIABLE_DT,U_BC_ZMAX_C1 

	IF (LES) THEN
	  WRITE(99,*) 'Model Type for LES Vel is', LES_MODEL_TYPE, 'and for Temperature is ', &
	  LES_MODEL_TYPE_TH  
	ENDIF

	DO N=1,N_TH
	  WRITE(6,*) 'Scalar number: ', N
	  WRITE(6,*) 'Final Buoyancy number: ', RI_FINAL
	  WRITE(6,*) 'Prandlt number: ', PR(N)
	  
	  WRITE(99,*) 'Scalar number: ', N
	  WRITE(99,*) 'Final buoyancy number: ', RI_FINAL
	  WRITE(99,*) 'Prandlt number: ', PR(N)
	ENDDO
	
	close(99)
      ENDIF
   
      NXM = NX-1
      NYM = NY-1
      NZM = NZ-1
      NXP_L = NX-(NXP+1)*rank-1
      NX2P_L = NKX+1-NX2P*rank-1
      write(6,'(a,i6,a,i6,a,i6)')'rank:',rank,', index bound:', min(NXP_L,NXP), &
      								'  index bound 2:',min(NX2P_L,NX2P)

      DO I=0,NXP
	GX(I) = (I*LX)/NX
	DX(I) = LX/NX
	IF (VERBOSITY .GT. 3) WRITE(6,*) 'GX(',I,') = ',GX(I)
      ENDDO
 
!          WRITE (6,*) 'Finite-difference in X'
!          OPEN (30,file='xgrid.txt',form='formatted',status='old')
!          READ (30,*) NX_T
! ! C Check to make sure that grid file is the correct dimensions
!          IF (NX_T.ne.NX) THEN
!            WRITE(6,*) 'NX, NX_T',NX,NX_T
!            STOP 'Error: xgrid.txt wrong dimensions'
!          END IF
!          DO I=1,NX+1
!            READ(30,*) GX(I)
!            IF (VERBOSITY .GT. 3) WRITE(6,*) 'GX(',I,') = ',GX(I)
!          END DO 
!          DO I=1,NX
!            READ(30,*) GXF(I)
!            IF (VERBOSITY .GT. 3) WRITE(6,*) 'GXF(',I,') = ',GXF(I)
!          END DO
! ! C Define ghost cells, if needed for this grid...
!          GXF(0)=2.d0*GXF(1)-GXF(2)
!          GXF(NX+1)=2.d0*GXF(NX)-GXF(NXM)
!          GX(0)=2.d0*GX(1)-GX(2)
! ! C Define the grid spacings 
!          DO I=1,NX+1
!            DX(I)=(GXF(I)-GXF(I-1))
!          END DO
!          DO I=1,NX
!            DXF(I)=(GX(I+1)-GX(I))
!          END DO
!          CLOSE(30)


! C   Initial Condition prescribed by the Zero's
 
      DO J=0,NY+1
       DO K=0,NZ+1
        U2_BAR(K,J) = 0.d0
        U1_BAR(J) = 0.d0
        U3_BAR(K,J) = 0.d0
       ENDDO 
      ENDDO
      
! C   Boundary condition prescribed by the Blasius for IC_TYPE = 5
! c      IF ( IC_TYPE .EQ. 5) THEN
! c      OPEN (80,file='./int_vel_xy.txt',form='formatted',status='old')     
! c       DO J=1,NY+2
! c        DO K=1,NZ+2
! c        read(80,*) U3_BAR(K,J), U2_BAR(K,J)
! c       ENDDO
! c       ENDDO
! 
! c       DO J=1,NY+1 
! c        DO K=2,NZ+1
! c          S1(0,K,J) = 0.5*(U3_BAR(K,J) + U3_BAR(K-1,J))
! c        ENDDO
! c       ENDDO
! 
! c       DO J=1,NY
! c        DO K=2,NZ
! c          U3_BAR(K,J) =S1(0,K,J)
! c        ENDDO
! c       ENDDO 
!        
! c       DO K=1,NZ
! c          U3_BAR(K,NY+1) =2.0*U3_BAR(K,NY)-U3_BAR(K,NY-1)
! c	  U2_BAR(K,NY+1) =2.0*U2_BAR(K,NY)-U2_BAR(K,NY-1)
! c       ENDDO
       
      DO J=0,NY+1
	DO K=0,NY+1
	    U3_BAR(K,J) = 0.0d0
	    U2_BAR(K,J) = 0.0d0
	ENDDO   
      ENDDO

!       CALL VEL_PROFILE              
      
      CALL allocation_ub (.true.)

      IF (LES) THEN
	CALL allocation_les_var
      ENDIF

      IF ( IC_TYPE .EQ. 7 ) THEN 
	CALL	JACOBI_TRANS
	CALL 	VEL_PROFILE
      ENDIF 


     CALL allocation_u (.true.)
     CALL allocation_v (.true.)
     CALL allocation_w (.true.)  
     CALL allocation_p (.true.)

     CALL allocation_R1 (.true.) 
     CALL allocation_R2 (.true.)
     CALL allocation_R3 (.true.) 
     CALL allocation_F1 (.true.)
     CALL allocation_F2 (.true.)
     CALL allocation_F3 (.true.)
 


      DO N=1,N_TH
       CALL allocation_th (.true.)
       CALL allocation_Rth (.true.)
       CALL allocation_Fth (.true.)
      ENDDO

! C Initialize storage arrays.
      DO K=0,NZV-1
        DO I=0,NXP
          DO J=0,NY+1
            U1X(I,K,J)=0.d0
            U3X(I,K,J)=0.d0
            U2X(I,K,J)=0.d0
            PX (I,K,J)=0.d0
            R1X(I,K,J)=0.d0
            R2X(I,K,J)=0.d0
            R3X(I,K,J)=0.d0
            F1X(I,K,J)=0.d0
            F2X(I,K,J)=0.d0
            F3X(I,K,J)=0.d0
            DO N=1,N_TH
             THX(I,K,J,N)=0.d0
             RTHX(I,K,J,N)=0.d0 
             FTHX(I,K,J,N)=0.d0
            END DO
            
! Array for LES subgrid model
            IF (LES) THEN
             NU_T(I,K,J)= 0.d0 !0.5 * NU  
             DO N=1,N_TH
               KAPPA_T(I,K,J,N)=0.d0
             ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO


      IF (LES) THEN
       DO K=0,NZ+1
        DO J=0,NY+1
          C_DYN(K,J)=0.d0
        ENDDO
       ENDDO
      ENDIF

! C Initialize FFT package (includes defining the wavenumber vectors).

!$OMP PARALLEL
      CALL INIT_FFT
!$OMP END PARALLEL

 
! C Initialize RKW3 parameters.
      H_BAR(1) = DELTA_T * (8.0d0/15.0d0)
      H_BAR(2) = DELTA_T * (2.0d0/15.0d0)
      H_BAR(3) = DELTA_T * (5.0d0/15.0d0)
      BETA_BAR(1) = 1.d0
      BETA_BAR(2) = 25.0d0/8.0d0
      BETA_BAR(3) = 9.0d0/4.0d0
      ZETA_BAR(1) = 0.d0
      ZETA_BAR(2) = -17.0d0/8.0d0
      ZETA_BAR(3) = -5.0d0/4.0d0
      
! C Initialize case-specific packages.
      CALL INIT_DUCT


!ccccccccccccccccccccccccccccccccccccccccccccccccc
! call open_MP_initialization
!ccccccccccccccccccccccccccccccccccccccccccccccccc


! C Initialize values for reading of scalars
      NUM_READ_TH=0
      
      DO N=1,N_TH
	IF (CREATE_NEW_TH(N)) THEN
	  NUM_READ_TH = NUM_READ_TH 
	ELSE
	  NUM_READ_TH=NUM_READ_TH + 1
	  READ_TH_INDEX(NUM_READ_TH)=N
	ENDIF

	CALL CREATE_TH_DUCT
      
      ENDDO

      
      IF (CREATE_NEW_FLOW) THEN
	
	CALL CREATE_FLOW_DUCT

	
	IF(rank.EQ.0) write(*,*) 'A new flowfield has been created'

	CALL SAVE_FLOW(.FALSE.)

      ELSE
	
  ! C 	Convert velocity back to Fourier space
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
      
	CU1X(:,:,:) = (0.d0,0.d0)
	CU2X(:,:,:) = (0.d0,0.d0)
	CU3X(:,:,:) = (0.d0,0.d0)
	CPX(:,:,:)  = (0.d0,0.d0)

	DO N=1,N_TH
	  IF (.NOT.(CREATE_NEW_TH(N))) THEN
	    CALL REAL_FOURIER_TRANS_TH (.true.)      
	    CALL REAL_FOURIER_TRANS_Rth (.true.)
	    CALL allocation_Fth(.false.)      

	    CTHX(:,:,:,N) = (0.d0,0.d0)
	  ENDIF
	ENDDO
	
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IF(rank.EQ.0) write(*,*) 'Reading flow...'
	  CALL READ_FLOW
	IF(rank.EQ.0) write(*,*) 'Done reading flow'

	IF ( INT_TREAT ) THEN
    ! TEMPORARY!! SCALE THE VELOCITY FLUCTIONS
	ENDIF

    ! C Initialize flow.
	IF (RESET_TIME .OR. CREATE_NEW_FLOW) THEN
	  PREVIOUS_TIME_STEP = 0
	  TIME_STEP = 0
	  TIME = 0
	ENDIF

	CALL SAVE_STATS_CURVI(.FALSE.)
	INIT_FLAG = .TRUE.       
	
	CALL SAVE_STATS_CURVI(.FALSE.)   

    ENDIF

  RETURN
END


! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE SAVE_FLOW2(FINAL)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
  USE mpi_var, 		ONLY : rank      

    IMPLICIT NONE

    CHARACTER*35 FNAME
    CHARACTER*35 FNAME_TH(N_TH)
    INTEGER      I, J, K, N, no_p,ny_p
    LOGICAL      FINAL, dir_exists
    PARAMETER  (no_p=4, ny_p=66)
      

    INQUIRE(DIRECTORY='./last_saved2/.', EXIST=dir_exists) 

    IF ( .NOT. dir_exists) THEN
      WRITE(6,*) 'last_saved2 directory does not exists'
      CALL system('mkdir last_saved2')
      WRITE(6,*) 'last_saved2 is created!'
    ENDIF    
      
      

!  write the latest file for restart
    I=1+RANK
    IF (FINAL) THEN
      FNAME='last_saved2/diablo_'      &
	  //CHAR(MOD(I,1000)/100+48)  & 
	  //CHAR(MOD(I,100)/10+48)    &
	  //CHAR(MOD(I,10)+48) // '.res'


      DO N=1,N_TH
	FNAME_TH(N)='last_saved2/diablo_th'  &
	  //CHAR(MOD(N,100)/10+48)          &
	  //CHAR(MOD(N,10)+48) //'_'        &
	  //CHAR(MOD(I,1000)/100+48)        & 
	  //CHAR(MOD(I,100)/10+48)          &
	  //CHAR(MOD(I,10)+48) // '.res'
      ENDDO

    ELSE
      FNAME='last_saved2/diablo_'       &
	  //CHAR(MOD(I,1000)/100+48)  &
	  //CHAR(MOD(I,100)/10+48)    &
	  //CHAR(MOD(I,10)+48) // '.saved'

      DO N=1,N_TH
	FNAME_TH(N)='last_saved2/diablo_th'  &
	  //CHAR(MOD(N,100)/10+48)          &
	  //CHAR(MOD(N,10)+48) //'_'        &
	  //CHAR(MOD(I,1000)/100+48)        & 
	  //CHAR(MOD(I,100)/10+48)          & 
	  //CHAR(MOD(I,10)+48) // '.saved'
      ENDDO
    ENDIF   

      
    WRITE(6,*) 'Writing flow to ',FNAME

    DO N=1,N_TH
      OPEN(UNIT=11,FILE=FNAME_TH(N),STATUS="UNKNOWN" &
	,FORM="UNFORMATTED")
      WRITE(11) NX, NY, NZ, TIME, TIME_STEP
      WRITE(11) (((CTHX(I,K,J,N),I=0,NX2P),K=0,NZ+1),J=0,NY+1)
      CLOSE(11)
    ENDDO

    CLOSE(11)

  RETURN
END


! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE SAVE_FLOW(FINAL)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
  USE mpi_var, 		ONLY : rank      

    IMPLICIT NONE

    CHARACTER*35 FNAME
    CHARACTER*35 FNAME_TH(N_TH)
    INTEGER      I, J, K, N, no_p,ny_p
    LOGICAL      FINAL, dir_exists
    PARAMETER  (no_p=4, ny_p=66)
    
    INQUIRE(DIRECTORY='./last_saved/.', EXIST=dir_exists) 

    IF ( (.NOT. dir_exists)) THEN
      WRITE(6,*) 'last_saved directory does not exists'
      CALL system('mkdir last_saved')
      WRITE(6,*) 'last_saved is created!'
    ENDIF    
    
    
      
    I=1+RANK
    IF (FINAL) THEN
      FNAME='last_saved/diablo_'      &
	  //CHAR(MOD(I,1000)/100+48)  & 
	  //CHAR(MOD(I,100)/10+48)    &
	  //CHAR(MOD(I,10)+48) // '.res'


      DO N=1,N_TH
	FNAME_TH(N)='last_saved/diablo_th'  &
	  //CHAR(MOD(N,100)/10+48)          &
	  //CHAR(MOD(N,10)+48) //'_'        &
	  //CHAR(MOD(I,1000)/100+48)        & 
	  //CHAR(MOD(I,100)/10+48)          &
	  //CHAR(MOD(I,10)+48) // '.res'
      END DO

    ELSE
      FNAME='last_saved/diablo_'       &
	  //CHAR(MOD(I,1000)/100+48)  &
	  //CHAR(MOD(I,100)/10+48)    &
	  //CHAR(MOD(I,10)+48) // '.saved'

      DO N=1,N_TH
	FNAME_TH(N)='last_saved/diablo_th'  &
	  //CHAR(MOD(N,100)/10+48)          &
	  //CHAR(MOD(N,10)+48) //'_'        &
	  //CHAR(MOD(I,1000)/100+48)        & 
	  //CHAR(MOD(I,100)/10+48)          & 
	  //CHAR(MOD(I,10)+48) // '.saved'
      ENDDO
    ENDIF

    
    WRITE(6,*) 'Writing flow to ',FNAME

    OPEN(UNIT=10,FILE=FNAME,STATUS="UNKNOWN",FORM="UNFORMATTED")
    WRITE(10) NX, NY, NZ, TIME, TIME_STEP


    WRITE(10) (((CU1X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1), &
	      (((CU2X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1),  &
	      (((CU3X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1), &
	      (((CPX(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1)
    DO N=1,N_TH
      OPEN(UNIT=11,FILE=FNAME_TH(N),STATUS="UNKNOWN" &
	,FORM="UNFORMATTED")
      WRITE(11) NX, NY, NZ, TIME, TIME_STEP
      WRITE(11) (((CTHX(I,K,J,N),I=0,NX2P),K=0,NZ+1),J=0,NY+1)
      CLOSE(11)
    ENDDO
      
    CLOSE(10)
    CLOSE(11)

  RETURN
END

!bob
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE BACKUP_FLOW  
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
  USE mpi_var, 		ONLY : rank      

    IMPLICIT NONE

    CHARACTER*45 FNAME
    CHARACTER*45 FNAME_TH(N_TH)
    CHARACTER*5 file_num
    INTEGER      I, J, K, N, m, no_p,ny_p
    LOGICAL      dir_exists
    PARAMETER  (no_p=4, ny_p=66)

            
	   m = time_step/BACKUP_FLOW_INT
            file_num = CHAR(MOD(m,10000)/1000+48) &
                  //CHAR(MOD(m,1000)/100+48)     &
                  //CHAR(MOD(m,100)/10+48)       &
                  //CHAR(MOD(m,10)+48)
    
    INQUIRE(DIRECTORY='./backup', EXIST=dir_exists) 

    IF ( (.NOT. dir_exists)) THEN
      WRITE(6,*) 'backup directory does not exists'
      CALL system('mkdir backup')
      WRITE(6,*) 'backup is created!'
    ENDIF    
    
    
      
    I=1+RANK
     FNAME='backup/diablo_'       &
	  //CHAR(MOD(I,1000)/100+48)  &
	  //CHAR(MOD(I,100)/10+48)    &
	  //CHAR(MOD(I,10)+48) // '.saved_' //file_num//''

     DO N=1,N_TH
	FNAME_TH(N)='backup/diablo_th'  &
	  //CHAR(MOD(N,100)/10+48)          &
	  //CHAR(MOD(N,10)+48) //'_'        &
	  //CHAR(MOD(I,1000)/100+48)        & 
	  //CHAR(MOD(I,100)/10+48)          & 
	  //CHAR(MOD(I,10)+48) // '.saved_'  //file_num//''
     ENDDO

    
    WRITE(6,*) 'Writing backup flow to ',FNAME

    OPEN(UNIT=10,FILE=FNAME,STATUS="UNKNOWN",FORM="UNFORMATTED")
    WRITE(10) NX, NY, NZ, TIME, TIME_STEP


    WRITE(10) (((CU1X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1), &
	      (((CU2X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1),  &
	      (((CU3X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1), &
	      (((CPX(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1)
    DO N=1,N_TH
      OPEN(UNIT=11,FILE=FNAME_TH(N),STATUS="UNKNOWN" &
	,FORM="UNFORMATTED")
      WRITE(11) NX, NY, NZ, TIME, TIME_STEP
      WRITE(11) (((CTHX(I,K,J,N),I=0,NX2P),K=0,NZ+1),J=0,NY+1)
      CLOSE(11)
    ENDDO
      
    CLOSE(10)
    CLOSE(11)

  RETURN
END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE READ_FLOW
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
  USE mpi_var, ONLY: rank      

    IMPLICIT NONE

    CHARACTER*35 FNAME
    CHARACTER*35 FNAME_TH(N_TH)
    INTEGER I, J, K, N
    REAL*8  RNUM1,RNUM2,RNUM3, DAMP_FACT
    LOGICAL KICK_AFTER_READ 
    
    KICK_AFTER_READ = .FALSE.
    
    I = RANK+1

    FNAME ='last_saved/diablo_'                 &
	  //CHAR(MOD(I,1000)/100+48)        & 
	  //CHAR(MOD(I,100)/10+48)              &
	  //CHAR(MOD(I,10)+48) // '.saved'
    DO N=1,N_TH
	FNAME_TH(N)='last_saved/diablo_th'       & 
	  //CHAR(MOD(N,100)/10+48)               & 
	  //CHAR(MOD(N,10)+48) //'_'             & 
	  //CHAR(MOD(I,1000)/100+48)             & 
	  //CHAR(MOD(I,100)/10+48)               &
	  //CHAR(MOD(I,10)+48) // '.saved'
    ENDDO

      
    WRITE(6,*)   'Reading flow from velocity data ',FNAME
    WRITE(6,*)   'Reading flow from theta data ', FNAME_TH(1)

    OPEN (UNIT=10, FILE=FNAME, STATUS="OLD", FORM="UNFORMATTED")
    READ (10) NX_T, NY_T, NZ_T, TIME, TIME_STEP


    WRITE(*,*) 'NX_T, NY_T, NZ_T: ', NX_T , NY_T, NZ_T, TIME, TIME_STEP
	
    WRITE(*,*) 'Size', size(CU1X)
	
    IF ((NX .NE. NX_T) .OR. (NY .NE. NY_T) .OR. (NZ .NE. NZ_T))THEN
	  WRITE(*,*) NX_T, NY_T, NZ_T
	  STOP 'Error: old flowfield wrong dimensions. '
    ENDIF
         

    WRITE(*,*) 'READING FLOW'
    READ (10) (((CU1X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1), &
	      (((CU2X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1),  &
	      (((CU3X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1),  &
	      (((CPX(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1) 
    WRITE(*,*) 'READ VEL, PRESSURE'
    DO N=1,NUM_READ_TH
! Specify in input.dat which scalars are to be read
      write(*,*) N, READ_TH_INDEX(N)
      OPEN(UNIT=11,FILE=FNAME_TH(READ_TH_INDEX(N)),STATUS="OLD" &
	    ,FORM="UNFORMATTED")
      READ (11) NX_T, NY_T, NZ_T, TIME, TIME_STEP
      READ (11) (((CTHX(I,K,J,READ_TH_INDEX(N)) &
	    ,I=0,NX2P),K=0,NZ+1),J=0,NY+1)
      CLOSE(11)
    ENDDO
    
    CLOSE(10)
    CLOSE(11)

!	OPEN(1231,file='velocity_data_input.dat',form='formatted', status='unknown')    
!	      
!	DO J=0,NY+1
!	  DO K=0,NZ+1
!	    write(1231,121) CU3X(1,K,J),CU1X(1,K,J)              
!	  ENDDO
!	ENDDO 
!          
!	CLOSE(1231)
!121   FORMAT(2f16.12)
!	WRITE(6,*) 'VELOCITY feild has been write from input file'

! if you wanted to  disturb the flow after reading    
    IF (KICK_AFTER_READ) THEN
      IF ( RANK.EQ.0 )  WRITE (6,*) '**WARNING: flow has been KICKED after reading**'

      DO I=1,MIN(NX2P,NX2P_L)
	DO J=1,NY
	  DO K=1,NZ
	!  Now, give the velocity field a random perturbation

	    CALL RANDOM_NUMBER(RNUM1)
	    CALL RANDOM_NUMBER(RNUM2)
	    CALL RANDOM_NUMBER(RNUM3)
! 	      Adding a random number to initial contdition based on this energy spectrum: E(k) = (K/K0)^4 * exp (-2*(K/K0)^2), it can be -5/3K/K0 as well
!	    CU1X(I,K,J)=CU1X(I,K,J)+(RNUM1-0.5d0) * KICK * (5.d0*I/NX)**4.d0 * EXP(-2.d0 * (5.d0*I/NX)**2.d0 )
!	    CU2X(I,K,J)=CU2X(I,K,J)+(RNUM2-0.5d0) * KICK * (5.d0*I/NX)**4.d0 * EXP(-2.d0 * (5.d0*I/NX)**2.d0 )
!	    CU3X(I,K,J)=CU3X(I,K,J)+(RNUM3-0.5d0) * KICK * (5.d0*I/NX)**4.d0 * EXP(-2.d0 * (5.d0*I/NX)**2.d0 )

!           random perturbation is mostly added to all fluctuations, but it can be
!           added like 5-10% to the mean flow (in peridoc directions i.e.  CU1 and CU3 in my cases)
             CU1X(I,K,J)=CU1X(I,K,J)+(RNUM2-0.5d0) * KICK
             CU2X(I,K,J)=CU2X(I,K,J)+(RNUM2-0.5d0) * KICK
             CU3X(I,K,J)=CU3X(I,K,J)+(RNUM3-0.5d0) * KICK

	    
	    IF (J .EQ. 0) THEN
	      DAMP_FACT = 1.0
	    ELSE
!	      DAMP_FACT = exp(-1.0d0*dble((J-1))/dble((NY-1)))
              DAMP_FACT = exp(-6 * dble((J)**2)/dble((NY)**2))
!	      IF (J.GE.NY/2) THEN 
!		DAMP_FACT = exp(-1.d0*dble((NY-J-1))/dble((NY-1)))
!	      ENDIF
	    ENDIF
	    
	     CU1X(I,K,J)=CU1X(I,K,J)*DAMP_FACT
	     CU2X(I,K,J)=CU2X(I,K,J)*DAMP_FACT
	     CU3X(I,K,J)=CU3X(I,K,J)*DAMP_FACT
   
	  ENDDO
	ENDDO
      ENDDO

!     random perturbation is mostly added to all fluctuations, but it can be
!     added like 5-10% to the mean flow (in peridoc directions i.e.  CU1 and CU3 in my cases)
!     Just 50% of grids in vertical driction
      DO J=0,0.5*NY
        DO K=1,NZ
          CALL RANDOM_NUMBER(RNUM1)
          CALL RANDOM_NUMBER(RNUM2)
          KICK=0.10
          !CU1X(0,K,J)=CU1X(0,K,J)*(1.0+(RNUM1-0.5d0) * KICK * exp(-10 *dble((J)**2)/dble((NY)**2)))
          !CU3X(0,K,J)=CU3X(0,K,J)*(1.0+(RNUM2-0.5d0) * KICK * exp(-10 *dble((J)**2)/dble((NY)**2)))
        ENDDO
      ENDDO 

    ENDIF


! C Apply initial boundary conditions, set ghost cells
    IF (IC_TYPE .EQ. 7)THEN
    
    ELSE 
      CALL APPLY_BC_VEL_LOWER
      CALL APPLY_BC_VEL_UPPER
      CALL APPLY_BC_VEL_LEFT
      CALL APPLY_BC_VEL_RIGHT
    ENDIF

  RETURN
END


! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE FILTER(n)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
USE TIME_STEP_VAR

      INTEGER n
      n = 1


  RETURN
END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE JACOBI_TRANS
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
  USE les_chan_var
  USE mpi_var, ONLY : rank
      
    IMPLICIT NONE

    INTEGER  I,J,K
    
    REAL*8  b,co_b,xx

    REAL(r8),ALLOCATABLE,DIMENSION(:,:,:) :: x_zeta, y_zeta, &
	x_eta, y_eta, zeta_x, zeta_y, eta_x, eta_y, IN_JACO

    ALLOCATE (x_zeta(NZ+1,NY+1,2), y_zeta(NZ+1,NY+1,2))
    ALLOCATE (x_eta(NZ+1,NY+1,2), y_eta(NZ+1,NY+1,2))
    ALLOCATE (zeta_x(NZ+1,NY+1,2), zeta_y(NZ+1,NY+1,2))
    ALLOCATE (eta_x(NZ+1,NY+1,2), eta_y(NZ+1,NY+1,2))
    ALLOCATE (IN_JACO(NZ+1,NY+1,2))
    ! allocate DYY just in case of wall model use
    IF ( W_BC_YMIN.EQ.3 ) THEN 
      ALLOCATE (DYY(0:NZ+1,0:NY+1))
    ENDIF

      
    IF(rank.EQ.0) WRITE(6,*) 'Jacobian tranformation is starting'
    open(21,file='GRID_YZ.dat',form='formatted',status='old') 
    open(31,file='density_data_input.dat',form='formatted', &
	status='old')
    
    DO J=0,NY+1
      DO K=0,NZ+1
	READ(21,*) xpoint(K,J),ypoint(K,J)
      ENDDO
    ENDDO
    
    CLOSE(21)
    
    DO J=1,NY+1
      DY(J) = ypoint(JSTART+2,J) - ypoint(JSTART+2,J-1)
    ENDDO

!DYY is a parameter used as DY in wall function
    IF ( W_BC_YMIN.EQ.3 ) THEN
      DO J=1,NY+1
	DO K=1,NZ+1
	  DYY(K,J) = ypoint(K,J) - ypoint(K,J-1)
	ENDDO
      ENDDO
    ENDIF
      
!       DO J=0,NY+1
!        DO K=0,NZ+1
!          TH_BAR(K,J)= (-ypoint(k,j)+ypoint(0,NY+1))
!        ENDDO
!       ENDDO

       
    DO J=0,NY+1
      DO K=0,NZ+1
	IF (CONT_STRAT) THEN 
	  TH_BAR(K,J) = (-ypoint(k,j)+ypoint(0,NY+1))
	ELSE
	  READ(31,*) TH_BAR(K,J)
	ENDIF         
      ENDDO
    ENDDO
      
    CLOSE(31)
    
111   format(2f16.8)
    
    DO J=1,NY
      DO K=1,NZ
	x_zeta(K,J,1) = xpoint(K+1,J) - xpoint(K,J)

	x_zeta(K,J,2) = 0.250d0 * ( xpoint(K+1,J+1) + xpoint(K+1,J) &
				  - xpoint(K-1,J+1) - xpoint(K-1,J) ) 

	y_zeta(K,J,1) = ypoint(K+1,J) - ypoint(K,J)
	
	y_zeta(K,J,2) = 0.250d0 * ( ypoint(K+1,J+1) + ypoint(K+1,J) &
				  - ypoint(K-1,J+1) - ypoint(K-1,J) )


	x_eta(K,J,2) = xpoint(K,J+1) - xpoint(K,J)

	x_eta(K,J,1) = 0.250d0 * ( xpoint(K,J+1) + xpoint(K+1,J+1) &
				  - xpoint(K,J-1) - xpoint(K+1,J-1) )
	
	y_eta(K,J,1) = 0.250d0 * ( ypoint(K,J+1) + ypoint(K+1,J+1)  &
				  - ypoint(K,J-1) - ypoint(K+1,J-1) )
	
	y_eta(K,J,2) = (ypoint(K,J+1)- ypoint(K,J)) 
	
      ENDDO
    ENDDO


    DO K=1,NZ+1
      x_eta(K,1,1)  = 2.0d0 * x_eta(K,2,1) - x_eta(K,3,1)
      DO I=1,2
      x_eta(K,NY+1,I)  = 2.0d0 * x_eta(K,NY,I) - x_eta(K,NY-1,I)
      x_zeta(K,NY+1,I) = 2.0d0 * x_zeta(K,NY,I) - x_zeta(K,NY-1,I) 
      ENDDO
      y_eta(K,1,1)  = 2.0d0 * y_eta(K,2,1) - y_eta(K,3,1)
      DO I=1,2
      y_eta(K,NY+1,I)  = 2.0d0 * y_eta(K,NY,I) - y_eta(K,NY-1,I)
      y_zeta(K,NY+1,I) = 2.0d0 * y_zeta(K,NY,I) - y_zeta(K,NY-1,I)
      ENDDO   
    ENDDO  

    DO J=1,NY+1
      x_zeta(1,J,2) = 2.0d0 * x_zeta(2,J,2) - x_zeta(3,J,2)
      DO I=1,2
	x_zeta(NZ+1,J,I) = 2.0d0*x_zeta(NZ,J,I)- x_zeta(NZ-1,J,I)
	x_eta(NZ+1,J,I) = 2.0d0*x_eta(NZ,J,I)- x_eta(NZ-1,J,I)
      ENDDO
      y_zeta(1,J,2)    = 2.0d0*y_zeta(2,J,2) - y_zeta(3,J,2)
      DO I=1,2
	y_zeta(NZ+1,J,I) = 2.0d0*y_zeta(NZ,J,I) - y_zeta(NZ-1,J,I)
	y_eta(NZ+1,J,I) = 2.0d0*y_eta(NZ,J,I) - y_eta(NZ-1,J,I)
      ENDDO 	 
    ENDDO

! C make smooth transformation
    DO J=1,NY+1
      DO K=1,NZ+1
	DO I=1,2
	  x_zeta(K,J,I) = dble(AINT(x_zeta(K,J,I)*10**8)) * dble(10.0d0**(-8))
	  x_eta(K,J,I)  = dble(AINT(x_eta(K,J,I)*10**8))  * dble(10.0d0**(-8))
	  
	  y_zeta(K,J,I) = dble(AINT(y_zeta(K,J,I)*10**8)) * dble(10.0d0**(-8))
	  y_eta(K,J,I)  = dble(AINT(y_eta(K,J,I)*10**8)) * dble(10.0d0**(-8))
	ENDDO
      ENDDO 	  
    ENDDO 	  
    
    
    DO J=1,NY+1
      DO K=1,NZ+1
      DO I = 1,2
	IN_JACO(K,J,I) = x_zeta(K,J,I) * y_eta(K,J,I) &
			- x_eta(K,J,I) * y_zeta(K,J,I)


	zeta_x(K,J,I) =   y_eta(K,J,I) / IN_JACO(K,J,I)  
	eta_x(K,J,I)  =  -y_zeta(K,J,I) / IN_JACO(K,J,I) 

	zeta_y(K,J,I) =  -x_eta(K,J,I) / IN_JACO(K,J,I) 
	eta_y(K,J,I)  =   x_zeta(K,J,I) / IN_JACO(K,J,I)
      ENDDO

      DO I = 1,2
	GMAT_11(K,J,I)   =  ( zeta_x(K,J,I) * zeta_x(K,J,I) &
			    + zeta_y(K,J,I) * zeta_y(K,J,I) ) * IN_JACO(K,J,I)
			    
	GMAT_12(K,J,I)   =  ( zeta_x(K,J,I) * eta_x(K,J,I) &
			    + zeta_y(K,J,I) * eta_y(K,J,I) ) * IN_JACO(K,J,I)
	
	GMAT_22(K,J,I)   =  ( eta_x(K,J,I) * eta_x(K,J,I) &
			    + eta_y(K,J,I) * eta_y(K,J,I) ) * IN_JACO(K,J,I) 

	CJOB_11(K,J,I)   = zeta_x(K,J,I) * IN_JACO(K,J,I)  
	CJOB_12(K,J,I)   = zeta_y(K,J,I) * IN_JACO(K,J,I)
	CJOB_21(K,J,I)   = eta_x(K,J,I) * IN_JACO(K,J,I)
	CJOB_22(K,J,I)   = eta_y(K,J,I) * IN_JACO(K,J,I)  
      ENDDO

      ENDDO
    ENDDO
    
    
    
    DO J=1,NY
      DO K=1,NZ
	INT_JACOB(K,J)= 0.250d0 * ( IN_JACO(K,J,1) + IN_JACO(K,J,2) &
				  + IN_JACO(K,J+1,1) + IN_JACO(K+1,J,2) )
      ENDDO
    ENDDO
    
    DO J=0,NY+1
! c       INT_JACOB(NZ,J)   = 2.0*INT_JACOB(NZ-1,J) - INT_JACOB(NZ-2,J)
      INT_JACOB(NZ+1,J) = 2.0d0 * INT_JACOB(NZ,J) - INT_JACOB(NZ-1,J)
      INT_JACOB(0,J) = 2.0d0 * INT_JACOB(1,J) - INT_JACOB(2,J)
    ENDDO
    
    DO K=0,NZ+1
      INT_JACOB(K,NY+1) = 2.0d0 * INT_JACOB(K,NY) - INT_JACOB(K,NY-1)
      INT_JACOB(K,0) = 2.0d0 * INT_JACOB(K,1) - INT_JACOB(K,2)
    ENDDO 
    
! c    Correction for 0.0
    INT_JACOB(0,0)  = 0.5d0 * ( 2.0d0 * INT_JACOB(1,0) - INT_JACOB(2,0) &
		    + 2.0d0 * INT_JACOB(0,1) - INT_JACOB(0,2) ) 
    
    INT_JACOB(NZ+1,NY+1)  = 0.5d0 * ( 2.0d0 * INT_JACOB(NZ,NY+1)  &
			    - INT_JACOB(NZ-1,NY+1) +  2.0d0 * INT_JACOB(NZ+1,NY) &
			    - INT_JACOB(NZ+1,NY-1) )   
    
    DO J=0,NY+1
      DO I=1,2
      CJOB_11(0,J,I) = 2.d0 * CJOB_11(1,J,I) - CJOB_11(2,J,I)
      CJOB_12(0,J,I) = 2.d0 * CJOB_12(1,J,I) - CJOB_12(2,J,I)
      CJOB_21(0,J,I) = 2.d0 * CJOB_21(1,J,I) - CJOB_21(2,J,I)
      CJOB_22(0,J,I) = 2.d0 * CJOB_22(1,J,I) - CJOB_22(2,J,I)
      
      GMAT_11(0,J,I) = 2.d0 * GMAT_11(1,J,I) - GMAT_11(2,J,I)
      GMAT_12(0,J,I) = 2.d0 * GMAT_12(1,J,I) - GMAT_12(2,J,I)
      GMAT_22(0,J,I) = 2.d0 * GMAT_22(1,J,I) - GMAT_22(2,J,I)
      ENDDO
    ENDDO
    
    DO K=0,NZ+1
      DO I=1,2
      CJOB_11(K,0,I) = 2.d0 * CJOB_11(K,1,I) - CJOB_11(K,2,I)
      CJOB_12(K,0,I) = 2.d0 * CJOB_12(K,1,I) - CJOB_12(K,2,I)
      CJOB_21(K,0,I) = 2.d0 * CJOB_21(K,1,I) - CJOB_21(K,2,I)
      CJOB_22(K,0,I) = 2.d0 * CJOB_22(K,1,I) - CJOB_22(K,2,I)
      
      GMAT_11(K,0,I) = 2.0d0 * GMAT_11(K,1,I) - GMAT_11(K,2,I)
      GMAT_12(K,0,I) = 2.0d0 * GMAT_12(K,1,I) - GMAT_12(K,2,I) 
      GMAT_22(K,0,I) = 2.0d0 * GMAT_22(K,1,I) - GMAT_22(K,2,I)
      ENDDO
    ENDDO
    
    DO I = 1,2
      GMAT_11(0,0,I) = 0.5d0 * ( 2.0d0 * GMAT_11(0,1,I) - GMAT_11(0,2,I) &
		      +  2.d0 * GMAT_11(1,0,I) - GMAT_11(2,0,I) )
		      
      GMAT_22(0,0,I) = 0.5d0 * ( 2.0d0 * GMAT_12(0,1,I) - GMAT_12(0,2,I) &
		      + 2.d0 * GMAT_12(1,0,I) - GMAT_12(2,0,I) )
		      
      GMAT_22(0,0,I) = 0.5d0 * ( 2.0d0 * GMAT_22(0,1,I) - GMAT_22(0,2,I) &
		      + 2.d0 * GMAT_22(1,0,I) - GMAT_22(2,0,I) )
	  
      CJOB_11(0,0,I) = 0.50d0 *( 2.d0 * CJOB_11(0,1,I) - CJOB_11(0,2,I) &
		      + 2.d0 * CJOB_11(1,0,I) - CJOB_11(2,0,I) )
	  
      CJOB_12(0,0,I) = 0.50d0 * (2.d0 * CJOB_12(0,1,I) - CJOB_12(0,2,I) &
		      + 2.d0 * CJOB_12(1,0,I) - CJOB_12(2,0,I) )
      
      CJOB_21(0,0,I) = 0.50d0 * (2.d0 * CJOB_21(0,1,I) - CJOB_21(0,2,I) &
		      + 2.d0 * CJOB_21(1,0,I) - CJOB_21(2,0,I) )
      
      CJOB_22(0,0,I) = 0.50d0 * (2.d0 * CJOB_22(0,1,I) - CJOB_22(0,2,I) &
		      + 2.d0 * CJOB_22(1,0,I) - CJOB_22(2,0,I) )
    ENDDO
    
! CC Additinal matrix required for MG  
    DO I = 1,2  
    
      DO K=0,NZ+1
	CJOB_11(K,NY+2,I) = CJOB_11(K,NY+1,I)
	CJOB_12(K,NY+2,I) = CJOB_12(K,NY+1,I)
	CJOB_21(K,NY+2,I) = CJOB_21(K,NY+1,I)
	CJOB_22(K,NY+2,I) = CJOB_22(K,NY+1,I)
      
	GMAT_11(K,NY+2,I) = GMAT_11(K,NY+1,I)
	GMAT_12(K,NY+2,I) = GMAT_12(K,NY+1,I)
	GMAT_22(K,NY+2,I) = GMAT_22(K,NY+1,I)
      ENDDO
      
      DO J=0,NY+1
	CJOB_11(NZ+2,J,I) = CJOB_11(NZ+1,J,I)
	CJOB_12(NZ+2,J,I) = CJOB_12(NZ+1,J,I)
	CJOB_21(NZ+2,J,I) = CJOB_21(NZ+1,J,I)
	CJOB_22(NZ+2,J,I) = CJOB_22(NZ+1,J,I)
      
	GMAT_11(NZ+2,J,I) = GMAT_11(NZ+1,J,I)
	GMAT_12(NZ+2,J,I) = GMAT_12(NZ+1,J,I)
	GMAT_22(NZ+2,J,I) = GMAT_22(NZ+1,J,I)
      ENDDO
      
      GMAT_11(NZ+2,NY+2,I) = GMAT_11(NZ+1,NY+1,I)
      GMAT_12(NZ+2,NY+2,I) = GMAT_12(NZ+1,NY+1,I)
      GMAT_22(NZ+2,NY+2,I) = GMAT_22(NZ+1,NY+1,I)
    
      CJOB_11(NZ+2,NY+2,I) = CJOB_11(NZ+1,NY+1,I)
      CJOB_12(NZ+2,NY+2,I) = CJOB_12(NZ+1,NY+1,I)
      CJOB_21(NZ+2,NY+2,I) = CJOB_21(NZ+1,NY+1,I)
      CJOB_22(NZ+2,NY+2,I) = CJOB_22(NZ+1,NY+1,I)
    ENDDO 	
    

    IF (LES) THEN

    DO J=1,NY+1
      DO K=1,NZ+1
      DO I=1,2
	GMAT_11_z(K,J,I)   =  ( 2.0d0 * zeta_x(K,J,I) * zeta_x(K,J,I) &
			      + zeta_y(K,J,I) * zeta_y(K,J,I) ) * IN_JACO(K,J,I)
		    
	GMAT_12_z(K,J,I)   =  ( 2.0d0 * zeta_x(K,J,I) * eta_x(K,J,I) &
			      + zeta_y(K,J,I) * eta_y(K,J,I) ) * IN_JACO(K,J,I)
		    
	GMAT_22_z(K,J,I)   =  ( 2.0d0 * eta_x(K,J,I) * eta_x(K,J,I) &
			      + eta_y(K,J,I) * eta_y(K,J,I) ) * IN_JACO(K,J,I)
		    
	GMAT_11_y(K,J,I)   =  ( zeta_x(K,J,I) * zeta_x(K,J,I) &
			      + 2.0d0 * zeta_y(K,J,I) * zeta_y(K,J,I) ) * IN_JACO(K,J,I)
		    
	GMAT_12_y(K,J,I)   =  ( zeta_x(K,J,I) * eta_x(K,J,I) &
			      + 2.0d0 * zeta_y(K,J,I) * eta_y(K,J,I) ) * IN_JACO(K,J,I)
		    
	GMAT_22_y(K,J,I)   =  ( eta_x(K,J,I) * eta_x(K,J,I) &
			      + 2.0d0 * eta_y(K,J,I) * eta_y(K,J,I) ) * IN_JACO(K,J,I)
      ENDDO
      ENDDO
    ENDDO

    DO J=0,NY+1
      DO I=1,2
      GMAT_11_z(0,J,I) = 2.d0 * GMAT_11_z(1,J,I) - GMAT_11_z(2,J,I)
      GMAT_12_z(0,J,I) = 2.d0 * GMAT_12_z(1,J,I) - GMAT_12_z(2,J,I)
      GMAT_22_z(0,J,I) = 2.d0 * GMAT_22_z(1,J,I) - GMAT_22_z(2,J,I)

      GMAT_11_y(0,J,I) = 2.d0 * GMAT_11_y(1,J,I) - GMAT_11_y(2,J,I)
      GMAT_12_y(0,J,I) = 2.d0 * GMAT_12_y(1,J,I) - GMAT_12_y(2,J,I)
      GMAT_22_y(0,J,I) = 2.d0 * GMAT_22_y(1,J,I) - GMAT_22_y(2,J,I)
      ENDDO
    ENDDO

    DO K=0,NZ+1
      DO I=1,2
      GMAT_11_z(K,0,I) = 2.0d0 * GMAT_11_z(K,1,I) - GMAT_11_z(K,2,I)
      GMAT_12_z(K,0,I) = 2.0d0 * GMAT_12_z(K,1,I) - GMAT_12_z(K,2,I)
      GMAT_22_z(K,0,I) = 2.0d0 * GMAT_22_z(K,1,I) - GMAT_22_z(K,2,I)

      GMAT_11_y(K,0,I) = 2.0d0 * GMAT_11_y(K,1,I) - GMAT_11_y(K,2,I)
      GMAT_12_y(K,0,I) = 2.0d0 * GMAT_12_y(K,1,I) - GMAT_12_y(K,2,I)
      GMAT_22_y(K,0,I) = 2.0d0 * GMAT_22_y(K,1,I) - GMAT_22_y(K,2,I)
      ENDDO
    ENDDO

    DO I = 1,2
      GMAT_11_z(0,0,I) = 0.5d0 * ( 2.0d0 * GMAT_11_z(0,1,I) - GMAT_11_z(0,2,I) &
				  + 2.d0 * GMAT_11_z(1,0,I) - GMAT_11_z(2,0,I) )
				  
      GMAT_22_z(0,0,I) = 0.5d0 * ( 2.0d0 * GMAT_12_z(0,1,I) - GMAT_12_z(0,2,I) &
				  + 2.d0 * GMAT_12_z(1,0,I) - GMAT_12_z(2,0,I) )
				  
      GMAT_22_z(0,0,I) = 0.5d0 * ( 2.0d0 * GMAT_22_z(0,1,I) - GMAT_22_z(0,2,I) &
				  + 2.d0 * GMAT_22_z(1,0,I) - GMAT_22_z(2,0,I) )

      GMAT_11_y(0,0,I) = 0.5d0 * ( 2.0d0 * GMAT_11_y(0,1,I) - GMAT_11_y(0,2,I) &
				  + 2.d0 * GMAT_11_y(1,0,I) - GMAT_11_y(2,0,I) )
		
      GMAT_22_y(0,0,I) = 0.5d0 * ( 2.0d0 * GMAT_12_y(0,1,I) - GMAT_12_y(0,2,I) &
				  + 2.d0 * GMAT_12_y(1,0,I) - GMAT_12_y(2,0,I) )
		
      GMAT_22_y(0,0,I) = 0.5d0 * ( 2.0d0 * GMAT_22_y(0,1,I) - GMAT_22_y(0,2,I) &
				  + 2.d0 * GMAT_22_y(1,0,I) - GMAT_22_y(2,0,I) )
    ENDDO

    DO I = 1,2
      DO K=0,NZ+1
	GMAT_11_z(K,NY+2,I) = GMAT_11_z(K,NY+1,I)
	GMAT_12_z(K,NY+2,I) = GMAT_12_z(K,NY+1,I)
	GMAT_22_z(K,NY+2,I) = GMAT_22_z(K,NY+1,I)

	GMAT_11_y(K,NY+2,I) = GMAT_11_y(K,NY+1,I)
	GMAT_12_y(K,NY+2,I) = GMAT_12_y(K,NY+1,I)
	GMAT_22_y(K,NY+2,I) = GMAT_22_y(K,NY+1,I)
      ENDDO
      
      DO J=0,NY+1
	GMAT_11_z(NZ+2,J,I) = GMAT_11_z(NZ+1,J,I)
	GMAT_12_z(NZ+2,J,I) = GMAT_12_z(NZ+1,J,I)
	GMAT_22_z(NZ+2,J,I) = GMAT_22_z(NZ+1,J,I)

	GMAT_11_y(NZ+2,J,I) = GMAT_11_y(NZ+1,J,I)
	GMAT_12_y(NZ+2,J,I) = GMAT_12_y(NZ+1,J,I)
	GMAT_22_y(NZ+2,J,I) = GMAT_22_y(NZ+1,J,I)
      ENDDO
      
      GMAT_11_z(NZ+2,NY+2,I) = GMAT_11_z(NZ+1,NY+1,I)
      GMAT_12_z(NZ+2,NY+2,I) = GMAT_12_z(NZ+1,NY+1,I)
      GMAT_22_z(NZ+2,NY+2,I) = GMAT_22_z(NZ+1,NY+1,I)

      GMAT_11_y(NZ+2,NY+2,I) = GMAT_11_y(NZ+1,NY+1,I)
      GMAT_12_y(NZ+2,NY+2,I) = GMAT_12_y(NZ+1,NY+1,I)
      GMAT_22_y(NZ+2,NY+2,I) = GMAT_22_y(NZ+1,NY+1,I)
    ENDDO

    ENDIF

      
! C111   format(2f12.8)      
	
112   format(f16.10)      
            
    WAVE_ABS = .FALSE.
    ART_VISCOSITY = 500.d0	!ART_VISCOSITY * VISCOSITY  

    SPONGE_SIGMA_OUT(:) = 0.d0

! C     Calculating sponge in Strea-Wise direction   

    IF(SPONGE_TYPE .EQ. 1 .AND. (NZ_S_R .LT. NZ+1))	THEN	!Parabolic Sponge

      b    = 1.d0 / REAL( xpoint(NZ+1,NY+1)-xpoint(NZ_S_R,NY+1) )
      
      DO K = 0, NZ+1
	IF ( K .GE. NZ_S_R ) THEN
		xx = ( xpoint(K,NY+1) - xpoint(NZ_S_R,NY+1) ) * b
		SPONGE_SIGMA_OUT(K) = xx**2		! 2*xx-xx**2
	ELSEIF( K .LE. NZ_S_L ) THEN
		xx = ( xpoint(NZ_S_L ,NY+1) - xpoint(K,NY+1) ) * b
		SPONGE_SIGMA_OUT(K) = xx**2		! 2*xx-xx**2
	ENDIF
      ENDDO

    ELSEIF((SPONGE_TYPE .EQ. 2) .AND. (NZ_S_R .LT. NZ+1) )	 THEN		!Exponential Sponge
!     co_b higher -> smoother sponge -> less effective 
      co_b  = 25.d0
!      b    = log(co_b) / ( REAL(xpoint(NZ+1,NY+1)-xpoint(NZ_S_R,NY+1)) )
      b = log(co_b)/(REAL(NZ+1-NZ_S_R))
      
      DO K = 0,NZ+1
	IF ( K .GE. NZ_S_R ) THEN
	!SPONGE_SIGMA_OUT(K) = ( exp(b*(DABS(xpoint(K,NY+1)-xpoint(NZ_S_R,NY+1)))) -1.d0 ) / (co_b-1.d0)
         SPONGE_SIGMA_OUT(K) = ( exp(b*ABS(K-NZ_S_R)) -1.0d0 ) / (co_b-1.d0)
	ELSEIF( K .LE. NZ_S_L ) THEN
	!SPONGE_SIGMA_OUT(K) = ( exp(b*(DABS(xpoint(K,NY+1)-xpoint(NZ_S_L ,NY+1)))) -1.d0 ) / (co_b-1.d0)
        SPONGE_SIGMA_OUT(K) = ( exp(b*ABS(NZ_S_L-K))-1.0d0 ) / (co_b-1.d0)
	ENDIF
      ENDDO
	
    ENDIF
    
! C	END OF SPONGE in stream wise direction
	  	      
!     IF(rank .eq. 0) THEN
! 	  DO K = 0, NZ+1
! 		  WRITE(666,*) xpoint(K,NY+1),SPONGE_SIGMA_OUT(K)
! 	  ENDDO
!     ENDIF
          
! C     Calculating sponge in Vertical direction
    SPONGE_SIGMA(:) = 0.d0
      
    IF(SPONGE_TYPE .EQ. 1)	THEN	!Parabolic Sponge
      
      b    = 1.d0 / REAL( ypoint(NZ+1,NY+1)-ypoint(NZ+1,NY_S_T) )
      
      DO J = NY_S_T, NY+1
 	xx = ( ypoint(NZ+1,J) - ypoint(NZ+1,NY_S_T) ) * b
 	SPONGE_SIGMA(J) = xx**2
      ENDDO
      
    ELSEIF(SPONGE_TYPE .EQ. 2)	THEN		!Exponential Sponge
!     co_b higher -> smoother sponge -> less effective 
      co_b  = 30.d0
      
      !b    = log(co_b)/(REAL(ypoint(NZ+1,NY+1)-ypoint(NZ+1,NY_S_T)))
      b    = log(co_b)/(REAL(NY+1-NY_S_T))
      
      DO J = NY_S_T, NY+1
 	!SPONGE_SIGMA(J) = ( exp(b*(DABS(ypoint(NZ+1,J) - ypoint(NZ+1,NY_S_T)))) -1.0d0 ) / (co_b-1.d0)
 	SPONGE_SIGMA(J) = ( exp(b*(J-NY_S_T)) -1.0d0 ) / (co_b-1.d0)
      ENDDO
    ENDIF
! C	END OF SPONGE in Vertical direction

!     IF(rank .eq. 0) THEN
!       DO J = 0, NY+1
! 	WRITE(667,*) ypoint(NZ+1,J),SPONGE_SIGMA(J)
!       ENDDO
!     ENDIF

    DEALLOCATE (x_zeta, y_zeta, x_eta, y_eta, zeta_x, zeta_y, eta_x, eta_y, IN_JACO)
      
    IF(rank.EQ.0) write(6,*) 'Jacobian Tranformation is Done'
      
  RETURN
END 

! !--------------------------------------------------------------
!             SUBROUTINE open_MP_initialization
! !--------------------------------------------------------------
! USE ntypes
! USE Domain
! USE Grid
! USE Fft_var
! USE TIME_STEP_VAR
! USE run_variable
! USE omp_lib
! USE mpi_var, ONLY : rank
!
! implicit none
! 
! !--.---------.---------.---------.---------.---------.-|-------|
! !$OMP PARALLEL PRIVATE(iam, nthreads, chunk)
! ! Compute the subset of iterations
! ! executed by each thread
! 
!       nthreads = omp_get_num_threads()
!       iam = omp_get_thread_num()
!       chunk = (NKX+1)/nthreads
! 
!       kstart = iam * chunk
!       kend = min((iam + 1) * chunk-1, NKX)
! 
!       if(iam == nthreads-1)then
!       kend = NKX
!       endif
! 
!       chunk = (NXM+1)/nthreads
! 
!       k_start = iam * chunk
!       k_end = min((iam + 1) * chunk-1, NXM)
! 
!       if(iam == nthreads-1)then
!       k_end = NXM
!       endif
! 
!       chunk = (JEND-JSTART+1)/nthreads
! 
!       j_start = iam * chunk +1 
!       j_end = min((iam + 1) * chunk, JEND)
! 
!       if(iam == nthreads-1)then
!       j_end = JEND
!       endif
! 
!       chunk = (NY+2)/nthreads
!          
!       jj_start = iam * chunk  
!       jj_end = min((iam + 1) * chunk-1, NY+1)
!             
!       if(iam == nthreads-1)then
!       jj_end = NY+1
!       endif 
! 
! !       write(6,*)'LOOPs', iam,'kstart',kstart,'kend',kend
! !      write(6,*)'LOOPs',iam,'k_start',k_start,'k_end',k_end
! !      write(6,*)'LOOPs',iam,'j_start',j_start,'j_end',j_end
!       if (rank .eq. 0) then 
!        write(6,*)'LOOPs',iam,'jj_start',jj_start,'jj_end',jj_end
!       endif
! !$omp end parallel
! 
!       write(6,*)'Threads are intialized'
! 
!       RETURN
!       END
 
