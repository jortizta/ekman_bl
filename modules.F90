subroutine global_allocation  

  RETURN
END

!@c
!DATA TYPES
MODULE ntypes

  INTEGER, PARAMETER :: r8=8 
  INTEGER, PARAMETER :: i4=4 

 END MODULE ntypes 

!DOMAIN
MODULE Domain
  USE ntypes

    INTEGER(i4)  :: 	nx, ny, nz, N_TH, NXM, NYM, NZM, TNKZ, TNKY, NZ_S_R, NZ_S_L, NY_S_T 
    INTEGER(i4)  :: 	NKX
    INTEGER(i4)  :: 	NP, NXP, NZP, NXV, NZV, NKXV, NKXP, NX2V, NX2P, NXP_L, NX2P_L

!c mudpack multigrid values 
 
    INTEGER(i4)  ::	iixp, jjyq, iiex, jjey, nnx, nny, llwork
    
    INCLUDE 	'grid_def'
    INCLUDE	'mpif.h'

END MODULE Domain

!GRID
MODULE Grid

  USE ntypes
    REAL(r8)   :: LX, LY, LZ, CSX, CSY, CSZ         !Length
    INTEGER    :: jstart, jend , ZSTART,  ZEND      ! 
    
    REAL(r8),ALLOCATABLE,DIMENSION(:) :: GX,  GY,  GZ,  DX,  DY,  DZ, &
					  GXF, GYF, GZF, DXF, DYF, DZF
    REAL(r8),ALLOCATABLE,DIMENSION(:,:) :: DYY, xpoint, ypoint
  
END MODULE Grid

!TIME_STEP_VAR
MODULE TIME_STEP_VAR

  USE ntypes
  
    REAL(r8) DELTA_T, DELTA_T_set
    INTEGER  N_TIME_STEPS
    

END MODULE TIME_STEP_VAR

!Fft
MODULE fft_var

  USE ntypes,	 ONLY : r8

!FFT plans
    INTEGER(8) :: FFTW_X_TO_P_PLAN, FFTW_X_TO_F_PLAN, &
		  FFTW_Y_TO_P_PLAN, FFTW_Y_TO_F_PLAN, &
		  FFTW_Z_TO_P_PLAN, FFTW_Z_TO_F_PLAN

    INTEGER    NKY, NKZ
    REAL(r8)   PI, EPS, RNX, RNY, RNZ
! PARAMETER  (NKX=NX/3) 

!FFT VARIABLES
!REAL FIELDS
    REAL(r8),ALLOCATABLE,DIMENSION(:,:,:)    :: FIELD_IN
    REAL(r8),ALLOCATABLE,DIMENSION(:,:)      :: PLANE_IN

!COMPLEX FIELDS
    COMPLEX(r8)                              :: CI  
    COMPLEX(r8),ALLOCATABLE,DIMENSION(:,:)   :: CZX_PLANE, CYZ_PLANE
    COMPLEX(r8),ALLOCATABLE,DIMENSION(:)     :: CIKX, CIKY, CIKZ,CIKXP

!WAVE INDICES  
    REAL(r8),ALLOCATABLE,DIMENSION(:)    :: kx, ky, kz,kxp

!WAVE NUMBERS 
    REAL(r8),ALLOCATABLE,DIMENSION(:)    :: rkx, rky, rkz

!MODIFIED WAVE NUMBERS
    REAL(r8),ALLOCATABLE,DIMENSION(:)    ::  kx2, ky2, kz2,kx2p
    REAL(r8) kx2_c
    
    common /fft_pln/ FFTW_X_TO_P_PLAN, FFTW_X_TO_F_PLAN
!$omp threadprivate(/fft_pln/)
END MODULE fft_var

!Run_vari
MODULE run_variable

  USE ntypes,	 ONLY : r8
  USE Domain
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Input parameters and runtime variables
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    LOGICAL  :: CONT_STRAT, MUDPACK_USE
    
    REAL(r8) :: NU, KICK, WBULK0, PX0,U_BC_XMIN_C1, U_BC_XMIN_C2, U_BC_XMIN_C3,       &
		V_BC_XMIN_C1,V_BC_XMIN_C2,V_BC_XMIN_C3, W_BC_XMIN_C1, W_BC_XMIN_C2,   &
		W_BC_XMIN_C3,U_BC_YMIN_C1,U_BC_YMIN_C2,U_BC_YMIN_C3, V_BC_YMIN_C1,    &
		V_BC_YMIN_C2, V_BC_YMIN_C3, W_BC_YMIN_C1, W_BC_YMIN_C2, W_BC_YMIN_C3, &
		U_BC_ZMIN_C1, U_BC_ZMIN_C2, U_BC_ZMIN_C3, V_BC_ZMIN_C1, V_BC_ZMIN_C2, &
		V_BC_ZMIN_C3, W_BC_ZMIN_C1, W_BC_ZMIN_C2, W_BC_ZMIN_C3, TH_BC_XMIN_C1(1:N_TH), &
		TH_BC_XMIN_C2(1:N_TH), TH_BC_XMIN_C3(1:N_TH), TH_BC_YMIN_C1(1:N_TH),  &
		TH_BC_YMIN_C2(1:N_TH), TH_BC_YMIN_C3(1:N_TH), TH_BC_ZMIN_C1(1:N_TH),  &
		TH_BC_ZMIN_C2(1:N_TH), TH_BC_ZMIN_C3(1:N_TH)

    REAL(r8)  ::  U_BC_XMAX_C1, U_BC_XMAX_C2, U_BC_XMAX_C3, &
		  V_BC_XMAX_C1, V_BC_XMAX_C2, V_BC_XMAX_C3, &
		  W_BC_XMAX_C1, W_BC_XMAX_C2, W_BC_XMAX_C3, &
		  U_BC_YMAX_C1, U_BC_YMAX_C2, U_BC_YMAX_C3, &
		  V_BC_YMAX_C1, V_BC_YMAX_C2, V_BC_YMAX_C3, &
		  W_BC_YMAX_C1, W_BC_YMAX_C2, W_BC_YMAX_C3, &
		  U_BC_ZMAX_C1, U_BC_ZMAX_C2, U_BC_ZMAX_C3, &
		  V_BC_ZMAX_C1, V_BC_ZMAX_C2, V_BC_ZMAX_C3, &
		  W_BC_ZMAX_C1, W_BC_ZMAX_C2, W_BC_ZMAX_C3, &
		  Qn1, Qn, Q0, WEIGHT, PXV

    REAL(r8)  ::  TH_BC_XMAX_C1(1:N_TH),TH_BC_XMAX_C2(1:N_TH), &
		  TH_BC_XMAX_C3(1:N_TH),TH_BC_YMAX_C1(1:N_TH), &
		  TH_BC_YMAX_C2(1:N_TH),TH_BC_YMAX_C3(1:N_TH), &
		  TH_BC_ZMAX_C1(1:N_TH), TH_BC_ZMAX_C2(1:N_TH), &
		  TH_BC_ZMAX_C3(1:N_TH), CFL,DFN,INT_PI,RI_FINAL, &
		  dtc,C_avg


    REAL(r8) :: ART_VISCOSITY
    
    REAL(r8),ALLOCATABLE,DIMENSION(:) :: SPONGE_SIGMA,SPONGE_SIGMA_OUT, &
					  C_int, C_int_le
    
    INTEGER :: NY_S, NX_T,NY_T,NZ_T, SPONGE_TYPE

    INTEGER :: 	VERBOSITY, &  
		SAVE_FLOW_INT, SAVE_STATS_INT, IC_TYPE, F_TYPE,COUNT_DATA, BACKUP_FLOW_INT,&
		U_BC_XMIN, V_BC_XMIN, W_BC_XMIN, TH_BC_XMIN(1:N_TH), &
		U_BC_XMAX, V_BC_XMAX, W_BC_XMAX, TH_BC_XMAX(1:N_TH), &
		U_BC_YMIN, V_BC_YMIN, W_BC_YMIN, TH_BC_YMIN(1:N_TH), &
		U_BC_YMAX, V_BC_YMAX, W_BC_YMAX, TH_BC_YMAX(1:N_TH), &
		U_BC_ZMIN, V_BC_ZMIN, W_BC_ZMIN, TH_BC_ZMIN(1:N_TH), &
		U_BC_ZMAX, V_BC_ZMAX, W_BC_ZMAX, TH_BC_ZMAX(1:N_TH)

    INTEGER :: PREVIOUS_TIME_STEP,SIGNAL_BC,UPDATE_DT

    LOGICAL :: VARIABLE_DT,FIRST_TIME, INT_TREAT, WAVE_ABS, STOCHASTIC_FORCING
    LOGICAL :: MOVIE,CREATE_NEW_FLOW,WRITE_VEL_TECPLOT

    REAL(r8) :: TIME, START_TIME,END_TIME, INT_TIME, INT_THX
    INTEGER  :: TIME_STEP, RK_STEP
    INTEGER  :: TIME_ARRAY(8)
    
    REAL(r8) :: RI_TAU(1:N_TH), PR(1:N_TH)
    LOGICAL  :: CREATE_NEW_TH(1:N_TH), BACKGROUND_GRAD(1:N_TH)
    INTEGER  :: NUM_READ_TH
    INTEGER  :: READ_TH_INDEX(1:N_TH)
    LOGICAL  :: FILTER_TH(1:N_TH)
    INTEGER  :: FILTER_INT(1:N_TH)
    INTEGER  :: JSTART_TH(1:N_TH),JEND_TH(1:N_TH),ZSTART_TH(1:N_TH), &
		ZEND_TH(1:N_TH) 
    REAL(r8) :: OMEGA0, AMP_OMEGA0, ANG_BETA,In_H0,Q_H0,H0,f_0

    REAL(r8),ALLOCATABLE,DIMENSION(:,:) :: U_BC_UPPER_NWM, W_BC_UPPER_NWM, & 
					    W_BC_LOWER_WALLMODEL, U_BC_LOWER_WALLMODEL, V_BC_LOWER_WALLMODEL
    
    REAL(r8),ALLOCATABLE,DIMENSION(:) ::     U1_bar
    REAL(r8),ALLOCATABLE,DIMENSION(:,:) ::   U2_bar, &
					      U3_bar, TH_BAR 
    
    COMPLEX(r8),ALLOCATABLE,DIMENSION(:,:) :: temp_mean 

    REAL(r8) :: UTAU_MEAN_LOWER,UTAU_MEAN_UPPER,TAUWALL_MEAN
    REAL(r8) :: UTAU_AVE
    
    LOGICAL :: FILTER_VEL
    INTEGER :: FILTER_VEL_INT

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|! Parameters for Large Eddy Simulation
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    LOGICAL :: LES
    INTEGER :: LES_MODEL_TYPE,LES_MODEL_TYPE_TH, WALL_MODEL_TYPE
    REAL(r8),ALLOCATABLE,DIMENSION(:,:,:) :: NU_T
    REAL(r8),ALLOCATABLE,DIMENSION(:,:,:,:) :: KAPPA_T
    INTEGER ::  J1,J2

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! RKW3 parameters
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    REAL(r8)  H_BAR(3), BETA_BAR(3), ZETA_BAR(3)

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Global variables
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  !     ARRAYS FOR WHEN DATA IS SPLITED IN X DIRECTION

  !      REAL(r8),target ::  U1 (0:NX+1,0:NZ+1,0:NY+1)
  !      REAL(r8),target ::  U2 (0:NX+1,0:NZ+1,0:NY+1)
  !      REAL(r8),target ::  U3 (0:NX+1,0:NZ+1,0:NY+1)
  !      REAL(r8),target ::  P  (0:NX+1,0:NZ+1,0:NY+1)
  !      REAL(r8),target ::  TH (0:NX+1,0:NZ+1,0:NY+1,1:N_TH)

  !      REAL(r8)        ::  U1b(0:NX+1,0:NZ+1,0:NY+1)
  !      REAL(r8)        ::  U2b(0:NX+1,0:NZ+1,0:NY+1)
  !      REAL(r8)        ::  U3b(0:NX+1,0:NZ+1,0:NY+1)

  !     ARRAYS FOR WHEN DATA IS SPLITED IN X DIRECTION

  !      complex(r8),target :: CU1(0:NX/2,0:NZ+1,0:NY+1)
  !      complex(r8),target :: CU2(0:NX/2,0:NZ+1,0:NY+1)
  !      complex(r8),target :: CU3(0:NX/2,0:NZ+1,0:NY+1)
  !      complex(r8),target :: CP (0:NX/2,0:NZ+1,0:NY+1)
  !      complex(r8),target :: CTH(0:NX/2,0:NZ+1,0:NY+1,1:N_TH)

  !      complex(r8)        ::  CU1b(0:NX/2,0:NZ+1,0:NY+1)
  !      complex(r8)        ::  CU2b(0:NX/2,0:NZ+1,0:NY+1) 
  !      complex(r8)        ::  CU3b(0:NX/2,0:NZ+1,0:NY+1) 

  !      EQUIVALENCE (U1,CU1), (U2,CU2), (U3,CU3), &
  !        (U1b,CU1b), (U2b,CU2b), (U3b,CU3b),     &
  !        (TH, CTH), (P,CP) 


  !     REAL(r8) R1 (0:NX+1,0:NZ+1,0:NY+1), R2 (0:NX+1,0:NZ+1,0:NY+1),    &
  !             R3 (0:NX+1,0:NZ+1,0:NY+1), F1 (0:NX+1,0:NZ+1,0:NY+1),     &
  !             F2 (0:NX+1,0:NZ+1,0:NY+1), F3 (0:NX+1,0:NZ+1,0:NY+1),     &
  !             S1 (0:NX+1,0:NZ+1,0:NY+1), S2 (0:NX+1,0:NZ+1,0:NY+1)

  !     REAL(r8)FTH (0:NX+1,0:NZ+1,0:NY+1,1:N_TH),                      &
  !             RTH (0:NX+1,0:NZ+1,0:NY+1,1:N_TH)



  !     complex(r8) CR1(0:NX/2,0:NZ+1,0:NY+1), CR2(0:NX/2,0:NZ+1,0:NY+1),  &
  !                 CR3(0:NX/2,0:NZ+1,0:NY+1), CF1(0:NX/2,0:NZ+1,0:NY+1),  &
  !                 CF2(0:NX/2,0:NZ+1,0:NY+1), CF3(0:NX/2,0:NZ+1,0:NY+1),  &
  !                 CS1(0:NX/2,0:NZ+1,0:NY+1), CS2(0:NX/2,0:NZ+1,0:NY+1)

  !     complex(r8) CFTH(0:NX/2,0:NZ+1,0:NY+1,1:N_TH), &
  !                CRTH(0:NX/2,0:NZ+1,0:NY+1,1:N_TH) 


  !     EQUIVALENCE (R1,CR1), (R2,CR2) &
  !        , (R3,CR3) , (F1,CF1), (F2,CF2), (F3,CF3), (S1,CS1) &
  !        , (S2,CS2),  (RTH,CRTH) , (FTH, CFTH)


  !      REAL(r8),DIMENSION(:,:,:),pointer      :: u_p
  !      complex(r8),DIMENSION(:,:,:),pointer    :: cu_p

  !       EQUIVALENCE (temp,ctemp)

    REAL(r8),ALLOCATABLE,DIMENSION(:,:) :: TH_BACK

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!       TRAIL VARIABLE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    REAL(r8),DIMENSION(0:NXP,0:NZV-1,0:NY+1)    	::  S1X,S2X
    COMPLEX(r8),DIMENSION(0:NX2P,0:NZV-1,0:NY+1)	::  CS1X, CS2X
    
    REAL(r8),DIMENSION(0:NXV-1,0:NZP,0:NY+1)    	::  S1Z,S2Z 
    COMPLEX(r8),DIMENSION(0:NX2V-1,0:NZP,0:NY+1)	::  CS1Z,CS2Z
    
    REAL(r8)    ::  VARP (0:NX+1,0:NZP,0:NY+1)
    COMPLEX(r8) ::  CVARP(0:NX/2,0:NZP,0:NY+1)

    EQUIVALENCE (S1X,CS1X,S1Z,CS1Z), (S2X,CS2X,S2Z,CS2Z)

    EQUIVALENCE (VARP,CVARP)

    REAL(r8),ALLOCATABLE,DIMENSION(:,:,:), target		::  U1X,U2X,U3X,PX,U2bX,U3bX,R1X,R2X,R3X,F1X,F2X,F3X
    COMPLEX(r8),ALLOCATABLE,DIMENSION(:,:,:), target		::  CU1X,CU2X,CU3X,CPX,CU2bX,CU3bX,CR1X,CR2X,CR3X,CF1X,CF2X,CF3X,CU1X_WALL_MODEL,CU2X_WALL_MODEL,CU3X_WALL_MODEL

    REAL(r8),ALLOCATABLE,DIMENSION(:,:,:,:), target 		::  THX,RTHX,FTHX
    COMPLEX(r8),ALLOCATABLE,DIMENSION(:,:,:,:), target 	::  CTHX,CRTHX,CFTHX, CTHX_WALL_MODEL

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! ADI variables
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

  !     real(r8),ALLOCATABLE,DIMENSION(:,:) ::  MATL, MATD, &
  !               MATU, VEC
  !     real(r8),ALLOCATABLE,DIMENSION(:,:) ::  MATL_Z, MATD_Z, MATU_Z, VEC_Z

  !      REAL(r8)     MATL (0:NX-1,0:NY+1), MATD(0:NX-1,0:NY+1), &
  !          MATU(0:NX-1,0:NY+1), VEC(0:NX-1,0:NY+1)
  !      REAL(r8)     MATL_Z (0:NX-1,0:NZ+1), MATD_Z(0:NX-1,0:NZ+1), &
  !          MATU_Z(0:NX-1,0:NZ+1), VEC_Z(0:NX-1,0:NZ+1)

    REAL*8 TIME_OLD,TIME_OLD_ENG
    
    INTEGER NSAMPLES

    REAL(r8),ALLOCATABLE,DIMENSION(:,:) :: INT_JACOB
    REAL(r8),ALLOCATABLE,DIMENSION(:,:,:) :: GMAT_11,&
	GMAT_12, GMAT_22, &
	CJOB_11, CJOB_12, & 
	CJOB_21, CJOB_22

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! openMP variables - ALL ARRAYS
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    INTEGER kstart, kend, chunk,nthreads,iam,& !omp_get_num_threads, &
	    k_start, k_end, j_start, j_end, jj_start, jj_end
    REAL(r8)  wtime, rktime !, omp_get_wtime
    common /bounds/ kstart, kend, k_start, k_end, j_start, j_end, &
    jj_start, jj_end

  !      common /THOMAS_TH/ MATL, MATD, MATU, VEC
  !!	$omp threadprivate(/bounds/)
END MODULE run_variable

MODULE ADI_var

  USE ntypes, ONLY : r8
  USE Domain

!      REAL(r8)     MATL (0:NX-1,0:NY+1), MATD(0:NX-1,0:NY+1), &
!          MATU(0:NX-1,0:NY+1), VEC(0:NX-1,0:NY+1)
!      REAL(r8)     MATL_Z (0:NX-1,0:NZ+1), MATD_Z(0:NX-1,0:NZ+1), &
!          MATU_Z(0:NX-1,0:NZ+1), VEC_Z(0:NX-1,0:NZ+1)
!     common /THOMAS_TH_Z/ MATL_Z, MATD_Z, MATU_Z, VEC_Z
!     common /THOMAS_TH/ MATL, MATD, MATU, VEC
!!$omp threadprivate(/THOMAS_TH_Z/, /THOMAS_TH/)
    REAL(r8),ALLOCATABLE,DIMENSION(:,:) :: MATLX,MATDX,MATUX,VECX
    REAL(r8),ALLOCATABLE,DIMENSION(:,:) :: MATLY,MATDY,MATUY,VECY
     
END MODULE ADI_var

MODULE variable_stat

  USE ntypes, 	ONLY : r8
  USE Domain,	ONLY : N_TH
  ! Variables for outputting statistics
    REAL(r8),ALLOCATABLE,DIMENSION(:) :: UBAR,VBAR,WBAR
    REAL(r8),ALLOCATABLE,DIMENSION(:,:) :: URMS,VRMS, &
	  WRMS, UV,UW, WV,PV,PU, &
	  DWDY, DUDY,DWDZ, DUDZ, &
	  SHEAR, OMEGA_X,OMEGA_Y,OMEGA_Z,TKE
	
    REAL(r8) ::     URMS_B,VRMS_B,WRMS_B,TKE_B

  ! Variables needed for SAVE_STATS_TH
    
    REAL(r8),ALLOCATABLE,DIMENSION(:,:) ::   THBAR, PE_DISS, Rig

    REAL(r8) THRMS_B(1:N_TH)

    REAL(r8),ALLOCATABLE,DIMENSION(:,:,:) :: THRMS, &
	      THV, THW, DTHDY, DTHDZ

  ! Variables for tkebudget
    REAL(r8),ALLOCATABLE,DIMENSION(:,:) :: epsilon, &
	    tke_mean,tke_mean_old,energy_mean,energy_mean_old


    REAL(r8),ALLOCATABLE,DIMENSION(:,:) :: tke_1,tke_2, &
	tke_2_1, tke_2_2, tke_3, tke_3_1, tke_3_2, tke_3_3, &
	tke_4, tke_6_1,tke_6_1_1, tke_6_1_2, tke_6_1_3, &
	tke_6_2, tke_6_2_1, tke_6_2_2, tke_6_2_3, tke_7, &
	S1_mean, p_mean, transport


    REAL(r8),ALLOCATABLE,DIMENSION(:,:,:) :: tke_5
     
    REAL(r8),ALLOCATABLE,DIMENSION(:,:) :: UU, VV, WW, UUU, VVV, WWW, UUUU, VVVV, WWWW, UUV, WWV, UUUV, WWWV
    REAL(r8),ALLOCATABLE,DIMENSION(:,:) :: dWdz1, dWdz2, dWdz3, dWdz4, dWdz5 
END MODULE variable_stat


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Multigrid variables - ALL ARRAYS 1D
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
MODULE mg_vari

  USE ntypes
      
    INTEGER ifail, istart, maxit, iprep, ndid, nout
    INTEGER MUDPACK_BC(4), iout(6)
    LOGICAL INIT_FLAG 

    REAL(r8),ALLOCATABLE,DIMENSION(:,:) :: V, VC, &
	    VB,VBC, RHS, RHSC, A, &
	    AC, LDU, LDUC, WORK, WA, WAC, WB, WBC

    REAL*8 TOL, RESNO

END 

MODULE forcing
  USE Domain
  USE ntypes
   
    REAL(r8),ALLOCATABLE,DIMENSION(:,:) :: f, g
    REAL(r8),ALLOCATABLE,DIMENSION(:,:) :: f_prime, g_prime
    COMPLEX(r8),ALLOCATABLE, DIMENSION(:) ::  u_prime,v_prime, beam_th
    COMPLEX*8 im

! real(r8) x0,y0,x_local,y_local, D, C

END 

MODULE pr_rem

  USE ntypes
  USE Domain
  
    REAL*8 P_TMP(0:NZ+1,0:NY+1), P_TMP2(0:NZ+1,0:NY+1), &
	  RHS_TMP(0:NZ+1,0:NY+1)

	  common /mg_rh/ P_TMP, P_TMP2,RHS_TMP
  !$omp threadprivate(/mg_rh/)

END MODULE pr_rem

MODULE IO

 USE ntypes, ONLY: i4, r8
 
  INTEGER(i4)             :: I_OUT,IOUT_MASTER,IOUT_SLAVE !Unit to write output 6 is screen
  CHARACTER(len=100)      :: resultDIR,tempDIR,runDIR,ext,penDIR,plnDIR,statDIR,flowDIR, gridDIR,MGDIR,relaxDIR
  
END MODULE IO


MODULE les_chan_var

  USE ntypes
  USE Domain

  ! Variables for dynamic Smagrinsky
    REAL(r8),ALLOCATABLE,DIMENSION(:,:,:)	:: numerator,denominator
    REAL(r8),ALLOCATABLE,DIMENSION(:,:)	:: denominator_sum,numerator_sum

    REAL(r8),ALLOCATABLE,DIMENSION(:,:)	:: NU_T_mean, DELTA_Y,U1_bar_les, U3_bar_les,U2_bar_les,C_DYN, KAPPA_T_MEAN
    REAL(r8),ALLOCATABLE,DIMENSION(:,:,:) 	:: C_DYN_TH
    REAL(r8),ALLOCATABLE,DIMENSION(:,:,:) 	:: Sij_mean, TAU_mean
    REAL(r8),ALLOCATABLE,DIMENSION(:,:) 	:: tke_sgs_diss,tke_sgs_p,tke_sgs_t
    REAL(r8),ALLOCATABLE,DIMENSION(:,:,:,:)	:: U_BAR_TIL,U_2BAR,U_4BAR
    REAL(r8),ALLOCATABLE,DIMENSION(:,:,:)   	:: S_2BAR
    REAL(r8),ALLOCATABLE,DIMENSION(:,:,:,:)	:: St_rij 
    REAL(r8),ALLOCATABLE,DIMENSION(:,:,:)	:: GMAT_11_z,GMAT_11_y,GMAT_12_z,GMAT_12_y,GMAT_22_z,GMAT_22_y


    INTEGER J1i, J2e
    

    REAL*8 C_DYN_H(0:NY+1),C_DYN_V(0:NY+1)
  ! Variables for plane-averaged momentum budget
    REAL*8 NU_U1(0:NY+1)
    REAL*8 NU_U3(0:NY+1)
    REAL*8 DELTA_YF(0:NY+1)

  ! For the TKE part
    REAL*8 tke_sgs_mm(0:NY+1),tke_sgs_evm(0:NY+1)

    REAL*8 TEMP(0:NXP,0:NZV-1,0:NY+1),Mij(0:NXP,0:NZV-1,0:NY+1)
    REAL*8 TEMP_1(0:NXP,0:NZV-1,0:NY+1)
    REAL*8 TEMP_2(0:NXP,0:NZV-1,0:NY+1)
    REAL*8 TEMP_3(0:NXP,0:NZV-1,0:NY+1)
    REAL*8  Sij(0:NXP,0:NZV-1,0:NY+1,1:6)

    REAL*8 cross

    COMPLEX*16 CSij(0:NX2P,0:NZV-1,0:NY+1,1:6)
    COMPLEX*16 CTEMP(0:NX2P,0:NZV-1,0:NY+1)
    COMPLEX*16 CTEMP_1(0:NX2P,0:NZV-1,0:NY+1)
    COMPLEX*16 CTEMP_2(0:NX2P,0:NZV-1,0:NY+1)
    COMPLEX*16 CTEMP_3(0:NX2P,0:NZV-1,0:NY+1)
    COMPLEX*16 CMij(0:NX2P,0:NZV-1,0:NY+1)

    EQUIVALENCE (TEMP,CTEMP)     &
	    ,(Mij,CMij)       &
	    ,(TEMP_1,CTEMP_1) &
	    ,(TEMP_2,CTEMP_2) &
	    ,(TEMP_3,CTEMP_3) 

END


MODULE mpi_var
  
  USE ntypes
  USE Domain
  
    INTEGER :: MPI_COMM_ROW, MPI_COMM_COL,MPI_COMM_CART
    INTEGER :: NPROCES, RANK, RANK_ROW, RANK_COL
    INTEGER :: status(MPI_STATUS_SIZE), IERROR


    INTEGER(i4),PARAMETER                   :: realtype  =MPI_DOUBLE_PRECISION
    INTEGER(i4),PARAMETER                   :: inttype   =MPI_INTEGER
    INTEGER(i4),PARAMETER                   :: chartype  =MPI_CHARACTER
    INTEGER(i4),PARAMETER                   :: logictype =MPI_LOGICAL
    INTEGER(i4),PARAMETER                   :: cmplxtype =MPI_DOUBLE_COMPLEX
    INTEGER(i4),PARAMETER                   :: commx1x2x3=MPI_COMM_WORLD

    ! complex(r8),ALLOCATABLE,DIMENSION(:,:,:)    :: IN_CZ,IN_CX
    ! complex(r8),ALLOCATABLE,DIMENSION(:,:,:)    :: OUT_CZ,OUT_CX


    REAL(r8),ALLOCATABLE,DIMENSION(:)        :: M_IN,M_OUT
    COMPLEX(r8),ALLOCATABLE,DIMENSION(:)     :: M_IN_C,M_OUT_C
    
END
