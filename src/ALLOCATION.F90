subroutine allocate_var

USE ntypes
USE Domain
USE Grid
USE Fft_var
USE TIME_STEP_VAR
USE run_variable
USE mg_vari
USE ADI_var
USE mpi_var
USE IO,     ONLY: I_OUT


implicit none

!Passed Variables
 INTEGER        :: ok

!Local Variables
 INTEGER           :: s_1

 !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 ! grid allocation
 ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

ALLOCATE (GX(0:NXP+1), stat=s_1) 
ALLOCATE (GY(0:NY+1), stat=s_1)
ALLOCATE (GZ(0:NZ+1), stat=s_1)
ALLOCATE (DX(0:NXP+1), stat=s_1)
ALLOCATE (DY(0:NY+1), stat=s_1)
ALLOCATE (DZ(0:NZ+1), stat=s_1)
ALLOCATE (GXF(0:NX), stat=s_1)
ALLOCATE (GYF(0:NY+1), stat=s_1) 
ALLOCATE (GZF(0:NZ), stat=s_1)
ALLOCATE (DXF(0:NX), stat=s_1)
ALLOCATE (DYF(0:NY), stat=s_1)
ALLOCATE (DZF(0:NZ), stat=s_1)
ALLOCATE (xpoint(0:NZ+1,0:NY+1), ypoint(0:NZ+1,0:NY+1))

  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
IF (s_1.NE.0) THEN
  write(I_OUT,*) "Error Allocating Grid Variables" 
  GOTO 1000
ENDIF


  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! fft allocation
  ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
ALLOCATE  ( KX (0:NX/3), stat=s_1)
allocate  ( KY(0:2*(NY/3)), stat=s_1)
ALLOCATE  ( KZ  (0:2*(NZ/3)), stat=s_1)
ALLOCATE  ( KX2 (0:NX/3), stat=s_1) 
ALLOCATE  ( KY2 (0:2*(NY/3)), stat=s_1)
ALLOCATE  ( KZ2 (0:2*(NZ/3)), stat=s_1)
ALLOCATE  (CIKX(0:NX/3), stat=s_1)
ALLOCATE  (CIKY(0:2*(NY/3)), stat=s_1)
ALLOCATE  ( CIKZ(0:2*(NZ/3)), stat=s_1)
ALLOCATE  ( CZX_PLANE(0:NZ,0:NX2P), stat=s_1)
ALLOCATE  ( CYZ_PLANE(0:NY,0:2*(NZ/3)), stat=s_1)

ALLOCATE  ( KXP(0:NX2P), stat=s_1)
ALLOCATE  ( KX2P(0:NX2P), stat=s_1)
ALLOCATE  (CIKXP(0:NX2P), stat=s_1)

IF (s_1.NE.0) THEN
  write(I_OUT,*) "Error Allocating FFT Variables" 
  GOTO 1000
ENDIF

! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! Input parameters and runtime variables
  ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

ALLOCATE (SPONGE_SIGMA(0:NY+1), stat=s_1)
ALLOCATE (SPONGE_SIGMA_OUT(0:NZ+1), stat=s_1)
ALLOCATE (C_int(0:NY+1),C_int_le(0:NY+1), stat=s_1)

ALLOCATE (U_BC_UPPER_NWM(0:NXP,0:NZ+1), stat=s_1)
ALLOCATE (W_BC_UPPER_NWM(0:NXP,0:NZ+1), stat=s_1)

ALLOCATE (W_BC_LOWER_WALLMODEL(0:NXP,0:NZ+1), stat=s_1)
ALLOCATE (V_BC_LOWER_WALLMODEL(0:NXP,0:NZ+1), stat=s_1)
ALLOCATE (U_BC_LOWER_WALLMODEL(0:NXP,0:NZ+1), stat=s_1)

ALLOCATE (U1_bar(0:NY+2), stat=s_1)
ALLOCATE (U2_bar(0:NZV-1,0:NY+1), stat=s_1)
ALLOCATE (U3_bar(0:NZV-1,0:NY+1), stat=s_1)
ALLOCATE (TH_BAR(0:NZ+2,0:NY+2), stat=s_1)
ALLOCATE (temp_mean(0:NZV-1,0:NY+1), stat=s_1) 


IF (s_1.NE.0) THEN
  write(I_OUT,*) "Error Allocating VEL_BAR Variables" 
  GOTO 1000
ENDIF

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! ! Global variables
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! ALLOCATE      (U1 (0:NX+1,0:NZ+1,0:NY+1), U2 (0:NX+1,0:NZ+1,0:NY+1), stat=s_1)
! ALLOCATE      (U3 (0:NX+1,0:NZ+1,0:NY+1), P  (0:NX+1,0:NZ+1,0:NY+1), stat=s_1)
! ALLOCATE      (R1 (0:NX+1,0:NZ+1,0:NY+1), R2 (0:NX+1,0:NZ+1,0:NY+1), stat=s_1)
! ALLOCATE      (R3 (0:NX+1,0:NZ+1,0:NY+1), F1 (0:NX+1,0:NZ+1,0:NY+1), stat=s_1)
! ALLOCATE      (F2 (0:NX+1,0:NZ+1,0:NY+1), F3 (0:NX+1,0:NZ+1,0:NY+1), stat=s_1)
! ALLOCATE      (S1 (0:NX+1,0:NZ+1,0:NY+1), U1b(0:NX+1,0:NZ+1,0:NY+1), stat=s_1)
! ALLOCATE      (U2b(0:NX+1,0:NZ+1,0:NY+1), U3b(0:NX+1,0:NZ+1,0:NY+1), stat=s_1)
! ALLOCATE      (S2 (0:NX+1,0:NZ+1,0:NY+1)) 
! ALLOCATE      (TH (0:NX+1,0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
! ALLOCATE      (TH_BACK(0:NY+1,1:N_TH), stat=s_1)
! ALLOCATE      (FTH (0:NX+1,0:NZ+1,0:NY+1,1:N_TH), stat=s_1) 
! ALLOCATE      (RTH (0:NX+1,0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
! IF (s_1.NE.0) THEN
!   write(I_OUT,*) "Error Allocating VEL_GOL_REL Variables" 
!   GOTO 1000
! ENDIF
! 
! ALLOCATE     (CU1(0:NX/2,0:NZ+1,0:NY+1), CU2(0:NX/2,0:NZ+1,0:NY+1), stat=s_1) 
! ALLOCATE     (CU3(0:NX/2,0:NZ+1,0:NY+1), CP (0:NX/2,0:NZ+1,0:NY+1), stat=s_1)
! ALLOCATE     (CR1(0:NX/2,0:NZ+1,0:NY+1), CR2(0:NX/2,0:NZ+1,0:NY+1), stat=s_1)
! ALLOCATE     (CR3(0:NX/2,0:NZ+1,0:NY+1), CF1(0:NX/2,0:NZ+1,0:NY+1), stat=s_1)
! ALLOCATE     (CF2(0:NX/2,0:NZ+1,0:NY+1), CF3(0:NX/2,0:NZ+1,0:NY+1), stat=s_1)
! ALLOCATE     (CS1(0:NX/2,0:NZ+1,0:NY+1), CU1b(0:NX/2,0:NZ+1,0:NY+1), stat=s_1)
! ALLOCATE     (CU2b(0:NX/2,0:NZ+1,0:NY+1), CU3b(0:NX/2,0:NZ+1,0:NY+1), stat=s_1) 
! ALLOCATE     (CS2(0:NX/2,0:NZ+1,0:NY+1), stat=s_1)
! ALLOCATE     (CTH(0:NX/2,0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
! ALLOCATE     (CFTH(0:NX/2,0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
! ALLOCATE     (CRTH(0:NX/2,0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
! 
! IF (s_1.NE.0) THEN
!   write(I_OUT,*) "Error Allocating VEL_GOL_COMPLX Variables" 
!   GOTO 1000
! ENDIF

ALLOCATE      (TH_BACK(0:NY+1,1:N_TH), stat=s_1)
! EQUIVALENCE (U1,CU1), (U2,CU2), (U3,CU3), (R1,CR1), (R2,CR2) &
!         , (R3,CR3) , (F1,CF1), (F2,CF2), (F3,CF3), (P,CP), (S1,CS1) &
!         , (S2,CS2), (U1b,CU1b), (U2b,CU2b), (U3b,CU3b), (RTH,CRTH)  &
!         , (TH, CTH), (FTH, CFTH)

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! ! ADI variables
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!ALLOCATE     (MATL (0:NX-1,0:NY+1), MATD(0:NX-1,0:NY+1), stat=s_1)
!ALLOCATE     (MATU(0:NX-1,0:NY+1), VEC(0:NX-1,0:NY+1), stat=s_1)
!ALLOCATE     (MATL_Z (0:NX-1,0:NZ+1), MATD_Z(0:NX-1,0:NZ+1), stat=s_1)
!ALLOCATE     (MATU_Z(0:NX-1,0:NZ+1), VEC_Z(0:NX-1,0:NZ+1), stat=s_1)

ALLOCATE     ( MATLX(0:NXP,0:NZ+1),MATDX(0:NXP,0:NZ+1), stat=s_1)
ALLOCATE     ( MATUX(0:NXP,0:NZ+1),VECX(0:NXP,0:NZ+1), stat=s_1)

ALLOCATE     ( MATLY(0:NXP,0:NY+1),MATDY(0:NXP,0:NY+1), stat=s_1)
ALLOCATE     ( MATUY(0:NXP,0:NY+1),VECY(0:NXP,0:NY+1), stat=s_1)

IF (s_1.NE.0) THEN
  write(I_OUT,*) "Error Allocating ADI Variables" 
  GOTO 1000
ENDIF

! JACOBIAN VAR
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
ALLOCATE     (INT_JACOB(0:NZ+2,0:NY+2), stat=s_1)
ALLOCATE     (GMAT_11(0:NZ+2,0:NY+2,2))
ALLOCATE     (GMAT_12(0:NZ+2,0:NY+2,2), GMAT_22(0:NZ+2,0:NY+2,2), stat=s_1)
ALLOCATE     (CJOB_11(0:NZ+2,0:NY+2,2), CJOB_12(0:NZ+2,0:NY+2,2), stat=s_1)
ALLOCATE     (CJOB_21(0:NZ+2,0:NY+2,2), CJOB_22(0:NZ+2,0:NY+2,2), stat=s_1)

IF (s_1.NE.0) THEN
  write(I_OUT,*) "Error Allocating JACOBIAN Variables for theta" 
  GOTO 1000
ENDIF
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Multigrid variables - ALL ARRAYS 1D
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
ALLOCATE      (V(1:(NY+2)*(NZ+2),0:NX2P), stat=s_1 )
ALLOCATE      (VC(1:NM,0:NX2P), stat=s_1)
ALLOCATE      (VB(1:(NY+2)*(NZ+2),0:NX2P), stat=s_1) 
ALLOCATE      (VBC(1:NM,0:NX2P), stat=s_1)
ALLOCATE      (RHS(1:(NY+2)*(NZ+2),0:NX2P), stat=s_1)
ALLOCATE      (RHSC(1:NM,0:NX2P), stat=s_1)
ALLOCATE      (A(1:(NY+2)*(NZ+2)*9,0:NX2P), stat=s_1)
ALLOCATE      (AC(1:9*NM,0:NX2P), stat=s_1)
ALLOCATE      (LDU(1:3*(NY+2)*(NZ+2),0:NX2P), stat=s_1)
ALLOCATE      (LDUC(1:3*NM,0:NX2P), stat=s_1)
ALLOCATE      (WORK(1:(NZ+2)*12,0:NX2P), stat=s_1)
ALLOCATE      (WA(1:(NY+2)*(NZ+2),0:NX2P), stat=s_1)
ALLOCATE      (WAC(NM,0:NX2P), stat=s_1)
ALLOCATE      (WB(1:(NY+2)*(NZ+2),0:NX2P), stat=s_1)
ALLOCATE      (WBC(1:NM,0:NX2P), stat=s_1)

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! MPI variables - ALL ARRAYS 1D
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!ALLOCATE       (OUT_CX(0:NX2P,0:NZV-1,0:NY+1),IN_CZ(0:NX2V-1,0:NZP,0:NY+1), stat=s_1)
!ALLOCATE       (IN_CX(0:NX2P,0:NZV-1,0:NY+1),OUT_CZ(0:NX2V-1,0:NZP,0:NY+1), stat=s_1)

ALLOCATE      (M_IN_C(1:NX2V*(NY+2)*NZV/NP),M_OUT_C(1:NX2V*(NY+2)*NZV/NP), stat=s_1)
ALLOCATE      (M_IN(1:NXV*(NY+2)*NZV/NP), M_OUT(1:NXV*(NY+2)*NZV/NP), stat=s_1)


IF (s_1.NE.0) THEN
  write(I_OUT,*) "Error Allocating MG Variables for theta" 
  GOTO 1000
ENDIF

write(I_OUT,'(a)') "ALLOCATION COMPLETED"
return
1000 continue 
 ok = s_1
 write(I_OUT,'(a)') "ALLOCATION FAILED"

return
END

subroutine allocation_u (FLAG)

  USE ntypes
  USE Domain
  USE run_variable
  USE IO,     ONLY: I_OUT

    IMPLICIT NONE
    
    LOGICAL FLAG
    INTEGER           :: s_1

    s_1 = 0

    IF (FLAG) THEN

      IF (.not. allocated(U1X) ) ALLOCATE( U1X(0:NXP,0:NZV-1,0:NY+1) )
      IF ( allocated(CU1X) ) THEN
	DEALLOCATE( CU1X, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating cu1"
      ENDIF

    ELSE
      IF (.not. allocated(CU1X) ) ALLOCATE( CU1X(0:NX2P,0:NZV-1,0:NY+1))
      IF ( allocated(U1X) ) THEN
	DEALLOCATE(U1X, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating u1"
      ENDIF

    ENDIF


  RETURN
END

subroutine allocation_v (FLAG)
  USE ntypes
  USE Domain
  USE run_variable
  USE IO,     ONLY: I_OUT

    IMPLICIT NONE
    
    LOGICAL FLAG
    INTEGER           :: s_1
    
    s_1 = 0

    IF (FLAG) THEN

      IF (.not. allocated(U2X) ) ALLOCATE( U2X(0:NXP,0:NZV-1,0:NY+1) )
      IF ( allocated(CU2X) ) THEN
	DEALLOCATE( CU2X, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating cu1"
      ENDIF

    ELSE

      IF (.not. allocated(CU2X) ) ALLOCATE( CU2X(0:NX2P,0:NZV-1,0:NY+1))
      IF ( allocated(U2X) ) THEN
	DEALLOCATE(U2X, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating u2"
      ENDIF

    ENDIF

  RETURN
END

subroutine allocation_w (FLAG)

  USE ntypes
  USE Domain
  USE run_variable
  USE IO,     ONLY: I_OUT

    IMPLICIT NONE
    
    LOGICAL FLAG
    INTEGER           :: s_1

    s_1 = 0
    
    IF (FLAG) THEN

      IF (.not. allocated(U3X) ) ALLOCATE( U3X(0:NXP,0:NZV-1,0:NY+1) )
      IF ( allocated(CU3X) ) THEN
	DEALLOCATE( CU3X, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating cu1"
      ENDIF

    ELSE

      IF (.not. allocated(CU3X) ) ALLOCATE( CU3X(0:NX2P,0:NZV-1,0:NY+1))
      IF ( allocated(U3X) ) THEN
	DEALLOCATE(U3X, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating u3"
      ENDIF

    ENDIF


RETURN
END

subroutine allocation_p (FLAG)

  USE ntypes
  USE Domain
  USE run_variable
  USE IO,     ONLY: I_OUT

  IMPLICIT NONE
  
    LOGICAL FLAG
    INTEGER           :: s_1

    s_1 = 0
      
    IF (FLAG) THEN

      IF (.not. allocated(PX) ) ALLOCATE( PX(0:NXP,0:NZV-1,0:NY+1) )
      IF ( allocated(CPX) ) THEN
	DEALLOCATE( CPX, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating CPX"
      ENDIF

    ELSE

      IF (.not. allocated(CPX) ) ALLOCATE( CPX(0:NX2P,0:NZV-1,0:NY+1))
      IF ( allocated(PX) ) THEN
	DEALLOCATE(PX, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating PX"
      ENDIF

      ENDIF

  RETURN
END

subroutine allocation_th (FLAG)
  USE ntypes
  USE Domain
  USE run_variable
  USE IO,     ONLY: I_OUT

    IMPLICIT NONE
    
    LOGICAL FLAG
    INTEGER           :: s_1
    
    s_1 = 0

    IF (FLAG) THEN
      
      IF (.not. allocated(THX) ) ALLOCATE( THX(0:NXP,0:NZV-1,0:NY+1,1:N_TH) ) 
      IF ( allocated(CTHX ) ) THEN
	DEALLOCATE( CTHX, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating CTHX"
      ENDIF
      
    ELSE

      IF (.not. allocated(CTHX ) ) ALLOCATE( CTHX(0:NX2P,0:NZV-1,0:NY+1,1:N_TH) ) 
      IF ( allocated(THX) ) THEN
	DEALLOCATE(THX, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating THX"
      ENDIF
    
    ENDIF

  return
END


subroutine allocation_ub (FLAG)

  USE ntypes
  USE Domain
  USE run_variable
  USE IO,     ONLY: I_OUT

    IMPLICIT NONE
    
    LOGICAL FLAG
    INTEGER           :: s_1

    s_1 = 0

    IF (FLAG) THEN

    IF (allocated(CU2bx) ) THEN
      DEALLOCATE( CU2bX )
      DEALLOCATE( CU3bX )
    ENDIF

    IF (.not. allocated(U2bX) ) THEN
      ALLOCATE( U2bX(0:NXP,0:NZV-1,0:NY+1) )
      allocate( U3bX(0:NXP,0:NZV-1,0:NY+1) )
    ENDIF

    ELSE

    IF (allocated(U2bX) ) THEN
    DEALLOCATE( U2bX )
    DEALLOCATE( U3bX )
    ENDIF

    IF (.not. allocated(CU2bX) ) THEN
    allocate( CU2bX(0:NX2P,0:NZV-1,0:NY+1) )
    allocate( CU3bX(0:NX2P,0:NZV-1,0:NY+1) )
    ENDIF

    ENDIF


  RETURN
END


subroutine allocation_R1(FLAG)

  USE ntypes
  USE Domain
  USE run_variable
  USE IO,     ONLY: I_OUT

    IMPLICIT NONE
    
    LOGICAL FLAG
    INTEGER           :: s_1
    
    s_1 = 0

    IF (FLAG) THEN

      IF (.not. allocated(R1X) ) ALLOCATE( R1X(0:NXP,0:NZV-1,0:NY+1) )
      IF ( allocated(CR1X) ) THEN
	DEALLOCATE( CR1X, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating CR1X"
      ENDIF

    ELSE

      IF (.not. allocated(CR1X) ) ALLOCATE( CR1X(0:NX2P,0:NZV-1,0:NY+1))
      IF ( allocated(R1X) ) THEN
	DEALLOCATE(R1X, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating R1"
      ENDIF

      ENDIF

  RETURN
END

subroutine allocation_R2(FLAG)
  USE ntypes
  USE Domain
  USE run_variable
  USE IO,     ONLY: I_OUT

    IMPLICIT NONE
    LOGICAL FLAG
    INTEGER           :: s_1

    s_1 = 0
    
    IF (FLAG) THEN

      IF (.not. allocated(R2X) ) ALLOCATE( R2X(0:NXP,0:NZV-1,0:NY+1) )
      IF ( allocated(CR2X) ) THEN
	DEALLOCATE( CR2X, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating CR2X"
      ENDIF

    ELSE

      IF (.not. allocated(CR2X) ) ALLOCATE( CR2X(0:NX2P,0:NZV-1,0:NY+1))
      IF (.not. allocated(CU1X_WALL_MODEL) .AND. (W_BC_YMIN.EQ.3) .AND.(WALL_MODEL_TYPE.EQ.2 .OR. WALL_MODEL_TYPE.EQ.4) ) ALLOCATE( CU1X_WALL_MODEL(0:NX2P,0:NZV-1,0:NY+1))
      IF (.not. allocated(CU2X_WALL_MODEL) .AND. (W_BC_YMIN.EQ.3) .AND. (WALL_MODEL_TYPE.EQ.2 .OR. WALL_MODEL_TYPE.EQ.4) ) ALLOCATE( CU2X_WALL_MODEL(0:NX2P,0:NZV-1,0:NY+1))
      IF (.not. allocated(CU3X_WALL_MODEL) .AND. (W_BC_YMIN.EQ.3) .AND. (WALL_MODEL_TYPE.EQ.2.OR. WALL_MODEL_TYPE.EQ.4) ) ALLOCATE( CU3X_WALL_MODEL(0:NX2P,0:NZV-1,0:NY+1))
      IF ( allocated(R2X) ) THEN
	DEALLOCATE(R2X, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating R2"
      ENDIF

      ENDIF

  RETURN
END

subroutine allocation_R3(FLAG)
  USE ntypes
  USE Domain
  USE run_variable
  USE IO,     ONLY: I_OUT
  
    IMPLICIT NONE
    
    LOGICAL FLAG
    INTEGER           :: s_1
    
    s_1 = 0
    
    IF (FLAG) THEN
    
      IF (.not. allocated(R3X) ) ALLOCATE( R3X(0:NXP,0:NZV-1,0:NY+1) )
      IF ( allocated(CR3X) ) THEN
	DEALLOCATE( CR3X, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating CR3X"
      ENDIF
    
    ELSE

      IF (.not. allocated(CR3X) ) ALLOCATE( CR3X(0:NX2P,0:NZV-1,0:NY+1))
      IF ( allocated(R3X) ) THEN
	DEALLOCATE(R3X, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating R3"
      ENDIF
    
    ENDIF

  RETURN
END

subroutine allocation_Rth (FLAG)
  USE ntypes
  USE Domain
  USE run_variable
  USE IO,     ONLY: I_OUT

    IMPLICIT NONE
    
    LOGICAL FLAG
    INTEGER           :: s_1
    
    s_1 = 0

    IF (FLAG) THEN
    
      IF (.not. allocated(RTHX) ) ALLOCATE( RTHX(0:NXP,0:NZV-1,0:NY+1,1:N_TH) ) 
      IF ( allocated(CRTHX ) ) THEN
	DEALLOCATE( CRTHX, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating CRTHX"
      ENDIF
    
    ELSE

      IF (.not. allocated(CRTHX ) ) ALLOCATE( CRTHX(0:NX2P,0:NZV-1,0:NY+1,1:N_TH) )
      IF (.NOT. allocated(CTHX_WALL_MODEL) .AND. (W_BC_YMIN.EQ.3) .AND. (WALL_MODEL_TYPE.EQ.2)) ALLOCATE( CTHX_WALL_MODEL(0:NX2P,0:NZV-1,0:NY+1,1:N_TH) ) 
      IF ( allocated(RTHX) ) THEN
	DEALLOCATE(RTHX, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating RTHX"
      ENDIF
    
    ENDIF

  RETURN
END


subroutine allocation_F1(FLAG)
  USE ntypes
  USE Domain
  USE run_variable
  USE IO,     ONLY: I_OUT

    IMPLICIT NONE
    LOGICAL FLAG
    INTEGER           :: s_1
    
    s_1 = 0

    IF (FLAG) THEN

      IF (.not. allocated(F1X) ) ALLOCATE( F1X(0:NXP,0:NZV-1,0:NY+1) )
      IF ( allocated(CF1X) ) THEN
	DEALLOCATE( CF1X, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating CF1X"
      ENDIF

    ELSE

      IF (.not. allocated(CF1X) ) ALLOCATE( CF1X(0:NX2P,0:NZV-1,0:NY+1))
      IF ( allocated(F1X) ) THEN
	DEALLOCATE(F1X, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating F1X"
      ENDIF

    ENDIF

  RETURN
END

subroutine allocation_F2(FLAG)
  USE ntypes
  USE Domain
  USE run_variable
  USE IO,     ONLY: I_OUT

    IMPLICIT NONE
    
    LOGICAL FLAG
    INTEGER           :: s_1

    s_1 = 0
    
    IF (FLAG) THEN

      IF (.not. allocated(F2X) ) ALLOCATE( F2X(0:NXP,0:NZV-1,0:NY+1) )
      IF ( allocated(CF2X) ) THEN
	DEALLOCATE( CF2X, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating CF2X"
      ENDIF

    ELSE

      IF (.not. allocated(CF2X) ) ALLOCATE( CF2X(0:NX2P,0:NZV-1,0:NY+1))
      IF ( allocated(F2X) ) THEN
	DEALLOCATE(F2X, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating F2X"
      ENDIF

    ENDIF

  RETURN
END

subroutine allocation_F3(FLAG)
  USE ntypes
  USE Domain 
  USE run_variable 
  USE IO,     ONLY: I_OUT
  
    IMPLICIT NONE
    
    LOGICAL FLAG 
    INTEGER           :: s_1 
    
    s_1 = 0
    
    IF (FLAG) THEN
    
      IF (.not. allocated(F3X) ) ALLOCATE( F3X(0:NXP,0:NZV-1,0:NY+1) )
      IF ( allocated(CF3X) ) THEN
	DEALLOCATE( CF3X, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating CF3X"
      ENDIF
    
    ELSE

      IF (.not. allocated(CF3X) ) ALLOCATE( CF3X(0:NX2P,0:NZV-1,0:NY+1))
      IF ( allocated(F3X) ) THEN
	DEALLOCATE(F3X, stat=s_1)
      ENDIF
      IF (s_1.NE.0) THEN
	write(I_OUT,*) "Error deallocating F3X"
      ENDIF

    ENDIF 
  
  RETURN 
END

SUBROUTINE allocation_Fth (FLAG)
  USE ntypes
  USE Domain
  USE run_variable
  USE IO,     ONLY: I_OUT

    IMPLICIT NONE
  
    LOGICAL FLAG
    INTEGER           :: s_1

    s_1 = 0

    IF (FLAG)THEN

    IF (.not. allocated(FTHX) ) ALLOCATE( FTHX(0:NXP,0:NZV-1,0:NY+1,1:N_TH) ) 
    IF ( allocated(CFTHX ) ) THEN
      DEALLOCATE( CFTHX, stat=s_1)
    ENDIF
    IF (s_1.NE.0) THEN
      write(I_OUT,*) "Error deallocating CRTHX"
    ENDIF

    ELSE

    IF (.not. allocated(CFTHX ) ) ALLOCATE( CFTHX(0:NX2P,0:NZV-1,0:NY+1,1:N_TH) ) 
    IF ( allocated(FTHX) ) THEN
      DEALLOCATE(FTHX, stat=s_1)
    ENDIF
    IF (s_1.NE.0) THEN
      write(I_OUT,*) "Error deallocating RTHX"
    ENDIF

    ENDIF

  RETURN
END


SUBROUTINE allocate_temps
  USE Domain
  USE Grid
  USE variable_stat

  USE IO,     ONLY: I_OUT
 
    IMPLICIT NONE
    !Passed Variables
    INTEGER           :: ok

    !Local Variables
    INTEGER                    :: s_1
    ok=0
    s_1=0

    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    ! Variables for outputting statistics
    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    ALLOCATE     (UBAR(0:NY+1),VBAR(0:NY+1),WBAR(0:NY+1), stat=s_1)
    ALLOCATE     (URMS(0:NZ+1,0:NY+1),VRMS(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (WRMS(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (UV(0:NZ+1,0:NY+1),UW(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (WV(0:NZ+1,0:NY+1),PV(0:NZ+1,0:NY+1), PU(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (DWDY(0:NZ+1,0:NY+1), DUDY(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (DWDZ(0:NZ+1,0:NY+1), DUDZ(0:NZ+1,0:NY+1), stat=s_1)

    ALLOCATE     (UU(0:NZ+1,0:NY+1),VV(0:NZ+1,0:NY+1), WW(0:NZ+1,0:NY+1),stat=s_1)
    ALLOCATE     (UUU(0:NZ+1,0:NY+1),VVV(0:NZ+1,0:NY+1), WWW(0:NZ+1,0:NY+1),stat=s_1)
    ALLOCATE     (UUV(0:NZ+1,0:NY+1),WWV(0:NZ+1,0:NY+1),stat=s_1)
    ALLOCATE     (UUUU(0:NZ+1,0:NY+1),VVVV(0:NZ+1,0:NY+1), WWWW(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (UUUV(0:NZ+1,0:NY+1),WWWV(0:NZ+1,0:NY+1),stat=s_1)

    ALLOCATE     (dWdz1(0:NZ+1,0:NY+1),dWdz2(0:NZ+1,0:NY+1),dWdz3(0:NZ+1,0:NY+1),stat=s_1)
    ALLOCATE     (dWdz4(0:NZ+1,0:NY+1),dWdz5(0:NZ+1,0:NY+1),stat=s_1)
    
    ALLOCATE     (SHEAR(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (OMEGA_X(0:NZ+1,0:NY+1),OMEGA_Y(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (OMEGA_Z(0:NZ+1,0:NY+1),TKE(0:NZ+1,0:NY+1), stat=s_1)

    IF (s_1.NE.0) THEN
      write(I_OUT,*) "Error Allocating statistics Variables" 
    ENDIF
    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    ! Variables needed for SAVE_STATS_TH
    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    ALLOCATE     (THBAR(0:NY+1,1:N_TH), PE_DISS(0:NY+1,1:N_TH), stat=s_1)
    ALLOCATE     (THRMS(0:NZ+1,0:NY+1,1:N_TH),Rig(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (THV(0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
    ALLOCATE     (THW(0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
    ALLOCATE     (DTHDY(0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
    ALLOCATE     (DTHDZ(0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
    IF (s_1.NE.0) THEN
      write(I_OUT,*) "Error Allocating statistics Variables for theta" 
    ENDIF
    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      
    ! Variables for tkebudget
    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    ALLOCATE     (EPSILON(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (tke_mean(0:NZ+1,0:NY+1),tke_mean_old(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (energy_mean(0:NZ+1,0:NY+1),energy_mean_old(0:NZ+1,0:NY+1), stat=s_1)
    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

! ALLOCATE tke variavles 

    ALLOCATE (tke_1(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (tke_2(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (tke_2_1(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (tke_2_2(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (tke_3(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (tke_3_1(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (tke_3_2(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (tke_3_3(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (tke_4(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (tke_5(0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
    ALLOCATE (tke_6_1(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (tke_6_1_1(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (tke_6_1_2(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (tke_6_1_3(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (tke_6_2(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (tke_6_2_1(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (tke_6_2_2(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (tke_6_2_3(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (tke_7(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (S1_mean(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (p_mean(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE (transport(0:NZ+1,0:NY+1), stat=s_1)
    
    IF (s_1.NE.0) THEN
      write(I_OUT,*) "Error Allocating statistics Variables for tke" 
    ENDIF      


    IF (s_1.NE.0) write(I_OUT,'(a,i5)') "ERROR ALLOCATING TEMPS: stat=",s_1
      ok=s_1

  RETURN
END 

SUBROUTINE deallocate_temps
  USE IO,     ONLY: I_OUT
  USE variable_stat

   IMPLICIT NONE
    !Passed Variables
    INTEGER       :: ok

    !Local Variables
    INTEGER                   :: s_1
    s_1=0
    ok=0

    ! Deallocate Variables for outputting statistics
	DEALLOCATE (UBAR,VBAR,WBAR)
	DEALLOCATE (URMS,VRMS, &
		WRMS, UV,UW, WV,PV,PU, &
		DWDY, DUDY,DWDZ, DUDZ, &
		SHEAR, OMEGA_X,OMEGA_Y,OMEGA_Z,TKE)
    ! deallocate variables needed for SAVE_STATS_TH
	
	  DEALLOCATE (THBAR, PE_DISS, Rig)
	  DEALLOCATE (THRMS,THV, THW, DTHDY, DTHDZ)

    ! Variables for tkebudget
	  DEALLOCATE (epsilon)
!		  tke_mean,tke_mean_old,energy_mean,energy_mean_old)     

    ! allocate tke variavles
	    DEALLOCATE (tke_1,tke_2,tke_2_1, &
		tke_2_2,tke_3,tke_3_1,tke_3_2,tke_3_3,tke_4,tke_5, &
		tke_6_1,tke_6_1_1,tke_6_1_2,tke_6_1_3,tke_6_2,     &
		tke_6_2_1,tke_6_2_2,tke_6_2_3,tke_7,S1_mean,p_mean, &
		transport)



    IF (s_1.NE.0) write(I_OUT,'(a31,i5)') "ERROR DEALLOCATING TEMPS: stat=",s_1
    ok=s_1

  RETURN
END 

SUBROUTINE allocation_les_var
  USE ntypes
  USE Domain
  USE run_variable, ONLY : NU_T, KAPPA_T
  USE les_chan_var
  USE IO,     ONLY: I_OUT


    IMPLICIT NONE

    !Passed Variables
    INTEGER        :: ok

    !Local Variables
    INTEGER           :: s_1
    !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    !  variables
    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    ALLOCATE     (tke_sgs_t(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (tke_sgs_p(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (tke_sgs_diss(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (Sij_mean(0:NZ+1,0:NY+1,1:6), stat=s_1)
    ALLOCATE     (TAU_mean(0:NZ+1,0:NY+1,1:6), stat=s_1)
    ALLOCATE     (NU_T_mean(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (NU_T(0:NXP,0:NZV-1,0:NY+1),  stat=s_1)
    ALLOCATE     (KAPPA_T(0:NXP,0:NZV-1,0:NY+1,1), stat=s_1)

    ALLOCATE     (U_2BAR(0:NXP,0:NZV-1,0:NY+1,1:3), stat=s_1)
    ALLOCATE     (U_BAR_TIL(0:NXP,0:NZV-1,0:NY+1,1:3), stat=s_1)
    ALLOCATE     (U_4BAR(0:NXP,0:NZV-1,0:NY+1,1:3), stat=s_1)
    ALLOCATE     (S_2BAR(0:NXP,0:NZV-1,0:NY+1), stat=s_1)
    ALLOCATE     (U1_bar_les(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (U2_bar_les(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (U3_bar_les(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (DELTA_Y(0:NZ+1,0:NY+1), stat=s_1)
    ALLOCATE     (C_DYN(0:NZ+1,0:NY+1), stat=s_1)

    ALLOCATE     (C_DYN_TH(0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
    ALLOCATE     (KAPPA_T_MEAN(0:NZ+1,0:NY+1), stat=s_1)

    ALLOCATE     (GMAT_11_z(0:NZ+2,0:NY+2,1:2), stat=s_1)
    ALLOCATE     (GMAT_11_y(0:NZ+2,0:NY+2,1:2), stat=s_1)
    ALLOCATE     (GMAT_22_z(0:NZ+2,0:NY+2,1:2), stat=s_1)
    ALLOCATE     (GMAT_22_y(0:NZ+2,0:NY+2,1:2), stat=s_1)
    ALLOCATE     (GMAT_12_z(0:NZ+2,0:NY+2,1:2), stat=s_1)
    ALLOCATE     (GMAT_12_y(0:NZ+2,0:NY+2,1:2), stat=s_1)

    ALLOCATE     (ST_rij(0:NXP,0:NZ+1,0:NY+1,1:6), stat=s_1)


    IF (s_1.NE.0) THEN
      write(I_OUT,*) "Error Allocating LES Variables "
      GOTO 1000
    ENDIF

    write(I_OUT,'(a)') "ALLOCATION FOR LES COMPLETED"
    return
    1000 continue
    ok = s_1
    write(I_OUT,'(a)') "ALLOCATION LES VAR FAILED"

  RETURN
END


SUBROUTINE allocate_les_tmp
  USE Domain
  USE Grid
  USE les_chan_var
  USE IO,     ONLY: I_OUT

    IMPLICIT NONE
    !Passed Variables
    INTEGER           :: ok

    !Local Variables
    INTEGER                    :: s_1
    ok=0
    s_1=0

    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    ! Variables for les  c_dyn
    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    ALLOCATE (numerator(0:NXP,0:NZ+1,0:NY+1), stat=s_1 )
    ALLOCATE (denominator(0:NXP,0:NZ+1,0:NY+1), stat=s_1 )
    ALLOCATE (denominator_sum(0:NZ+1,0:NY+1), stat=s_1 )
    ALLOCATE (numerator_sum(0:NZ+1,0:NY+1), stat=s_1 )

  RETURN
END 


SUBROUTINE deallocate_les_tmp
  USE Domain
  USE Grid
  USE les_chan_var
  USE IO,     ONLY: I_OUT

    IMPLICIT NONE
    !Passed Variables
    INTEGER           :: ok

    !Local Variables
    INTEGER                    :: s_1
    ok=0
    s_1=0

    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    ! Variables for les  c_dyn
    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    DEALLOCATE (numerator, stat=s_1 )
    DEALLOCATE (denominator, stat=s_1 )
    DEALLOCATE (denominator_sum, stat=s_1 )
    DEALLOCATE (numerator_sum, stat=s_1 )
    
  RETURN
end 
