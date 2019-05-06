!*|-------------------------------------------------------|
SUBROUTINE INT_MPI 
!*|-------------------------------------------------------|
  USE ntypes
  USE Domain
  USE Grid
  USE   mpi_var    
  
    IMPLICIT NONE
      
    CHARACTER*35 MPI_IO_NUM


! This subroutine initializes all mpi variables

    CALL MPI_INIT(IERROR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCES,IERROR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,RANK,IERROR)

    write(6,'(a,2i6,i12)') 'MPI INITIALIZED: NPROCS, RANK, MPI_COMM_WORLD= ',NPROCES,RANK,MPI_COMM_WORLD

! Set a string to determine which input/output files to USE
! When MPI is USEd, each process will read/write to files with the
! number of their rank (+1) appended to the end of the file.
! The string MPI_IO_NUM will be USEd to define the RANK+1 for each process
      IF (NPROCES.le.10) THEN
	MPI_IO_NUM=CHAR(MOD(RANK+1,10)+48)
      ELSE IF (NPROCES.le.100) THEN
	MPI_IO_NUM=CHAR(MOD(RANK+1,100)/10+48)   &
		//CHAR(MOD(RANK+1,10)+48)
      ELSE IF (NPROCES.le.1000) THEN
	MPI_IO_NUM=CHAR(MOD(RANK+1,1000)/100+48) &
		//CHAR(MOD(RANK+1,100)/10+48)   &
		//CHAR(MOD(RANK+1,10)+48)
      ELSE IF (NPROCES.le.10000) THEN
	MPI_IO_NUM=CHAR(MOD(RANK+1,10000)/1000+48) &
		//CHAR(MOD(RANK+1,1000)/100+48)   &
		//CHAR(MOD(RANK+1,100)/10+48)     &
		//CHAR(MOD(RANK+1,10)+48)
      ELSE
	  WRITE(6,*) 'ERROR, NPROCES>10,000, Unsupported problem size'
      END IF

  RETURN
END


SUBROUTINE MPI_TRANSPOSE_REAL_X_TO_Z(IN_R,OUT_R)

! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite dIFferences over the local gridpoints.
! We need to update the contents of the ghost zero cells at the start of
! les process

  USE ntypes
  USE Domain
  USE Grid
  USE   mpi_var

    IMPLICIT NONE

      INTEGER i,j,k,N
      INTEGER(r8) :: M, SEND_NO 
      REAL(r8),DIMENSION(0:NXP,0:NZV-1,0:NY+1) :: IN_R
      REAL(r8),DIMENSION(0:NXV-1,0:NZP,0:NY+1) :: OUT_R  
!      REAL(r8),DIMENSION(1:NXV*(NY+2)*NZV/NP)  :: M_IN,M_OUT

      DO N=0,NP-1
       DO J=0,NY+1
        DO I=0,NXP
         DO K=N*(NZP+1)+1,(N+1)*(NZP+1) 
          M= (K-N*(NZP+1)) + (NZP+1)*I + (NXP+1)*(NZP+1)*J + (NXP+1)*(NZP+1)*(NY+2)*N
          M_IN(M) = IN_R(I,K-1,J) 
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      SEND_NO = (NXP+1)*(NZP+1)*(NY+2)
!      write(6,*)'The Size of M', SEND_NO,NXV*(NY+2)*NZV/NP,shape(M_IN), shape(M_OUT),MPI_COMM_WORLD 

      CALL MPI_ALLTOALL(M_IN,SEND_NO,MPI_DOUBLE_PRECISION,M_OUT,SEND_NO,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR) 

      OUT_R(:,:,:)=0.d0

      DO N=0,NP-1
       DO J=0,NY+1
        DO I=N*(NXP+1)+1,(N+1)*(NXP+1)
         Do K=0,NZP
          M= (NZP+1)*(I-N*(NXP+1)-1) + (K+1) + (NXP+1)*(NZP+1)*J + (NXP+1)*(NZP+1)*(NY+2)*N
          OUT_R(I-1,K,J) = M_OUT(M)
         ENDDO
        ENDDO
       ENDDO
      ENDDO


  RETURN
END

!--|*-----------------------------------------------------
SUBROUTINE MPI_TRANSPOSE_REAL_Z_TO_X(IN_R,OUT_R)
!--|*-----------------------------------------------------

! This subroutine is part of the MPI package for the duct flow
! Diablo package.
! Here, we define a set of ghost cells on each process

  USE ntypes
  USE Domain
  USE Grid
  USE   mpi_var

    IMPLICIT NONE

      INTEGER i,j,k,N
      INTEGER(r8) :: M, SEND_NO 
      REAL(r8),DIMENSION(0:NXP,0:NZV-1,0:NY+1) :: OUT_R
      REAL(r8),DIMENSION(0:NXV-1,0:NZP,0:NY+1) :: IN_R  
!      REAL(r8),DIMENSION(1:NXV*(NY+2)*NZV/NP)  :: M_IN,M_OUT

      DO N=0,NP-1
       DO J=0,NY+1
        DO I=N*(NXP+1)+1,(N+1)*(NXP+1)
         Do K=0,NZP
          M= (NZP+1)*(I-N*(NXP+1)-1) + (K+1) + (NXP+1)*(NZP+1)*J + (NXP+1)*(NZP+1)*(NY+2)*N  ! Bishakh
!          M = (I-N*(NXP+1)) + (NXP+1)*K + (NXP+1)*(NZP+1)*J + (NXP+1)*(NZP+1)*(NY+2)*N      ! Eric 
          M_IN(M)=IN_R(I-1,K,J) 
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      SEND_NO = (NXP+1)*(NZP+1)*(NY+2)
!      write(6,*)'The Size of M', SEND_NO,NXV*(NY+2)*NZV/NP,shape(M_IN), shape(M_OUT),MPI_COMM_WORLD 

      CALL MPI_ALLTOALL(M_IN,SEND_NO,MPI_DOUBLE_PRECISION,M_OUT,SEND_NO,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR) 

      OUT_R(:,:,:)=0.d0

      DO N=0,NP-1
       DO J=0,NY+1
        DO I=0,NXP
         DO K=N*(NZP+1)+1,(N+1)*(NZP+1)
          M= (K-N*(NZP+1)) + (NZP+1)*I + (NXP+1)*(NZP+1)*J + (NXP+1)*(NZP+1)*(NY+2)*N   ! Bishakh
!           M = (I+1) + (NXP+1)*(K-N*(NZP+1)-1) + (NXP+1)*(NZP+1)*J + (NXP+1)*(NZP+1)*(NY+2)*N ! Eric
          OUT_R(I,K-1,J)=M_OUT(M) 
         ENDDO
        ENDDO
       ENDDO
      ENDDO

  RETURN
END

!--|*-----------------------------------------------------
SUBROUTINE MPI_TRANSPOSE_COMPLEX_Z_TO_X(IN_CZ,OUT_CX)
!--|*-----------------------------------------------------

! This subroutine is part of the MPI package for the duct flow
! Diablo package.
! Here, we define a set of ghost cells on each process

  USE ntypes
  USE Domain
  USE Grid
  USE   mpi_var

    IMPLICIT NONE

      INTEGER i,j,k,N
      INTEGER(r8) :: M, SEND_NO 
      complex(r8),DIMENSION(0:NX2P,0:NZV-1,0:NY+1) :: OUT_CX
      complex(r8),DIMENSION(0:NX2V-1,0:NZP,0:NY+1) :: IN_CZ  
!      complex(r8),DIMENSION(1:NX2V*(NY+2)*NZV/NP)  :: M_IN_C,M_OUT_C

      DO N=0,NP-1
       DO J=0,NY+1
        DO I=N*(NX2P+1)+1,(N+1)*(NX2P+1)
         Do K=0,NZP
          M= (NZP+1)*(I-N*(NX2P+1)-1) + (K+1) + (NX2P+1)*(NZP+1)*J + (NX2P+1)*(NZP+1)*(NY+2)*N  ! Bishakh
!          M = (I-N*(NX2P+1)) + (NX2P+1)*K + (NX2P+1)*(NZP+1)*J + (NX2P+1)*(NZP+1)*(NY+2)*N      ! Eric 
          M_IN_C(M)=IN_CZ(I-1,K,J) 
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      SEND_NO = (NX2P+1)*(NZP+1)*(NY+2)
!      write(6,*)'The Size of M', SEND_NO,NXV*(NY+2)*NZV/NP,shape(M_IN), shape(M_OUT),MPI_COMM_WORLD 

      CALL MPI_ALLTOALL(M_IN_C,SEND_NO,MPI_DOUBLE_COMPLEX,M_OUT_C,SEND_NO,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,IERROR) 

      OUT_CX(:,:,:)=0.d0

      DO N=0,NP-1
       DO J=0,NY+1
        DO I=0,NX2P
         DO K=N*(NZP+1)+1,(N+1)*(NZP+1)
          M= (K-N*(NZP+1)) + (NZP+1)*I + (NX2P+1)*(NZP+1)*J + (NX2P+1)*(NZP+1)*(NY+2)*N   ! Bishakh
!           M = (I+1) + (NX2P+1)*(K-N*(NZP+1)-1) + (NX2P+1)*(NZP+1)*J + (NX2P+1)*(NZP+1)*(NY+2)*N ! Eric
          OUT_CX(I,K-1,J)=M_OUT_C(M) 
         ENDDO
        ENDDO
       ENDDO
      ENDDO

  RETURN
END


SUBROUTINE MPI_TRANSPOSE_COMPLEX_X_TO_Z(IN_C,OUT_C)

! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost zero cells at the start of
! les process

  USE ntypes
  USE Domain
  USE Grid
  USE   mpi_var

    IMPLICIT NONE

      INTEGER i,j,k,N
      INTEGER(r8) :: M, SEND_NO 
      complex(r8),DIMENSION(0:NX2P,0:NZV-1,0:NY+1) :: IN_C
      complex(r8),DIMENSION(0:NX2V-1,0:NZP,0:NY+1) :: OUT_C  
!      complex(r8),DIMENSION(1:NX2V*(NY+2)*NZV/NP)  :: M_IN_C,M_OUT_C

      DO N=0,NP-1
       DO J=0,NY+1
        DO I=0,NX2P
         DO K=N*(NZP+1)+1,(N+1)*(NZP+1) 
          M= (K-N*(NZP+1)) + (NZP+1)*I + (NX2P+1)*(NZP+1)*J + (NX2P+1)*(NZP+1)*(NY+2)*N
          M_IN_C(M) = IN_C(I,K-1,J) 
         ENDDO
        ENDDO
       ENDDO
      ENDDO


      SEND_NO = (NX2P+1)*(NZP+1)*(NY+2)
!      write(6,*)'The Size of M', SEND_NO,NXV*(NY+2)*NZV/NP,shape(M_IN), shape(M_OUT),MPI_COMM_WORLD 

      CALL MPI_ALLTOALL(M_IN_C,SEND_NO,MPI_DOUBLE_COMPLEX,M_OUT_C,SEND_NO,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,IERROR) 

      OUT_C(:,:,:)=(0.d0,0.d0)

      DO N=0,NP-1
       DO J=0,NY+1
        DO I=N*(NX2P+1)+1,(N+1)*(NX2P+1)
         Do K=0,NZP
          M= (NZP+1)*(I-N*(NX2P+1)-1) + (K+1) + (NX2P+1)*(NZP+1)*J + (NX2P+1)*(NZP+1)*(NY+2)*N
          OUT_C(I-1,K,J) = M_OUT_C(M)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

  RETURN
END
   
SUBROUTINE REAL_FOURIER_TRANS_U1 (flag) 
  USE ntypes 
  USE Domain 
  USE Grid
  USE TIME_STEP_VAR
  USE run_variable 
  USE mpi_var
    
    IMPLICIT NONE
    INTEGER :: i
    LOGICAL :: flag              

    IF (flag) THEN
      S1X=U1X
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
      CALL allocation_u(.false.)
      CU1X=CS1X

    ELSE 

      CS1X=CU1X
      CALL allocation_u(.true.)
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
      U1X=S1X
    ENDIF

  RETURN 
END  

SUBROUTINE REAL_FOURIER_TRANS_U2 (flag) 
  USE ntypes 
  USE Domain 
  USE Grid
  USE TIME_STEP_VAR
  USE run_variable 
  USE mpi_var
    
    IMPLICIT NONE
    INTEGER :: i
    LOGICAL :: flag            

     
    IF (flag) THEN
      S1X=U2X
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
      CALL allocation_v(.false.)
      CU2X=CS1X

    ELSE 
  
      CS1X=CU2X
      CALL allocation_v(.true.)
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
      U2X=S1X
    ENDIF

  RETURN 
END  
 
SUBROUTINE REAL_FOURIER_TRANS_U3 (flag) 
  USE ntypes 
  USE Domain 
  USE Grid
  USE TIME_STEP_VAR
  USE run_variable 
  USE mpi_var
    
    IMPLICIT NONE
    INTEGER :: i
    LOGICAL :: flag             
    
    IF (flag) THEN
      S1X=U3X
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
      CALL allocation_w(.false.)
      CU3X=CS1X

    ELSE 
    
      CS1X=CU3X
      CALL allocation_w(.true.)
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
      U3X=S1X
    ENDIF

  RETURN 
END  

SUBROUTINE REAL_FOURIER_TRANS_P (flag) 
  USE ntypes 
  USE Domain 
  USE Grid
  USE TIME_STEP_VAR
  USE run_variable 
  USE mpi_var
    
    IMPLICIT NONE
    INTEGER :: i
    LOGICAL :: flag              
  
    IF (flag) THEN
      S1X=PX
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
      CALL allocation_p(.false.)
      CPX=CS1X

    ELSE 

      CS1X=CPX
      CALL allocation_p(.true.)
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
      PX=S1X
    ENDIF

  RETURN 
END  
 
 
SUBROUTINE REAL_FOURIER_TRANS_TH (flag) 
  USE ntypes 
  USE Domain 
  USE Grid
  USE TIME_STEP_VAR
  USE run_variable 
  USE mpi_var
    
    IMPLICIT NONE
    INTEGER :: i
    LOGICAL :: flag            

    IF (flag) THEN
      S1X=THX(:,:,:,1)
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
      CALL allocation_th(.false.)
      CTHX(:,:,:,1)=CS1X

    ELSE 
      
      CS1X=CTHX(:,:,:,1)
      CALL allocation_th(.true.)
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
      THX(:,:,:,1)=S1X
    ENDIF

  RETURN 
END  
 
SUBROUTINE REAL_FOURIER_TRANS_R1 (flag) 
  USE ntypes 
  USE Domain 
  USE Grid
  USE TIME_STEP_VAR
  USE run_variable 
  USE mpi_var
    
    IMPLICIT NONE
    INTEGER :: i
    LOGICAL :: flag             
    
    IF (flag) THEN
      S1X=R1X
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
      CALL allocation_R1 (.false.)
      CR1X=CS1X

    ELSE 

      CS1X=CR1X
      CALL allocation_R1(.true.)
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
      R1X=S1X
    ENDIF

  RETURN 
END  


SUBROUTINE REAL_FOURIER_TRANS_R2 (flag) 
  USE ntypes 
  USE Domain 
  USE Grid
  USE TIME_STEP_VAR
  USE run_variable 
  USE mpi_var
    
    IMPLICIT NONE
    INTEGER :: i
    LOGICAL :: flag            

     
    IF (flag) THEN
      S1X=R2X
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
      CALL allocation_R2 (.false.)
      CR2X=CS1X

    ELSE 

      CS1X=CR2X
      CALL allocation_R2(.true.)
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
      R2X=S1X
    ENDIF

  RETURN 
END  

SUBROUTINE REAL_FOURIER_TRANS_R3 (flag) 
  USE ntypes 
  USE Domain 
  USE Grid
  USE TIME_STEP_VAR
  USE run_variable 
  USE mpi_var
    
    IMPLICIT NONE
    INTEGER :: i
    LOGICAL :: flag             

    
    IF (flag) THEN
      S1X=R3X
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
      CALL allocation_R3 (.false.)
      CR3X=CS1X

    ELSE 

      CS1X=CR3X
      CALL allocation_R3(.true.)
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
      R3X=S1X
    ENDIF

  RETURN 
END  

SUBROUTINE REAL_FOURIER_TRANS_Rth (flag) 
  USE ntypes 
  USE Domain 
  USE Grid
  USE TIME_STEP_VAR
  USE run_variable 
  USE mpi_var
      
    IMPLICIT NONE
    INTEGER :: i
    LOGICAL :: flag             

    
    IF (flag) THEN
      S1X=RTHX(:,:,:,1)
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
      CALL allocation_Rth(.FALSE.)
      CRTHX(:,:,:,1)=CS1X

    ELSE 

      CS1X=CRTHX(:,:,:,1)
      CALL allocation_Rth(.TRUE.)
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
      RTHX(:,:,:,1)=S1X
      
    ENDIF

  RETURN 
END  
 
 
      
SUBROUTINE MPI_BCAST_REAL(IN_R,Z_SIZE,Y_SIZE)

! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost zero cells at the start of
! les process

  USE ntypes
  USE Domain
  USE Grid
  USE   mpi_var
    
    IMPLICIT NONE

    INTEGER send_rank,Z_SIZE,Y_SIZE
    REAL(r8),DIMENSION(0:Z_SIZE*Y_SIZE-1) :: IN_R

    send_rank = 0


    CALL MPI_BCAST(IN_R,Z_SIZE*Y_SIZE,MPI_DOUBLE_PRECISION,send_rank,MPI_COMM_WORLD,IERROR) 

  RETURN
END


SUBROUTINE MPI_BCAST_COMPLEX(IN_R,Z_SIZE,Y_SIZE)

! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost zero cells at the start of
! les process

  USE ntypes
  USE Domain
  USE Grid
  USE   mpi_var

    IMPLICIT NONE

    INTEGER send_rank,Z_SIZE,Y_SIZE
    complex(r8),DIMENSION(0:Z_SIZE*Y_SIZE-1) :: IN_R

    send_rank = 0


    CALL MPI_BCAST(IN_R,Z_SIZE*Y_SIZE,MPI_DOUBLE_COMPLEX,send_rank,MPI_COMM_WORLD,IERROR)

  RETURN
END

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE MPI_COMBINE_STATS(STAT,Z_SIZE,Y_SIZE)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  USE ntypes
  USE Domain
  USE Grid
  USE   mpi_var
    
    IMPLICIT NONE    
  
    INTEGER :: I, J, Z_SIZE,Y_SIZE,send_rank
    REAL(r8),DIMENSION(0:Z_SIZE*Y_SIZE-1)      :: STAT
    REAL(r8),DIMENSION(0:Z_SIZE*Y_SIZE-1,1:NP) :: STAT_TMP

    send_rank = 0
    CALL MPI_GATHER(STAT,Z_SIZE*Y_SIZE,MPI_DOUBLE_PRECISION,STAT_TMP,Z_SIZE*Y_SIZE,MPI_DOUBLE_PRECISION,send_rank,MPI_COMM_WORLD,IERROR)
    IF ( RANK .EQ. 0) THEN
      STAT(:) = 0D0
      DO I = 1,NP
	DO J =0, Z_SIZE*Y_SIZE-1
	STAT(J) = STAT(J) + STAT_TMP(J,I)
	ENDDO
      ENDDO
    ENDIF
 
  RETURN
END 


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE MPI_COURANT(DT)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  USE ntypes
  USE Domain
  USE Grid
  USE mpi_var
  
    IMPLICIT NONE

    REAL(r8),INTENT(INOUT) :: DT
    REAL(r8),DIMENSION(1:NP) :: IPACK
    REAL(r8) :: OPACK

  ! READ IN LOCAL DT BASED ON CFL
    OPACK = DT

  ! SEND ALL LOCAL DT'S TO EVERY PROCESS
    CALL MPI_ALLGATHER(OPACK,1,MPI_DOUBLE_PRECISION,IPACK,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)

    DT = MINVAL(IPACK)

  RETURN
END


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE FILTER_VAR(filter_type)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  USE ntypes
  USE Domain
  USE Grid
  USE mpi_var
  USE run_variable, ONLY: S1X,S1Z

    IMPLICIT NONE

    INTEGER :: filter_type

    CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
    CALL les_filter_chan (S1Z,0,NZP,0,NY+1,filter_type)
    CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)


  RETURN
END

     !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE MPI_AVG(IN)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  USE ntypes
  USE Domain
  USE Grid
  USE mpi_var
  
    IMPLICIT NONE
  
    REAL(r8),INTENT(INOUT) :: IN
    REAL(r8),DIMENSION(1:NP) :: IPACK
    REAL(r8) :: OPACK

    OPACK = IN

  ! SEND ALL LOCAL IN'S TO EVERY PROCESS
    CALL MPI_ALLGATHER(OPACK,1,MPI_DOUBLE_PRECISION,IPACK,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)
   
    IN = SUM(IPACK(1:NP))/NP
  
  RETURN
END
