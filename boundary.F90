
      
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_REAL_SP(A,B,C,G,D,K,NY,NX)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! C Uses the Thomas algorithm to solve Ax=b for tridiagonal A
! C The RHS vector and solution are real
! C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
! C Returns solution in x
! C The indexing should be done by ROW, ie.
! C [ b1  c1   0   0   0 ...
! C [ a2  b2  c2   0   0 ...
! C [  0  a3  b3   c3  0 ...

      INTEGER I, J, NX, NY,k
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY), G(0:NX,0:NY)
      REAL*8 D

       J= K

       DO I=0,NX
! c          D=-D/A(I,J+1)
          A(I,J+1)=A(I,J+1)-B(I,J)*D/A(I,J)
          B(I,J+1)=B(I,J+1)-C(I,J)*D/A(I,J)
          G(I,J+1)=G(I,J+1)-G(I,J)*D/A(I,J)
       END DO



      DO J=0,NY-1
        DO I=0,NX
          A(I,J+1)=-A(I,J+1)/B(I,J)
          B(I,J+1)=B(I,J+1)+A(I,J+1)*C(I,J)
          G(I,J+1)=G(I,J+1)+A(I,J+1)*G(I,J)
        END DO
      END DO
      DO I=0,NX
        G(I,NY)=G(I,NY)/B(I,NY)
      END DO
      DO J=NY-1,0,-1
        DO I=0,NX
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        END DO
      END DO

      RETURN
      END




! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_REAL(A,B,C,G,NY,NX)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! C Uses the Thomas algorithm to solve Ax=b for tridiagonal A
! C The RHS vector and solution are real
! C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
! C Returns solution in x
! C The indexing should be done by ROW, ie.
! C [ b1  c1   0   0   0 ...
! C [ a2  b2  c2   0   0 ...
! C [  0  a3  b3   c3  0 ...

      INTEGER I, J, NX, NY
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY), G(0:NX,0:NY)

      DO J=0,NY-1
        DO I=0,NX
          A(I,J+1)=-A(I,J+1)/B(I,J)
          B(I,J+1)=B(I,J+1)+A(I,J+1)*C(I,J)
          G(I,J+1)=G(I,J+1)+A(I,J+1)*G(I,J)
        END DO
      END DO
      DO I=0,NX
        G(I,NY)=G(I,NY)/B(I,NY)
      END DO
      DO J=NY-1,0,-1
        DO I=0,NX
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        END DO
      END DO

      RETURN
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|    
      SUBROUTINE THOMAS_COMPLEX(A,B,C,G,NY,NX)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

! C Uses the Thomas algorithm to solve Ax=b for tridiagonal A
! C The RHS vector and solution is complex
! C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
! C Returns solution in x
! C The indexing should be done by ROW, ie.
! C [ b1  c1   0   0   0 ...
! C [ a2  b2  c2   0   0 ...
! C [  0  a3  b3   c3  0 ...

      INTEGER I, J, NY, NX
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY)
      COMPLEX*16 G(0:NX,0:NY)

      DO J=0,NY-1
        DO I=0,NX
          A(I,J+1)=-A(I,J+1)/B(I,J)
          B(I,J+1)=B(I,J+1)+A(I,J+1)*C(I,J)
          G(I,J+1)=G(I,J+1)+A(I,J+1)*G(I,J)
        END DO
      END DO
      DO I=0,NX
        G(I,NY)=G(I,NY)/B(I,NY)
      END DO
      DO I=0,NX
        DO J=NY-1,0,-1
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        END DO
      END DO

      RETURN
      END


! C----*|--.---------.---------.---------.---------.---------.---------.-|------
SUBROUTINE APPLY_BC_TH_LOWER(N,K)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----
  USE ntypes
  USE Domain
  USE Grid
  ! USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
  USE ADI_var
  USE mg_vari, ONLY : INIT_FLAG
	
    IMPLICIT NONE
    INTEGER I,N,K
    REAL*8 TH_BC_YMIN_var

! C Bottom Wall:
     IF ( TH_BC_YMIN(N) .EQ. 0) THEN
! C Dirichlet
        DO I=0,NXP
!          MATLY(I,0)=0. 
!          MATDY(I,0)=1.
!          MATUY(I,0)=0.                   
!          VECY(I,0)=0.

          MATLY(I,0)=0. 
          MATDY(I,0)=1.
          MATUY(I,0)=0.                   
          VECY(I,0)=TH_BC_YMIN_C1(N) 
        ENDDO
     ELSEIF (TH_BC_YMIN(N) .EQ. 1) THEN
! C Neumann
!        DO I=0,NXP
!          MATLY(I,1)=0.
!          MATDY(I,1)=-1.
!          MATUY(I,1)=1.
!          VECY(I,1)=0.5d0*TH_BC_YMIN_C1(N)*  &
!               (INT_JACOB(K,1)+INT_JACOB(K,2))/CJOB_22(K,2,1)
!        ENDDO

        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=-1.
          MATUY(I,0)=1.
          VECY(I,0)=TH_BC_YMIN_C1(N)/(0.5*(CJOB_22(K,1,2)+CJOB_22(K+1,1,1)))*INT_JACOB(K,1)
        ENDDO 


     ELSEIF (TH_BC_YMIN(N) .EQ. 5) THEN
! C Neumann for finite time (check diablo.F90)

        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=-1.
          MATUY(I,0)=1.
          VECY(I,0)=TH_BC_YMIN_C1(N)/(0.5*(CJOB_22(K,1,2)+CJOB_22(K+1,1,1)))*INT_JACOB(K,1)
        ENDDO

     ELSEIF ( TH_BC_YMIN(N) .EQ. 2) THEN
! C dTHdt = TH_BC_YMIN_C1
        DO I=0,NXP
          TH_BC_YMIN_var = TH_BC_YMIN_C1(N)*(TIME-INT_TIME) 
          !TH_BC_YMIN_var=TH_BC_YMIN_C1(N)*(DELTA_T)
          MATLY(I,0)=0. 
          MATDY(I,0)=1.
          MATUY(I,0)=0.                   
          VECY(I,0) = INT_THX + TH_BC_YMIN_var

!          IF ((W_BC_YMIN.EQ.3) .AND. (IC_TYPE == 7) .AND. (WALL_MODEL_TYPE.NE.2) .AND. (WALL_MODEL_TYPE.NE.4)) THEN 
!          ! in case of wall model, we need to change this becasUSE surface is now between nodes 1 and 2
!	    MATLY(I,1)=0. 
!	    MATDY(I,1)=1.
!	    MATUY(I,1)=0.                   
!	    VECY(I,1) = INT_THX + TH_BC_YMIN_var
!          ENDIF
        ENDDO
      
      ELSEIF (TH_BC_YMIN(N) .EQ. 6) THEN
! C Periodic
        DO I=0,NXP
          MATLY(I,0)=0. 
          MATDY(I,0)=1.
          MATUY(I,0)=0.                   
          VECY(I,0)=THX(I,K,NY,N)
	ENDDO        	 
      ELSEIF (TH_BC_YMIN(N) .EQ. 7) THEN
! C Mixed boundary condition for w at the bottom wall 
       IF ( K .lE. 80 ) THEN
        DO I=0,NXP
          MATLY(I,1)=0.
          MATDY(I,1)=-1.
          MATUY(I,1)=1.
          VECY(I,1)=DY(2)*TH_BC_YMIN_C1(N)          

          MATLY(I,0)=0.
          MATDY(I,0)=-1.
          MATUY(I,0)=1.
          VECY(I,0)=DY(1)*TH_BC_YMIN_C1(N)
        END DO
       ELSE
        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=1.
          MATUY(I,0)=0.
          VECY(I,0)=0.

          MATLY(I,1)=0.
          MATDY(I,1)=1.
          MATUY(I,1)=0.
          VECY(I,1)=TH_BC_YMIN_C1(N)
        ENDDO   
       ENDIF   
      ELSEIF (TH_BC_YMIN(N) .eq. 8) THEN
! C Neumann when sponge is USEd at left and right boundary
        DO I=0,NXP
          MATLY(I,1)=0.
          MATDY(I,1)=-1.
          MATUY(I,1)=1.
! c          VECY(I,1)=(0.5d0*TH_BC_YMIN_C1(N)*
! c     &          (INT_JACOB(K,1)+INT_JACOB(K,2))/CJOB_22(K,2,1))
! c     &          /(1.d0 + SPONGE_SIGMA_OUT(K))
          VECY(I,1)=-(TH_BAR(K,2)-TH_BAR(K,1)) &
               /(1.d0 + SPONGE_SIGMA_OUT(K))
        ENDDO

        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=-1.
          MATUY(I,0)=1.
! c          VECY(I,0)=(0.5d0*TH_BC_YMIN_C1(N)*
! c     &          (INT_JACOB(K,0)+INT_JACOB(K,1))/CJOB_22(K,1,1))
! c     &          /(1.d0 + SPONGE_SIGMA_OUT(K)) 
          VECY(I,0)=-(TH_BAR(K,1)-TH_BAR(K,0)) &
               /(1.d0 + SPONGE_SIGMA_OUT(K)) 
        ENDDO
      ELSE
        write(*,*) 'WARNING: TH_BC_LOWER is of unknown type'
      ENDIF

  RETURN 
END



! C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_TH_UPPER(N,K)
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
USE ntypes
USE Domain
USE Grid
! USE Fft_var
USE TIME_STEP_VAR
USE run_variable
USE ADI_var
USE mg_vari, ONLY : INIT_FLAG
      
IMPLICIT NONE

      INTEGER I, N,K

! C Top wall
      IF (TH_BC_YMAX(N) .EQ. 0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=TH_BC_YMAX_C1(N)

!         MATLY(I,NY)=0.
!         MATDY(I,NY)=1.
!         MATUY(I,NY)=0.
!         VECY(I,NY)=TH_BC_YMAX_C1(N)
        END DO
      ELSE IF (TH_BC_YMAX(N) .eq. 1) THEN
! C Neumann
        DO I=0,NXP
          MATLY(I,NY)=-1.
          MATDY(I,NY)=1.
          MATUY(I,NY)=0.
          VECY(I,NY)=0.5d0*TH_BC_YMAX_C1(N)* &
               (INT_JACOB(K,NY)+INT_JACOB(K,NY-1))/CJOB_22(K,NY-1,1)
        END DO
        DO I=0,NXP
          MATLY(I,NY+1)=-1
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=0.5d0*TH_BC_YMAX_C1(N)* &
               (INT_JACOB(K,NY+1)+INT_JACOB(K,NY))/CJOB_22(K,NY,1)
        END DO 
      ELSE IF (TH_BC_YMAX(N) .eq. 6) THEN	
! c periodic 
       DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=THX(I,K,1,N)
       ENDDO	  
      
      END IF

      RETURN
      END   

! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_TH_LEFT(N,J)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----
USE ntypes
USE Domain
USE Grid
! USE Fft_var
USE TIME_STEP_VAR
USE run_variable
USE ADI_var
USE mg_vari, ONLY : INIT_FLAG
      
IMPLICIT NONE
      INTEGER I,N,J

! C Left Wall:
      IF (TH_BC_ZMIN(N) .EQ.0) THEN
! C Dirichlet
       IF (IC_TYPE .EQ. 5) THEN
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=1.
          MATUX(I,0)=0.
          VECX(I,0)=0.

          MATLX(I,1)=0.
          MATDX(I,1)=1.
          MATUX(I,1)=0.
          VECX(I,1)=U3_BAR(1,J)
        END DO
       ElSE  
        DO I=0,NXP          
          MATLX(I,0)=0. 
          MATDX(I,0)=1.
          MATUX(I,0)=0.                   
          VECX(I,0)=0.

          MATLX(I,1)=0. 
          MATDX(I,1)=1.
          MATUX(I,1)=0.                   
          VECX(I,1)=TH_BC_ZMIN_C1(N) 
        END DO
       ENDIF
      ELSE IF (TH_BC_ZMIN(N).eq.1) THEN
! C Neumann
! c        DO I=0,NXP
! c          MATLX(I,1)=0.
! c          MATDX(I,1)=-1.
! c          MATUX(I,1)=1.
! c          VECX(I,1)=TH_BC_ZMIN_C1(N)
! c        END DO

        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=-1.
          MATUX(I,0)=1.
          VECX(I,0)=TH_BC_ZMIN_C1(N)
        END DO
      ELSE IF (TH_BC_ZMIN(N).eq.6) THEN
! C Periodic
       DO I=0,NXP          
          MATLX(I,0)=0. 
          MATDX(I,0)=1.
          MATUX(I,0)=0.                   
          VECX(I,0)=THX(I,NZ,J,N)
! 
! c          MATLX(I,1)=0. 
! c          MATDX(I,1)=1.
! c          MATUX(I,1)=0.                   
! c          VECX(I,1)=TH_BC_ZMIN_C1(N) 
        END DO

      ELSE
        write(*,*) 'WARNING: TH_BC_LEFT is of unknown type'
      END IF

      RETURN 
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_TH_RIGHT(N,J)
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
USE ntypes
USE Domain
USE Grid
! USE Fft_var
USE TIME_STEP_VAR
USE run_variable
USE ADI_var
USE mg_vari, ONLY : INIT_FLAG
      
IMPLICIT NONE
      INTEGER I,N,J

! C Right wall
      IF (TH_BC_ZMAX(N).EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=0.

          MATLX(I,NZ)=0.
          MATDX(I,NZ)=1.
          MATUX(I,NZ)=0.
          VECX(I,NZ)=TH_BC_ZMAX_C1(N)
        END DO
      ELSE IF (TH_BC_ZMAX(N) .eq. 1) THEN
! C Neumann
! c        DO I=0,NXP
! c          MATLX(I,NZ)=-1.
! c          MATDX(I,NZ)=1.
! c          MATUX(I,NZ)=0.
! c          VECX(I,NZ)=TH_BC_ZMAX_C1(N)
! c        END DO
        DO I=0,NXP
          MATLX(I,NZ+1)=-1.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=TH_BC_ZMAX_C1(N)
        END DO
      ELSE IF (TH_BC_ZMAX(N) .eq. 5) THEN
        DO I=0,NXP
          MATLX(I,NZ+1)=-2.0
          MATDX(I,NZ+1)=1.d0
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1) =0.
        END DO
      ELSE IF (TH_BC_ZMAX(N) .eq. 6) THEN
! C periodic condition      	
	DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=THX(I,1,J,N)
	ENDDO  
      END IF

      RETURN
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_1_LOWER(K)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----
USE ntypes
USE Domain
USE Grid
! USE Fft_var
USE TIME_STEP_VAR
USE run_variable
USE ADI_var
USE mg_vari, ONLY : INIT_FLAG
      
IMPLICIT NONE
      INTEGER I,K

! C Bottom Wall:
      IF (U_BC_YMIN.EQ.0) THEN
! C Dirichlet
       If (F_TYPE.EQ.5)THEN
        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=1.
          MATUY(I,0)=0.
          VECY(I,0)=0.

          MATLY(I,1)=0.
          MATDY(I,1)=1.
          MATUY(I,1)=0.
          VECY(I,1)=-Q_H0*(1.d0/cos(ANG_BETA)+GX(NX/2)*tan(ANG_BETA) &
                  *In_H0)*sin(OMEGA0*TIME)           
! c           VECY(I,1)=0.0
        END DO
       ELSE
        DO I=0,NXP
          MATLY(I,0)=0. 
          MATDY(I,0)=1.
          MATUY(I,0)=0.                   
          VECY(I,0)=U_BC_YMIN_C1

!          MATLY(I,1)=0. 
!          MATDY(I,1)=1.
!          MATUY(I,1)=0.                   
!          VECY(I,1)=U_BC_YMIN_C1 
        END DO
       ENDIF
      ELSE IF (U_BC_YMIN.eq.1) THEN
! C Neumann
!be aware of the location of this BC, it depends on Jsatrt which is solely a function of U_BC_YMIN
	IF (IC_TYPE == 7) THEN
	  DO I=0,NXP
	    MATLY(I,1)=0.
	    MATDY(I,1)=-1.
	    MATUY(I,1)=1.
	    VECY(I,1)=U_BC_YMIN_C1
	  ENDDO
        ELSE
	  DO I=0,NXP
	    MATLY(I,1)=0.
	    MATDY(I,1)=-1.
	    MATUY(I,1)=1.
	    VECY(I,1)=DY(1)*U_BC_YMIN_C1
	  ENDDO
        ENDIF
      ELSE IF (U_BC_YMIN.eq.3) THEN
! C Wall model artificial slip boundary condition
! C Neumann
	IF ( IC_TYPE == 7) THEN
	 IF (WALL_MODEL_TYPE.EQ.2 .OR. WALL_MODEL_TYPE.EQ.1) THEN
!	    DO I=0,NXP
!	      MATLY(I,1)=0.
!	      MATDY(I,1)=1.
!	      MATUY(I,1)=0.
!	      VECY(I,1)=U_BC_LOWER_WALLMODEL(I,K)
!	    ENDDO
    ! C This gridpoint is not USEd here
	    DO I=0,NXP
	      MATLY(I,0)=0.
	      MATDY(I,0)=1.
	      MATUY(I,0)=0.
	      VECY(I,0)=U_BC_LOWER_WALLMODEL(I,K)
	    ENDDO
	  ELSE
	  DO I=0,NXP
	      MATLY(I,0)=0.
	      MATDY(I,0)=-1.
	      MATUY(I,0)=1.
	      VECY(I,0)=U_BC_LOWER_WALLMODEL(I,K)
	    ENDDO
    ! C This gridpoint is not USEd here
!	    DO I=0,NXP
!	      MATLY(I,0)=0.
!	      MATDY(I,0)=1.
!	      MATUY(I,0)=0.
!	      VECY(I,0)=0.
!	    ENDDO
	  ENDIF
        ELSE
	  DO I=0,NXP
	    MATLY(I,1)=0.
	    MATDY(I,1)=-1.
	    MATUY(I,1)=1.
	    VECY(I,1)=DY(2)*U_BC_LOWER_WALLMODEL(I,K)
	  ENDDO
  ! C This gridpoint is not used here
	  DO I=0,NXP
	    MATLY(I,0)=0.
	    MATDY(I,0)=1.
	    MATUY(I,0)=0.
	    VECY(I,0)=0.
	  ENDDO
	ENDIF
      ElSE IF (U_BC_YMIN .EQ. 6) THEN
! C Periodic
        DO I=0,NXP
          MATLY(I,0)=0. 
          MATDY(I,0)=1.
          MATUY(I,0)=0.                   
          VECY(I,0)=U1X(I,K,NY)
	ENDDO	
      ELSE
        write(*,*) 'WARNING: U_BC_LOWER is of unknown type'
      END IF

      RETURN 
      END



! C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_1_UPPER(K)
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
USE ntypes
USE Domain
USE Grid
! USE Fft_var
USE TIME_STEP_VAR
USE run_variable
USE ADI_var
USE mg_vari, ONLY : INIT_FLAG
      
IMPLICIT NONE
      INTEGER I,K

! C Top wall
      IF (U_BC_YMAX.EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=U_BC_YMAX_C1

!          MATLY(I,NY)=0.
!          MATDY(I,NY)=1.
!          MATUY(I,NY)=0.
!          VECY(I,NY)=U_BC_YMAX_C1
        END DO
      ELSE IF (U_BC_YMAX.eq.1) THEN
! C Neumann
        DO I=0,NXP
          MATLY(I,NY+1)=-1.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=U_BC_YMAX_C1
        END DO
        
      ELSE IF (U_BC_YMAX.EQ.3) THEN
! C Use wall model artificial boundary condition
! C Neumann
        DO I=0,NXP
          MATLY(I,NY)=-1.
          MATDY(I,NY)=1.
          MATUY(I,NY)=0.
          VECY(I,NY)=DY(NY)*U_BC_UPPER_NWM(I,K)
        END DO
! C This gridpoint is not USEd here
        DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=U1X(I,K,NY)
        END DO
      ELSE IF (U_BC_YMAX .eq. 6) THEN	
! c periodic 
       DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=U1X(I,K,1)
       ENDDO	
       
      END IF

      RETURN
      END   

! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_1_LEFT(J)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----
USE ntypes
USE Domain
USE Grid
! USE Fft_var
USE TIME_STEP_VAR
USE run_variable
USE ADI_var
USE mg_vari, ONLY : INIT_FLAG
      
IMPLICIT NONE

      INTEGER I,J

! C Left Wall:
      IF (U_BC_ZMIN.EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLX(I,0)=0. 
          MATDX(I,0)=1.
          MATUX(I,0)=0.                   
          VECX(I,0)=U_BC_ZMIN_C1
        IF (IC_TYPE .EQ. 5) THEN
          MATLX(I,1)=0.
          MATDX(I,1)=1.
          MATUX(I,1)=0.
          VECX(I,1)= U1_BAR(J)
        ELSE  
          MATLX(I,1)=0. 
          MATDX(I,1)=1.
          MATUX(I,1)=0.                   
          VECX(I,1)=U_BC_ZMIN_C1
         ENDIF 
        END DO
      ELSE IF (U_BC_ZMIN.eq.1) THEN
! C Neumann
       IF (IC_TYPE .EQ. 7) THEN
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=-1.
          MATUX(I,0)=1.
          VECX(I,0)=U_BC_ZMIN_C1
        END DO
       ELSE
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=-1.
          MATUX(I,0)=1.
          VECX(I,0)=DZ(1)*U_BC_ZMIN_C1
        END DO	  
       ENDIF	
      ELSE IF (U_BC_ZMIN.eq.3) THEN
! C Wall model artificial slip boundary condition
! C Neumann
        DO I=0,NXP
          MATLX(I,1)=0.
          MATDX(I,1)=-1.
          MATUX(I,1)=1.
          VECX(I,1)=DZ(2)*U_BC_LOWER_WALLMODEL(I,J)
        END DO
! C This gridpoint is not USEd here
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=1.
          MATUX(I,0)=0.
          VECX(I,0)=0.
        END DO
      ELSE IF (U_BC_ZMIN .eq. 6) THEN
! C Periodic
       DO I=0,NXP          
          MATLX(I,0)=0. 
          MATDX(I,0)=1.
          MATUX(I,0)=0.                   
          VECX(I,0)=U1X(I,NZ,J)
        END DO	
      ELSE
        write(*,*) 'WARNING: U_BC_LEFT is of unknown type'
      END IF

      RETURN 
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_1_RIGHT(J)
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
USE ntypes
USE Domain
USE Grid
! USE Fft_var
USE TIME_STEP_VAR
USE run_variable
USE ADI_var
USE mg_vari, ONLY : INIT_FLAG
      
IMPLICIT NONE

      INTEGER I,J

! C Right wall
      IF (U_BC_ZMAX.EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=0.

          MATLX(I,NZ)=0.
          MATDX(I,NZ)=1.
          MATUX(I,NZ)=0.
          VECX(I,NZ)=U_BC_ZMAX_C1
        END DO
      ELSE IF (U_BC_ZMAX.eq.1) THEN
! C Neumann
        IF (IC_TYPE == 7) THEN
	 DO I=0,NXP
           MATLX(I,NZ+1)=-1.
           MATDX(I,NZ+1)=1.
           MATUX(I,NZ+1)=0.
           VECX(I,NZ+1)=U_BC_ZMAX_C1
         END DO
	ELSE
         DO I=0,NXP
           MATLX(I,NZ)=-1.
           MATDX(I,NZ)=1.
           MATUX(I,NZ)=0.
           VECX(I,NZ)=DZ(NZ)*U_BC_ZMAX_C1
         END DO
	ENDIF
        DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=0.
        END DO
      ELSE IF (U_BC_ZMAX.EQ.3) THEN
! C Use wall model artificial boundary condition
! C Neumann
        DO I=0,NXP
          MATLX(I,NZ)=-1.
          MATDX(I,NZ)=1.
          MATUX(I,NZ)=0.
          VECX(I,NZ)=DZ(NZ)*U_BC_UPPER_NWM(I,J)
        END DO
! C This gridpoint is not USEd here
        DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=U1X(I,J,NZ)
        END DO
	
      ELSE IF (U_BC_ZMAX .eq. 6) THEN
! C periodic condition      	
	DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=U1X(I,1,J)
	ENDDO	
	
      END IF

      RETURN
      END   

! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_2_LOWER(K)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----
USE ntypes
USE Domain
USE Grid
! USE Fft_var
USE TIME_STEP_VAR
USE run_variable
USE ADI_var
USE mg_vari, ONLY : INIT_FLAG
      
IMPLICIT NONE

      INTEGER I,K

! C Bottom Wall:
      IF (V_BC_YMIN.EQ.0) THEN
! C Dirichlet
       IF ( IC_TYPE .EQ. 7) THEN
        DO I=0,NXP
          MATLY(I,0)=0. 
          MATDY(I,0)=1.
          MATUY(I,0)=0.                   
          VECY(I,0)=V_BC_YMIN_C1 

!          MATLY(I,1)=0. 
!          MATDY(I,1)=1.
!          MATUY(I,1)=0.                   
!          VECY(I,1)=V_BC_YMIN_C1 
        END DO 
       ELSE
        DO I=0,NXP
          MATLY(I,1)=0. 
          MATDY(I,1)=1.
          MATUY(I,1)=0.                   
          VECY(I,1)=V_BC_YMIN_C1

          MATLY(I,2)=0. 
          MATDY(I,2)=1.
          MATUY(I,2)=0.                   
          VECY(I,2)=V_BC_YMIN_C1 
        END DO
       ENDIF	
      ELSE IF (V_BC_YMIN.eq.1) THEN
! C Neumann
	IF (IC_TYPE == 7) THEN
	 DO I=0,NXP
          MATLY(I,1)=0.
          MATDY(I,1)=-1.
          MATUY(I,1)=1.
          VECY(I,1)=V_BC_YMIN_C1
        ENDDO
       ELSE
        DO I=0,NXP
          MATLY(I,1)=0.
          MATDY(I,1)=-1.
          MATUY(I,1)=1.
          VECY(I,1)=DYF(1)*V_BC_YMIN_C1
        ENDDO
       ENDIF
      ELSE IF (V_BC_YMIN.eq.3) THEN
! C Wall model artificial slip boundary condition
! C Neumann
	IF ( IC_TYPE == 7) THEN
	 IF (WALL_MODEL_TYPE.EQ.2 .OR. WALL_MODEL_TYPE.EQ.1) THEN
!	    DO I=0,NXP
!	      MATLY(I,1)=0.
!	      MATDY(I,1)=1.
!	      MATUY(I,1)=0.
!	      VECY(I,1)=V_BC_LOWER_WALLMODEL(I,K)
!	    ENDDO
    ! C This gridpoint is not USEd here
	    DO I=0,NXP
	      MATLY(I,0)=0.
	      MATDY(I,0)=1.
	      MATUY(I,0)=0.
	      VECY(I,0)=V_BC_LOWER_WALLMODEL(I,K)+W_BC_LOWER_WALLMODEL(I,K)*SIN(ATAN(-CJOB_21(K,0,1)/CJOB_22(K,0,1)))
	    ENDDO
	  ELSE
	  DO I=0,NXP
	      MATLY(I,0)=0.
	      MATDY(I,0)=-1.
	      MATUY(I,0)=1.
	      VECY(I,0)=V_BC_LOWER_WALLMODEL(I,K)
	    ENDDO
    ! C This gridpoint is not USEd here
!	    DO I=0,NXP
!	      MATLY(I,0)=0.
!	      MATDY(I,0)=1.
!	      MATUY(I,0)=0.
!	      VECY(I,0)=0.
!	    ENDDO
	  ENDIF
        ELSE
	  DO I=0,NXP
	    MATLY(I,1)=0.
	    MATDY(I,1)=-1.
	    MATUY(I,1)=1.
	    VECY(I,1)=DY(2)*V_BC_LOWER_WALLMODEL(I,K)
	  ENDDO
  ! C This gridpoint is not used here
	  DO I=0,NXP
	    MATLY(I,0)=0.
	    MATDY(I,0)=1.
	    MATUY(I,0)=0.
	    VECY(I,0)=0.
	  ENDDO
	ENDIF
      ElSE IF (V_BC_YMIN .EQ. 6) THEN
! C Periodic
        DO I=0,NXP
          MATLY(I,0)=0. 
          MATDY(I,0)=1.
          MATUY(I,0)=0.                   
          VECY(I,0)=U2X(I,K,NY)
	ENDDO	
      ELSE
        write(*,*) 'WARNING: V_BC_LOWER is of unknown type'
      END IF
! 
! C The following is ONLY a placeholder, this row is USEd for U1 and U3
! c      DO I=0,NXP
! c        MATLY(I,0) = 0.
! c        MATDY(I,0) = 1.
! c        MATUY(I,0) = 0.
! c        VECY(I,0) = 0.
! c      END DO



      RETURN 
      END



! C----*|--.---------.---------.---------.---------.---------.---------.-|----
SUBROUTINE APPLY_BC_2_UPPER(K)
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
  USE ntypes
  USE Domain
  USE Grid
  ! USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
  USE ADI_var
  USE mg_vari, ONLY : INIT_FLAG
      
    IMPLICIT NONE

    INTEGER I,K

! C Top wall
      IF (V_BC_YMAX.EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=V_BC_YMAX_C1

          MATLY(I,NY)=0.
          MATDY(I,NY)=1.
          MATUY(I,NY)=0.
          VECY(I,NY)=V_BC_YMAX_C1
        END DO
      ELSE IF (V_BC_YMAX.eq.1) THEN
! C Neumann
      IF (IC_TYPE ==7 ) THEN
       DO I=0,NXP
          MATLY(I,NY)=-1.
          MATDY(I,NY)=1.
          MATUY(I,NY)=0.
          VECY(I,NY)=V_BC_YMAX_C1
        END DO 
      ELSE	 
        DO I=0,NXP
          MATLY(I,NY+1)=-1.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=DYF(NY)*V_BC_YMAX_C1
        END DO
       ENDIF	
      ELSE IF (V_BC_YMAX.EQ.3) THEN
! C Use wall model artificial boundary condition
! C Neumann
        DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=0.
        END DO
! C This gridpoint is not USEd here
        DO I=0,NXP
          MATLY(I,NY)=0.
          MATDY(I,NY)=1.
          MATUY(I,NY)=0.
          VECY(I,NY)=0.
        END DO
      ELSE IF (V_BC_YMAX .eq. 6) THEN	
! c periodic 
       DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=U2X(I,K,1)
       ENDDO	
      END IF

  RETURN
END   

! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_2_LEFT(J)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----
USE ntypes
USE Domain
USE Grid
! USE Fft_var
USE TIME_STEP_VAR
USE run_variable
USE ADI_var
USE mg_vari, ONLY : INIT_FLAG
      
IMPLICIT NONE

      INTEGER I,J

! C Left Wall:
      IF (V_BC_ZMIN.EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
	IF (IC_TYPE .EQ. 7) THEN
          MATLX(I,0)=0. 
          MATDX(I,0)=1.
          MATUX(I,0)=0.                   
          VECX(I,0)=U2_BAR(1,j)
         
          MATLX(I,1)=0.
          MATDX(I,1)=1.
          MATUX(I,1)=0.
          VECX(I,1)=U2_BAR(1,j) 
         ELSE
	  MATLX(I,0)=0. 
          MATDX(I,0)=1.
          MATUX(I,0)=0.                   
          VECX(I,0)=0.
	 
          MATLX(I,1)=0. 
          MATDX(I,1)=1.
          MATUX(I,1)=0.                   
          VECX(I,1)=V_BC_ZMIN_C1
         ENDIF 
        END DO
      ELSE IF (V_BC_ZMIN.eq.1) THEN
! C Neumann
      IF (IC_TYPE .EQ. 7) THEN
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=-1.
          MATUX(I,0)=1.
          VECX(I,0)=V_BC_ZMIN_C1
        END DO
      ELSE
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=-1.
          MATUX(I,0)=1.
          VECX(I,0)=DZ(1)*V_BC_ZMIN_C1
        END DO 
      ENDIF 
      ELSE IF (V_BC_ZMIN.eq.3) THEN
! C Wall model artificial slip boundary condition
! C Neumann
! c        DO I=0,NXP
! c          MATLX(I,1)=0.
! c          MATDX(I,1)=-1.
! c          MATUX(I,1)=1.
! c          VECX(I,1)=0.
! c        END DO
! C This gridpoint is not USEd here
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=1.
          MATUX(I,0)=0.
          VECX(I,0)=DZ(1)*V_BC_ZMIN_C1
        END DO
      ELSE IF (V_BC_ZMIN .eq. 6) THEN
! C Periodic
       DO I=0,NXP          
          MATLX(I,0)=0. 
          MATDX(I,0)=1.
          MATUX(I,0)=0.                   
          VECX(I,0)=U2X(I,NZ,J)
        END DO	
      ELSE IF (V_BC_ZMIN.EQ.7) THEN
	 DO I=0,NXP
	  MATLX(I,0)=0.
          MATDX(I,0)=(1. - dtc*C_int_le(J)*0.25*  &
          ( CJOB_11(0,J,2) + CJOB_11(1,J,2)        &
          + CJOB_11(0,J,1) + CJOB_11(0,J+1,1))	   &
     	  /INT_JACOB(0,J))
          MATUX(I,0)=(dtc*C_int_le(J)*0.25*       &
         ( CJOB_11(0,J,2) + CJOB_11(1,J,2)         &
          + CJOB_11(0,J,1) + CJOB_11(0,J+1,1))	   &
     	  /INT_JACOB(0,J))
          VECX(I,0) =U2X(I,0,J)
         END DO	
      ELSE
        write(*,*) 'WARNING: V_BC_LEFT is of unknown type'
      END IF

      RETURN 
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_2_RIGHT(J)
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
USE ntypes
USE Domain
USE Grid
! USE Fft_var
USE TIME_STEP_VAR
USE run_variable
USE ADI_var
USE mg_vari, ONLY : INIT_FLAG
      
IMPLICIT NONE

      INTEGER I,J

! C Right wall
      IF (V_BC_ZMAX.EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=0.

          MATLX(I,NZ)=0.
          MATDX(I,NZ)=1.
          MATUX(I,NZ)=0.
          VECX(I,NZ)=V_BC_ZMAX_C1
        END DO
      ELSE IF (V_BC_ZMAX.eq.1) THEN
! C Neumann
       IF (IC_TYPE .EQ. 7) THEN
        DO I=0,NXP
          MATLX(I,NZ+1)=-1. 
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=V_BC_ZMAX_C1
        END DO
       ELSE
        DO I=0,NXP
          MATLX(I,NZ)=-1.
          MATDX(I,NZ)=1.
          MATUX(I,NZ)=0.
          VECX(I,NZ)=DZ(NZ)*V_BC_ZMAX_C1
        END DO
        DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=0.
        END DO
       ENDIF
      ELSE IF (V_BC_ZMAX.EQ.4) THEN
! C  outlet Boundary condition
        DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=U2X(I,NZ-1,J)+(U2X(I,NZ,J)-U2X(I,NZ-1,J)) &
                   *DZ(NZ+1)/DZ(NZ)

          MATLX(I,NZ)=-1.d0/DZ(NZ)
          MATDX(I,NZ)=1.d0/DZ(NZ)+1.d0/DZ(NZ+1)
          MATUX(I,NZ)=-1.d0/DZ(NZ+1)
          VECX(I,NZ) = 0.
        END DO

! c        DO I=0,NXP
! c          MATLX(I,NZ+1)=0.
! c          MATDX(I,NZ+1)=1.
! c          MATUX(I,NZ+1)=0.
! c          VECX(I,NZ+1)=U2(I,NZ-1,J)+(U2(I,NZ,J)-U2(I,NZ-1,J))
! c     +              *DZ(NZ+1)/DZ(NZ)
! 
! c          MATLX(I,NZ)=-(1.d0/DZ(NZ-1)+1.d0/DZ(NZ))
! c          MATDX(I,NZ)=1.d0/DZ(NZ)
! c          MATUX(I,NZ)=0.
! c          VECX(I,NZ) =0.
! c        END DO
      ELSE IF ( V_BC_ZMAX.EQ.5 ) THEN
       IF ( IC_TYPE == 7) THEN
        DO I=0,NXP
          MATLX(I,NZ+1)=-2.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1) =0.
        END DO
       ELSE
        DO I=0,NXP
! c         MATLX(I,NZ)=-(1.d0/DZ(NZ-1)+1.d0/DZ(NZ))
! c          MATDX(I,NZ)=1.d0/DZ(NZ)
! c          MATUX(I,NZ)=0.
! c          VECX(I,NZ) =0.

          MATLX(I,NZ+1)=-(1.d0/DZ(NZ)+1.d0/DZ(NZ+1))
          MATDX(I,NZ+1)=1.d0/DZ(NZ+1)
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1) =0.
        END DO
       ENDIF
      ELSE IF (V_BC_ZMAX.EQ.3) THEN
! C Use wall model artificial boundary condition
! C Neumann
        DO I=0,NXP
          MATLX(I,NZ)=-1.
          MATDX(I,NZ)=1.
          MATUX(I,NZ)=0.
          VECX(I,NZ)=0.0
        END DO
! C This gridpoint is not USEd here
        DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=U2X(I,NZ,J)
        END DO
      ELSE IF (V_BC_ZMAX .eq. 6) THEN
! C periodic condition      	
	DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=U2X(I,1,J)
	ENDDO
      ELSE IF (V_BC_ZMAX.EQ.7) THEN
	 DO I=0,NXP
	  MATLX(I,NZ+1)=-(dtc*C_int(J)*0.25*             &
         ( CJOB_11(NZ+1,J,2) + CJOB_11(NZ+2,J,2)          &
          + CJOB_11(NZ+1,J,1) + CJOB_11(NZ+1,J+1,1))	  & 
     	  /INT_JACOB(NZ+1,J))
          MATDX(I,NZ+1)=(1. + dtc*C_int(J)*0.25*       &
          ( CJOB_11(NZ+1,J,2) + CJOB_11(NZ+2,J,2)       &
          + CJOB_11(NZ+1,J,1) + CJOB_11(NZ+1,J+1,1))	&  
     	  /INT_JACOB(NZ+1,J))
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1) =U2X(I,NZ+1,J)
         END DO			
      END IF

      RETURN
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|------
SUBROUTINE APPLY_BC_3_LOWER(K)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----
  USE ntypes
  USE Domain
  USE Grid
  ! USE Fft_var
  USE TIME_STEP_VAR
  USE run_variable
  USE ADI_var
  USE mg_vari, ONLY : INIT_FLAG
	
    IMPLICIT NONE

    INTEGER I,K

! C Bottom Wall:
      IF (W_BC_YMIN.EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLY(I,0)=0. 
          MATDY(I,0)=1.
          MATUY(I,0)=0.                   
          VECY(I,0)=W_BC_YMIN_C1

          !MATLY(I,1)=0. 
          !MATDY(I,1)=1.
          !MATUY(I,1)=0.                   
          !VECY(I,1)=W_BC_YMIN_C1 
        END DO
      ELSE IF (W_BC_YMIN.eq.1) THEN
! C Neumann
       IF (IC_TYPE == 7) THEN
        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=-1.
          MATUY(I,0)=1.
          VECY(I,0)=W_BC_YMIN_C1
        END DO
       ELSE
        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=-1.
          MATUY(I,0)=1.
          VECY(I,0)=DY(1)*W_BC_YMIN_C1
         ENDDO
       ENDIF
      ELSE IF (W_BC_YMIN.eq.3) THEN
! C Wall model artificial slip boundary condition
! C Neumann
	IF ( IC_TYPE == 7) THEN
	 IF (WALL_MODEL_TYPE.EQ.2 .OR. WALL_MODEL_TYPE.EQ.1) THEN
!	    DO I=0,NXP
!	      MATLY(I,1)=0.
!	      MATDY(I,1)=1.
!	      MATUY(I,1)=0.
!	      VECY(I,1)=W_BC_LOWER_WALLMODEL(I,K)
!	    ENDDO
    ! C This gridpoint is not USEd here
	    DO I=0,NXP
	      MATLY(I,0)=0.
	      MATDY(I,0)=1.
	      MATUY(I,0)=0.
	      VECY(I,0)=W_BC_LOWER_WALLMODEL(I,K)*COS(ATAN(-CJOB_21(K,0,1)/CJOB_22(K,0,1)))
	    ENDDO
	  ELSE
	   DO I=0,NXP
	      MATLY(I,0)=0.
	      MATDY(I,0)=-1.
	      MATUY(I,0)=1.
	      VECY(I,0)=W_BC_LOWER_WALLMODEL(I,K)
	    ENDDO
    ! C This gridpoint is not USEd here
!	    DO I=0,NXP
!	      MATLY(I,0)=0.
!	      MATDY(I,0)=1.
!	      MATUY(I,0)=0.
!	      VECY(I,0)=0.
!	    ENDDO
	  ENDIF
        ELSE
	  DO I=0,NXP
	    MATLY(I,1)=0.
	    MATDY(I,1)=-1.
	    MATUY(I,1)=1.
	    VECY(I,1)=DY(2)*W_BC_LOWER_WALLMODEL(I,K)
	  ENDDO
  ! C This gridpoint is not used here
	  DO I=0,NXP
	    MATLY(I,0)=0.
	    MATDY(I,0)=1.
	    MATUY(I,0)=0.
	    VECY(I,0)=0.
	  ENDDO
	ENDIF
        
      ElSE IF (W_BC_YMIN .EQ. 6) THEN
! C Periodic
        DO I=0,NXP
          MATLY(I,0)=0. 
          MATDY(I,0)=1.
          MATUY(I,0)=0.                   
          VECY(I,0)=U3X(I,K,NY)
	ENDDO	
      ElSE IF (W_BC_YMIN .EQ. 7) THEN
! C Mixed boundary condition for w at the bottom wall 
       IF ( K .lE. 80 ) THEN
        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=-1.
          MATUY(I,0)=1.
          VECY(I,0)=DY(1)*W_BC_YMIN_C1
        END DO
       ELSE
        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=1.
          MATUY(I,0)=0.
          VECY(I,0)=0.

          MATLY(I,1)=0.
          MATDY(I,1)=1.
          MATUY(I,1)=0.
          VECY(I,1)=W_BC_YMIN_C1
        END DO
       ENDIF
 
      ELSE
        write(*,*) 'WARNING: W_BC_LOWER is of unknown type'
      END IF

  RETURN 
END



! C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_3_UPPER(K)
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
USE ntypes
USE Domain
USE Grid
! USE Fft_var
USE TIME_STEP_VAR
USE run_variable
USE ADI_var
USE mg_vari, ONLY : INIT_FLAG
      
IMPLICIT NONE

      INTEGER I,K

! C Top wall
      IF (W_BC_YMAX.EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=W_BC_YMAX_C1

          MATLY(I,NY)=0.
          MATDY(I,NY)=1.
          MATUY(I,NY)=0.
          VECY(I,NY)=W_BC_YMAX_C1
        END DO
      ELSE IF (W_BC_YMAX.eq.1) THEN
! C Neumann
       IF (IC_TYPE == 7) THEN
        DO I=0,NXP
          MATLY(I,NY+1)=-1.0
          MATDY(I,NY+1)=1.0
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=W_BC_YMAX_C1
        END DO
       ELSE
        DO I=0,NXP
          MATLY(I,NY)=-1.
          MATDY(I,NY)=1.
          MATUY(I,NY)=0.
          VECY(I,NY)=DY(NY)*W_BC_YMAX_C1
        END DO
        DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=0.
        END DO
       ENDIF	
      ELSE IF (W_BC_YMAX.EQ.3) THEN
! C Use wall model artificial boundary condition
! C Neumann
        DO I=0,NXP
          MATLY(I,NY)=-1.
          MATDY(I,NY)=1.
          MATUY(I,NY)=0.
          VECY(I,NY)=DY(NY)*W_BC_UPPER_NWM(I,K)
        END DO
! C This gridpoint is not USEd here
        DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=U3X(I,K,NY)
        END DO
      ELSE IF (W_BC_YMAX .eq. 6) THEN	
! c periodic 
       DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=U3X(I,K,1)
       ENDDO	
      END IF

      RETURN
      END   

! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_3_LEFT(J)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----
USE ntypes
USE Domain
USE Grid
! USE Fft_var
USE TIME_STEP_VAR
USE run_variable
USE ADI_var
USE mg_vari, ONLY : INIT_FLAG
      
IMPLICIT NONE

      INTEGER I,J

! C Bottom Wall:
      IF (W_BC_ZMIN.EQ.0) THEN
! C Dirichlet
       
        DO I=0,NXP
         IF (IC_TYPE .EQ. 5) THEN
          MATLX(I,1)=0.
          MATDX(I,1)=1.
          MATUX(I,1)=0.
          VECX(I,1)= U3_BAR(1,J)

          MATLX(I,2)=0.
          MATDX(I,2)=1.
          MATUX(I,2)=0.
          VECX(I,2)=U3_BAR(1,J)
        ELSEIF (IC_TYPE .EQ. 7) THEN
	  MATLX(I,0)=0. 
          MATDX(I,0)=1.
          MATUX(I,0)=0.                   
          VECX(I,0)=U3_BAR(1,J)
	  
	  MATLX(I,1)=0. 
          MATDX(I,1)=1.
          MATUX(I,1)=0.                   
          VECX(I,1)= U3_BAR(1,J) !W_BC_ZMIN_C1
	ELSE
          MATLX(I,1)=0. 
          MATDX(I,1)=1.
          MATUX(I,1)=0.                   
          VECX(I,1)=W_BC_ZMIN_C1

          MATLX(I,2)=0. 
          MATDX(I,2)=1.
          MATUX(I,2)=0.                   
          VECX(I,2)=W_BC_ZMIN_C1
         ENDIF 
        END DO
      ELSE IF (W_BC_ZMIN.eq.1) THEN
! C Neumann
      IF (IC_TYPE .EQ. 7) THEN 
	DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=-1.
          MATUX(I,0)=1.
          VECX(I,0)=W_BC_ZMIN_C1
        END DO
       ELSE
        DO I=0,NXP
          MATLX(I,1)=0.
          MATDX(I,1)=-1.
          MATUX(I,1)=1.
          VECX(I,1)=DZF(1)*W_BC_ZMIN_C1
        END DO
       ENDIF	
	
      ELSE IF (W_BC_ZMIN.eq.3) THEN
! C Wall model artificial slip boundary condition
! C Neumann
        DO I=0,NXP
          MATLX(I,1)=0.
          MATDX(I,1)=1.
          MATUX(I,1)=0.
          VECX(I,1)=0.
  
          MATLX(I,2)=0.
          MATDX(I,2)=1.
          MATUX(I,2)=0.
          VECX(I,2)=0.
        END DO
! C This gridpoint is not USEd here
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=1.
          MATUX(I,0)=0.
          VECX(I,0)=0.
        END DO
      ELSE IF (W_BC_ZMIN .eq. 6) THEN
! C Periodic
       DO I=0,NXP          
          MATLX(I,0)=0. 
          MATDX(I,0)=1.
          MATUX(I,0)=0.                   
          VECX(I,0)=U3X(I,NZ,J)
        END DO
      ELSE IF (W_BC_ZMIN.EQ.7) THEN
	 DO I=0,NXP
	  MATLX(I,0)=0.
          MATDX(I,0)=(1. - dtc*C_int_le(J)*0.25* &
          ( CJOB_11(0,J,2) + CJOB_11(1,J,2)      &
          + CJOB_11(0,J,1) + CJOB_11(0,J+1,1))	   &
     	  /INT_JACOB(0,J))
          MATUX(I,0)=(dtc*C_int_le(J)*0.25*       &
         ( CJOB_11(0,J,2) + CJOB_11(1,J,2)        &
          + CJOB_11(0,J,1) + CJOB_11(0,J+1,1))	  &
     	  /INT_JACOB(0,J))
          VECX(I,0) =U3X(I,0,J)
         END DO		
      ELSE
        write(*,*) 'WARNING: W_BC_LEFT is of unknown type'
      END IF

! C The following is ONLY a placeholder, this row is USEd for U1 and U2
! c      DO I=0,NXP
! c        MATLX(I,0) = 0.
! c        MATDX(I,0) = 1.
! c        MATUX(I,0) = 0.
! c        VECX(I,0) = 0.
! c      END DO

      RETURN 
      END



! C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_3_RIGHT(J)
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
USE ntypes
USE Domain
USE Grid
! USE Fft_var, ONLY : NX2P
USE TIME_STEP_VAR
USE run_variable
USE ADI_var
USE mg_vari, ONLY : INIT_FLAG
      
IMPLICIT NONE
      INTEGER I,K,J

! C Top wall
      IF (W_BC_ZMAX.EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=W_BC_ZMAX_C1

          MATLX(I,NZ)=0.
          MATDX(I,NZ)=1.
          MATUX(I,NZ)=0.
          VECX(I,NZ)=W_BC_ZMAX_C1
        END DO
      ELSE IF (W_BC_ZMAX.eq.1) THEN
! C Neumann
      IF (IC_TYPE == 7) THEN
        DO I=0,NXP
          MATLX(I,NZ+1)=-1.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=W_BC_ZMAX_C1
        END DO
      ELSE  
	DO I=0,NXP
          MATLX(I,NZ+1)=-1.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=DZF(NZ)*W_BC_ZMAX_C1
        END DO
       ENDIF	
      ELSE IF (W_BC_ZMAX.EQ.4) THEN
! C  outlet Boundary condition
         DO I=0,NXP
           MATLX(I,NZ+1)=0.
           MATDX(I,NZ+1)=1.
           MATUX(I,NZ+1)=0.
           VECX(I,NZ+1)=U3X(I,NZ-1,J)+(U3X(I,NZ,J)-U3X(I,NZ-1,J)) &
                   *DZF(NZ)/DZF(NZ-1)

           MATLX(I,NZ)=-1.d0/DZF(NZ-1)
           MATDX(I,NZ)=1.d0/DZF(NZ)+1.d0/DZF(NZ-1)
           MATUX(I,NZ)=-1.d0/DZF(NZ)
           VECX(I,NZ) = 0.
          END DO 
	  
       ELSE IF (W_BC_ZMAX.EQ.5) THEN
        IF (IC_TYPE ==7 )THEN
	 DO I=0,NXP
	  MATLX(I,NZ+1)=-2.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1) =0.
         END DO
	ELSE
         DO I=0,NXP
          MATLX(I,NZ+1)=-(1.d0/DZF(NZ-1)+1.d0/DZF(NZ))
          MATDX(I,NZ+1)=1.d0/DZF(NZ)
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1) =0.

! c          MATLX(I,NZ)=-(1.d0/DZF(NZ-1)+1.d0/DZF(NZ-2))
! c          MATDX(I,NZ)=1.d0/DZF(NZ-1)
! c          MATUX(I,NZ)=0.
! c          VECX(I,NZ) =0.
        END DO
	ENDIF
      
      ELSE IF (W_BC_ZMAX.EQ.3) THEN
! C Use wall model artificial boundary condition
! C Neumann
        DO I=0,NXP
          MATLX(I,NZ+1)=-1.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=DZF(NZ)*W_BC_UPPER_NWM(I,K)
        END DO
! C This gridpoint is not USEd here
        DO I=0,NXP
          MATLX(I,NZ)=0.
          MATDX(I,NZ)=1.
          MATUX(I,NZ)=0.
          VECX(I,NZ)=0.
        END DO
      ELSE IF (W_BC_ZMAX .eq. 6) THEN
! C periodic condition      	
	DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=U3X(I,1,J)
	ENDDO
      ELSE IF (W_BC_ZMAX.EQ.7) THEN
	 DO I=0,NXP
	  MATLX(I,NZ+1)=-(dtc*C_int(J)*0.25*  &
         ( CJOB_11(NZ+1,J,2) + CJOB_11(NZ+2,J,2) &
          + CJOB_11(NZ+1,J,1) + CJOB_11(NZ+1,J+1,1))&
     	  /INT_JACOB(NZ+1,J))
          MATDX(I,NZ+1)=(1. + dtc*C_int(J)*0.25* &
         ( CJOB_11(NZ+1,J,2) + CJOB_11(NZ+2,J,2) &
          + CJOB_11(NZ+1,J,1) + CJOB_11(NZ+1,J+1,1))	&  
     	  /INT_JACOB(NZ+1,J))
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1) =U3X(I,NZ+1,J)
         END DO		
      END IF

      RETURN
      END
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_VEL_LOWER
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
! C This subroutine is called after initializing the flow
! C It sets the appropriate boundary conditions including ghost cell values
! C  on the velocity field in Fourier space
USE ntypes
USE Domain
USE Grid
!USE Fft_var, ONLY : NX2P
USE TIME_STEP_VAR
USE run_variable
USE ADI_var
USE mg_vari, ONLY : INIT_FLAG
      
IMPLICIT NONE
   


      RETURN
      END

     
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE  APPLY_BC_VEL_UPPER
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! C This subroutine is called after initializing the flow
! C It sets the appropriate boundary conditions including ghost cell values
! C  on the velocity field in Fourier space
USE ntypes
USE Domain
USE Grid
!USE Fft_var, ONLY : NX2P
USE TIME_STEP_VAR
USE run_variable
USE mg_vari, ONLY : INIT_FLAG
      
IMPLICIT NONE

   
      RETURN
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_VEL_LEFT
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
! C This subroutine is called after initializing the flow
! C It sets the appropriate boundary conditions including ghost cell values
! C  on the velocity field in Fourier space
USE ntypes
USE Domain
USE Grid
!USE Fft_var, ONLY : NX2P
USE TIME_STEP_VAR
USE run_variable
USE mg_vari, ONLY : INIT_FLAG
      
IMPLICIT NONE


      RETURN
      END

    
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE  APPLY_BC_VEL_RIGHT
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! C This subroutine is called after initializing the flow
! C It sets the appropriate boundary conditions including ghost cell values
! C  on the velocity field in Fourier space
USE ntypes
USE Domain
USE Grid
!USE Fft_var, ONLY : NX2P
USE TIME_STEP_VAR
USE run_variable
USE mg_vari, ONLY : INIT_FLAG
      
IMPLICIT NONE

     
   
      RETURN
      END
