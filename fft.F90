!C******************************************************************************|
!C fft.f, the FFT package for diablo.                               VERSION 1.1
!C
!C This file isolates all calls to the FFTW package (available at: www.fftw.org)
!C These wrapper routines were written by T. Bewley (spring 2001).
!C******************************************************************************|

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!C The arrangement of the significant real numbers in the arrays (denoted by +)
!C in physical space, in Fourier space, and in Fourier space after packing are
!C shown below for the 2D (X-Z) plane.  The third direction (Y) is handled in
!C an identical matter as the Z direction shown here.
!C
!C      oooooooooooooooooo         oooooooooooooooooo         oooooooooooooooooo
!C      oooooooooooooooooo         oooooooooooooooooo         oooooooooooooooooo
!C NZ-1 ++++++++++++++++oo     -1  ++++++++++++oooooo         oooooooooooooooooo
!C      ++++++++++++++++oo     -2  ++++++++++++oooooo         oooooooooooooooooo
!C      ++++++++++++++++oo     -3  ++++++++++++oooooo         oooooooooooooooooo
!C      ++++++++++++++++oo         ++++++++++++oooooo         oooooooooooooooooo
!C      ++++++++++++++++oo    -NKZ ++++++++++++oooooo         oooooooooooooooooo
!C      ++++++++++++++++oo         oooooooooooooooooo     -1  ++++++++++++oooooo
!C      ++++++++++++++++oo         oooooooooooooooooo     -2  ++++++++++++oooooo
!C      ++++++++++++++++oo         oooooooooooooooooo     -3  ++++++++++++oooooo
!C      ++++++++++++++++oo         oooooooooooooooooo         ++++++++++++oooooo
!C      ++++++++++++++++oo         oooooooooooooooooo    -NKZ ++++++++++++oooooo
!C      ++++++++++++++++oo     NKZ ++++++++++++oooooo     NKZ ++++++++++++oooooo
!C      ++++++++++++++++oo         ++++++++++++oooooo         ++++++++++++oooooo
!C   3  ++++++++++++++++oo      3  ++++++++++++oooooo      3  ++++++++++++oooooo
!C   2  ++++++++++++++++oo      2  ++++++++++++oooooo      2  ++++++++++++oooooo
!C   1  ++++++++++++++++oo      1  ++++++++++++oooooo      1  ++++++++++++oooooo
!C   0  ++++++++++++++++oo      0  +o++++++++++oooooo      0  +o++++++++++oooooo
!C      ^^^^           ^           ^ ^ ^     ^                ^ ^ ^     ^
!C      0123           NX-1        0 1 2     NKX              0 1 2     NKX
!C
!C       PHYSICAL SPACE              FOURIER SPACE         FOURIER SPACE (PACKED)
!C
!C After the Real->Fourier transform, the significant coefficients are put next
!C to each other in the array, so a loop such as
!C
!C        DO K=0,TNKZ           [where TNKZ = 2*NKZ = 2*(NZ/3) ]
!C          DO I=0,NKX          [where  NKX = NX/3             ]
!C            CP(I,K,J)= ...
!C          END DO
!C        END DO
!C
!C includes all the Fourier coefficients of interest.  The subsequent loops in
!C Fourier space just work on these coefficients in the matrix.
!C  
!C Before a Fourier->Real transform, the significant coefficients are unpacked
!C and the higher wavenumbers are SET TO ZERO before the inverse transform.
!C This has the effect of doing the required dealiasing.
!C
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE INIT_FFT
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use mpi_var
!use run_variable

      IMPLICIT NONE

      INTEGER I

      INTEGER         FFTW_FORWARD,      FFTW_BACKWARD,&
                     FFTW_ESTIMATE,     FFTW_MEASURE,&
                     FFTW_OUT_OF_PLACE, FFTW_IN_PLACE,&
                     FFTW_USE_WISDOM,   FFTW_THREADSAFE
      PARAMETER(      FFTW_FORWARD=-1,      FFTW_BACKWARD=1, &
                     FFTW_ESTIMATE=0,      FFTW_MEASURE=1, &
                     FFTW_OUT_OF_PLACE=0,  FFTW_IN_PLACE=8,&
                     FFTW_USE_WISDOM=16,   FFTW_THREADSAFE=128 )

      WRITE(6,*) 'Initializing FFTW package.'
      
      PI = 4. * ATAN(1.0)
      CI = CMPLX(0.0,1.0)
      EPS = 0.000000001


      CALL RFFTWND_F77_CREATE_PLAN (FFTW_X_TO_F_PLAN, 1, NX, &
				    FFTW_FORWARD,  FFTW_MEASURE + FFTW_IN_PLACE )
      
      CALL RFFTWND_F77_CREATE_PLAN(FFTW_X_TO_P_PLAN, 1, NX, &
				    FFTW_BACKWARD, FFTW_MEASURE + FFTW_IN_PLACE )
!        NKX=NX/3
      RNX=1.0*NX
      DO I=0,NKX
	KX(I) = I * (2.*PI)/LX
	KX2(I) = KX(I)*KX(I)
	CIKX(I) = CI*KX(I)
      END DO

      DO I=0,NX2P
	KXP(I) = (I+rank*(NX2P+1)) * (2.*PI)/LX
	KX2P(I) = KXP(I) * KXP(I)
	CIKXP(I) = CI * KXP(I)
      END DO

    WRITE(6,*)
    WRITE(6,*) 'FFTW package Successfully initialized along Span-Wise direction'


  RETURN
END

!C******************************************************************************|
!C-------------> The transform routines for the duct flow follow. <-------------|
!C******************************************************************************|

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE FFT_X_TO_FOURIER(U,CU,JMIN,JMAX,KMIN,KMAX)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!C This routine transforms (in 1 direction) planes JMIN-JMAX to Fourier space.
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var
  USE TIME_STEP_VAR

    IMPLICIT NONE

      INTEGER JMIN, JMAX, KMIN, KMAX, I, J, K
      REAL*8     U (0:NX+1,0:NZ+1,0:NY+1)
      COMPLEX*16 CU(0:NX/2,0:NZ+1,0:NY+1)

!C Looping over the planes of interest, simply perform a real -> complex
!C transform in place in the big storage array, scaling appropriately.

      DO J=JMIN,JMAX
	CALL RFFTWND_F77_REAL_TO_COMPLEX( FFTW_X_TO_F_PLAN,(KMAX-KMIN+1),&
					   U(0,KMIN,J), 1, NX+2, CU(0,KMIN,J), 1, NX/2+1) 

	DO K=KMIN,KMAX
	  DO I=0,NKX
	    CU(I,K,J) = CU(I,K,J)/RNX
	  ENDDO
	ENDDO

      ENDDO
    
  RETURN
END

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE FFT_X_TO_PHYSICAL(CU,U,JMIN,JMAX,KMIN,KMAX)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!C This routine transforms (in 1 direction) planes JMIN-JMAX to physical space.
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
USE ntypes
USE Domain
USE Grid
USE Fft_var
USE TIME_STEP_VAR
!use run_variable

  IMPLICIT NONE

    INTEGER JMIN, JMAX, KMIN, KMAX, I, J, K
    REAL*8     U (0:NX+1,0:NZ+1,0:NY+1)
    COMPLEX*16 CU(0:NX/2,0:NZ+1,0:NY+1)

!C Looping over the planes of interest, simply set the higher wavenumbers to
!C zero and then perform a complex -> real transform in place in the big
!C storage array.

    DO J=JMIN,JMAX
      DO K=KMIN,KMAX
	DO I=NKX+1,NX/2
	  CU(I,K,J)=0.d0
	ENDDO
      ENDDO
      CALL RFFTWND_F77_COMPLEX_TO_REAL( FFTW_X_TO_P_PLAN,(KMAX-KMIN+1),&
					 CU(0,KMIN,J), 1, NX/2+1, U(0,KMIN,J), 1, NX+2) 
    ENDDO

  RETURN
END


!C******************************************************************************|

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE FFT_X_TO_FOURIER_OP(U,CU,JMIN,JMAX,KMIN,KMAX)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!C This routine transforms (in 1 direction) planes JMIN-JMAX to Fourier space.
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
USE ntypes
USE Domain
USE Grid
USE Fft_var
USE TIME_STEP_VAR
!use run_variable

  IMPLICIT NONE

    INTEGER JMIN, JMAX, KMIN, KMAX, I, J, K
    REAL*8     U (0:NX+1,0:NZP,0:NY+1)
    COMPLEX*16 CU(0:NX/2,0:NZP,0:NY+1)

!C Looping over the planes of interest, simply perform a real -> complex
!C transform in place in the big storage array, scaling appropriately.

    DO J=JMIN,JMAX
      CALL RFFTWND_F77_REAL_TO_COMPLEX( FFTW_X_TO_F_PLAN,(KMAX-KMIN+1),&
				      U(0,KMIN,J), 1, NX+2, CU(0,KMIN,J), 1, NX/2+1)
      DO K=KMIN,KMAX
	DO I=0,NKX
	  CU(I,K,J)=CU(I,K,J)/RNX
	ENDDO
      ENDDO
    ENDDO 

  RETURN
END

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE FFT_X_TO_PHYSICAL_OP(CU,U,JMIN,JMAX,KMIN,KMAX)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!C This routine transforms (in 1 direction) planes JMIN-JMAX to physical space.
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  USE ntypes
  USE Domain
  USE Grid
  USE Fft_var
  USE TIME_STEP_VAR

  IMPLICIT NONE

    INTEGER JMIN, JMAX, KMIN, KMAX, I, J, K
    REAL*8     U (0:NX+1,0:NZP,0:NY+1)
    COMPLEX*16 CU(0:NX/2,0:NZP,0:NY+1)

!C Looping over the planes of interest, simply set the higher wavenumbers to
!C zero and then perform a complex -> real transform in place in the big
!C storage array.

    DO J=JMIN,JMAX
      DO K=KMIN,KMAX
	DO I=NKX+1,NX/2
	  CU(I,K,J)=0.d0
	ENDDO
      ENDDO
      
      CALL RFFTWND_F77_COMPLEX_TO_REAL( FFTW_X_TO_P_PLAN,(KMAX-KMIN+1),&
					 CU(0,KMIN,J), 1, NX/2+1, U(0,KMIN,J), 1, NX+2)
    ENDDO


  RETURN
END



