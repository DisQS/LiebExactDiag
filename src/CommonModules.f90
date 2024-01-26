!-----------------------------------------------------------------
!
! LEDdiag
!
!-----------------------------------------------------------------
! exact diagonalization of 2D and 3D extended Lieb models
! see https://github.com/DisQS/LiebExactDiag
!-----------------------------------------------------------------

MODULE MyNumbers     
  IMPLICIT NONE

  INTEGER, PARAMETER :: IKIND = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(15,307)

  REAL(KIND=RKIND) :: PI, ONEPLS, ONEMNS

  REAL(KIND=RKIND), PARAMETER :: ZERO = 0.0, ONE = 1.0 ,TWO = 2.0, THREE = 3.0, FOUR = 4.0
  COMPLEX(KIND=RKIND), PARAMETER :: CZERO = (0.0d0,0.0d0), CONE = (1.0d0,0.0d0), &
       CIMAGONE= (0.0d0,1.0d0)

  REAL (KIND=RKIND), PARAMETER :: HALF = 0.5D0, QUARTER = 0.25D0, EIGHTH = 0.125D0

  REAL(KIND=RKIND) :: TINY= 1.0D-9

CONTAINS
  SUBROUTINE INIT_NUMBERS
    PI = 4.0D0* ATAN(1.0D0)
    ONEMNS = SQRT(EPSILON(ONEMNS))
    ONEPLS = ONE + ONEMNS
    ONEMNS = ONE - ONEMNS
  END SUBROUTINE INIT_NUMBERS

  FUNCTION ARG(X,Y)
    
    REAL(KIND=RKIND) ARG, X, Y
    
    IF( X > 0. ) THEN 
       ARG= ATAN(Y/X)
    ELSE IF ( (X == 0.) .and. (Y > 0. )) THEN 
       ARG = PI/2.0D0
    ELSE IF ( (X == 0.) .and. (Y < 0. )) THEN 
       ARG = -PI/2.0D0
    ELSE IF ( (X < 0. ) .and. (Y >= 0.)) THEN 
       ARG = PI + ATAN(Y/X)
    ELSE IF ( (X < 0. ) .and. (Y < 0. )) THEN 
       ARG = -PI + ATAN(Y/X)
    ENDIF
    
    RETURN
  END FUNCTION ARG

END MODULE MyNumbers



