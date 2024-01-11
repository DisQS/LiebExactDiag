!-----------------------------------------------------------------
!
! LEDdiag
!
!-----------------------------------------------------------------
! exact diagonalization of 2D and 3D extended Lieb models
! see https://github.com/DisQS/LiebExactDiag
!-----------------------------------------------------------------

!!--------------------------------------------------------------------
MODULE CConstants
  CHARACTER*18, PARAMETER :: RStr= "$Revision: 1.3 $ "
  CHARACTER*30, PARAMETER :: DStr= "$Date: 2007/09/25 20:45:08 $ "
  CHARACTER*16, PARAMETER :: AStr= "$Author: phrfar $ "
END MODULE CConstants

!! MAXGamma needs to be equal to MAXWidth, as we need to find ALL
!! Lyapunov exponents, so do not change!
MODULE IConstants
  INTEGER, PARAMETER :: MAXWidth= 1000, MAXGamma= MAXWidth, MAXIter=2147483646
  INTEGER, PARAMETER :: MAXKeepFlag= 3, MAXWriteFlag= 4, MAXFluxFlag= 3, MAXRNGFlag=2
  INTEGER, PARAMETER :: MAXSortFlag=1, MAXBCFlag=2, MAXLevelFlag=0, MAXConvFlag=3
  INTEGER, PARAMETER :: MAXStripeFlag=3, MINDimenFlag=2,MAXDimenFlag=3
  INTEGER, PARAMETER :: MAXFiles= 5, MINIter=3
END MODULE IConstants

!!--------------------------------------------------------------------
MODULE IPara
  USE MyNumbers
  INTEGER(KIND=IKIND) :: ISeed, NSeed
  INTEGER(KIND=IKIND) :: Dim, Nx
  INTEGER(KIND=IKIND) :: Width0, Width1, dWidth, IKeepFlag, IWriteFlag, ISortFlag
  INTEGER(KIND=IKIND) :: IFluxFlag, IBCFlag, IRNGFlag, ILevelFlag, IConvFlag
  INTEGER(KIND=IKIND) :: IStripeFlag, IDimenFlag, IStateFlag
END MODULE IPara

!!--------------------------------------------------------------------
MODULE DPara
  USE MyNumbers
  REAL(KIND=RKIND) :: CubeDis0,CubeDis1,dCubeDis
  REAL(KIND=RKIND) :: CubeConPot0,CubeConPot1,dCubeConPot
  REAL(KIND=RKIND) :: CubeDis, LiebDis
  REAL(KIND=RKIND) :: CubeConPot, LiebConPot
  REAL(KIND=RKIND) :: CubeShiftOD, CubeDisOD, CubeDisOD0,CubeDisOD1,dCubeDisOD
  REAL(KIND=RKIND) :: LiebShiftOD, LiebDisOD, LiebDisOD0,LiebDisOD1,dLiebDisOD
  REAL(KIND=RKIND) :: flux, flux0,flux1,dflux
  REAL(KIND=RKIND) :: Kappa, MagFlux
  REAL(KIND=RKIND) :: MyEpsilon

  REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE :: PsiA, PsiB
END MODULE DPara

!!--------------------------------------------------------------------
!!      Input- and Outputchannels
MODULE IChannels
  USE MyNumbers
  INTEGER(KIND=IKIND), PARAMETER :: IChInp= 40, IChOut= 41, IChOutGam= 42, &
       IChOutPsi= 43, IChOutRHO= 44, IChOutAvgRHO= 45, &
       IChOutAvgRHO1= 46, IChOutAvgRHOL= 47, &
       ICHtmp= 48, IChOutHAV= 49, &
       IChEVal=50, IChEVec=51
END MODULE IChannels






