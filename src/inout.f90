!-----------------------------------------------------------------
!
! LEDdiag
!
!-----------------------------------------------------------------
! exact diagonalization of 2D and 3D extended Lieb models
! see https://github.com/DisQS/LiebExactDiag
!-----------------------------------------------------------------

! --------------------------------------------------------------------
! Input:
!
! IErr	error code
!----------------------------------------------------------------------

SUBROUTINE Input(IErr)

  USE MyNumbers
  
  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara
  
  USE IChannels
  
  INTEGER IErr, ILine
!  REAL(KIND=RKIND) RIter
  
  !	PRINT*,"DBG: Input()"
  
  IErr = 0
  ILine= 0
  
!!$  OPEN(UNIT= IChInp, ERR=120, FILE= "LEDdiag.inp",STATUS= 'OLD')
  OPEN(UNIT= IChInp, ERR=120, FILE= "/dev/stdin",STATUS= 'OLD')

  ILine= ILine+1
  READ(IChInp,10,ERR=20) ISeed
  PRINT*,"ISeed        = ",ISeed

  ILine= ILine+1
  READ(IChInp,10,ERR=20) NSeed
  PRINT*,"NSeed        = ",NSeed

  ILine= ILine+1
  READ(IChInp,10,ERR=20) Dim
  PRINT*,"Dim          = ",Dim

  ILine= ILine+1
  READ(IChInp,10,ERR=20) Nx
  PRINT*,"Nx           = ",Nx
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IBCFlag
  PRINT*,"IBCFlag      = ", IBCFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IRNGFlag
  PRINT*,"IRNGFlag     = ", IRNGFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IKeepFlag
  PRINT*,"IKeepFlag    = ", IKeepFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IWriteFlag
  PRINT*,"IWriteFlag   = ", IWriteFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IStateFlag
  PRINT*,"IWriteFlag   = ", IStateFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IFluxFlag
  PRINT*,"IFluxFlag    = ", IFluxFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) Width0
  PRINT*,"Width0       = ",Width0
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) Width1
  PRINT*,"Width1       = ", Width1
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) dWidth
  PRINT*,"dWidth       = ", dWidth

  ILine= ILine+1
  READ(IChInp,15,ERR=20) CubeConPot0
  PRINT*,"CubeConPot0  = ", CubeConPot0

  ILine= ILine+1
  READ(IChInp,15,ERR=20) CubeConPot1
  PRINT*,"CubeConPot1  = ", CubeConPot1

  ILine= ILine+1
  READ(IChInp,15,ERR=20) dCubeConPot
  PRINT*,"dCubeConPot  = ", dCubeConPot

  ILine= ILine+1
  READ(IChInp,15,ERR=20) CubeDis0
  PRINT*,"CubeDis0     = ", CubeDis0
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20) CubeDis1
  PRINT*,"CubeDis1     = ", CubeDis1
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20) dCubeDis
  PRINT*,"dCubeDis     = ", dCubeDis

  ILine= ILine+1
  READ(IChInp,15,ERR=20) LiebConPot
  PRINT*,"LiebConPot   = ", LiebConPot

  ILine= ILine+1
  READ(IChInp,15,ERR=20) LiebDis
  PRINT*,"LiebDis      = ", LiebDis
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20) OffDShift0
  PRINT*,"OffDShift0   = ", OffDShift0
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20) OffDShift1
  PRINT*,"OffDShift1   = ", OffDShift1
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20) dOffDShift
  PRINT*,"dOffDShift   = ", dOffDShift
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20) OffDDis0
  PRINT*,"OffDDis0     = ", OffDDis0
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20) OffDDis1
  PRINT*,"OffDDis1     = ", OffDDis1
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20) dOffDDis
  PRINT*,"dOffDDis     = ", dOffDDis
  
10 FORMAT(16X,I15.1)
  ! 10	FORMAT("IMAXIteration= ",I15.1)
15 FORMAT(16X,F18.9)
  ! 15     FORMAT("IMAXIteration= ",F18.9)

  ! check the parameters for validity
  
  IF(IWriteFlag.GE.2) THEN
     PRINT*,"ISeed          = ", ISeed
     PRINT*,"NSeed          = ", NSeed
     PRINT*,"Dim            = ", Dim
     PRINT*,"Nx             = ", Nx
     PRINT*,"IBCFlag        = ", IBCFlag
     PRINT*,"IRNGFlag       = ", IRNGFlag
     PRINT*,"IKeepFlag      = ", IKeepFlag
     PRINT*,"IWriteFlag     = ", IWriteFlag
     PRINT*,"IStateFlag     = ", IStateFlag
     PRINT*,"Width0         = ", Width0
     PRINT*,"Width1         = ", Width1
     PRINT*,"dWidth         = ", dWidth
     PRINT*,"0 CubeConPot0  = ", CubeConPot0
     PRINT*,"0 CubeConPot1  = ", CubeConPot1
     PRINT*,"0 dCubeConPot  = ", dCubeConPot
     PRINT*,"1 CubeDis0     = ", CubeDis0
     PRINT*,"1 CubeDis1     = ", CubeDis1
     PRINT*,"1 dCubeDis     = ", dCubeDis
     PRINT*,"2 LiebConPot   = ", LiebConPot
     PRINT*,"3 LiebDis      = ", LiebDis
     PRINT*,"4 OffDShift0   = ", OffDShift0
     PRINT*,"4 OffDShift1   = ", OffDShift1
     PRINT*,"4 dOffDShift   = ", dOffDShift
     PRINT*,"5 OffDDis0     = ", OffDDis0
     PRINT*,"5 OffDDis1     = ", OffDDis1
     PRINT*,"5 dOffDDis     = ", dOffDDis

  ENDIF

  CLOSE(IChInp)
  RETURN

120 PRINT*,"Input(): ERR in OPEN(), line=", ILine  
  IErr= 1
  RETURN

20 PRINT*,"Input(): ERR in READ(), line=", ILine  
  IErr= 1
  RETURN
  
END SUBROUTINE Input

! --------------------------------------------------------------------
! GetFileName
!
! IErr	error code
!----------------------------------------------------------------------

!!$FUNCTION GetFileName(fnamestr,vdata,n)
!!$!write(*,*) trim(GetFileName('/Output/QH_L(I4)_NL(I1)_NS(I1)_B(F5.3)_C(F4.2)_S(I5).dat',&
!!$!     (/L_INPUT,NLevels,NSpins,BField,Coul,Seed/),6))
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(in)         :: n
!!$  REAL(8),INTENT(in)         :: vdata(n)
!!$  CHARACTER(len=*),INTENT(in) :: fnamestr
!!$
!!$  CHARACTER(len=*)   :: GetFileName
!!$  CHARACTER(len=500) :: str,fstr,vstr
!!$  
!!$  INTEGER           :: i
!!$
!!$  str = fnamestr
!!$  
!!$  DO i=1,n
!!$     
!!$     fstr = str(INDEX(str,'('):INDEX(str,')'))
!!$     IF(fstr(2:2).EQ.'I')THEN
!!$        WRITE(vstr,fstr) INT(vdata(i))
!!$     ELSE
!!$        WRITE(vstr,fstr) vdata(i)
!!$     END IF
!!$  
!!$     str = str(1:INDEX(str,'(')-1) // TRIM(ADJUSTL(vstr)) // str(INDEX(str,')')+1:LEN(str)) 
!!$     
!!$  END DO
!!$  
!!$  GetFileName = TRIM(ADJUSTL(str))
!!$
!!$  RETURN
!!$END FUNCTION GetFileName

!--------------------------------------------------------------------
! MakeMiddleName
!
! IErr	error code

FUNCTION MakeMiddleName(IWidth, IErr)
!write(*,*) trim(GetFileName('/Output/QH_L(I4)_NL(I1)_NS(I1)_B(F5.3)_C(F4.2)_S(I5).dat',&
!     (/L_INPUT,NLevels,NSpins,BField,Coul,Seed/),6))

  USE MyNumbers
  USE IPara
  USE DPara

  IMPLICIT NONE

  CHARACTER*200          :: MakeMiddleName
  INTEGER, INTENT(in)    :: IWidth
  INTEGER, INTENT(inout) :: IErr
  CHARACTER*200 middlenamestr, mnstr1, mnstr2, mnstr3
  CHARACTER*1 SymbolCP,symbolLP
  
  PRINT*,"DBG: MakeMiddleName()"

  IErr= 0

  IF(CubeConPot.GE.0.0D0) THEN
     SymbolCP="+"
  ELSE
     SymbolCP="-"
  END IF

  IF(LiebConPot.GE.0.0D0) THEN
     SymbolLP="+"
  ELSE
     SymbolLP="-"
  END IF

  !   WRITE out the input parameter
  WRITE(mnstr1, &
   '(A1,I1,I1,A2,I4.4,A3)') &
       "L", Dim, Nx, &
       "_M", IWidth, &
       "_Es"
  
  WRITE(mnstr2, &
   '(A3,A1,I6.6,A3,I6.6,A3,A1,I6.6,A3,I6.6)') &
       "_CP", SymbolCP,NINT(100.*ABS(CubeConPot)), &
       "_CD", NINT(100.*ABS(CubeDis)), &
       "_LP", SymbolLP,NINT(100.*ABS(LiebConPot)), &
       "_LD", NINT(100.*ABS(LiebDis)) 
  
  WRITE(mnstr3, &
   '(A3,I6.6,A3,I6.6)') &
       "_oS", NINT(100.*ABS(OffDShift)), &
       "_oD", NINT(100.*ABS(OffDDis))
  
!!$  PRINT*, TRIM(mnstr1)
!!$  PRINT*, TRIM(mnstr2)
!!$  PRINT*, TRIM(mnstr3)
  middlenamestr= TRIM(mnstr1) // TRIM(mnstr2) // TRIM(mnstr3)
  
  MakeMiddleName = TRIM(ADJUSTL(middlenamestr))

  RETURN
END FUNCTION MakeMiddleName

!--------------------------------------------------------------------
! MakeDataDir:
!
! IErr	error code

!!$SUBROUTINE MakeDataDir(Dim, Nx, Width, CubeDis, LiebDis, CubeConPot, LiebConPot, str)
SUBROUTINE MakeDataDir(IWidth, dirname, IErr)

  USE MyNumbers
  USE IPara
  USE DPara

  INTEGER IWidth, IErr
  CHARACTER(len=*) dirname
  !CHARACTER*200 MakeDataDir

  LOGICAL*4  ierr1
  EXTERNAL SYSTEM
  
  PRINT*, "MakeDataDir(): checking for ", TRIM(dirname)

#ifdef ifort
  INQUIRE(directory=TRIM(ADJUSTL(dirname)), Exist=ierr1) ! ifort
#else
  INQUIRE(file=TRIM(ADJUSTL(dirname)), Exist=ierr1) ! gfortran
#endif
  IF(ierr1)THEN
     PRINT*,"MakeDataDir(): using existing directory ", TRIM(ADJUSTL(dirname))
     !WRITE(*,'(/)')
  ELSE
     PRINT*,"MakeDataDir(): creating NEW directory ", TRIM(ADJUSTL(dirname))
     !WRITE(*,'(/)')
     CALL System("mkdir -p "//TRIM(ADJUSTL(dirname)) )
  END IF

  RETURN
  
END SUBROUTINE MakeDataDir

!--------------------------------------------------------------------
! CheckOutput:
!
! IErr	error code

SUBROUTINE CheckOutput( IWidth, PreSeed, DirName, MiddleName, IErr )
  !Dim, Nx, IWidth, CubeDis, LiebDis, CubeConPot, LiebConPot, PreSeed, str, IErr )

  USE MyNumbers 
  USE IChannels
  USE IPara
  USE DPara

  INTEGER(KIND=IKIND) IWidth, PreSeed, IErr
  CHARACTER(len=*) DirName, MiddleName
  
  CHARACTER*50 Prefix, Postfix
  CHARACTER*200 FileName

  PRINT*,"DBG: CheckOutput()"
!  IErr= 0

  WRITE(Prefix, '(A5)') "Eval_"
  WRITE(Postfix, '(A2,I5.5,A4)') "-c",  PreSeed, ".raw" 
  !PRINT*, Prefix, Postfix
  
  FileName= TRIM(Prefix)//TRIM(MiddleName)//TRIM(Postfix)
  !PRINT*, FileName
  
  OPEN(UNIT= IChOut, ERR= 10, STATUS= 'NEW', FILE= TRIM(ADJUSTL(DirName))//"/"//FileName)
  !PRINT*, "CheckOutput(): ", TRIM(FileName), "DOES NOT exist -- proceeding!"
  WRITE(*,'(A16,A54,A31)') "CheckOutput(): ", &
       TRIM(FileName), " DOES NOT exist -- proceeding!"

  IErr= 0
  
20 CLOSE(UNIT= IChOut, ERR= 100)
  
  RETURN
  
10 WRITE(*,'(A16,A54,A22)') "CheckOutput(): ", &
        TRIM(FileName), " exists -- skipping!"

  !PRINT*, "CheckOutput(): ", TRIM(FileName), " exists -- skipping!"
  IErr= 2
  GOTO 20
  
  !  ERR in CLOSE detected
100 &
  PRINT*,"CheckOutput(): ERR in CLOSE()"
  IErr= 1
  RETURN
  
END SUBROUTINE CheckOutput


!--------------------------------------------------------------------
! WriteOutputEVal:
!
! IErr	error code

SUBROUTINE WriteOutputEVal(NEVals,EIGS, IWidth,PreSeed, DirName,MiddleName, IErr)

  USE MyNumbers
  USE IChannels
  USE DPara
  USE IPara

  INTEGER(KIND=IKIND) NEVals, IWidth, IErr, PreSeed, i
  REAL(KIND=RKIND) EIGS(NEVals)

  CHARACTER(len=*) DirName, MiddleName
  CHARACTER*50 Prefix, Postfix
  CHARACTER*200 FileName

!!$  PRINT*,"DBG: WriteOutputEVal()"
  IErr= 0

  WRITE(Prefix, '(A5)') "Eval_"
  WRITE(Postfix, '(A2,I5.5,A4)') "-c",  PreSeed, ".raw" 
  !PRINT*, Prefix, Postfix
  
  FileName= TRIM(Prefix)//TRIM(MiddleName)//TRIM(Postfix)
  !PRINT*, FileName

  IF(IWriteFlag.GE.2) THEN
     PRINT*, "WriteOutputEVal(): ", FileName
  ENDIF
  
  OPEN(UNIT= IChEVal, ERR= 10, STATUS='UNKNOWN', FILE= TRIM(ADJUSTL(DirName))//"/"//FileName)
  
  IF(NEVals .GT. 0)THEN
     DO i=1,NEVals
        WRITE(IChEVal, FMT=15, ERR=20) EIGS(i)
15      FORMAT(f30.20)
     ENDDO
  END IF
  
  CLOSE(UNIT=IChEVal, ERR= 30)
  
  RETURN
  
  !	error in OPEN detected
10 PRINT*, "WriteOutputEVal(): ERR in OPEN()"
  IErr= 1
  RETURN
  
  !	error in WRITE detected
20 PRINT*,"WriteOutputEVal(): ERR in WRITE()"
  IErr= 1
  RETURN

  ! ERR in CLOSE detected
30 PRINT*,"OutputEVal(): ERR in CLOSE()"
  IErr= 1
  RETURN
  
END SUBROUTINE WriteOutputEVal


!--------------------------------------------------------------------
! WriteOutputEVec:
!
! IErr	error code

SUBROUTINE WriteOutputEVec( Inum, NEVals, Lsize, VECS, VECS_size, &
     IWidth, PreSeed, DirName, MiddleName, IErr)

  USE MyNumbers
  USE IChannels
  USE DPara
  USE IPara

  INTEGER(KIND=IKIND) Inum, PreSeed, ISSeed, IWidth, IErr, Lsize, VECS_size, NEVals, i

  REAL(KIND=RKIND) VECS(VECS_size,VECS_size)

  CHARACTER(len=*) DirName, MiddleName
  CHARACTER*50 Prefix, Postfix
  CHARACTER*200 FileName
  
!!$  PRINT*,"DBG: WriteOutputEvec()"
  IErr= 0

  WRITE(Prefix, '(A5)') "Evec_"
  WRITE(Postfix, '(A2,I5.5A2,I5.5,A4)') "-c",  PreSeed, "-N", Inum, ".raw" 
  !PRINT*, Prefix, Postfix
  
  FileName= TRIM(Prefix)//TRIM(MiddleName)//TRIM(Postfix)
  !PRINT*, FileName

  IF(IWriteFlag.GE.2) THEN
     PRINT*, "WriteOutputEVec(): ", FileName
  END IF

  OPEN(UNIT= IChEVec, ERR= 40, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(DirName))//"/"//FileName)

  !DO i= 1+( Lsize*( Inum -1) ), 1+( Lsize*( Inum -1) ) + Lsize
  DO i=1,LSize
     WRITE(UNIT=IChEVec, FMT=45, ERR=50) VECS(Inum,i)
  ENDDO

  CLOSE(UNIT= IChEVec, ERR= 60)
  
  RETURN

45 FORMAT(f30.20)

  !	error in OPEN detected
40 PRINT*,"WriteOutputEVec(): ERR in OPEN()"
  IErr= 1
  RETURN

  !	error in WRITE detected
50 PRINT*,"WriteOutputEVec(): ERR in WRITE()"
  IErr= 1
  RETURN

  ! ERR in CLOSE detected
60 PRINT*,"OutputEVec(): ERR in CLOSE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteOutputEVec

!--------------------------------------------------------------------
! WriteOutputEVecBULK:
!
! IErr	error code

SUBROUTINE WriteOutputEVecBULK(NEVals, Lsize, VECS, VECS_size, &
     IWidth, PreSeed, DirName, Middlename, IErr)

  USE MyNumbers
  USE IChannels
  USE DPara
  USE IPara

  INTEGER(KIND=IKIND) PreSeed, ISSeed, IWidth, IErr, Lsize, VECS_size, NEVals, i,j

  REAL(KIND=RKIND) VECS(VECS_size,VECS_size)

  CHARACTER(len=*) DirName, MiddleName
  CHARACTER*50 Prefix, Postfix
  CHARACTER*200 FileName

!!$  PRINT*,"DBG: WriteOutputEvecBULK()"

  IErr= 0

  WRITE(Prefix, '(A6)') "EvecB_"
  WRITE(Postfix, '(A2,I5.5,A4)') "-c",  PreSeed, ".raw" 
  !PRINT*, Prefix, Postfix
  
  FileName= TRIM(Prefix)//TRIM(MiddleName)//TRIM(Postfix)
  !PRINT*, FileName

  IF(IWriteFlag.GE.2) THEN  
     PRINT*, "WriteOutputEVecBULK(): ", FileName
  END IF
  
  OPEN(UNIT= IChEVec, ERR= 40, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(DirName))//"/"//FileName)

  !DO i= 1+( Lsize*( Inum -1) ), 1+( Lsize*( Inum -1) ) + Lsize
  DO j=1,NEVals
     DO i=1,LSize
        WRITE(UNIT=IChEVec, FMT=45, ERR=50) j, i, VECS(j,i)
     ENDDO
  END DO

  CLOSE(UNIT= IChEVec, ERR= 60)
  
  RETURN

45 FORMAT(I6,I6,f30.20)

  !	error in OPEN detected
40 PRINT*,"WriteOutputEVec(): ERR in OPEN()"
  IErr= 1
  RETURN

  !	error in WRITE detected
50 PRINT*,"WriteOutputEVec(): ERR in WRITE()"
  IErr= 1
  RETURN

  ! ERR in CLOSE detected
60 PRINT*,"OutputEVec(): ERR in CLOSE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteOutputEVecBULK

!--------------------------------------------------------------------
! WriteOutputEVecProj:
!
! IErr	error code

SUBROUTINE WriteOutputEVecProj( Inum, NEVals, &
     EIGS, LSize, &
     CubeProb, CubePart, Cube_size, &
     LiebProb, LiebPart, Lieb_size, &
     FullPart, Full_size, &
     IWidth, PreSeed, DirName, Middlename, IErr)

  USE MyNumbers
  USE IChannels
  USE DPara
  USE IPara
  USE iso_fortran_env, ONLY: output_unit

  INTEGER(KIND=IKIND) Inum, PreSeed, ISSeed, IWidth, IErr
  INTEGER(KIND=IKIND) Lsize, Cube_size,Lieb_size,Full_size, NEVals, i,j

  REAL(KIND=RKIND) EIGS(LSize), &
       CubeProb(Cube_size), CubePart(Cube_size), &
       LiebProb(Lieb_size), LiebPart(Lieb_size), &
       FullPart(Full_size)

  CHARACTER(len=*) DirName, MiddleName
  CHARACTER*50 Prefix, Postfix
  CHARACTER*200 FileName
  
!!$  PRINT*,"DBG: WriteOutputEvecProj()"
  IErr= 0
  
  WRITE(Prefix, '(A5)') "Eval_"
  WRITE(Postfix, '(A2,I5.5,A4)') "-c",  PreSeed, ".raw" 
  !PRINT*, Prefix, Postfix
  
  FileName= TRIM(Prefix)//TRIM(MiddleName)//TRIM(Postfix)

  IF(IWriteFlag.GE.2) THEN
     PRINT*, "WriteOutputEVecProj(): ", FileName
  END IF
  
  OPEN(UNIT= IChEVec, ERR= 40, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(DirName))//"/"//FileName)

  !DO i= 1+( Lsize*( Inum -1) ), 1+( Lsize*( Inum -1) ) + Lsize
  DO Inum=1,Cube_size
     WRITE(UNIT=IChEVec, FMT=45, ERR=50) Inum, EIGS(Inum), &
          CubeProb(Inum),LiebProb(Inum), CubePart(Inum),LiebPart(Inum), FullPart(Inum)
!!$     IF( MOD(Inum,Cube_size/10)==0 ) THEN
!!$        WRITE( output_unit, '(I4)', advance = 'no') NINT(REAL(Inum)/REAL(Cube_size)*100)
!!$     END IF
  END DO
!!$  PRINT*, " "
  
  CLOSE(UNIT= IChEVec, ERR= 60)
  
  RETURN

45 FORMAT(I6,6f30.20)

  !	error in OPEN detected
40 PRINT*,"WriteOutputEVec(): ERR in OPEN()"
  IErr= 1
  RETURN

  !	error in WRITE detected
50 PRINT*,"WriteOutputEVec(): ERR in WRITE()"
  IErr= 1
  RETURN

  ! ERR in CLOSE detected
60 PRINT*,"OutputEVec(): ERR in CLOSE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteOutputEVecProj




