! ********************************************************************
!       
! For Lieb model, call the function DSYEV() to calculate eigenvalues
! and eigenvectors
!
! ********************************************************************
      

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
  
  !OPEN(UNIT= IChInp, ERR=120, FILE= "LEDdiag.inp",STATUS= 'OLD')
OPEN(UNIT= IChInp, ERR=120, FILE= "/dev/stdin",STATUS= 'OLD')

  ILine= ILine+1
  READ(IChInp,10,ERR=20) ISeed
  !PRINT*,"ISeed        = ",ISeed

  ILine= ILine+1
  READ(IChInp,10,ERR=20) NSeed
  !PRINT*,"NSeed        = ",NSeed

  ILine= ILine+1
  READ(IChInp,10,ERR=20) Dim
  !PRINT*,"Dim         = ",Dim

  ILine= ILine+1
  READ(IChInp,10,ERR=20) Nx
  !PRINT*,"Nx          = ",Nx
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IBCFlag
  !PRINT*,"IBCFlag      = ", IBCFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IRNGFlag
  !PRINT*,"IRNGFlag     = ", IRNGFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IKeepFlag
  !PRINT*,"IKeepFlag    = ", IKeepFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IWriteFlag
  !PRINT*,"IWriteFlag   = ", IWriteFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IStateFlag
  !PRINT*,"IWriteFlag   = ", IStateFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) Width0
  !PRINT*,"Width0       = ",Width0
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) Width1
  !PRINT*,"Width1       = ", Width1
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) dWidth
  !PRINT*,"dWidth       = ", dWidth
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20) HubDis0
  !PRINT*,"HubDis0     = ", HubDis0
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20) HubDis1
  !PRINT*,"HubDis1     = ", HubDis1
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20) dHubDis
  !PRINT*,"dHubDis     = ", dHubDis

  ILine= ILine+1
  READ(IChInp,15,ERR=20) RimDis
  !PRINT*,"RimDis      = ", RimDis

  
10 FORMAT(16X,I15.1)
  ! 10	FORMAT("IMAXIteration= ",I15.1)
15 FORMAT(16X,F18.9)
  ! 15     FORMAT("IMAXIteration= ",F18.9)

  ! check the parameters for validity
  
  IF(IWriteFlag.GE.2) THEN
     PRINT*,"ISeed        = ", ISeed
     PRINT*,"NSeed        = ", NSeed
     PRINT*,"Dim          = ", Dim
     PRINT*,"Nx           = ", Nx
     PRINT*,"IBCFlag      = ", IBCFlag
     PRINT*,"IRNGFlag     = ", IRNGFlag
     PRINT*,"IKeepFlag    = ", IKeepFlag
     PRINT*,"IWriteFlag   = ", IWriteFlag
     PRINT*,"IStateFlag   = ", IStateFlag
     PRINT*,"Width0       = ", Width0
     PRINT*,"Width1       = ", Width1
     PRINT*,"dWidth       = ", dWidth
     PRINT*,"HubDis0     = ", HubDis0
     PRINT*,"HubDis1     = ", HubDis1
     PRINT*,"dHubDis     = ", dHubDis
     PRINT*,"RimDis   = ", RimDis
    
  ENDIF

  CLOSE(IChInp)
  RETURN

120 PRINT*,"Input(): ERR in OPEN()"  
  IErr= 1
  RETURN

20 PRINT*,"Input(): ERR in WRITE()"  
  IErr= 1
  RETURN
  
END SUBROUTINE Input

! --------------------------------------------------------------------
! GetFileName
!
! IErr	error code
!----------------------------------------------------------------------

FUNCTION GetFileName(fnamestr,vdata,n)
!write(*,*) trim(GetFileName('/Output/QH_L(I4)_NL(I1)_NS(I1)_B(F5.3)_C(F4.2)_S(I5).dat',&
!     (/L_INPUT,NLevels,NSpins,BField,Coul,Seed/),6))

  IMPLICIT NONE

  INTEGER,INTENT(in)         :: n
  REAL(8),INTENT(in)         :: vdata(n)
  CHARACTER(len=*),INTENT(in) :: fnamestr

  CHARACTER(len=*)   :: GetFileName
  CHARACTER(len=500) :: str,fstr,vstr
  
  INTEGER           :: i

  str = fnamestr
  
  DO i=1,n
     
     fstr = str(INDEX(str,'('):INDEX(str,')'))
     IF(fstr(2:2).EQ.'I')THEN
        WRITE(vstr,fstr) INT(vdata(i))
     ELSE
        WRITE(vstr,fstr) vdata(i)
     END IF
  
     str = str(1:INDEX(str,'(')-1) // TRIM(ADJUSTL(vstr)) // str(INDEX(str,')')+1:LEN(str)) 
     
  END DO
  
  GetFileName = TRIM(ADJUSTL(str))

  RETURN
END FUNCTION GetFileName

!--------------------------------------------------------------------
! CheckOutput:
!
! IErr	error code

SUBROUTINE CheckOutput( Dim, Nx, IWidth, Energy, HubDis, RimDis, PreSeed, str, IErr )

  USE MyNumbers 
  USE IChannels
!  USE DPara
!  USE IPara

  
  INTEGER(KIND=IKIND) Dim, Nx, IWidth, IErr, PreSeed, ISSeed
  REAL(KIND=RKIND) HubDis, RimDis, Energy
  
  CHARACTER*100 FileName, str
  
  PRINT*,"DBG: CheckOutput()"
  
  IErr= 0
  
  !   WRITE out the input parameter
  IF(Energy.GE.0.0D0) THEN
     WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A5,A3,I6.6,A3,I6.6,A2,I5.5,A4)') &
          "Eval-", "L", Dim, Nx, &
          "-M",IWidth, &
          "-Espec", &! NINT(100.*ABS(Energy)), &
          "-hD", NINT(100.*ABS(HubDis)), &
          "-rD", NINT(100.*ABS(RimDis)), "-c",& 
          PreSeed, ".raw" !"_s", ISSeed, 
  ELSE
     WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A6,A3,I6.6,A3,I6.6,A2,I5.5,A4)') &
          "Eval-","L",Dim, Nx, &
          "-M",IWidth, &
          "-Espec", & !NINT(100.*ABS(Energy)), &
          "-hD",NINT(100.*ABS(HubDis)), &
          "-rD", NINT(100.*ABS(RimDis)), "-c",&
          PreSeed, ".raw" !"_s", ISSeed, 
  ENDIF
  
  OPEN(UNIT= IChOut, ERR= 10, STATUS= 'NEW', FILE= TRIM(ADJUSTL(str))//"/"//FileName)
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

SUBROUTINE WriteOutputEVal(Dim, Nx, NEVals, EIGS, IWidth, Energy, HubDis, RimDis, PreSeed, str, IErr)

  USE MyNumbers
  USE IChannels
!  USE DPara
!  USE IPara

  INTEGER(KIND=IKIND) Dim, Nx
  INTEGER(KIND=IKIND) IWidth, IErr, NEVals, PreSeed, ISSeed, i
  REAL(KIND=RKIND) HubDis, RimDis, Energy
  REAL(KIND=RKIND) EIGS(NEVals)
  
  CHARACTER*100 FileName, str

  PRINT*,"DBG: WriteOutputEVal()"
  
  IErr= 0
  
  !   WRITE out the input parameter
  IF(Energy.GE.0.0D0) THEN
     WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A6,A3,I6.6,A3,I6.6,A2,I5.5,A4)') &
          "Eval-", "L", Dim, Nx, &
          "-M",IWidth, &
          "-Espec", & !NINT(100.*ABS(Energy)), &
          "-hD", NINT(100.*ABS(HubDis)), &
          "-rD", NINT(100.*ABS(RimDis)), "-c",& 
          PreSeed, ".raw" !"_s", ISSeed, 
  ELSE
     WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A6,A3,I6.6,A3,I6.6,A2,I5.5,A4)') &
          "Eval-","L",Dim, Nx, &
          "-M",IWidth, &
          "-Espec", & !NINT(100.*ABS(Energy)), &
          "-hD",NINT(100.*ABS(HubDis)), &
          "-rD", NINT(100.*ABS(RimDis)), "-c",&
          PreSeed, ".raw" !"_s", ISSeed, 
  ENDIF
  
!  IF(IWriteFlag.GE.2) THEN
     PRINT*, "WriteOutputEVal(): ", FileName
!  ENDIF
  
!!$  OPEN(UNIT= IChEVal, ERR= 10, STATUS='UNKNOWN', FILE=Trim(str)//"/"//FileName)  
  
  OPEN(UNIT= IChEVal, ERR= 10, STATUS='UNKNOWN', FILE= TRIM(ADJUSTL(str))//"/"//FileName)
  
  IF(NEVals .GT. 0)THEN
     DO i=1,NEVals
        WRITE(IChEVal, FMT=15, ERR=20) EIGS(i)
     ENDDO
  END IF
  
  CLOSE(UNIT=IChEVal, ERR= 30)
  
  RETURN
  
15 FORMAT(f30.20)

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

SUBROUTINE WriteOutputEVec( Dim, Nx, Inum, NEVals, Lsize, VECS, VECS_size, &
                   IWidth, Energy, HubDis, RimDis, PreSeed, str, IErr)

  USE MyNumbers
  USE IChannels
!  USE DPara
!  USE IPara

  INTEGER(KIND=IKIND) Dim, Nx
  INTEGER(KIND=IKIND) Inum, PreSeed, ISSeed, IWidth, IErr, Lsize, VECS_size, NEVals, i
  REAL(KIND=RKIND) HubDis, RimDis, Energy

  REAL(KIND=RKIND) VECS(VECS_size,VECS_size)

  CHARACTER*100 FileName, str

  PRINT*,"DBG: WriteOutputEvec()"

  IErr= 0

  !   WRITE out the input parameter
  IF(Energy.GE.0.0D0) THEN
     WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A6,A3,I6.6,A3,I6.6,A2,I5.5,A2,I4.4,A4)') &
     !WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A2,A5,I9.9,A7,I7.7,A7,I7.7,A2,I4.4,A1,I5.5,A4)') &
          "Evec-","L", Dim, Nx, &
          "-M", IWidth, &
          "-Espec", & !NINT(100.0D0*ABS(Energy)), &
          "-hD", NINT(100.0D0*ABS(HubDis)), &
          "-rD", NINT(100.0D0*ABS(RimDis)), &
          "-c", PreSeed, "-N", Inum, & !"_s", ISSeed, 
          ".raw"
  ELSE
     WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A6,A3,I6.6,A3,I6.6,A2,I5.5,A2,I4.4,A4)') &
     !WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A2,A5,I9.9,A7,I7.7,A7,I7.7,A2,I4.4,A1,I5.5,A4)') &
          "Evec-","L",Dim, Nx, &
          "-M", IWidth, &
          "-Espec", & !NINT(100.0D0*ABS(Energy)), &
          "-dD", NINT(100.0D0*ABS(HubDis)), &
          "-rD", NINT(100.0D0*ABS(RimDis)), &
          "-c", PreSeed, "-N", Inum, & !"_s", ISSeed, 
          ".raw"
  ENDIF

  PRINT*, "WriteOutputEVec(): ", FileName

  OPEN(UNIT= IChEVec, ERR= 40, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(str))//"/"//FileName)

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

SUBROUTINE WriteOutputEVecBULK( Dim, Nx, Inum, NEVals, Lsize, VECS, VECS_size, &
                   IWidth, Energy, HubDis, RimDis, PreSeed, str, IErr)

  USE MyNumbers
  USE IChannels
!  USE DPara
!  USE IPara

  INTEGER(KIND=IKIND) Dim, Nx
  INTEGER(KIND=IKIND) Inum, PreSeed, ISSeed, IWidth, IErr, Lsize, VECS_size, NEVals, i,j
  REAL(KIND=RKIND) HubDis, RimDis, Energy

  REAL(KIND=RKIND) VECS(VECS_size,VECS_size)

  CHARACTER*100 FileName, str

  PRINT*,"DBG: WriteOutputEvecBULK()"

  IErr= 0

  !   WRITE out the input parameter
  IF(Energy.GE.0.0D0) THEN
     WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A6,A3,I6.6,A3,I6.6,A2,I5.5,A2,I4.4,A4)') &
     !WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A2,A5,I9.9,A7,I7.7,A7,I7.7,A2,I4.4,A1,I5.5,A4)') &
          "Evec-","L", Dim, Nx, &
          "-M", IWidth, &
          "-EFULL", & !NINT(100.0D0*ABS(Energy)), &
          "-hD", NINT(100.0D0*ABS(HubDis)), &
          "-rD", NINT(100.0D0*ABS(RimDis)), &
          "-c", PreSeed, "-N", Inum, & !"_s", ISSeed, 
          ".raw"
  ELSE
     WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A6,A3,I6.6,A3,I6.6,A2,I5.5,A2,I4.4,A4)') &
     !WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A2,A5,I9.9,A7,I7.7,A7,I7.7,A2,I4.4,A1,I5.5,A4)') &
          "Evec-","L",Dim, Nx, &
          "-M", IWidth, &
          "-EFULL", & !NINT(100.0D0*ABS(Energy)), &
          "-dD", NINT(100.0D0*ABS(HubDis)), &
          "-rD", NINT(100.0D0*ABS(RimDis)), &
          "-c", PreSeed, "-N", Inum, & !"_s", ISSeed, 
          ".raw"
  ENDIF

  PRINT*, "WriteOutputEVecBULK(): ", FileName

  OPEN(UNIT= IChEVec, ERR= 40, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(str))//"/"//FileName)

  !DO i= 1+( Lsize*( Inum -1) ), 1+( Lsize*( Inum -1) ) + Lsize
  DO j=1,Inum
     DO i=1,LSize
        WRITE(UNIT=IChEVec, FMT=45, ERR=50) j, i, VECS(Inum,i)
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Create Folder !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GetDirec(Dim, Nx, Width, HubDis, RimDis, Energy, str)
  USE MyNumbers

  INTEGER*4 Dim, Nx, Width, Seed
  REAL*8 HubDis, RimDis, Energy
  CHARACTER(len=100) str
  CHARACTER(len=10) fid1, fid2, fid3, fid4, fid5, fid6
  LOGICAL*4 ierr1

  PRINT*, "GetDirec(): ", Dim, Nx, Width, HubDis, RimDis, Energy

  WRITE(fid1,'(I1)') Dim; fid1=TRIM(ADJUSTL(fid1))
  WRITE(fid2,'(I1)') Nx; fid2=TRIM(ADJUSTL(fid2))
  WRITE(fid3,'(I3)') Width; fid3=TRIM(ADJUSTL(fid3))
  WRITE(fid4,'(I6.6)') NINT(HubDis*100.); fid4=TRIM(ADJUSTL(fid4))
  WRITE(fid5,'(I6.6)') NINT(RimDis*100.); fid5=TRIM(ADJUSTL(fid5))
  !WRITE(fid6,'(I4.4)') Seed
  !WRITE(fid6,'(I6.6)') NINT(ABS(Energy)*100.); fid6=TRIM(ADJUSTL(fid6))
  
  IF(Energy.GE.ZERO) THEN
     str='L'//TRIM(fid1)//TRIM(fid2)//'_M'//TRIM(fid3)//'_hD'//TRIM(fid4) &
          //'_rD'//TRIM(fid5)!//"_E"//TRIM(fid6)!//'_DATA'
  ELSE
     str='L'//TRIM(fid1)//TRIM(fid2)//'_M'//TRIM(fid3)//'_hD'//TRIM(fid4) &
          //'_rD'//TRIM(fid5)!//"_E-"//TRIM(fid6)!//'_DATA'
  END IF

!  Write(str,'(A1,I1,I1,A2,I3.1,A7,f6.1,A7,f6.1,A6)') &
!       "L", Dim, Nx, "_M", Width, "_HubDis", HubDis, &
!       "_RimDis", RimDis, "_.DATA"

!!$  PRINT*,str
  PRINT*, "GetDirec(): checking for ", str, Dim, Nx, Width, HubDis, RimDis, Energy

#ifdef ifort
  INQUIRE(directory=TRIM(ADJUSTL(str)), Exist=ierr1) ! ifort
#else
  INQUIRE(file=TRIM(ADJUSTL(str)), Exist=ierr1) ! gfortran
#endif
  IF(ierr1)THEN
     PRINT*,"GetDirec(): directory exists and doesn't need to be created!"
     WRITE(*,'(/)')
  ELSE
     PRINT*,"GetDirec(): directory doesn't exist and is NOW being createed!"
     WRITE(*,'(/)')
     CALL System("mkdir -p "//TRIM(ADJUSTL(str)) )
  END IF
  
  RETURN
  
END SUBROUTINE GetDirec



