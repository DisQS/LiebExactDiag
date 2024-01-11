!-----------------------------------------------------------------
!
! LEDdiag
!
!-----------------------------------------------------------------
! exact diagonalization of 2D and 3D extended Lieb models
! see https://github.com/DisQS/LiebExactDiag
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!
!   Call function DSYEV() to calculate eigenvalues and eigenvectors
!   for Lieb matrix and its extendsions(2D and 3D)
!
!   IBCFlag   Not finish in this version, we calculate periodic boundary case
!
!   IRNGFlag  0 CubeConPot on Cube sites.
!             1 CubeConPot and -CubeConPot randomly on Cube sites.
!             2 CubeConPot and -CubeConPot on Cube sites with checkboard pattern   
!
!--------------------------------------------------------------------------------------

PROGRAM LiebExactDiag

!!$  use, intrinsic :: iso_c_binding
  USE MyNumbers  
  USE CConstants
  USE IConstants
  USE IPara
  USE DPara
  USE IChannels
  USE RNG_MT
  USE mt95

  IMPLICIT NONE

  !----------------------------------------
  ! Variable declaration
  !----------------------------------------

  ! paramters for Lieb matrix
  INTEGER(KIND=IKIND) IWidth, ISSeed(5)

  INTEGER(KIND=IKIND) &
       ucl, &   ! the number of atoms in a unit cell
       n_uc, &  ! the number of unit cell
       LSize, & ! the whole number of atoms in system
       NEIG     ! #eigenvalues = LSize

  REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE:: HAMMAT0

  ! Parameters for call function DSYEV()

  INTEGER(KIND=IKIND) LDA, LWMAX, INFO, LWORK

  INTEGER(KIND=IKIND), DIMENSION(:), ALLOCATABLE:: LiebSites, CubeSites
  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE:: EIGS, WORK
  REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE:: HAMMAT

  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE:: &
       CubeProb, CubePart, LiebProb, LiebPart, FullPart

  INTRINSIC        INT, MIN
  EXTERNAL         DSYEV

  ! Parameters for eigenverctor, participation numbers

  INTEGER(KIND=IKIND) Seed, i, j, Inum, IErr, col, row, len
  REAL(KIND=RKIND) drandval,SUMHUBrandval,SUMRIMrandval, CubeNorm,LiebNorm
  REAL(KIND=RKIND),ALLOCATABLE :: norm(:), part_nr(:)

  CHARACTER*100 str

  ! ----------------------------------------------------------
  ! start of main code
  ! ----------------------------------------------------------

  ! ----------------------------------------------------------
  ! protocol feature via git
  ! set: git tag -a v0.0 -m 'Version 0.0'; git push --tags
  ! ----------------------------------------------------------
#ifdef git
  PRINT*,"LiebExactDiag ", TRIM("GITVERSION"), ", ", &
       TRIM("GITBRANCH"), ", compiled: ", TRIM("COMPILED")
#else
  PRINT*,"LiebExactDiag()"
#endif

  ! ----------------------------------------------------------
  ! inout handling
  ! ----------------------------------------------------------

  CALL  Input(IErr)
  IF(IErr.NE.0) THEN
     PRINT*,"main: Input() finds IErr=", IErr
     STOP
  ENDIF

  ! ----------------------------------------------------------
  ! start of main IWidth loop
  ! ----------------------------------------------------------

  DO IWidth= Width0, Width1, dWidth

     ! ----------------------------------------------------------
     IF(IWriteFlag.GE.0) THEN
        PRINT*, "START@ IWidth=", IWidth, " Config=", ISeed
     ENDIF

     !--------------------------------------------------------------------------
     ! Setting the parameters size passed to function of generating Lieb Matrix
     !--------------------------------------------------------------------------

     ucl = (Dim * Nx) + 1 ! number of elements in unit cell
     n_uc = IWidth**Dim   ! number of unit cells
     LSize = ucl * n_uc   ! number of sites

     PRINT*,"#unit cells=", n_uc, " elements/unit cell=", ucl, " total sites=", LSize

     !--------------------------------------------------------------------------
     ! Setting the parameters size passed to function DSYEV
     !--------------------------------------------------------------------------

     LDA = LSize
     LWMAX = 100000

     ! ----------------------------------------------------------
     ! ALLOCATing memory
     ! ----------------------------------------------------------

     ALLOCATE( HAMMAT0(LSize, LSize) )
     ALLOCATE( EIGS( LSize ) )
     ALLOCATE( WORK( LWMAX ) )
     ALLOCATE( HAMMAT( LSize, LSize ) )
     ALLOCATE( CubeSites( n_uc ) )
     ALLOCATE( LiebSites( LSize - n_uc ) )
     ALLOCATE( CubeProb( LSize ) )
     ALLOCATE( LiebProb( LSize ) )
     ALLOCATE( CubePart( LSize ) )
     ALLOCATE( LiebPart( LSize ) )
     ALLOCATE( FullPart( LSize ) )
     !ALLOCATE ( norm(LSize) )
     !ALLOCATE ( part_nr(LSize) )

     HAMMAT0(:,:) = 0.0D0

     ! ----------------------------------------------------------
     IF(IWriteFlag.GE.1) THEN
        PRINT*, "---------- starting matrix build -----------"
     ENDIF

     CALL MakeLiebMatrixStructrue(Dim, Nx, IWidth, ucl, n_uc, LSize, HAMMAT0, CubeSites, LiebSites)

     !print*,"CubeSites=", CubeSites
     !print*,"LiebSites=", LiebSites

     ! ----------------------------------------------------------
     ! Determining the "flux" to vary
     ! ----------------------------------------------------------

     CubeConPot=CubeConPot0
     CubeDis=CubeDis0
     CubeDisOD=CubeDisOD0
     LiebDisOD=LiebDisOD0

     SELECT CASE(IFluxFlag)
     CASE(0) ! CubeConPot
        PRINT*,"main: (0) CubeConPot-loop"
        flux0= CubeConPot0; flux1= CubeConPot1; dflux= dCubeConPot
     CASE(1) ! CubeDis
        PRINT*,"main: (1) CubeDis-loop"
        flux0= CubeDis0; flux1= CubeDis1; dflux= dCubeDis
     CASE(2) ! CubeDis
        PRINT*,"main: (2) LiebDis-loop"
        PRINT*,"main: NOT implemented yet --- ABORTING!"; STOP
     CASE(3) ! CubeShiftOD
        PRINT*,"main: (3) CubeShiftOD-loop"
        PRINT*,"main: NOT implemented yet --- ABORTING!"; STOP
     CASE(4) ! CubeDisOD
        PRINT*,"main: (4) CubeDisOD-loop"
        flux0= CubeDisOD0; flux1= CubeDisOD1; dflux= dCubeDisOD
     CASE(5) ! LiebShiftOD
        PRINT*,"main: (5) LiebShiftOD-loop"
        PRINT*,"main: NOT implemented yet --- ABORTING!"; STOP
     CASE(6) ! LiebDisOD
        PRINT*,"main: (6) LiebDisOD-loop"
        flux0= LiebDisOD0; flux1= LiebDisOD1; dflux= dLiebDisOD
     END SELECT

     ! ----------------------------------------------------------
     ! main continuous parameter loop
     ! ----------------------------------------------------------
     DO flux = flux0, flux1, dflux

        ! ----------------------------------------------------------
        ! Setting the "flux" to vary
        ! ----------------------------------------------------------
        SELECT CASE(IFluxFlag)
        CASE(0)
           CubeConPot= flux
        CASE(1)
           CubeDis= flux
        CASE(4)
           CubeDisOD= flux
        CASE(6)
           LiebDisOD=flux
        END SELECT
        
        CALL GetDirec(Dim, Nx, IWidth, CubeDis, LiebDis, CubeConPot, LiebConPot, str)

        DO Seed=ISeed, ISeed+NSeed-1

           ! ----------------------------------------------------------
           ! Compute actual seed
           ! ----------------------------------------------------------

           !CALL SRANDOM(Seed)

           ISSeed(1)= Seed
           ISSeed(2)= IWidth
           ISSeed(3)= NINT(CubeDis*1000.) ! in lieu of "Energy"
           ISSeed(4)= NINT(CubeDis*1000.)
           ISSeed(5)= NINT(LiebDis*1000.)

           !           CALL genrand_int31(ISSeed) ! MT95 with 5 seeds

           SELECT CASE(IWriteFlag)
           CASE(1,2)
              PRINT*, "--> Seed=", Seed
              PRINT*, "--> ISSeed=", ISSeed
           CASE(3,4)
!!$              PRINT*, "IS: IW=", IWidth, "hD=", NINT(CubeDis*1000.), "E=", NINT(Energy*1000.), &
!!$                   "S=", Seed, "IS=", ISSeed
              CALL genrand_int31(ISSeed) ! MT95 with 5 seeds
              CALL genrand_real1(drandval)
              CALL SRANDOM5(ISSeed)
              drandval=DRANDOM5(ISSeed)
              WRITE(*, '(A7,I3,A4,F6.3,A4,F5.3,A3,I5,A4,F16.10)') &
                   "IS: IW=", IWidth, " CD=", CubeDis, " LD=", LiebDis, &
                   " S=", Seed, " R=", drandval
              PRINT*, "ISSeed=", ISSeed
           CASE DEFAULT
              PRINT*,"main: Seed=", Seed
           END SELECT

           ! ----------------------------------------------------------
           ! CHECK if same exists and can be overwritten
           ! ----------------------------------------------------------

           SELECT CASE(IKeepFlag)
           CASE(1)
              CALL CheckOutput( Dim,Nx, IWidth, CubeDis, LiebDis, &
                   CubeConPot, LiebConPot, Seed, str, IErr )
              IF(IErr.EQ.2) CYCLE
           END SELECT

           !CALL genrand_int31(ISSeed) ! MT95 with 5 seeds, before: CALL SRANDOM(ISSeed
           CALL SRANDOM5(ISSeed) ! MT95 with 5 seeds, before: CALL SRANDOM(ISSeed)

           ! ----------------------------------------------------------
           ! LOG which order is being used
           ! ----------------------------------------------------------

           SELECT CASE(IRNGFlag)
           CASE(0)
              PRINT*,"--- constant CubeConPot on each cube site"
           CASE(1)
              PRINT*,"--- +/- CubeConPot on RANDOM cube sites"
           CASE(2)
              PRINT*,"--- CHECKERBOARD +/- CubeConPot on each cube site"
           CASE DEFAULT
              PRINT*,"--- this IRNGFlag value is NOT implemented --- ABORTING"
           END SELECT

           ! ----------------------------------------------------------
           ! ENTER random values into matrix
           ! ----------------------------------------------------------

           HAMMAT(:,:) = HAMMAT0(:,:)
           SUMHUBrandval= 0.0D0; SUMRIMrandval= 0.0D0

           ! Give the Lieb matrix different onsite potentials
           DO i=1, n_uc

              ! CUBE onsite potentials
              drandval= DRANDOM5(ISSeed)

              SUMHUBrandval=SUMHUBrandval + CubeDis*(drandval - 0.5D0)

              SELECT CASE (IRNGFlag)
              CASE(0) ! constant CubeConPot on each cube site
                 HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = CubeConPot + CubeDis*(drandval - 0.5D0)
              CASE(1) ! +/- CubeConPot on random cube sites

                 IF(MOD(n_uc,2) .NE. 0) THEN
                    PRINT*, "main: WRNG, cube size are odd, so we can not achieve 0 potential!"
                    PRINT*, "main: WRNG, calculation will proceed, but output is questionable."
                 END IF

                 IF(DRANDOM5(ISSeed)>=0.5) THEN
                    HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = &
                         CubeConPot + CubeDis*(drandval - 0.5D0)
                 ELSE
                    HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = &
                         -CubeConPot + CubeDis*(drandval - 0.5D0)
                 END IF

              CASE(2) ! checkerboard +/- CubeConPot on each cube site

                 IF(MOD(n_uc,2) .NE. 0) THEN
                    PRINT*, "main: WRNG, cube size are odd, so we can not achieve 0 potential!"
                    PRINT*, "main: WRNG, calculation will proceed, but output is questionable."
                 END IF

                 IF(Dim==2)THEN

                    len=1
                    IF(MOD(i,IWidth)==0)THEN
                       col=IWidth
                       row=INT(i/IWidth)
                    ELSE
                       col=MOD(i,IWidth)
                       row=INT(i/IWidth)+1
                    END IF

                 ELSE IF(Dim==3)THEN

                    IF(MOD(i,IWidth**2)==0)THEN
                       col=IWidth
                       row=IWidth
                       len=INT(i/IWidth**2)
                    ELSE
                       IF(MOD(MOD(i,IWidth**2),IWidth)==0)THEN
                          col=IWidth
                          row=INT(MOD(i,IWidth**2)/IWidth)
                       ELSE
                          col=MOD(MOD(i,IWidth**2),IWidth)
                          row=INT(MOD(i,IWidth**2)/IWidth)+1
                       END IF
                       len=INT(i/IWidth**2)+1
                    END IF

                 ELSE
                    PRINT*,"Not finished yet!"
                    STOP
                 END IF

                 drandval= DRANDOM5(ISSeed)
                 SUMHUBrandval=SUMHUBrandval + CubeDis*(drandval - 0.5D0)

                 IF(MOD(len,2)==1)THEN

                    IF((MOD(col,2)==1 .AND. MOD(row,2)==1) .OR. (MOD(col,2)==0 .AND. MOD(row,2)==0))THEN
                       HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = CubeConPot + CubeDis*(drandval - 0.5D0)
                    ELSE
                       HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = -1.0*CubeConPot + CubeDis*(drandval - 0.5D0)
                    END IF

                 ELSE
                    IF((MOD(col,2)==1 .AND. MOD(row,2)==1) .OR. (MOD(col,2)==0 .AND. MOD(row,2)==0))THEN
                       HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = -1.0*CubeConPot + CubeDis*(drandval - 0.5D0)
                    ELSE
                       HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = CubeConPot + CubeDis*(drandval - 0.5D0)
                    END IF

                 END IF

              END SELECT

              ! LIEB onsite potentials
              DO j=1, ucl-1

                 drandval= DRANDOM5(ISSeed)
                 SUMRIMrandval=SUMRIMrandval + LiebDis*(drandval - 0.5D0)
                 HAMMAT((i-1)*ucl + j + 1 , (i-1)*ucl + j + 1) = LiebConPot + LiebDis*(drandval - 0.5D0)

              END DO

           END DO

!!$           Do i=1,Lsize
!!$              Print*,i, HAMMAT(i,i)
!!$           END Do
!!$           Pause

           ! ----------------------------------------------------------
           ! WRITE SUMrandval to allow identification of accidental states
           ! ----------------------------------------------------------

           PRINT*,"main(): Seed=", Seed, ", SHrv=", SUMHUBrandval/n_uc, &
                ", SRrv=", SUMRIMrandval/n_uc

           ! ----------------------------------------------------------
           ! START the diagonalizstion process
           ! ----------------------------------------------------------

           PRINT*, "STARTing the diagonalization process"

           !PRINT*, "DSYEV: computing WORK space"
           LWORK =  -1  !3*LSize
           CALL DSYEV( 'V', 'Upper', LSize, HAMMAT, LSize, EIGS, WORK, LWORK, INFO )

           !PRINT*, "DSYEV: using WORK space"
           LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
           CALL DSYEV( 'V', 'Upper', LSize, HAMMAT, LSize, EIGS, WORK, LWORK, INFO ) 

           IF(INFO.NE.0) THEN
              PRINT*,"main: DYSEV(INFO)=", INFO
              STOP
           ENDIF

           ! ----------------------------------------------------------
           ! WRITE the eigenvalues and -vectors
           ! ----------------------------------------------------------

           !CALL WriteEvals(Dim, Nx, IWidth, LSize, CubeDis, LiebDis, EIGS, HAMMAT, norm, part_nr, Seed, INFO)

           NEIG=LSize ! this is complete diagonalization

           CALL WriteOutputEVal( Dim, Nx, NEIG, EIGS, IWidth, &
                CubeDis, LiebDis, CubeConPot, LiebConPot, Seed, str, IErr)

           SELECT CASE(IStateFlag)
           CASE(0)
              CONTINUE
           CASE(1)
              PRINT*,"main: DYSEV() eigenvectors will now be saved into individual files"
              DO Inum= 1,NEIG
                 CALL WriteOutputEVec(Dim, Nx, Inum, NEIG, Lsize, HAMMAT, LSize, &
                      IWidth, CubeDis, LiebDis, CubeConPot, LiebConPot,&
                      Seed, str, IErr)
              END DO
           CASE(2)
              PRINT*,"main: DYSEV() eigenvectors will now be saved into single BULK file"

              CALL WriteOutputEVecBULK(Dim, Nx, Lsize, NEIG, Lsize, EIGS, LSize, &
                   IWidth, CubeDis, LiebDis, CubeConPot, LiebConPot, &
                   Seed, str, IErr)
           CASE(-1)
              PRINT*,"main: Cube/Lieb site projections of DYSEV() eigenvectors"

              ! compute projected probabilities 
              DO Inum= 1,NEIG

                 CubeProb(Inum)= 0.0
                 DO i=1,n_uc
                    CubeProb(Inum)= CubeProb(Inum) + &
                         HAMMAT(CubeSites(i),Inum) * HAMMAT(CubeSites(i),Inum)
                 END DO

                 LiebProb(Inum)= 0.0
                 DO i=1,LSize-n_uc
                    LiebProb(Inum)= LiebProb(Inum) + &
                         HAMMAT(LiebSites(i),Inum) * HAMMAT(LiebSites(i),Inum)
                 END DO

              ENDDO

              ! compute FULL participation numbers (already normalized) 
              DO Inum= 1,NEIG

                 FullPart(Inum)= 0.0
                 DO i=1,LSize
                    FullPart(Inum)= FullPart(Inum) + &
                         HAMMAT(i,Inum) * HAMMAT(i,Inum) * &
                         HAMMAT(i,Inum) * HAMMAT(i,Inum) 
                 END DO
                 FullPart(Inum)=1.0D0/FullPart(Inum) /LSize
                 !PRINT*,Inum, FullPart(Inum)

              END DO

              ! compute the projected normalizations
              DO Inum= 1,NEIG

                 CubeNorm= 0.0
                 DO i=1,n_uc
                    CubeNorm= CubeNorm + &
                         HAMMAT(CubeSites(i),Inum) * HAMMAT(CubeSites(i),Inum) 
                 END DO

                 LiebNorm= 0.0
                 DO i=1,LSize-n_uc
                    LiebNorm= LiebNorm + &
                         HAMMAT(LiebSites(i),Inum) * HAMMAT(LiebSites(i),Inum)
                 END DO

                 !renormalize
                 DO i=1,n_uc
                    HAMMAT(CubeSites(i),Inum) =  HAMMAT(CubeSites(i),Inum) /SQRT(CubeNorm)
                 END DO

                 DO i=1,LSize-n_uc
                    HAMMAT(LiebSites(i),Inum) = HAMMAT(LiebSites(i),Inum) /SQRT(LiebNorm)
                 END DO

              ENDDO

              ! compute normalized participation numbers
              DO Inum= 1,NEIG

                 CubePart(Inum)= 0.0
                 DO i=1,n_uc
                    CubePart(Inum)= CubePart(Inum) + &
                         HAMMAT(CubeSites(i),Inum) * HAMMAT(CubeSites(i),Inum) * &
                         HAMMAT(CubeSites(i),Inum) * HAMMAT(CubeSites(i),Inum) 
                 END DO
                 CubePart(Inum)=1/CubePart(Inum) /n_uc

                 LiebPart(Inum)= 0.0
                 DO i=1,LSize-n_uc
                    LiebPart(Inum)= LiebPart(Inum) + &
                         HAMMAT(LiebSites(i),Inum) * HAMMAT(LiebSites(i),Inum) * &
                         HAMMAT(LiebSites(i),Inum) * HAMMAT(LiebSites(i),Inum)
                 END DO
                 LiebPart(Inum)=1/LiebPart(Inum) /(LSize-n_uc)

              ENDDO

              CALL WriteOutputEVecProj(Dim, Nx, Inum, NEIG, &
                   EIGS, LSize, &
                   CubeProb, CubePart, Lsize, &
                   LiebProb, LiebPart, LSize, &
                   FullPart, LSize, &
                   IWidth, CubeDis, LiebDis, &
                   CubeConPot, LiebConPot, &
                   Seed, str, IErr)

           END SELECT

        END DO ! ISeed cycle

     END DO  ! flux cycle

     DEALLOCATE ( HAMMAT0, EIGS, WORK, HAMMAT, CubeSites, LiebSites )

  END DO ! IWidth cycle

  STOP 'LiebExactDiag'

END PROGRAM





