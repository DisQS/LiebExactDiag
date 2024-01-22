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

  REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE:: HAMMAT

  ! Parameters for call function DSYEV()

  INTEGER(KIND=IKIND) LDA, LWMAX, INFO, LWORK

  INTEGER(KIND=IKIND), DIMENSION(:), ALLOCATABLE:: LiebSites, CubeSites
  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE:: EIGS, WORK
  !REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE:: HAMMAT

  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE:: &
       CubeProb, CubePart, LiebProb, LiebPart, FullPart

  INTRINSIC        INT, MIN
  EXTERNAL         DSYEV

  ! Parameters for eigenverctor, participation numbers

  INTEGER(KIND=IKIND) Seed, i, j, Inum, IErr, col, row, len
  REAL(KIND=RKIND) drandval,SUMHUBrandval,SUMRIMrandval, CubeNorm,LiebNorm
  REAL(KIND=RKIND),ALLOCATABLE :: norm(:), part_nr(:)

  CHARACTER*200 middlename, MakeMiddleName, dirname

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

  CALL Input(IErr)
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

     IF(IWriteFlag.GE.1) THEN
        PRINT*,"main: total sites=", LSize
        PRINT*,"main: #unit cells=", n_uc, " elements/unit cell=", ucl
     ENDIF
     
     !--------------------------------------------------------------------------
     ! Setting the parameters passed to function DSYEV
     !--------------------------------------------------------------------------

     LDA = LSize
     LWMAX = 100000

     ! ----------------------------------------------------------
     ! ALLOCATing memory and initializing
     ! ----------------------------------------------------------

     ALLOCATE( HAMMAT(LSize, LSize) )
     ALLOCATE( EIGS( LSize ) )
     ALLOCATE( WORK( LWMAX ) )
     !ALLOCATE( HAMMAT( LSize, LSize ) )
     ALLOCATE( CubeSites( n_uc ) )
     ALLOCATE( LiebSites( LSize - n_uc ) )
     ALLOCATE( CubeProb( LSize ) )
     ALLOCATE( LiebProb( LSize ) )
     ALLOCATE( CubePart( LSize ) )
     ALLOCATE( LiebPart( LSize ) )
     ALLOCATE( FullPart( LSize ) )
     !ALLOCATE ( norm(LSize) )
     !ALLOCATE ( part_nr(LSize) )

     HAMMAT(:,:) = 0.0D0

     ! ----------------------------------------------------------
     IF(IWriteFlag.GE.1) THEN
        PRINT*, "---------- starting matrix build -----------"
     ENDIF

     CALL MakeLiebMatrixStructrue(Dim, Nx, IWidth, ucl, n_uc, LSize, HAMMAT, CubeSites, LiebSites)

     !print*,"CubeSites=", CubeSites
     !print*,"LiebSites=", LiebSites

     ! ----------------------------------------------------------
     ! Determining the "flux" to vary
     ! ----------------------------------------------------------

     CubeConPot=CubeConPot0
     CubeDis=CubeDis0
     OffDShift=OffDShift0
     OffDDis=OffDDis0

     SELECT CASE(IFluxFlag)
     CASE(0) ! CubeConPot
        PRINT*,"main: (0) CubeConPot-loop"
        flux0= CubeConPot0; flux1= CubeConPot1; dflux= dCubeConPot
     CASE(1) ! CubeDis
        PRINT*,"main: (1) CubeDis-loop"
        flux0= CubeDis0; flux1= CubeDis1; dflux= dCubeDis
     CASE(2) ! LiebConPot
        PRINT*,"main: (2) LiebConPot-loop"
        PRINT*,"main: NOT implemented yet --- ABORTING!"; STOP
     CASE(3) ! LiebDis
        PRINT*,"main: (3) LiebDis-loop"
        PRINT*,"main: NOT implemented yet --- ABORTING!"; STOP
     CASE(4) ! OffDShift
        PRINT*,"main: (4) OffDShift-loop"
        flux0= OffDShift0; flux1= OffDShift1; dflux= dOffDShift
     CASE(5) ! OffDDis
        PRINT*,"main: (5) OffDDis-loop"
        flux0= OffDDis0; flux1= OffDDis1; dflux= dOffDDis
     CASE DEFAULT
        PRINT*,"main: NOT implemented --- ABORTING!"; STOP
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
           OffDShift= flux
        CASE(5)
           OffDDis=flux
        END SELECT

        ! ----------------------------------------------------------
        ! generate filename templates and data directory
        ! ----------------------------------------------------------

        middlename= TRIM(MakeMiddleName( IWidth, IErr ) )
        dirname= middlename
        !PRINT*, middlename
        
        CALL MakeDataDir(IWidth, dirname, IErr)

        ! ----------------------------------------------------------
        ! main loop over samples
        ! ----------------------------------------------------------
        
        DO Seed=ISeed, ISeed+NSeed-1

           ! ----------------------------------------------------------
           ! Compute seed specific for chosen parameters
           ! ----------------------------------------------------------

           !CALL SRANDOM(Seed)

           ISSeed(1)= Seed
           ISSeed(2)= IWidth
           ISSeed(3)= NINT(CubeConPot*1000.) ! in lieu of "Energy"
           ISSeed(4)= NINT(CubeDis*1000.) + NINT(OffDShift*10000.)
           ISSeed(5)= NINT(LiebDis*1000.) + NINT(OffDDis*10000.)

           SELECT CASE(IWriteFlag)
           CASE(1,2 )! showing seeds constructed
              PRINT*, "--> Seed=", Seed
              PRINT*, "--> ISSeed=", ISSeed
           CASE(3,4) ! for testing seed construction
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
              CALL CheckOutput( IWidth,Seed, dirname,middlename, IErr )
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
           ! ENTER random values into matrix DIAGONALS

           CALL MakeLiebOnsiteDisorder(Dim, Nx, IWidth, ucl, n_uc, LSize, HAMMAT, LSize, ISSeed)

           ! ----------------------------------------------------------
           ! ENTER random values into matrix OFF-diagonals

           CALL MakeLiebHoppingDisorder(Dim, Nx, IWidth, ucl, n_uc, LSize, HAMMAT, CubeSites, LiebSites)
           
!!$           HAMMAT(:,:) = HAMMAT0(:,:)
!!$           SUMHUBrandval= 0.0D0; SUMRIMrandval= 0.0D0
!!$
!!$           ! Give the Lieb matrix different onsite potentials
!!$           DO i=1, n_uc
!!$
!!$              ! CUBE onsite potentials
!!$              drandval= DRANDOM5(ISSeed)
!!$
!!$              SUMHUBrandval=SUMHUBrandval + CubeDis*(drandval - 0.5D0)
!!$
!!$              SELECT CASE (IRNGFlag)
!!$              CASE(0) ! constant CubeConPot on each cube site
!!$                 HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = CubeConPot + CubeDis*(drandval - 0.5D0)
!!$              CASE(1) ! +/- CubeConPot on random cube sites
!!$
!!$                 IF(MOD(n_uc,2) .NE. 0) THEN
!!$                    PRINT*, "main: WRNG, cube size are odd, so we cannot achieve 0 potential!"
!!$                    PRINT*, "main: WRNG, calculation will proceed, but output is questionable."
!!$                 END IF
!!$
!!$                 IF(DRANDOM5(ISSeed)>=0.5) THEN
!!$                    HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = &
!!$                         CubeConPot + CubeDis*(drandval - 0.5D0)
!!$                 ELSE
!!$                    HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = &
!!$                         -CubeConPot + CubeDis*(drandval - 0.5D0)
!!$                 END IF
!!$
!!$              CASE(2) ! checkerboard +/- CubeConPot on each cube site
!!$
!!$                 IF(MOD(n_uc,2) .NE. 0) THEN
!!$                    PRINT*, "main: WRNG, cube size are odd, so we can not achieve 0 potential!"
!!$                    PRINT*, "main: WRNG, calculation will proceed, but output is questionable."
!!$                 END IF
!!$
!!$                 IF(Dim==2)THEN
!!$
!!$                    len=1
!!$                    IF(MOD(i,IWidth)==0)THEN
!!$                       col=IWidth
!!$                       row=INT(i/IWidth)
!!$                    ELSE
!!$                       col=MOD(i,IWidth)
!!$                       row=INT(i/IWidth)+1
!!$                    END IF
!!$
!!$                 ELSE IF(Dim==3)THEN
!!$
!!$                    IF(MOD(i,IWidth**2)==0)THEN
!!$                       col=IWidth
!!$                       row=IWidth
!!$                       len=INT(i/IWidth**2)
!!$                    ELSE
!!$                       IF(MOD(MOD(i,IWidth**2),IWidth)==0)THEN
!!$                          col=IWidth
!!$                          row=INT(MOD(i,IWidth**2)/IWidth)
!!$                       ELSE
!!$                          col=MOD(MOD(i,IWidth**2),IWidth)
!!$                          row=INT(MOD(i,IWidth**2)/IWidth)+1
!!$                       END IF
!!$                       len=INT(i/IWidth**2)+1
!!$                    END IF
!!$
!!$                 ELSE
!!$                    PRINT*,"Not finished yet!"
!!$                    STOP
!!$                 END IF
!!$
!!$                 drandval= DRANDOM5(ISSeed)
!!$                 SUMHUBrandval=SUMHUBrandval + CubeDis*(drandval - 0.5D0)
!!$
!!$                 IF(MOD(len,2)==1)THEN
!!$
!!$                    IF((MOD(col,2)==1 .AND. MOD(row,2)==1) .OR. (MOD(col,2)==0 .AND. MOD(row,2)==0))THEN
!!$                       HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = CubeConPot + CubeDis*(drandval - 0.5D0)
!!$                    ELSE
!!$                       HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = -1.0*CubeConPot + CubeDis*(drandval - 0.5D0)
!!$                    END IF
!!$
!!$                 ELSE
!!$                    IF((MOD(col,2)==1 .AND. MOD(row,2)==1) .OR. (MOD(col,2)==0 .AND. MOD(row,2)==0))THEN
!!$                       HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = -1.0*CubeConPot + CubeDis*(drandval - 0.5D0)
!!$                    ELSE
!!$                       HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = CubeConPot + CubeDis*(drandval - 0.5D0)
!!$                    END IF
!!$
!!$                 END IF
!!$
!!$              END SELECT
!!$
!!$              ! LIEB onsite potentials
!!$              DO j=1, ucl-1
!!$
!!$                 drandval= DRANDOM5(ISSeed)
!!$                 SUMRIMrandval=SUMRIMrandval + LiebDis*(drandval - 0.5D0)
!!$                 HAMMAT((i-1)*ucl + j + 1 , (i-1)*ucl + j + 1) = LiebConPot + LiebDis*(drandval - 0.5D0)
!!$
!!$              END DO
!!$
!!$           END DO

           
           DO i=1,Lsize
              PRINT*,i, HAMMAT(i,i)
           END DO
           !PAUSE

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

           CALL WriteOutputEVal(NEIG, EIGS, IWidth, Seed, dirname, middlename, IErr)

           SELECT CASE(IStateFlag)
           CASE(0)
              CONTINUE
           CASE(1)
              PRINT*,"main: DYSEV() eigenvectors will now be saved into individual files"
              DO Inum= 1,NEIG
                 CALL WriteOutputEVec(Inum, NEIG, Lsize, HAMMAT, LSize, &
                      IWidth, Seed, dirname, middlename, IErr)
              END DO
           CASE(2)
              PRINT*,"main: DYSEV() eigenvectors will now be saved into single BULK file"

              CALL WriteOutputEVecBULK(NEIG, Lsize, EIGS, LSize, &
                   IWidth, Seed, dirname, middlename, IErr)
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

              CALL WriteOutputEVecProj(Inum, NEIG, &
                   EIGS, LSize, &
                   CubeProb, CubePart, Lsize, &
                   LiebProb, LiebPart, LSize, &
                   FullPart, LSize, &
                   IWidth, Seed, dirname, middlename, IErr)

           END SELECT

        END DO ! ISeed cycle

     END DO  ! flux cycle

     DEALLOCATE ( HAMMAT, EIGS, WORK, CubeSites, LiebSites )

  END DO ! IWidth cycle

  STOP 'LiebExactDiag()'

END PROGRAM





