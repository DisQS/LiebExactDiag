!
!
!   Call function DSYEV() to calculate eigenvalues and eigenvectors
!   for Lieb matrix and its extendsions(2D and 3D)
!
!
!--------------------------------------------------------------------------------------

PROGRAM Lieb

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
  
  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE:: CubeProb, CubePart, LiebProb, LiebPart

  INTRINSIC        INT, MIN
  EXTERNAL         DSYEV

  ! Parameters for eigenverctor, participation numbers
  
  INTEGER(KIND=IKIND) Seed, i, j, Inum, IErr
  REAL(KIND=RKIND) drandval
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
  PRINT*,"LiebExactDiag (", TRIM("GITVERSION"), ", ", TRIM("GITBRANCH"), ", compiled: ", TRIM("COMPILED"), ")"
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
     
     ALLOCATE ( HAMMAT0(LSize, LSize) )
     ALLOCATE ( EIGS( LSize ) )
     ALLOCATE ( WORK( LWMAX ) )
     ALLOCATE ( HAMMAT( LSize, LSize ) )
     ALLOCATE ( CubeSites( n_uc ) )
     ALLOCATE ( LiebSites( LSize - n_uc ) )
     ALLOCATE( CubeProb( LSize ) )
     ALLOCATE( LiebProb( LSize ) )
     ALLOCATE( CubePart( LSize ) )
     ALLOCATE( LiebPart( LSize ) )
     !ALLOCATE ( norm(LSize) )
     !ALLOCATE ( part_nr(LSize) )
 
     HAMMAT0(:,:) = 0.0D0

     ! ----------------------------------------------------------
     IF(IWriteFlag.GE.1) THEN
        PRINT*, "--- starting matrix build"
     ENDIF

     CALL MakeLiebMatrixStructrue(Dim, Nx, IWidth, ucl, n_uc, LSize, HAMMAT0, CubeSites, LiebSites)

     !print*,"CubeSites=", CubeSites
     !print*,"LiebSites=", LiebSites

!!$
!!$     END DO
!!$     DO i= 1, LSize
!!$        write(*,108) i
!!$        Do j= 1, LSize 
!!$           IF( matr(i,j).ne.(0.0) )Then
!!$              Write(*,108) j
!!$           END IF
!!$           
!!$        END DO
!!$        write(*,*)" "
!!$     END DO
!!$108  format(1x,1I3\)
     
     DO HubDis= HubDis0,HubDis1,dHubDis

        PRINT*,"main: HubDis=", HubDis

        CALL GetDirec(Dim, Nx, IWidth, HubDis, RimDis, 0.0, str)

        DO Seed=ISeed, ISeed+NSeed-1

           ! ----------------------------------------------------------
           ! Compute actual seed
           ! ----------------------------------------------------------
           
           !CALL SRANDOM(Seed)
           
           ISSeed(1)= Seed
           ISSeed(2)= IWidth
           ISSeed(3)= NINT(HubDis*1000.) ! in lieu of "Energy"
           ISSeed(4)= NINT(HubDis*1000.)
           ISSeed(5)= NINT(RimDis*1000.)
           
!           CALL genrand_int31(ISSeed) ! MT95 with 5 seeds

           SELECT CASE(IWriteFlag)
           CASE(1,2)
              PRINT*, "-- Seed=", Seed
              PRINT*, "-> ISSeed=", ISSeed
           CASE(3,4)
!!$                 PRINT*, "IS: IW=", IWidth, "hD=", NINT(HubDis*1000.), "E=", NINT(Energy*1000.), &
!!$                      "S=", Seed, "IS=", ISSeed
              CALL genrand_int31(ISSeed) ! MT95 with 5 seeds
              CALL genrand_real1(drandval)
              CALL SRANDOM5(ISSeed)
              drandval=DRANDOM5(ISSeed)
              WRITE(*, '(A7,I3,A4,F6.3,A4,F5.3,A3,I5,A4,F16.10)') &
                   "IS: IW=", IWidth, " hD=", HubDis, " rD=", RimDis, &
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
              CALL CheckOutput( Dim,Nx, IWidth, 0.0, HubDis, RimDis, &
                   Seed, str, IErr )
              IF(IErr.EQ.2) CYCLE
           END SELECT
           
           !CALL genrand_int31(ISSeed) ! MT95 with 5 seeds, before: CALL SRANDOM(ISSeed
           CALL SRANDOM5(ISSeed) ! MT95 with 5 seeds, before: CALL SRANDOM(ISSeed)

           ! ----------------------------------------------------------
           ! ENTER random values into matrix
           ! ----------------------------------------------------------
              
           !norm(:) = 0d0
           !part_nr(:) = 0d0

           HAMMAT(:,:) = HAMMAT0(:,:)

           ! Give the Lieb matrix different onsite potensial
           DO i=1, n_uc

              drandval= DRANDOM5(ISSeed)
              HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = HubDis*(drandval - 0.5D0)

              DO j=1, ucl-1

                 drandval= DRANDOM5(ISSeed)
                 HAMMAT((i-1)*ucl + j + 1 , (i-1)*ucl + j + 1) = RimDis*(drandval - 0.5D0)

              END DO

           END DO

           ! ----------------------------------------------------------
           ! START the diagonalizstion process
           ! ----------------------------------------------------------

           PRINT*, "STARTing the diagonalizstion process"
              
           LWORK =  -1  !3*LSize

           CALL DSYEV( 'V', 'Upper', LSize, HAMMAT, LSize, EIGS, WORK, LWORK, INFO )

           LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

           CALL DSYEV( 'V', 'Upper', LSize, HAMMAT, LSize, EIGS, WORK, LWORK, INFO ) 

           IF(INFO.NE.0) THEN
              PRINT*,"main: DYSEV(INFO)=", INFO
              STOP
           ENDIF

           ! ----------------------------------------------------------
           ! WRITE the eigenvalues and -vectors
           ! ----------------------------------------------------------

           !CALL WriteEvals(Dim, Nx, IWidth, LSize, HubDis, RimDis, EIGS, HAMMAT, norm, part_nr, Seed, INFO)

           NEIG=LSize ! this is complete diagonalization

           CALL WriteOutputEVal( Dim, Nx, NEIG, EIGS, &
                IWidth, 0., HubDis, RimDis, Seed, str, IErr)
           SELECT CASE(IStateFlag)
           CASE(0)
              CONTINUE
           CASE(1)
              PRINT*,"main: DYSEV() eigenvectors will now be saved into individual files"
              DO Inum= 1,NEIG
                 Call WriteOutputEVec(Dim, Nx, Inum, NEIG, Lsize, &
                      HAMMAT, LSize, IWidth, 0.0, HubDis, & 
                      RimDis, Seed, str, IErr)
              END DO
           CASE(2)
              PRINT*,"main: DYSEV() eigenvectors will now be saved into single BULK file"

              Call WriteOutputEVecBULK(Dim, Nx, Lsize, NEIG, Lsize, &
                   EIGS, LSize, IWidth, 0.0, HubDis, & 
                   RimDis, Seed, str, IErr)
           CASE(-1)
              PRINT*,"main: Cube/Lieb site projections of DYSEV() eigenvectors"

              ! compute projection and participation numbers for Cube sites
              DO Inum= 1,NEIG

                 CubeProb(Inum)= 0.0
                 CubePart(Inum)= 0.0

                 DO i=1,n_uc
                    CubeProb(Inum)= CubeProb(Inum) + &
                         HAMMAT(Inum, CubeSites(i)) * HAMMAT(Inum, CubeSites(i))
                    CubePart(Inum)= CubePart(Inum) + &
                         HAMMAT(Inum, CubeSites(i)) * HAMMAT(Inum, CubeSites(i)) * &
                         HAMMAT(Inum, CubeSites(i)) * HAMMAT(Inum, CubeSites(i)) 
                 END DO
                 CubePart(Inum)=n_uc*n_uc/CubePart(Inum)/LSize/LSize

                 LiebProb(Inum)= 0.0
                 LiebPart(Inum)= 0.0

                 DO i=1,LSize-n_uc
                    LiebProb(Inum)= LiebProb(Inum) + &
                         HAMMAT(Inum, LiebSites(i)) * HAMMAT(Inum, LiebSites(i))
                    LiebPart(Inum)= LiebPart(Inum) + &
                         HAMMAT(Inum, LiebSites(i)) * HAMMAT(Inum, LiebSites(i)) * &
                         HAMMAT(Inum, LiebSites(i)) * HAMMAT(Inum, LiebSites(i))
                 END DO
                 LiebPart(Inum)=(LSize-n_uc)*(LSize-n_uc)/LiebPart(Inum)/LSize/LSize

              ENDDO

              Call WriteOutputEVecProj(Dim, Nx, Inum, NEIG, &
                   EIGS, LSize, &
                   CubeProb, CubePart, Lsize, &
                   LiebProb, LiebPart, LSize, &
                   IWidth, 0.0, HubDis, & 
                   RimDis, Seed, str, IErr)
              
           END SELECT

        END DO ! ISeed cycle

     END DO  ! Disorder cycle

     DEALLOCATE ( HAMMAT0, EIGS, WORK, HAMMAT, CubeSites, LiebSites )

  END DO ! IWidth cycle

END PROGRAM Lieb





