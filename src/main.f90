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

  REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE:: MATRIX0

  ! Parameters for call function DSYEV()
  
  INTEGER(KIND=IKIND) LDA, LWMAX, INFO, LWORK

  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE:: W, WORK
  REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE:: MATRIX

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
  ! set: git tag -a v0.0 -m 'Version 0.0'
  ! ----------------------------------------------------------
#ifdef git
  PRINT*,"LiebExactDiag (", TRIM("GITVERSION"), ")"
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
     
     ucl = (Dim * Nx) + 1
     n_uc = IWidth**Dim
     LSize = ucl * n_uc

     !--------------------------------------------------------------------------
     ! Setting the parameters size passed to function DSYEV
     !--------------------------------------------------------------------------

     LDA = LSize
     LWMAX = 100000
     
     ! ----------------------------------------------------------
     ! ALLOCATing memory
     ! ----------------------------------------------------------
     
     ALLOCATE ( MATRIX0(LSize, LSize) )
     ALLOCATE ( w( LSize ) )
     ALLOCATE ( WORK( LWMAX ) )
     ALLOCATE ( MATRIX( LSize, LSize ) )
     ALLOCATE ( norm(LSize) )
     ALLOCATE ( part_nr(LSize) )
 

     MATRIX0(:,:) = 0.0D0
     MATRIX(:,:) = 0.0D0

     ! ----------------------------------------------------------
     IF(IWriteFlag.GE.1) THEN
        PRINT*, "--- starting matrix build"
     ENDIF

     CALL MakeLiebMatrixStructrue(Dim, Nx, IWidth, ucl, n_uc, LSize, MATRIX0)

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
     
     norm(:) = 0d0
     part_nr(:) = 0d0

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
              CALL CheckOutput( Dim,Nx, IWidth, HubDis, RimDis, &
                   Seed, str, IErr )
              IF(IErr.EQ.2) CYCLE
           END SELECT
           
           !CALL genrand_int31(ISSeed) ! MT95 with 5 seeds, before: CALL SRANDOM(ISSeed
           CALL SRANDOM5(ISSeed) ! MT95 with 5 seeds, before: CALL SRANDOM(ISSeed)

           ! ----------------------------------------------------------
           ! ENTER random values into matrix
           ! ----------------------------------------------------------
              
           MATRIX(:,:) = MATRIX0(:,:)

           ! Give the Lieb matrix different onsite potensial
           DO i=1, n_uc

              drandval= DRANDOM5(ISSeed)
              MATRIX( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = HubDis*(drandval - 0.5D0)

              DO j=1, ucl-1

                 drandval= DRANDOM5(ISSeed)
                 MATRIX((i-1)*ucl + j + 1 , (i-1)*ucl + j + 1) = RimDis*(drandval - 0.5D0)

              END DO

           END DO

           ! ----------------------------------------------------------
           ! START the diagonalizstion process
           ! ----------------------------------------------------------
              
           LWORK =  -1  !3*LSize

           CALL DSYEV( 'V', 'Upper', LSize, MATRIX, LSize, W, WORK, LWORK, INFO )

           LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

           CALL DSYEV( 'V', 'Upper', LSize, MATRIX, LSize, W, WORK, LWORK, INFO )      

           ! ----------------------------------------------------------
           ! WRITE the eigenvalues and -vectors
           ! ----------------------------------------------------------

           !CALL WriteEvals(Dim, Nx, IWidth, LSize, HubDis, RimDis, W, MATRIX, norm, part_nr, Seed, INFO)

           NEIG=LSize ! this is complete diagonalization

           CALL WriteOutputEVal( Dim, Nx, NEIG, MATRIX, &
                IWidth, 0., HubDis, RimDis, Seed, str, IErr)
           IF(IStateFlag.NE.0) THEN
              PRINT*,"main: DYSEV() found eigenvectors, these will now be saved into file"
              DO Inum= 1,NEIG
                 Call WriteOutputEVec(Dim, Nx, Inum, NEIG, Lsize, &
                      MATRIX, LSize, IWidth, HubDis, & 
                      RimDis, Seed, str, IErr)
              END DO
           END IF !IStateFlag IF

           MATRIX(:,:) = 0d0
           norm(:) = 0d0
           part_nr(:) = 0d0

        END DO ! ISeed cycle

     END DO  ! Disorder cycle

     DEALLOCATE ( MATRIX0, w, WORK, MATRIX, norm, part_nr )

  END DO ! IWidth cycle

END PROGRAM Lieb





