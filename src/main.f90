!
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
  Integer(KIND=IKIND), DIMENSION(:), ALLOCATABLE:: CRA

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
     ALLOCATE( CRA( n_uc ) )
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
     
     DO CubeDis= CubeDis0,CubeDis1,dCubeDis

        PRINT*,"main: Cubedis-loop, CubeDis=", CubeDis

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
              PRINT*, "-- Seed=", Seed
              PRINT*, "-> ISSeed=", ISSeed
           CASE(3,4)
!!$                 PRINT*, "IS: IW=", IWidth, "hD=", NINT(CubeDis*1000.), "E=", NINT(Energy*1000.), &
!!$                      "S=", Seed, "IS=", ISSeed
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
           ! ENTER random values into matrix
           ! ----------------------------------------------------------
              
           !norm(:) = 0d0
           !part_nr(:) = 0d0

           HAMMAT(:,:) = HAMMAT0(:,:)
           SUMHUBrandval= 0.0D0; SUMRIMrandval= 0.0D0
           
           ! Give the Lieb matrix different onsite potensial

           If(IRNGFlag==0)then ! constant CubeConPot on each cube site
              DO i=1, n_uc
                 
                 drandval= DRANDOM5(ISSeed)
                 SUMHUBrandval=SUMHUBrandval + CubeDis*(drandval - 0.5D0)
                 
                 HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = CubeConPot + CubeDis*(drandval - 0.5D0)
                 
                 DO j=1, ucl-1
                    
                    drandval= DRANDOM5(ISSeed)
                    SUMRIMrandval=SUMRIMrandval + LiebDis*(drandval - 0.5D0)
                    HAMMAT((i-1)*ucl + j + 1 , (i-1)*ucl + j + 1) = LiebConPot + LiebDis*(drandval - 0.5D0)
                    
                 END DO
                 
              END DO
              
           Else if(IRNGFlag==1)then ! +/- CubeConPot on random cube site
              DO i=1, n_uc
                 
                 drandval= DRANDOM5(ISSeed)
                 SUMHUBrandval=SUMHUBrandval + CubeDis*(drandval - 0.5D0)
                 
                 HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = CubeConPot + CubeDis*(drandval - 0.5D0)
                 
                 DO j=1, ucl-1
                    
                    drandval= DRANDOM5(ISSeed)
                    SUMRIMrandval=SUMRIMrandval + LiebDis*(drandval - 0.5D0)
                    HAMMAT((i-1)*ucl + j + 1 , (i-1)*ucl + j + 1) = LiebConPot + LiebDis*(drandval - 0.5D0)
                    
                 END DO
              END DO
                            
              IF(Mod(n_uc,2) .ne. 0) Then
                 Print*, "Cube size are odd can not achieve 0 potential! Stop!!!!"
                 Stop
              End IF

              ! Add the checkerboard pattern to the cube sites
              Call www_fcode_cn( CRA, n_uc)              
              DO i=1, n_uc/2
                 drandval= DRANDOM5(ISSeed)
                 HAMMAT( (CRA(i)-1)*ucl + 1 , (CRA(i)-1)*ucl + 1 ) = -1.0D0*CubeConPot + CubeDis*(drandval - 0.5D0)
              END DO

           Else if(IRNGFlag==2) then ! checkerboard +/- CubeConPot on each cube site
              DO i=1, n_uc

                 If(Dim==2)Then

                    len=1
                    If(Mod(i,IWidth)==0)Then
                       col=IWidth
                       row=Int(i/IWidth)
                    Else
                       col=Mod(i,IWidth)
                       row=Int(i/IWidth)+1
                    End if

                 Else if(Dim==3)Then
                    
                    If(Mod(i,IWidth**2)==0)Then
                       col=IWidth
                       row=IWidth
                       len=Int(i/IWidth**2)
                    Else
                       If(Mod(Mod(i,IWidth**2),IWidth)==0)Then
                          col=IWidth
                          row=Int(Mod(i,IWidth**2)/IWidth)
                       Else
                          col=Mod(Mod(i,IWidth**2),IWidth)
                          row=Int(Mod(i,IWidth**2)/IWidth)+1
                       End if                      
                       len=Int(i/IWidth**2)+1
                    End If

                 Else
                    Print*,"Not finish yet!"
                    Stop
                 End If
                 
                 drandval= DRANDOM5(ISSeed)
                 SUMHUBrandval=SUMHUBrandval + CubeDis*(drandval - 0.5D0)

                 If(Mod(len,2)==1)Then

                    If((Mod(col,2)==1 .and. Mod(row,2)==1) .or. (Mod(col,2)==0 .and. Mod(row,2)==0))Then
                       HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = CubeConPot + CubeDis*(drandval - 0.5D0)
                    Else
                       HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = -1.0*CubeConPot + CubeDis*(drandval - 0.5D0)
                    End If

                 Else
                    If((Mod(col,2)==1 .and. Mod(row,2)==1) .or. (Mod(col,2)==0 .and. Mod(row,2)==0))Then
                       HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = -1.0*CubeConPot + CubeDis*(drandval - 0.5D0)
                    Else
                       HAMMAT( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = CubeConPot + CubeDis*(drandval - 0.5D0)
                    End If

                 End If
                 
                 DO j=1, ucl-1
                    
                    drandval= DRANDOM5(ISSeed)
                    SUMRIMrandval=SUMRIMrandval + LiebDis*(drandval - 0.5D0)
                    HAMMAT((i-1)*ucl + j + 1 , (i-1)*ucl + j + 1) = LiebConPot + LiebDis*(drandval - 0.5D0)
                    
                 END DO
                 
!!$                 If(Int(n_uc/IWidth**(Dim-1)) .le. 1))
                 
              END DO
              
              
           Else
              Print*,"Not finish yet!!!"
              Stop
           End If

           
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
                 Call WriteOutputEVec(Dim, Nx, Inum, NEIG, Lsize, HAMMAT, LSize, &
                      IWidth, CubeDis, LiebDis, CubeConPot, LiebConPot,&
                      Seed, str, IErr)
              END DO
           CASE(2)
              PRINT*,"main: DYSEV() eigenvectors will now be saved into single BULK file"

              Call WriteOutputEVecBULK(Dim, Nx, Lsize, NEIG, Lsize, EIGS, LSize, &
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

              Call WriteOutputEVecProj(Dim, Nx, Inum, NEIG, &
                   EIGS, LSize, &
                   CubeProb, CubePart, Lsize, &
                   LiebProb, LiebPart, LSize, &
                   FullPart, LSize, &
                   IWidth, CubeDis, LiebDis, &
                   CubeConPot, LiebConPot, &
                   Seed, str, IErr)
              
           END SELECT

        END DO ! ISeed cycle

     END DO  ! Disorder cycle

     DEALLOCATE ( HAMMAT0, CRA, EIGS, WORK, HAMMAT, CubeSites, LiebSites )

  END DO ! IWidth cycle
  
  STOP 'LiebExactDiag'

END PROGRAM LiebExactDiag





