!-----------------------------------------------------------------
!
! LEDdiag
!
!-----------------------------------------------------------------
! exact diagonalization of 2D and 3D extended Lieb models
! see https://github.com/DisQS/LiebExactDiag
!-----------------------------------------------------------------

!-----------------------------------------------------------------
SUBROUTINE MakeLiebMatrixStructrue(dm, nu, n, ucl, n_uc, nt, matr, cubesites, liebsites)

  USE MyNumbers
  USE IChannels
  
  IMPLICIT NONE

  INTEGER(KIND=IKIND) &
       dm, & ! the dimension
       n, &  ! the number of unit cell in each dimension
       nu, & ! the number of site in each edge
       ucl, & ! the number of atoms in a unit cell
       nt, &  ! the whole number of atoms in system
       n_uc   ! the number of unit cells

  INTEGER(KIND=IKIND) i, j, k, ind, cubecount,liebcount

  LOGICAL(KIND=8) Flag

  INTEGER(KIND=IKIND), ALLOCATABLE :: ucl_d(:) 
  REAL(KIND=RKIND) matr(nt, nt)! , matr_W( nt, nt )
  INTEGER(KIND=IKIND) cubesites(n_uc), liebsites(nt-n_uc)
  
  PRINT*,"MakeLiebMatrixStructure()"

  matr(:,:) = 0.0D0

  SELECT CASE(dm)
  CASE(2)
     ALLOCATE(ucl_d(dm))
     ucl_d(1) = 2          ! The first Lieb atoms
     ucl_d(2) = nu + 2
  CASE(3)
     ALLOCATE(ucl_d(dm))
     ucl_d(1) = 2
     ucl_d(2) = nu + 2
     ucl_d(3) = 2 * nu + 2
  CASE DEFAULT
     PRINT*, "We only finished the 2D and 3D cases for the Lieb models --- ABORTING!"
     STOP
  END SELECT
    
  ! ----------------------------------------------------------
  ! Start construction of Lieb Matrix
  ! ----------------------------------------------------------
  
  cubecount=0
  liebcount=0

  DO i=1,n_uc
     ! ----------------------- Cube atoms ------------------------!
     ind = (i-1)*ucl + 1

     matr(ind,ind) = 0.0D0
     !print*,"Cube site", ind
     cubecount= cubecount+1
     cubesites(cubecount)=ind

     ! hopping terms for Cube atoms 
     DO j=1,dm

        IF(j==1)THEN
           Flag = ( MOD(i,n)==1 )
        ELSE IF(j==2)THEN
           Flag = ( MOD(i,n**2)>0 .AND. MOD(i,n**2)<=n )
        ELSE
           Flag = ( i<=n**2 )
        END IF

        matr(ind, (i-1)*ucl + ucl_d(j)) = 1.0D0
        matr((i-1)*ucl + ucl_d(j), ind) = 1.0D0

        IF(Flag)THEN
           matr(ind, (i-1)*ucl + ucl_d(j) + (nu-1) - ucl*(n)**(j-1) + ucl*(n)**j) = 1.0D0
           matr((i-1)*ucl + ucl_d(j) + (nu-1) - ucl*(n)**(j-1) + ucl*(n)**j, ind) = 1.0D0
        ELSE
           matr(ind, (i-1)*ucl + ucl_d(j) + (nu-1) - ucl*(n)**(j-1)) = 1.0D0
           matr((i-1)*ucl + ucl_d(j) + (nu-1) - ucl*(n)**(j-1), ind) = 1.0D0
        END IF

     END DO
     ! ----------------------- Lieb atoms ------------------------!

     ! ------- For Lieb atoms except close to Cube atom of other unit cell-----------

     IF(nu>1)THEN
        DO k=1, nu-1

           DO j=1, dm

              ind = (i-1)*ucl + ucl_d(j) + k - 1
              !print*,"Lieb site", ind
              liebcount= liebcount+1
              liebsites(liebcount)=ind

              matr(ind,ind) = 0.0D0
              !Hopping term
              IF(k==1)THEN
                 matr(ind, ind - (j-1)*nu -1) = 1.0D0
                 matr(ind - (j-1)*nu -1, ind) = 1.0D0                
              ELSE
                 matr(ind, ind - 1) = 1.0D0
                 matr(ind - 1, ind) = 1.0D0
              END IF
              matr(ind, ind + 1) = 1.0D0
              matr(ind + 1, ind) = 1.0D0

           END DO

        END DO
     END IF

     ! For Lieb atoms close to Cube atoms of other unit cells

     DO j=1,dm
        IF(j==1)THEN
           Flag = ( MOD(i,n**j)==0 )
        ELSE IF(j==2)THEN
           Flag = ( MOD(i,n**j)==0 .OR. (MOD(i,n**j).GT.(n**j-n)) )
        ELSE
           Flag = ( i.GT.(n**j-n**2) )
        END IF

        ind = (i-1)*ucl + ucl_d(j) + nu - 1
        !print*,"Lieb site", ind
        liebcount= liebcount+1
        liebsites(liebcount)=ind

        matr(ind, ind) = 0.0D0 

        IF(nu==1) THEN
           matr(ind, ind - ucl_d(j) +1) = 1.0D0
           matr(ind - ucl_d(j) +1, ind) = 1.0D0
        ELSE
           matr(ind, ind - 1) = 1.0D0
           matr(ind - 1, ind) = 1.0D0
        END IF

        IF(Flag)THEN
           matr(ind, (i-1)*ucl + 1 - (n-1)*ucl*(n)**(j-1)) = 1.0D0
           matr((i-1)*ucl + 1 - (n-1)*ucl*(n)**(j-1), ind) = 1.0D0
        ELSE
           matr(ind, (i-1)*ucl + 1 + ucl*(n)**(j-1)) = 1.0D0
           matr((i-1)*ucl + 1 + ucl*(n)**(j-1), ind) = 1.0D0
        END IF
     END DO

  END DO
!!!!! End construct Lieb Matrix

  RETURN

END SUBROUTINE MakeLiebMatrixStructrue

!-----------------------------------------------------------------
SUBROUTINE MakeLiebOnsiteDisorderORIGINAL(dm, nu, n, ucl, n_uc, nt, matr, cubesites, liebsites)

  USE MyNumbers
  USE IChannels
  
  IMPLICIT NONE

  INTEGER(KIND=IKIND) &
       dm, & ! the dimension
       n, &  ! the number of unit cell in each dimension
       nu, & ! the number of site in each edge
       ucl, & ! the number of atoms in a unit cell
       nt, &  ! the whole number of atoms in system
       n_uc   ! the number of unit cells

  INTEGER(KIND=IKIND) i, j, k, ind, cubecount,liebcount

  LOGICAL(KIND=8) Flag

  INTEGER(KIND=IKIND), ALLOCATABLE :: ucl_d(:) 
  REAL(KIND=RKIND) matr(nt, nt)! , matr_W( nt, nt )
  INTEGER(KIND=IKIND) cubesites(n_uc), liebsites(nt-n_uc)
  
  PRINT*,"MakeLiebOnsiteDisorderORIGINAL()"

  !matr(:,:) = 0.0D0

  ! ----------------------------------------------------------
  ! Start construction for the DIAGONAL elements of the Lieb Matrix
  ! ----------------------------------------------------------
  
  DO i=1,n_uc

     ! ----------------------- Cube atoms ------------------------!
     ind = (i-1)*ucl + 1

     matr(ind,ind) = 0.0D0

     !print*,"Cube site", ind
     cubecount= cubecount+1
     cubesites(cubecount)=ind

     ! ----------------------- Lieb atoms ------------------------!

     ! ------- For Lieb atoms except close to Cube atom of other unit cell-----------
     IF(nu>1)THEN
        DO k=1, nu-1

           DO j=1, dm

              ind = (i-1)*ucl + ucl_d(j) + k - 1
              !print*,"Lieb site", ind
              liebcount= liebcount+1
              liebsites(liebcount)=ind

              matr(ind,ind) = 0.0D0
           END DO

        END DO
     END IF

     ! For Lieb atoms close to Cube atoms of other unit cells

     DO j=1,dm
        IF(j==1)THEN
           Flag = ( MOD(i,n**j)==0 )
        ELSE IF(j==2)THEN
           Flag = ( MOD(i,n**j)==0 .OR. (MOD(i,n**j).GT.(n**j-n)) )
        ELSE
           Flag = ( i.GT.(n**j-n**2) )
        END IF

        ind = (i-1)*ucl + ucl_d(j) + nu - 1
        !print*,"Lieb site", ind
        liebcount= liebcount+1
        liebsites(liebcount)=ind

        matr(ind, ind) = 0.0D0 
     END DO

  END DO
!!!!! End construct Lieb Matrix

  RETURN

END SUBROUTINE MakeLiebOnsiteDisorderORIGINAL

!-----------------------------------------------------------------
SUBROUTINE MakeLiebOnsiteDisorder(dm, nu, n, ucl, n_uc, nt, matr, size, seed)

  USE MyNumbers
  USE IChannels
  USE IPara
  USE DPara

  USE RNG_MT

  IMPLICIT NONE

  INTEGER(KIND=IKIND) &
       dm, & ! the dimension
       n, &  ! the number of unit cell in each dimension
       nu, & ! the number of site in each edge
       ucl, & ! the number of atoms in a unit cell
       nt, &  ! the whole number of atoms in system
       n_uc   ! the number of unit cells

  INTEGER(KIND=IKIND) i, j, k, ind, cubecount,liebcount

  LOGICAL(KIND=8) Flag

  INTEGER(KIND=IKIND), ALLOCATABLE :: ucl_d(:) 
  REAL(KIND=RKIND) matr(nt, nt)! , matr_W( nt, nt )
  !INTEGER(KIND=IKIND) cubesites(n_uc), liebsites(nt-n_uc)

  INTEGER seed(5), size, len,col,row 
  REAL(KIND=RKIND) drandval, SUMCUBErandval,SUMLIEBrandval
  !EXTERNAL DRANDOM5
  
  PRINT*,"MakeLiebOnsiteDisorder()"

  !matr(:,:) = 0.0D0

  ! ----------------------------------------------------------
  ! Start construction for the DIAGONAL elements of the Lieb Matrix
  ! ----------------------------------------------------------

  SUMCUBErandval= 0.0D0; SUMLIEBrandval= 0.0D0

  ! Give the Lieb matrix different onsite potentials
  DO i=1, n_uc

     ! CUBE onsite potentials
     drandval= DRANDOM5(seed)

     SUMCUBErandval=SUMCUBErandval + CubeDis*(drandval - 0.5D0)

     SELECT CASE (IRNGFlag)
     CASE(0) ! constant CubeConPot on each cube site
        matr( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = CubeConPot + CubeDis*(drandval - 0.5D0)
     CASE(1) ! +/- CubeConPot on random cube sites

        IF(MOD(n_uc,2) .NE. 0) THEN
           PRINT*, "main: WRNG, cube size is odd, so we cannot achieve 0 SUMMED potential!"
           PRINT*, "main: WRNG, calculation will proceed, but output is questionable."
        END IF

        IF(DRANDOM5(seed)>=0.5) THEN
           matr( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = &
                CubeConPot + CubeDis*(drandval - 0.5D0)
        ELSE
           matr( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = &
                -CubeConPot + CubeDis*(drandval - 0.5D0)
        END IF

     CASE(2) ! checkerboard +/- CubeConPot on each cube site

        IF(MOD(n_uc,2) .NE. 0) THEN
           PRINT*, "main: WRNG, cube size are odd, so we can not achieve 0 potential!"
           PRINT*, "main: WRNG, calculation will proceed, but output is questionable."
        END IF

        IF(dm==2)THEN

           len=1
           IF(MOD(i,size)==0)THEN
              col=size
              row=INT(i/size)
           ELSE
              col=MOD(i,size)
              row=INT(i/size)+1
           END IF

        ELSE IF(Dim==3)THEN

           IF(MOD(i,size**2)==0)THEN
              col=size
              row=size
              len=INT(i/size**2)
           ELSE
              IF(MOD(MOD(i,size**2),size)==0)THEN
                 col=size
                 row=INT(MOD(i,size**2)/size)
              ELSE
                 col=MOD(MOD(i,size**2),size)
                 row=INT(MOD(i,size**2)/size)+1
              END IF
              len=INT(i/size**2)+1
           END IF

        ELSE
           PRINT*,"Not finished yet!"
           STOP
        END IF

        drandval= DRANDOM5(seed)
        SUMCUBErandval=SUMCUBErandval + CubeDis*(drandval - 0.5D0)

        IF(MOD(len,2)==1)THEN

           IF((MOD(col,2)==1 .AND. MOD(row,2)==1) .OR. (MOD(col,2)==0 .AND. MOD(row,2)==0))THEN
              matr( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = CubeConPot + CubeDis*(drandval - 0.5D0)
           ELSE
              matr( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = -1.0*CubeConPot + CubeDis*(drandval - 0.5D0)
           END IF

        ELSE
           IF((MOD(col,2)==1 .AND. MOD(row,2)==1) .OR. (MOD(col,2)==0 .AND. MOD(row,2)==0))THEN
              matr( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = -1.0*CubeConPot + CubeDis*(drandval - 0.5D0)
           ELSE
              matr( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = CubeConPot + CubeDis*(drandval - 0.5D0)
           END IF

        END IF

     END SELECT

     ! LIEB onsite potentials
     DO j=1, ucl-1

        drandval= DRANDOM5(seed)
        SUMLIEBrandval=SUMLIEBrandval + LiebDis*(drandval - 0.5D0)
        matr((i-1)*ucl + j + 1 , (i-1)*ucl + j + 1) = LiebConPot + LiebDis*(drandval - 0.5D0)

     END DO

  END DO

  ! ----------------------------------------------------------
  ! WRITE SUMrandval to allow identification of accidental states
  ! ----------------------------------------------------------

  PRINT*,"MakeLiebOnsiteDisorder(): Seed=", seed, ", SHrv=", SUMCUBErandval/n_uc, &
       ", SRrv=", SUMLIEBrandval/n_uc

  RETURN

END SUBROUTINE MakeLiebOnsiteDisorder

!-----------------------------------------------------------------
SUBROUTINE MakeLiebHoppingDisorder(dm, nu, n, ucl, n_uc, nt, matr, cubesites, liebsites)

  USE MyNumbers
  USE IChannels
  
  IMPLICIT NONE

  INTEGER(KIND=IKIND) &
       dm, & ! the dimension
       n, &  ! the number of unit cell in each dimension
       nu, & ! the number of site in each edge
       ucl, & ! the number of atoms in a unit cell
       nt, &  ! the whole number of atoms in system
       n_uc   ! the number of unit cells

  INTEGER(KIND=IKIND) i, j, k, ind, cubecount,liebcount

  LOGICAL(KIND=8) Flag

  INTEGER(KIND=IKIND), ALLOCATABLE :: ucl_d(:) 
  REAL(KIND=RKIND) matr(nt, nt)! , matr_W( nt, nt )
  INTEGER(KIND=IKIND) cubesites(n_uc), liebsites(nt-n_uc)
  
  PRINT*,"MakeLiebMatrixStructure()"

  !matr(:,:) = 0.0D0

  ! ----------------------------------------------------------
  ! Start construction for the OFF-DIAGONAL elements of the Lieb Matrix
  ! ----------------------------------------------------------
  
  cubecount=0
  liebcount=0

  DO i=1,n_uc
     ! ----------------------- Cube atoms ------------------------!
     ind = (i-1)*ucl + 1

     cubecount= cubecount+1
     cubesites(cubecount)=ind

     ! hopping terms for Cube atoms 
     DO j=1,dm

        IF(j==1)THEN
           Flag = ( MOD(i,n)==1 )
        ELSE IF(j==2)THEN
           Flag = ( MOD(i,n**2)>0 .AND. MOD(i,n**2)<=n )
        ELSE
           Flag = ( i<=n**2 )
        END IF

        matr(ind, (i-1)*ucl + ucl_d(j)) = 1.0D0
        matr((i-1)*ucl + ucl_d(j), ind) = 1.0D0

        IF(Flag)THEN
           matr(ind, (i-1)*ucl + ucl_d(j) + (nu-1) - ucl*(n)**(j-1) + ucl*(n)**j) = 1.0D0
           matr((i-1)*ucl + ucl_d(j) + (nu-1) - ucl*(n)**(j-1) + ucl*(n)**j, ind) = 1.0D0
        ELSE
           matr(ind, (i-1)*ucl + ucl_d(j) + (nu-1) - ucl*(n)**(j-1)) = 1.0D0
           matr((i-1)*ucl + ucl_d(j) + (nu-1) - ucl*(n)**(j-1), ind) = 1.0D0
        END IF

     END DO
     ! ----------------------- Lieb atoms ------------------------!

     ! ------- For Lieb atoms except close to Cube atom of other unit cell-----------

     IF(nu>1)THEN
        DO k=1, nu-1

           DO j=1, dm

              ind = (i-1)*ucl + ucl_d(j) + k - 1
              !print*,"Lieb site", ind
              liebcount= liebcount+1
              liebsites(liebcount)=ind

              !Hopping term
              IF(k==1)THEN
                 matr(ind, ind - (j-1)*nu -1) = 1.0D0
                 matr(ind - (j-1)*nu -1, ind) = 1.0D0                
              ELSE
                 matr(ind, ind - 1) = 1.0D0
                 matr(ind - 1, ind) = 1.0D0
              END IF
              matr(ind, ind + 1) = 1.0D0
              matr(ind + 1, ind) = 1.0D0

           END DO

        END DO
     END IF

     ! For Lieb atoms close to Cube atoms of other unit cells

     DO j=1,dm
        IF(j==1)THEN
           Flag = ( MOD(i,n**j)==0 )
        ELSE IF(j==2)THEN
           Flag = ( MOD(i,n**j)==0 .OR. (MOD(i,n**j).GT.(n**j-n)) )
        ELSE
           Flag = ( i.GT.(n**j-n**2) )
        END IF

        ind = (i-1)*ucl + ucl_d(j) + nu - 1
        !print*,"Lieb site", ind
        liebcount= liebcount+1
        liebsites(liebcount)=ind

        IF(nu==1) THEN
           matr(ind, ind - ucl_d(j) +1) = 1.0D0
           matr(ind - ucl_d(j) +1, ind) = 1.0D0
        ELSE
           matr(ind, ind - 1) = 1.0D0
           matr(ind - 1, ind) = 1.0D0
        END IF

        IF(Flag)THEN
           matr(ind, (i-1)*ucl + 1 - (n-1)*ucl*(n)**(j-1)) = 1.0D0
           matr((i-1)*ucl + 1 - (n-1)*ucl*(n)**(j-1), ind) = 1.0D0
        ELSE
           matr(ind, (i-1)*ucl + 1 + ucl*(n)**(j-1)) = 1.0D0
           matr((i-1)*ucl + 1 + ucl*(n)**(j-1), ind) = 1.0D0
        END IF
     END DO

  END DO
!!!!! End construct Lieb Matrix

  RETURN

END SUBROUTINE MakeLiebHoppingDisorder
