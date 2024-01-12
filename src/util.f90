!-----------------------------------------------------------------
!
! LEDdiag
!
!-----------------------------------------------------------------
! exact diagonalization of 2D and 3D extended Lieb models
! see https://github.com/DisQS/LiebExactDiag
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

SUBROUTINE MakeLiebOnsiteDisorder(dm, nu, n, ucl, n_uc, nt, matr, cubesites, liebsites)

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
  
  PRINT*,"MakeLiebOnsiteDisorder()"

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

END SUBROUTINE MakeLiebOnsiteDisorder

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
