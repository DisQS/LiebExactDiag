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

  IF(dm==2)THEN
     ALLOCATE(ucl_d(dm))
     ucl_d(1) = 2          ! The first Rim atoms
     ucl_d(2) = nu + 2
  ELSE IF(dm==3)THEN
     ALLOCATE(ucl_d(dm))
     ucl_d(1) = 2
     ucl_d(2) = nu + 2
     ucl_d(3) = 2 * nu + 2
  ELSE
     PRINT*, "We Only Finished the 2D and 3D cases for Lieb model"
     STOP
  END IF
    
  ! ----------------------------------------------------------
  ! Start construction of Lieb Matrix
  ! ----------------------------------------------------------
  
  cubecount=0
  liebcount=0

  DO i=1,n_uc
     ! ----------------------- hub atoms ------------------------!
     ind = (i-1)*ucl + 1

     matr(ind,ind) = 0.0D0
     !print*,"Cube site", ind
     cubecount= cubecount+1
     cubesites(cubecount)=ind

     ! hopping terms for hub atoms 
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
     ! ----------------------- Rim atoms ------------------------!

     ! ------- For rim atoms except close to hub atom of orther unit cell-----------

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


     ! For rim atoms close to hub atoms of other unit cells

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

!!$SUBROUTINE WriteEvals(dm, nu, n, nt, HubDiagDis, RimDiagDis, W, matr_W, norm, part_nr, ISample, INFO)
!!$
!!$  IMPLICIT NONE
!!$  
!!$  INTEGER, PARAMETER :: IKIND = SELECTED_INT_KIND(9)
!!$  INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(15,307)
!!$
!!$  INTEGER(KIND=IKIND) dm, nu, n, nt, ISample, INFO
!!$  INTEGER(KIND=IKIND) i, j
!!$  REAL(KIND=RKIND) HubDiagDis, RimDiagDis
!!$  REAL(KIND=RKIND) W( nt ), matr_W( nt, nt ), norm( nt ), part_nr( nt )
!!$  
!!$  CHARACTER*100 FileName
!!$
!!$  WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A1,A2,I4.4,A1,A2,I4.4,A1,I4.4,A4)')&
!!$       "Eval-","L",dm,nu,&
!!$       "-M",n, "-",&
!!$       "WH", NINT(100.D0*ABS(HubDiagDis)),&
!!$       "-","WR", NINT(100.D0*ABS(RimDiagDis)),&
!!$       "-",ISample,".raw"
!!$
!!$  PRINT*, "FileName: ", FileName
!!$
!!$  OPEN(Unit=9, FILE=FileName)
!!$
!!$  IF(INFO==0)THEN
!!$
!!$     DO i=1,nt
!!$        DO j=1,nt	
!!$           norm(i) = norm(i) + matr_W(i,j)**2
!!$           part_nr(i) = part_nr(i) + matr_W(i,j)**4
!!$        END DO
!!$     END DO
!!$
!!$     DO i=1,nt
!!$        WRITE(9,'(2f30.20)') W(i) , (norm(i)**2) / part_nr(i)
!!$     END DO
!!$
!!$  ELSE
!!$     PRINT*, "ERROR IN CALL DSYEV()" 
!!$  END IF
!!$
!!$  CLOSE(9)
!!$
!!$  RETURN
!!$       
!!$END SUBROUTINE WRITEEVALS


