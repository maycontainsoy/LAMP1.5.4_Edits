subroutine complex_inverse(M,A,Ainv,info)
	Implicit none
	integer, intent(IN)  :: M   
	complex (kind=8),dimension(1:M,1:M), intent(IN) ::A
        complex (kind=8),dimension(1:M,1:M), intent(OUT) ::Ainv
	complex (kind=8),dimension(1:M,1:M) ::WORK
	integer,dimension(1:M)::IPIV
	integer :: i,j,error
	integer,intent(OUT) :: info
	complex (kind=8) :: zdet
	
	Ainv = A

	call ZGETRF(M,M,Ainv,M,IPIV,info)
	
	if(info .eq. 0) then
!	  write(*,*)"LU factorization succeded"
	else
	 write(*,*)"complex LU factorization failed"

	 
         return
	end if
 zdet = (1.0,0.0)
 do i = 1,M
	 zdet = zdet*Ainv(i,i)
 end do
! print*,' det = ',zdet
	
	call ZGETRI(M,Ainv,M,IPIV,WORK,M,info)
	if(info .eq. 0) then
!	  write(*,*)"inverse succeded"
	else
	 write(*,*)"complex inverse failed"
         return
	end if

!	call ComplexMatrixReform(M,Ainv)
	return

end subroutine complex_inverse

subroutine ComplexMatrixReform(dim, Mat)

	IMPLICIT NONE
	
	INTEGER, INTENT(IN) :: dim
	COMPLEX (kind=8), DIMENSION(1:dim,1:dim), INTENT(INOUT) :: Mat
	
	INTEGER :: i, j
	REAL (kind=8) :: eps = 1.0D-10, max, min
	COMPLEX (kind=8), PARAMETER :: IUNIT = (0.0D0, 1.0D0)
	
	max = MAXVAL(ABS(Mat(1:dim,1:dim)))
	min = MINVAL(ABS(Mat(1:dim,1:dim)))
	
	DO i = 1, dim
	DO j = 1, dim
	
	   IF( ABS(Mat(i,j)) .LT. eps .OR. ABS(Mat(i,j)) .LT. max*1.0D-10 ) THEN
	
	      Mat(i,j) = (0.0D0,0.0D0)
	
	   ELSE IF (ABS(AIMAG(Mat(i,j))) .LT. eps ) THEN
	
	      Mat(i,j) = DBLE(Mat(i,j))
	
	   ELSE IF (ABS(DBLE(Mat(i,j))) .LT. eps ) THEN
	
	      Mat(i,j) = IUNIT * AIMAG(Mat(i,j))
	
	   END IF
	
	END DO
	END DO
	
	RETURN

end subroutine ComplexMatrixReform


subroutine double_inverse(M,A,Ainv,info)
	Implicit none
	integer, intent(IN)  :: M   
	real (kind=8),dimension(1:M,1:M), intent(IN) ::A
    real (kind=8),dimension(1:M,1:M), intent(OUT) ::Ainv
	real (kind=8),dimension(1:M,1:M) ::WORK
	integer,dimension(1:M)::IPIV
	integer :: i,j,error
	integer,intent(OUT) :: info

	Ainv = A

	call DGETRF(M,M,Ainv,M,IPIV,info)
	if(info .eq. 0) then
	  write(*,*)"LU factorization succeded"
	else
	 write(*,*)"double LU factorization failed"
         return
	end if
	
	call DGETRI(M,Ainv,M,IPIV,WORK,M,info)
	if(info .eq. 0) then
	  write(*,*)"inverse succeded"
	else
	 write(*,*)"inverse failed"
         return
	end if

	return

end subroutine double_inverse
