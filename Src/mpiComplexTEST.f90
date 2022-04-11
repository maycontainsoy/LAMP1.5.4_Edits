!.....................................................................!
! This program is not associated with LAMP but was created to test the
! passing of complex arrays and the differences between MPI types on 
! the broadcast function. This did not solve the issue I am currently 
! experiencing with allocateSlaterDet so there is something else I am 
! missing. 
!
! In this test, there does not appear to be a difference (at least no
! errors are produced during compilation and the results match at run)
! between COMPLEX8 and COMPLEX16 in MPI when the vector is defined as 
! a COMPLEX8 in Fortran.
!.....................................................................!

PROGRAM mpiComplexTEST
  IMPLICIT NONE

	INCLUDE 'mpif.h'

	INTEGER :: nMPIranks, myMPIrank, icomm 
	INTEGER :: root 
	INTEGER :: ierr

	INTEGER :: vecSize, i

	REAL(KIND=4) :: a, b

	COMPLEX(KIND=4), ALLOCATABLE :: pvec4(:)
	COMPLEX(KIND=8), ALLOCATABLE :: pvec8(:)

	CALL MPI_INIT(icomm)
	icomm = MPI_COMM_WORLD
	CALL MPI_COMM_RANK(icomm, myMPIrank, ierr)
	CALL MPI_COMM_SIZE(icomm, nMPIranks, ierr)
	root = 0

	IF (myMPIrank == root) THEN 
		PRINT*, ' Enter vecSize: '
		READ*, vecSize 
	END IF 

	CALL MPI_BCAST(vecSize,1,MPI_INTEGER,root,icomm,ierr)

	ALLOCATE(pvec4(vecSize))
	ALLOCATE(pvec8(vecSize))

	IF (myMPIrank == root) THEN 
		DO i = 1, vecSize
			CALL random_number(a)
			CALL random_number(b) 
			pvec4(i) = CMPLX(a,b)
			pvec8(i) = CMPLX(a,b)
		END DO ! i
	END IF ! myMPIrank == root

	CALL MPI_BARRIER(icomm,ierr)
	CALL MPI_BCAST(pvec4,vecSize,MPI_COMPLEX16,root,icomm,ierr)
	CALL MPI_BARRIER(icomm,ierr)

	CALL MPI_BARRIER(icomm,ierr)
	CALL MPI_BCAST(pvec8,vecSize,MPI_COMPLEX16,root,icomm,ierr)
	CALL MPI_BARRIER(icomm,ierr)

	PRINT*, 'Node # ', myMPIrank, ' vec4: ', pvec4(:)
	PRINT*, 'Node # ', myMPIrank, ' vec8: ', pvec8(:)

	CALL MPI_FINALIZE(icomm,ierr)

END PROGRAM mpiComplexTEST