
!.....................................................................!
! Contains modules: 
! 	sporbit - single particle orbit information 
! 	spstate - single particle states from obit
! Both modules contain subroutines that pull from the .spo/.sps 
! file to populate derived type arrays.
!.....................................................................!

!.....................................................................!
! Module: sporbit
!.....................................................................!
MODULE sporbit 
	! Single-particle orbit information
	USE nodeinfo 
	IMPLICIT NONE

	INTEGER :: numorb 	! # of s.p. orbits 
	INTEGER :: numruns 

!.... CREATE A DEFINED TYPE 
	TYPE orb
		INTEGER :: nr 			! radial quantum number
		INTEGER :: j 				! 2xj 
		INTEGER :: l 				! L 
		INTEGER :: w 				! excitation 
		INTEGER :: spstart 	! where these orbits correspond to ! start in spqn 
		INTEGER :: par 			! parity 
	END TYPE orb 

	LOGICAL :: allsameparity ! NOTE in 1.5.1 allsameparity replaces PAIRTEST 
	TYPE (orb), ALLOCATABLE, TARGET :: orbqn(:)

CONTAINS 

!.....................................................................!
! Subroutine: get_orb_info
! Reads in s.p. space information;
! fills in quantum numbers into orbqn and spsqn
!
! Reads in either .spo or .sps files 
! it first looks for .spo file; if that fails, 
! automatically looks for .sps file of same name 
!.....................................................................!
SUBROUTINE get_orb_info 
	IMPLICIT NONE 

!.... FILE CONTROL
	CHARACTER(LEN=25) :: filename 
	CHARACTER(LEN=1)  :: achar 

	CHARACTER(LEN=70) :: title 
	INTEGER :: ilast 

!.... TEMP 
	REAL :: xn, xl, xj, xw ! orbital quantum numbers 

!.... DUMMY COUNTERS
	INTEGER :: i, j, k, m 
	LOGICAL :: success 
	INTEGER :: isp 	! counter for single particle states?
	INTEGER :: ierr ! for MPI
	INTEGER :: aerr ! Allocation error

!.... OPEN A FILE 
	success = .FALSE. 
	IF (myMPIrank == root) THEN 
		DO WHILE (.NOT.success) 
			PRINT*,' Enter file with s.p. orbit information (.spo/.sps)'
			READ(5,'(A)') filename 
			ilast = index(filename,' ') - 1 
!.... ATTEMPT TO OPEN .spo FILE 
			OPEN(UNIT=1,FILE=filename(1:ilast)//'.spo',STATUS='OLD',ERR=101)
			success = .TRUE. 
			CYCLE 
101		CONTINUE
!.... ATEMPT TO OPEN .sps FILE 
			OPEN(UNIT=1,FILE=filename(1:ilast)//'.sps',STATUS='OLD',ERR=102)
			success = .TRUE. 
			CYCLE 
102		CONTINUE 
			PRINT*,filename(1:ilast),'.spo/.sps file does not exist'
		END DO 

!.... READ PAST TITLE CARDS
		success = .FALSE. 
		DO WHILE(.NOT.success)
			READ(1,'(A)') achar 
			BACKSPACE(1)
			IF (achar /= '#' .AND. achar /= '!') THEN 
				success = .TRUE. 
			ELSE 
				READ(1,'(A)') title 
				WRITE(6,*) title 
			END IF ! searching for title card 
		END DO ! while .not.success
!.... READ PAST POSSIBLE LABEL OF ISO/PN
		READ(1,'(A)') achar 
		IF (achar == 'p' .OR. achar == 'P') THEN 
			PRINT*, ' .sps file in pn formalism, cannot handle'
			PRINT*, ' program ended in file: LAMPsplib.f90'
			PRINT*, ' subroutine: get_orb_info'
			STOP 
		ELSEIF (achar /= 'i' .AND. achar /= 'I') THEN 
			BACKSPACE(1)	
		END IF ! checking for format 

!.... READ NUMBER OF ORBITS 
		READ(1,*) numorb 
	END IF ! myMPIrank == root

	CALL MPI_BARRIER(icomm,ierr)
	CALL MPI_Bcast(numorb,1,MPI_INT,root,icomm,ierr)
	! CALL MPI_Barrier(icomm,ierr)

	! DEBUGGING, REMOVE WHEN FINISHED
	! PRINT*, 'My MPI rank = ', ' number of orbits = ', numorb
	! DEBUGGING, REMOVE WHEN FINISHED 

!.... ALLOCATE MEMORY FOR QUANTUM NUMBERS turbo
	!IF (myMPIrank == root) THEN
	IF (numruns .NE. 0) THEN 
		DEALLOCATE(orbqn)
	END IF 
	ALLOCATE(orbqn(numorb), stat=aerr)
	IF (aerr /= 0) PRINT*, 'Error allocating array'
	!END IF ! myMPIrank == root
	!CALL MPI_Barrier(icomm,ierr)

!.... READ IN 
	isp = 1 
	allsameparity = .TRUE. 
	DO i = 1, numorb 
		IF (myMPIrank == root) THEN 
			READ(1,*,END=2001) xn, xl, xj, xw
			! PRINT*, xn, xl, xj, xw
			orbqn(i)%nr      = INT(xn)
			orbqn(i)%l       = INT(xl)
			orbqn(i)%j       = INT(2*xj)
			orbqn(i)%w       = INT(xw)
			orbqn(i)%spstart = isp 
			isp = isp + orbqn(i)%j + 1 
			orbqn(i)%par     = (-1)**orbqn(i)%l
		END IF ! myMPIrank == root

		CALL MPI_BARRIER(icomm,ierr)
		CALL MPI_BCAST(orbqn(i)%nr,1,MPI_INT,root,icomm,ierr)
		CALL MPI_BCAST(orbqn(i)%l, 1,MPI_INT,root,icomm,ierr)
		CALL MPI_BCAST(orbqn(i)%j, 1,MPI_INT,root,icomm,ierr)
		CALL MPI_BCAST(orbqn(i)%w, 1,MPI_INT,root,icomm,ierr)
		CALL MPI_BCAST(orbqn(i)%spstart,1,MPI_int,root,icomm,ierr)
		CALL MPI_BCAST(isp, 1, MPI_INT,root,icomm,ierr)
		! Note sure if I need to broadcast isp or just the array information
		CALL MPI_BCAST(orbqn(i)%par,1,MPI_INT,root,icomm,ierr)

		IF (myMPIrank == root) THEN 
			IF (i > 1 .AND. orbqn(i)%par /= orbqn(1)%par) allsameparity = .FALSE.
		END IF ! myMPIrank
		CALL MPI_Bcast(allsameparity,1,MPI_logical,root,icomm,ierr)
	END DO ! i = 1, numorb 
	
	IF (myMPIrank == root) CLOSE(UNIT=1) ! .spo/.sps file 
	RETURN

!.... ERROR TRAP FOR END OF FILE 
2001	CONTINUE 
	IF (myMPIrank == root) THEN
		PRINT*, ' Sudden end of file in ', filename(1:ilast)
		PRINT*, ' Error trap occurs in LAMPsplib.f90'
		PRINT*, ' Subroutine: get_orb_info'
		CALL MPI_ABORT(icomm,101,ierr)
		STOP 1234 
	END IF ! myMPIrank == root 

	RETURN 
END SUBROUTINE get_orb_info 
!.....................................................................!

END MODULE sporbit 

!.....................................................................!
! Module: spstate 
! Single particle state information 
!.....................................................................!
MODULE spstate 
	USE sporbit ! module defined above 
	USE nodeinfo
	IMPLICIT NONE 

	INTEGER :: nsps ! number of s.p. states 

!.... CREATE A DEFINED TYPE 
	TYPE spst
		INTEGER :: nr 		! radial quantum number 
		INTEGER :: j 			! 2xj 
		INTEGER :: m 			! 2xjz 
		INTEGER :: l 			! L 
		INTEGER :: w 			! excitation 
		INTEGER :: par 		! parity 
		INTEGER :: orb 		! orbit label 
	END TYPE spst

	TYPE (spst), ALLOCATABLE :: spsqn(:)

CONTAINS 

!.....................................................................!
! Subroutine: unfold_spstate 
! Creates single particle states from single particle orbits 
!.....................................................................!
SUBROUTINE unfold_spstates 
	IMPLICIT NONE 

	INTEGER :: i, j, k, m 
	INTEGER :: ierr ! for MPI communication error

	nsps = 0 
	DO i = 1, numorb 
		IF (myMPIrank == root) THEN 
			nsps = nsps + orbqn(i)%j+1
		END IF 
	END DO 

	CALL MPI_BARRIER(icomm,ierr)
	CALL MPI_BCAST(nsps,1,MPI_INT,root,icomm,ierr)
	!PRINT*, 'Node: ', myMPIrank, ' nsps: ', nsps
	!CALL MPI_BARRIER(icomm,ierr)

	ALLOCATE(spsqn(nsps))
	k = 0 
	DO i = 1, numorb 
		j = orbqn(i)%j 
		DO m = -j,j,2 
			k = k + 1 
			IF (myMPIrank == root) THEN 
				spsqn(k)%nr  = orbqn(i)%nr 
				spsqn(k)%l   = orbqn(i)%l 
				spsqn(k)%par = (-1)**(orbqn(i)%l)
				spsqn(k)%j   = orbqn(i)%j 
				spsqn(k)%w   = orbqn(i)%j 
				spsqn(k)%m   = m 
				spsqn(k)%orb = i 
				spsqn(k)%par = orbqn(i)%par 
			END IF ! myMPIrank

			! Update other nodes 
			CALL MPI_BARRIER(icomm,ierr)
			CALL MPI_BCAST(spsqn(k)%nr,  1,MPI_INT,root,icomm,ierr)
			CALL MPI_BCAST(spsqn(k)%l,   1,MPI_INT,root,icomm,ierr)
			CALL MPI_BCAST(spsqn(k)%par, 1,MPI_INT,root,icomm,ierr)
			CALL MPI_BCAST(spsqn(k)%j,   1,MPI_INT,root,icomm,ierr)
			CALL MPI_BCAST(spsqn(k)%w,   1,MPI_INT,root,icomm,ierr)
			CALL MPI_BCAST(spsqn(k)%m,   1,MPI_INT,root,icomm,ierr)
			CALL MPI_BCAST(spsqn(k)%orb, 1,MPI_INT,root,icomm,ierr)
		END DO ! m 
	END DO ! i 

	! Checking that nodes are receiving the same information
	! PRINT*, 'node, k_last = ', myMPIrank, k, ' nr  = ', spsqn(:)%nr 
	! PRINT*, 'node, k_last = ', myMPIrank, k, ' l   = ', spsqn(:)%l 
	! PRINT*, 'node, k_last = ', myMPIrank, k, ' par = ', spsqn(:)%par 
	! PRINT*, 'node, k_last = ', myMPIrank, k, ' j   = ', spsqn(:)%j 
	! PRINT*, 'node, k_last = ', myMPIrank, k, ' w   = ', spsqn(:)%w
	! PRINT*, 'node, k_last = ', myMPIrank, k, ' m   = ', spsqn(:)%m 
	! PRINT*, 'node, k_last = ', myMPIrank, k, ' orb = ', spsqn(:)%orb 

		RETURN 
END SUBROUTINE unfold_spstates 
!.....................................................................!

END MODULE spstate 