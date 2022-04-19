!.....................................................................!
! Revised main routine for LAMP 
! (Linear algebra Angular Momentum Projection)
! Current version 1.5.4 Oct 2021 (with edits in Nov 2021)
! Starting with version 1.3.0 Nov 2016
! by CWJ based upon prior work by Joshua Staker and Kevin O'Mara
!.....................................................................!
USE psis 
USE sporbit 
USE spstate 
USE lamplight
USE managers 
USE hamlib 
USE wignersfriend 
USE eigenpackage 
USE dens1body 
USE outerputer 

! NOTE FOR SML
! Need to revise timing module for MPI vs non-MPI use

! Added by SML
! Contains variables for MPI environment including: myMPIrank, nMPIranks, icomm & root
USE nodeinfo ! LAMPmodule_parallel.f90 (modelled after Sherpa parallel module)
IMPLICIT NONE 

CHARACTER :: achar, ychar, menu_char, eval_char 

INTEGER :: ierr 																							! Error message variable for MPI 
INTEGER :: ij
INTEGER :: num_threads, mythread 															! for OpenMP
INTEGER :: omp_get_thread_num, omp_get_num_threads 						! for OpenMP
INTEGER(KIND=4) :: numsdused 																	! allows user to reduce number of SDs used
INTEGER(KIND=4) :: nprint 																		! how many to print out 
INTEGER(KIND=8) :: nlevelmax 																	! # of levels possible 
INTEGER(KIND=8) :: nftotal 																		! actual number of levels

REAL(KIND=4) :: tolerance, elapsedTime, Jtolerance 
REAL(KIND=4) :: newJmax 
REAL(KIND=8) :: clock_start, clock_stop 											! timing variables
REAL(KIND=8) :: norm_start, norm_stop, ham_start, ham_stop 		! timing variables
REAL(KIND=8) :: proj_start, proj_stop
REAL(KIND=8) :: x 
REAL(KIND=8) :: normSum, hamSum 

usepnform = .FALSE. 

! Set up MPI environment 
! Variables below defined in LAMP_modulesParallel.f90 (module nodeinfo) except ierr
! I might try adding ierr to that module so I don't have to define it in 
! each module that makes use of MPI
CALL MPI_INIT(icomm)
icomm = MPI_COMM_WORLD 
CALL MPI_COMM_RANK(icomm, myMPIrank, ierr)
CALL MPI_COMM_SIZE(icomm, nMPIranks, ierr)
root = 0 

IF (myMPIrank == root) THEN 
	PRINT*, '  WELCOME to LAMP'
	PRINT*, '  Version 1.5.4, Oct 2021'
	PRINT*, '  Authors: C. W. Johnson, S. M. Lauber, C. F. Jiao'
	PRINT*, '  based upon work by J. Staker and K. O Mara'
	PRINT*, ''
END IF  

! Initialize OpenMP information 
!$OMP PARALLEL SHARED(num_threads)
num_threads = omp_get_num_threads()
mythread = omp_get_thread_num() 
num_global_threads = num_threads
IF (mythread == 0 .AND. myMPIrank == root) THEN  
	WRITE(6,*) ' NUM_THREADS = ', num_threads 
	WRITE(6,*) ' NUM MPI RANKS = ', nMPIranks
END IF 
!$OMP END PARALLEL 

! Start to count calculating time 
IF (myMPIrank == root) CALL CPU_TIME(clock_start)
CALL get_orb_info							! located in LAMPsplib.f90
CALL unfold_spstates					! located in LAMPsplib.f90
CALL set_nuclide 							! located in LAMPpsilib.f90
CALL factorial(-70.d0,x)			! Just to initialize
CALL allocateSlaterDet 				! located in LAMPpsilib.f90 
CALL inputJmax 								! located in LAMPutils.f90  
CALL default_Jmesh 						! likely run independently, located in LAMPutils.f90 SIMPLE!

numOfBeta = numOfJ 
doHam = .FALSE. 

IF (myMPIrank == root) CALL CPU_TIME(norm_start)			! might need to modify or protect for MPI runs
CALL projectorator 												! located in LAMPmanagelib.f90, keep an eye on some subroutines
CALL tracemaster 													! located in LAMPoutput.f90
IF (myMPIrank == root) CALL printNorm(6) 	! located in LAMPoutput.f90
IF (myMPIrank == root) CALL CPU_TIME(norm_stop)			! might need to modify or protect for MPI runs

IF (myMPIrank == root) THEN 
	PRINT*, ''
	PRINT*, ' Choose one of the following: '
	PRINT*, ' (e) Compute energies or similar operator ONLY'
	PRINT*, ' (d) 1-body densities (will compute energies first)'
	READ(5,'(A)') menu_char 
END IF 

CALL MPI_BARRIER(icomm,ierr)
CALL MPI_BCAST(menu_char,1,MPI_CHARACTER,root,icomm,ierr)

SELECT CASE(menu_char)
CASE DEFAULT ! menu_char
	IF (myMPIrank == root) THEN 
		PRINT*, ' That choice is not an option'
		PRINT*, ' Please choose (e) eneriges or (d) 1-body densities'
		STOP 
	END IF ! myMPIrank 

!.... SOLVE FOR ENERGIES 
CASE('e','E') ! menu_char
	! Getting stuck somewhere in this subroutine
	CALL ham_boss ! Read in Hamiltonian information, in LAMP_hamlib.f90
	IF (myMPIrank == root) THEN 
		PRINT*, ' Enter the tolerance for accept Js (typically ~0.00001 to 0.0001)'
		READ*, Jtolerance 
	END IF ! myMPIrank == root 
	CALL MPI_BARRIER(icomm,ierr)
	CALL MPI_BCAST(Jtolerance,1,MPI_INTEGER,root,icomm,ierr)
	CALL findNewJmax(Jtolerance,newJmax)	! LAMPmanagelib.f90
	CALL default_Jmesh 										! LAMPutils.f90
	numOfBeta = numOfJ 
	IF (myMPIrank == root) CALL CPU_TIME(ham_start)

!.... CALCULATE NORM AND HAMILTONIAN MATRICES 
	doHam = .TRUE. 
	CALL deallocator 			! LAMPutils.f90
	CALL hmultMPIdistro		! LAMP_hamlib,f90 

	! Added timing around Hamiltonian computation
	IF (myMPIrank == root) CALL CPU_TIME(proj_start) 
	CALL projectorator 		! LAMPmanagelib.f90
	IF (myMPIrank == root) CALL CPU_TIME(proj_stop)

	CALL tracemaster 			! LAMPoutput.f90 

	IF (ychar == 'y' .OR. ychar == 'Y' .AND. myMPIrank == root) THEN 
		CALL printNorm(6)		! LAMPoutput.f90
	END IF 

	IF (myMPIrank == root) CALL CPU_TIME(clock_stop)
	ham_stop = clock_stop 

!.... ADDED IN 1.4.4, JULY 2019 BY CWJ 
!     OPTIONS TO WRITE TO FILE / COMPUTE EXPECTATION VALUE 
	IF (myMPIrank == root) THEN 
		PRINT*, '' 
		PRINT*, ' Evaluation options: '
		PRINT*, ' (e) (default) compute energies'
		PRINT*, ' (w) write matrices to file'
		PRINT*, ' (x) read matrices and compute expectation values'
		PRINT*, '' 
		READ(5,'(A)') eval_char 
	END IF ! myMPIrank == root 
	CALL MPI_BARRIER(icomm,ierr)
	CALL MPI_BCAST(eval_char,1,MPI_CHARACTER,root,icomm,ierr)

	SELECT CASE(eval_char)
	CASE('w','W') ! eval_char
		PRINT*, ''
		PRINT*, ' Writing matrices to file'
		PRINT*, '' 
		CALL write_mat2file 

		PRINT*, ''
		PRINT*, ' Do you want to continue (y/n)?'
		READ(5,'(A)') ychar 

		IF (ychar == 'n' .OR. ychar == 'N') THEN 
			PRINT*, '' 
			WRITE(*,*) '    TIME of total LAMP calculation ', clock_stop - clock_start, ' seconds'
			WRITE(*,*) '    TIME of norm calculation: ', norm_stop - norm_start, ' seconds'
			WRITE(*,*) '    TIME of ham calculation: ', ham_stop - ham_start, ' seconds'
			PRINT*, '' 
			PRINT*, ' Goodbye!'
			STOP 
		END IF 

	CASE('x','X') ! eval_char
		PRINT*, '' 
		PRINT*, ' Evalulation of expectation values'
		PRINT*, '' 

		compute_expect = .TRUE. 
		CALL read_matfromfile 
	CASE DEFAULT ! e
		IF (myMPIrank == root) THEN 
			PRINT*, '' 
			PRINT*, ' Computing energies and eigenvectors'
			PRINT*, ''
		END IF ! myMPIrank == root 

		compute_expect = .FALSE. 
	END SELECT ! eval_char

!.... SOLVE GENERALIZED EIGENVALUE PROBLEM 
!     FIRST, count up how many levels we could have 
	nlevelmax = 0 
	DO ij =  1, numOfJ 
		nlevelmax = nlevelmax + NINT(2 * xJlist(ij) + 1)
	END DO 
	nlevelmax = nlevelmax * numsd 
	ALLOCATE(jall(nlevelmax),pallPair(2,nlevelmax),obsall(2,nlevelmax))
	ALLOCATE(problist(2,numOfJ),hamlist(2,numOfJ))
	
	IF (myMPIrank == root) THEN 
		PRINT*, ' How many energies to print to screen? (all will be written to file)'
		READ*, nprint 
	END IF ! myMPIrank == root 
	CALL MPI_BARRIER(icomm,ierr)
	CALL MPI_BCAST(nprint,1,MPI_INT,root,icomm,ierr)

!.... ADDED IN 1.4.2: ABILITY TO DO CHANGES IN TOLERANCE 
	ychar = 'y' 
	CALL inputTolerance(tolerance) ! located in managelilbs
	CALL MPI_BARRIER(icomm,ierr)
	CALL MPI_BCAST(tolerance,1,MPI_REAL,root,icomm,ierr)

	jtarget = -1. 
	numsdused = numsd 
	!--------------------------------------------------
	! MIGHT NEED TO ROOT PROTECT THIS ENTIRE WHILE LOOP
	!--------------------------------------------------
	DO WHILE (ychar == 'y' .OR. ychar == 'Y')		
		CALL EigenSolverPackage(tolerance, nlevelmax, numsdused, nftotal, jall, pallPair, obsall, normSum, hamSum, problist, hamlist)	! in LAMPeigen.f90 
		PRINT*, '' 
		IF (myMPIrank == root) THEN 
			PRINT*, ' Sum of norms = ', DBLE(normSum)
			PRINT*, ' Sum of trace(H) = ', DBLE(hamSum)

			CALL J_WriteResults(tolerance, nlevelmax, nftotal, jall, pallPair, obsall, problist, hamlist, .NOT.allsameparity, nprint)		! in LAMPoutput.f90
			PRINT*, ''
			PRINT*, ' Do you want to run with other parameters (y/n)?'
			READ(5,'(A)') ychar 
		END IF ! myMPIrank == root 
		CALL MPI_BARRIER(icomm,ierr)
		CALL MPI_BCAST(ychar,1,MPI_CHARACTER,root,icomm,ierr)

		IF ( ychar == 'y' .OR. ychar == 'Y' ) THEN 
			CALL inputTolerance(tolerance)
			IF ( numsd > 1 ) THEN 
				PRINT*, ' Enter number of SDs to use (max', numsd, ')'
				READ*, numsdused 
				IF (numsdused > numsd) numsdused = numsd 
			END IF 

			PRINT*, ' Choose a target j to print out (y/n)?'
			READ(5,'(A)') achar 
			IF (achar == 'y' .OR. achar == 'Y') THEN 
				PRINT*, ' Enter target J '
				READ*, jtarget 
			ELSE 
				jtarget = -1 
			END IF ! achar
		END IF ! ychar 
	END DO ! ychar 

	IF (myMPIrank == root) THEN
		WRITE(*,*) '    TIME of total LAMP calculation: ', clock_stop - clock_start, ' seconds'
		WRITE(*,*) '    TIME of norm calculation: ', norm_stop - norm_start, ' seconds'
		WRITE(*,*) '    TIME of ham calculation: ', ham_stop - ham_start, ' seconds'
		WRITE(*,*) '    TIME of projectorator calculation: ', proj_stop - proj_start, ' seconds'
	END IF ! myMPIrank == root 

CASE('d','D') ! menu_char
	PRINT*, ' Must first find energies'
!.... Read in Hamiltonian information 
	CALL ham_boss 
	! REPLACED BELOW BY CWJ RECOMMENDATION
	! PRINT*, ' Enter the tolerance for accept Js (typically ~0.00001 to 0.0001)'
	! READ*, Jtolerance 

	! CWJ RECOMMENDED CORRECTION
	IF (myMPIrank == root) THEN 
		PRINT*, ' Enter the tolerance for accept Js (typically ~0.00001 to 0.0001)'
		READ*,  Jtolerance 
	END IF ! myMPIrank == root 
	CALL MPI_BARRIER(icomm,ierr)
	CALL MPI_BCAST(Jtolerance,1,MPI_INTEGER,root,icomm,ierr)
	! CWJ RECOMMENDED CORRECTION 

	CALL findNewJmax(Jtolerance,newJmax)
	CALL default_Jmesh 
	numOfBeta = numOfJ 
	CALL CPU_TIME(ham_start)
!.... Calculate norm and Hamiltonian matrices 
	doHam = .TRUE. 
	CALL deallocator 
	CALL projectorator 
	CALL tracemaster 

	IF (ychar == 'y' .OR. ychar == 'Y') THEN 
		CALL printNorm(6) 
	END IF 

	CALL CPU_TIME(clock_stop)
	ham_stop = clock_stop 

	CALL write_mat2file 

	finddense = .TRUE. 

	CALL default_Jmesh 
	numOfBeta = numOfJ 
	CALL deallocator 
	CALL projectorator_dens 
END SELECT ! menu_char

CALL MPI_BARRIER(icomm,ierr)
CALL MPI_FINALIZE(icomm,ierr) 						! TEMPORARY LOCATION FOR TESTING

END ! of main program 
