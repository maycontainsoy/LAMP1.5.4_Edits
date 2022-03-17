module outerputer
	USE nodeinfo
	implicit none
    real (kind=8), allocatable :: jall(:), pallPair(:,:),obsall(:,:)!!
    real (kind=8),  allocatable  :: jlist(:)!!
    real (kind=8), allocatable  :: problist(:,:), hamlist(:,:),obslist(:,:)!!
contains

!==============================================================================
!
! prints formatted problist to screen or file
!
!  INPUT:
!   numOfJ = # of beta values -> int(Jmax + 1.)
!   isOdd = logical flag for odd A nuclides
!   ParityTest = logical, .true. if even/odd parity
!   nlist = length of problist, hamlist, and jlist; int(Jmax-Jmin)+1
!   problist(2,:) = fraction of original HF state in the norm for each J
!   printTo = target of write() statements
!               6 == to screen (standard output)
!               anything >=10 is an external file  
!
!  INTERNAL:
!   jlist(:) = list of floating point J values (i.e. 0.0, 1.0, ...)
!   probsum = sum of all elements in problist(2,:) should be 1.0 (or 2.0 if ParityTest)
!
	subroutine printNorm(printTo)
		use lamplight
		use tracy
		implicit none

		integer (kind=4) :: printTo

!............INTERNAL....................
		integer (kind=4) :: intJ
		real(kind=8) :: probsum,negnormsum
		character*70 :: filename
		integer      :: ilast 
	
		logical      :: print_fM = .false.
		real(kind=4) :: xJmax
	
		xJmax = xJlist(numOfJ)
	
!	print*,' local conversion ',xJmax+1
	
		if(printTo /= 6)then  ! open a file to write to
			print*, ' '
			print*,' Enter FULL name of output file for norm data '
			print*,' (Will overwrite existing file )'
			read(5,'(a)')filename
			ilast=index(filename,' ')-1
			open(unit=printTo,file=filename(1:ilast),status='unknown')
		end if

		write(printTo,*)' '
		write(printTo,*)' Fraction in original HF state: Norm'
		write(printTo,*)' '
		probsum = 0.d0
		negnormsum = 0.d0
		if(.not.allsameparity)then
			write(printTo,*)'   J    frac(+)   frac(-)'	
			write(printTo,*)'-------------------------'
			do intJ = 1,numOfJ
				write(printTo,2001)xJlist(intJ), normtrace(intJ),pnormtrace(intJ)
!            write(printTo,2001)jlist(intJ),problist(1,intJ),problist(2,intJ)
				probsum = probsum+normtrace(intJ)+pnormtrace(intJ)
				if(normtrace(intJ) < 0.d0)negnormsum = negnormsum-normtrace(intJ)
				if(pnormtrace(intJ) < 0.d0)negnormsum = negnormsum-pnormtrace(intJ)
			end do
2001 format(f6.1,2f10.6)
		
			if(print_fM)then
				write(printTo,*)' '
				write(printTo,*)'   M    frac(+)   frac(-)'
				write(printTo,*)'-------------------------'			
				do intJ = 1,J2max1
					write(printTo,2001)-xJmax+intJ-1,mnormtrace(intJ),pmnormtrace(intJ)
				end do			
			end if ! print_fM
		else ! .not.allsameparity
			write(printTo,*)'   J    frac'
			write(printTo,*)'-----------------------'
			do intJ = 1,numOfJ
				write(printTo,2001)xJlist(intJ),normtrace(intJ)
				probsum = probsum+normtrace(intJ)
				if(normtrace(intJ) < 0.d0)negnormsum = negnormsum-normtrace(intJ)
			end do ! intJ
		
			if(print_fM)then
				write(printTo,*)' '
				write(printTo,*)'  M/K    frac(M)     frac(K)'
				write(printTo,*)'-----------------------'			
				do intJ = 1,numOfM
!				print*,intJ-mshift,intJ-xJmax-1
					write(printTo,2001)-xJmax+intJ-1,mnormtrace(intJ),mknorm(intJ)
				end do ! intJ
			end if ! print_fM
		end if ! .not.allsameparity
		write(printTo,*)' Total of sum of trace(N) = ',probsum

		if(negnormsum > 0.001)then
			print*,' '
			print*,' ********** '
			print*,' ATTENTION:  significant fraction of negative norms ',negnormsum
			print*,' Try again with smaller and/or larger max J '
			print*,' ********** '
			print*,'   '
		end if

		if(printTo/=6)close(printTo)

		return 
	end subroutine printNorm
!==============================================================================
	subroutine tracemaster
		use nodeinfo 
		use lamplight
		use tracy
		use psis,only:numsd
	
		implicit none
		interface
	   	SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,  INFO ) ! INTERFACE
				CHARACTER          JOBZ, UPLO
				INTEGER            INFO, LDA, LWORK, N
				DOUBLE PRECISION   RWORK( * ), W( * )
				COMPLEX*16         A( LDA, * ), WORK( * )
			 END SUBROUTINE zheev      ! INTERFACE
		end interface	

		real(kind=4) :: floatK,floatM,xJ
		integer :: iJ,ii,kk,localdim
		real(kind=8) :: trace,ptrace,traceall,htrace
		real(kind=8) :: htraceall
		complex*16,allocatable :: overlapmat(:,:)
		real (kind = 8),ALLOCATABLE :: rwork(:),overeig(:)
		COMPLEX*16, ALLOCATABLE :: work(:)
		integer :: lwork,info
	
		real(kind=8) :: xJmax
		integer :: localoffset,im,ik
		integer :: sdX,sdY ! indices for multiple SDs

		INTEGER :: ierr ! For MPI messages

		traceall = 0.d0
		htraceall = 0.d0
	
		allocate(normtrace(numOfJ),hamtrace(numOfJ))
		allocate(overlapmat(numSD,numSD))
		overlapmat= 0.d0
	
		if(.not.allsameparity)allocate(pnormtrace(numOfJ))

		!PRINT*, ' Node = ', myMPIrank, ' makes it into tracemaster 1'	! TESTING, REMOVE
	
!........ SET UP ARRAYS mnormtrace,pmnormtrace	FOR f_M ......

		xJmax = xJlist(numOfJ)
		numOfM = nint(2*xJmax+1)
		!PRINT*, ' Node = ', myMPIrank, ' xJmax = ', xJmax, ' numofM = ', numOfM, ' tracemaster 2' ! TESTING, REMOVE
		if(numOfM/=J2max1 .AND. myMPIrank == root)then
			print*,' Huh, some mismatch ',J2max1,numOfM
!			stop
		end if ! includes myMPIrank == root

		if(allocated(mnormtrace))deallocate(mnormtrace)
		if(allocated(pmnormtrace))deallocate(pmnormtrace)
		allocate(mnormtrace(numOfM))
		if(allocated(mknorm))deallocate(mknorm)
		allocate(mknorm(numOfM))

		mnormtrace = 0.d0
		mknorm = 0.d0
		if(.not.allsameparity)then
			allocate(pmnormtrace(numOfM))
			if(allocated(pmknorm))deallocate(pmknorm)
			allocate(pmknorm(numofM))
			pmnormtrace= 0.d0
			pmknorm=0.d0
		end if
	
		do 	iJ = 1,numOfJ
			xJ = xJlist(iJ)
			localdim = nint(2*xJ+1)
			trace  = 0.d0
			ptrace = 0.d0
			htrace = 0.d0
			localoffset = nint(xJmax-xJ)
			do sdX= 1,numSD
				do ii = 1,localdim
					im = ii + localoffset
					if(.not.allsameparity)then
						trace  = trace  + 0.5* real(N_J_MK(sdX,sdX,iJ)%MK(ii,ii)+ PN_J_MK(sdX,sdX,iJ)%MK(ii,ii),kind=8)
						ptrace = ptrace + 0.5* real(N_J_MK(sdX,sdX,iJ)%MK(ii,ii)- PN_J_MK(sdX,sdX,iJ)%MK(ii,ii),kind=8)
						mnormtrace(im)  = mnormtrace(im)  + 0.5* real(N_J_MK(sdX,sdX,iJ)%MK(ii,ii)+ PN_J_MK(sdX,sdX,iJ)%MK(ii,ii),kind=8)
						pmnormtrace(im) = pmnormtrace(im) + 0.5* real(N_J_MK(sdX,sdX,iJ)%MK(ii,ii)- PN_J_MK(sdX,sdX,iJ)%MK(ii,ii),kind=8)
					else
!				print*,sdX,ij,ii,real(N_J_MK(sdX,sdX,iJ)%MK(ii,ii),kind=8)
						!PRINT*, ' Node = ', myMPIrank, ' tracemaster 3' ! TESTING, REMOVE
						trace = trace + real(N_J_MK(sdX,sdX,iJ)%MK(ii,ii),kind=8)
						mnormtrace(im)= mnormtrace(im)+ real(N_J_MK(sdX,sdX,iJ)%MK(ii,ii),kind=8)
						do kk=1,localdim
							mknorm(im)=mknorm(im)+  abs(N_J_MK(sdX,sdX,iJ)%MK(ii,kk))
						end do
						if(doHam)htrace = htrace + real(H_J_MK(sdX,sdX,iJ)%MK(ii,ii),kind=8)
					endif  ! .not.allsameparity
					if(doHam)then
						htraceall = htraceall + real(H_J_MK(sdX,sdX,iJ)%MK(ii,ii),kind=8)				
					end if
					if(.not.doHam .and. numSD >1)then
						do sdY = sdX,numSD
							overlapmat(sdX,sdY)= overlapmat(sdX,sdY)+N_J_MK(sdX,sdY,iJ)%MK(ii,ii)					
!					if(sdY> sdX)overlapmat(sdY,sdY)= overlapmat(sdY,sdX)+dconjg(N_J_MK(sdX,sdY,iJ)%MK(ii,ii)	)				
						end do
					end if
				end do ! ii 
			end do ! sdX
			normtrace(iJ)=trace
			if(doHam)hamtrace(iJ) = htrace
			if(.not.allsameparity)pnormtrace(iJ)=ptrace
!		print*,' f_J: ',xJ,trace
			traceall = traceall+trace+ptrace
!		if(doHam)print*,xJ,trace,htrace
		end do ! iJ 
	
		IF (myMPIrank == root) THEN
			if(doHam)print*,' Trace Hamiltonian = ',htraceall
			print*,' sum rule ',traceall
		END IF ! myMPIrank == root
	
		if(.not.doHam .and. numsd > 1)then  ! DIAGONALIZE overlapmat TO SEE IF THERE ARE ANY SDs WHICH OBVIOUSLY OUGHT TO BE THROWN OUT
			do sdX = 1,numsd-1
				do sdY = sdX+1,numsd
					overlapmat(sdY,sdX)= conjg(overlapmat(sdX,sdY))
				end do
			end do
!			print*,overlapmat
			lwork = 2*numsd - 1
			ALLOCATE(work(lwork),overeig(numsd))
			ALLOCATE(rwork(3*numsd - 2))
!		print*,overlapmat(1,1)
			CALL zheev('V','U',numsd,overlapmat,numsd,overeig,work,lwork,rwork,info)
!		print*,info
			IF (myMPIrank == root) THEN
				print*,' SD overlap eigenvalues = ',overeig
			END IF ! myMPIrank == root
!		print*,overlapmat
		end if ! .not.doHam and numsd > 1

		!CALL MPI_BARRIER(icomm,ierr)

		return
	end subroutine tracemaster

!====================================================================
!
!  revised CWJ June 2015;
!  also prints out probabilities
!  revsied again CWJ April 2018:
!  more logical output; sorts on energies
! 
!  INPUT:
!     np      :  master declared dimension of arrays
!     nstates : = # of states
!     jval(:) : array of jvalues
!     eval(:) : array of energies
!     pval(:) : array of parities (+,-1)
!     obsvals(:) :: array of observables (if compute_expect)
!     parityflag : logical flag if true then possible to have different parities
!     nprint: how many to print to screen
!

SUBROUTINE J_WriteResults(tol,np,nf,jvals,pevals,obsvals,problist,hamlist,parityflag,nprint)

!use errortests
	use eigenpackage
	use lamplight
	implicit none

	real,    intent(in) :: tol
	integer (kind=8), intent(in) :: np, nf
	real    (kind=8), intent(in) :: jvals(np), pevals(2,np),obsvals(2,np)
	real    (kind=8), intent(in) :: problist(2,numOfJ),hamlist(2,numOfJ)
	logical, intent(in) :: parityflag
	integer(4),intent(in) :: nprint
	integer iprint,icount
	real(kind=8) :: probsum,hamSum,ex
	INTEGER :: i,iostatus
	CHARACTER (LEN = 1) :: choice
	CHARACTER (LEN = 4) :: shell
	CHARACTER (LEN = 26) :: name
	CHARACTER (LEN = 100) :: filename!, path

!path = '/home/Jtstaker/Desktop/JohnsonResearch/ProjectedData/'//TRIM(shell)//'_data/'

	call sortresults

	write(6,*)' '
	iprint = 0
	if(parityflag)then
		if(compute_expect)then
			write(6,*)' State    E     Ex      J     parity      <obs>'
		else
			write(6,*)' State    E     Ex      J     parity'
		end if ! compute_expect
		write(6,*)' ----------------------------------------------'
	
		do i = 1,nlevels !min(nlevels,nprint)
			ex= energylevels(i)-energylevels(1)
			if(jtarget < 0 .or. jlevels(i)==jtarget)then
				if(compute_expect)then
					write(6,3001)i,energylevels(i),ex,jlevels(i),paritylevels(i),obslevels(i)
				else ! compute_expect
					write(6,3001)i,energylevels(i),ex,jlevels(i),paritylevels(i)
				endif ! compute_expect
				iprint= iprint+1
				if(iprint ==nprint)exit
			end if ! jtarget
		end do
	else ! parityFlag
		if(compute_expect)then
			write(6,*)' State    E     Ex      J    <obs>  '
		else ! compute_expect
			write(6,*)' State    E     Ex      J     '
		end if ! compute_expect
		write(6,*)' ----------------------------------------'
		do i = 1,nlevels !min(nlevels,nprint)
			ex= energylevels(i)-energylevels(1)
			if(jtarget < 0 .or. jlevels(i)==jtarget)then
				if(compute_expect)then
					write(6,3002)i,energylevels(i),ex,jlevels(i),obslevels(i)
				else
					write(6,3002)i,energylevels(i),ex,jlevels(i)
				end if
				iprint= iprint+1
				if(iprint ==nprint)exit
			end if ! jtarget
		end do ! i 
	end if	! parityflag


1000 FORMAT(I3,4(G15.8))
	!PRINT*, ' Makes it here in LAMPoutput.f90 (~355)' ! TESTING, REMOVE
	WRITE(*,*) 'Write output to file? (Y or N)'
	!PRINT*, ' Makes it here in LAMPoutput.f90 (~357)' ! TESTING, REMOVE

	DO
		READ(*,*) choice
		IF ((choice == 'y').OR.(choice == 'Y')) THEN
			EXIT
		ELSEIF ((choice == 'n').OR.(choice == 'N')) THEN
			RETURN
		ELSE
			WRITE(*,*) 'Y or N please.'
		ENDIF
	END DO

	DO 
		WRITE(*,*) 'Enter file name (without extention -- .res added): '
		READ(*,*) name
		filename = TRIM(name)//'.res'
		OPEN(UNIT = 10, FILE = filename, STATUS = 'NEW',IOSTAT = iostatus)
		IF (iostatus > 0) THEN
			WRITE(*,*) "File already exists: Overwrite (o), append (a), new file name (n)?"
			READ(*,*) choice
			DO icount = 1,5
				IF ((choice == 'o').OR.(choice == 'O')) THEN
					CLOSE(UNIT = 10)
					OPEN(UNIT = 10, FILE = filename, STATUS = 'REPLACE')
					EXIT
				ELSEIF ((choice == 'a').OR.(choice == 'A')) THEN
					CLOSE(UNIT = 10)
					OPEN(UNIT = 10,FILE = filename, STATUS = 'OLD', POSITION = 'APPEND')
					EXIT
				ELSEIF ((choice == 'n').OR.(choice == 'N')) THEN
					CLOSE(UNIT = 10)
					EXIT
				ELSE
					WRITE(*,*) 'Incorrect choice.  Please select again.'
				END IF ! choice flag 
			END DO ! icount
			IF ((choice.NE.'n').AND.(choice.NE.'N')) EXIT
		ELSE
			EXIT
		END IF ! iostatus 
	END DO

	call sortresults

	write(10,*)' LAMP Results '
	write(10,*)' (cutoff criterion = ',tol,')'
	write(10,*)' '
	if(parityflag)then
		if(compute_expect)then
			write(10,*)' State    E     Ex      J     parity    < obs >'
		else
			write(10,*)' State    E     Ex      J     parity'
		end if
		write(10,*)' ------------------------------------------------'
!	DO i = 1, nf
!		WRITE(10,1001) i,pevals(1,i),pevals(2,i),jvals(i)		
!	END DO
1001 format(i3,2x,2F13.5,2x,f4.1)		
	
		do i = 1,nlevels
			ex= energylevels(i)-energylevels(1)
			if(compute_expect)then
				write(10,3001)i,energylevels(i),ex,jlevels(i),paritylevels(i),obslevels(i)
			else
				write(10,3001)i,energylevels(i),ex,jlevels(i),paritylevels(i)			
			end if
		end do
3001 format(i3,2x,2f13.5,2x,f4.1,2x,a1,2x,f13.5)

	else ! parity flag 
		if(compute_expect)then
			write(10,*)' State       E        Ex         J       <obs >  '
		else
			write(10,*)' State       E        Ex         J     '
		end if 
		write(10,*)' ------------------------------------'
!	DO i = 1, nf
!		WRITE(10,1002) i,pevals(1,i),jvals(i)		
!	END DO
1002 format(i3,2x,F12.5,2x,f4.1)	
		do i = 1,nlevels
			ex= energylevels(i)-energylevels(1)
			if(compute_expect)then
				write(10,3002)i,energylevels(i),ex,jlevels(i),obslevels(i)
			else
				write(10,3002)i,energylevels(i),ex,jlevels(i)
			end if
		end do

3002 format(i3,2x,2f13.5,2x,f4.1,2x,f13.5)
	end if	! parityflag 
!DO i = 1, nf
!	WRITE(10,1000) i,pevals(1,i),pevals(2,i),jvals(i)
!END DO
	probsum = 0.d0 
	write(10,*)' '
	write(10,*)' Fraction in original HF state: Norm'
	write(10,*)' '
	if(parityflag)then
		write(10,*)' J    frac(+)   frac(-)'
		write(10,*)'-----------------------'
		do i = 1,numOfJ
			write(10,2001)xjlist(i),problist(1,i),problist(2,i)
			probsum = probsum+problist(1,i)+problist(2,i)
		end do
2001 format(f4.1,2f10.6)
	else
		write(10,*)' J    frac'
		write(10,*)'-----------------------'
		do i = 1,numOfJ
			write(10,2001)xjlist(i),problist(1,i)
			probsum = probsum+problist(1,i)
		end do
	end if ! parityflag 
	write(10,*)' Total of HF state = ',probsum

	hamsum = 0.d0
	write(10,*)' '
	write(10,*)' Fraction in original HF state: Hamiltonian'
	write(10,*)' '
	if(parityflag)then
		write(10,*)' J    frac(+)   frac(-)'
		write(10,*)'-----------------------'
		do i = 1,numOfJ
			write(10,2002)xjlist(i),hamlist(1,i),hamlist(2,i)
			hamsum = hamsum+hamlist(1,i)+hamlist(2,i)
		end do
2002 format(f4.1,2(1X,f10.6))
	else
		write(10,*)' J    frac'
		write(10,*)'-----------------------'
		do i = 1,numOfJ
			write(10,2002)xjlist(i),hamlist(1,i)
			hamsum = hamsum+hamlist(1,i)
		end do
	end if
	write(10,*)' Total of sum of trace(H) = ',hamSum

	CLOSE(UNIT = 10)

	WRITE(*,*) ''
	WRITE(*,*) 'Data written to:',filename
	return

	
END SUBROUTINE J_WriteResults
!===================================================================
!---- ADDED IN 1.3.6-----------
!
subroutine sortresults
	use eigenpackage
	implicit none
	integer :: i,j,k
	character ::parswap
	real(8) :: eswap,jswap,etest,oswap
	
	do i = 1,nlevels -1
		etest = energylevels(i)

		k = i
		do j = i+1,nlevels
			if(energylevels(j)< etest)then
				k = j
				etest = energylevels(j)
			end if
		
		end do
		if(k > i)then ! swap
			eswap = energylevels(i)
			jswap = jlevels(i)
			parswap=paritylevels(i)
			oswap = obslevels(i)
			energylevels(i)=energylevels(k)
			jlevels(i) = jlevels(k)
			paritylevels(i)=paritylevels(k)
			obslevels(i)=obslevels(k)
			energylevels(k)=eswap
			jlevels(k)     =jswap
			paritylevels(k)=parswap
			obslevels(k)=oswap
			
			
		end if
	
	end do
	
	return
	
end subroutine sortresults

end module outerputer
