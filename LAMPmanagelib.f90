
module managers
	use nodeinfo
	implicit none
	
contains
!============================================================	
subroutine projectorator
	use lamplight
!	use sampler
!	use projections
	implicit none

	integer :: intM,intK
	real(kind=4) :: floatK,floatM,xJ
	integer :: iJ,ii,localdim
	real(kind=8) :: trace,ptrace,traceall

	REAL :: rate 
	INTEGER :: c1, c2, cr

	! Need to redo timing for MPI
	! CALL SYSTEM_CLOCK(count_rate=cr)
	! rate = REAL(cr)

	! Need to redo timing for MPI
	! CALL SYSTEM_CLOCK(c1)
	call set_euler_meshes		! in LAMPutils.f90
	! CALL SYSTEM_CLOCK(c2)
	! WRITE(*,*) '**************'
	! WRITE(*,*) 'Set euler meshes system clock: ', (c2-c1)/rate
	! WRITE(*,*) '**************'
	
!..... SET UP FINAL MATRICES.....................
	CALL allocate_delta														! in LAMPutils.f90, no MPI protection needed
	CALL allocate_arrays4sampling									! in LAMPutils.f90, no MPI protection but need to locate doHam
	! doHam might not be needed at this step, I believe this is just initialization
	
!.... sample over Euler angles, e.g.	< Psi | R(alpha_i, beta_j, gamma_k) | Psi> = N_ijk
! 	
	! Need to redo timing for MPI
	! CALL SYSTEM_CLOCK(c1)
	!PRINT*, ' Node = ', myMPIrank, ' before basic_sampling'
	call basic_sampling			! in LAMPsampler.f90
	!PRINT*, ' Node = ', myMPIrank, ' after basic_sampling'
	! CALL SYSTEM_CLOCK(c2)
	! WRITE(*,*) '**************'
	! WRITE(*,*) 'basic sampling system clock: ', (c2-c1)/rate
	! WRITE(*,*) '**************'
	
!
!...... solve for M,K from alpha and gamma => N_i,MK.....................	
	CALL populateZmat															! in LAMPutils.f90, no MPI protection needed
	CALL invert_alpha_and_gamma										! in LAMPutils.f90, ok for now?
																								! may become an issue as numSD >> 1

	CALL allocate_tilded_matrices									! in LAMPutils.f90, no MPI protection needed
	CALL allocate_X_J_MK													! in LAMPutils.f90, no MPI protection needed
	DO intM = 1, J2max1 !-Jmax to Jmax
!		    if(.not.chooseMval(intM))cycle
		DO intK = 1, J2max1 ! intM !-Jmax to M; M>=J
!			    if(.not.chooseMval(intK))cycle
			CALL compute_delta(intM,intK)							! in LAMPutils.f90, no MPI protection needed
			CALL compute_tilded_matrices(intM,intK)		! in LAMPutils.f90, ok for now? numSD issue
!................. NOW SOLVE FOR BETA -> J..............				
!                 obtain N^J_MK etc
			CALL solve4projected_matrices(intM,intK)	! in LAMPutils.f90, ok for now? numSD issues
		END DO ! intK
	END DO ! int M
	
	RETURN
	
end subroutine projectorator
!------------------------ MAIN DRIVER ----------
!
!  SUBROUTINES CALLED:
!    set_euler_meshes
!    allocate_delta
!    perturb_beta_mesh
!    testdeltas
!    allocate_arrays4sampling
!    Projection_with_Parity
!    populateZmat, invert_alpha_and_gamma, allocate_tilded_matrices, allocate_X_J_MK
!    compute_delta
!    solve4projected_matrices
!
subroutine projectorator_dens
	use lamplight
	use project1bdense
	use dens1body
	implicit none

	integer :: intM,intK
	integer :: eyeroll
	real(kind=4) :: floatK,floatM,xJ
	integer :: iJ,ii,localdim
	real(kind=8) :: trace,ptrace,traceall

	REAL :: rate 
	INTEGER :: c1, c2, cr

	CALL SYSTEM_CLOCK(count_rate=cr)
	rate = REAL(cr)

	WRITE(*,*) "*********"
	WRITE(*,*) "THIS IS EVEN WORKING?!"
	WRITE(*,*) "*********"
	
	call allocate_delta
	call allocate_arrays4sampling
	call setup_den1b_orig_array(numsd,J2max1**2*numOfJ)

	call density_sampling

    print*,' DENSITIES SAMPLED; inverting alpha, gamma '
!...... solve for M,K from alpha and gamma => N_i,MK.....................	
    call populateZmat    

	! Timing by SML
	CALL SYSTEM_CLOCK(c1)
    call invert_alpha_and_gamma4dense
	CALL SYSTEM_CLOCK(c2)
	WRITE(*,*) "*********"
	WRITE(*,*) "Invert alpha and gamma system clock: ", (c2-c1)/rate
	WRITE(*,*) "*********"

	call allocate_tilded_density
	call allocate_rho_J_MK

    print*,' inverting beta '

    do intM = 1, J2max1 !-Jmax to Jmax
           do intK = 1, J2max1 ! intM !-Jmax to M; M>=J
				
				call compute_tilded_density(intM,intK)
!................. NOW SOLVE FOR BETA -> J..............				
!                 obtain N^J_MK etc

				call solve4projected_density(intM,intK)				
						
			end do ! intK
	end do !int M
	print*,' projected '
	call setup_transformed_rhos

	! Timing by SML
	CALL SYSTEM_CLOCK(c1)
	call transform_rhos_step1
	CALL SYSTEM_CLOCK(c2)
	WRITE(*,*) "*********"
	WRITE(*,*) "Transform rhos step1 system clock: ", (c2-c1)/rate
	WRITE(*,*) "*********"

    print*,' all transformed '
	call writeout_dens
	
	return
end subroutine projectorator_dens
!====================================================
subroutine findNewJmax(Jtol,newJmax)
	use lamplight
	use tracy
	implicit none
	real(4) :: Jtol
	real(4) :: newJmax
	integer :: intJ,newnumOfJ
	INTEGER :: ierr ! For MPI communication
	
	newJmax = 0.0
	newnumOfJ=1
    do intJ = 1,numOfJ
		if(.not.allsameparity)then
			if(normtrace(intJ)+pnormtrace(intJ) > Jtol)then
				newJmax = xJlist(intJ)
				newnumOfJ  = intJ
			endif
		elseif(normtrace(intJ)> Jtol)then
				newJmax = xJlist(intJ)
				newnumOfJ  = intJ
		endif
    end do
	Jmax = newJmax
    J2max1 = int(2.*newJmax)+1
	numOfJ=newnumOfJ
	
	numOfJout = numOfJ
	J2max1out = J2max1   ! used for output; may be modified later
	
	IF (myMPIrank == root) THEN 
		print*,' '
		print*,' New Jmax = ',newJmax,numOfJ
!............ ADD ABILITY TO RESTRICT OUTPUT............	
		print*,' '
		print*,' NEW OPTION: '
		print*,' Enter a maximum J for output (enter -1 to output all )'
		print*,' (Large values of J, especially for a large number of reference states, '
		print*,'  can take a long time to diagonalize/write out)'
		read*,newJmax
	END IF ! 
	CALL MPI_BARRIER(icomm,ierr)
	CALL MPI_BCAST(newJmax,1,MPI_REAL,root,icomm,ierr)
	!PRINT*, 'Node = ', myMPIrank, ' newJmax = ', newJmax ! TESTING, REMOVE
	if(newJmax < 0)return
	
	if(newJmax >= Jmax)return
	
	J2max1out = int(2.*newJmax)+1
	
	do intJ = 1,numOfJ
		if(xJlist(intJ)==newJmax)then
			numOfJout = intJ
			exit
		end if
		if(xJlist(intJ) > newJmax)then
			numOfJout = intJ-1
			exit
		end if
		
	end do
	
	
	return
	
end subroutine findNewJmax

!==============================================================================
!
! user enters values for the Norm matrix tolerance and the number of points for inversion.
!
! OUTPUT:
!   tolerance = cuttoff point for eliminating zeroes in the Norm matrix (during zhegvdiag)
!
subroutine inputTolerance(tolerance)

    implicit none
!............OUTPUT......................
    real (kind=4), intent(out) :: tolerance
		INTEGER :: ierr ! for MPI communication

!---------------------- USER ENTERS VALUES ------------------------
		IF (myMPIrank == root) THEN 
    	print *, ' Enter tolerance for norm (typical = 0.0001 to 0.001) '
    	read(*,*) tolerance
		END IF ! myMPIrank == root 
		CALL MPI_BARRIER(icomm,ierr)
		CALL MPI_BCAST(tolerance,1,MPI_REAL,root,icomm,ierr)

return 
end subroutine inputTolerance

!==============================================================================	


!
!  ROUTINE TO WRITE GENERALIZED MATRICES TO FILE
!  added 1.4.4  July 2019 by CWJ @ SDSU
!
subroutine write_mat2file
!	use system_parameters
	use sporbit
	use spstate
	use lamplight
	use project1bdense,only:binary_dens
	implicit none
	
	character*80 ::filename
	integer :: ilast
	
	integer :: intJ,matdim
	real :: floatJ
	integer :: ii,jj,sdX, sdY
	
	print*,' '
	print*,' Enter name of .mat file (do not enter extension)'
	print*,' (will overwrite existing data if file already exists)'
	
	if(binary_dens)then
		print*,' (Writing as unformatted binary file )'
	end if
	read(5,'(a)')filename
	ilast = index(filename,' ')-1
	
	if(binary_dens)then
		open(unit=21,file=filename(1:ilast)//".mat",status='unknown',form='unformatted')
		
	else
		open(unit=2,file=filename(1:ilast)//".mat",status='unknown')
		
	end if
	
	if(binary_dens)then
		write(21)numorb,nsps
		write(21)numprot,numneut,numsd
		write(21)J2max1out,numOfJout,allsameparity
	else
		write(2,*)numorb,nsps
		write(2,*)numprot,numneut,numsd
		write(2,*)J2max1out,numOfJout,allsameparity
	end if

    do intJ = 1, numOfJout
		
        floatJ = xJlist(intJ)
		matdim = N_J_MK(1,1,intJ)%dimMK
		
		do sdX = 1,numsd
			do sdY = 1,numsd
				do ii = 1,matdim
					if(binary_dens)then
						write(21)(N_J_MK(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
						
					else
						write(2,*)(N_J_MK(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
					end if
				end do !ii
			end do !sdY
		end do ! sdX
		do sdX = 1,numsd
			do sdY = 1,numsd
				do ii = 1,matdim
					if(binary_dens)then
						write(21)(H_J_MK(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
					else
						write(2,*)(H_J_MK(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
					end if
				end do !ii
			end do !sdY
		end do ! sdX
		if(allsameparity)cycle
		do sdX = 1,numsd
			do sdY = 1,numsd
				do ii = 1,matdim
					if(binary_dens)then
						write(21)(PN_J_MK(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
					else
						write(2,*)(PN_J_MK(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
					end if
				end do !ii
			end do !sdY
		end do ! sdX
		do sdX = 1,numsd
			do sdY = 1,numsd
				do ii = 1,matdim
					if(binary_dens)then
						write(21)(PH_J_MK(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
					else
						write(2,*)(PH_J_MK(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
					end if
				end do !ii
			end do !sdY
		end do ! sdX		
		
	end do  ! int J

    if(binary_dens)then
		close(21)
	else
 	    close(2)
	end if
	return
end subroutine write_mat2file
!
!  ROUTINE TO WRITE GENERALIZED MATRICES TO FILE
!  added 1.4.4  July 2019 by CWJ @ SDSU
!
!  NOTE this version assumes some files have already been created
!
!  need to modify to allow for reading binary files
!
subroutine read_matfromfile
!	use system_parameters
	use sporbit
	use spstate
	use lamplight
	use project1bdense,only:binary_dens
	
	implicit none
	
	character*80 ::filename
	integer :: ilast
	
	integer :: intJ,matdim
	real :: floatJ
	integer :: ii,jj,sdX, sdY
	
	integer :: dum1,dum2,dum3
	logical :: ptest
	
	if(binary_dens)then
		print*,' '
		print*,' Warning, subroutine read_matfromfile needs to be modified'
		print*,'  to read unformatted files '
		stop
	end if
	
	print*,' '
	print*,' Enter name of .mat file (do not enter extension)'
	read(5,'(a)')filename
	ilast = index(filename,' ')-1
	open(unit=3,file=filename(1:ilast)//".mat",status='old')
	
	read(3,*)dum1,dum2
	if(dum1 /= numorb)then
		print*,' Error, mismatch in num of orbits ',numorb,dum1
		stop
	end if
	if(dum2 /= nsps)then
		print*,' Error, mismatch in num of s.p. states ',nsps,dum2
		stop
	end if
!	write(2,*)numorb,nsps
    read(3,*)dum1,dum2,dum3
	if(dum1 /= numprot)then
		print*,' Error, mismatch in num prot',numprot,dum1
		stop
	end if
	if(dum2 /= numneut)then
		print*,' Error, mismatch in num neutrons',numneut,dum2
		stop
	end if
	if(dum3 /= numsd)then
		print*,' Error, mismatch in num SDs',numsd,dum3
		stop
	end if
!	write(2,*)numprot,numneut,numsd
    read(3,*)dum1,numOfJObs,ptest
    if(ptest .neqv. allsameparity)then
		print*,' Error, mismatch in parity ',allsameparity,ptest
		stop
	end if
	if(numOfJObs > numOfJ)then
		print*,' Need to increase maximum J ',numOfJ,numOfJObs
		stop
	end if
!................... NOW ALLOCATE.....................


    allocate(Nobs_J_MK(numsd,numsd,numOfJobs))	
    allocate(Hobs_J_MK(numsd,numsd,numOfJobs))	
	if(.not.allsameparity) allocate(PNobs_J_MK(numsd,numsd,numOfJobs),PHobs_J_MK(numsd,numsd,numOfJobs))	

	do intJ = 1,numOfJobs
		
        floatJ = xJlist(intJ)
		matdim = N_J_MK(1,1,intJ)%dimMK
		do sdX =1,numsd
			do sdy = 1,numsd
				allocate(Nobs_J_MK(sdx,sdy,intJ)%MK(matdim,matdim))
				allocate(Hobs_J_MK(sdx,sdy,intJ)%MK(matdim,matdim))
				if(.not.allsameparity)allocate(PNobs_J_MK(sdx,sdy,intJ)%MK(matdim,matdim))
				if(.not.allsameparity)allocate(PHobs_J_MK(sdx,sdy,intJ)%MK(matdim,matdim))	
			end do
		end do
		
	end do
!	write(2,*)J2max1,numOfJ,ParityTest
    do intJ = 1, numOfJobs
		
        floatJ = xJlist(intJ)
		matdim = N_J_MK(1,1,intJ)%dimMK
		
		do sdX = 1,numsd
			do sdY = 1,numsd
				do ii = 1,matdim
					read(3,*)(Nobs_J_MK(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
				end do !ii
			end do !sdY
		end do ! sdX
		do sdX = 1,numsd
			do sdY = 1,numsd
				do ii = 1,matdim
					read(3,*)(Hobs_J_MK(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
				end do !ii
			end do !sdY
		end do ! sdX
		if(allsameparity)cycle
		do sdX = 1,numsd
			do sdY = 1,numsd
				do ii = 1,matdim
					read(3,*)(PNobs_J_MK(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
				end do !ii
			end do !sdY
		end do ! sdX
		do sdX = 1,numsd
			do sdY = 1,numsd
				do ii = 1,matdim
					read(3,*)(PHobs_J_MK(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
				end do !ii
			end do !sdY
		end do ! sdX		
		
	end do  ! int J

	close(3)
	return
end subroutine read_matfromfile


end module managers