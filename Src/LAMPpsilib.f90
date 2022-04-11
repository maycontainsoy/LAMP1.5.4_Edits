module psis
!
!  wave functions--slater determinants
!
! subsumes old module system_parameters 
	USE nodeinfo
  implicit none
	  
  integer numprot,numneut,numsd ! from old module system_parameters

  logical special_override_p, special_override_n
	integer :: p_override, n_override
  type sdmaster
  complex(kind=8),pointer :: psd(:,:),nsd(:,:)
  end type sdmaster

  type(sdmaster), pointer :: psir(:),psil(:)  ! of it = time
		  
	complex (kind=8), allocatable :: psd(:,:),nsd(:,:)
	complex (kind=8), allocatable :: psdf(:,:,:),nsdf(:,:,:)

contains

! SML QUESTION: has this subroutine been replaced? Doesn't appeared to be called and 
! I suspect it has been replaced.
  subroutine readinsd(nsps,n,ifile,psi,errflag)
!
!  subroutine to read in SD from file written by SHERPA
!  NB: probably will have other routines as well
!
!  INPUT:
!  nsps   = # of single-particle states
!  n      = # of particles
!  ifile  = logical number of the file
!  
!  OUTPUT
!   psi(i,j) = slater determinant (real)
!   errflag  = true if there was a problem in reading
!


    implicit none
    integer nsps,n
    real psi(nsps,n) ! slater determinant

    integer i,j
    integer ifile
    logical errflag

    errflag=.false.

    do i = 1,n
      read(ifile,err=103,end=103)(psi(j,i),j=1,nsps)
    enddo
    return
103 continue
    errflag=.true.
    return

  end subroutine readinsd  

!=====================================================================
! MPI added by SML in Jan 2022
  subroutine set_nuclide
		implicit none
		INTEGER :: ierr ! for MPI communication

		special_override_p =.false.
		special_override_n =.false.
		
    IF (myMPIrank == root) THEN
    	PRINT*,' Enter Z, N '
	    READ*,numprot,numneut			
		END IF 

		! Update other ranks
		!CALL MPI_BARRIER(icomm,ierr) ! REMOVE
		CALL MPI_BCAST(numprot,1,MPI_INT,root,icomm,ierr)
		CALL MPI_BCAST(numneut,1,MPI_INT,root,icomm,ierr)
		CALL MPI_BARRIER(icomm,ierr) ! REMOVE?
		
		return
	end subroutine set_nuclide      	

!=====================================================================

  subroutine GetSherpaSD(pSD,nSD)
!
!  master routine for retrieving SD written out by SHERPA
!  must also convert them from real*4 to complex*8
!
!      use system_parameters
    use spstate
    implicit none

    integer ifile
    complex(kind = 8) ::pSD(nsps,numprot),nSD(nsps,numneut)
    real(kind=4),allocatable :: sdtmp(:,:)

    ifile = 73

    call OpenSherpaSD(ifile,nsps,numprot,numneut)
	  
    if(allocated(sdtmp))deallocate(sdtmp)
    if(numprot > 0 .or. special_override_p)then
			if(special_override_p)then
	      allocate(sdtmp(nsps,p_override))
	      call ReadSherpa(ifile,nsps,p_override,sdtmp)  
			  ! DO NOT CONVERT
		  else ! special_override_p 
        allocate(sdtmp(nsps,numprot))
        call ReadSherpa(ifile,nsps,numprot,sdtmp)
        call ConvertSD(nsps,numprot,sdtmp,pSD)
	    end if ! special_override_p
      deallocate(sdtmp)
    endif ! numprot > 0 or special_override_p
      
		if(numneut > 0 .or. special_override_n)then
			if(special_override_n)then
	      allocate(sdtmp(nsps,n_override))
	      call ReadSherpa(ifile,nsps,n_override,sdtmp)  
			  ! DO NOT CONVERT
		  else ! special_override_n
        allocate(sdtmp(nsps,numneut))
        call ReadSherpa(ifile,nsps,numneut,sdtmp)
        call ConvertSD(nsps,numneut,sdtmp,nSD)
	    end if ! special_override_n
      deallocate(sdtmp)
    endif ! numneut > 0 or special_override_n
    close(ifile)

    return
  end subroutine GetSherpaSD

!=====================================================================      

  subroutine OpenSherpaSD(ifile,nsps,Z,N)
!	use system_parameters

    use sporbit
    implicit none

!.......... 
    integer N,Z
    integer nsps
    integer ifile

!..............FILE HANDLING..........................
    character ychar*1
    character filename*25  		! 
    integer ilast
    integer tempfile			! location of temporary file
    data tempfile/99/
    logical errflag
    character title*60
    logical openflag,failflag

!-------------- DUMMIES 
    integer i,ii,j
    integer norb
  	integer zz,nn
		INTEGER :: ierr ! for MPI communication

!-------------- OPEN FILE ---------------
1  continue

    openflag = .false.
		!CALL MPI_BARRIER(icomm,ierr) ! REMOVE?

		!IF ( myMPIrank == root ) THEN 
    	do while(.not.openflag)
	      write(6,*)' Enter input filename (.sd)'
  	    read(5,'(a)')filename
    	  ilast=index(filename,' ')
    	  if(ilast.ne.0)then
    	    ilast=ilast-1
    	  else
    	    ilast=15
    	  endif
    	  open(unit=ifile,file=filename(1:ilast)//'.sd',status='old', err=33, form='unformatted')
    	  openflag = .true.
33  continue
        IF (.not.openflag) THEN
          write(6,*)' That file does not exist; do you wish to try another file (y/n)?'
          read(5,'(a)')ychar
          if(ychar == 'N' .or. ychar == 'n')stop ! might need to add MPI abort here instead
				END IF 
    	END DO   ! while on openflag
		!END IF 

		!CALL MPI_BARRIER(icomm,ierr) ! REMOVE?
!..............READ IN HEADER INFO...............
!              CHECK THAT MATCHES SYSTEMS
! SML NOTES: what is norb in the SD file? It appears to be read once to 
! check the value against numorb from the .sps file but is read in again 
! and appears to be dependent on i (1,numorb)?
		! Seems to be consistently making it here
		!PRINT*, ' node # = ', myMPIrank, ' test allocateSD here 1 (202)' ! REMOVE
		!IF ( myMPIrank == root ) THEN 
	    read(ifile)norb
	    failflag = .false.
  	  if(norb /= numorb)then
				write(6,*)' # of orbits mismatch ',norb,numorb
      	failflag = .true.
    	endif
		!END IF ! myMPIrank == root 

		! This is the bcast causing the error
		!CALL MPI_BCAST(norb,1,MPI_INT,root,icomm,ierr)

		! Seems to be consistently making it here 
		!PRINT*, ' node # = ', myMPIrank, ' test allocateSD here 2 (211) norb = ', norb
    
		!IF ( myMPIrank == root ) THEN 
			do i = 1,numorb
				read(ifile)ii,norb,j
				if(norb /= orbqn(i)%nr .or. j /= orbqn(i)%j)then
	  			write(6,*)' mismatch n,j:',norb,orbqn(i)%nr,j,orbqn(i)%j
	  			failflag = .true.
      	endif ! norb
    	enddo ! i 
		!END IF ! myMPIrank == root, moved further down

    	if(failflag)then
    	  write(6,*)' The single particle space does not match ', &
     							' with that previously chosen.'
      	write(6,*)' Exit (x) or choose another file (c)?'
      	read(5,'(a)')ychar
      	if(ychar.eq.'x' .or. ychar.eq.'X')stop
      	close(ifile)
      	goto 1
    	endif
		!END IF ! myMPIrank == root 

		!CALL MPI_BARRIER(icomm,ierr) ! REMOVE?
		!PRINT*, ' node #', myMPIrank, ' makes it here 2'

!...............CHECK IF N,Z match........
! SML NOTES: it seems reasonable to restrict the majority of this 
! subroutine to root 
		!IF ( myMPIrank == root ) THEN
	    read(ifile)zz,nn
!      if(n.ne.nn .or. z.ne.zz)then
		  
!.......... ADDED in 1.5.1..... special override
  	  if( n /= nn)then
			  if(special_override_n)then
				  print*,' not reading in neutron SD '
				  n_override = nn
			  elseif(n==0)then
				  print*,' Do you mean to ignore the neutron SD?'
				  read(5,'(a)')ychar
				  if(ychar=='y'.or. ychar=='Y')then
					  special_override_n=.true.
					  n_override = nn
				  end if
			  else
		   	  write(6,*)' mismatch N: Old: ',n,', new: ',nn
		      write(6,*)' Exit (x) or choose another file (c)?'
		      read(5,'(a)')ychar
		      if(ychar.eq.'x' .or. ychar.eq.'X')stop
		      close(ifile)
		      goto 1
			  end if ! special_override_n
		  end if ! n /= nn
    
			if( z/= zz )then
			  if(special_override_p)then
				  print*,' not reading in proton SD '
				  p_override = zz
			  elseif(z==0)then ! special_override_p
				  print*,' Do you mean to ignore the proton SD?'
				  read(5,'(a)')ychar
			  
					if(ychar=='y'.or. ychar=='Y')then
					  special_override_p=.true.
					  p_override = zz
				  end if
			  
					if(special_override_p .and. special_override_n)then
					  print*,'cannot have an empty nucleus!'
					  stop
				  end if
			  else ! special_override_p
	  	 	  write(6,*)' mismatch Z: Old: ',z,', new: ',zz
	  	    write(6,*)' Exit (x) or choose another file (c)?'
	  	    read(5,'(a)')ychar
	  	    if(ychar.eq.'x' .or. ychar.eq.'X')stop
	  	    close(ifile)
	  	    goto 1
			  end if ! special_override_p
	  	end if
		

			read(ifile)title
    	write(6,'(60a)')'Title card: "',title,'"'
		!END IF ! myMPIrank == root 

		!CALL MPI_BARRIER(icomm,ierr)
		!PRINT*, ' node #', myMPIrank, ' makes it here 3 (OpenSherpaSD end)'
    return
  end subroutine OpenSherpaSD

!=====================================================================
  subroutine ReadSherpa(ifile,nsps,np,sd)

!....... NOTE THIS ASSUMES real*4 
    implicit none

    integer nsps,np
    real(kind=4),intent(OUT) :: sd(nsps,np) ! slater determinant

    integer i,j
    integer ifile
    logical errflag

    errflag=.false.

    do i = 1,np
      read(ifile,err=103,end=103)(sd(j,i),j=1,nsps)
    enddo
    return
103 continue
    errflag=.true.
    return

    return
  end subroutine ReadSherpa

!=====================================================================
! Converts temporary proton and neutron Slater determinants from real4 
! (from Sherpa) to double complex.
  subroutine ConvertSD(nsps,np,sd,zsd)

    implicit none
    integer nsps,np

    real(kind=4) :: sd(nsps,np)
    complex(kind = 8) :: zsd(nsps,np)

    integer i,j

    do i = 1,nsps
      do j = 1,np
        zsd(i,j) = dcmplx( sd(i,j),0.0)
      enddo	! loop over j
		enddo	! loop over i

    return
  end subroutine ConvertSD

!=====================================================================
  subroutine read_sd_txt(psdtrial,nsdtrial)
!=====================================================================
! Reads a (time-reversal) Slater determinant
! in text format
!====================================================================
!      use system_parameters
  	use spstate

  	implicit none

  	integer  :: i,j,k,n,indx,ilast  
  	complex(kind=8) :: psdtrial(nsps,numprot),nsdtrial(nsps,numneut)
  	real vec(nsps)
  	character :: ychar
  	character*25 :: filename

311 continue
  	write(6,*)' Enter input filename (.tsd)'
  	read(5,'(a)')filename
  	ilast=index(filename,' ')
  	if(ilast.ne.0)then
  	  ilast=ilast-1
  	else
  	  ilast=15
  	endif

  	open(unit=10,file=filename(1:ilast)//'.tsd',status='old',err=33,form='formatted')
  	goto 44
33  continue
  	write(6,*)' That file does not exist; do you wish to try another file (y/n)?'
    read(5,'(a)')ychar
    if(ychar.eq.'n' .or.ychar.eq.'N')then
      return
    else
      goto 311
    endif
44  continue

    psdtrial(:,:)=dcmplx(0.d0,0.d0)
    if(numneut/=0)nsdtrial(:,:)=dcmplx(0.d0,0.d0)

    do n=1,numprot
      read(10,*)(vec(i),i=1,nsps)
      do i=1,nsps
        psdtrial(i,n)=dcmplx(dble(vec(i)),0.d0)
      enddo
    enddo
    
		do n=1,numneut
      read(10,*)(vec(i),i=1,nsps)
      do i=1,nsps
        nsdtrial(i,n)=dcmplx(dble(vec(i)),0.d0)
      enddo
    enddo
122 continue
    close(10)

    return
  end subroutine read_sd_txt
!
!==============================================================================
!
! allocate the slater determinants psd / nsd & psdf / nsdf
!
! OUTPUT:
!   psd, nsd = proton/neutron slater determinants
!   psdf, nsdf = ^^ with a third dimension for multiple slater deterimants (i.e. numsd > 1)
!   ParityTest = checks for proper parity [?]
!
! SUBROUTINES CALLED:
!   PairLog = checks for parity and flags PairLog appropriately
!
subroutine allocateSlaterDet 
	use spstate
	implicit none

!............INTERNAL....................
  integer :: ii,jj,iii
	integer :: ifile
	logical :: errflag
	integer :: nSDs,isd,jsd,ksd
!	real(kind=4),allocatable :: psdtmp(nsps,numprot),nsdtmp(nsps,numneut)

	real(kind=4),allocatable :: psdtmp(:,:),nsdtmp(:,:)
  complex(kind=8) :: pSDx(nsps,numprot),nSDx(nsps,numneut) 
  !COMPLEX(KIND=8), ALLOCATABLE :: pvec(:), nvec(:) ! For MPI communication
	COMPLEX(KIND=8), ALLOCATABLE :: pvec(:), nvec(:)
	real :: eee

	INTEGER :: ierr ! for MPI communication
	INTEGER :: i, j, k, a ! integers for broadcasting SDs
	REAL :: b, c
	
	ifile  = 73

!----------------- ALLOCATE SLATER DETERMINANTS -------------------
  allocate(psd(nsps,numprot),nsd(nsps,numneut))
	allocate(psdtmp(nsps,numprot))
	allocate(nsdtmp(nsps,numneut))
	
  IF (myMPIrank == root) THEN 
		print*,' Enter number of Slater determinants '
		read*,numsd
	END IF ! myMPIrank 

	CALL MPI_BARRIER(icomm,ierr)
	CALL MPI_BCAST(numsd,1,MPI_INT,root,icomm,ierr)
	!CALL MPI_BARRIER(icomm,ierr)

!	print*,' testing ',numsd,nsps,numprot,numneut, ' TESTING'
	! Check that these are correctly allocated on all ranks
  ALLOCATE (psdf(numsd,nsps,numprot),nsdf(numsd,nsps,numneut))
	ALLOCATE (pvec(nsps), nvec(nsps))

	jj=0

	! ~*~*~ Notes to SML from SML 
	! special override is not commonly used, but should be tested with MPI
	! I don't think this is the source of the error but it is something that 
	! we may need to implement in the future.
	! ADD ROOT PROTECTION
	CALL MPI_BARRIER(icomm,ierr) ! TESTING, REMOVE?
	IF (myMPIrank == root) THEN 
  	do ii = 1,numsd
			if(jj==numsd)return

    	call OpenSherpaSD(ifile,nsps,numprot,numneut) ! located in this file

			! NOT CURRENTLY USED 
		  if(special_override_p)then
				deallocate(psdtmp)
		   	allocate(psdtmp(nsps,p_override))
		  end if	
	  
			! NOT CURRENTLY USED 
			if(special_override_n)then
				deallocate(nsdtmp)
		   	allocate(nsdtmp(nsps,n_override))
		  end if

!... special override is for separating proton and neutron SDs....
			! NOT CURRENTLY USED 
    	if(special_override_p)call count_SDs_in_file(ifile,nsps,p_override,numneut,nSDs)
    	if(special_override_n)call count_SDs_in_file(ifile,nsps,numprot,n_override,nSDs)

!........... 'NORMAL' CASE.....................	   
	  	if(.not.special_override_p .and. .not.special_override_n)then
			  call count_SDs_in_file(ifile,nsps,numprot,numneut,nSDs) ! located in this file
	  	end if 

	  	rewind(ifile)
	  	call read_SD_header(ifile,.false.,errflag) ! located in this file
	   
	  	if(nSDs==1)then
	  	  if(numprot > 0)then			   
	    	  call ReadSherpa(ifile,nsps,numprot,psdtmp) 	! located in this file, reads into array from file j = 1, nsps and i = 1, numProt
	    	  call ConvertSD(nsps,numprot,psdtmp,pSDx) 		! located in this file, coverts psdtmp to DCMPLX type
				end if ! numprot > 0 
		  
				if(special_override_p)then			   
	  	    call ReadSherpa(ifile,nsps,p_override,psdtmp)
			  end if
			   
	  	  if(numneut > 0)then
	  	    call ReadSherpa(ifile,nsps,numneut,nsdtmp)	! located in this file, reads into array from file j = 1, nsps and i = 1, numProt
	  	    call ConvertSD(nsps,numneut,nsdtmp,nSDx)		! located in this file, coverts dsdtmp to DCMPLX type
			  end if ! numneut > 0 
		  
				if(special_override_n)then
	  	    call ReadSherpa(ifile,nsps,n_override,nsdtmp)
			  end if
		  
				! Adds columns to total proton and neutron Slater determinant arrays
				! Need to figure out where psdf and nsdf are allocated, I think these will need 
				! to be broadcast to the other ranks? Or if they are allocated on each rank
				! I can broadcast the psdx and nsdz arrays as columns onto other ranks?
				jj = jj + 1
      	psdf(jj,:,:) = psdx
      	nsdf(jj,:,:) = nsdx
		  	close(ifile)
		   
	  	else ! nSDs == 1
		  	print*,' There are ',nSDs,' Slater determinants in this file '
		  	print*,' Enter a Slater determinant index (0 or -1 to stop )'
		  	print*,' To read from an ordered list, enter -999 '
				
		  	read*,isd
		  	if(isd==-999)then
				  open(unit=33,file="orderedSDlist.dat",status='old')
				  do iii=1,nSDs
					  read(33,*,end=9909)ksd,eee
					  print*,ksd,eee
					  rewind(ifile)
					  call read_SD_header(ifile,.false.,errflag)
					  do jsd = 1,ksd
					    if(numprot > 0)then
					      call ReadSherpa(ifile,nsps,numprot,psdtmp)
					      call ConvertSD(nsps,numprot,psdtmp,pSDx)
						  end if
					  
							! NOT CURRENTLY USED 
							if(special_override_p)then
					      call ReadSherpa(ifile,nsps,p_override,psdtmp)
						  end if
				    
							if(numneut > 0)then
					      call ReadSherpa(ifile,nsps,numneut,nsdtmp)
					      call ConvertSD(nsps,numneut,nsdtmp,nSDx)
						  end if
					  
							! NOT CURRENTLY USED 
							if(special_override_n)then
					      call ReadSherpa(ifile,nsps,n_override,nsdtmp)
						  end if
					  end do ! jsd
					  jj = jj + 1
		    	  if (numprot>0) psdf(jj,:,:) = psdx
		    	  if (numneut>0) nsdf(jj,:,:) = nsdx
			   
				  	if ( jj == numsd ) then
						  close(ifile)		   
						  return
				  	end if 
				  end do ! iii
9909  continue
		  	end if ! isd = -999
		   
		  	do while(isd > 0 .and. isd <= nSDs)
				  rewind(ifile)
				  call read_SD_header(ifile,.false.,errflag)
				  do jsd = 1,isd
				    if(numprot > 0)then
				      call ReadSherpa(ifile,nsps,numprot,psdtmp)
				      call ConvertSD(nsps,numprot,psdtmp,pSDx) 
					  end if ! numprot > 0 
				  
						! NOT CURRENTLY USED 
						if(special_override_p)then
			  	    call ReadSherpa(ifile,nsps,p_override,psdtmp)
					  end if
			    
						if(numneut > 0)then
			  	    call ReadSherpa(ifile,nsps,numneut,nsdtmp)
			  	    call ConvertSD(nsps,numneut,nsdtmp,nSDx)
					  end if
				  
						! NOT CURRENTLY USED 
						if ( special_override_n ) then
			  	  	call ReadSherpa(ifile,nsps,n_override,nsdtmp)
					  end if
				  end do ! jsd 
				  jj = jj+1
	    	  if (numprot>0) psdf(jj,:,:) = psdx
	    	  if (numneut>0) nsdf(jj,:,:) = nsdx
			   
				  if ( jj == numsd ) then
					  close(ifile)		   
					  !return
					  goto 201 ! replaces old return statement (above)
				  end if
				  print*,' Enter a Slater determinant index (0 or -1 to stop )'
				  read*,isd
			  end do ! isd > 0 and <= nSDs
	  	end if ! nSDs == 1
  	end do ! ii 
		201 continue 
	END IF ! myMPIrank == root 
	! End MPI root protection 

	!
	!
	!
	! The location of the zeroing out of pvec and nvec has been tested 
	! at various points in this subroutine.
	pvec(:) = (0.d0,0.d0)
	nvec(:) = (0.d0,0.d0)
	!
	!
	!
	! I was attemping to see if 'i' was causing the issue and if different 
	! ranks were getting incorrect information; this came from the discovery 
	! of the data 'offset'.
	i = 0

	! This portion of the code is functioning but can probably be improved
	! and optimized. This should probably been tested for speed up in two or 
	! three different ways to ensure efficiency. 
	CALL MPI_BARRIER(icomm,ierr) ! TESTING
	DO isd = 1, numsd 
		PRINT*, ' Node = ', myMPIrank, ' isd = ', isd
		DO a = 1, numprot 
			IF (myMPIrank == root) THEN
				DO i = 1, nsps 
					pvec(i) = psdf(isd,i,a)
				END DO ! i 
			END IF ! myMPIrank 
			CALL MPI_BARRIER(icomm,ierr)
			CALL MPI_BCAST(pvec(:),nsps,MPI_COMPLEX8,root,icomm,ierr)

			IF (myMPIrank /= root) THEN
				DO i = 1, nsps 
					psdf(isd,i,a) = pvec(i)
				END DO 
			END IF ! myMPIrank /= root
		END DO ! a
	END DO ! isd 

	DO isd = 1, numsd 
		DO a = 1, numneut 
			IF (myMPIrank == root) THEN 
				DO i = 1, nsps 
					nvec(i) = nsdf(isd,i,a)
				END DO ! i 
			END IF ! myMPIrank 
			CALL MPI_BARRIER(icomm,ierr)
			CALL MPI_BCAST(nvec(:),nsps,MPI_COMPLEX8,root,icomm,ierr)

			IF (myMPIrank /= root) THEN 
				DO i = 1, nsps 
					nsdf(isd,i,a) = nvec(i)
				END DO ! i 
			END IF ! myMPIrank /= root 
		END DO ! a 
	END DO ! isd

	!CALL MPI_BARRIER(icomm,ierr)

	! ~*~*~ Notes from SML to SML
	! Need to test this subroutine to make sure that the information is being 
	! broadcast correctly and that the code runs for different model spaces, different 
	! number of Slater determinants and reading information in from an ordered list. 
	
	! Need to check the type of data for the Slater determinants before broadcast. 
	! I think it comes in as real valued from Sherpa but gets converted to a double
	! complex in the convertSD subroutine (in this file). 
	RETURN
END SUBROUTINE allocateSlaterDet
!======================================================
subroutine read_SD_header(ifile,verbose,failflag)
	
	use sporbit
!	use system_parameters
	
	implicit none
	integer,intent(IN) :: ifile
	logical,intent(IN) :: verbose
	logical,intent(OUT) :: failflag
	
	integer :: norb,i,ii,j
	integer :: nn,zz
  character title*60
	
	failflag = .false.
	
  read(ifile)norb
 
  if(norb /= numorb)then
    write(6,*)' # of orbits mismatch ',norb,numorb
    failflag = .true.
  endif
  
	do i = 1,numorb
    read(ifile)ii,norb,j
    if(norb /= orbqn(i)%nr .or. j /= orbqn(i)%j)then
      write(6,*)' mismatch n,j:',norb,orbqn(i)%nr,j,orbqn(i)%j
      failflag = .true.
    endif
  enddo
    
	if(failflag)then
		close(ifile)
		return
	end if

!...............CHECK IF N,Z match........
  read(ifile)zz,nn
  if (((numneut.ne.nn) .and. .not. special_override_n) .or. & 
	   ((numprot.ne.zz)  .and. .not. special_override_p))then
	  write(6,*)' mismatch Z,N. Old: ',numprot,numneut, ', new: ',zz,nn
		failflag = .true.
		close(ifile)
		return
  endif

  read(ifile)title

  if(verbose)write(6,'(60a)')'Title card: "',title,'"'
		
	return
end subroutine read_SD_header

!===================================================================

subroutine count_SDs_in_file(ifile,nsps,z,n,nSDs)
	implicit none
	integer,intent(IN) :: ifile
	integer,intent(IN) :: nsps,z,n
	integer,intent(OUT) :: nSDs
	integer :: isd
	logical :: errflag

	INTEGER :: ierr ! for MPI communication
	
	real(kind=4) :: psdtmp(nsps,z),nsdtmp(nsps,n)
	
	nSDs = 0
	
	do isd = 1,100000
		errflag=.false.
   	if(z > 0)call ReadSherpa2(ifile,nsps,z,psdtmp,errflag)
		if(errflag)return
   	if(n > 0)call ReadSherpa2(ifile,nsps,n,nsdtmp,errflag)
		if(errflag)return
		nSDs = nSDs+1
	end do
	print*,' Wait should not have gotten here ',nSDs
	stop

end subroutine count_SDs_in_file
!===================================================================
subroutine ReadSherpa2(ifile,nsps,np,sd,errflag)

!....... NOTE THIS ASSUMES real*4 
!use spstate, only:spsqn
  implicit none

  integer nsps,np
  real(kind=4),intent(OUT) :: sd(nsps,np) ! slater determinant

  integer i,j
  integer ifile
  logical errflag

	INTEGER :: ierr ! for MPI communication

  errflag=.false.

  do i = 1,np
    read(ifile,err=103,end=103)(sd(j,i),j=1,nsps)
!		print*,' sd ',sd(:,1)
  enddo

  return
103 continue
  errflag=.true.
  return

  return
end subroutine ReadSherpa2
	
end module psis
