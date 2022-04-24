!
!  Program REORDERINT
!
!  reorders BIGSTICK interaction files from one sps file to another
!  also capability to truncate and (optionally) include induced s.p.e. from core
!
!  by CWJ @ SDSU Dec 2015
!

module information
    implicit none
	integer :: norbiti,norbitf  ! number of orbits in initial and final spaces
	                            ! must have norbiti >= norbitf
!------------ CREATE A DEFINED TYPE-----------------------------------------
	type orb
		integer :: nr            ! radial quantum number
		integer :: j             ! 2 x j
		integer :: l             ! L
		integer :: par           ! parity = +/- 1
		integer :: map           ! maps to destination
		integer :: w             ! weighting; needed for adding in single particle energies
		logical :: core          ! if in "core"
    end type orb
	type (orb),allocatable,target :: orbqni(:), orbqnf(:) ! initial, final orbqn(iorb); assume single species

	real, allocatable :: spe(:)
	integer :: ntbmei,ntbmef  ! initial, final # of TBMEs
    integer, allocatable :: qlist(:,:)   ! labels for each TBME
	real, allocatable :: vtbme(:)
	logical :: inducedspes  ! whether or not to created induced spes from sum over "core"
	real :: ecore

end module information

!================= MAIN PROGRAM ==============================================

print*,' '
print*,' Program to reorder single-particle orbits in .int files '
print*,' '

call open_sps_file('i')
call open_sps_file('f')
call mapspaces
call openhamiltonian
call readhamiltonian
call writehamiltonian

end

!================= END MAIN PROGRAM ==========================================

subroutine open_sps_file(ichar)
	
	use information
	implicit none
	character :: ichar  ! determines if initial or final space
	
	type (orb),pointer :: curorbqn(:)
	integer :: iorb,norb
    character*20 :: filename
	character*3 :: dummy
	integer  :: ilast
	logical  :: success
	integer :: ierr
	real    :: xn,xl,xj
	integer :: w
	
	
	success = .false.
	
	do while(.not.success)
	
  	   if(ichar=='i')then
	       print*,' Enter name of INITIAL .sps file (leave off extension )'	 
	   else
	       print*,' Enter name of FINAL .sps file (leave off extension )'	 
  	   end if
	   read(5,'(a)')filename
	   ilast = index(filename,' ')-1
	   open(unit=1,file=filename(1:ilast)//".sps",status='old',iostat=ierr)
	   if(ierr/=0)then
		   print*,filename(1:ilast)//".sps",' does not exist '
	   else
		   success=.true.
	   end if
    end do
   	read(1,*)dummy   ! SHOULD BE 'iso'
	read(1,*)norb
	if(ichar=='i')then
		allocate(orbqni(norb))
		curorbqn=> orbqni
		norbiti = norb
	else
		allocate(orbqnf(norb))
		curorbqn=> orbqnf
		norbitf = norb
	end if
	do iorb =1,norb
		read(1,*)xn,xl,xj,w
		curorbqn(iorb)%nr = nint(xn)
		curorbqn(iorb)%l  =  nint(xl)
		curorbqn(iorb)%j = nint(xj*2)
        curorbqn(iorb)%map = 0         ! initialize
		curorbqn(iorb)%w   = w        ! used to determine core
		curorbqn(iorb)%core = .false.   ! not always used
	end do
	close(1)	

	return
	
end subroutine open_sps_file

!=============================================================================
!
!  maps initial <--> final s.p. orbits
!  also: determines "core"
!
subroutine mapspaces
	use information
	implicit none
	
	integer orbi,orbf
	logical :: foundmap
	logical :: allsameW
	integer :: w0,maxW
	character :: ychar
	

	do orbf = 1,norbitf
		foundmap = .false.
		do orbi = 1,norbiti
			if(orbqni(orbi)%map/=0)cycle  ! skip those already mapped
			if(orbqni(orbi)%nr ==orbqnf(orbf)%nr  .and.  & 
  			     orbqni(orbi)%l ==orbqnf(orbf)%l .and.   &
			      orbqni(orbi)%j ==orbqnf(orbf)%j )then
				  foundmap = .true.
				  orbqni(orbi)%map = orbf
				  orbqnf(orbf)%map = orbi
				  exit
		    end if
		end do
		if(.not.foundmap)then
			print*,' Could not find a matching state '
			print*,orbf, orbqnf(orbf)%nr, orbqnf(orbf)%l,orbqnf(orbf)%j
			stop
		end if
	end do ! orbi

	w0 = orbqni(1)%w
	allsameW = .true.
	do orbi = 1,norbiti
	   if(orbqni(orbi)%w/=w0)allsameW = .false.
   end do
	
	if(allsameW)then
		print*,' All initial s.p. orbits have same weighting'
		print*,' Cannot determine a core for induced s.p.e.s '
	else
		print*,' Do you want to sum over core for induced s.p.e.s (y/n)?'
		read(5,'(a)')ychar
		if(ychar=='y' .or. ychar=='Y')then
			inducedspes = .true.
			print*,' Enter max W for inclusion in core '
			print*,' (You may have to check the initial .sps file)'
			
!.......... SEARCH FOR LIKELY MAXW............................
            maxW = 10000
			do orbf=1,norbitf
				orbi = orbqnf(orbf)%map
                      !      print*,orbf,orbi,orbqni(orbi)%w
				maxW = min(maxW,orbqni(orbi)%w)
			end do
			print*,' Suggested max W = ',maxW-1				
			read*,maxW
			print*,' '
			print*,' These orbits are included in the core '
			print*,'  N   L  2J '
			do orbi = 1,norbiti
				if(orbqni(orbi)%w <= maxW)then
					orbqni(orbi)%core=.true.
					write(6,'(3i4)')orbqni(orbi)%nr,orbqni(orbi)%l,orbqni(orbi)%j
					
				else
					orbqni(orbi)%core = .false.
					
				end if
				
			end do
			
		else
			inducedspes = .false.
		end if
		
	end if
	return
end subroutine mapspaces

!=============================================================================
!
!  open Hamiltonian file; check for compatibility
!

subroutine openhamiltonian
	implicit none
    character*60 :: filename
	character*3 :: dummy
	integer  :: ilast
	logical  :: success
	integer :: ierr
	
	success=.false.
	do while(.not.success)
		print*,' Enter name of initial .int interaction file '
		read(5,'(a)')filename
		ilast = index(filename,' ')-1
		open(unit=3,file=filename(1:ilast)//".int",status='old',iostat=ierr)
		if(ierr/=0)then
			print*,filename(1:ilast)//".int",' does not exist '
		else
			success=.true.
		end if
	end do
        print*,' successfully opened '
	return
	
end subroutine openhamiltonian

!=============================================================================
!
!  count up how many TBMES are in the new Hamiltonian
!
subroutine readhamiltonian
	use information
	implicit none
	character*70 :: dummy
	logical :: endofheader
	integer :: i,n,jj,tt
	real    :: v
	integer(4) :: nlab(4)
	logical :: keep
	integer ai,bi,af
	integer :: delta
	integer :: k
	real    :: tbmescale
	
	tbmescale = 1.0
	
	print*,' '
	print*,' Enter overall scaling for TBMEs (choose 1 if you dont want to scale) '
	read*,tbmescale
	
    print*,' '
!--------- READ PAST HEADER-------	
	endofheader=.false.
	do while(.not.endofheader)
		read(3,'(a)')dummy
		if(dummy(1:1)/='!' .and. dummy(1:1)/='#')then
			backspace(3)
			endofheader=.true.
		else
			print*,dummy
		end if
	end do
	allocate(spe(norbiti))
	read(3,*)ntbmei,(spe(i),i=1,min(norbiti,10))
	
	if(norbiti>10)then
		k = 10
		do while(k < norbiti)
			read(3,*)(spe(k+i),i=1,min(norbiti-k,10))
			k = k+10
			
		end do
		
	end if
	
	ntbmef = 0
	allocate(vtbme(ntbmei))
	allocate(qlist(ntbmei,6))
!  .... SET UP CORE ENERGY........
	ecore = 0.0
	if(inducedspes)then
       do i = 1,norbiti
		   if(orbqni(i)%core)ecore = ecore+ spe(i)*2*(orbqni(i)%j+1)		
	   end do
	end if
	print*,' Ecore so far ',ecore
!.... LOOP OVER MATRIX ELEMENTS
	do i = 1,ntbmei
		read(3,*)(nlab(n),n=1,4),jj,tt,v
		v=v*tbmescale
!----------- CHECK CONSISTENCY ----------------		
        if ( jj > (orbqni(nlab(1))%j + orbqni(nlab(2))%j)/2 .or. & 
		jj > (orbqni(nlab(3))%j + orbqni(nlab(4))%j)/2 .or. & 
		jj < abs(orbqni(nlab(1))%j - orbqni(nlab(2))%j)/2 .or. & 
		jj < abs(orbqni(nlab(3))%j - orbqni(nlab(4))%j)/2 )then
		print*,' mismatch in quantum numbers # ',i
		print*,nlab(:),jj,tt
		stop
  	    end if
        if ( (-1)**( orbqni(nlab(1))%l + orbqni(nlab(2))%l + orbqni(nlab(3))%l + orbqni(nlab(4))%l   )/=1)then
			print*,' mismatch in parity # ',i
			print*,nlab(:),jj,tt			
			stop
		end if

!----------- CHECK MAPPING ------------------
        keep = .true.
		do n =1,4
			if( orbqni(nlab(n))%map==0)keep=.false.
		end do
		if(keep)then
			ntbmef = ntbmef+1
			vtbme(ntbmef)=v
			do n = 1,4
				qlist(ntbmef,n)=orbqni(nlab(n))%map
			end do
			qlist(ntbmef,5)=jj
			qlist(ntbmef,6)=tt
		end if
!------------ OPTIONAL: ADD TO INDUCED SINGLE-PARTICLE-ENERGIES -------
        if(inducedspes)then
			if(nlab(1)==nlab(3) .and. nlab(2)==nlab(4) )then
!	----- INDUCED SPEs --------------
				if( (orbqni(nlab(1))%core .and. .not.  orbqni(nlab(2))%core .and. orbqni(nlab(2))%map > 0) .or.  & 
				    (orbqni(nlab(2))%core .and. .not.  orbqni(nlab(1))%core .and. orbqni(nlab(1))%map > 0 ))then
				     if(orbqni(nlab(1))%core )then
						 bi = nlab(1)
						 ai = nlab(2)
						 af = orbqni(nlab(2))%map

					 else
						 bi = nlab(2)
						 ai = nlab(1)
						 af = orbqni(nlab(1))%map						 
						 
					 end if    
					 spe(ai)=spe(ai)+ 0.5*(2*tt+1)*(2*JJ+1)*v/float( orbqnf(af)%j+1)
				
			    end if
!------------- INDUCED CORE ENERGY --------------
                if( orbqni(nlab(1))%core .and. orbqni(nlab(2))%core )then
				  bi = nlab(2)
				  ai = nlab(1)				
				  delta=1
				  if(ai==bi)delta=2
				  ecore = ecore + (2*jj+1)*(2*tt+1)*v*float(delta)!/2
				  	
				end if
				
				
			end if
			
		end if
		
	end do
	
	print*,' Kept ',ntbmef,' matrix elements out of ',ntbmei
	if(inducedspes)then
		print*,' Energy of exclude core = ',ecore
	end if
	close(3)
	return
end subroutine readhamiltonian

!=============================================================================

subroutine writehamiltonian
	use information
	implicit none
    character*70 :: filename
	character*1 :: ychar
	integer  :: ilast
	logical  :: success
	integer :: ierr
	integer :: iorbf
	integer :: i,j
    integer nline,nres
	
	success=.false.
	do while(.not.success)
		print*,' Enter name of final .int interaction file '
		read(5,'(a)')filename
		ilast = index(filename,' ')-1
		open(unit=33,file=filename(1:ilast)//".int",status='new',iostat=ierr)
		if(ierr/=0)then
			print*,filename(1:ilast)//".int",' already exists; overwrite (y/n)? '
			read(5,'(a)')ychar
			if(ychar=='y' .or. ychar=='Y')then
				open(unit=33,file=filename(1:ilast)//".int",status='old',iostat=ierr)
				success=.true.
			end if
		else
			success=.true.
		end if
	end do	
	
	write(33,'("!SP ORB  N  L  J" )')
	do iorbf = 1, norbitf
		write(33,'("!SP ",4i3,"/2")')iorbf,orbqnf(iorbf)%nr,orbqnf(iorbf)%l,orbqnf(iorbf)%j
	end do  ! iorbf 
	if(inducedspes)then
		write(33,'("! CORE ENERGY =",f10.4)')ecore
	end if
        print*,orbqnf(:)%map
        do j = 1,norbitf
           print*,spe(orbqnf(j)%map)
        end do
	write(33,'(i8,10f9.4)')ntbmef,(spe(orbqnf(j)%map),j=1,min(10,norbitf))
	print*,norbitf,' final s.p. orbits '
	if(norbitf > 10)then
        nline = norbitf/10
        nres = norbitf-nline*10
		print*,nline,nres
        if(nline > 1)then
           do j = 2,nline
               write(33,'(10f9.4)')(spe(orbqnf(i+(j-1)*10)%map),i=1,10)
            end do
        end if
        if(nres > 0) write(33,'(10f9.4)')(spe(orbqnf(i+(nline)*10)%map),i=1,nres)
		print*, ' Testing sum rule ', nres+nline*10,norbitf
		
	end if
	do i =1,ntbmef
		write(33,'(4i4,2x,2i4,2x,f10.5)')(qlist(i,j),j=1,6),vtbme(i)
	end do
			
	close(33)
	
	return
end subroutine writehamiltonian

