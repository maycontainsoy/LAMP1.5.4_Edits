module hamlib
	use coupledinteraction
	use nodeinfo ! Added by SML
	implicit none
	
!!---------- uncoupled TBMEs ------------------
	integer :: nmatpp,nmatnn
	integer :: nmatpn
	integer :: nmatXX	! does this replace nmatpp/nn?

	real,allocatable  :: hmatxx(:,:)
	integer, allocatable :: hmatorbxx(:,:,:)

	real, allocatable :: hmatpn(:)
	integer, allocatable :: hmatorbPN(:,:)

!------------ coding of the pairs -------------
!  take states i,j,k,l  from list in hspsqn
!
!  We encode a_i a_j, or a^+_i a^+_j, as follows:
!
!  assume i >= j
!
!  then  (i,j) =>  i(i-1)/2 +j 
!!-----------------------------

!!------ uncoupled SPEs
	real, allocatable :: speunX(:,:,:),speunXshift(:,:,:)

!------------ SPE shifts
	real, allocatable :: speshift(:)  !ULTIMATELY ALLOW FOR different for p,n

	contains

!========================================================================
!======================= SUBROUTINES ===========================


!.....................................................................!
! Subroutine: ham_boss 
!.....................................................................!
	subroutine ham_boss
!	use interaction,only:usenewham
!	use hamlib
		IMPLICIT NONE
		character :: choicechar
		INTEGER :: ierr 

		if(usenewham)then
			choicechar = 'n'
			IF (myMPIrank == root) THEN 
				do while(choicechar /= 'i' .and. choicechar/='x')
					print*,' '
					print*,' Do you want to read in files in (i) isospin conserving format (old format) '
					print*," or (x) explicit proton-neutron (xpn) format?  (enter i/x)"
					read(5,'(a)')choicechar
					if(choicechar=='I')choicechar='i'
					if(choicechar=='X')choicechar='x'
					if(choicechar /= 'i' .and. choicechar/='x')then
						print*,choicechar,' is not an acceptable option '
					else
						exit
					end if
				end do ! choicechar
			END IF ! myMPIrank == root 
			CALL MPI_BARRIER(icomm,ierr)
			CALL MPI_BCAST(choicechar,1,MPI_CHARACTER,root,icomm,ierr)
	
			select case (choicechar)
			case('i')
				call setup4tbmes			! in LAMP_interact.f90, no MPI
				!IF (myMPIrank == root) 
				call readvtbme				! in LAMP_interact.f90, currently stuck here
				call uncoupleXXmaster	! in this file
				call undoSPE					! in this file 
				call uncouplePNmaster	! in this file
			case('x')
				call make_hvec_xpn
				call uncoupleXXmaster
				call undoSPE
				call uncouplePNmaster		
			CASE DEFAULT 
				print*,' should not have gotten here '
				print*,choicechar
				stop
			end select
		else   ! usenewham ! OLD WAY
			call setup4tbmes
			call readvtbme
			call uncoupleXXmaster
			call undoSPE
			call uncouplePNmaster
		end if ! usenewham
	  return
	end subroutine ham_boss

!.....................................................................!
! Subroutine: make_hvec_xpn
!.....................................................................!
	subroutine make_hvec_xpn
!	use hamiltonian
!	use observables
		implicit none
		
		character*55 :: intfilename 
		character*1  :: ychar
		integer :: ilast 
		integer(4) :: ierr 
		logical :: finished
		logical :: success,normalizedpn
		logical :: emptyham,pnnormalized
	
!	call old_orbinfo2new
		usepnform = .true.
		print*,' '            
		print*,' Note: interaction files must be in explicit proton-neutron formalism '
		print*,' In this format proton and neutron orbits are sequential and do not overlap, '
		print*,' E.g., proton orbits are 1,2,3 and neutron orbits are 4,5,6. '
		print*,' '
		call setup_tbme_XXcouples(1,.false.) 
		call setup_tbme_XXcouples(1,.true.) 
		call set_tbme_XXarray(1) 
		call setup_tbme_XXcouples(2,.false.) 
		call setup_tbme_XXcouples(2,.true.) 
		call set_tbme_XXarray(2) 
		call setup_tbme_PNcouples(.false.) 
		call setup_tbme_PNcouples(.true.) 
		call set_tbme_PNarray
		finished = .false.
		emptyham = .true.

		do while (.not.finished)
			success = .false.
			do while(.not.success) 
				write(6,*)' Enter interaction file name  (leave off .int extension)'
				write(6,*)' (Enter "end" to stop reading in file )'
				read(5,'(a)')intfilename 
				ilast = index(intfilename,' ')-1         
				if(intfilename(1:3)=='end' .or. intfilename=='END')then
					finished = .true.
					exit
				end if
				
				inquire(file=intfilename(1:ilast)//'.int',exist=success) 
				if(success)then 
					open(unit=1,file=intfilename(1:ilast)//'.int',status='old') 
					print*,' Are the pn matrix element normalized (default) or unnormalize ? (n/u)'
					read(5,'(a)')ychar
					if(ychar.eq.'u'.or.ychar=='U')then
						normalizedpn=.false.
					else
						normalizedpn = .true.
					end if
					call readin_tbmes(normalizedpn,emptyham)
					close(1) 
				else ! success
					print*,' file ',intfilename(1:ilast),'.int does not appear to exist '
					cycle
				end if ! success
			end do ! not success
		end do ! not finished
		
		if(emptyham)then
			print*,' WARNING! No interaction matrix elements found! '
		end if
!	call readexternalfields
!	call count_create_XXuncoupled(1,hamflag,.false.)
!	call count_create_XXuncoupled(1,hamflag,.true.)
!	call count_create_XXuncoupled(2,hamflag,.false.)
!	call count_create_XXuncoupled(2,hamflag,.true.)	
!	call count_create_PNuncoupled(hamflag,.false.)
!	call count_create_PNuncoupled(hamflag,.true.)

		RETURN
	END SUBROUTINE make_hvec_xpn

!.....................................................................!
! Subroutine: uncoupleXXmaster
!.....................................................................!
	subroutine uncoupleXXmaster
!		  use system_parameters
!		  use interaction
		use psis,only:numneut,numprot
		implicit none

		integer mmax
		integer nuncpairs
		integer, pointer :: sppair(:,:),mvpair(:),parpair(:)
		integer nutbmeXX
		IF (myMPIrank == root) PRINT*,' about to decouple '
		
		call countcreateuncoupledpairs(.true.,mmax,nuncpairs, sppair,mvpair,parpair)		! in this file
		call countcreateuncoupledpairs(.false.,mmax,nuncpairs, sppair,mvpair,parpair)		! in this file 
		call countuncoupledtbmeXX(nuncpairs,mvpair,parpair)															! in this file

		if(usepnform)then
			if(numprot> 1)call untbmeXXnew(1,mmax,nuncpairs,sppair, mvpair,parpair)				! in this file
			if(numneut> 1)call untbmeXXnew(2,mmax,nuncpairs,sppair, mvpair,parpair)				! in this file
		else ! usepnform
			call untbmeXX(mmax,nuncpairs,sppair,mvpair,parpair)														! in this file
		end if ! usepnform
		deallocate(sppair,mvpair)
		IF (myMPIrank == root) PRINT*,' decoupled '
		RETURN
	END SUBROUTINE uncoupleXXmaster

!.....................................................................!
! Subroutine: countcreateuncoupledpairs
!.....................................................................!
	subroutine countcreateuncoupledpairs(countflag,jmax,nuncpairs, sppair,mvpair,parpair)
		use spstate
		use sporbit
!		use interaction

		implicit none

		integer :: jmax
		integer :: j
		integer :: itbme
		integer :: asp,bsp
		integer :: ma,mb
		integer :: nuncpairs
		real :: tol
		logical :: countflag
		integer, pointer :: sppair(:,:),mvpair(:),parpair(:)

		tol = 1.0e-5

!-------------------- find maxJ = max M
		if(countflag)then
			jmax = -99
			do asp = 1,numorb
				jmax = max(jmax,orbqn(asp)%j)
			end do ! asp
		endif ! countflag 

!........... COUNT UP # OF UNCOUPLED PAIRS 
		nuncpairs = 0
		do asp = 2,nsps 
			ma = spsqn(asp)%m
			do bsp = 1,asp-1
				mb = spsqn(bsp)%m
				if(abs(ma+mb) <= 2*jmax)then
					nuncpairs = nuncpairs + 1
					if(.not.countflag)then
						sppair(nuncpairs,1) = asp
						sppair(nuncpairs,2) = bsp
						mvpair(nuncpairs) = (ma+mb)/2
						parpair(nuncpairs)= (-1)**(spsqn(asp)%l + spsqn(bsp)%l)
					endif ! not countflag 
				endif ! abs
			enddo  ! loop over bsp
		enddo  ! loop over asp

		IF (myMPIrank == root) PRINT*,' there are ',nuncpairs,' uncoupled pairs '
		if(countflag)then
			allocate(sppair(nuncpairs,2),mvpair(nuncpairs),parpair(nuncpairs))
		endif

		RETURN
	END SUBROUTINE countcreateuncoupledpairs

!.....................................................................!
! Subroutine: countuncoupledtbmeXX
!.....................................................................!
	subroutine countuncoupledtbmeXX(nuncpairs,mvpair,parpair)
		use spstate
		use sporbit
!		use interaction

		implicit none

		integer :: nuncpairs
		integer :: m
		integer :: i,j
		integer,pointer :: mvpair(:),parpair(:)
		integer :: parity
		integer :: a,b,c,d

		nmatXX = 0

		do i = 1,nuncpairs
			m = mvpair(i)
			parity=parpair(i)
!		write(99,*)i,mvpair(i),parpair(i)
			do j = 1,i
				if(parity/=parpair(j))cycle
				if(m == mvpair(j))nmatXX = nmatXX +1
			enddo ! j 
		enddo ! i 

		IF (myMPIrank == root) print*,' There are ',nmatxx,' PP/NN uncoupled TBMEs '
		if (allocated(hmatxx).or.allocated(hmatorbxx)) then
			deallocate(hmatxx,hmatorbxx)
		end if
		if(usepnform)then
			allocate(hmatxx(2,nmatxx),hmatorbxx(2,nmatxx,4))
		else
			allocate(hmatxx(1,nmatxx),hmatorbxx(1,nmatxx,4))
		end if
		hmatXX = 0.d0

		return
	end subroutine countuncoupledtbmeXX

!.....................................................................!
! Subroutine: untbmeXX
!.....................................................................!
	subroutine untbmeXX(jmax,nuncpairs,sppair,mvpair,parpair)
!		use interaction
		use spstate

		implicit none
		integer :: jmax
		integer :: nuncpairs
		integer,pointer :: sppair(:,:),mvpair(:),parpair(:)
		integer :: a,b,c,d
		integer :: asp,bsp,csp,dsp
		integer :: ma,mb,mc,md
		integer :: ja,jb,jc,jd
		integer :: m
		integer :: imatXX
		integer :: i,j
		integer :: ipair,jpair
		integer :: indx
		integer :: JJ
		real :: vtmp
		real :: cleb

		imatXX = 0
		do i = 1,nuncpairs
			m = mvpair(i)
			asp = sppair(i,1) 
			bsp = sppair(i,2)

			a = spsqn(asp)%orb 
			b = spsqn(bsp)%orb
			if(b > a)then
				print*,' wrong order ',a,b
				stop
			endif

			ja = spsqn(asp)%j
			jb = spsqn(bsp)%j
			ma = spsqn(asp)%m
			mb = spsqn(bsp)%m
			ipair = a*(a-1)/2 + b
			do j = 1,i
				if(parpair(i)/=parpair(j))cycle
				if(m == mvpair(j))then
					imatXX = imatXX + 1

!----------------- EXTRACT THE CORRESPONDING STATES
					csp = sppair(j,1)  
					dsp = sppair(j,2)

					c = spsqn(csp)%orb
					d = spsqn(dsp)%orb
					if(d > c)then
						print*,' wrong order c d ',c,d
						stop
					endif

					jc = spsqn(csp)%j
					jd = spsqn(dsp)%j
					mc = spsqn(csp)%m
					md = spsqn(dsp)%m
					jpair = c*(c-1)/2+d

					if(jpair <= ipair)then 
						indx = ipair*(ipair-1)/2+jpair
					else
						indx = jpair*(jpair-1)/2+ipair
					endif

					vtmp = 0.0
					do JJ = vtbme(indx)%jmin,min(jmax,vtbme(indx)%jmax)
						vtmp = vtmp+ vtbme(indx)%v(jj,1) &
						     *zeta(a,b)*zeta(c,d) &
						     *cleb(ja,ma,jb,mb,2*jj,2*m)*cleb(jc,mc,jd,md,2*jj,2*m)
					enddo  ! loop over j

!---------------- KEEP FROM DOUBLE COUNTING FOR DIAGONAL 
					if(asp == csp .and. bsp == dsp)vtmp = vtmp*.5
					hmatXX(1,imatXX) = vtmp
					if(usepnform)hmatXX(2,imatXX)=vtmp
					hmatorbXX(1,imatxx,1) = asp
					hmatorbXX(1,imatxx,2) = bsp
					hmatorbXX(1,imatxx,3) = csp
					hmatorbXX(1,imatxx,4) = dsp
				endif
			enddo ! loop over j
		enddo    ! loop over i

		return
	end subroutine untbmeXX

!.....................................................................!
! Subroutine: untbmeXXnew
!.....................................................................!
	subroutine untbmeXXnew(it,jmax,nuncpairs,sppair,mvpair,parpair)
!		use interaction
		use spstate
!		use coupledinteraction

		implicit none
		integer :: it  !species
		integer :: jmax
		integer :: nuncpairs
		integer,pointer :: sppair(:,:),mvpair(:),parpair(:)
		integer :: a,b,c,d
		integer :: asp,bsp,csp,dsp
		integer :: ma,mb,mc,md
		integer :: ja,jb,jc,jd
		integer :: m
		integer :: imatXX
		integer :: i,j
		integer :: ipair,jpair
		integer :: indx
		integer :: JJ,JJmin,JJmax
		real :: vtmp
		real :: cleb

		imatXX = 0
		do i = 1,nuncpairs
			m = mvpair(i)
			asp = sppair(i,1) 
			bsp = sppair(i,2)

			a = spsqn(asp)%orb 
			b = spsqn(bsp)%orb

			if(b > a)then
				print*,' wrong order ',a,b
				stop
			endif

			ja = spsqn(asp)%j
			jb = spsqn(bsp)%j
			ma = spsqn(asp)%m
			mb = spsqn(bsp)%m
			ipair = a*(a-1)/2 + b
			do j = 1,i
				if(parpair(i)/=parpair(j))cycle
! check parity
				if(m == mvpair(j))then

!----------------- EXTRACT THE CORRESPONDING STATES
					csp = sppair(j,1)  
					dsp = sppair(j,2)

					c = spsqn(csp)%orb
					d = spsqn(dsp)%orb
			  ! check parity
					if((-1)**(orbqn(a)%l+orbqn(b)%l+orbqn(c)%l+orbqn(d)%l)/= 1)cycle
					imatXX = imatXX + 1
					if(asp==0 .or. bsp==0 .or. csp== 0 .or. dsp ==0)then
						print*,' bad sp indices ',asp,bsp,csp,dsp
						stop
					end if
					if(d > c)then
						print*,' wrong order c d ',c,d
						stop
					endif

					jc = spsqn(csp)%j
					jd = spsqn(dsp)%j
					mc = spsqn(csp)%m
					md = spsqn(dsp)%m
					vtmp = 0.0
					jjmin = max(abs(ja-jb),abs(jc-jd))/2
					jjmax = min(ja+jb,jc+jd)/2
					do JJ = jjmin,jjmax
						vtmp = vtmp+ twobodyme(a,b,c,d,JJ,it)*zeta(a,b)*zeta(c,d) &
							   *cleb(ja,ma,jb,mb,2*jj,2*m)*cleb(jc,mc,jd,md,2*jj,2*m)
					enddo  ! loop over JJ 

!---------------- KEEP FROM DOUBLE COUNTING FOR DIAGONAL 
					if(asp == csp .and. bsp == dsp)vtmp = vtmp*.5
					hmatXX(it,imatXX) = vtmp
			  
!			  if(it==1)print*,imatXX, vtmp
					hmatorbXX(it,imatxx,1) = asp
					hmatorbXX(it,imatxx,2) = bsp
					hmatorbXX(it,imatxx,3) = csp
					hmatorbXX(it,imatxx,4) = dsp
				endif
			enddo ! loop over j
		enddo    ! loop over i
		nmatXX = imatxx
		! print*,' testing XX count ',nmatXX

		return
	end subroutine untbmeXXnew

!.....................................................................!
! Subroutine: uncouplePNmaster
!.....................................................................!
	subroutine uncouplePNmaster
		use spstate
!		use interaction
		implicit none

		if(usepnform)then
			call countcreatetbmePN_xpn(.true.)
			call countcreatetbmePN_xpn(.false.)
		else
			call countcreatetbmePN(.true.)
			call countcreatetbmePN(.false.)
		end if

		return
	end subroutine uncouplePNmaster

!.....................................................................!
! Subroutine: countcreatetbmePN
!.....................................................................!
	subroutine countcreatetbmePN(countflag)
		use spstate
!		use interaction
		implicit none

		logical :: countflag
		integer :: nct
		integer :: pa,pc,nb,nd
		integer :: a,b,c,d
		integer :: jpa,jpc,mpa,mpc
		integer :: jnb,jnd,mnb,mnd
		integer :: Mp, Mn
		integer :: parP,parN
		integer :: Jmin,Jmax,J,m
		logical :: phaseab,phasecd
		integer :: fact0,fact1
		integer :: pair1,pair2,indx
		real :: vtmp

		INTEGER :: ierr ! for MPI communication

!------------ FUNCTIONS CALL
		real :: cleb !,zeta
		nct = 0

!------------- LOOP OVER PROTON OPS
		do pa = 1,nsps
			a   = spsqn(pa)%orb
			jpa = spsqn(pa)%j
			mpa = spsqn(pa)%m
			do pc = 1, pa !
				c    = spsqn(pc)%orb
				jpc  = spsqn(pc)%j
				mpc  = spsqn(pc)%m                        
				Mp   = mpa - mpc
				parP = spsqn(pa)%par*spsqn(pc)%par

!------------- loop over NEUTRON OPS
				do nb = 1,nsps
					b   = spsqn(nb)%orb
					jnb = spsqn(nb)%j
					mnb = spsqn(nb)%m
					do nd = 1,nsps
						if(pc ==pa .and. nd > nb)cycle  
						d = spsqn(nd)%orb
						jnd =spsqn(nd)%j
						mnd = spsqn(nd)%m
						Mn = mnb - mnd
						if( Mp + Mn /= 0)cycle
						parN = spsqn(nb)%par*spsqn(nd)%par
						if( parN /= parP)cycle

!----------------------- FIND INDEX OF MATRIX ELEMENT --------
!                        plus phases
!               standard is a >= b, c >= d
						if(a >= b)then
							pair1 = a*(a-1)/2 + b
							phaseab = .false.
						else
							pair1 = b*(b-1)/2 + a			
							phaseab = .true.
						endif
						if(c >= d)then
							pair2 = c*(c-1)/2 + d
							phasecd = .false.
						else
							pair2 = d*(d-1)/2 + c
							phasecd = .true.
						endif
						if(pair1 >= pair2)then
							indx = pair1*(pair1-1)/2 + pair2
						else
							indx = pair2*(pair2-1)/2 + pair1
						endif

!------------------------ DECOUPLE MATRIX ELEMENT
						vtmp = 0.0
						m = (mpa + mnb)/2
						jmin = abs(m)
						jmin = max(jmin,vtbme(indx)%jmin)
						Jmax = min( jpa + jnb, jpc+jnd)/2
						jmax = min(jmax,vtbme(indx)%jmax)
						if(jmin > jmax)cycle
						do J = Jmin,Jmax
							fact0 = 1
							fact1 = 1
							if(phaseab)then
								fact0 = (-1)**((jpa+jnb)/2+J)
								fact1 = - fact0
							endif
							if(phasecd)then
								fact0 = fact0*(-1)**((jpc+jnd)/2+J)
								fact1 = -fact1*(-1)**((jpc+jnd)/2+J)
							endif
							vtmp = vtmp+ 0.5*(fact0*vtbme(indx)%v(j,0) + &
							       fact1*vtbme(indx)%v(j,1)) *zeta(a,b)*zeta(c,d) & 
							      *cleb(jpa,mpa,jnb,mnb,2*j,2*m)*cleb(jpc,mpc,jnd,mnd,2*j,2*m)
						enddo  ! loop over J
						if(abs(vtmp) < 0.00001)cycle
						if(pa==pc .and. nb==nd)vtmp =vtmp*0.5  ! prevent double-counting
						nct = nct + 1
						if(.not. countflag)then
							hmatPN(nct) = vtmp
							hmatorbPN(1,nct) = pa
							hmatorbPN(2,nct) = nb
							hmatorbPN(3,nct) = pc
							hmatorbPN(4,nct) = nd
						endif
					enddo ! loop over nd
				enddo  ! loop over nb
			enddo  ! loop over pc
		enddo  ! loop over pa
		if(countflag)then
			nmatPN = nct
			IF (myMPIrank == root) print*,' There are ',nmatPN,' PN matrix elements '
			if(allocated(hmatPN))then
				deallocate(hmatPN)
			end if
			allocate(hmatPN(nct))
			hmatPN = 0.0
			if(allocated(hmatorbPN))then
				deallocate(hmatorbPN)
			end if
			allocate(hmatorbPN(4,nct))
		endif

		! TEMPORARY FIX, MISMATCH IN NMATPN COULD CAUSE ERRORS ELSEWHERE
		! CALL MPI_BARRIER(icomm,ierr)
		! CALL MPI_BCAST(nmatPN,1,MPI_INTEGER,root,icomm,ierr)

		return
	end subroutine countcreatetbmePN

!.....................................................................!
! Subroutine: countcreatetbmePN_xpn
! For use with XPN format
!.....................................................................!
	subroutine countcreatetbmePN_xpn(countflag)
		use spstate
!		use interaction
!		use coupledinteraction
		implicit none

		logical :: countflag
		integer :: nct
		integer :: pa,pc,nb,nd
		integer :: a,b,c,d
		integer :: jpa,jpc,mpa,mpc
		integer :: jnb,jnd,mnb,mnd
		integer :: Mp, Mn
		integer :: parP,parN
		integer :: Jmin,Jmax,J,m
		logical :: phaseab,phasecd
		integer :: fact0,fact1
		integer :: pair1,pair2,indx
		real :: vtmp

!------------ FUNCTIONS CALL
		real cleb !zeta
		nct = 0

!------------- LOOP OVER PROTON OPS
		do pa = 1,nsps
			a = spsqn(pa)%orb
			jpa = spsqn(pa)%j
			mpa = spsqn(pa)%m
			do pc = 1,pa !nsps
				c = spsqn(pc)%orb
				jpc = spsqn(pc)%j
				mpc = spsqn(pc)%m                        
				Mp = mpa - mpc
				parP = spsqn(pa)%par*spsqn(pc)%par

!------------- loop over NEUTRON OPS
				do nb = 1,nsps
					b = spsqn(nb)%orb
					jnb = spsqn(nb)%j
					mnb = spsqn(nb)%m
					do nd = 1,nsps
						if(pc ==pa .and. nd > nb)cycle  
						d = spsqn(nd)%orb
						if((-1)**(orbqn(a)%l+orbqn(b)%l+orbqn(c)%l+orbqn(d)%l)/= 1)cycle
						jnd =spsqn(nd)%j
						mnd = spsqn(nd)%m
						Mn = mnb - mnd
						if( Mp + Mn /= 0)cycle
						parN = spsqn(nb)%par*spsqn(nd)%par
						if( parN /= parP)cycle

!----------------------- FIND INDEX OF MATRIX ELEMENT --------
!                        plus phases
!               standard is a >= b, c >= d
						if(a >= b)then
							pair1 = a*(a-1)/2 + b
							phaseab = .false.
						else
							pair1 = b*(b-1)/2 + a			
							phaseab = .true.
						endif

						if(c >= d)then
							pair2 = c*(c-1)/2 + d
							phasecd = .false.
						else
							pair2 = d*(d-1)/2 + c
							phasecd = .true.
						endif

!------------------------ DECOUPLE MATRIX ELEMENT
						vtmp = 0.0
						m = (mpa + mnb)/2
						jmin = abs(m)
						jmin = max(jmin,abs(jpa-jnb)/2)
						jmin = max(jmin,abs(jpc-jnd)/2)
						Jmax = min( jpa + jnb, jpc+jnd)/2
 !                   jmax = min(jmax,vtbme(indx)%jmax)
						if(jmin > jmax)cycle
						do J = Jmin,Jmax
							fact0 = 1
							fact1 = 1
							if(phaseab)then
								fact0 = (-1)**((jpa+jnb)/2+J)
								fact1 = - fact0
							endif
							if(phasecd)then
								fact0 = fact0*(-1)**((jpc+jnd)/2+J)
								fact1 = -fact1*(-1)**((jpc+jnd)/2+J)
							endif
							vtmp = vtmp+ twobodyme(a,b,c,d,J,0) &
							      *cleb(jpa,mpa,jnb,mnb,2*j,2*m)*cleb(jpc,mpc,jnd,mnd,2*j,2*m)
!     &     *zeta(a,b)*zeta(c,d)  !DON"T DOUBLECOUNT" old version
						enddo  ! loop over J
						if(abs(vtmp) < 0.00001)cycle
						if(pa==pc .and. nb==nd)vtmp = vtmp*0.5  ! prevent double-counting
						nct = nct + 1
						if(.not. countflag)then
							hmatPN(nct) = vtmp  
							hmatorbPN(1,nct) = pa
							hmatorbPN(2,nct) = nb
							hmatorbPN(3,nct) = pc
							hmatorbPN(4,nct) = nd
						endif
					enddo ! loop over nd
				enddo  ! loop over nb
			enddo  ! loop over pc
		enddo  ! loop over pa

		if(countflag)then
			nmatPN = nct
			IF (myMPIrank == root) print*,' There are ',nmatPN,' PN matrix elements '
			if(allocated(hmatPN))then
				deallocate(hmatPN)
			end if
			allocate(hmatPN(nct))
			if(allocated(hmatorbPN))then
				deallocate(hmatorbPN)
			end if
			allocate(hmatorbPN(4,nct))
		endif

		return
	end subroutine countcreatetbmePN_xpn

!.....................................................................!
! Function: zeta
!.....................................................................!
	real function zeta(i,j)
		implicit none

		integer i,j

		zeta = 1.0
		if(i ==j)zeta = sqrt(2.)
		return
	end function zeta

!.....................................................................!
! Subroutine: undpSPE
! At this time assume spes diagonal however can generalize
! allow assume neutron, proton part the same
!.....................................................................!
	subroutine undoSPE
		use sporbit
		use spstate
!		use interaction
!		use coupledinteraction

		integer asp,a
		integer it

		if(allocated(speunX)) then
			deallocate(speunX)
		end if

		allocate(speunX(nsps,nsps,2))
		speunX = 0.0

		if(usepnform)then
			do asp = 1,nsps
				a = spsqn(asp)%orb
				speunX(asp,asp,1) = pspe(a)
				speunX(asp,asp,2) = nspe(a)
			end do
		else
			do it = 1,2
				do asp = 1,nsps
					a = spsqn(asp)%orb
					speunX(asp,asp,it) = spe(a)
				enddo
			enddo
		end if

		return
	end subroutine undoSPE

!.....................................................................!
! Subroutine: undoSPEshifted
! At this time assume SPEs diagonal however can generalize
! allow assume neutron, proton part the same
!.....................................................................!
	subroutine undoSPEshifted
		use sporbit
		use spstate
!		use interaction

		integer :: asp,a
		integer :: it

		allocate(speunXshift(nsps,nsps,2))

		speunXshift = 0.0
		do it = 1,2
			do asp = 1,nsps
				a = spsqn(asp)%orb
				speunXshift(asp,asp,it) = spe(a)+speshift(a)
			enddo
		enddo

		return
	end subroutine undoSPEshifted


	!.....................................................................!
! Subroutine: hmultMPIdistro
!.....................................................................!
	SUBROUTINE hmultMPIdistro 
		!USE nodeinfo
		!USE hamlib 
		IMPLICIT NONE 

		INTEGER :: irank 
		INTEGER :: ierr ! for MPI communication
		INTEGER(8) :: rankchunkXX, rankchunkPN 
		!INTEGER(8) :: rankchunkPP, rankchunkNN, rankchunkPN
		INTEGER, ALLOCATABLE :: nmatXXstart(:), nmatXXstop(:)
		!INTEGER, ALLOCATABLE :: nmatPPstart(:), nmatPPstop(:)
		!INTEGER, ALLOCATABLE :: nmatNNstart(:), nmatNNstop(:) 
		INTEGER, ALLOCATABLE :: nmatPNstart(:), nmatPNstop(:)

		!ALLOCATE(nmatPPstart(0:nMPIranks-1),nmatPPstop(0:nMPIranks-1))
		!ALLOCATE(nmatNNstart(0:nMPIranks-1),nmatNNstop(0:nMPIranks-1))
		ALLOCATE(nmatXXstart(0:nMPIranks-1),nmatXXstop(0:nMPIranks-1))
		ALLOCATE(nmatPNstart(0:nMPIranks-1),nmatPNstop(0:nMPIranks-1))

		IF (nMPIranks == 1) THEN 
			! nmatPPstart(0) = 1 
			! nmatPPstop(0)  = nmatpp 
			! nmatNNstart(0) = 1 
			! nmatNNstop(0)  = nmatnn 
			nmatXXstart(0) = 1 
			nmatXXstop(0) = nmatxx 
			nmatPNstart(0) = 1 
			nmatPNstop(0)  = nmatpn 
			RETURN 
		END IF 

		! rankchunkPP = nmatpp/nMPIranks 
		! rankchunkNN = nmatnn/nMPIranks 
		rankchunkXX = nmatxx/nMPIranks 
		rankchunkPN = nmatpn/nMPIranks 

		DO irank = 0, nMPIranks - 1
			! nmatPPstart(irank) = 1 + irank*rankchunkPP 
			! nmatPPstop(irank)  = (irank + 1)*rankchunkPP 
			! nmatNNstart(irank) = 1 + irank*rankchunkNN 
			! nmatNNstop(irank)  = (irank + 1)*rankchunkNN 
			nmatXXstart(irank) = 1 + irank*rankchunkXX 
			nmatXXstop(irank)  = (irank + 1)*rankchunkXX  
			nmatPNstart(irank) = 1 + irank*rankchunkPN
			nmatPNstop(irank)  = (irank + 1)*rankchunkPN
		END DO 
		! This is to account for division by # MPI ranks not exactly even 
		! nmatPPstop(nMPIranks-1) = nmatpp 
		! nmatNNstop(nMPIranks-1) = nmatnn 
		nmatXXstop(nMPIranks-1) = nmatxx
		nmatPNstop(nMPIranks-1) = nmatpn 

		RETURN 
	END SUBROUTINE hmultMPIdistro 

!.....................................................................!
end module hamlib
!.....................................................................!