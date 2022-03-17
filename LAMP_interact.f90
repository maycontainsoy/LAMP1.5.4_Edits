!
!  data and codes for the proton-neutron part of the interaction
!  many of the routines and ideas borrowed from BIGSTICK
!
!  for efficient mapping, first set up "PNcouples" which 
!  couple up proton and neutron orbits, as well as where their 
!  quantum numbers start and stop --subroutine setup_tbme_PNcouples
!
!  based upon that information, set up and allocate derived type array
!  pnme in set_tbme_PNarray 
!
!  finally matrix elements in explicit pn formalism are read in, 
!  skipping over pp and nn matrix elements
!
!  subroutines in this file:
!     setup_tbme_PNcouples       
!     set_tbme_PNarray 
!     readin_pntbmes
!     master_readin_pntbmes
!
!    * - * - * - * - * - * - * - * - * - *
!  CRUCIAL: How matrix elements are mapped:
!  for a matrix element < ab: J | V | cd: J>
!  where a and c are proton labels, b and d neutron labels
!  form: pair1 =  numorb(2)*(a-1) + b 
!        pair2 = numorb(2)*(c-1) + d 
!   require pair1 >= pair2 (else swap)
!   pair1 and pair2 must have the same parity (else skip)
!   from the derived type PNcouples get iref and istart for a given parity
!  
!    then finally indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart    
!
!  the matrix element is store in pnme(indx)%v(J)
!
!  pn(indx) also has jmin,jmax
!    * - * - * - * - * - * - * - * - * - *
!
!  ALSO:  in the input file, must be in explicit proton neutron format
!  so if originally we have 3 protons orbits and 5 neutron orbits,
!  then labels 1-3 belong to protons and 4-8 to neutrons.
!  Hence in the interaction file
!   1 2 2 3 J T xxxx is a pp matrix element
!   4 5 5 7 J T xxxxx is a nn matrix element
!   1 5 2 7 J T xxxxx is a pn matrix element
!  
!  the code cvtiso2pn.f90 will read in standard isospin format files 
!  and convert to this explicit p-n formalism
!
!  NOTE some files for NuSHell have "unnormalized" pn-matrix elements
!                                           --has to do with a sqrt(2)
!  however in general we will default to "normalized" pn matrix elements.
!  (This means  state | a_p b_n J > in  2-body matrix element is normalized to 1)
!

module coupledinteraction
	use sporbit
	use nodeinfo ! added by SML
	implicit none    
!
!  information on COUPLED two-body pairs 
!
	type coupair_qn
		integer :: par
		integer :: indx
		integer :: ia,ib
		end type coupair_qn

		type coupairinfo
			type (coupair_qn), allocatable :: pairc(:)
			integer :: meref(-1:1), mestart(-1:1)  ! i parity
		end type coupairinfo
		integer,target :: ncouplesPN,ncouplesXX(2)
		type (coupairinfo),target :: PNcouples,XXcouples(2)
        
!----------- coupled  TBMES -----------------
		integer :: npnVmes,nv2bmedim(2)
		type vjs
			integer  :: Jmin,Jmax
			real,pointer :: v(:)
		end type vjs
		type (vjs), allocatable,target :: pnme(:),ppme(:),nnme(:)    
		real,allocatable,target :: pspe(:),nspe(:) 
		integer, allocatable,target :: PNcouplemap(:),PPcouplemap(:),NNcouplemap(:)       
		logical :: usepnform					! added in 1.4.5  July 2019 by CWJ
		logical :: usenewham  =.true.	! added in 1.4.5  July 2019 by CWJ
                          ! to aid in adding XPN format 
		real,allocatable :: spe(:)
		type pair_qn
			integer :: M
			integer :: par
			integer :: W
			integer :: indx
			integer :: ia,ib
		end type pair_qn
		integer :: nvdim			! dimension
		integer :: ntbme			! # of two-body matrix elements
		integer :: norb_allow	! # of allowed orbits
		integer,allocatable :: orblist(:)  

!C----------- coupled  TBMES -----------------
		type vjts
			integer  :: Jmin,Jmax
			real,allocatable :: v(:,:)
			integer  :: orb(4)
		end type vjts
		type (vjts),allocatable :: vtbme(:)  

CONTAINS

!.....................................................................!
! Subroutine: setup4tbmes
! V(a b, c d; J T)
! assume: a >= b, c >= d, a >= d etc. 
!.....................................................................!
	subroutine setup4tbmes
		use sporbit
!		use interaction
		implicit none

		integer :: iorb
		integer :: icount
		integer :: npair
		integer :: pair1,pair2
		integer :: indx
		integer :: jmin,jmax
		integer :: ia,ib,ic,id
		integer :: na,nb,nc,nd
		integer :: ja,jb,jc,jd
		integer :: dstart

!     SOME INTRONS
		!PRINT*, ' Node = ', myMPIrank, ' gets to start of setup4tbmes'
		norb_allow = numorb
		if (numruns.ne.0) then
			deallocate(orblist,spe,vtbme)
		end if
		allocate(orblist(norb_allow))
		allocate(spe(numorb))

		icount = 0

		do iorb = 1,numorb
			icount = icount+1
			orblist(icount) = iorb
		enddo

!------------ figure out approx dimension of unique matrix elements
		npair = norb_allow*(norb_allow+1)/2  ! # of pairs
		nvdim = npair*(npair+1)/2
		allocate(vtbme(nvdim))

!---------------FIND MIN,MAX J
		do ia = 1,norb_allow
			na = orblist(ia)
			ja = orbqn(na)%j
			do ib = 1,ia
				nb = orblist(ib)
				jb = orbqn(nb)%j
				pair1 = ia*(ia-1)/2 + ib
				do ic = 1,ia
					nc = orblist(ic)
					jc = orbqn(nc)%j
					if(ia ==ic)then
						dstart = ib
					else
						dstart = ic
					endif
					do id = 1,dstart 
						nd = orblist(id)
						jd = orbqn(nd)%j
						pair2 = ic*(ic-1)/2+id
						indx = pair1*(pair1-1)/2+pair2
!                write(18,202)indx,ia,ib,ic,id,pair1,pair2
202         format(i3,3x,4i2,3x,2i2)
						jmax = min((ja+jb)/2, (jc+jd)/2)
						jmin = max(abs(ja-jb)/2, abs(jc-jd)/2)
						vtbme(indx)%jmax = jmax
						vtbme(indx)%jmin = jmin

						if(jmin <= jmax)then
							allocate(vtbme(indx)%v(jmin:jmax,0:1))
							vtbme(indx)%v = 0.0
						endif

!------------------ STORE ORBITAL INDICES 
						vtbme(indx)%orb(1) = ia
						vtbme(indx)%orb(2) = ib
						vtbme(indx)%orb(3) = ic
						vtbme(indx)%orb(4) = id
					enddo  ! loop over id
				enddo  ! loop over ic
			enddo  ! loop over ia
		enddo  ! loop over ib

		!PRINT*, ' Node = ', myMPIrank, ' gets to end of setup4tbmes'

	end subroutine setup4tbmes

!.....................................................................!
! Subroutine: readvtbme 
! Reads in TBME from interaction file. Scales and orders the matrix 
! elements accordingly and populates vtbme(indx)%v(j,t)
!.....................................................................!
	SUBROUTINE readvtbme 
		USE sporbit 
		USE psis,ONLY:numprot,numneut 
		IMPLICIT NONE 
	
	!---- FILE CONTROL ----------------------------------------------------
		CHARACTER(LEN=25) :: filename 
		CHARACTER(LEN=1)  :: achar 
		CHARACTER(LEN=70) :: title 
		INTEGER :: ilast 
	
	!---- INTERNAL VARIABLES ----------------------------------------------
		INTEGER :: ia, ib, ic, id 
		INTEGER :: na, nb, nc, nd 
		INTEGER :: pair1, pair2, indx 
		INTEGER :: j,t 
		INTEGER :: phase 
		INTEGER :: nme 			! Number of MEs in .int file
		INTEGER :: dw 
		INTEGER :: ierr 		! For MPI communications
	
		REAL :: V 
		REAL :: spscale, a, b, x, vscale 	! For scaling interactions
		REAL,ALLOCATABLE :: spetmp(:)
	
	!---- DUMMY COUNTERS --------------------------------------------------
		INTEGER :: i, m, mmax 
		INTEGER :: L 
		INTEGER :: numorbplus 
	
		LOGICAL :: success 
		LOGICAL :: smint 		! successor to .int 
		LOGICAL :: finished, autoscale 
	
	!---- BEGIN -----------------------------------------------------------
		spe = 0. 
		finished = .FALSE. 
		DO WHILE (.NOT.finished)
	!---- OPEN A FILE -----------------------------------------------------
			success = .FALSE. 
			DO WHILE (.NOT.success)
				IF (myMPIrank == root) THEN 
					PRINT*, ' Enter interaction file name (.smint/.int)'
					PRINT*, ' (Enter END to stop)'
					READ(5,'(A)') filename 
				!END IF 
				CALL MPI_BARRIER(icomm,ierr)
				CALL MPI_BCAST(filename,25,MPI_CHARACTER,root,icomm,ierr)
					IF (filename == 'END' .OR. filename == 'end') THEN 
						finished = .TRUE. 
						RETURN 
					END IF ! filename == end or END 
					ilast = INDEX(filename,' ') - 1
				END IF ! myMPI
	!---- ATTEMPT TO OPEN .smint FILE -------------------------------------
				IF (myMPIrank == root) THEN 
					OPEN(UNIT=1,file=filename(1:ilast)//'.smint',STATUS='OLD',ERR=101)
					success = .TRUE. 
					smint  = .TRUE. 
				END IF ! myMPIrank == root
				CALL MPI_BARRIER(icomm,ierr)
				CALL MPI_BCAST(success,1,MPI_LOGICAL,root,icomm,ierr)
				CALL MPI_BCAST(smint,  1,MPI_LOGICAL,root,icomm,ierr)
				CYCLE 
	101		CONTINUE 
	
	!---- ATTEMPT TO OPEN .int FILE ---------------------------------------
				IF (myMPIrank == root) THEN 
					OPEN(UNIT=1,file=filename(1:ilast)//'.int',STATUS='OLD',ERR=102)
					success = .TRUE. 
					smint   = .FALSE. 
				END IF ! myMPIrank == root
				CALL MPI_BARRIER(icomm,ierr)
				CALL MPI_BCAST(success,1,MPI_LOGICAL,root,icomm,ierr)
				CALL MPI_BCAST(smint,  1,MPI_LOGICAL,root,icomm,ierr)
				CYCLE 
	102		CONTINUE 
			END DO ! while .NOT.success 
	
	!---- READ PAST TITLE CARDS -------------------------------------------
			success = .FALSE.
			DO WHILE (.NOT.success)
				IF (myMPIrank == root) THEN 
					READ(1,'(A)') achar 
					BACKSPACE(1)
					IF (achar /= '#' .AND. achar /= '!') THEN 
						success = .TRUE. 
					ELSE 
						READ(1,'(A)') title 
						WRITE(6,*) title 
					END IF ! achar 
				END IF ! myMPIrank == root 
				CALL MPI_BARRIER(icomm,ierr)
				CALL MPI_BCAST(success,1,MPI_LOGICAL,root,icomm,ierr)
			END DO ! while .NOT.success

			IF (myMPIrank == root) READ(1,*) nme 
			CALL MPI_BARRIER(icomm,ierr)
			CALL MPI_BCAST(nme,1,MPI_INTEGER,root,icomm,ierr)
	
			IF (nme < 0) THEN 
				autoscale = .TRUE. 
				IF (myMPIrank == root) THEN 
					PRINT*, ' '
					PRINT*, ' AUTOSCALING two-body matrix elements'
					PRINT*, ' '
				END IF ! myMPIrank == root 
				spscale = 1.
				numorbplus = numorb + 3 
			ELSE ! nme < 0 
				autoscale = .FALSE. 
				numorbplus = numorb 
			END IF ! nme < 0 
			IF (myMPIrank == root) BACKSPACE(1)
	
	!---- ENTER SCALING ---------------------------------------------------
			IF (.NOT.autoscale) THEN 
				IF (myMPIrank == root) THEN 
					PRINT*, ' Enter scaling for spes, A0, A, X ((A0/A)^X) for TBMEs'
					PRINT*, ' (If A or X = 0, then scale by A0)'
					PRINT*, ' (Typically A0 = Acore + 2)'
					READ*, spscale, a, b, x 
				END IF ! myMPIrank == root 

				CALL MPI_BARRIER(icomm,ierr)
				CALL MPI_BCAST(spscale,1,MPI_REAL,root,icomm,ierr)
				CALL MPI_BCAST(a,      1,MPI_REAL,root,icomm,ierr)
				CALL MPI_BCAST(b,      1,MPI_REAL,root,icomm,ierr)
				CALL MPI_BCAST(x,      1,MPI_REAL,root,icomm,ierr)
	
				IF (b == 0. .OR. x == 0.0) THEN 
					vscale = a 
				ELSE 
					IF (b == 0. .OR. x == 0.) THEN 
						vscale = 1 
						IF (myMPIrank == root) PRINT*, ' SCALING set = 1'
					ELSE 
						vscale = (a/b)**x 
					END IF 
				END IF 
			END IF ! .NOT.autoscale 
	
	!---- READ IN SPEs ----------------------------------------------------
			ALLOCATE(spetmp(numorbplus))
			IF (myMPIrank == root) THEN 
				READ(1,*)nme,(spetmp(i),i=1,min(10,numorbplus))
	
				IF (numorbplus > 10) THEN 
					DO m = 10, numorbplus, 10 
						mmax = min(10+m,numorbplus)
						READ(1,*)(spetmp(i),i=1+m,mmax)
					END DO ! m 
				END IF ! numorbplus > 10 
	
				IF (autoscale) THEN 
					spscale = 1. 
					a = spetmp(numorb+2)
					b = spetmp(numorb+1) + numprot + numneut 
					x = spetmp(numorb+3)
					vscale = (a/b)**x 
					PRINT*, ' Autoscale parameters ', a, b, x, vscale 
					nme = ABS(nme)
				END IF ! autoscale
	
				DO i  = 1, numorb 
					spe(i) = spe(i) + spscale*spetmp(i)
				END DO 
			END IF ! myMPIrank == root 

			CALL MPI_BARRIER(icomm,ierr)
			! MIGHT NOT NEED THESE BROADCAST IF SCALING IS NOT NEEDED LATER
			CALL MPI_BCAST(a,1,MPI_REAL,root,icomm,ierr)
			CALL MPI_BCAST(b,1,MPI_REAL,root,icomm,ierr)
			CALL MPI_BCAST(x,1,MPI_REAL,root,icomm,ierr)
			CALL MPI_BCAST(vscale,1,MPI_REAL,root,icomm,ierr)
			! MIGHT NOT NEED THESE BROADCAST IF SCALING IS NOT NEEDED LATER
			CALL MPI_BCAST(spe(:),numorb,MPI_REAL,root,icomm,ierr)
			CALL MPI_BCAST(nme,1,MPI_INTEGER,root,icomm,ierr)

			DEALLOCATE(spetmp)
	
	!---- READ IN TBMEs ---------------------------------------------------
			DO i = 1, nme 
				IF (myMPIrank == root) THEN 
					READ(1,*) ia, ib, ic, id, j, t, v 
	
	!---- CHECK PARITY
					L = orbqn(ia)%l + orbqn(ib)%l + orbqn(ic)%l + orbqn(id)%l 
					IF ((-1)**(L) == -1) THEN 
						PRINT*, ' Error in parity'
						PRINT*, ia, ib, ic, id, j, t, v 
						STOP 
					END IF 
	
	!---- PUT INTO CORRECT ORDER; PICK UP PHASES --------------------------
	!			"CORRECT" ORDER: a >= b, c >= d
					phase = 1
					IF (ia < ib) THEN 
						na = ib 
						ib = ia 
						ia = na 
						phase = (-1)**(J + T + (orbqn(ia)%j + orbqn(ib)%j)/2)
					END IF 
	
					IF (ic < id) THEN 
						nc = id 
						id = ic 
						ic = nc 
						phase = phase*(-1)**(J + T + (orbqn(ic)%j + orbqn(id)%j)/2)
					END IF 
	
					IF (ia < ic .OR. (ia == ic .AND. ib < id)) THEN 
						na = ic 
						nb = id 
						ic = ia 
						id = ib 
						ia = na 
						ib = nb 
					END IF 
	
	!---- CONVERT ---------------------------------------------------------
					na = orblist(ia)
					nb = orblist(ib)
					nc = orblist(ic)
					nd = orblist(id)
					pair1 = ia*(ia-1)/2 + ib 
					pair2 = ic*(ic-1)/2 + id 
					indx = pair1*(pair1-1)/2+pair2
	
	!---- ERROR TRAP ------------------------------------------------------
					IF (j > vtbme(indx)%jmax .OR. j < vtbme(indx)%jmin) THEN 
						PRINT*, ' error in Js ', pair1, pair2, indx 
						PRINT*, ia, ib, ic, id, J, T, v 
						PRINT*, orbqn(ia)%j, orbqn(ib)%j, orbqn(ic)%j, orbqn(id)%j 
						PRINT*, vtbme(indx)%jmax, vtbme(indx)%jmin 
						STOP 
					END IF 
					vtbme(indx)%v(j,t) = vtbme(indx)%v(j,t) + v*vscale*phase
				END IF ! myMPIrank == root 

				CALL MPI_BARRIER(icomm,ierr)
				CALL MPI_BCAST(indx,1,MPI_INT,root,icomm,ierr)
				CALL MPI_BCAST(j   ,1,MPI_INT,root,icomm,ierr)
				CALL MPI_BCAST(t   ,1,MPI_INT,root,icomm,ierr)
				CALL MPI_BCAST(vtbme(indx)%v(j,t),1,MPI_REAL,root,icomm,ierr)
			END DO 
			IF (myMPIrank == root) CLOSE(1) 
		END DO ! while .NOT.finished 
	END SUBROUTINE readvtbme


!---------------------------------------------------------------------!
! Subroutine: setup_tbme_XXcouples
!---------------------------------------------------------------------!
	subroutine setup_tbme_XXcouples(it, create) 
		use sporbit 

		implicit none 
		integer :: it ! = species of XX, either pp or nn 
		logical :: create  ! if count or just create 
		integer :: nrawpairs 
		integer :: a,b 
		integer :: ncount 
		integer :: W 
		integer, pointer :: cmap(:) 
		integer :: indx 
		integer :: par 
		integer :: neven,ieven 
		integer :: i,j,itmp 
		type (coupair_qn) :: cptemp 
		integer :: nv 
		integer :: aerr 
		type (orb),pointer :: xorbqn(:)
!		logical :: allsameparity
!		nrawpairs = numorb(it)*(numorb(it) +1)/2 
		nrawpairs = numorb*(numorb +1)/2 

		if(create)then 
			if( ncouplesXX(it) == 0)return 
			if(.not.allocated(XXcouples(it)%pairc))then
				allocate( XXcouples(it)%pairc( ncouplesXX(it) ) , stat=aerr) 
			endif
			if ( it == 1 ) then 
				if(.not.allocated(PPcouplemap))allocate( PPcouplemap(nrawpairs) , stat=aerr) 
				cmap => PPcouplemap 
				xorbqn => orbqn
!				allsameparity = pallsameparity
			else ! it == 1
				if(.not.allocated(NNcouplemap))allocate( NNcouplemap(nrawpairs) , stat=aerr) 
				cmap => NNcouplemap 
				xorbqn => orbqn
!				allsameparity = nallsameparity
			end if ! it == 1 
			cmap(:) = -1 
		else ! create  
			if ( it == 1 ) then 
				cmap => PPcouplemap 
				xorbqn => orbqn
!				allsameparity = pallsameparity
			else ! it == 1 
				cmap => NNcouplemap 
				xorbqn => orbqn
!				allsameparity = nallsameparity
			end if ! it == 1 
		endif  ! create
 
		ncount = 0 
		do a = 1,numorb !(it) 
			do b = 1,a 
				ncount = ncount + 1 
				if(create)then 
					indx = a*(a-1)/2 + b 
					cmap(indx) = ncount 
					XXcouples(it)%pairc(ncount)%ia = a 
					XXcouples(it)%pairc(ncount)%ib = b 
					XXcouples(it)%pairc(ncount)%indx = indx 
					par = orbqn(a)%par*orbqn(b)%par 
					XXcouples(it)%pairc(ncount)%par = par 
!			print*,'parity ',a,b,par
				end if ! create
			end do ! b 
		end do ! a 

		if(.not. create)then 
			ncouplesXX(it) = ncount 
		end if 
 
!---------------- SORT ! on parity only ------------------- 
		if(create)then 
			if(allsameparity)then 
				XXcouples(it)%mestart(1) = 0 
				XXcouples(it)%meref(1) = 0 
				XXcouples(it)%mestart(-1) =-1 
				XXcouples(it)%meref(-1) =ncouplesXX(it) 
			else 
!............ CLEVERSORT on parity............ 
!             FIRST count up # of even parities 
				neven = 0 
				do i = 1,ncouplesXX(it) 
					if( XXcouples(it)%pairc(i)%par ==1)neven = neven+1 
				end do 
				if(neven == 0 .or. neven == ncouplesXX(it))then 
					if(neven == ncouplesXX(it))then 
						XXcouples(it)%mestart(1) = 0 
						XXcouples(it)%meref(1) = 0 
						XXcouples(it)%mestart(-1) =-1 
						XXcouples(it)%meref(-1) =ncouplesXX(it) 
					else  ! neven 
						XXcouples(it)%mestart(1) =-1 
						XXcouples(it)%meref(1) = 0 
						XXcouples(it)%mestart(-1) =0 
						XXcouples(it)%meref(-1) =0 
					endif ! neven
				else 
					ieven = 0 
					i = 1 
					j = ncouplesXX(it) 
					do while( ieven /= neven) 
						if( XXcouples(it)%pairc(i)%par == -1)then 
							do while( XXcouples(it)%pairc(j)%par ==-1 ) 
								j = j -1 
								if(j <= i)then 
									print*,' whoops missed ' 
									stop 
								endif 
							end do 
!----------------------- SWAP 
							cptemp = XXcouples(it)%pairc(i) 
							XXcouples(it)%pairc(i) = XXcouples(it)%pairc(j) 
							XXcouples(it)%pairc(j) = cptemp 
							cmap(  XXcouples(it)%pairc(i)%indx )  =i 
							cmap(  XXcouples(it)%pairc(j)%indx )  =j 
						endif ! XXcouples 
						ieven = ieven + 1 
						i = i + 1 
					end do ! ieven /= neven
!..................... FIND START, FINISH................ 
					XXcouples(it)%mestart(1) = 0 
					XXcouples(it)%meref(1) = 0    
					XXcouples(it)%mestart(-1) = neven*(neven+1)/2 
					XXcouples(it)%meref(-1) = neven 
				endif 
			endif 
!	 print*,it,XXcouples(it)%mestart(1),XXcouples(it)%mestart(-1)
!	 print*,XXcouples(it)%meref(1),XXcouples(it)%meref(-1)
		endif 
 
		return 
	end subroutine setup_tbme_XXcouples 
!============================================================ 
! 
! set up storage of coupled PP,NN TBMEs 
!  
!  started 7/2010 by CWJ @SDSU 
! 
 
subroutine set_tbme_XXarray(it) 
 
  use sporbit 
!  use coupledmatrixelements 
!  use butil_mod
  implicit none 
 
  integer it  ! species, 1 = P, 2 =N 
 
  integer cpair,dpair,ddpair,ccpair 
  integer pair1,pair2,tmp 
  integer dparity 
  integer iref,istart 
  integer itbme 
  integer a,b,c,d 
  integer ja,jb,jc,jd 
  integer jmin,jmax 
  integer nv 
  integer ncpair 
  type (vjs), pointer :: xxme(:) 
  integer, pointer :: cmap(:) 
  integer :: aerr 
!  type (orb),pointer :: orbqn(:)
 
!.................... COMPUTE # OF MATRIX ELEMENTS 
  nv = 0 
  ncpair = XXcouples(it)%meref(-1) 
  nv = ncpair*(ncpair+1)/2 
  ncpair = ncouplesXX(it) - XXcouples(it)%meref(-1) 
  nv = nv + ncpair*(ncpair+1)/2 
  nv2bmedim(it) = nv 
  if(nv == 0)return 
  if(it == 1)then 
     allocate( ppme( nv2bmedim(it) ), stat=aerr) 
!     if(aerr /= 0) call memerror("set_tbme_XXarray 1") 
     xxme => ppme 
     cmap => PPcouplemap 
!	 orbqn => orbqn
  else 
     allocate( nnme( nv2bmedim(it) ), stat=aerr) 
!     if(aerr /= 0) call memerror("set_tbme_XXarray 2") 
     xxme => nnme 
     cmap => NNcouplemap 
!	 orbqn => norbqn
	 
  endif 
 
  do dpair = 1,ncouplesXX(it) 
 
     dparity = XXcouples(it)%pairc(dpair)%par 
     iref = XXcouples(it)%meref(dparity) 
     istart = XXcouples(it)%mestart(dparity) 
     c = XXcouples(it)%pairc(dpair)%ia 
     d = XXcouples(it)%pairc(dpair)%ib 
     jc = orbqn(c)%j 
     jd = orbqn(d)%j 
 
     do cpair = 1, dpair 
 
          if( XXcouples(it)%pairc(cpair)%par /= dparity )cycle 
          a = XXcouples(it)%pairc(cpair)%ia 
          b = XXcouples(it)%pairc(cpair)%ib 
 
          ja = orbqn(a)%j 
          jb = orbqn(b)%j 
          jmin = MAX( abs(ja-jb), abs(jc-jd) )/2 
          jmax = MIN(ja+jb,jc+jd)/2 
           
          itbme = istart + (dpair-iref)*(dpair-iref-1)/2+cpair-iref 
		  
          xxme(itbme)%jmin = jmin 
          xxme(itbme)%jmax = jmax 
          if(jmax >= jmin) then 
 
              allocate( xxme(itbme)%v(jmin:jmax), stat=aerr )  
!              if(aerr /= 0) call memerror("set_tbme_XXarray 3") 
              xxme(itbme)%v(:) = 0.0 
 
          endif       
     end do  ! cpair 
 
  end do ! dpair 
 
 
  return 
end subroutine set_tbme_XXarray 
 
!============================================================ 

! 
!  subroutine to compute minimal set of coupled pairs for PN TBMEs 
!  uses information from restricted quantum numbers to determine max W of pairs 
!
!  NOTE: If we have particle-hole conjugation allow for all possible pairs
! 
! INPUT:  
!  create : if F then just count; if T then fill in 
! 
! CALLED BY master_readin_pntbmes
 
subroutine setup_tbme_PNcouples(create) 
 
!  use ntuple_info 
!  use W_info 

  implicit none 
  logical create  ! if count or just create 
   
  integer nrawpairs 
  integer a,b 
  integer ncount 
  integer W 
  integer indx 
  integer par 
  integer neven,ieven 
  integer i,j,itmp 
  type (coupair_qn) :: cptemp 
  integer :: aerr 
   
 
  nrawpairs = numorb*numorb
  if(nrawpairs == 0)return 
 
!................. 
 
  if(create)then 
     if( ncouplesPN == 0)return 
     if(.not.allocated(PNcouples%pairc))then
       allocate( PNcouples%pairc( ncouplesPN ),Pncouplemap(ncouplesPN) , stat=aerr) 
       if(aerr /= 0)print*,"setup_tbme_PNcouples 1",aerr
     end if	 
  endif 
 
  ncount = 0 
  do a = 1,numorb 
     do b = 1,numorb

         ncount = ncount + 1 
         if(create)then 
            
            indx = (a-1)*numorb + b 
            PNcouplemap(indx) = ncount 
            PNcouples%pairc(ncount)%ia = a 
            PNcouples%pairc(ncount)%ib = b 
            PNcouples%pairc(ncount)%indx = indx 
            par = orbqn(a)%par*orbqn(b)%par 
			
!            par = porbqn(a)%par*norbqn(b)%par 
!            par = (3-par)/2 
            PNcouples%pairc(ncount)%par = par 
         end if 
 
     end do ! b 
 
  end do ! a 
 
  if(.not. create)then 
    ncouplesPN = ncount 
!    print*,ncouplesPN,' coupled PN pairs kept ',it 
  end if 
 
!---------------- SORT ! on parity only ------------------- 
  if(create)then 
      if(allsameparity)then 
	  
!     if(pallsameparity.and. nallsameparity)then 
         PNcouples%mestart(1) = 0 
         PNcouples%meref(1) = 0 
         PNcouples%mestart(-1) =-1 
         PNcouples%meref(-1) =ncouplesPN 
     else 
!............ CLEVERSORT on parity............ 
!             FIRST count up # of even parities 
         neven = 0 
         do i = 1,ncouplesPN 
            if( PNcouples%pairc(i)%par ==1)neven = neven+1 
         end do 
         if(neven == 0 .or. neven == ncouplesPN )then 
           if(neven == ncouplesPN )then 
              PNcouples%mestart(1) = 0 
              PNcouples%meref(1) = 0 
              PNcouples%mestart(-1) =-1 
              PNcouples%meref(-1) =ncouplesPN 
           else 
              PNcouples%mestart(1) =-1 
              PNcouples%meref(1) = 0 
              PNcouples%mestart(-1) =0 
              PNcouples%meref(-1) =0 
           endif 
         else 
            ieven = 0 
            i = 1 
            j = ncouplesPN 
            do while( ieven /= neven) 
               if( PNcouples%pairc(i)%par == -1)then 
                   do while( PNcouples%pairc(j)%par ==-1 ) 
                       j = j -1 
                       if(j <= i)then 
                           print*,' whoops missed in pn ' 
                           stop 
                       endif 
                   end do 
!----------------------- SWAP 
                   cptemp = PNcouples%pairc(i) 
                   PNcouples%pairc(i) = PNcouples%pairc(j) 
                   PNcouples%pairc(j) = cptemp 
                   PNcouplemap(  PNcouples%pairc(i)%indx )  =i 
                   PNcouplemap(  PNcouples%pairc(j)%indx )  =j 
                    
               endif 
               ieven = ieven + 1 
               i = i + 1 
            end do 
!..................... FIND START, FINISH................ 
            PNcouples%mestart(1) = 0 
            PNcouples%meref(1) = 0    
            PNcouples%mestart(-1) = neven*(neven+1)/2 
            PNcouples%meref(-1) = neven            
 
         endif 
     endif 
 
  endif 
 
 
  return 
end subroutine setup_tbme_PNcouples 
 
!============================================================ 
! 
! set up storage of coupled PN TBMEs 
! 
! CALLED BY master_readin_pntbmes
! 
subroutine set_tbme_PNarray 

  implicit none 
 
  integer cpair,dpair,ccpair,ddpair,pair1,pair2 
  integer dparity 
  integer iref,istart 
  integer itbme 
  integer a,b,c,d 
  integer ja,jb,jc,jd 
  integer jmin,jmax 
  integer nv 
  integer ncpair 
  integer :: aerr 
 
!.................... COMPUTE # OF MATRIX ELEMENTS 
  nv = 0 
  ncpair = PNcouples%meref(-1) 
  nv = ncpair*(ncpair+1)/2 
  ncpair = ncouplesPN - PNcouples%meref(-1) 
  nv = nv + ncpair*(ncpair+1)/2 
  npnvmes = nv 
  if(npnvmes == 0)then 
     return 
  else 
     if(.not.allocated(pnme))allocate( pnme( npnvmes ) , stat=aerr) 
     if(aerr /= 0) print*,"set_tbme_PNarray 1"
  endif 
 
  do dpair = 1,ncouplesPN 
     dparity = PNcouples%pairc(dpair)%par 
     iref = PNcouples%meref(dparity) 
     istart = PNcouples%mestart(dparity) 
     c = PNcouples%pairc(dpair)%ia 
     d = PNcouples%pairc(dpair)%ib 
!     jc = porbqn(c)%j 
!     jd = norbqn(d)%j 
     jc = orbqn(c)%j 
     jd = orbqn(d)%j 
     do cpair = 1, dpair 
          if( PNcouples%pairc(cpair)%par /= dparity )cycle 
          a = PNcouples%pairc(cpair)%ia 
          b = PNcouples%pairc(cpair)%ib 
!          ja = porbqn(a)%j 
!          jb = norbqn(b)%j
          ja = orbqn(a)%j 
          jb = orbqn(b)%j  
          jmin = MAX( abs(ja-jb), abs(jc-jd) )/2 
          jmax = MIN(ja+jb,jc+jd)/2 
          itbme = istart + (dpair-iref)*(dpair-iref-1)/2+cpair-iref 
          pnme(itbme)%jmin = jmin !j, not 2xj
          pnme(itbme)%jmax = jmax 
          if( jmax >= jmin)then 
             allocate( pnme(itbme)%v(jmin:jmax) , stat=aerr) 
             if(aerr /= 0) print*,"set_tbme_PNarray 2"
             pnme(itbme)%v(:) = 0.0 
          endif 
 
     end do  ! cpair 
 
  end do ! dpair 
 
  return 
end subroutine set_tbme_PNarray 

!============================================================
!
! called by
!  master_readin_pntbmes
!
!  reads in p-n part of the interaction,
!  skips pp and nn
!
! file must be in explict PN formalism
!
!
subroutine readin_tbmes(normalizedpn,emptyham)
    use psis,only:numprot,numneut
    implicit none
    logical :: emptyham,normalizedpn
    
    integer :: nme
    integer :: i,m,mmax
    logical :: success
!----------------------------------------------- 
    character*1 ychar  
    character*70 title 
    real:: ax,bx,x,vscale,spscale
    real, allocatable :: spetmp(:) 
    integer :: ia,ib,ic,id,jj,tt
    integer :: a,b,c,d
    real    :: vv,factor
    integer :: L
    integer :: pair1,pair2,tmp
    integer :: ipar,iref,indx,istart
    
    integer :: aerr
	integer :: neutron_offset
	integer :: it,phase
	integer :: na,nb,nc,nd
	
	logical :: autoscale
	integer :: nmetmp
	
    neutron_offset = numorb
    
!.......... PRINT OUT MAPPING................  
  
  write(6,*)' Here is my understanding of the mapping of states '
          write(6,'(a)')' label   N    L    2J    PROTONS  '
          do i = 1,numorb
              write(6,'(2x,4i5,10x,i3)')i,orbqn(i)%nr,orbqn(i)%l, orbqn(i)%j
          end do
          write(6,'(a)')' label   N    L    2J    NEUTRONS  '
          do i = 1,numorb
              write(6,'(2x,4i5,10x,i3)')i+numorb,orbqn(i)%nr,orbqn(i)%l,orbqn(i)%j
          end do
 
!-------------- READ PAST TITLE CARDS --------------------------- 
  success = .false. 
  do while(.not.success) 
     read(1,'(a)')ychar 
     backspace(1) 
     if(ychar /= '#' .and. ychar /= '!')then 
        success = .true. 
     else 
        read(1,'(a)')title 
        write(6,'(a)')title 
     end if 
  end do 


!--------------- CHECK FOR AUTOSCALING----------------

allocate(spetmp(numorb*2), stat=aerr) 
if(aerr /= 0) print*,"readv2bme_xpn 2"
read(1,*,err=1899,end=1899)nme,(spetmp(i),i= 1,MIN(10,numorb*2) )


if(numorb*2 > 10)then 
   do m = 10,numorb*2-1,10 
      mmax = MIN(10+m,numorb*2) 
      read(1,*,err=1899,end=1899)(spetmp(i),i=1+m,mmax) 
   end do 
endif 

autoscale = .false.
if(nme < 0 )then 
	autoscale = .true.
	nme = abs(nme)
	
	
	print*,' '
	print*,' Attempting to autoscale-- do NOT enter scaling '
	print*,' '
	deallocate(spetmp)
	backspace(1)
	allocate(spetmp(numorb*2+3), stat=aerr) 
	if(aerr /= 0) print*,"readv2bme_xpn 2b"
	read(1,*,err=1899,end=1899)nmetmp,(spetmp(i),i= 1,numorb*2+3) 
	if(nme/=abs(nmetmp))then
		print*,' Sorry, I failed to autoscale, please remove - sign from # matrix elements '
		print*,' and enter scaling by hand '
		stop
	end if
	
	
	
end if


!-------------- ENTER SCALING ----------------------------- 


if(autoscale)then
	spscale = 1.0
	ax = spetmp(numorb*2+2)
	bx = spetmp(numorb*2+1)+numprot+numneut
	x = spetmp(numorb*2+3)
!	print*,ax,bx,x
    if(bx==0. .or. x ==0.)then
		print*,' Scaling = 1'
		vscale=1.0
	else

       vscale = (ax/bx)**x 
     end if
else

  print*,' Enter global scaling for spes, A,B,X ( (A/B)^X ) for TBMEs '
  print*,' (If B or X = 0, then scale by A ) ' 
  read*,spscale,ax,bx,x 
  if ( bx.eq.0.0 .or. x.eq.0.0 ) then 
    vscale = ax 
  else 
     vscale = (ax/bx)**x 
  end if  

end if
print*,' Scaling two-body matrix elements by ',vscale,' s.p.e.s by ',spscale
!-------------- READ IN SPEs --------------
  print*,' * * '
  print*,' * * NOTICE: I expect single-particles space with ',numorb*2,' orbits '
  print*,' * * ' 

  
  
  if(allocated(pspe)) deallocate(pspe)
  if(allocated(nspe)) deallocate(nspe)
  allocate(pspe(numorb),nspe(numorb))
  do i = 1,numorb
	  pspe(i)=spetmp(i)*spscale
  end do
  
  do i = 1,numorb
	  nspe(i)=spetmp(i+numorb)*spscale
  end do
 
  deallocate( spetmp ) 
  
  do i = 1,nme 
!---------- ERROR TRAP ADDED July 2011 CWJ ------------------' 
     
      read(1,*,err=1899,end=1899)ia,ib,ic,id,jj,tt,vv 
      
 
!---------- ERROR TRAP ADDED July 2011 CWJ ------------------' 
 
      if( tt /= 0 .and. tt /= 1)then 
            print*,' bad TBME matrix element; T value is bad ',tt
            print*,' matrix element # ',i 
            print*,ia,ib,ic,id,jj,tt,vv 
            print*,' Check that single particle space matches hamiltonian ' 
            print*,' # of single-particle orbits = ',2*numorb !numorb(1) +numorb(2)
            stop 
       end if 
!-------------- READ IN PP --------------------	   
     if( ia > neutron_offset .or. ib > neutron_offset .or. ic > neutron_offset .or. id > neutron_offset)goto 1001 
     a = ia 
     b = ib 
     c = ic  
     d = id 
	 it = 1
!----------- CHECK PARITY ------------------- 
!     L = porbqn(a)%l+ porbqn(b)%l+porbqn(c)%l+porbqn(d)%l 
     L = orbqn(a)%l+ orbqn(b)%l+orbqn(c)%l+orbqn(d)%l 
 
     if( (-1)**(L) ==-1)then 
        print*,' Oops! Error in parity in (pp) matrix element # ',i
		print*,' Most likely this is a mismatch between single-particle orbits ' 
		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,jj,tt,vv 
        write(6,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))') ! & 
!		     porbqn(a)%j,porbqn(a)%l,  porbqn(b)%j,porbqn(b)%l,  porbqn(c)%j,porbqn(c)%l, & 
!			  porbqn(d)%j,porbqn(d)%l
	 
        stop 
     endif 
 
     phase = 1 
 
 
     if(a < b)then 
        na = b 
        b = a 
        a = na 
        phase = (-1)**( JJ+TT+(orbqn(a)%j+orbqn(b)%j)/2)  ! check 
!        phase = (-1)**( JJ+TT+(porbqn(a)%j+porbqn(b)%j)/2)  ! check 
     endif 
 
     if(c < d)then 
        nc = d 
        d = c 
        c = nc 
        phase = phase*(-1)**( JJ+TT+(orbqn(c)%j+orbqn(d)%j)/2)  ! check 
!        phase = phase*(-1)**( JJ+TT+(porbqn(c)%j+porbqn(d)%j)/2)  ! check 

     endif 
 
     if(a < c .or. (a==c .and. b < d))then 
        na = c 
        nb = d 
        c = a 
        d = b 
        a = na 
        b = nb 
     endif 
 
!---------- CONVERT ------------------------- 
 
     pair1 = a*(a-1)/2 + b 
     pair2 = c*(c-1)/2 + d 
     if( PPcouplemap(pair1) == -1) cycle !goto 1001 
     if( PPcouplemap(pair2) == -1) cycle !goto 1001 
     pair1 = PPcouplemap(pair1) 
     pair2 = PPcouplemap(pair2) 
 
     ipar = XXcouples(it)%pairc(pair1)%par 
     if(ipar /= XXcouples(it)%pairc(pair2)%par)then 
        print*,' problem with parity, boss  (PP) ',a,b,c,d 
        stop 
     endif 
     iref = XXcouples(it)%meref(ipar) 
     istart = XXcouples(it)%mestart(ipar) 
     if(pair1 < pair2)then 
         tmp = pair1 
         pair1 = pair2 
         pair2 = tmp 
     endif 
      
     indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart    
!--------------- ERROR TRAP --------------- 
 
     if(jj > ppme(indx)%jmax .or. jj < ppme(indx)%jmin)then 
		 
        print*,' Oops! Error in adding Js in (pp) matrix element # ',i
 		print*,' Most likely this is a mismatch between single-particle orbits ' 
 		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,jj,tt,vv 
        write(6,'(" I think orbital 2xJ are (respectively) ",4(i4,","))') ! & 
! 		     porbqn(a)%j, porbqn(b)%j,  porbqn(c)%j, porbqn(d)%j
		write(6,*)' See log file for additional information '	 
	
         
        stop 
     endif 
     ppme(indx)%v(jj)=ppme(indx)%v(jj)+vv*phase  *vscale
	 cycle    ! here we can only have pp, OR nn, OR pn
1001 continue 
 
!.......................................................................	 
!---------------------------------- NN INTERACTION ------------- 
!.......................................................................	 
    if( ia <= neutron_offset .or.  ic <= neutron_offset )goto 1002
! print*,i,' reading NN '

    if(ic <= neutron_offset .or. id <= neutron_offset)then
		print*,' Some error should not have reached here '
		print*,' NN ?'
		print*,ia,ib,ic,id
		stop
	end if
 
     a = ia -neutron_offset
     b = ib -neutron_offset
     c = ic  -neutron_offset
     d = id -neutron_offset
	 it = 2
!----------- CHECK PARITY ------------------- 
!     L = norbqn(a)%l+ norbqn(b)%l+norbqn(c)%l+norbqn(d)%l 
     L = orbqn(a)%l+ orbqn(b)%l+orbqn(c)%l+orbqn(d)%l 
 
     if( (-1)**(L) ==-1)then 
        print*,' Oops! Error in parity in matrix element # ',i
 		print*,' Most likely this is a mismatch between single-particle orbits ' 
 		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,jj,tt,vv 
         write(6,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))') ! & 
 !		     norbqn(a)%j,norbqn(a)%l,  norbqn(b)%j,norbqn(b)%l,  norbqn(c)%j,norbqn(c)%l, & 
!			  norbqn(d)%j,norbqn(d)%l
	
        stop 
     endif 
 
     phase = 1 
      
     if(a < b)then 
        na = b 
        b = a 
        a = na 
        phase = (-1)**( Jj+Tt+(orbqn(a)%j+orbqn(b)%j)/2)  ! check 
     endif 
 
     if(c < d)then 
        nc = d 
        d = c 
        c = nc 
        phase = phase*(-1)**( Jj+Tt+(orbqn(c)%j+orbqn(d)%j)/2)  ! check 
     endif 
 
     if(a < c .or. (a==c .and. b < d))then 
        na = c 
        nb = d 
        c = a 
        d = b 
        a = na 
        b = nb 
     endif 
 
     pair1 = a*(a-1)/2 + b 
     pair2 = c*(c-1)/2 + d 
     if( NNcouplemap(pair1) == -1) cycle !goto 1002 
     if( NNcouplemap(pair2) == -1) cycle !goto 1002 
     pair1 = NNcouplemap(pair1) 
     pair2 = NNcouplemap(pair2) 
 
     ipar = XXcouples(it)%pairc(pair1)%par 
     if(ipar /= XXcouples(it)%pairc(pair2)%par)then 
        print*,' problem with parity, boss (NN in xpn) ',a,b,c,d 
        stop 
     endif 
     iref = XXcouples(it)%meref(ipar) 
     istart = XXcouples(it)%mestart(ipar) 
      if(pair1 < pair2)then 
         tmp = pair1 
         pair1 = pair2 
         pair2 = tmp 
     endif  
      
     indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart    
!--------------- ERROR TRAP --------------- 
 
     if(jj > nnme(indx)%jmax .or. jj < nnme(indx)%jmin)then  
        print*,' Oops! Error in adding Js in (nn) matrix element # ',i
  		print*,' Most likely this is a mismatch between single-particle orbits ' 
  		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,jj,tt,vv 
        write(6,'(" I think orbital 2xJ are (respectively) ",4(i4,","))') ! & 
!  		     norbqn(a)%j, norbqn(b)%j,  norbqn(c)%j, norbqn(d)%j

        stop 
     endif 
     nnme(indx)%v(jj)=nnme(indx)%v(jj)+vv*phase *vscale
 
!............................................................... 
!     if ( nproc > 1 ) then 
!        n_nn = n_nn + 1 
!        nn_indx_J(1,n_nn) = indx 
!        nn_indx_J(2,n_nn) = j 
!        nn_bcast(n_nn) = nnme(indx)%v(j) 
!     end if 
!............................................................... 
     cycle
1002 continue 
!.......................................................................	 
!--- PN --- PN --- PN  --- PN --- PN ---- PN ----- PN ---- PN ------
!.......................................................................	  
!----------- PUT INTO CORRECT ORDER; PICK UP PHASES ------------- 
!           "CORRECT" ORDER: a >= b, c >= d 
 

       if( ia > numorb .or. ib <= numorb .or. ic >numorb .or. id <= numorb)cycle
       a = ia 
       b = ib - numorb
       c = ic  
       d = id -numorb
     
     
       ! factor due to identical label 
       ! factor should ALREADY be accounted for, if "normalized"
       ! if "unnormalized" must take into account
       !OG test: force flag = false
       !normalizedpn = .false.
       if(normalizedpn)then
             factor=1.0
       else
          print*,'normalizedpn is false'
          factor = 0.5 
          
          if(a == b)factor = factor*sqrt(2.) 
          if(c == d) factor = factor*sqrt(2.) 
       end if     
	
  !     print*,i,a,b,c,d,factor
 
  !----------- CHECK PARITY ------------------- 
       L = orbqn(a)%l+ orbqn(b)%l+orbqn(c)%l+orbqn(d)%l 
 
!       L = porbqn(a)%l+ norbqn(b)%l+porbqn(c)%l+norbqn(d)%l 
 
       if( (-1)**(L) ==-1)then 
          print*,' Oops! Error in parity in (pn 1) matrix element # ',i
           print*,' Most likely this is a mismatch between single-particle orbits ' 
           write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,jj,tt,vv 
           write(6,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))') ! & 
 !               porbqn(a)%j,porbqn(a)%l,  norbqn(b)%j,norbqn(b)%l,  porbqn(c)%j,porbqn(c)%l,  norbqn(d)%j,norbqn(d)%l
             
     
          stop 
       endif 
 
!-------- NO PHASE, ALREADY IN CORRECT ORDER 
!       phase = 1 

       
!       pair1 = numorb(2)*(a-1) + b 
!       pair2 = numorb(2)*(c-1) + d 
       pair1 = numorb*(a-1) + b 
       pair2 = numorb*(c-1) + d 
       if( PNcouplemap(pair1) == -1) cycle
       if( PNcouplemap(pair2) == -1) cycle
       pair1 = PNcouplemap(pair1) 
       pair2 = PNcouplemap(pair2) 
 
       ipar = PNcouples%pairc(pair1)%par 
       if(ipar /= PNcouples%pairc(pair2)%par)then 
          print*,' problem with parity, boss (pn 1)  ' 
		  print*,ipar, PNcouples%pairc(pair2)%par
          stop 
       endif 
       if(pair1 < pair2)then 
           tmp = pair1 
           pair1 = pair2 
           pair2 = tmp 
       endif 
       iref = PNcouples%meref(ipar) 
       istart = PNcouples%mestart(ipar) 
      
       indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart    
 
  !--------------- ERROR TRAP --------------- 
 
       if(jj > pnme(indx)%jmax .or. jj < pnme(indx)%jmin)then 
         
           print*,' Oops! Error in adding Js in (pn1) matrix element # ',i
            print*,' Most likely this is a mismatch between single-particle orbits ' 
            write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,jj,tt,vv 
           write(6,'(" I think orbital 2xJ are (respectively) ",4(i4,","))') ! & 
 !             porbqn(a)%j,  norbqn(b)%j,  porbqn(c)%j, norbqn(d)%j
 
          stop 
       endif 
       
       ! print*,a,b,c,d,jj,vv ,'indx',indx,'factor,vscale',factor,vscale
       
       pnme(indx)%v(jj)=pnme(indx)%v(jj)+vv*factor*vscale
	!print*,pnme(indx)%v(jj)
       if(vv /=0)emptyham = .false.
     
     
 end do
 
 print*,' '
 print*,' All finished reading in '
 print*,' '
    
return
1899 continue
    print*,' Error in reading file; problem may be incommensurate size of single-particle space '
    print*,' I expect single-particles space with ',numorb*2,' orbits '

!    print*,' I expect single-particles space with ',numorb(1)+numorb(2),' orbits '
    print*,nme,i
    backspace(1)
    read(1,'(a)')title
    write(6,*)title
end subroutine readin_tbmes


!============================================================
!
!  a function to get the two-body matrix element with the proper phase
!  for arbitrary ordering of indices
!
!  NOTE DURING DEVELOPMENT: worry about phase!
!
real function twobodyme(a,b,c,d,J,iT)
implicit none

integer :: a,b,c,d   ! orbit labels
integer :: J         ! ang momentum
integer :: iT		! = 1 is PP, =2 is NN, = 0 is PN

integer :: pair1,pair2
integer :: iref,istart,indx
integer :: phase
real    :: vvv
integer :: na,nb,nc,nd,tmp,ipar
type (vjs), pointer :: xme(:)
type (coupairinfo),pointer :: xcouples
integer,pointer :: Xcouplemap(:)	
!type (orb), pointer :: xorbqn(:)
	
 select case (iT)
 
 case (0)      !PN
 
 xme => pnme
 xcouples => PNcouples
 Xcouplemap => PNcouplemap

 case (1)      ! PP
 
 xme => ppme
 xcouples => XXcouples(1)
 Xcouplemap=> PPcouplemap
! xorbqn   => porbqn 
 
 case (2)      !NN
 
 xme => nnme
 xcouples => XXcouples(2)
 Xcouplemap=> NNcouplemap
! xorbqn   => norbqn
 
 case default
 
 twobodyme = -9999.00
 return
	 
 end select
 
 phase = 1
 
 if(iT == 0)then
!       pair1 = numorb(2)*(a-1) + b 
!       pair2 = numorb(2)*(c-1) + d 
       pair1 = numorb*(a-1) + b 
       pair2 = numorb*(c-1) + d 
 else
	 
  if(a < b)then 
    na = b 
    b = a 
    a = na 
    phase = phase*(-1)**( J+1+(orbqn(a)%j+orbqn(b)%j)/2)  ! check 
  endif 

  if(c < d)then 
    nc = d 
    d = c 
    c = nc 
    phase = phase*(-1)**( J+1+(orbqn(c)%j+orbqn(d)%j)/2)  ! check 
  endif 
  pair1 = a*(a-1)/2 + b 
  pair2 = c*(c-1)/2 + d  
 
end if


!print*,pair1,pair2,a,b,c,d
pair1 = Xcouplemap(pair1) 
pair2 = Xcouplemap(pair2) 
if( Xcouplemap(pair1) == -1 .or. Xcouplemap(pair2) == -1) then
	print*,' Some wrong couples '
	stop
end if

ipar = Xcouples%pairc(pair1)%par 
if(ipar /= Xcouples%pairc(pair2)%par)then 
   print*,' problem with parity, boss  twobodyme  ',it
  print*,ipar, Xcouples%pairc(pair2)%par
  print*,pair1,pair2
  print*,a,b,c,d
    
   stop 
endif 
if(pair1 < pair2)then 
    tmp = pair1 
    pair1 = pair2 
    pair2 = tmp 
endif 
iref = Xcouples%meref(ipar) 
istart = Xcouples%mestart(ipar) 

indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart    

twobodyme = xme(indx)%v(j)*phase

return

end function twobodyme
!=========================================================================
!
!  routine to get the index and phase 
!  for arbitrary ordering of indices
!
!  NOTE DURING DEVELOPMENT: worry about phase!
!
subroutine twobodyindexphase(a,b,c,d,J,iT,indx,phase)
implicit none

integer :: a,b,c,d   ! orbit labels
integer :: J         ! ang momentum
integer :: iT		! = 1 is PP, =2 is NN, = 0 is PN

integer :: pair1,pair2
integer :: iref,istart,indx
integer :: phase
real    :: vvv
integer :: na,nb,nc,nd,tmp,ipar
type (coupairinfo),pointer :: xcouples
integer,pointer :: Xcouplemap(:)	
!type (orb), pointer :: xorbqn(:)

    phase = 1
	indx = -1
	
 select case (iT)
 
 case (0)      !PN
 
 xcouples => PNcouples
 Xcouplemap => PNcouplemap

 case (1)      ! PP
 
 xcouples => XXcouples(1)
 Xcouplemap=> PPcouplemap
! xorbqn   => porbqn 
 
 case (2)      !NN

 xcouples => XXcouples(2)
 Xcouplemap=> NNcouplemap
! xorbqn   => norbqn
 
 case default
 
 indx = -1
 return
	 
 end select
 
 phase = 1
 
 if(iT == 0)then
!       pair1 = numorb(2)*(a-1) + b 
!       pair2 = numorb(2)*(c-1) + d 
       pair1 = numorb*(a-1) + b 
       pair2 = numorb*(c-1) + d 
 else
	 
  if(a < b)then 
    na = b 
    b = a 
    a = na 
    phase = phase*(-1)**( J+1+(orbqn(a)%j+orbqn(b)%j)/2)  ! check 
  endif 

  if(c < d)then 
    nc = d 
    d = c 
    c = nc 
    phase = phase*(-1)**( J+1+(orbqn(c)%j+orbqn(d)%j)/2)  ! check 
  endif 
  
  pair1 = a*(a-1)/2 + b 
  pair2 = c*(c-1)/2 + d  
end if

if( Xcouplemap(pair1) == -1 .or. Xcouplemap(pair2) == -1) then
	print*,' Some wrong couples '
	stop
end if
pair1 = Xcouplemap(pair1) 
pair2 = Xcouplemap(pair2) 

ipar = Xcouples%pairc(pair1)%par 
if(ipar /= Xcouples%pairc(pair2)%par)then 
   print*,' problem with parity, boss (pn twobody indexphase)  ' 
   stop 
endif 
if(pair1 < pair2)then 
    tmp = pair1 
    pair1 = pair2 
    pair2 = tmp 
endif 
iref = Xcouples%meref(ipar) 
istart = Xcouples%mestart(ipar) 

indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart    


return

end subroutine twobodyindexphase
!=========================================================================    
    
end module coupledinteraction
