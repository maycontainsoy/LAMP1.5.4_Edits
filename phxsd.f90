!
!  reads in HF basis file from SHERPA and
!  creates particle-hole excited Slater deteminants
!
      module spspace
         implicit none      
!-------------SINGLE-PARTICLE STATES---------------------------   

         integer norb(2)		! # of single-particle j-orbits
         				! 1 = p,  2 = n
         integer nsps(2)		! # of single-particle m-states
      		   		        ! 1 = p, 2 = n
      
         integer,allocatable  :: jorb(:)	! 2 x J of orbit (note: p,n different)
         integer,allocatable  :: nrorb(:)	! radial quantum number N of j-orbit
         integer,allocatable  :: torb(:)	! 2 x Tz of orbit ( p/n = +/-1 )
      
         real,allocatable     :: orb_qn(:,:)	! quantum #'s of j-orbits  -- from WE0
      				! note: p,n have different orbits
      				! orb_qn(i,1) = N of ith orbit
      				! orb_qn(i,2) = L of ith orbit
      				! orb_qn(i,3) = J of ith orbit

         integer,allocatable  :: spsqn(:,:,:)
      				! quantum #'s of s.p. m-states - from WEO
      				! spsqnp(it,i,j)  it = 1,2 proton/neutron  
      				! i = label of s.p. m-state
      				! j = 1 -> label of j-orbit (note p,n different)
      				! j = 2 -> radial quantum number N
      				! j = 3 -> L 
      				! j = 4 -> 2 x J of state
      				! j = 5 -> 2 x M of state
      				! j = 6 -> 2 x Tz of state ( p/n=+/-1 )
        logical spinless
		
	contains 
		
! reads from sps information and reconstruct orbital information
		
		subroutine extract_orbital_info
			implicit none
			
			integer iorb
			integer :: isps,jsps
			
			integer :: it
			
			
			
!			do it = 1,2
             it = 1
			 
			 norb(it)=1
			 do isps = 1,nsps(it)
				 norb(it)=max(norb(it),spsqn(it,isps,1))
			 end do
             print*,' space contains ',norb(1),' unique orbits '
             norb(2)=norb(1)
			 
			 allocate(nrorb(norb(it)),jorb(norb(it)))
			 nrorb=-1
			 jorb =-2
				
			do isps = 1,nsps(it)
				iorb = spsqn(it,isps,1)
!				print*,isps,iorb,jorb(iorb)
				if( jorb(iorb) < 0 )then
					nrorb(iorb)=spsqn(it,isps,2)
					jorb(iorb)=spsqn(it,isps,4)
				end if

			end do

			
			return
			
		end subroutine extract_orbital_info
		
      end module spspace
!===================================================================
	  module nuclide
	     implicit none
	     integer                          :: np,nn   ! # of protons/neutrons
	  end module nuclide

!===================================================================

	        module wfn
	           implicit none
	           real,allocatable,dimension(:,:)  :: pSD,nSD ! proton/neutron SD's      
	            				         ! columns = particles, rows = m-states
	  		 integer tempfile			! location of temporary file for storing input sds
			 
			 logical writesdout    !  write out sd
								 
	  	end module wfn
	   
	  
module hfbasis
   implicit none
   real, allocatable,target :: ehfp(:),ehfn(:)
   real,allocatable,target :: basishfp(:,:),basishfn(:,:)
   
   real :: avgDE,maxE,minE
   integer :: nphstates

   contains

subroutine read_in_hf_basis
	use spspace
	use nuclide
   implicit none
   integer :: isps
   integer :: i,ii
   real    :: ex

   open(unit=1,file='hf.basis',status='old',err=111)
   
   read(1,*)nsps(1)
   nsps(2)=nsps(1)
   
   allocate(spsqn(2,nsps(1),6))
   spsqn=0
   do i = 1,nsps(1)
	   read(1,*)ii,spsqn(1,i,1),spsqn(1,i,2),spsqn(1,i,3),spsqn(1,i,4), & 
	      spsqn(1,i,5)
!	   print*,ii,spsqn(1,i,2),spsqn(1,i,3),spsqn(1,i,4),spsqn(1,i,5)
	   if(ii/=i)then
		   print*,' Error in reading sps q#s ',i,ii
		   stop
	   end if
	   do ii = 1,5
		   spsqn(2,i,ii)=spsqn(1,i,ii)
	   end do
   end do
   print*,' there are ',nsps(1),' single-particle states '
   read(1,*)np,nn
   print*,' Valence Z, N  = ',np,nn

   call extract_orbital_info

!...... READ IN HF ENERGIES ......................
   
   allocate(ehfp(nsps(1)),ehfn(nsps(2)))
   
   do i = 1,nsps(1)
	   read(1,*)ii,ehfp(i)
	   if(ii/=i)then
		   print*,' Error in reading proton E ',i,ii
		   stop
	   end if   
   end do
   do i = 1,nsps(2)
	   read(1,*)ii,ehfn(i)
	   if(ii/=i)then
		   print*,' Error in reading neutron E ',i,ii
		   stop
	   end if   
   end do
!.............. READ IN HF BASIS VECTORS..................

allocate(basishfp(nsps(1),nsps(1)))   
do i = 1,nsps(1)
	read(1,*)ii,ex
    if(ii/=i)then
	   print*,' Error in reading proton sp wf ',i,ii
	   stop
    end if   	
	if(abs(ex-ehfp(i))> 0.00001)then
		print*,' Error in comparing proton sp energies ',ex,ehfp(i)
		stop
		
	end if
	read(1,*)(basishfp(ii,i),ii=1,nsps(1))
	
	
end do
   
allocate(basishfn(nsps(2),nsps(2)))   
do i = 1,nsps(2)
	read(1,*)ii,ex
    if(ii/=i)then
	   print*,' Error in reading neutron sp wf ',i,ii
	   stop
    end if   	
	if(abs(ex-ehfn(i))> 0.00001)then
		print*,' Error in comparing neutron sp energies ',ex,ehfn(i)
		stop
		
	end if
	read(1,*)(basishfn(ii,i),ii=1,nsps(2))
	
	
end do
close(1)
print*,' file read !'
   return
111 continue
print*,' file hf.basis does not seem to exist '
stop

end subroutine read_in_hf_basis

!
!  plist and nlist are lists of occupied proton, neutron states in the hf basis
!  subroutine default_selection fills the lowest energies
!

subroutine default_selection(np,nn,plist,nlist)
	use spspace
	implicit none
	integer,intent(in) :: np,nn ! Z,N valence
	integer :: plist(np),nlist(nn)
	integer i,ip,isps
	
	if(np > 0)then
		do i = 1,np
			plist(i)=i
		end do
	end if
	if(nn > 0)then
		do i = 1,nn
			nlist(i)=i
		end do
	end if
	return
end subroutine default_selection

!  plist and nlist are lists of occupied proton, neutron states in the hf basis
! subroutine fill_selected_states fills the Slater determinant arrays psd, nsd
! with occupied proton, neutron s.p. wavefunctions from the hf basis
!

subroutine fill_selected_states(np,nn,plist,nlist)
	use wfn
	use spspace
	implicit none
	integer,intent(in) :: np,nn ! Z,N valence
	integer :: plist(np),nlist(nn)
	integer i,ip,isps
	
	if(np > 0)then
		do i= 1,np
			ip = plist(i)
			do isps  = 1,nsps(1)
				psd(isps,i)=basishfp(isps,ip)
			end do
		end do		
	end if
	if(nn > 0)then
		do i= 1,nn
			ip = nlist(i)
			do isps  = 1,nsps(2)
				nsd(isps,i)=basishfn(isps,ip)
			end do
		end do		
	end if
	return
	
end subroutine fill_selected_states

end module hfbasis


!=========================================================
      subroutine save_sd
!
!  saves slater determinant wfn to a file
!  also saves s.p. information, plus a header
!  has option to write out multiple SDs from tempfile
!

!  SUBROUTINES CALLED
!	writeoutsd
!	readinsd
!
      use spspace
      use wfn
      use nuclide

      implicit none

      real,allocatable           :: sdtemp(:,:)
!..............FILE HANDLING..........................

      character*1 ychar
      character filename*45  		! 
      integer ilast
      logical append			! controls appending
      logical errflag

!..............MISC........................

      integer i,ii,j,n,z
      character title*60

311    continue
      write(6,*)' Enter output filename (.sd)'
	  write(6,*)' (Enter "none" to only compute statistics )'
      read(5,'(a)')filename
      ilast=index(filename,' ')
      if(ilast.ne.0)then
         ilast=ilast-1
      else
         ilast=15
      endif
	  
	  if(filename(1:4)=='none')then
		  writesdout= .false.
		  return
	  else
		  writesdout=.true.
	  end if

      append = .false.
      open(unit=2,file=filename(1:ilast)//'.sd',status='new',err=33,  form='unformatted')
      goto 44
33    continue
      write(6,*)' That file already exists; overwrite (o) or append (a) or new (n) file?(o/n/a) '
!
! IF A FILE ALREADY EXISTS YOU CAN OVERWRITE OR APPEND

      read(5,'(a)')ychar
      if(ychar.eq.'O' .or. ychar.eq.'o')then
	open(unit=2,file=filename(1:ilast)//'.sd',status='unknown' ,  form='unformatted')
      elseif(ychar.eq.'A' .or. ychar.eq.'a')then
        append = .true.
	open(unit=2,file=filename(1:ilast)//'.sd',status='old',  form='unformatted')
      else
	goto 311
      endif
44    continue

!
! AS LONG AS NOT APPENDING, WRITE DOWN S.P. INFORMATION
! ASSUME PROTON, NEUTRON SPACES THE SAME
!
      if(.not.append)then
      	write(6,*)'(I assume the proton, neutron spaces the same)'
	write(2)norb(1)
	do i =1,norb(1)
	   write(2)i,nrorb(i),jorb(i)
        enddo
	write(2)np,nn
        write(6,*)' Enter a comment '
        read(5,'(60a)')title
        write(2)title
      else
!
!  IF APPEND, THEN CHECK TO MAKE SURE INFO MATCHES
!
 	read(2)n
	if(n.ne.norb(1))then
	    write(6,*)' inconsistent orbit info ',  norb(1),n
	    close(unit=2)
	    goto 311	    
        endif
	do i =1,norb(1)
	   read(2)ii,n,j
	   if(n.ne.nrorb(i).or. j.ne.jorb(i))then
		write(6,*)' inconsistent orbits ', i,nrorb(i),n,jorb(i),j
		close(unit=2)
		goto 311
	   endif
        enddo
	read(2)z,n
        if(z.ne.np .or. n.ne.nn)then
		write(6,*)' wrong nuclide ',np,z,nn,n
		close(unit=2)
		goto 311
	endif
        read(2)title
        write(6,*)title
!............ now readout to dummy array.............
	do
          allocate(sdtemp(nsps(1),np))
	  call readinsd(nsps(1),np,2,sdtemp,errflag)
          deallocate(sdtemp)
	  if(errflag)goto 443
          allocate(sdtemp(nsps(1),nn))
	  call readinsd(nsps(2),nn,2,sdtemp,errflag)
          deallocate(sdtemp)
	  if(errflag)goto 443
        enddo
443     continue

      endif
      if(tempfile.eq.0)then
	call writeoutsd(nsps(1),np,2,psd)
	call writeoutsd(nsps(2),nn,2,nsd)
      else
        rewind(tempfile)
	do i = 1,1000
	  call readinsd(nsps(1),np,tempfile,psd,errflag)
	  if(errflag)goto 444
 	  call writeoutsd(nsps(1),np,2,psd)
	  call readinsd(nsps(2),nn,tempfile,nsd,errflag)
	  if(errflag)goto 444
 	  call writeoutsd(nsps(2),nn,2,nsd)

        enddo
444     continue
      endif

!      close(unit=2)
      return
      end
	  
!======================================================================

      subroutine readinsd(nsps,n,ifile,psi,errflag)

      implicit none

      integer,intent(IN)      :: nsps,n
      real,intent(OUT)        :: psi(nsps,n) ! slater determinant
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

      end   
!
!=========================================================
!
      subroutine writeoutsd(nsps,n,ifile,psi) 
		  
		  use wfn,only:writesdout

      implicit none
      
      integer,intent(IN)   :: nsps,n
      real,intent(IN)      :: psi(nsps,n) ! slater determinant
      integer i,j
      integer ifile
      
	  if(.not.writesdout)return
      do i = 1,n
        write(ifile)(psi(j,i),j=1,nsps)
      enddo
100   format(20(F12.6,1x))
      return
      end
      
!======================================================================
!
!  MAIN PROGRAM

   
use hfbasis

use nuclide
use wfn
use spspace
implicit none
integer,allocatable :: plist(:),nlist(:)
real :: ecut

integer :: nphchoice

call read_in_hf_basis

if(np > 0)then
	allocate(plist(np))
	allocate(psd(nsps(1),np))
end if
if(nn > 0)then
	allocate(nlist(nn))
	allocate(nsd(nsps(2),nn))
end if

call default_selection(np,nn,plist,nlist)
call fill_selected_states(np,nn,plist,nlist)

call save_sd

avgDE = 0.0
maxE = 0.0
minE = 100000.0
nphstates = 0

print*,' Do you want 1p-1h or 2p-2h (enter 1/2, or 3 for both)'
read*,nphchoice

print*,' Enter cutoff if desired '
print*,' (enter  < 0 if no cut off )'
read*,ecut

select case (nphchoice)

case(1)
print*,' Proton particle-hole states '
call select_1p1h_boss(1,np,nn,plist,nlist,ecut)
print*,' Neutron particle-hole states '

call select_1p1h_boss(2,np,nn,plist,nlist,ecut)

case(2)
call select_2p2h_boss(1,np,nn,plist,nlist,ecut)
call select_2p2h_boss(2,np,nn,plist,nlist,ecut)
call select_pn_2p2h_boss(np,nn,plist,nlist,ecut)

case (3)
call select_1p1h_boss(1,np,nn,plist,nlist,ecut)

call select_1p1h_boss(2,np,nn,plist,nlist,ecut)
call select_2p2h_boss(1,np,nn,plist,nlist,ecut)
call select_2p2h_boss(2,np,nn,plist,nlist,ecut)
call select_pn_2p2h_boss(np,nn,plist,nlist,ecut)

case default

print*,' bad choice ',nphchoice
stop

end select

print*,' '
print*,' Avg excitation = ',avgDE/float(nphstates),' out of ',nphstates,' states '
print*,' Min, max excitation = ',minE,maxE

end
   
!=========================================================

subroutine select_1p1h_boss(it,np,nn,plist,nlist,ecut)
	
	use spspace
	use hfbasis
	use wfn
	implicit none
	integer,intent(in) :: it  ! species
	integer,intent(in) :: np,nn
	integer,intent(out) :: plist(np),nlist(nn)
	integer :: nph   ! # of particle-hole
	integer :: npart,nhole
	integer :: ipart,jhole
	real    :: eph
	real    :: ecut
	
	
	if(it==1)then
		nhole = np
		npart = nsps(it)-np
	else
		nhole = nn
		npart = nsps(it)-nn		
		
	end if
	
	nph = npart*nhole
	if(nph < 1)return
	if(it==1)print*,nph,' proton  1 particle- 1 hole states '
	if(it==2)print*,nph,' neutron 1 particle- 1 hole states '
	
	do ipart = 1,npart   ! loop over unoccupied states  (note, technically these are called "particle" states)
		do jhole = 1,nhole   ! loop over occupied states (note, technically called "hole" states)
			call default_selection(np,nn,plist,nlist)
			if(it==1)then
				plist(jhole)=ipart+nhole
				eph = ehfp(ipart+nhole)-ehfp(jhole)
!				print*,ehfp(ipart+nhole),ehfp(jhole)
			else
				nlist(jhole)=ipart+nhole
				eph = ehfn(ipart+nhole)-ehfn(jhole)
!				print*,ehfn(ipart+nhole),ehfn(jhole)
			
			end if
			if(ecut > 0. .and. eph > ecut )cycle
			
!			print*,eph,' particle-hole energy ',ipart+nhole,jhole
			nphstates = nphstates+1
			avgDe = avgDE+ eph
			maxe = max(maxE,eph)
			mine = min(mine,eph)
			call fill_selected_states(np,nn,plist,nlist)
			call writeoutsd(nsps(1),np,2,psd)
			call writeoutsd(nsps(2),nn,2,nsd)
		

		end do
	end do
	return


end subroutine select_1p1h_boss

!=========================================================


! for PP or NN excitations
subroutine select_2p2h_boss(it,np,nn,plist,nlist,ecut)
	
	use spspace
	use hfbasis
	use wfn
	implicit none
	integer,intent(in) :: it  ! species
	integer,intent(in) :: np,nn
	integer,intent(out) :: plist(np),nlist(nn)
	integer :: nph   ! # of particle-hole
	integer :: npart,nhole
	integer :: ipart,jpart,ihole,jhole
	real    :: eph
	real    :: ecut
	
	
	if(it==1)then
		nhole = np
		npart = nsps(it)-np
	else
		nhole = nn
		npart = nsps(it)-nn		
		
	end if
	
	nph = npart*nhole
	if(nph < 1)return
	if(it==1)print*,nph,' proton  2 particle-2 hole states '
	if(it==2)print*,nph,' neutron 2 particle-2 hole states '
	
	do ipart = 1,npart-1   ! loop over unoccupied states  (note, technically these are called "particle" states)
		do jpart = ipart+1,npart
			
		do ihole = 1,nhole-1   ! loop over occupied states (note, technically called "hole" states)
			do jhole = ihole+1,nhole
			call default_selection(np,nn,plist,nlist)
			if(it==1)then
				plist(ihole)=ipart+nhole
				plist(jhole)=jpart+nhole
				eph = ehfp(ipart+nhole)+ehfp(jpart+nhole)-ehfp(jhole)-ehfp(ihole)
!				print*,ehfp(ipart+nhole),ehfp(jhole)
			else
				nlist(ihole)=ipart+nhole
				nlist(jhole)=jpart+nhole
				eph = ehfn(ipart+nhole)+ehfn(jpart+nhole)-ehfn(jhole)-ehfn(ihole)
!				print*,ehfn(ipart+nhole),ehfn(jhole)
			
			end if
			if(ecut > 0. .and. eph > ecut )cycle
			
!			print*,eph,' particle-hole energy ',ipart+nhole,jhole
			nphstates = nphstates+1
			avgDe = avgDE+ eph
			maxe = max(maxE,eph)
			mine = min(mine,eph)
			call fill_selected_states(np,nn,plist,nlist)
			call writeoutsd(nsps(1),np,2,psd)
			call writeoutsd(nsps(2),nn,2,nsd)
		end do
	end do
		

		end do
	end do
	return


end subroutine select_2p2h_boss

subroutine select_pn_2p2h_boss(np,nn,plist,nlist,ecut)
	
	use spspace
	use hfbasis
	use wfn
	implicit none
	integer,intent(in) :: np,nn
	integer,intent(out) :: plist(np),nlist(nn)
	integer :: nph   ! # of particle-hole
	integer :: npartp,nholep,npartn,nholen
	integer :: ipartp,jpartn,iholep,jholen
	real    :: eph
	real    :: ecut
	
	
		nholep = np
		npartp = nsps(1)-np
		nholen = nn
		npartn = nsps(2)-nn		
		
	
	nph = npartp*nholep*npartn*nholen
	if(nph < 1)return
	print*,nph,' proton-neutron 2 particle-2 hole states '
	
	do ipartp = 1,npartp   ! loop over unoccupied states  (note, technically these are called "particle" states)
		do jpartn = 1,npartn
			
		do iholep = 1,nholep   ! loop over occupied states (note, technically called "hole" states)
			do jholen = 1,nholen
			call default_selection(np,nn,plist,nlist)
				plist(iholep)=ipartp+nholep
				nlist(jholen)=jpartn+nholen
				eph = ehfp(ipartp+nholep)+ehfn(jpartn+nholen)-ehfn(jholen)-ehfp(iholep)
!				print*,ehfp(ipart+nhole),ehfp(jhole)

			if(ecut > 0. .and. eph > ecut )cycle
			
!			print*,eph,' particle-hole energy ',ipart+nhole,jhole
			nphstates = nphstates+1
			avgDe = avgDE+ eph
			maxe = max(maxE,eph)
			mine = min(mine,eph)
			call fill_selected_states(np,nn,plist,nlist)
			call writeoutsd(nsps(1),np,2,psd)
			call writeoutsd(nsps(2),nn,2,nsd)
		end do
	end do
		

		end do
	end do
	return


end subroutine select_pn_2p2h_boss






	  