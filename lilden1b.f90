module dens1b
	
	use little
	use sporbit
	implicit none
	
    logical :: binary_dens = .false.  ! if true, write out as binary file
	
	
	type rhobase
		complex(kind=8), allocatable :: prho(:,:),nrho(:,:)
		complex(kind=8), allocatable :: Pprho(:,:),Pnrho(:,:)   ! parity-inverted on initial
		
	end type rhobase	
	type rhobaseJt
		integer :: Jtmin,Jtmax
		type (rhobase), allocatable :: rhoJt(:)
	end type rhobaseJt
	
	
	type rhobaseJiJf
		type(rhobase), allocatable :: mat(:,:)
		
	end type rhobaseJiJf

		

!  rhotr(sdf,sdi,Jf,Ji,Kf,Ki)	
	type (rhobaseJt), allocatable,target  :: rhotr(:,:,:,:,:,:)	
	type (rhobaseJiJf), allocatable,target ::  rhofinal(:,:,:,:) 	
		

	   
   logical :: finddense	   
   
contains
	
!================================================================

subroutine readin_densitymat
	implicit none
	
	
	character*80 ::filename
	integer :: ilast
	integer :: a,b
	integer :: iJ,fJ,iMi,iMf
	integer :: iJshift,fJshift
	real :: xJi,xJf,xMi,xMf
	integer :: Jtmin,Jtmax,Jt
	integer :: sdx,sdy
	integer :: N,Z,numsdused
	integer :: n1,n2,n3,n4,n5,n6,n7
	logical :: test
	
	integer :: numorbdum,nspsdum
	
	integer :: i,Jspmax
	logical :: density_file
	density_file=.false.
	print*,' Enter name of .denmat file (do not include extension)'
	
	read(5,'(a)')filename
	ilast = index(filename,' ')-1
	open(unit=33,file=filename(1:ilast)//".denmat",status='unknown',err=111)
	
	read(33,*,err=111)numorbdum,nspsdum
	goto 222

111 continue

    print*,' May be unformattted (binary); trying that '
	close(33)
	open(unit=33,file=filename(1:ilast)//".denmat",status='unknown',form='unformatted')
	read(33)numorbdum,nspsdum
    density_file=.true.

222 continue	
	
	if(numorbdum /= numorb)then
		print*,' Mimatch on # of orbits '
		print*,numorbdum,' expect ',numorb
		stop
		
	end if
	if(nspsdum /= nsps)then
		print*,' Mimatch on # of s.p. states '
		print*,nspsdum,' expect ',nsps
		stop
		
	end if
	if(density_file)then
		read(33)Z,N,numsdused
		
	else
		read(33,*)Z,N,numsdused
		
	endif
	
	if(Z/=numprot)then
		print*,'Mismatch on Z',z, numprot
		stop
		
	end if
	if(N/=numneut)then
		print*,'Mismatch on N',N, numneut
		stop
		
	end if
	if(numsd/=numsdused)then
		print*,'Mismatch on # SDs',numsdused,numsd
		stop
		
	end if
!	print*,numOfJ
	  bigJmax = nint(2*xJlist(numOfJ)+1)
	  
	  if(density_file)then
	  	read(33)n1,n2,test
		  
	  else
	  	read(33,*)n1,n2,test
		  
	  end if
	  
	
	if(n1 /=bigJmax)then
		print*,' Mismatch bigJmax ',n1,bigJmax
		print*,' You may not have set your maximum J consistently  in generating your files '
		stop
	end if
	
	if(n2 /= numOfJ)then
		print*,' Mismatch numOfJ ',n2,numOfJ
		stop
	end if
	
	if(test .eqv. Paritytest)then
		print*,' Mismatch parity flag ',test,ParityTest
		stop
		
	end if
			
!.......... READ IN ORBITS.......

    allocate(orbqn(numorb))
	allsameparity=.true.
	do i = 1,numorb
		if(density_file)then
			read(33)orbqn(i)%j,orbqn(i)%l,orbqn(i)%nr
			
		else
			read(33,*)orbqn(i)%j,orbqn(i)%l,orbqn(i)%nr
			
		end if
		orbqn(i)%par = (-1)**orbqn(i)%l
		if(i > 1)then
		  if(orbqn(i)%par/= orbqn(1)%par)allsameparity=.false.
	    end if
		
	end do	
	print*,' '
	
	Jspmax = 0
	do i = 1,numorb
		Jspmax = max(Jspmax,orbqn(i)%j)
	end do
		
	allocate(rhotr(numSD,numSD,numOfJ,bigJmax,numofJ,bigJmax))
	
	rhotr%Jtmin = -999
	rhotr%Jtmax = -999

	Jmax = 0.5*(bigJmax-1)
		
    do iJ = 1, numOfJ
		
        xJi = xJlist(iJ)
		iJshift = nint(Jmax-xJi)
		do fJ = 1,numOfJ
			xJf = xJList(fJ)
			fJshift = nint(Jmax-xJf)
			
			do iMi = 1,bigJmax ! this is M_beta in my notes
				xMi = float(iMi)-(Jmax+1.) ! this is K_mu (associate with beta= initial ) in my notes
				if(xJi < abs(xMi))cycle
				
				do iMf = 1, bigJmax ! this is M_alpha in my notes
					xMf = float(iMf)-(Jmax+1.) ! this is K_lambda (associated with alpha=final) in my notes
					if(xJf < abs(xMf))cycle
				
					Jtmin = nint(abs(xJi-xJf))
					Jtmax = nint(xJi+xjF)!rhotr(1,1,fJ,iMf-fJshift,iJ,iMi-iJshift)%Jtmax
					Jtmax = min(Jtmax,Jspmax)
					
					if(Jtmin > Jtmax)cycle
					if(Jtmin==-999 .or. Jtmax ==-999)cycle
										do sdX = 1,numsd
						do sdy = 1,numsd
							rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%Jtmin=Jtmin
							rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%Jtmax=Jtmax

							allocate(rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jtmin:Jtmax) )
							
							do Jt = Jtmin,Jtmax
								if(numprot> 0)allocate(rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%prho(numorb,numorb))
								if(numneut> 0)allocate(rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%nrho(numorb,numorb))
								if(.not.allsameparity)then
									if(numprot> 0)allocate(rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%Pprho(numorb,numorb))
									if(numneut> 0)allocate(rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%Pnrho(numorb,numorb))									
									
								end if
								
							end do
						end do
					end do

					do Jt = Jtmin,Jtmax
						do sdX = 1,numsd
							do sdY = 1,numsd
								if(density_file)then
									read(33)n1,n2,n3,n4,n5,n6,n7 !sdx,sdy,fJ,iMf,iJ,iMf,Jt
									
								else
									read(33,*)n1,n2,n3,n4,n5,n6,n7 !sdx,sdy,fJ,iMf,iJ,iMf,Jt
									
								end if

								if(n1/=sdx .or. n2 /= sdy .or. n3 /=fj .or. n5 /=iJ)then
									print*,' Mismatch (1) ',n1,sdx,n2,sdy,n3,fj,n5,ij
									print*,xji,xmi,xjf,xmf
									stop
								end if
								if(n4/=iMf .or. n6/= iMi .or. n7 /=Jt)then
									print*,' Mismatch (2) ',n4,iMf,n6,iMi,n7,Jt
									stop
								end if
								
								if(numprot> 0)then
									if(density_file)then
										read(33)rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%prho(:,:)
									
									else
										read(33,*)rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%prho(:,:)
									
									end if


								end if
								if(numneut> 0)then
									if(density_file)then
										read(33)rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%nrho(:,:)
									
									else
										read(33,*)rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%nrho(:,:)
									
									end if
								end if
								if(.not.allsameparity)then
								
									if(numprot> 0)then
										if(density_file)then
											read(33)rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%Pprho(:,:)
									
										else
											read(33,*)rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%Pprho(:,:)
									
										end if


									end if
									if(numneut> 0)then
										if(density_file)then
											read(33)rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%Pnrho(:,:)
									
										else
											read(33,*)rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%Pnrho(:,:)
									
										end if
									end if

								end if
								
								
							end do ! sdY
						end do ! sdZ

					end do ! Jt
				end do ! imF
			end do ! iMi
		end do ! fJ
	end do ! iJ	
			


	close(33)
	
	

end subroutine readin_densitymat	
	
	
	
end module dens1b