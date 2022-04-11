!
!  revision projection in 1.4.7
!

module project1bdense
!	use sampler
	use wignersfriend
	
	implicit none
	
	logical :: binary_dens = .true.  ! if true, write out as binary file
	
contains


!==================================================================================
!
! transforms N_ijk -> N_jMK etc
!
! here N_jMK = sum Z_mi Z_Kk N_ijk  etc
!
!  CALLED BY: projectorator_dens
!
subroutine invert_alpha_and_gamma4dense
!	use system_parameters
	use dens1body
	use spstate !,only:nsps
	use lamplight
 	implicit none

 	integer ialp, kgam, iM,iK,J, indexerDim 
 	real(kind=8):: xM, xK

 	REAL :: timeStart, timeStop ! Timing variables
	
	complex(kind=8) :: ptmp(nsps,nsps),ntmp(nsps,nsps),Pptmp(nsps,nsps),Pntmp(nsps,nsps)
	
	integer :: sdX,sdY   ! indices for multiple SDs	
	
	integer, allocatable :: indexer(:,:)
	integer :: i
	integer :: i1,i2
	
	
	allocate(rho1b_jMK(numsd,numsd,numOfJ,J2max1,J2max1))
	do sdX = 1,numsd
	   do sdY = 1,numsd
 	      do iM = 1,J2max1		
 		     do iK = 1,J2max1
			    do J = 1,numofJ
	!				if(sdX ==2 .and. sdY == 1 .and. j ==9 .and. iM ==1)print*,sdX,sdY,j,iM,iK,' huh '
					
			  	   if(numprot> 0)then
					   allocate(rho1b_jMK(sdX,sdY,J,iM,IK)%prho(nsps,nsps))
				   end if
			 	   if(numneut>0) then
					   allocate(rho1b_jMK(sdX,sdY,J,iM,IK)%nrho(nsps,nsps))
				   end if
				   if(.not.allsameparity)then
				  	   if(numprot> 0)then
						   allocate(rho1b_jMK(sdX,sdY,J,iM,IK)%Pprho(nsps,nsps))
						   
					   end if
				 	   if(numneut>0) then
						   allocate(rho1b_jMK(sdX,sdY,J,iM,IK)%Pnrho(nsps,nsps))
						   
					   end if					   					   
				   end if
				   
			    end do
		     end do
	     end do
	 end do
   end do



!...... SET UP UNROLLED LOOP
indexerDim = (J2Max1**2)*(numsd**2)*(numofJ)
allocate(indexer(indexerDim,5))
i = 0
DO sdX = 1, numsd
	DO sdY = 1, numsd 
		DO iM = 1,J2max1		 
			DO iK = 1,J2max1
				DO j = 1,numofJ
					i = i + 1
					indexer(i,1) = sdX 
					indexer(i,2) = sdY 
					indexer(i,3) = im
					indexer(i,4) = ik
					indexer(i,5) = j
				END DO ! j 
			END DO ! iK
		END DO ! iM
	END DO ! sdY
END DO ! sdX


!$OMP PARALLEL DO private(im,ik,j,ptmp,ntmp,Pptmp,Pntmp,ialp,kgam,sdX,sdY) &
!$OMP shared(rho1b_jMK,Zmatinv) 
DO i = 1, indexerDim 
	sdX = indexer(i,1)
	sdY = indexer(i,2)
	im  = indexer(i,3)
	ik  = indexer(i,4)
	j   = indexer(i,5)
			    if(numprot >0)then
						ptmp = (0.d0,0.d0)
					end if

			    if(numneut >0)then
						ntmp = (0.d0,0.d0)
					end if
				if(.not.allsameparity)then
				    if(numprot >0)Pptmp=(0.d0,0.d0)			
				    if(numneut >0)Pntmp=(0.d0,0.d0)												
				end if

 			do ialp = 1,J2max1
 				do kgam = 1,J2max1	

 				   if(numprot>0)ptmp=ptmp+ Zmatinv(iM,kgam)*Zmatinv(iK,ialp)*Xrho1b(sdX,sdY,ialp,j,kgam)%prho
 				   if(numneut>0)ntmp=ntmp+ Zmatinv(iM,kgam)*Zmatinv(iK,ialp)*Xrho1b(sdX,sdY,ialp,j,kgam)%nrho
				   if(.not.allsameparity)then
					   if(numprot>0)then
						   Pptmp=Pptmp+ Zmatinv(iM,kgam)*Zmatinv(iK,ialp)*Xrho1b(sdX,sdY,ialp,j,kgam)%Pprho
						   
					   end if
					   if(numneut>0)then
						   Pntmp=Pntmp+ Zmatinv(iM,kgam)*Zmatinv(iK,ialp)*Xrho1b(sdX,sdY,ialp,j,kgam)%Pnrho
						   
					   end if  
					   
				   end if


 				end do ! kgam
 			end do  ! ialp
			

				if(numprot>0) rho1b_jMK(sdX,sdY,j,iM,iK)%prho = ptmp
				if(numneut>0) rho1b_jMK(sdX,sdY,j,iM,iK)%nrho = ntmp
				if(.not.allsameparity .and. numprot> 0)rho1b_jMK(sdX,sdY,j,iM,iK)%Pprho = Pptmp
				if(.not.allsameparity .and. numneut> 0)rho1b_jMK(sdX,sdY,j,iM,iK)%Pnrho = Pntmp
				

 	end do  ! i
!$OMP end parallel do

print*,' done inverting '
	
 	return
	
 end subroutine invert_alpha_and_gamma4dense

!==================================================
!==================================================
!
! cALLED BY: projectorator_dens
!
subroutine allocate_tilded_density
!	use system_parameters
	use spstate,only:nsps
	use dens1body
	implicit none
    integer :: isd,jsd,iJ,iM,iK
	
	integer ::Jstart
	real(4) :: Jmin
	real(4) :: xM,xK

    allocate(rhotild_Jprime_MK(numsd,numsd,numOfJ,J2max1,J2max1))

    do isd = 1,numsd
		do jsd = 1,numsd
				do iM = 1,J2max1
					do iK = 1,J2max1
					    xM = float(iM)-(Jmax+1.)
					    xK = float(iK)-(Jmax+1.)
	
						Jmin = max(abs(xM),abs(xK))
						Jstart = startJ2list(nint(2*Jmin))
						
						do iJ = jstart,numOfJ
						if(numprot> 0)then
							allocate(rhotild_Jprime_MK(isd,jsd,iJ,iM,iK)%prho(1:nsps,1:nsps))
							rhotild_Jprime_MK(isd,jsd,iJ,iM,iK)%prho=(0.d0,0.d0)
						end if
						if(numneut> 0)then
								allocate(rhotild_Jprime_MK(isd,jsd,iJ,iM,iK)%nrho(nsps,nsps))
								rhotild_Jprime_MK(isd,jsd,iJ,iM,iK)%nrho=(0.d0,0.d0)
						end if
						if(.not.allsameparity)then
							if(numprot> 0)allocate(rhotild_Jprime_MK(isd,jsd,iJ,iM,iK)%Pprho(nsps,nsps))
							if(numneut> 0)allocate(rhotild_Jprime_MK(isd,jsd,iJ,iM,iK)%Pnrho(nsps,nsps))							
							
						end if
					end do ! iJ
						
					end do ! iK
					
				end do ! im
				
			
		end do ! jsd
		
	end do ! isd

	return
end subroutine allocate_tilded_density
!==================================================
!==================================================
!
! subroutine to compute \tilde{N}^Jprime_MK
! inspired by SVD
!
!  here ~N^J'_MK = sum d^J'_MK (beta_j)  N_jMK
! 
!  for fixed M,K  (important!)
!
!  CALLED BY  projectorator_dens
!
subroutine compute_tilded_density(intM,intK)
!	use system_parameters
	use dens1body
	use spstate,only:nsps
!	use lamplight
	implicit none
	

	
	integer :: intK, intM
	integer :: aJ,bJ,ij
	real(8) :: xM,xK,xJ,xJp
	real(4) :: Jmin
	integer :: localNofJs
	complex(kind=8),allocatable :: pztmp(:,:),nztmp(:,:),Ppztmp(:,:),Pnztmp(:,:)

	
	real(4) :: wigtmp
	integer :: jstart
	integer :: sdX,sdY   ! indices for multiple SDs
		
	if(.not.allocated(pztmp))allocate(pztmp(1:nsps,1:nsps))
	if(.not.allocated(nztmp))allocate(nztmp(nsps,nsps))
	if(.not.allocated(ppztmp))allocate(Ppztmp(nsps,nsps))
	if(.not.allocated(pnztmp))allocate(pnztmp(nsps,nsps))


!........... RECONSTRUCT LIMITS...............
    xM = real(intM,kind=8)-(Jmax+1.)
    xK = real(intK,kind=8)-(Jmax+1.)
	
	Jmin = max(abs(xM),abs(xK))
	Jstart = startJ2list(nint(2*Jmin))
	do sdX = 1,numsd
	do sdY = 1,numsd
	do aJ = Jstart,numOfJ
		xJp = real(xJlist(aJ),kind=8)
		pztmp =(0.d0,0.d0)
		nztmp =(0.d0,0.d0)
		Ppztmp =(0.d0,0.d0)
		Pnztmp =(0.d0,0.d0)

		do iJ = 1,numOfJ
			wigtmp = wigner_d(xJp,xM,xK,beta_j(ij))
			
			if(numprot> 0)pztmp = pztmp + wigtmp * rho1b_jMK(sdX,sdY,ij,intM,intK)%prho
			if(numneut> 0)nztmp = nztmp + wigtmp * rho1b_jMK(sdX,sdY,ij,intM,intK)%nrho
			if(.not.allsameparity)then
				if(numprot> 0)Ppztmp = Ppztmp + wigtmp * rho1b_jMK(sdX,sdY,ij,intM,intK)%Pprho
				if(numneut> 0)Pnztmp = Pnztmp + wigtmp * rho1b_jMK(sdX,sdY,ij,intM,intK)%Pnrho									
				
			end if

		end do

		if(numprot> 0)rhotild_Jprime_MK(sdX,sdY,aJ,intM,intK)%prho= pztmp
		if(numneut> 0)rhotild_Jprime_MK(sdX,sdY,aJ,intM,intK)%nrho= nztmp
		if(.not.allsameparity)then
			if(numprot> 0)rhotild_Jprime_MK(sdX,sdY,aJ,intM,intK)%Pprho= Ppztmp
			if(numneut> 0)rhotild_Jprime_MK(sdX,sdY,aJ,intM,intK)%Pnrho= Pnztmp
	
		end if


	end do 
    end do  
    end do
	
	return	
	
end subroutine compute_tilded_density

!==================================================
!
!  set up the final arrays
!  uses derived types ('jagged arrays')
!
!  want: N^J_MK, H^J_MK for both even and odd parities
!        
!  CALLED BY  projectorator_dens
!
subroutine allocate_rho_J_MK
!	use system_parameters
	use spstate,only:nsps
	use dens1body
	implicit none
	
	integer :: ij,iM,iK,Jstart
	real(4) :: xM,xK,Jmin
	real    :: xJ
	
	integer :: sizeOfMK
	integer :: sdX,sdY   ! indices for multiple SDs
	integer :: Jtmax,Jtmin
	
    allocate(rho1b_J_MK(numsd,numsd,numOfJ))

	
	do ij = 1,numOfJ
		xJ = xJlist(ij)
		sizeOfMK = nint(2*xJ)+1		
		do sdX = 1,numsd
			do sdY = 1,numsd
				rho1b_J_MK(sdX,sdY,ij)%dimMK = sizeOfMK
				allocate(rho1b_J_MK(sdX,sdY, ij)%MK(sizeOfMK,sizeOfMK))
				do iM = 1,sizeOfMK
					do iK = 1,sizeOfMK
						if(numprot > 0)then
							allocate(rho1b_J_MK(sdX,sdY,iJ)%MK(iM,iK)%prho(nsps,nsps))
							rho1b_J_MK(sdX,sdY,iJ)%MK(iM,iK)%prho=(0.d0,0.d0)
						end if
						if(numneut > 0)then
							allocate(rho1b_J_MK(sdX,sdY,iJ)%MK(iM,iK)%nrho(nsps,nsps))
							rho1b_J_MK(sdX,sdY,iJ)%MK(iM,iK)%nrho=(0.d0,0.d0)
							
						end if
						if(.not.allsameparity)then
							if(numprot > 0)then
								allocate(rho1b_J_MK(sdX,sdY,iJ)%MK(iM,iK)%Pprho(nsps,nsps))
								rho1b_J_MK(sdX,sdY,iJ)%MK(iM,iK)%Pprho=(0.d0,0.d0)

							end if
							if(numneut > 0)then
								allocate(rho1b_J_MK(sdX,sdY,iJ)%MK(iM,iK)%Pnrho(nsps,nsps))
								rho1b_J_MK(sdX,sdY,iJ)%MK(iM,iK)%Pnrho=(0.d0,0.d0)

							end if
							
						endif
						
					end do ! ik
				end do ! im
			end do ! sdY
		end do !sdZ
	end do ! iJ

	return
		  		  	  			  
end subroutine allocate_rho_J_MK
!==================================================
!
!  solve sum_J Delta^J'J N^J_MK = ~N^J'_MK
!  for fixed M,K
!  also find H^J_MK etc
!
! CALLS lapack routine ZGESV to carry out LU decomposition
! and solve linear equation; if more than one solution needed,
! addition solutions invoke lapack routine ZGETRS  
! (so don't need to redo LU decomposition)
!
! CALLED BY:  projectorator_dens
!
subroutine solve4projected_density(intM,intK)
!	    use system_parameters
		use lamplight
		use dens1body
		implicit none

	
!............... INPUT...............................

		integer :: intK, intM

!............... INTERNAL................		
        integer :: sizeOfNorm
		integer :: info
		complex(kind=8),allocatable :: ztmpvec(:)
		real(8) :: Jmin,xJ
		real(8) :: xM,xK
		integer(4) :: Jstart,Jshift
		complex(kind=8), allocatable :: ztmpmat(:,:)
!	    real (kind = 8),ALLOCATABLE :: rwork(:),bb(:)
        COMPLEX (kind = 8), ALLOCATABLE :: zwork(:)
		integer :: lwork		
		
		integer :: i,j,k
		integer :: sdX, sdY
!		complex(kind=8),allocatable :: pztmp(:,:),nztmp(:,:),Ppztmp(:,:),Pnztmp(:,:)
		complex(kind=8) :: zme
		real(kind=8) :: dtmp
		integer :: aJ,bJ
		integer :: localNofJs
		real(8) :: xJp
	

		
!......... COMPUTE DIMENSION...................		
        xM = float(intM)-(Jmax+1.)
        xK = float(intK)-(Jmax+1.)
        Jmin = max(abs(xM),abs(xK))
        Jstart = startJ2list(nint(2*Jmin))
		sizeOfNorm = numOfJ -Jstart +1

!.......... COPY OVER TO ARRAYS
		
		allocate(ztmpvec(sizeOfNorm),ztmpmat(sizeOfNorm,sizeOfNorm))
		allocate(zwork(3*sizeOfNorm))
		if(allocated(IPIVOT))deallocate(IPIVOT)
		allocate(IPIVOT(sizeOfNorm))

		deltaJpJ(:,:)=0.d0

		localNofJs = numOfJ-Jstart+1

		do aJ = Jstart,numOfJ
			xJp = xJlist(aJ)
			do bJ = Jstart,numOfJ
				xJ = xJlist(bJ)
				dtmp = 0.d0
	!			if(xM==0.0 .and. xK == 0.0)print*,aJ,xjp,bj,xj
			
				do j = 1, numOfJ
					dtmp = dtmp + wigner_d(xJp,xM,xK,beta_j(j))* &
	                            wigner_d(xJ,xM,xK,beta_j(j))
				end do
				deltaJpJ(aJ,bJ)=dtmp
		
			end do
		end do 


		do i = 1,sizeOfNorm

			do j = 1,sizeOfNorm
				ztmpmat(i,j)=cmplx(deltaJpJ(i+Jstart-1,j+Jstart-1),0.d0)
			end do
		end do
		
!...........................		

        call zgetrf(sizeOfNorm,sizeOfNorm,ztmpmat,sizeOfNorm,IPIVOT,info)
        call errorZGESV(INFO)
		
		if(info/=0)then
			print*,' error in LU decomposition in solve4projected_density '
			stop
		end if
		call zgetri(sizeOfNorm,ztmpmat,sizeOfNorm,IPIVOT,zwork,3*sizeOfNorm,info)
        call errorZGESV(INFO)
		
		do i = 1,sizeOfNorm
!....... MUST FIGURE OUT SHIFT...........
!
!  normally intM,intK run from 1 to J2max1 
!  which corresponds to -Jmax to +Jmax
!  here, however, it runs from -xJ to +xJ
!  so the shift must be Jmax-xJ
!
            xJ = xJlist(i+Jstart-1)
			Jshift = nint(Jmax-xJ)
			do k =1,sizeOfNorm
				zme = ztmpmat(i,k)
				
			
			   do sdX = 1,numSD
			  	  do sdY = 1,numSD
					  if(numprot>0)then
						  rho1b_J_MK(sdX,sdY,i+Jstart-1)%MK(intM-Jshift,intK-Jshift)%prho = & 
						   rho1b_J_MK(sdX,sdY,i+Jstart-1)%MK(intM-Jshift,intK-Jshift)%prho + & 
						  zme * rhotild_Jprime_MK(sdX,sdY,k+Jstart-1,intM,intK)%prho 
					  end if
					  if(numneut>0)then
						  rho1b_J_MK(sdX,sdY,i+Jstart-1)%MK(intM-Jshift,intK-Jshift)%nrho = & 
						  rho1b_J_MK(sdX,sdY,i+Jstart-1)%MK(intM-Jshift,intK-Jshift)%nrho + & 
						  zme * rhotild_Jprime_MK(sdX,sdY,k+Jstart-1,intM,intK)%nrho 
					  end if
					  if(.not.allsameparity)then
						  if(numprot>0)then
							  rho1b_J_MK(sdX,sdY,i+Jstart-1)%MK(intM-Jshift,intK-Jshift)%Pprho = & 
							   rho1b_J_MK(sdX,sdY,i+Jstart-1)%MK(intM-Jshift,intK-Jshift)%Pprho + & 
							  zme * rhotild_Jprime_MK(sdX,sdY,k+Jstart-1,intM,intK)%Pprho 

						  end if
						  if(numneut>0)then
							  rho1b_J_MK(sdX,sdY,i+Jstart-1)%MK(intM-Jshift,intK-Jshift)%Pnrho = & 
							  rho1b_J_MK(sdX,sdY,i+Jstart-1)%MK(intM-Jshift,intK-Jshift)%Pnrho + & 
							  zme * rhotild_Jprime_MK(sdX,sdY,k+Jstart-1,intM,intK)%Pnrho 
							  
						  end if				  

					  end if
					
				  end do  !sdy
				
			   end do ! sdX
		   end do  !k

		end do  !u

deallocate(IPIVOT,ztmpvec,ztmpmat)

		
end subroutine solve4projected_density

!==============================================================================

!==============================================================================

subroutine setup_transformed_rhos
	use sporbit
	use dens1body
	use lamplight
!	use system_parameters
	implicit none
	integer :: sdf,sdi,iJ,fJ,iM,fM
	real :: xJi,xJf,xMi,xMf,xK
	integer :: Jtmin,Jtmax
	integer :: Jt
	integer :: i
	integer :: Jspmax
	integer :: iJshift,fJshift
	
	Jspmax = 0
	do i = 1,numorb
		Jspmax = max(Jspmax,orbqn(i)%j)
	end do
	
	allocate(rhotr(numSD,numSD,numOfJ,J2max1,numofJ,J2max1))
	
	rhotr%Jtmin = -999
	rhotr%Jtmax = -999
	
	do sdf = 1,numSD
		do sdi = 1, numSD
			
			do iJ= 1,numOfJ
				xJi = xJlist(iJ)
				iJshift = nint(Jmax-xJi)
				
				do fJ = 1,numOfJ
					xJf = xJlist(fJ)
					fJshift = nint(Jmax-xJf)
					
					Jtmax = nint(xJi+xJf)
					Jtmax = min(Jtmax,jspmax)
					Jtmin = nint(abs(xJi-xJf))
					if(Jtmax < Jtmin)cycle
					do iM = 1,J2max1
				        xMi = float(iM)-(Jmax+1.)
						if(xJi < abs(xMi))cycle
						do fM = 1,J2max1
					        xMf = float(fM)-(Jmax+1.)
							if(xJf < abs(xMf))cycle
							
							rhotr(sdf,sdi,fJ,fM-fJshift,iJ,iM-iJshift)%Jtmin = Jtmin
							rhotr(sdf,sdi,fJ,fM-fJshift,iJ,iM-iJshift)%Jtmax = Jtmax
							
							
							allocate(rhotr(sdf,sdi,fJ,fM-fJshift,iJ,iM-iJshift)%rhoJt(Jtmin:Jtmax))
							do Jt = Jtmin,Jtmax
								allocate(rhotr(sdf,sdi,fJ,fM-fJshift,iJ,iM-iJshift)%rhoJt(Jt)%prho(numorb,numorb) )
								allocate(rhotr(sdf,sdi,fJ,fM-fJshift,iJ,iM-iJshift)%rhoJt(Jt)%nrho(numorb,numorb) )
								rhotr(sdf,sdi,fJ,fM-fJshift,iJ,iM-iJshift)%rhoJt(Jt)%prho=(0.d0,0.d0)
								rhotr(sdf,sdi,fJ,fM-fJshift,iJ,iM-iJshift)%rhoJt(Jt)%nrho=(0.d0,0.d0)

								if(.not.allsameparity)then
									allocate(rhotr(sdf,sdi,fJ,fM-fJshift,iJ,iM-iJshift)%rhoJt(Jt)%Pprho(numorb,numorb) )
									allocate(rhotr(sdf,sdi,fJ,fM-fJshift,iJ,iM-iJshift)%rhoJt(Jt)%Pnrho(numorb,numorb) )								
									rhotr(sdf,sdi,fJ,fM-fJshift,iJ,iM-iJshift)%rhoJt(Jt)%Pprho=(0.d0,0.d0)
									rhotr(sdf,sdi,fJ,fM-fJshift,iJ,iM-iJshift)%rhoJt(Jt)%Pnrho=(0.d0,0.d0)
									
								end if
			
							end do
							
						end do ! ik
					end do ! iM

				end do ! fJ
				
			end do ! iJ
		
		end do ! sdi
		
	end do  ! sdf
	
	return
end subroutine setup_transformed_rhos

!==============================================================================
!
subroutine transform_rhos_step1
	use dens1body
	use lamplight
!	use system_parameters
	use spstate
	implicit none
	integer :: sdf,sdi,iJ,fJ
	integer :: iMi,iMf,iKp,Mt
	real :: xJi,xJf,xMi,xMf,xKp,xMt
	integer :: Jtmin,Jtmax
	integer :: a,b,ia,ib
	integer :: Jt
	real    :: fact,fact2
	real :: clebr,cleb
	integer :: Jspmax
	integer :: iabindx,abindx
	integer :: ja,jb,ma,mb
	complex(8):: sum
	integer :: Jshift,iJshift,fJshift

	INTEGER :: i, indexerDim
	INTEGER, ALLOCATABLE :: indexer(:,:)
	
	
	type rhotrans
		real(4), allocatable :: trans(:,:)
	end type rhotrans
	
	type (rhotrans), allocatable :: transrhoJM(:,:)
	
	
	Jspmax = 0
	do a = 1,numorb
		Jspmax = max(Jspmax,orbqn(a)%j)
	end do
	
	allocate(transrhoJM(0:Jspmax,-Jspmax:Jspmax))
	
	do Jt = 0,Jspmax
		do Mt = -Jt,Jt
			allocate(transrhoJM(Jt,Mt)%trans(numorb**2,nsps**2))
			
			transrhoJM(Jt,Mt)%trans= 0.d0
			
			do ia = 1,nsps
				ja = spsqn(ia)%j
				ma = spsqn(ia)%m
				a = spsqn(ia)%orb
				
				do ib = 1,nsps
					jb = spsqn(ib)%j
					mb = spsqn(ib)%m
					b = spsqn(ib)%orb
					iabindx = (ia-1)*nsps + ib
					abindx = (a-1)*numorb+b
					
					transrhoJM(Jt,Mt)%trans(abindx,iabindx)=cleb(ja,ma,jb,-mb,2*Jt,2*Mt)*(-1)**( (jb-mb)/2)
					
				end do
			end do
			
		end do
		
	end do

! Dimension of new index: numOfJ**2 * J2max1 ** 2
! Unrolling fJ, iJ, iMi, iMf 
! Added sdf and sdi to unrolled loop on Sept 10th
! for first pass 
! MODIFIED REGION
indexerDim = (J2Max1**2)*(numOfJ**2)*(numsd**2)
ALLOCATE(indexer(indexerDim,6))

i = 1
DO sdf = 1, numsd 
	DO sdi = 1, numsd 
		DO fJ = 1, numOfJ
			DO iJ = 1, numOfJ
				DO iMi = 1, J2max1
					DO iMf = 1, J2max1
						indexer(i,1) = sdf
						indexer(i,2) = sdi 
						indexer(i,3) = fJ  
						indexer(i,4) = iJ 
						indexer(i,5) = iMi 
						indexer(i,6) = iMf 
						i = i + 1 
					END DO ! iJ
				END DO ! fJ
			END DO ! iMf
		END DO ! iMi
	END DO ! sdY
END DO ! sdX

			! MODIFIED REGION
!$OMP PARALLEL DO private(fJ, iJ, iMi, iMf, iKp, Jt, ia, ib, a, b, xJf, fJshift, xJi, Jtmax, Jtmin, Jshift, iJshift, & 
!$OMP xMi, xMf, xKp, xMt, Mt, ma, mb, iabindx, abindx, fact, fact2, sdf, sdi) &
!$OMP default(SHARED) 

!xJf, fJshift, xJi, Jtmax, Jtmin, Jshift, iJshift, xMi, xMf, xKp, xMt, Mt, fact, ma, mb, iabindx, abindx, fact2) 
!
			DO i = 1, indexerDim
			! MODIFIED REGION
				sdf = indexer(i,1)
				sdi = indexer(i,2)
				fJ  = indexer(i,3)
				iJ  = indexer(i,4)
				iMi = indexer(i,5)
				iMf = indexer(i,6)
				!do fJ= 1,numOfJ 
			
				xJf = xJlist(fJ) ! this is J_alpha in my notes
				fJshift = nint(Jmax-xJf)
				!do iJ = 1,numOfJ   
				
				xJi = xJlist(iJ) ! this is J_beta in my notes
				Jtmax = nint(xJi+xJf)
				Jtmin = nint(abs(xJi-xJf))
				Jshift = nint(Jmax-xJi)
				iJshift = Jshift

				!do iMi = 1,J2max1 ! this is M_beta in my notes
				xMi = float(iMi)-(Jmax+1.) ! this is K_mu (associate with beta= initial ) in my notes
				if(xJi < abs(xMi))cycle

				!do iMf = 1, J2max1 ! this is M_alpha in my notes
				xMf = float(iMf)-(Jmax+1.) ! this is K_lambda (associated with alpha=final) in my notes
				if(xJf < abs(xMf))cycle
							
				do iKp = 1,J2max1
					xKp = float(iKp)-(Jmax+1.) ! this is Klambda - Mt in my notes
					if(xJi < abs(xKp))cycle 
					xMt = xMf - xKp
					Mt = int(xMt) ! NEED TO CONSTRAIN Mt

					Jtmin = rhotr(sdf,sdi,fJ,iMf-fJshift,iJ,iMi-iJshift)%Jtmin 
					Jtmax = rhotr(sdf,sdi,fJ,iMf-fJshift,iJ,iMi-iJshift)%Jtmax

					if(Jtmin < 0 .or. Jtmax < 0)cycle

					if(abs(Mt) > Jtmax)cycle
					Jtmin = max(Jtmin,abs(Mt))

					do Jt = Jtmin,Jtmax
						fact = sqrt(2*xJf+1.0)/sqrt(2.0*Jt+1.0)
						fact = fact*clebr(xJi,xKp,float(Jt),xMt,xJf,xMf)
						!....... NOW I HAVE TO PULL OUT (a,b,Jt,Mt)
						! here a, b are orbit labels, ia, ib s.p. state labels
						! this can probably be done more efficiently and compactly, but start with this

						do ia = 1,nsps
							ma = spsqn(ia)%m
							do ib = 1,nsps
								mb = spsqn(ib)%m
								if(ma-mb /= 2*Mt)cycle
								iabindx = (ia-1)*nsps+ib
															
								do a = 1,numorb
									do b = 1,numorb
										abindx = (a-1)*numorb +b
										fact2 = transrhoJM(Jt,Mt)%trans(abindx,iabindx)

										if(fact2 == 0.0)cycle

								    if(numprot>0) rhotr(sdf,sdi,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%prho(a,b) = & 
									                rhotr(sdf,sdi,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%prho(a,b) + & 
									                fact*fact2*rho1b_J_MK(sdf,sdi,iJ)%MK(iKp-Jshift,iMi-Jshift)%prho(ia,ib)

						    		if(numneut>0) rhotr(sdf,sdi,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%nrho(a,b) = & 
									                rhotr(sdf,sdi,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%nrho(a,b) + & 
									                fact*fact2* rho1b_J_MK(sdf,sdi,iJ)%MK(iKp-Jshift,iMi-Jshift)%nrho(ia,ib) 
													
									if(.not.allsameparity)then
									    if(numprot>0) rhotr(sdf,sdi,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%Pprho(a,b) = & 
										                rhotr(sdf,sdi,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%Pprho(a,b) + & 
										                fact*fact2*rho1b_J_MK(sdf,sdi,iJ)%MK(iKp-Jshift,iMi-Jshift)%Pprho(ia,ib)

							    		if(numneut>0) rhotr(sdf,sdi,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%Pnrho(a,b) = & 
										                rhotr(sdf,sdi,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%Pnrho(a,b) + & 
										                fact*fact2* rho1b_J_MK(sdf,sdi,iJ)%MK(iKp-Jshift,iMi-Jshift)%Pnrho(ia,ib) 

																											
									end if				
										end do ! b
									end do ! a
								end do ! ib
							end do ! ia
						end do ! jt
					end do ! ikp

				
			end do ! i
!$OMP END PARALLEL DO
	
	return
	
end subroutine transform_rhos_step1
!----------------------------------------------------------------------------------------------	
!
! write to a file for post-processing
!
! key flag:  binary_dens
!
subroutine writeout_dens
	use dens1body
	use lamplight
!	use phf_vals
!	use system_parameters
	implicit none
	integer :: sdf,sdi,iJ,fJ,iM,fM

	integer :: ii,ff,a,b,kk,ll,jj


	character*80 ::filename
	integer :: ilast
	
	integer :: intJ,matdim
	real :: floatJ
	integer :: sdX, sdY
	integer :: iMi,iMf,iJshift,fJshift
	real :: xJi,xJf,xMi,xMf
	integer :: Jt,Jtmin,Jtmax
	real(8) :: sum
	
	integer :: delJshift ! needed to account for difference between J2max1 and J2max1out
!	integer :: a,b
	
	delJshift = (J2max1-J2max1out)/2
	
	print*,' Writing out density matrices to file for postprocessing with lilLAMP '
	if(binary_dens)then
		write(6,*)' (Writing out as a binary file; if you want to change, '
		write(6,*)' then set binary_dens = .false. in LAMPproject2lib2.f90 )'
	else
		
		write(6,*)' (Writing out as ASCII file; if you want to write as binary, '
		write(6,*)' then set binary_dens = .true. in LAMPproject2lib2.f90 )'	
		
	end if
	
	print*,' '
	print*,' Enter name of .denmat file (do not enter extension)'
	print*,' (will overwrite existing data if file already exists)'
	read(5,'(a)')filename
	ilast = index(filename,' ')-1
	if(binary_dens)then
       open(unit=31,file=filename(1:ilast)//".denmat",status='unknown',form='unformatted')
		
	else
   	   open(unit=3,file=filename(1:ilast)//".denmat",status='unknown')
    end if
	
	if(binary_dens)then
 	   write(31)numorb,nsps
 	   write(31)numprot,numneut,numsd
 	   write(31)J2max1out,numOfJout,allsameparity		
		
	else
	   write(3,*)numorb,nsps
	   write(3,*)numprot,numneut,numsd
	   write(3,*)J2max1out,numOfJout,allsameparity
   end if
!.......... NEED TO WRITE OUT ORBITS.........

    do a = 1,numorb
		if(binary_dens)then
			write(31)orbqn(a)%j,orbqn(a)%l,orbqn(a)%nr		
			
		else
			write(3,*)orbqn(a)%j,orbqn(a)%l,orbqn(a)%nr		
			
		end if
	end do	
    do iJ = 1, numOfJout
		
		
        xJi = xJlist(iJ)
		iJshift = nint(Jmax-xJi)
		
		do fJ = 1,numOfJout
			xJf = xJList(fJ)
			fJshift = nint(Jmax-xJf)
!			print*,xJi,xJf,delJshift
			do iMi = 1,J2max1 ! this is M_beta in my notes
!			do iMi = 1,J2max1out ! this is M_beta in my notes

				xMi = float(iMi)-(Jmax+1.) ! this is K_mu (associate with beta= initial ) in my notes
				if(xJi < abs(xMi))cycle
				
				do iMf = 1, J2max1 ! this is M_alpha in my notes
!				do iMf = 1,J2max1out ! this is M_alpha in my notes

					xMf = float(iMf)-(Jmax+1.) ! this is K_lambda (associated with alpha=final) in my notes
					if(xJf < abs(xMf))cycle
				
					Jtmin = rhotr(1,1,fJ,iMf-fJshift,iJ,iMi-iJshift)%Jtmin 
					Jtmax = rhotr(1,1,fJ,iMf-fJshift,iJ,iMi-iJshift)%Jtmax
					
					if(Jtmin > Jtmax)cycle
					if(Jtmin==-999 .or. Jtmax ==-999)cycle
					
					do Jt = Jtmin,Jtmax
						do sdX = 1,numsd
							do sdY = 1,numsd
								if(binary_dens)then
									write(31)sdx,sdy,fJ,iMf-delJshift,iJ,iMi-delJshift,Jt
									if(numprot> 0)then
											write(31)rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%prho(:,:)

									end if
									if(numneut> 0)then

											write(31)rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%nrho(:,:)
	
									end if
									if(.not.allsameparity)then
										if(numprot> 0)then
												write(31)rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%Pprho(:,:)

										end if
										if(numneut> 0)then

												write(31)rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%Pnrho(:,:)
	
										end if		
										
									end if
			
								else
									write(3,*)sdx,sdy,fJ,iMf-delJshift,iJ,iMi-delJshift,Jt
									if(numprot> 0)then
											write(3,*)rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%prho(:,:)

									end if
									if(numneut> 0)then

											write(3,*)rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%nrho(:,:)
	
									end if
									if(.not.allsameparity)then
										if(numprot> 0)then
												write(3,*)rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%Pprho(:,:)

										end if
										if(numneut> 0)then

												write(3,*)rhotr(sdx,sdy,fJ,iMf-fJshift,iJ,iMi-iJshift)%rhoJt(Jt)%Pnrho(:,:)
	
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
			

	if(binary_dens)then
		close(31)
	else
		close(3)
	end if
	return

	
end subroutine writeout_dens



!----------------------------------------------------------------------------------------------	

end module project1bdense
