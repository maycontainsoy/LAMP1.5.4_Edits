

!
! central sampling routine
! rewritten afresh May 2021 version 1.5.1
!
! NOTE: the euler angles are referenced two different ways:
!   first, each angle has its own list: alpha_i, beta_j, gamma_k, stored in module lamplight.
!   each alpha_i etc is unique; these are set in routine set_euler_meshes in module lamplight (in LAMPutils.f90)
!   However for sampling, we have a list of meshpoints and alpha_list, beta_list, gamma_list
!   In this case you may have alpha_list constant (for some ways) while beta_list changes etc.
!   The reason for this is to make it easier to divide up across MPI processes.

subroutine basic_sampling
	
	use lamplight
	use wignersfriend
	use densitylib
	use applyh
	implicit none
	
	logical :: errflag

	logical, ALLOCATABLE :: Check(:,:,:)  !Added by Changfeng Jiao
	real(kind=8) :: Amass !added by Changfeng Jiao
	integer :: iphase  ! phase needed for symmetries
	integer(kind=4) :: nu, icount      !added by Changfeng Jiao

	integer,allocatable :: alpha_list(:),beta_list(:),gamma_list(:)	 ! the Euler angles from a serial list of meshpoints
	logical,allocatable :: useSymmetry(:)   ! if true, then use symmetry to find values
	integer :: ialp,jbet,kgam
	integer :: ilist,numEuler 
	complex (KIND = 8), ALLOCATABLE ::  parOp(:,:)  ! matrix for making parity transformation
	
	COMPLEX (KIND = 8), allocatable :: RotMat(:,:)
	COMPLEX (KIND = 8),ALLOCATABLE :: psdr(:,:),nsdr(:,:),psdpr(:,:),nsdpr(:,:) 
	COMPLEX (KIND = 8),ALLOCATABLE :: rhopij(:,:),rhonij(:,:)
	COMPLEX (KIND = 8) :: ovlpp,ovlpn
	COMPLEX (KIND = 8) :: vme
	
	complex(kind=8) :: zone = (1.d0,0.d0)
	complex(kind=8) :: zzero = (0.d0,0.d0)
	
	integer :: isd,fsd  ! indices for Slater determinants
	
	integer :: i,j

	INTEGER :: ierr ! for MPI communication
		
!.... DETERMINE WHAT SYMMETRY PATTERN WE USE...SEE PAPER
    if (MOD(J2max1,4) .eq. 1) then
       nu = INT(Jmax)
    else if(MOD(J2max1,4) .eq. 2) then
       nu = INT(Jmax - 0.5)
    else if (MOD(J2max1, 4) .eq. 3) then
       nu = INT(Jmax - 1.)
    else if (MOD(J2max1, 4) .eq. 0) then
       nu = INT(Jmax + 0.5)
    end if	
	
!... SET UP ARRAYS FOR LATTICE OF EULER ANGLES.....	
	allocate(alpha_list(J2max1**2*numOfJ),beta_list(J2max1**2*numOfJ),gamma_list(J2max1**2*numOfJ) )
!....SET UP ARRAYS FOR USING SYMMETRIES TO REDUCE WORK....	
	allocate(useSymmetry( J2max1**2*numOfJ))
	ALLOCATE (Check(0:J2max1,numofJ,0:J2max1))   !Added by Changfeng Jiao
	Check( 0:J2max1, 1:numofJ, 0:J2max1) = .false.
	Amass = dble(numprot+numneut)  ! needed for computing phases while using symmetries
	
!............ NOW SET UP FULL 3-D EULER MESH, USING SYMMETRIES...............
!  the individual meshes alpha_i, beta_j, gamma_k, as well as the  array prime (which determines symmetries)
!  have already been set up in set_euler_meshes in LAMPutils.f90
!	
	ilist = 0
	useSymmetry=.false.
	
	
!	beta_j = 0.d0
	do ialp= 1,J2max1
		do jbet=1,numOfJ
			do kgam=1,J2max1
				if(check(prime(kgam),jbet,prime(ialp)))cycle
				ilist =ilist+1
				alpha_list(ilist)=ialp
				beta_list(ilist)=jbet			
				gamma_list(ilist)=kgam
				if(prime(ialp) .ne. 0 .and. prime(kgam) .ne. 0 .and. .not. check(ialp,jbet,kgam) )then
					check(ialp,jbet,kgam)=.true.
					useSymmetry(ilist)=.true.
				end if
			end do
		end do
	end do		
	numEuler= ilist ! this is the number of mesh points actually sampled; the rest are filled in by symmetries
	
	IF (myMPIrank == root) PRINT*,' Use only ',numEuler,' Euler angles out of ',J2max1*numOfJ*J2max1
	
	allocate(parOp(nsps,nsps))
	call MakeParityOp(parOP)

	allocate(RotMat(nsps,nsps))
	ALLOCATE (psdr(nsps,numprot),nsdr(nsps,numneut))
	if(.not.allsameparity)ALLOCATE (psdpr(nsps,numprot),nsdpr(nsps,numneut))
	ALLOCATE(rhopij(nsps,nsps),rhonij(nsps,nsps))

!.......... NOW LOOP OVER ALL EULER MESH POINTS except those accesses via symmetry....	
	! distribute with MPI 
	do ilist = 1,numEuler
		ialp = alpha_list(ilist)
		jbet = beta_list(ilist)
		kgam = gamma_list(ilist)
		
		call make_rotmat(alpha_i(ialp),beta_j(jbet),gamma_k(kgam),RotMat)

		do isd = 1,numsd

			psdr=zzero
			nsdr=zzero
			if(numprot> 0)call zgemm('N','N',nsps,numprot,nsps,zone,RotMat,nsps,psdf(isd,:,:),nsps,zzero,psdr,nsps)
			if(numneut> 0)call zgemm('N','N',nsps,numneut,nsps,zone,RotMat,nsps,nsdf(isd,:,:),nsps,zzero,nsdr,nsps)

			if(.not.allsameparity)then
				if(numprot> 0)call zgemm('N','N',nsps,numprot,nsps,zone,parOp,nsps,psdr,nsps,zzero,psdpr,nsps)
				if(numneut> 0)call zgemm('N','N',nsps,numneut,nsps,zone,parOp,nsps,nsdr,nsps,zzero,nsdpr,nsps)				
			end if

			do fsd = 1, numsd
				CALL makerhoij(1,numprot,psdf(fsd,:,:),psdr,ovlpp,rhopij,errflag)	! in LAMPdenslib.f90
				if(errflag)then
					print*,' Error between (proton) states ',isd,fsd
					stop
				end if
				CALL makerhoij(2,numneut,nsdf(fsd,:,:),nsdr,ovlpn,rhonij,errflag)
				
				if(errflag)then
					print*,' Error between (neutron) states ',isd,fsd
					stop
				end if

				N_ijk(fsd,isd,ialp,jbet,kgam) = ovlpp*ovlpn
				
				if(doHam)then
!					print*,ovlpp,ovlpn
!					print*,rhopij
					!CALL hmultMPIdistro ! DONT PUT THIS HERE
					CALL TBMEmaster(nsps,rhopij,rhonij,ovlpp,ovlpn,vme)	! in LAMPapplyh.f90
					H_ijk(fsd,isd,ialp,jbet,kgam) = vme
!					print*,vme
					PRINT*, ' Node = ', myMPIrank, ' end of TBMEmaster (in sampler)'
				end if
!.................... IF WE PROJECT PARITY.....
				if(.not.allsameparity)then
					CALL makerhoij(1,numprot,psdf(fsd,:,:),psdpr,ovlpp,rhopij,errflag)
					if(errflag)then
						print*,' Error between (proton) states ',isd,fsd
						stop
					end if
					CALL makerhoij(2,numneut,nsdf(fsd,:,:),nsdpr,ovlpn,rhonij,errflag)
					if(errflag)then
						print*,' Error between (neutron) states ',isd,fsd
						stop
					end if
					PN_ijk(fsd,isd,ialp,jbet,kgam) = ovlpp*ovlpn
					if(doHam)then
						CALL TBMEmaster(nsps,rhopij,rhonij,ovlpp,ovlpn,vme)
						PH_ijk(fsd,isd,ialp,jbet,kgam) = vme
					end if					
				end if				! END FINDING PARITY TRANSFORMED
!.......................... USE SYMMETRIES TO REDUCE SAMPLING................				
				if (useSymmetry(ilist) ) then
					if((ialp .LE. nu .and. kgam .LE. nu) .or. (ialp .GT. nu .and. kgam .GT. nu)) then
						iphase = (-1)**Amass
					else
						iphase =1
					end if

					N_ijk(isd,fsd,prime(kgam),jbet,prime(ialp)) =  iphase * conjg(N_ijk(fsd,isd,ialp,jbet,kgam))
					IF (doHam) THEN
						H_ijk(isd,fsd,prime(kgam),jbet,prime(ialp)) =  iphase * conjg(H_ijk(fsd,isd,ialp,jbet,kgam))
					END IF
					
					IF (.not.allsameparity) THEN
						PN_ijk(isd,fsd,prime(kgam),jbet,prime(ialp)) =  iphase * conjg(PN_ijk(fsd,isd,ialp,jbet,kgam))
						IF (doHam) THEN
							PH_ijk(isd,fsd,prime(kgam),jbet,prime(ialp))  =  iphase * conjg(PH_ijk(fsd,isd,ialp,jbet,kgam) )
						END IF
					END IF
				end if
			end do ! jsd
		end do ! isd
	end do ! ilist

	PRINT*, ' Node = ', myMPIrank, ' end of basic_sampling'
	
	deallocate(alpha_list,beta_list,gamma_list)
	return
	
end subroutine basic_sampling


subroutine density_sampling
	use lamplight
	use wignersfriend
	use densitylib
	use dens1body
	implicit none
	
	logical :: errflag

	logical, ALLOCATABLE :: Check(:,:,:)  !Added by Changfeng Jiao
	real(kind=8) :: Amass !added by Changfeng Jiao
	integer :: iphase  ! phase needed for symmetries
	integer(kind=4) :: nu, icount      !added by Changfeng Jiao

	integer,allocatable :: alpha_list(:),beta_list(:),gamma_list(:)	 ! the Euler angles from a serial list of meshpoints
	logical,allocatable :: useSymmetry(:)   ! if true, then use symmetry to find values
	integer :: ialp,jbet,kgam
	integer :: ilist,numEuler 
	COMPLEX (KIND = 8), ALLOCATABLE ::  parOp(:,:)  ! matrix for making parity transformation
	
	COMPLEX (KIND = 8), allocatable :: RotMat(:,:)
	COMPLEX (KIND = 8),ALLOCATABLE :: psdr(:,:),nsdr(:,:),psdpr(:,:),nsdpr(:,:)     ! right (initial) and parity reflection
	COMPLEX (KIND = 8),ALLOCATABLE :: rhopij(:,:),rhonij(:,:)
	COMPLEX (KIND = 8) :: ovlpp,ovlpn
	COMPLEX (KIND = 8) :: vme
	
	complex(kind=8) :: zone = (1.d0,0.d0)
	complex(kind=8) :: zzero = (0.d0,0.d0)
	
	integer :: isd,fsd  ! indices for Slater determinants
	
	integer :: i,j
		
!.... DETERMINE WHAT SYMMETRY PATTERN WE USE...SEE PAPER
    if (MOD(J2max1,4) .eq. 1) then
       nu = INT(Jmax)
    else if(MOD(J2max1,4) .eq. 2) then
       nu = INT(Jmax - 0.5)
    else if (MOD(J2max1, 4) .eq. 3) then
       nu = INT(Jmax - 1.)
    else if (MOD(J2max1, 4) .eq. 0) then
       nu = INT(Jmax + 0.5)
    end if	
	
!... SET UP ARRAYS FOR LATTICE OF EULER ANGLES.....	
	allocate(alpha_list(J2max1**2*numOfJ),beta_list(J2max1**2*numOfJ),gamma_list(J2max1**2*numOfJ) )
!....SET UP ARRAYS FOR USING SYMMETRIES TO REDUCE WORK....	
	allocate(useSymmetry( J2max1**2*numOfJ))
	ALLOCATE (Check(0:J2max1,numofJ,0:J2max1))   !Added by Changfeng Jiao
	Check( 0:J2max1, 1:numofJ, 0:J2max1) = .false.
	Amass = dble(numprot+numneut)  ! needed for computing phases while using symmetries
	
!............ NOW SET UP FULL 3-D EULER MESH, USING SYMMETRIES...............
!  the individual meshes alpha_i, beta_j, gamma_k, as well as the  array prime (which determines symmetries)
!  have already been set up in set_euler_meshes in LAMPutils.f90
!	
	ilist = 0
	useSymmetry=.false.

	!useSymmetry=.false.
	do ialp= 1,J2max1
		do jbet=1,numOfJ
			do kgam=1,J2max1
				   ilist =ilist+1
				   alpha_list(ilist)=ialp
				   beta_list(ilist)=jbet			
				   gamma_list(ilist)=kgam

			end do
		end do
	end do		
	numEuler= ilist
	print*,' Use only ',numEuler,' Euler angles out of ',J2max1*numOfJ*J2max1	

	
	allocate(parOp(nsps,nsps))
	call MakeParityOp(parOP)
	
	allocate(RotMat(nsps,nsps))
	ALLOCATE (psdr(nsps,numprot),nsdr(nsps,numneut))
    if(.not.allsameparity)ALLOCATE (psdpr(nsps,numprot),nsdpr(nsps,numneut))
	
	ALLOCATE(rhopij(nsps,nsps),rhonij(nsps,nsps))

!.......... NOW LOOP OVER ALL EULER MESH POINTS except those accesses via symmetry....	
	do ilist = 1,numEuler
		
		ialp = alpha_list(ilist)
		jbet = beta_list(ilist)
		kgam = gamma_list(ilist)
		
		call make_rotmat(alpha_i(ialp),beta_j(jbet),gamma_k(kgam),RotMat)
		
		do isd = 1,numsd
			psdr=zzero
			nsdr=zzero
			if(numprot> 0)call zgemm('N','N',nsps,numprot,nsps,zone,RotMat,nsps,psdf(isd,:,:),nsps,zzero,psdr,nsps)
			if(numneut> 0)call zgemm('N','N',nsps,numneut,nsps,zone,RotMat,nsps,nsdf(isd,:,:),nsps,zzero,nsdr,nsps)
			
			if(.not.allsameparity)then
				psdpr=zzero
				nsdpr=zzero
				if(numprot> 0)call zgemm('N','N',nsps,numprot,nsps,zone,parOp,nsps,psdr,nsps,zzero,psdpr,nsps)
				if(numneut> 0)call zgemm('N','N',nsps,numneut,nsps,zone,parOp,nsps,nsdr,nsps,zzero,nsdpr,nsps)	
							
			end if

			do fsd = 1, numsd


	            CALL makerhoij(1,numprot,psdf(fsd,:,:),psdr,ovlpp,rhopij,errflag)
				if(errflag)then
					print*,' Error between (proton) states ',isd,fsd
					stop
				end if
	            CALL makerhoij(2,numneut,nsdf(fsd,:,:),nsdr,ovlpn,rhonij,errflag)
				
				if(errflag)then
					print*,' Error between (neutron) states ',isd,fsd
					stop
				end if

	            N_ijk(fsd,isd,ialp,jbet,kgam) = ovlpp*ovlpn
					if(numprot > 0)then
					    Xrho1b(fsd,isd,ialp,jbet,kgam)%prho=rhopij*ovlpp*ovlpn

					end if

					if(numneut > 0)then
					    Xrho1b(fsd,isd,ialp,jbet,kgam)%nrho=rhonij*ovlpp*ovlpn
					end if
!.................... IF WE PROJECT PARITY.....
                if(.not.allsameparity)then
		            CALL makerhoij(1,numprot,psdf(fsd,:,:),psdpr,ovlpp,rhopij,errflag)
					if(errflag)then
						print*,' Error between (proton) states ',isd,fsd
						stop
					end if
		            CALL makerhoij(2,numneut,nsdf(fsd,:,:),nsdpr,ovlpn,rhonij,errflag)
					if(errflag)then
						print*,' Error between (neutron) states ',isd,fsd
						stop
					end if
		            PN_ijk(fsd,isd,ialp,jbet,kgam) = ovlpp*ovlpn
					if(numprot > 0)then
					    Xrho1b(fsd,isd,ialp,jbet,kgam)%Pprho=rhopij*ovlpp*ovlpn

					end if

					if(numneut > 0)then
					    Xrho1b(fsd,isd,ialp,jbet,kgam)%Pnrho=rhonij*ovlpp*ovlpn
					end if

				
				end if				! END FINDING PARITY TRANSFORMED

				
				
			end do ! jsd
			
			
		end do ! isd

	
	end do ! ilist
	
	deallocate(alpha_list,beta_list,gamma_list)
	
	
	return
	
	
end subroutine density_sampling
!end module sampler