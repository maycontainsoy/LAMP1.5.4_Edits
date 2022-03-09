
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
module jaggedArrayType

  type jaggedArray
		complex(kind=8), allocatable, dimension(:,:) :: MK
		integer(kind=4) :: dimMK
  end type jaggedArray

end module jaggedArrayType
	  
!========================================================
!
!  variables needed for solving the PHF equations
!
!  IMPORTANT NOTE ON MAPPING OF INDICES
!
!  MAPPING OF M,K:
!   index(M) = M - Jmax +1
!
!  DEFAULT MAPPING OF J
!   if( even )
!     index(J)= J+1
!  else
!     index(J)= J+1/2
!
! SOME IMPORTANT VARIABLES
!
!  real :: Jmax
!  J2max1 = 2Jmax+1, replaces bigJmax 
!  numOfJ = int(Jmax +1)
	
module lamplight
	USE nodeinfo ! added by SML
	use jaggedArrayType
	use psis
	use sporbit,only:allsameparity
	use spstate
	IMPLICIT NONE
		  
	logical :: doHam                      ! if computing the Hamiltonian
	real (kind=4) :: Jmax
	integer (kind=4) :: J2max1, J2max1out ! replace bigJmax  ! !2*Jmax+1
	logical :: isOdd
	integer (kind=4) :: numOfJ,numOfJout  !default is int(Jmax+1.)
	integer (kind=4) :: numOfJObs         ! for scalar observables
	integer (kind=4) :: numOfM
		  
!.............. IMPORTANT ADDITION TO 1.3.0: CHOOSING THE POSSIBLE Jvalues
	real, allocatable :: xJlist(:)         !list of chosen 2xJ values
	integer, allocatable :: mapJ2indx(:)   ! maps 2xJ to the index
	integer, allocatable :: startJ2list(:) ! list of starting indices

!.................. MESHES OF EULER ANGLES..............................		  
	integer(kind=4) :: numOfBeta
	real (kind=8), allocatable :: alpha_i(:), gamma_k(:), beta_j(:)  ! meshes of Euler angles
	integer(kind=4), allocatable :: prime(:)    !For symmetry, added by Changfeng Jiao
		  
!........ FUNDAMENTAL SAMPLES OF MESHES
	complex (kind=8), allocatable :: N_ijk(:,:,:,:,:), H_ijk(:,:,:,:,:), PN_ijk(:,:,:,:,:), PH_ijk(:,:,:,:,:)

!........... NEW ARRAYS SET UP IN 1.4.4 MAY 2019..................		  
	complex (kind=8), allocatable :: Nijk(:,:,:), Hijk(:,:,:), PNijk(:,:,:), PHijk(:,:,:)

!............... ARRAYS TO INVERT ALPHA, GAMMA SAMPLING
	complex(kind=8), allocatable :: Zmat(:,:)  ! used for inverting alpha, gamma		  
	complex(kind=8), allocatable :: Zmatinv(:,:) ! added by Changfeng Jiao

!.......... PARTIALLY TRANSFORMED ARRAYS........
	complex (kind=8), allocatable :: N_jMK(:,:,:,:,:), H_jMK(:,:,:,:,:), PN_jMK(:,:,:,:,:), PH_jMK(:,:,:,:,:)

!........... ARRAY TO INVERT ON BETA.................
	real (kind=8), allocatable :: deltaJpJ(:,:) !  USED FOR TRANSFORMING 

!.............. MODIFIED ARRAYS FOR INVERSION..............
	complex (kind=8), allocatable :: Ntild_Jprime_MK(:,:,:,:,:), Htild_Jprime_MK(:,:,:,:,:)
	complex(kind=8),allocatable :: PNtild_Jprime_MK(:,:,:,:,:), PHtild_Jprime_MK(:,:,:,:,:)

!.................. FINAL ARRAYS...........................
	type (jaggedArray), allocatable,target :: N_J_MK(:,:,:), H_J_MK(:,:,:), PN_J_MK(:,:,:), PH_J_MK(:,:,:)   ! FINAL ARRAYS OF ARRAYS

!..................  ARRAYS OF SCALAR OBS...........................
	logical :: compute_expect
	type (jaggedArray), allocatable,target :: Nobs_J_MK(:,:,:), Hobs_J_MK(:,:,:), PNobs_J_MK(:,:,:), PHobs_J_MK(:,:,:)   ! FINAL ARRAYS OF observables		  

	!............. NOTE ON INDEXING.....
!  for a given J, N_J_MK(iJ)%MK(iM,iK)
!  here iM,iK starts from 1 and goes to 2J+1
!			  

!............. OTHER ARRAYS............................
!	      complex (kind=8), allocatable :: allOvlpp(:,:,:,:), allOvlpn(:,:,:,:)
	integer (kind=4), allocatable :: IPIVOT(:)	
		  
!............ INFORMATION ON DELETED J - VALUES .................
	integer :: ndelJ		  
	logical, allocatable :: chooseme(:)  ! for J-values
!	logical, allocatable :: chooseMval(:) ! for M-values
	integer :: nproblems,npossible
	real(4) :: eigcrit      ! criterion number for improving eigenvalues

	contains

!-----------------------------------------------------------------------
!
! user enters a value for Jmax -> the highest value to project
! also determines J2max1, isOdd, and numOfJ
!
! OUTPUT:
!   Jmax = for projection [user input]
!   J2max1 = 2*Jmax + 1
!   isOdd = logical flag for odd A nuclides
!   numOfJ = # of beta values -> int(Jmax + 1.)
!-----------------------------------------------------------------------
	subroutine inputJmax 
  	implicit none

!............INTERNAL....................
		integer :: AA
		real (kind=8) :: check
		character (len = 12) :: jstring

		INTEGER :: ierr ! for MPI communication

!---------------------- USER ENTERS "Jmax" ------------------------
		AA = numprot + numneut
		if(mod(AA,2) == 0)then
			check = 0.0d0; jstring = 'integer'; isOdd = .false.
		else
			check = 1.0d0; jstring = 'half-integer'; isOdd = .true.
		endif ! MOD(AA,2) == 0

		IF (myMPIrank == root) THEN 
			write(*,*) 'Enter J-max to project (',trim(jstring),' value):'
			do
				read(*,*) Jmax
				if(mod(2.0d0*Jmax,2.0d0) == check)then
					exit
				else
					write(*,*) 'Incorrect j.  A = ', AA
					write(*,*) 'Expecting ', trim(jstring), ' values for J'
					write(*,*) 'Please re-enter J-max.'
				end if ! mod
			end do !
		END IF ! myMPIrank == root

		CALL MPI_BARRIER(icomm,ierr)
		CALL MPI_BCAST(Jmax,1,MPI_INT,root,icomm,ierr)

!--------------------- SET "numOfJ" & "J2max1" ----------------------
		numOfJ = int(Jmax + 1.)
		J2max1 = int(2.*Jmax)+1
		numOfJout = numOfJ
		J2max1out = J2max1   ! used for output; may be modified later

		RETURN
	END SUBROUTINE inputJmax

!-----------------------------------------------------------------------
! set up default "mesh" of J,M,K values:
!    J runs from either 0 or 1/2 up to some Jmax;
!    M,K both run from -Jmax to + Jmax
!-----------------------------------------------------------------------
	subroutine default_Jmesh
		implicit none
		real :: xJ
		integer :: ij
	
		if(allocated(xJlist))deallocate(xJlist)
		if(allocated(mapJ2indx))deallocate(mapJ2indx)
		if(allocated(startJ2list))deallocate(startJ2list)
		allocate(xJlist(numOfJ),mapJ2indx(0:numOfJ*2),startJ2list(0:numOfJ*2))
	
		if(isOdd)then
			xJ = 0.5
		else
			xJ = 0.
		end if
	
		xJ = xJ-1.0
		do ij = 1,numOfJ
			xJ = xJ +1.0
			xJlist(ij)= xJ
			mapJ2indx(nint(2.0*xJ))=ij
			startJ2list(nint(2.0*xJ))=ij
		end do ! ij
	
		if(xJ /= Jmax .AND. myMPIrank == root)then 
			print*,' Huh. Somehow did not set up J mesh correctly '
			print*,Jmax,xJ
			stop
		end if
	
		if(allocated(chooseme))deallocate(chooseme)
		allocate(chooseme(numOfJ))
		chooseme =.true.

		RETURN
	END SUBROUTINE default_Jmesh

!-----------------------------------------------------------------------
!  routine to set up mesh of points for sampling Euler angles alpha, beta, gamma
!  default is to set them equally spaced
!-----------------------------------------------------------------------
	subroutine set_euler_meshes
		implicit none
		integer :: ii,iii, nu
		real (kind=4), parameter :: pi = acos(-1.)
		logical, allocatable :: skip(:)
		real :: rv
	
		if(allocated(alpha_i))deallocate(alpha_i)
		if(allocated(beta_j))deallocate(beta_j)
		if(allocated(gamma_k))deallocate(gamma_k)
		if(allocated(prime))deallocate(prime)
	
!-------------------- POPULATE ALPHA(i) & GAMMA(k) ----------------
		allocate(alpha_i(1:J2max1))
		allocate(gamma_k(1:J2max1))
		allocate(prime(1:J2max1))

! Edited by Changfeng Jiao
		if (MOD(J2max1,4) .eq. 1) then
			nu = INT(Jmax)
			do ii = 1, nu
				gamma_k(ii) = float(ii) * pi / (float(nu) + 1.)
				prime(ii) = nu + 1 - ii
			end do ! ii 

			do ii = nu+1, 2*nu
				gamma_k(ii) = (float(ii) - float(nu))*pi/ (float(nu) + 1.) + pi
				prime(ii) = 3*nu + 1 -ii
			end do ! ii 
			gamma_k(2*nu+1) = pi
			prime(2*nu+1) = 0
		else if(MOD(J2max1,4) .eq. 2) then
			nu = INT(Jmax - 0.5)
			do ii = 1, nu
				gamma_k(ii) = float(ii) * pi / (float(nu) + 1.)
				prime(ii) = nu + 1 - ii
			end do ! ii 

			do ii = nu+1, 2*nu
				gamma_k(ii) = (float(ii) - float(nu)) * pi/ (float(nu) + 1. ) + pi
				prime(ii) = 3*nu + 1 - ii
			end do ! ii 
			gamma_k(2*nu+1) = 0.0
			gamma_k(2*nu+2) = pi
			prime(2*nu+1) = 0 !2*nu+1
			prime(2*nu+2) = 0
		else if (MOD(J2max1, 4) .eq. 3) then
			nu = INT(Jmax - 1.)
			do ii = 1, nu
				gamma_k(ii) = float(ii) * pi / (float(nu) + 1.)
				prime(ii) = nu + 1 - ii
			end do

			do ii = nu+1, 2*nu
				gamma_k(ii) = (float(ii) - float(nu)) * pi/ (float(nu) + 1. ) + pi
				prime(ii) = 3*nu + 1 - ii
			end do

			gamma_k(2*nu+1) = 0.0
			gamma_k(2*nu+2) = pi / 2.0
			gamma_k(2*nu+3) = pi
			prime(2*nu+1) = 0
			prime(2*nu+2) = 2*nu+2
			prime(2*nu+3) = 0
		else if (MOD(J2max1, 4) .eq. 0) then
			nu = INT(Jmax + 0.5)
			do ii = 1, nu
				gamma_k(ii) = float(ii) * pi / (float(nu) + 1.)
				prime(ii) = nu + 1 - ii
			end do

			do ii = nu+1, 2*nu
				gamma_k(ii) = (float(ii) - float(nu)) * pi / (float(nu) + 1. ) + pi
				prime(ii) = 3*nu + 1 - ii
			end do
		else
			stop "Something wrong with the Jmax"
		end if

		do ii = 1, J2max1
			alpha_i(ii) = gamma_k(ii)
		end do

!---------------------- POPULATE BETA(j) --------------------------
		allocate(beta_j(1:numOfJ))

!.............. 	SET UP ARRAY FOR SKIPPING ANGLES .......
		allocate(skip(numOfJ+ndelJ)) 
		skip(:)= .false.

		ii = 0
		do iii = 1, numOfJ+ndelJ
			if(ndelJ > 0 .and. .not.chooseme(iii))cycle
			ii = ii + 1
			if(ii > numOfJ)cycle
			beta_j(ii) = (float(iii) -.5)*pi/(float(numOfJ+ndelJ))
		end do

		if(ii /= numOfJ)print*,' WHOOPS DID NOT ADD UP ',ii,numOfJ
		print*,' '
		IF (myMPIrank == root) PRINT*,' Number of ANGLES in Euler mesh = ',numOfJ*J2max1*J2max1
!	 print*,beta_j
		RETURN
	END SUBROUTINE set_euler_meshes		

!-----------------------------------------------------------------------
! this routine follows SET_EULER_MESHES
! but unrolls the loops into one big list
!-----------------------------------------------------------------------
	SUBROUTINE set_euler_list
		IMPLICIT NONE 

		RETURN
	END SUBROUTINE set_euler_list  

!-----------------------------------------------------------------------
! set up arrays for basic sampling
!
! N_ijk = < Psi | R(alpha_i,beta_j, gamma_k) | Psi >
!  where R is the standard rotation matrix
!
! H_ijk = < Psi | R(alpha_i,beta_j, gamma_k) H | Psi >
!
! PN_ijk = < Psi | R(alpha_i,beta_j, gamma_k) P | Psi >
! where P = parity operator
!
! PH_ijk = < Psi | R(alpha_i,beta_j, gamma_k) H P | Psi >
!-----------------------------------------------------------------------
	SUBROUTINE allocate_arrays4sampling
		IMPLICIT NONE

		INTEGER :: temperr 
		
		allocate(N_ijk(numsd,numsd,J2max1,numOfJ,J2max1),STAT=temperr)
		! PRINT*, ' Node = ', myMPIrank, ' allsameparity = ', allsameparity ! TESTING, REMOVE
		! PRINT*, ' Node = ', myMPIrank, ' temperr 4 sampling = ', temperr 	! TESTING, REMOVE
    if(.not.allsameparity)allocate(PN_ijk(numsd,numsd,J2max1,numOfJ,J2max1))
		if(doHam)then
			allocate(H_ijk(numsd,numsd,J2max1,numOfJ,J2max1))	
			if(.not.allsameparity)allocate(PH_ijk(numsd,numsd,J2max1,numOfJ,J2max1))
		end if

		RETURN
	END SUBROUTINE allocate_arrays4sampling		  

!-----------------------------------------------------------------------
!  set up the final arrays
!  uses derived types ('jagged arrays')
!
!  want: N^J_MK, H^J_MK for both even and odd parities
!-----------------------------------------------------------------------
	subroutine allocate_X_J_MK
		implicit none
	
		integer :: ij
		real    :: xJ
		integer :: sizeOfMK
		integer :: sdX,sdY   ! indices for multiple SDs
	
		allocate(N_J_MK(numsd,numsd,numOfJ))
		if(.not.allsameparity)	allocate(PN_J_MK(numsd,numsd,numOfJ))
		if(doHam)then
			allocate(H_J_MK(numsd,numsd,numOfJ))
			if(.not.allsameparity)	allocate(PH_J_MK(numsd,numsd,numOfJ))
		end if
	
		do ij = 1,numOfJ
			xJ = xJlist(ij)
			sizeOfMK = nint(2*xJ)+1		
			do sdX = 1,numsd
				do sdY = 1,numsd
					N_J_MK(sdX,sdY,ij)%dimMK = sizeOfMK
					allocate(N_J_MK(sdX,sdY, ij)%MK(sizeOfMK,sizeOfMK))
					N_J_MK(sdX,sdY, ij)%MK= (0.d0,0.d0)
					if(.not.allsameparity)then
						PN_J_MK(sdX,sdY,ij)%dimMK = sizeOfMK
						allocate(PN_J_MK(sdX,sdY,ij)%MK(sizeOfMK,sizeOfMK))
						PN_J_MK(sdX,sdY,ij)%MK=(0.d0,0.d0)
					end if
					if(doHam)then
						H_J_MK(sdX,sdY,ij)%dimMK = sizeOfMK
						allocate(H_J_MK(sdX,sdY,ij)%MK(sizeOfMK,sizeOfMK))
						H_J_MK(sdX,sdY,ij)%MK=(0.d0,0.d0)
						if(.not.allsameparity)then
							PH_J_MK(sdX,sdY,ij)%dimMK = sizeOfMK
							allocate(PH_J_MK(sdX,sdY,ij)%MK(sizeOfMK,sizeOfMK))	
							PH_J_MK(sdX,sdY,ij)%MK=(0.d0,0.d0)
						end if		
					end if
				end do !sdY
			end do ! sdX
		end do !ij
	
		RETURN
	END SUBROUTINE allocate_X_J_MK

!-----------------------------------------------------------------------
! creates the Z matrix for inverting alpha, gamma indices
!  replaces subroutine populateZeta which had only vectors
!
! [equation 17]
!  NOTE: assume distribution of alpha, gamma meshes are the same! very important
!
! Z_mi = 1/(2Jmax+1) exp(-i M_m alpha_i)
!-----------------------------------------------------------------------
	subroutine populateZmat
		implicit none
!............INTERNAL....................
		integer (kind=4) :: ik,jk, info
		real(kind=8)  :: xJ,xK

		INTEGER :: temperr 
	
!	integer :: kk
!	complex(kind=8) :: ztmp
!    complex (kind=8),allocatable :: Amat(:,:) 

		if(allocated(Zmat))deallocate(Zmat)
		allocate(Zmat(J2max1,J2max1),STAT=temperr)
		!PRINT*, ' Node = ', myMPIrank, ' Zmat allocation temperr = ', temperr 		! TESTING, REMOVE
		if(allocated(Zmatinv))deallocate(Zmatinv)
		allocate(Zmatinv(J2max1,J2max1),STAT=temperr)
		!PRINT*, ' Node = ', myMPIrank, ' Zmatinv allocation temperr = ', temperr	! TESTING, REMOVE

		xJ = dble(J2max1-1)*0.5d0
		xK = -xJ

		! PRINT*, ' Node = ', myMPIrank, ' alpha_i(1) = ', alpha_i(1)
		do ik = 1, J2max1
			do jk = 1,J2max1			 
				Zmat(jk,ik) = exp((0.d0,1.d0)*xK*alpha_i(jk))
			end do ! jk
			xK = xK + 1.d0
		end do ! ik

		! PRINT*, ' Node = ', myMPIrank, ' Zmat(1,1) = ', Zmat(1,1)

		call complex_inverse(J2max1,Zmat,Zmatinv,info)

!			ztmp = (0.d0,0.d0)
!			do kk = 1,j2max1
!				ztmp = ztmp + zmatinv(ik,kk)*zmat(kk,jk)
!			end do
!			print*,ik,jk,ztmp			
!		end do
!	end do
		return
	end subroutine populateZmat

!-----------------------------------------------------------------------
! transforms N_ijk -> N_jMK etc
!
! here N_jMK = sum Z_mi Z_Kk N_ijk  etc
!-----------------------------------------------------------------------
	subroutine invert_alpha_and_gamma
		implicit none

		integer iJ,jJ,iM,iK,J
		integer :: ialp,kgam
		real(kind=8):: xM, xK
		complex(kind=8), allocatable :: ztmp1(:),ztmp2(:),ztmp3(:),ztmp4(:)
		integer :: sdX,sdY   ! indices for multiple SDs	
	
		if(allocated(N_jMK))deallocate(N_jMK)
		if(allocated(PN_jMK))deallocate(PN_jMK)
		if(allocated(H_jMK))deallocate(H_jMK)
		if(allocated(PH_jMK))deallocate(PH_jMK)
	
		allocate(N_jMK(numsd,numsd,numOfJ,J2max1,J2max1))
		if(.not.allsameparity)allocate(PN_jMK(numsd,numsd,numOfJ,J2max1,J2max1))
		if(doHam)then
			allocate(H_jMK(numsd,numsd,numOfJ,J2max1,J2max1))
			if(.not.allsameparity)allocate(PH_jMK(numsd,numsd,numOfJ,J2max1,J2max1))
		end if
	
		allocate(ztmp1(numOfJ),ztmp2(numOfJ),ztmp3(numOfJ),ztmp4(numOfJ))    
	
		do sdX = 1,numsd
			do sdY = sdX,numsd
				do iM = 1,J2max1		
					do iK = 1,J2max1
						ztmp1(:) = (0.d0,0.d0)
						ztmp2(:) = (0.d0,0.d0)
						ztmp3(:) = (0.d0,0.d0)
						ztmp4(:) = (0.d0,0.d0)
						do ialp = 1,J2max1
							do kgam = 1,J2max1
								do J = 1,numOfJ
! 						ztmp1(j)=ztmp1(j)+ Zmatinv(kgam,im)*Zmatinv(ialp,ik)*N_ijk(sdX,sdY,ialp,j,kgam)
									ztmp1(j)=ztmp1(j)+ Zmatinv(iM,kgam)*Zmatinv(iK,ialp)*N_ijk(sdX,sdY,ialp,j,kgam)
								end do ! J 
								if(.not.allsameparity)then
									do J = 1,numOfJ
										ztmp2(j)=ztmp2(j)+ Zmatinv(iM,kgam)*Zmatinv(iK,ialp)*PN_ijk(sdX,sdY,ialp,j,kgam)
									end do ! J 
								end if ! .not.allsameparity
								if(doHam)then
									do J = 1,numOfJ
										ztmp3(j)=ztmp3(j)+ Zmatinv(iM,kgam)*Zmatinv(iK,ialp)*H_ijk(sdX,sdY,ialp,j,kgam)
									end do
									if(.not.allsameparity)then
										do J = 1,numOfJ
											ztmp4(j)=ztmp4(j)+ Zmatinv(iM,kgam)*Zmatinv(iK,ialp)*PH_ijk(sdX,sdY,ialp,j,kgam)
										end do
									end if						
								end if
							end do ! kgam
						end do ! ialp
						do J = 1,NumOfJ
							N_jMK(sdX,sdY,j,iM,iK)= ztmp1(j)
							if(.not.allsameparity)then
								PN_jMK(sdX,sdY,j,iM,iK)= ztmp2(j)
							end if
 							if(doHam) THEN
								H_jMK(sdX,sdY,j,iM,iK)= ztmp3(j)
								if(.not.allsameparity)then
									PH_jMK(sdX,sdY,j,iM,iK)= ztmp4(j)
								end if
							end if
						end do ! J 
					end do ! iK 
				end do ! iM
			end do ! sdY 
		end do ! sdX
		deallocate(ztmp1,ztmp2,ztmp3,ztmp4)

		return
	end subroutine invert_alpha_and_gamma

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
	subroutine allocate_delta
		INTEGER :: temperr 
		INTEGER :: ierr
		if(allocated(deltaJpJ))deallocate(deltaJpJ)
		allocate(deltaJpJ(numOfJ,numOfJ), STAT=temperr)

		deltaJpJ(:,:) = 0.d0

		return
	end subroutine allocate_delta

!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------

	subroutine allocate_tilded_matrices

		if(allocated(Ntild_Jprime_MK))deallocate(Ntild_Jprime_MK)
		if(allocated(PNtild_Jprime_MK))deallocate(PNtild_Jprime_MK)
		if(allocated(Htild_Jprime_MK))deallocate(Htild_Jprime_MK)
		if(allocated(PHtild_Jprime_MK))deallocate(PHtild_Jprime_MK)

		allocate(Ntild_Jprime_MK(numsd,numsd,numOfJ,J2max1,J2max1))
		if(.not.allsameparity)allocate(PNtild_Jprime_MK(numsd,numsd,numOfJ,J2max1,J2max1))
		if(doHam)then
			allocate(Htild_Jprime_MK(numsd,numsd,numOfJ,J2max1,J2max1))
			if(.not.allsameparity)allocate(PHtild_Jprime_MK(numsd,numsd,numOfJ,J2max1,J2max1))
		end if ! doHam

		return
	end subroutine allocate_tilded_matrices

!-----------------------------------------------------------------------
!
! computes, for fixed K,M, deltaJ'J
!  = sum_j  d^J'_MK(beta_j)*d^J_MK(beta_j)
!
! both J and J' go from max(|K|,|M|) to maxJ
!
! where d is the Wigner little-d function
!
!  CALLED BY:
!   projectorator
!   testdeltas
!-----------------------------------------------------------------------
	subroutine compute_delta(intM,intK)	
		use wignersfriend
		implicit none

		integer :: intK, intM
		integer :: aJ,bJ,j
		real(8) :: xM,xK,xJ,xJp
		real(4) :: Jmin
		integer :: localNofJs
		real(kind=8) :: dtmp
		integer :: jstart
		complex(kind=8), allocatable :: ztmpmat(:,:)
		real (kind = 8),ALLOCATABLE :: rwork(:),bb(:)
		COMPLEX (kind = 8), ALLOCATABLE :: work(:)
		integer :: lwork,info
		real, parameter :: etol = 0.001
	
		deltaJpJ(:,:)=0.d0

!........... RECONSTRUCT LIMITS...............
		xM = float(intM)-(Jmax+1.)
		xK = float(intK)-(Jmax+1.)
		
		Jmin = max(abs(xM),abs(xK))
		Jstart = startJ2list(nint(2*Jmin))
		localNofJs = numOfJ-Jstart+1

		allocate(ztmpmat(localNofJs,localNofJs))
		do aJ = Jstart,numOfJ
			xJp = xJlist(aJ)
			do bJ = Jstart,numOfJ
				xJ = xJlist(bJ)
				dtmp = 0.d0			
				do j = 1, numOfJ
					dtmp = dtmp + wigner_d(xJp,xM,xK,real(beta_j(j),kind=8))* &
                            wigner_d(xJ,xM,xK,real(beta_j(j),kind=8))
				end do ! j 
				deltaJpJ(aJ,bJ)=dtmp
				ztmpmat(aJ-Jstart+1,bJ-Jstart+1)=cmplx(dtmp)
			end do ! bJ 
		end do  ! aJ 

!....... CHECK EIGENVALUES.......a bit kludgy
		lwork = 2*localNofJs - 1
		ALLOCATE(work(lwork),bb(localNofJs))
		ALLOCATE(rwork(3*localNofJs - 2))
		CALL zheev('n','u',localNofJs,ztmpmat,localNofJs,bb,work,lwork,rwork,info)
	
		do j = 1,localNofJs
			if(bb(j)< etol)then
				nproblems =nproblems+1
				eigcrit = eigcrit+bb(j)
			end if
			npossible = npossible+1
		end do

		deallocate(bb,ztmpmat,work,rwork)
	
		return
	end subroutine compute_delta

!-----------------------------------------------------------------------
! subroutine to compute \tilde{N}^Jprime_MK
! inspired by SVD
!
!  here ~N^J'_MK = sum d^J'_MK (beta_j)  N_jMK
! 
!  for fixed M,K  (important!)
!-----------------------------------------------------------------------
	subroutine compute_tilded_matrices(intM,intK)
		use wignersfriend

		implicit none
	
		integer :: intK, intM
		integer :: aJ,bJ,ij
		real(8) :: xM,xK,xJ,xJp
		real(4) :: Jmin
		integer :: localNofJs
		complex(kind=8) :: ztmp1,ztmp2,ztmp3,ztmp4
		real(4) :: wigtmp
		integer :: jstart
		integer :: sdX,sdY   ! indices for multiple SDs
	
!........... RECONSTRUCT LIMITS...............
		xM = float(intM)-(Jmax+1.)
		xK = float(intK)-(Jmax+1.)
	
		Jmin = max(abs(xM),abs(xK))
		Jstart = startJ2list(nint(2*Jmin))
		do sdX = 1,numsd
			do sdY = sdX,numsd
				do aJ = Jstart,numOfJ
					xJp = xJlist(aJ)
					ztmp1 =(0.d0,0.d0)
					ztmp2 =(0.d0,0.d0)
					ztmp3 =(0.d0,0.d0)
					ztmp4 =(0.d0,0.d0)
					do iJ = 1,numOfJ
						wigtmp = wigner_d(xJp,xM,xK,real(beta_j(ij),kind=8))
						ztmp1 = ztmp1 + wigtmp * N_jMK(sdX,sdY,ij,intM,intK)
						if(.not.allsameparity)ztmp2 = ztmp2 + wigtmp * PN_jMK(sdX,sdY,ij,intM,intK)
						if(doHam)then
							ztmp3 = ztmp3 + wigtmp * H_jMK(sdX,sdY,ij,intM,intK)
							if(.not.allsameparity)ztmp4 = ztmp4 + wigtmp * PH_jMK(sdX,sdY,ij,intM,intK)
						end if ! doHam
					end do ! iJ
!		if(xJp==1.0)print*,xJp,intM,intK,ztmp1
					Ntild_Jprime_MK(sdX,sdY,aJ,intM,intK)= ztmp1
					if(.not.allsameparity)PNtild_Jprime_MK(sdX,sdY,aJ,intM,intK)= ztmp2
					if(doHam)then
						Htild_Jprime_MK(sdX,sdY,aJ,intM,intK)= ztmp3
						if(.not.allsameparity)PHtild_Jprime_MK(sdX,sdY,aJ,intM,intK)= ztmp4		
					end if
				end do	! aJ
			end do	! sdY
		end do	!sdX
	
		return
	end subroutine compute_tilded_matrices

!-----------------------------------------------------------------------
!  solve sum_J Delta^J'J N^J_MK = ~N^J'_MK
!  for fixed M,K
!  also find H^J_MK etc
!
! CALLS lapack routine ZGESV to carry out LU decomposition
! and solve linear equation; if more than one solution needed,
! addition solutions invoke lapack routine ZGETRS  
! (so don't need to redo LU decomposition)
!-----------------------------------------------------------------------
	subroutine solve4projected_matrices(intM,intK)
!	    use system_parameters,only:numsd
		implicit none
!............... INPUT...............................

		integer :: intK, intM

!............... INTERNAL................		
		integer :: sizeOfNorm
		integer :: info
		complex(kind=8),allocatable :: ztmpvec(:)
		real(4) :: Jmin,xJ
		real(4) :: xM,xK
		integer(4) :: Jstart,Jshift
		complex(kind=8), allocatable :: ztmpmat(:,:)
		real (kind = 8),ALLOCATABLE :: rwork(:),bb(:)
		COMPLEX (kind = 8), ALLOCATABLE :: work(:)
		integer :: lwork		
		integer :: i,j
		integer :: sdX, sdY

!......... COMPUTE DIMENSION...................		
		xM = float(intM)-(Jmax+1.)
		xK = float(intK)-(Jmax+1.)
		Jmin = max(abs(xM),abs(xK))
		Jstart = startJ2list(nint(2*Jmin))
		sizeOfNorm = numOfJ -Jstart +1
!		print*,xM,xK,sizeOfNorm,' size of norm ',Jmin,Jstart,numOfJ

!.......... COPY OVER TO ARRAYS
		allocate(ztmpvec(sizeOfNorm),ztmpmat(sizeOfNorm,sizeOfNorm))
		if(allocated(IPIVOT))deallocate(IPIVOT)
		allocate(IPIVOT(sizeOfNorm))
		xM = float(intM)-(Jmax+1.)
		xK = float(intK)-(Jmax+1.)
!		print*,xM,xK

		do sdX = 1,numSD
			do sdY = sdX,numSD
				do i = 1,sizeOfNorm
					ztmpvec(i)=Ntild_Jprime_MK(sdX,sdY,i+Jstart-1,intM,intK)
					do j = 1,sizeOfNorm
						ztmpmat(i,j)=cmplx(deltaJpJ(i+Jstart-1,j+Jstart-1),0.d0)
					end do ! j 
				end do ! i 
				call ZGESV(sizeOfNorm,1,ztmpmat,sizeOfNorm,IPIVOT,ztmpvec,sizeOfNorm,INFO)
				call errorZGESV(INFO)
				if(info/=0)then
					print*,' error in LU decomposition in solve4projected_matrices '
					stop
				end if
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
					N_J_MK(sdX,sdY,i+Jstart-1)%MK(intM-Jshift,intK-Jshift)=ztmpvec(i)
				end do ! i 
				if(.not.allsameparity)then
					do i = 1,sizeOfNorm
						ztmpvec(i)=PNtild_Jprime_MK(sdX,sdY,i+Jstart-1,intM,intK)
					end do	
					call ZGETRS('No transpose',sizeOfNorm,1,ztmpmat,sizeOfNorm,IPIVOT,ztmpvec,sizeOfNorm,INFO)
					call errorZGESV(INFO)
					do i = 1,sizeOfNorm
						xJ = xJlist(i+Jstart-1)
						Jshift = nint(Jmax-xJ)
						PN_J_MK(sdX,sdY,i+Jstart-1)%MK(intM-Jshift,intK-Jshift)=ztmpvec(i)
					end do ! i 
				end if ! .not.allsameparity
			end do ! sdY
		end do ! sdX

		if(.not.doHam)return

!................. SOLVE FOR HAMILTONIAN MATRICES.............
		do sdX = 1,numSD
			do sdY = sdX,numSD
				do i = 1,sizeOfNorm
					ztmpvec(i)=Htild_Jprime_MK(sdX,sdY,i+Jstart-1,intM,intK)
				end do		
				call ZGETRS('No transpose',sizeOfNorm,1,ztmpmat,sizeOfNorm,IPIVOT,ztmpvec,sizeOfNorm,INFO)
				call errorZGESV(INFO)
				do i = 1,sizeOfNorm
					xJ = xJlist(i+Jstart-1)
					Jshift = nint(Jmax-xJ)
					H_J_MK(sdX,sdY,i+Jstart-1)%MK(intM-Jshift,intK-Jshift)=ztmpvec(i)
				end do
				if(.not.allsameparity)then
					do i = 1,sizeOfNorm
						ztmpvec(i)=PHtild_Jprime_MK(sdX,sdY,i+Jstart-1,intM,intK)
					end do	
					call ZGETRS('No transpose',sizeOfNorm,1,ztmpmat,sizeOfNorm,IPIVOT,ztmpvec,sizeOfNorm,INFO)
					call errorZGESV(INFO)
					do i = 1,sizeOfNorm
						xJ = xJlist(i+Jstart-1)
						Jshift = nint(Jmax-xJ)
						PH_J_MK(sdX,sdY,i+Jstart-1)%MK(intM-Jshift,intK-Jshift)=ztmpvec(i)
					end do			
				end if		
			end do   ! sdY
		end do   ! sdX
		deallocate(IPIVOT,ztmpvec,ztmpmat)
		return
	end subroutine solve4projected_matrices

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
	subroutine deallocator
		use tracy
		implicit none
		if(allocated(N_ijk))deallocate(N_ijk)
		if(allocated(PN_ijk))deallocate(PN_ijk)
		if(allocated(H_ijk))deallocate(H_ijk)
		if(allocated(PH_ijk))deallocate(PH_ijk)
	
!	if(allocated(allovlpp))deallocate(allovlpp)
!	if(allocated(allovlpn))deallocate(allovlpn)
	
		if(allocated(N_J_MK))deallocate(N_J_MK)
		if(allocated(PN_J_MK))deallocate(PN_J_MK)
		if(allocated(H_J_MK))deallocate(H_J_MK)
		if(allocated(PH_J_MK))deallocate(PH_J_MK)
	
		if(allocated(normtrace))deallocate(normtrace)
		if(allocated(pnormtrace))deallocate(pnormtrace)
		if(allocated(hamtrace))deallocate(hamtrace)
		if(allocated(phamtrace))deallocate(phamtrace)
	
		return
	end subroutine deallocator

!-----------------------------------------------------------------------
!  creates a matrix with the parities of each s.p. state
!
! OUTPUT:  parOp(i) = +/-1 for each s.p. state i
!
! CALLED BY
!-----------------------------------------------------------------------
	SUBROUTINE MakeParityOp(parOP)

		USE spstate
		USE psis

		IMPLICIT NONE

		COMPLEX (KIND = 8), DIMENSION(nsps,nsps), INTENT(OUT) :: parOP
		INTEGER :: i

		parOP = 0.0d0
		DO i = 1, nsps
			parOP(i,i) = (-1.0d0)**(spsqn(i)%l)*(1.d0,0.d0)
		END DO

		return
END SUBROUTINE MakeParityOp

!-----------------------------------------------------------------------
!   Checks for errors in ZGESV subroutine.
!
!  INPUT:
!   INFO = contains error flag i.e. if INFO /= 0 then there is an issue
!-----------------------------------------------------------------------
	subroutine errorZGESV(INFO)
		implicit none

!........... INPUT........................
		integer, intent(in) :: INFO

!----------------------- CHECK FOR ERORS --------------------------
		if (INFO > 0) then
			write(*,*) 'ERROR, factorization completed, but the factor U is '
			write(*,*) 'exactly singular, so the solution could not be computed.'
		else if (INFO < 0) then
			write(*,'(A,I3.1)') 'ERROR, illegal value at element ', INFO
		else
			continue
		end if

		return
	end subroutine errorZGESV		  
!==================================================
end module lamplight	 
	

