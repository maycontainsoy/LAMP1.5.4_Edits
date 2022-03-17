module eigenpackage
	USE nodeinfo

	implicit none

!..... ADDED IN 1.3.6...
    integer :: nlevels
	real(kind=8), allocatable :: energylevels(:),jlevels(:)
!.... ADDED IN 1.4.4.....	
	real(kind=8), allocatable :: obslevels(:)
	character, allocatable :: paritylevels(:)
		
	real :: jtarget
	
	integer:: nparity_soln ! # of parity solutions

contains

!==============================================================================
!
! takes the solved Norm and Hamiltonian and runs them through "geneigsolverdiag" to eigensolve.
! after, prints out results to the screen.
!
!  INPUT:
!
!   tolerance = cuttoff point for eliminating zeroes in the Norm matrix (during zhegvdiag)
!   npmax = # of levels possible
!   numsdused = # of sds used in calculation--can be changed
!   nlist = length of problist, hamlist, and jlist; int(Jmax-Jmin)+1; usually numOfJ
!            that is, it is the number of possible J values
!
! OUTPUT:
!   nftotal = # of total states FOUND 
!   jall
!   pallPair
!   
!   problist(2,:) = fraction of original HF state in the norm for each J
!   hamlist
!   jlist
!   normsum
!   hamsum
!   
!
!  INTERNAL:
!   ParityEvals(i=1,numberFound) = energy eigenvalues found
!   normvals = not used
!   parityOUT = counts parity state (i.e. 1 or 2)
!   npair = number of parity states (i.e. 1 or 2) [via phf_vals module]
!   TSDHmat = temporary Slater determinant Hamiltonian matrix (holds values from H_J_MK)
!   TSDNmat = temporary Slater determinant Norm matrix (holds values from N_J_MK)
!   nf = number of states found
!
subroutine EigenSolverPackage(tolerance,npmax, numsdused, nftotal,jfound,Efound,obsfound,normsum,hamsum,problist,hamlist)
	use jaggedArrayType
	use lamplight
	use tracy
	use psis,only: numsd
	use dens1body

	implicit none

	real (kind=4), intent(in) :: tolerance
	integer (kind=8) :: npmax!!
	integer,intent(in) :: numsdused

!............OUTPUT......................
	integer (kind=8), intent(out) :: nftotal!!
	real (kind=8), intent(inout) :: jfound(npmax), Efound(2,npmax),obsfound(2,npmax) !jall(npmax), pallPair(2,npmax)!!
	real (kind = 8) :: normSum, hamSum
	real (kind=8) :: problist(2,numOfJ), hamlist(2,numOfJ)!!

!............INTERNAL....................
	integer (kind=4) :: intJ, kk,ll,ii,jj
	real (kind=4) :: floatJ
	real (kind=8), allocatable :: ParityEvals(:,:), normvals(:,:),obsvals(:,:)
	integer (kind=4) :: parityOUT
	complex(kind=8), allocatable, dimension(:,:) :: TSDHmat, TSDNmat, obsSDHmat
	integer (kind=8) :: nf, nkeep
	integer (kind=4) :: stateFound
	complex(kind=8) :: localtrace
	integer :: sdX,sdY  ! indices for muliple Slater Determinants
	integer :: localdim,matdim  ,a,b
	
!................. ADDED IN 1.3.6.... allocate memory ....................
	nlevels = 0
	do intJ = 1, numOfJ
		floatJ = xJlist(intJ)
		nlevels = nlevels + 2*nint(floatJ)+1
	end do
	nlevels = nlevels*numsdused
	if(.not.allsameparity)nlevels = 2*nlevels
	if(.not. allocated(energylevels))allocate(energylevels(nlevels),jlevels(nlevels),paritylevels(nlevels),obslevels(nlevels))
	
	energylevels = 0.0
	jlevels = 0.0
	paritylevels='+'
	obslevels = 0.0
	nlevels = 0

!....................................
	nftotal = 0!!
	jfound  = 0.d0
	Efound = 0.0d0!!
	normSum = 0.0d0
	hamSum = 0.0d0
	problist = 0.0d0
	hamlist = 0.0d0
	
	if(allsameparity)then
		nparity_soln=1
	else
		nparity_soln=2
	end if
	
!......... ADDED in 1.4.8 for densities
	if(finddense)then
		allocate(solutions(nparity_soln,numOfJ))	
	end if

!...........................	
	do intJ = 1, numOfJ
		floatJ = xJlist(intJ)
		matdim = N_J_MK(1,1,intJ)%dimMK
		localdim = matdim*numsdused
		allocate(ParityEvals(nparity_soln,localdim))
		allocate(normvals(nparity_soln,localdim))
		allocate(TSDHmat(localdim,localdim))
		allocate(TSDNmat(localdim,localdim))
		if(compute_expect)allocate(obsSDHmat(localdim,localdim),obsvals(nparity_soln,localdim))
		nkeep = 0
		ParityEvals = 0.0d0

		do parityOUT = 1, nparity_soln
			TSDHmat = cmplx(0.0d0,0.0d0)
			TSDNmat = cmplx(0.0d0,0.0d0)
			do sdX = 1,numsdused
				do ii = 1,matdim
					a = ii+ (sdX-1)*matdim
					do sdY = sdX,numsdused
						do jj = 1,matdim
							b = jj + (sdY-1)*matdim
							if (.not.allsameparity) then
								TSDHmat(a,b) = & 
								0.5d0*(H_J_MK(sdX,sdY,intJ)%MK(ii,jj) + (-1.0d0)**(parityOUT+1)*PH_J_MK(sdX,sdY,intJ)%MK(ii,jj))
								TSDNmat(a,b) = & 
								0.5d0*(N_J_MK(sdX,sdY,intJ)%MK(ii,jj) + (-1.0d0)**(parityOUT+1)*PN_J_MK(sdX,sdY,intJ)%MK(ii,jj))
								if(compute_expect)then
									obsSDHmat(a,b)=0.5d0*(Hobs_J_MK(sdX,sdY,intJ)%MK(ii,jj)+(-1.0d0)**(parityOUT+1)*PHobs_J_MK(sdX,sdY,intJ)%MK(ii,jj))
								end if
							else ! .not.allsameparity 
								TSDHmat(a,b) = H_J_MK(sdX,sdY,intJ)%MK(ii,jj)
								TSDNmat(a,b) = N_J_MK(sdX,sdY,intJ)%MK(ii,jj) 	
								if(compute_expect)obsSDHmat(a,b)=Hobs_J_MK(sdX,sdY,intJ)%MK(ii,jj)					
							end if ! .not.allsameparity 
!------------ ENFORCE HERMITICIY ------------------------------------------							
							if(sdY > sdX)then
								TSDNmat(b,a)=dconjg(TSDNmat(a,b))
								TSDHmat(b,a)=dconjg(TSDHmat(a,b))
								if(compute_expect)obsSDHmat(b,a)= dconjg(obsSDHmat(a,b))
							end if ! sdy > sdx 
						end do !jj
					end do ! sdy
				end do ! ii
			end do ! sdx

			do kk = 1, localdim
				problist(parityOUT,intJ) = problist(parityOUT,intJ) + dble(TSDNmat(kk,kk))
				hamlist(parityOUT,intJ)  = hamlist(parityOUT,intJ)  + dble(TSDHmat(kk,kk))
			end do
			
			! PROBABLY NEED TO ADD MPI ROOT PROTECTION (BELOW)
			if(localdim < 9)then
				write(66,*)floatJ,' matrices '
				write(66,*)' norm ',localdim
				do ii = 1, localdim
					write(66,'(10f8.3)')(TSDNmat(ii,kk),kk =1,localdim )
				end do
				write(66,*)' ham '
				do ii = 1, localdim
					write(66,'(10f8.3)')(TSDHmat(ii,kk),kk =1,localdim )
				end do
			end if ! localdim < 9 
			! PROBABLY NEED TO ADD MPI ROOT PROTECTION (ABOVE)


!-----------------
			if(compute_expect)then
				call geneigsolverdiag_wexpect(localdim,nparity_soln,parityOUT,TSDNmat,TSDHmat,ObsSDHmat, & 
					tolerance,nf,ParityEvals,normvals,obsvals)
			else
				if(finddense)then   !setting up for transformations for density matrices
					myarrayindex=intJ
					myparity=parityOUT
				end if
				call geneigsolverdiag(localdim,nparity_soln,parityOUT,TSDNmat,TSDHmat,tolerance,nf,ParityEvals,normvals)
			end if

			if (nf > nkeep) nkeep = nf
		end do ! parityOUT

		normsum = normsum + problist(1,intJ)
		hamsum = hamsum + hamlist(1,intJ)
		if (.not.allsameparity) then
			normsum = normsum + problist(2,intJ)
			hamsum = hamsum + hamlist(2,intJ)
		end if

		do stateFound = 1, int(nkeep)
			nftotal = nftotal + 1
			jfound(nftotal) = floatJ!!
			Efound(1,nftotal) = ParityEvals(1,stateFound)!!
			if(compute_expect)obsfound(1,nftotal)=obsvals(1,stateFound)
			if(compute_expect .and. .not.allsameparity)obsfound(2,nftotal)=obsvals(2,stateFound)
			if (.not.allsameparity) Efound(2,nftotal) = ParityEvals(2,stateFound)
			if(.not.allsameparity)then
				if(Efound(1,nftotal)/=0.0)then
					nlevels = nlevels + 1
					energylevels(nlevels)=Efound(1,nftotal)
					jlevels(nlevels)=jfound(nftotal)
					paritylevels(nlevels)='+'
					obslevels(nlevels)=obsfound(1,nftotal)
				end if ! Efound

				if(Efound(2,nftotal)/=0.0)then
					nlevels = nlevels + 1
					energylevels(nlevels)=Efound(2,nftotal)
					jlevels(nlevels)=jfound(nftotal)
					paritylevels(nlevels)='-'
					obslevels(nlevels)=obsfound(2,nftotal)
				end if ! Efound
			else ! .not. allsameparity 
				if(Efound(1,nftotal)/= 0.0)then
					nlevels = nlevels + 1
					energylevels(nlevels)=Efound(1,nftotal)
					jlevels(nlevels)=jfound(nftotal)
					paritylevels(nlevels)='+'
					obslevels(nlevels)=obsfound(1,nftotal)
				end if ! Efound
			end if ! .not. allsameparity
		end do ! stateFound
    
		IF (myMPIrank == root) THEN 
			if (nkeep<1) write(*,*)' No states for J = ',floatJ
		END IF ! myMPIrank == root

		deallocate(ParityEvals,normvals,TSDHmat,TSDNmat)
		if(compute_expect)deallocate(obsSDHmat,obsvals)
	end do ! intJ

	PRINT*, ' Node = ', myMPIrank, ' local dim = ', localdim

1001 format(I3,2X,2F12.5,2X,F4.1)
1002 format(I3,2X,F12.5,2X,F4.1)
	
	return
end subroutine EigenSolverPackage

!==============================================================================
!
! original generalized eigensolution routines
! using diagonalization
!
!  INPUT:
!    ndim = dimension of matrices
!    nparity = # of parity states considered
!    ipar = which parity (1 or 2)
!   Nmat(i,j) = ndim x ndim norm matrix
!   Hmat(i,j)= ndim x ndim hamiltonian matrix
!   tol = real variable, tolerance for eliminating zeros
!
! OUTPUT
!   nf = # of states found
!   pevals(i=1,nf) = energy eigenvalues found
!   normvals(i=1,ndim)
!
! subroutines called:
!
!  zheev  -- lapack routine for diagonalizing a hermitian, complex*16 matrix
!
subroutine geneigsolverdiag(ndim,nparity,ipar,NMat,HMat,tol,nf,pevals,normvals)
!	use errortests
	implicit none
	integer (kind=4) :: ndim, nparity   ! dimension of matrices

    integer ipar
	COMPLEX (KIND = 8):: Hmat(ndim,ndim),NMat(ndim,ndim)   !Input hamiltonian and norm matrices
		REAL :: tol
	
	INTEGER (KIND = 8), INTENT(OUT) :: nf
	REAL (KIND = 8), INTENT(OUT) :: pevals(nparity,ndim)
	real (kind=8) :: normvals(nparity,ndim)

	integer (kind=4) :: nftmp	
		
	
	LOGICAL :: invFlag = .true.
	integer :: i,it,k,jtt
	integer :: info,lwork,pairP,ii,jj
	integer loop1

	real (kind = 8), ALLOCATABLE :: e(:),normies(:)
	complex*16, allocatable :: vec(:,:)

	allocate(e(ndim),normies(ndim),vec(ndim,ndim))
    CALL zhegvdiag('n',ndim,ndim,Hmat,Nmat,tol,e,normies,vec,nftmp,info)

	if(info/=0)then
		print*,' Problem running zhegvdiag '
		return
	end if
	nf = nftmp
	pevals(ipar,:)=e(:)
	normvals(ipar,:)=normies(:)
!chfjiao
	deallocate(e,normies,vec)
!chfjiao
	return
	
end	subroutine geneigsolverdiag
!=========================================================
!
!  zhegvdiag
!
!  routine for solving generalized eigenvalue problem A vec = aa B vec
!  by diagonalizing B
!
!  calls lapack subroutine ZHEEV
!
! INPUT: jobz = 'N' or 'n' for eigenvalues only, 'V' or 'v' for eigenvectors
!        uplo = 'U' or 'L' for stored in upper or lower triangular
!        np   = declared dimension of the arrays
!        ndim = dimension used by the arrays
!        A(:,:),B(:,:) Hermitian, double complex input matrices
!        tol : real valued tolerance for discarding eigenvalues of B
! OUTPUT:
!      aa(:) : final eigenvalues from A vec = aa B vec
!      bb(:) : eigenvalues of matrx B; used as check
!      vec(:,:) : if jobz = 'v', filled afterwards with eigenvectors
!      nf : # of final eigenvalues kept ( = # of eigenvalues of B > tol )
!      info: if =0, then results okay
!

subroutine zhegvdiag(jobz,np,ndim,A,B,tol,aa,bb,vec,nf,info)
    use dens1body
	implicit none
	interface
	   SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,  INFO ) ! INTERFACE
	        CHARACTER          JOBZ, UPLO
	        INTEGER            INFO, LDA, LWORK, N
	        DOUBLE PRECISION   RWORK( * ), W( * )
	        COMPLEX*16         A( LDA, * ), WORK( * )
		end subroutine zheev      ! INTERFACE
    end interface	
!------------ INPUTS -----------------	
    CHARACTER          JOBZ
    INTEGER            NP,NDIM
	complex*16    ::  a(np,ndim),b(np,ndim)
	real :: tol
!------------ OUTPUTS ------------------
	
	real(kind=8)   :: aa(ndim),bb(ndim)
	complex*16     :: vec(np,ndim)
	integer :: nf
	integer :: info

!-------- INTERMEDIATE----------------
    complex*16, allocatable :: c(:,:)
    real (kind = 8),ALLOCATABLE :: rwork(:)
    COMPLEX (kind = 8), ALLOCATABLE :: work(:),tmpinvrt(:,:)
	integer :: lwork
	integer :: i,j,k,l,ii,jj
	complex*16 :: ztmp
	real(kind=8) :: facti,factj
	
!.......... FIRST STEP..... DIAGONALIZE B..................
!chfjiao
    aa = 0.0d0
    bb = 0.0d0
!chfjiao
    allocate(c(ndim,ndim))
	do i = 1,ndim
		do j = 1,ndim
			c(i,j)=b(i,j)
		end do
	end do
    lwork = 2*ndim - 1
    ALLOCATE(work(lwork))
    ALLOCATE(rwork(3*ndim - 2))
	CALL zheev('v','u',ndim,c,ndim,bb,work,lwork,rwork,info)
	if(info/=0)then
		deallocate(c,work,rwork)
		return
	end if
	nf = 0
	do i = 1,ndim
		if(bb(i) > tol)then
			nf = nf + 1
		end if
	end do
	if(finddense)then
		allocate(tmpinvrt(ndim,nf))
		solutions(myparity,myarrayindex)%nkeep=nf
		allocate(solutions(myparity,myarrayindex)%eval(nf))
		allocate(solutions(myparity,myarrayindex)%gvec(nf,nf))
		allocate(solutions(myparity,myarrayindex)%ninvrt(ndim,nf))

	end if
	
!	print*,nf,' states '
	if(nf ==0)then
		deallocate(c,work,rwork)
		return
	end if
!................. NOW CONSTRUCT NEW MATRIX..... V = 
    ii = 0
    do i = 1,ndim
		if(bb(i)> tol)then
			ii = ii +1
			facti = 1.d0/dsqrt(bb(i))
		else
			cycle
		end if
		
		if(finddense)then
			do j = 1,ndim
				tmpinvrt(j,ii)=c(j,i)*facti
			end do
		end if
		
		jj = 0
		do j = 1 ,ndim
			if(bb(j)> tol)then
				jj = jj +1
				factj = 1.d0/dsqrt(bb(j))
			else
				cycle
			end if
			ztmp = (0.d0,0.d0)
			do k = 1,ndim
				do l = 1,ndim
  				    ztmp = ztmp + dconjg(c(k,i))*a(k,l)*c(l,j)
				end do
			end do
			vec(ii,jj)= ztmp*facti*factj		
		end do		
	end do
	
	if(nf< 9)then
		do ii = 1, nf
			write(66,'(10f9.3)')(vec(ii,jj),jj =1,nf )
		end do		
	end if
		

!... NOW DIAGONALIZE....
	CALL zheev('v','u',nf,vec,np,aa,work,lwork,rwork,info)
!.... TRANSFORM EIGENVECTORS IF NEEDED...TO BE ADDED LATER...	

    if(finddense)then  ! store
			do i = 1,nf
				solutions(myparity,myarrayindex)%eval(i)=aa(i)
				do j = 1,nf
					solutions(myparity,myarrayindex)%gvec(j,i)=vec(j,i)
				end do
				do j = 1,ndim
					solutions(myparity,myarrayindex)%ninvrt(j,i)=tmpinvrt(j,i)
				end do
			end do

		
	end if
	
	deallocate(c,work,rwork)
	if(finddense)deallocate(tmpinvrt)
	
	return
	
end subroutine zhegvdiag
!!=========================================================
!==============================================================================
!
! original generalized eigensolution routines
!  adding expectation values
! using diagonalization
!
!  INPUT:
!    ndim = dimension of matrices
!    nparity = # of parity states considered
!    ipar = which parity (1 or 2)
!   Nmat(i,j) = ndim x ndim norm matrix
!   Hmat(i,j)= ndim x ndim hamiltonian matrix
!   tol = real variable, tolerance for eliminating zeros
!
! OUTPUT
!   nf = # of states found
!   pevals(i=1,nf) = energy eigenvalues found
!   normvals(i=1,ndim)
!
! subroutines called:
!
!  zheev  -- lapack routine for diagonalizing a hermitian, complex*16 matrix
!
subroutine geneigsolverdiag_wexpect(ndim,nparity,ipar,NMat,HMat,ObsMat,tol,nf,pevals,normvals,obsvals)
!	use errortests
	implicit none
	integer (kind=4) :: ndim, nparity   ! dimension of matrices

    integer ipar
	COMPLEX (KIND = 8):: Hmat(ndim,ndim),NMat(ndim,ndim),Obsmat(ndim,ndim)   !Input hamiltonian and norm matrices
		REAL :: tol
	
	INTEGER (KIND = 8), INTENT(OUT) :: nf
	REAL (KIND = 8), INTENT(OUT) :: pevals(nparity,ndim),obsvals(nparity,ndim)
	real (kind=8) :: normvals(nparity,ndim)

	integer (kind=4) :: nftmp	
		
	
	LOGICAL :: invFlag = .true.
	integer :: i,it,k,jtt
	integer :: info,lwork,pairP,ii,jj
	integer loop1

	real (kind = 8), ALLOCATABLE :: e(:),normies(:),obsies(:)
	complex*16, allocatable :: vec(:,:)

	allocate(e(ndim),normies(ndim),vec(ndim,ndim),obsies(ndim))
    CALL zhegvdiag_wexpect('n',ndim,ndim,Hmat,Nmat,Obsmat,tol,e,normies,obsies,vec,nftmp,info)

	if(info/=0)then
		print*,' Problem running zhegvdiag '
		return
	end if
	nf = nftmp
	pevals(ipar,:)=e(:)
	normvals(ipar,:)=normies(:)
	obsvals(ipar,:)=obsies(:)
!chfjiao
	deallocate(e,normies,vec,obsies)
!chfjiao
	return
	
end	subroutine geneigsolverdiag_wexpect
!=========================================================
!
!  zhegvdiag
!
!  routine for solving generalized eigenvalue problem A vec = aa B vec
!  by diagonalizing B
!
!  calls lapack subroutine ZHEEV
!
! INPUT: jobz = 'N' or 'n' for eigenvalues only, 'V' or 'v' for eigenvectors
!        uplo = 'U' or 'L' for stored in upper or lower triangular
!        np   = declared dimension of the arrays
!        ndim = dimension used by the arrays
!        A(:,:),B(:,:) Hermitian, double complex input matrices
!        tol : real valued tolerance for discarding eigenvalues of B
! OUTPUT:
!      aa(:) : final eigenvalues from A vec = aa B vec
!      bb(:) : eigenvalues of matrx B; used as check
!      vec(:,:) : if jobz = 'v', filled afterwards with eigenvectors
!      nf : # of final eigenvalues kept ( = # of eigenvalues of B > tol )
!      info: if =0, then results okay
!

subroutine zhegvdiag_wexpect(jobz,np,ndim,A,B,D,tol,aa,bb,dd,vec,nf,info)

	implicit none
	interface
	   SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,  INFO ) ! INTERFACE
	        CHARACTER          JOBZ, UPLO
	        INTEGER            INFO, LDA, LWORK, N
	        DOUBLE PRECISION   RWORK( * ), W( * )
	        COMPLEX*16         A( LDA, * ), WORK( * )
		end subroutine zheev      ! INTERFACE
    end interface	
!------------ INPUTS -----------------	
    CHARACTER          JOBZ
    INTEGER            NP,NDIM
	complex*16    ::  a(np,ndim),b(np,ndim),D(np,ndim)
	real :: tol
!------------ OUTPUTS ------------------
	
	real(kind=8)   :: aa(ndim),bb(ndim),dd(ndim)
	complex*16     :: vec(np,ndim),obs(np,ndim)
	integer :: nf
	integer :: info

!-------- INTERMEDIATE----------------
    complex*16, allocatable :: c(:,:)
    real (kind = 8),ALLOCATABLE :: rwork(:)
    COMPLEX (kind = 8), ALLOCATABLE :: work(:)
	integer :: lwork
	integer :: i,j,k,l,ii,jj
	complex*16 :: ztmp
	real(kind=8) :: facti,factj
	
!.......... FIRST STEP..... DIAGONALIZE B..................
!chfjiao
    aa = 0.0d0
    bb = 0.0d0
	obs = (0.d0,0.d0)
!chfjiao
    allocate(c(ndim,ndim))
	do i = 1,ndim
		do j = 1,ndim
			c(i,j)=b(i,j)
		end do
	end do
    lwork = 2*ndim - 1
    ALLOCATE(work(lwork))
    ALLOCATE(rwork(3*ndim - 2))
	CALL zheev('v','u',ndim,c,ndim,bb,work,lwork,rwork,info)
	if(info/=0)then
		deallocate(c,work,rwork)
		return
	end if
	nf = 0
	do i = 1,ndim
		if(bb(i) > tol)then
			nf = nf + 1
		end if
	end do
!	print*,nf,' states '
	if(nf ==0)then
		deallocate(c,work,rwork)
		return
	end if
!................. NOW CONSTRUCT NEW MATRIX..... V = 
    ii = 0
    do i = 1,ndim
		if(bb(i)> tol)then
			ii = ii +1
			facti = 1.d0/dsqrt(bb(i))
		else
			cycle
		end if
		jj = 0
		do j = 1 ,ndim
			if(bb(j)> tol)then
				jj = jj +1
				factj = 1.d0/dsqrt(bb(j))
			else
				cycle
			end if
!.............. TRANSFORM HAM.......................			
			ztmp = (0.d0,0.d0)
			do k = 1,ndim
				do l = 1,ndim
  				    ztmp = ztmp + dconjg(c(k,i))*a(k,l)*c(l,j)
				end do
			end do
			vec(ii,jj)= ztmp*facti*factj
!................ TRANSFORM OBS.............			
			ztmp = (0.d0,0.d0)
			do k = 1,ndim
				do l = 1,ndim
  				    ztmp = ztmp + dconjg(c(k,i))*d(k,l)*c(l,j)
				end do
			end do
			obs(ii,jj)= ztmp*facti*factj
		end do		
	end do
	
	if(nf< 9)then
		do ii = 1, nf
			write(66,'(10f9.3)')(vec(ii,jj),jj =1,nf )
		end do		
	end if
		

!... NOW DIAGONALIZE....
	CALL zheev('v','u',nf,vec,np,aa,work,lwork,rwork,info)

!....... COMPUTE OBSERVABLE EXPECTATION VALUES........
    do ii = 1,ndim
		ztmp = (0.d0,0.d0)
		ztmp = (0.d0,0.d0)
		do k = 1,ndim
			do l = 1,ndim
				ztmp = ztmp +dconjg(vec(k,ii))*obs(k,l)*vec(l,ii)
			end do
		end do		
		
		dd(ii)=real(ztmp)
	end do  ! ii
	deallocate(c,work,rwork)
	
	return
	
end subroutine zhegvdiag_wexpect





end module eigenpackage
