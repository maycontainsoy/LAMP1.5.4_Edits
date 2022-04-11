module lileigenpackage

	implicit none

!..... ADDED IN 1.3.6...
    integer :: nlevels
	real(kind=8), allocatable :: energylevels(:),jlevels(:)
!.... ADDED IN 1.4.4.....	
	real(kind=8), allocatable :: obslevels(:)
	character, allocatable :: paritylevels(:)
	
	logical :: stdeigen
	
	real :: jtarget
	
	
!	
! NEED TO STORE INFORMATION ON SOLUTIONS

   type solnstore
	   integer :: localdim ! # bare dimension 
	   integer :: mkdim
	   integer :: nkeep  ! # of non-null norm eigenvalues
	   complex(kind=8), allocatable :: ninvrt(:,:)
	   complex(kind=8), allocatable :: gvec(:,:)
	   real(kind=8), allocatable :: eval(:)
	   
   end type solnstore		
   
   integer :: myarrayindex  ! used to communicate which set
   integer :: myparity,myJ
   
   type (solnstore), allocatable,target :: solutions(:,:)

   logical :: savesolutions
   
   integer, allocatable :: mapstates(:,:)
   	   
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

    use little

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
	if(paritytest)nlevels = 2*nlevels
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

    if(savesolutions)allocate(solutions(numOfJ,npair))

    do intJ = 1, numOfJ
		
        floatJ = xJlist(intJ)
		matdim = N_J_MK(1,1,intJ)%dimMK
		localdim = matdim*numsdused
        allocate(ParityEvals(npair,localdim))
        allocate(normvals(npair,localdim))
        allocate(TSDHmat(localdim,localdim))
        allocate(TSDNmat(localdim,localdim))
		if(compute_expect)allocate(obsSDHmat(localdim,localdim),obsvals(npair,localdim))

        nkeep = 0
        ParityEvals = 0.0d0
		
!		print*,intJ,localdim

        do parityOUT = 1, npair
            TSDHmat = cmplx(0.0d0,0.0d0)
            TSDNmat = cmplx(0.0d0,0.0d0)
			if(savesolutions)then
				solutions(intJ,parityOUT)%localdim=localdim
				allocate(solutions(intJ,parityOUT)%eval(localdim))
				allocate(solutions(intJ,parityOUT)%gvec(localdim,localdim))
				allocate(solutions(intJ,parityOUT)%ninvrt(localdim,localdim))
				myJ = intJ
				myparity=parityOUT
				solutions(intJ,parityOUT)%mkdim = matdim
			
			end if
		
			do sdX = 1,numsdused
				do ii = 1,matdim
					a = ii+ (sdX-1)*matdim  ! this is the mapping
					do sdY = sdX,numsdused
						do jj = 1,matdim
							b = jj + (sdY-1)*matdim
				            if (ParityTest) then
             TSDHmat(a,b) = 0.5d0*(H_J_MK(sdX,sdY,intJ)%MK(ii,jj) + (-1.0d0)**(parityOUT+1)*PH_J_MK(sdX,sdY,intJ)%MK(ii,jj))
             TSDNmat(a,b) = 0.5d0*(N_J_MK(sdX,sdY,intJ)%MK(ii,jj) + (-1.0d0)**(parityOUT+1)*PN_J_MK(sdX,sdY,intJ)%MK(ii,jj))
			 if(compute_expect)then
	             obsSDHmat(a,b)=0.5d0*(Hobs_J_MK(sdX,sdY,intJ)%MK(ii,jj)+(-1.0d0)**(parityOUT+1)*PHobs_J_MK(sdX,sdY,intJ)%MK(ii,jj))
			 end if
			 
							else
				                TSDHmat(a,b) = H_J_MK(sdX,sdY,intJ)%MK(ii,jj)
				                TSDNmat(a,b) = N_J_MK(sdX,sdY,intJ)%MK(ii,jj) 	
								if(compute_expect)obsSDHmat(a,b)=Hobs_J_MK(sdX,sdY,intJ)%MK(ii,jj)					
							end if
!------------ ENFORCE HERMITICIY ------------------------------------------							
							if(sdY > sdX)then
								TSDNmat(b,a)=dconjg(TSDNmat(a,b))
								TSDHmat(b,a)=dconjg(TSDHmat(a,b))
								if(compute_expect)obsSDHmat(b,a)= dconjg(obsSDHmat(a,b))
								
							end if
							
						end do !jj
					end do ! sdY
				end do ! ii
			end do !sdX
			
			

            do kk = 1, localdim
                problist(parityOUT,intJ) = problist(parityOUT,intJ) + dble(TSDNmat(kk,kk))
                hamlist(parityOUT,intJ) = hamlist(parityOUT,intJ) + dble(TSDHmat(kk,kk))
            end do
			
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
			end if


!-----------------
            if(compute_expect)then
               call geneigsolverdiag_wexpect(localdim,npair,parityOUT,TSDNmat,TSDHmat,ObsSDHmat, & 
			       tolerance,nf,ParityEvals,normvals,obsvals)
		     else
               call geneigsolverdiag(localdim,npair,parityOUT,TSDNmat,TSDHmat,tolerance,nf,ParityEvals,normvals)
   		    end if
			if(savesolutions)solutions(intJ,parityOUT)%nkeep = nf ! this is wrong
			
            if (nf > nkeep) nkeep = nf
        end do

        normsum = normsum + problist(1,intJ)
        hamsum = hamsum + hamlist(1,intJ)
        if (parityTest) then
            normsum = normsum + problist(2,intJ)
            hamsum = hamsum + hamlist(2,intJ)
        end if
		

        do stateFound = 1, int(nkeep)
            nftotal = nftotal + 1
            jfound(nftotal) = floatJ!!
            Efound(1,nftotal) = ParityEvals(1,stateFound)!!
			if(compute_expect)obsfound(1,nftotal)=obsvals(1,stateFound)
			if(compute_expect .and. ParityTest)obsfound(2,nftotal)=obsvals(2,stateFound)
			
            if (ParityTest) Efound(2,nftotal) = ParityEvals(2,stateFound)!!
			
			if(ParityTest)then
				if(Efound(1,nftotal)/=0.0)then
					nlevels = nlevels + 1
					energylevels(nlevels)=Efound(1,nftotal)
					jlevels(nlevels)=jfound(nftotal)
					paritylevels(nlevels)='+'
					obslevels(nlevels)=obsfound(1,nftotal)

					
				end if
				if(Efound(2,nftotal)/=0.0)then
					nlevels = nlevels + 1
					energylevels(nlevels)=Efound(2,nftotal)
					jlevels(nlevels)=jfound(nftotal)
					paritylevels(nlevels)='-'
					obslevels(nlevels)=obsfound(2,nftotal)
					
				end if
			else
				if(Efound(1,nftotal)/= 0.0)then
					nlevels = nlevels + 1
					energylevels(nlevels)=Efound(1,nftotal)
					jlevels(nlevels)=jfound(nftotal)
					paritylevels(nlevels)='+'
					obslevels(nlevels)=obsfound(1,nftotal)
					
				end if
			end if
        end do
        if (nkeep<1) write(*,*)' No states for J = ',floatJ

        deallocate(ParityEvals,normvals,TSDHmat,TSDNmat)
		if(compute_expect)deallocate(obsSDHmat,obsvals)
    end do
    if(savesolutions)allocate(mapstates(nlevels,3))
	

    1001 format(I3,2X,2F12.5,2X,F4.1)
    1002 format(I3,2X,F12.5,2X,F4.1)
	
return; end subroutine EigenSolverPackage

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
	use little
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
	if(stdeigen)then
        CALL zhegvdiag('n',ndim,ndim,Hmat,Nmat,tol,e,normies,vec,nftmp,info)
	else
	    CALL zhegvdiag_alt('n',ndim,ndim,Hmat,Nmat,tol,e,normies,vec,nftmp,info)
	end if
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
    COMPLEX (kind = 8), ALLOCATABLE :: work(:)
	integer :: lwork
	integer :: i,j,k,l,ii,jj
	complex*16 :: ztmp
	real(kind=8) :: facti,factj
	
	
!.......... FIRST STEP..... DIAGONALIZE B (norm matrix)..................
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
!	print*,nf,' states '
	if(nf ==0)then
		deallocate(c,work,rwork)
		return
	end if
!................. NOW CONSTRUCT NEW MATRIX..... V = 
    ii = 0
	
	if(savesolutions)then
		solutions(myJ,myparity)%ninvrt=(0.d0,0.d0)
		solutions(myJ,myparity)%gvec=(0.d0,0.d0)
	end if
	
    do i = 1,ndim
		if(bb(i)> tol)then
			ii = ii +1
			facti = 1.d0/dsqrt(bb(i))
		else
			cycle
		end if
		
		if(savesolutions)then
			do k = 1,ndim
				solutions(myJ,myparity)%ninvrt(k,ii)=c(k,i)*facti
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
			
		end do		! j
	end do !i
	
	if(nf< 9)then
		do ii = 1, nf
			write(66,'(10f9.3)')(vec(ii,jj),jj =1,nf )
		end do		
	end if
		

!... NOW DIAGONALIZE....transformed Hamiltonian
	CALL zheev('v','u',nf,vec,np,aa,work,lwork,rwork,info)
!.... TRANSFORM EIGENVECTORS IF NEEDED...TO BE ADDED LATER...	

if(savesolutions)then
	solutions(myJ,myparity)%eval(1:nf)=aa(1:nf)
	do i = 1,nf
		solutions(myJ,myparity)%gvec(1:nf,i)=vec(1:nf,i)
	end do 
	
end if
	
	deallocate(c,work,rwork)
	
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
	use little
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
!!=========================================================
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

subroutine zhegvdiag_alt(jobz,np,ndim,A,B,tol,aa,bb,vec,nf,info)

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
	
	real(kind=8)   :: aa(ndim),bb(ndim),ee(ndim)
	complex*16     :: vec(np,ndim)
	integer :: nf
	integer :: info

!-------- INTERMEDIATE----------------
    complex*16, allocatable :: c(:,:)
    real (kind = 8),ALLOCATABLE :: rwork(:)
    COMPLEX (kind = 8), ALLOCATABLE :: work(:)
	integer :: lwork
	integer :: i,j,k,l,ii,jj
	complex*16 :: ztmp
	real(kind=8) :: facti,factj,dtmp
	
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
		print*,' some problem with norm '
		return
	end if
	
!......... NEXT, DIAGONALIZE A ...................	
	
    do i = 1,ndim
	   do j = 1,ndim
		  vec(i,j)=a(i,j)
	   end do
    end do
	CALL zheev('v','u',ndim,vec,ndim,aa,work,lwork,rwork,info)
	
!......... TAKE EXPECTATION VALUES OF B..........

    nf = 0
    do i = 1,ndim	
		dtmp = 0.d0
		do j = 1,ndim
			do k = 1,ndim
				dtmp = dtmp+vec(j,i)*b(j,k)*vec(k,i)
			end do
		end do
			
			if(dtmp > tol)then
				
				aa(i)=aa(i)/dtmp
				nf = nf+1
			else
				aa(i)=0.00
			end if
			
		
	end do

	
	deallocate(c,work,rwork)
	
	return
	
end subroutine zhegvdiag_alt
!!=========================================================
!====================================================================
!
!  revised CWJ June 2015;
!  also prints out probabilities
!  revised again CWJ April 2018:
!  more logical output; sorts on energies
! 
!  INPUT:
!     np      :  master declared dimension of arrays
!     nstates : = # of states
!     jval(:) : array of jvalues
!     eval(:) : array of energies
!     pval(:) : array of parities (+,-1)
!     obsvals(:) :: array of observables (if compute_expect)
!     parityflag : logical flag if true then possible to have different parities
!     nprint: how many to print to screen; if negative, do not ask to write to file
!

SUBROUTINE J_WriteResults(tol,np,nf,jvals,pevals,obsvals,problist,hamlist,parityflag,nprint)

use little
implicit none

real,    intent(in) :: tol
integer (kind=8), intent(in) :: np, nf
real    (kind=8), intent(in) :: jvals(np), pevals(2,np),obsvals(2,np)
real    (kind=8), intent(in) :: problist(2,numOfJ),hamlist(2,numOfJ)
logical, intent(in) :: parityflag
integer(4),intent(in) :: nprint  ! # to print to screen; if < 0, do not ask to write to file
integer iprint,icount
real(kind=8) :: probsum,hamSum,ex
INTEGER :: i,iostatus
CHARACTER (LEN = 1) :: choice
CHARACTER (LEN = 4) :: shell
CHARACTER (LEN = 26) :: name
CHARACTER (LEN = 100) :: filename!, path



call sortresults


write(6,*)' '
iprint = 0
if(parityflag)then
    if(compute_expect)then
		
	write(6,*)' State    E     Ex      J     parity      <obs>'
else
	write(6,*)' State    E     Ex      J     parity'
end if
	write(6,*)' ----------------------------------------------'
		
	
	do i = 1,min(nlevels,abs(nprint))
		ex= energylevels(i)-energylevels(1)
		if(jtarget < 0 .or. jlevels(i)==jtarget)then
			if(compute_expect)then
		write(6,3001)i,energylevels(i),ex,jlevels(i),paritylevels(i),obslevels(i)
	else
		
		write(6,3001)i,energylevels(i),ex,jlevels(i),paritylevels(i)
		
	endif
		
		iprint= iprint+1
		if(iprint ==abs(nprint))exit
     	end if
		
	end do
else
	if(compute_expect)then
	write(6,*)' State      E           Ex        J         <obs>  '
else
	write(6,*)' State      E           Ex        J     '
end if
	write(6,*)' ----------------------------------------'
do i = 1,nlevels !min(nlevels,nprint)
	ex= energylevels(i)-energylevels(1)
	if(jtarget < 0 .or. jlevels(i)==jtarget)then
		if(compute_expect)then
			write(6,3002)i,energylevels(i),ex,jlevels(i),obslevels(i)
		else
	        write(6,3002)i,energylevels(i),ex,jlevels(i)
        end if
	    iprint= iprint+1
	    if(iprint ==nprint)exit
    end if
	
end do

end if	


1000 FORMAT(I3,4(G15.8))

if(nprint < 0)return

WRITE(*,*) 'Write spectral output to .res file? (Y or N)'

DO
	READ(*,*) choice
	IF ((choice == 'y').OR.(choice == 'Y')) THEN
		EXIT
	ELSEIF ((choice == 'n').OR.(choice == 'N')) THEN
		RETURN
	ELSE
		WRITE(*,*) 'Y or N please.'
	ENDIF
END DO

DO 
	WRITE(*,*) 'Enter file name (without extention -- .res added): '
	READ(*,*) name
	filename = TRIM(name)//'.res'
	OPEN(UNIT = 10, FILE = filename, STATUS = 'NEW',IOSTAT = iostatus)
	IF (iostatus > 0) THEN
		WRITE(*,*) "File already exists: Overwrite (o), append (a), new file name (n)?"
		READ(*,*) choice
		DO icount = 1,5
			IF ((choice == 'o').OR.(choice == 'O')) THEN
				CLOSE(UNIT = 10)
				OPEN(UNIT = 10, FILE = filename, STATUS = 'REPLACE')
				EXIT
			ELSEIF ((choice == 'a').OR.(choice == 'A')) THEN
				CLOSE(UNIT = 10)
				OPEN(UNIT = 10,FILE = filename, STATUS = 'OLD', POSITION = 'APPEND')
				EXIT
			ELSEIF ((choice == 'n').OR.(choice == 'N')) THEN
				CLOSE(UNIT = 10)
				EXIT
			ELSE
				WRITE(*,*) 'Incorrect choice.  Please select again.'
			END IF
		END DO
		IF ((choice.NE.'n').AND.(choice.NE.'N')) EXIT
	ELSE
		EXIT
	END IF
END DO

call sortresults

write(10,*)' LAMP Results '
write(10,*)' (cutoff criterion = ',tol,')'
write(10,*)' '
if(parityflag)then
	if(compute_expect)then
	   write(10,*)' State    E     Ex      J     parity    < obs >'
   else
	   write(10,*)' State    E     Ex      J     parity'
    end if
	write(10,*)' ------------------------------------------------'
!	DO i = 1, nf
!		WRITE(10,1001) i,pevals(1,i),pevals(2,i),jvals(i)		
!	END DO
	1001 format(i3,2x,2F13.5,2x,f4.1)		
	
	do i = 1,nlevels
		ex= energylevels(i)-energylevels(1)
		if(compute_expect)then
			write(10,3001)i,energylevels(i),ex,jlevels(i),paritylevels(i),obslevels(i)
		else
			write(10,3001)i,energylevels(i),ex,jlevels(i),paritylevels(i)			
		end if
		
	end do
	3001 format(i3,2x,2f13.5,2x,f4.1,2x,a1,2x,f13.5)

else
	if(compute_expect)then
	write(10,*)' State       E        Ex         J       <obs >  '
    else
	write(10,*)' State       E        Ex         J     '
    end if 

	write(10,*)' ------------------------------------'
!	DO i = 1, nf
!		WRITE(10,1002) i,pevals(1,i),jvals(i)		
!	END DO
1002 format(i3,2x,F12.5,2x,f4.1)	
do i = 1,nlevels
	ex= energylevels(i)-energylevels(1)
	if(compute_expect)then
	write(10,3002)i,energylevels(i),ex,jlevels(i),obslevels(i)
    else
	write(10,3002)i,energylevels(i),ex,jlevels(i)
     end if
	
end do
3002 format(i3,2x,2f13.5,2x,f4.1,2x,f13.5)

end if	
!DO i = 1, nf
!	WRITE(10,1000) i,pevals(1,i),pevals(2,i),jvals(i)
!END DO
probsum = 0.d0 
write(10,*)' '
write(10,*)' Fraction in original HF state: Norm'
write(10,*)' '
if(parityflag)then
	write(10,*)' J    frac(+)   frac(-)'
	write(10,*)'-----------------------'
	do i = 1,numOfJ
		write(10,2001)xjlist(i),problist(1,i),problist(2,i)
		probsum = probsum+problist(1,i)+problist(2,i)
	end do
	2001 format(f4.1,2f10.6)
else
	write(10,*)' J    frac'
	write(10,*)'-----------------------'
	do i = 1,numOfJ
		write(10,2001)xjlist(i),problist(1,i)
		probsum = probsum+problist(1,i)
	end do
end if
write(10,*)' Total of HF state = ',probsum

hamsum = 0.d0
write(10,*)' '
write(10,*)' Fraction in original HF state: Hamiltonian'
write(10,*)' '
if(parityflag)then
	write(10,*)' J    frac(+)   frac(-)'
	write(10,*)'-----------------------'
	do i = 1,numOfJ
		write(10,2002)xjlist(i),hamlist(1,i),hamlist(2,i)
		hamsum = hamsum+hamlist(1,i)+hamlist(2,i)
	end do
	2002 format(f4.1,2(1X,f10.6))
else
	write(10,*)' J    frac'
	write(10,*)'-----------------------'
	do i = 1,numOfJ
		write(10,2002)xjlist(i),hamlist(1,i)
		hamsum = hamsum+hamlist(1,i)
	end do
end if
write(10,*)' Total of sum of trace(H) = ',hamSum

CLOSE(UNIT = 10)

WRITE(*,*) ''
WRITE(*,*) 'Data written to:',filename
return

END SUBROUTINE J_WriteResults
!===================================================================
!---- ADDED IN 1.3.6-----------
!
subroutine sortresults
	implicit none
	integer :: i,j,k
	character ::parswap
	real(8) :: eswap,jswap,etest,oswap
	
	do i = 1,nlevels -1
		etest = energylevels(i)

		k = i
		do j = i+1,nlevels
			if(energylevels(j)< etest)then
				k = j
				etest = energylevels(j)
			end if
		
		end do
		if(k > i)then ! swap
			eswap = energylevels(i)
			jswap = jlevels(i)
			parswap=paritylevels(i)
			oswap = obslevels(i)
			energylevels(i)=energylevels(k)
			jlevels(i) = jlevels(k)
			paritylevels(i)=paritylevels(k)
			obslevels(i)=obslevels(k)
			energylevels(k)=eswap
			jlevels(k)     =jswap
			paritylevels(k)=parswap
			obslevels(k)=oswap
			
			
		end if
	
	end do
	
	return
	
end subroutine sortresults
!===================================================================

subroutine writeoutlevels(outfile,nlevels)
	use sporbit,only:allsameparity
	implicit none
	
	integer :: outfile
	integer :: nlevels
	integer :: ilevel
	real :: ex
	if(allsameparity)then
	   write(outfile,'("  State     E            Ex       J     T (dummy)")')
    else 
 	   write(outfile,'("  State     E            Ex       J     T (dummy)  parity")')
		
	end if
	   
	
	do ilevel = 1,nlevels
		ex= energylevels(ilevel)-energylevels(1)
		if(allsameparity)then

			write(outfile,3002)ilevel,energylevels(ilevel),ex,jlevels(ilevel),0.000
		else
			write(outfile,3003)ilevel,energylevels(ilevel),ex,jlevels(ilevel),0.000,paritylevels(ilevel)
			
		end if
	end do ! ilevel
	3002 format(i3,2x,2f13.5,2x,f4.1,2x,f4.1)
	3003 format(i3,2x,2f13.5,2x,f4.1,2x,f4.1,10x,a1)
	
	write(outfile,*)' '
	
	return
end subroutine writeoutlevels


!===================================================================

subroutine writeoutorbits(outfile)
	use sporbit
	implicit none
	integer :: outfile
	integer :: a
	write(outfile,'("  Single particle state quantum numbers ")')
	write(outfile,'("ORBIT     N    L   2 x J ")')
	do a = 1,numorb
		write(outfile,'(4i5)')a,orbqn(a)%nr,orbqn(a)%l,orbqn(a)%j
		
	end do
	write(outfile,*)' '
	return
	
	
end subroutine writeoutorbits


!===================================================================
subroutine mapresults
	use little
	implicit none
	
	integer :: i,j,k,kstate
	real :: jshift
	real(8) :: ee
	integer :: iJ, iPar
	
	if(.not.savesolutions)return
	
	
	jshift = xJlist(1)
	do i = 1,nlevels
		iJ = nint(jlevels(i)-jshift+1)
		if(paritylevels(i)=='-')then
			iPar = 2
		else
			iPar = 1
		end if
		ee = energylevels(i)
		
		kstate = -1
		
!		print*,' testing ',i,jlevels(i),ij,nlevels
		do k = 1,solutions(ij,ipar)%nkeep
			if(abs(solutions(ij,ipar)%eval(k)-ee)< 0.00001)then
				kstate = k
				exit
			end if
			
		end do
		if(kstate < 0)then
			print*,' could not find state ',i,jlevels(i),paritylevels(i),energylevels(i)
			stop
		end if
		mapstates(i,1)=iJ
		mapstates(i,2)=ipar
		mapstates(i,3)=kstate
		
		
	end do
	
	
	return
	
end subroutine mapresults
!===================================================================

subroutine write_density_results(resultfile)
	use sporbit
	use dens1b
	implicit none
	integer :: resultfile
	
	
	integer :: maxstates
	
	integer :: istate,fstate
	integer :: iJ, fJ
	integer :: ipar, fpar,iloc,floc
	
	integer :: nilevels, nflevels
	integer :: ilevel, flevel
	integer :: ikeep,fkeep,nikeep,nfkeep
	
	complex(8) :: prho(numorb,numorb),nrho(numorb,numorb)
	
	complex(8),allocatable :: transi(:),transf(:)
	
	integer :: imatdim,fmatdim
	integer :: sdi,sdf,iM,fM
	
	integer :: Jtmin,Jtmax,Jt
	
	integer :: a,b
	integer :: apar, bpar
	real :: sumprot,sumneut
	
	complex(8) :: trnorm
	
	integer :: maxlevel
	real :: ei,ef
	
	character :: fparchar,iparchar
	
	

	print*,' Enter max # of states to include in output of density matrices '
	print*,' (max of ',nlevels,')'
	read*,maxstates
	if(maxstates > nlevels) maxstates=nlevels
	maxlevel = 0
    do istate = 1,maxstates
		iJ = mapstates(istate,1)
		ipar = mapstates(istate,2)
	 	maxlevel = max(maxlevel, solutions(iJ,ipar)%localdim)	
!	    print*,iJ,ipar,solutions(iJ,ipar)%localdim
	end do
	allocate(transi(maxlevel),transf(maxlevel))
	
	
!.............. WRITE OUT ENERGIES ETC

call writeoutlevels(33,maxstates)
!............ WRITE OUT ORBITS FOR REFERENCR.....
call writeoutorbits(33)

!.............. LOOP OVER DENSITY MATRICES........

    do istate = 1,maxstates
		iJ = mapstates(istate,1)
		ipar = mapstates(istate,2)
		
		if(ipar==1)then
			iparchar='+'
		else
			iparchar='-'
		end if
		iloc = mapstates(istate,3)
		nilevels = solutions(iJ,ipar)%localdim
		nikeep   = solutions(iJ,ipar)%nkeep
		
!............. CONSTRUCT TRANSFORMATION.......
        transi = (0.d0, 0.d0)
		
		do ikeep = 1,nikeep
			do ilevel = 1,nilevels
				transi(ilevel)= transi(ilevel) +  solutions(iJ,ipar)%ninvrt(ilevel,ikeep)*solutions(iJ,ipar)%gvec(ikeep,iloc)
			end do
			
		end do
		
		imatdim = solutions(iJ,ipar)%MKdim
		
		do fstate = 1,maxstates
			fJ = mapstates(fstate,1)
			fpar = mapstates(fstate,2)
			if(fpar==1)then
				fparchar='+'
			else
				fparchar='-'
			end if
			floc = mapstates(fstate,3)
			nflevels = solutions(fJ,fpar)%localdim
			nfkeep   = solutions(fJ,fpar)%nkeep
			
!................... transform.................
!............. CONSTRUCT TRANSFORMATION.......
            transf = (0.d0, 0.d0)
		
		    do fkeep = 1,nfkeep
			    do flevel = 1,nflevels
				    transf(flevel)= transf(flevel) +  solutions(fJ,fpar)%ninvrt(flevel,fkeep)*solutions(fJ,fpar)%gvec(fkeep,floc)
			    end do
			
  		    end do
			fmatdim = solutions(fJ,fpar)%MKdim
!.......... LOOP OVER Jt................
            Jtmax = nint(xJlist(iJ)+xJlist(fJ))
            Jtmin = abs(nint(xJlist(iJ)-xJlist(fJ)))
			ei =energylevels(istate) !solutions(iJ,ipar)%eval(ikeep)
			ef =energylevels(fstate) !solutions(fJ,fpar)%eval(fkeep)
			
            write(resultfile,*)' '
			
			if(allsameparity)then

                write(resultfile,333)istate,ei,int(2*xJlist(iJ)),abs(numneut-numprot)
                write(resultfile,334)fstate,ef,int(2*xJlist(fJ)),abs(numneut-numprot)
				
			else
                write(resultfile,1333)istate,ei,int(2*xJlist(iJ)),abs(numneut-numprot),iparchar
                write(resultfile,1334)fstate,ef,int(2*xJlist(fJ)),abs(numneut-numprot),fparchar
			end if

			do Jt = Jtmin,Jtmax

	            prho =(0.d0, 0.d0)
	            nrho = (0.d0, 0.d0)			
				
				trnorm=(0.d0,0.d0)
			
				do ilevel = 1,nilevels
					!a = ii+ (sdX-1)*matdim  ! this is the mapping
					iM = mod(ilevel-1,imatdim)+1
					sdi = (ilevel-iM)/imatdim+1
									
					do flevel = 1,nflevels
						fM = mod(flevel-1,fmatdim)+1
						sdf = (flevel-fM)/fmatdim+1
						
						
						if(Jt < rhotr(sdf,sdi,fJ,fM,iJ,iM)%Jtmin)cycle
						if(Jt > rhotr(sdf,sdi,fJ,fM,iJ,iM)%Jtmax)cycle
						if( rhotr(sdf,sdi,fJ,fM,iJ,iM)%Jtmin ==-999)cycle
						if( rhotr(sdf,sdi,fJ,fM,iJ,iM)%Jtmax ==-999)cycle
						
						if(allsameparity)then

						if(numprot > 0)prho = prho + dconjg(transf(flevel)) * transi(ilevel) * rhotr(sdf,sdi,fJ,fM,iJ,iM)%rhoJt(Jt)%prho
						if(numneut > 0)nrho = nrho + dconjg(transf(flevel)) * transi(ilevel) * rhotr(sdf,sdi,fJ,fM,iJ,iM)%rhoJt(Jt)%nrho
					    
						else  !this projects out the initial state (righthand side) with correct parity (ipar)
							
							if(numprot > 0)prho = prho + dconjg(transf(flevel)) * transi(ilevel) *  & 
							 0.5* ( rhotr(sdf,sdi,fJ,fM,iJ,iM)%rhoJt(Jt)%prho + (-1)**(ipar+1)*rhotr(sdf,sdi,fJ,fM,iJ,iM)%rhoJt(Jt)%Pprho  )
							  
							if(numneut > 0)nrho = nrho + dconjg(transf(flevel)) * transi(ilevel) *  & 
							0.5*  (rhotr(sdf,sdi,fJ,fM,iJ,iM)%rhoJt(Jt)%nrho + 	(-1)**(ipar+1)*rhotr(sdf,sdi,fJ,fM,iJ,iM)%rhoJt(Jt)%Pnrho ) 						
	
						end if
					
					end do
					
				end do
!................... write out real part of density matrices .................			



                 write(resultfile,32)Jt
32               format(' Jt = ',i3,', proton      neutron ')
				do a = 1,numorb
					apar = orbqn(a)%par
 					do b = 1,numorb
						bpar = orbqn(b)%par
						if(allsameparity)then
						if(abs(real(prho(a,b)))> 0.000001 .or.abs(real(nrho(a,b)))> 0.000001 )then
							write(resultfile,111)a,b,real(prho(a,b)),real(nrho(a,b))
						end if
						
						
					    else  ! check parity; by only writing out the combination of single-particle orbits
							! with the correct change in parity, we guarantee we 'connect' only to 
							! the final state (lefthand side) with the desired parity (fpar)
							if(apar*bpar /= (-1)**(ipar+fpar))cycle  
							if(abs(real(prho(a,b)))> 0.000001 .or.abs(real(nrho(a,b)))> 0.000001 )then
								write(resultfile,111)a,b,real(prho(a,b)),real(nrho(a,b))
							end if							
							
							
						end if
				 
				 
					end do
				end do
				111 format(2i5,2f10.5)
!...............  check sum rules

                  if(jt ==0 .and. istate==fstate)then
					  sumprot = 0.0
					  sumneut = 0.0
					  do a = 1,numorb
						  sumprot = sumprot + real(prho(a,a))*sqrt( orbqn(a)%j + 1.)
						  sumneut = sumneut + real(nrho(a,a))*sqrt( orbqn(a)%j + 1.)
						  
					  end do
					  
					  print*,xJlist(iJ),'testing sum rule ',istate,sumprot/(sqrt(2.*xJlist(iJ)+1.0)),sumneut/(sqrt(2.*xJlist(iJ)+1.0))
					  
				  end if			
		    end do ! jt

		
		end do ! fstate
		

		
	end do ! istate
	
	
	433     format(' Initial state #',i5,' E = ',f10.5,' 2xJ   = ',i4) 
	434     format(' Final state   #',i5,' E = ',f10.5,' 2xJ   = ',i4) 
	333     format(' Initial state #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
	334     format(' Final state   #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
	1333     format(' Initial state #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4,5x,a1) 
	1334     format(' Final state   #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4,5x,a1) 		
    write(resultfile,*)'+++++++++++++++++++++++++++++++++++++++++++++'
    write(resultfile,*)' '
    write(resultfile,*)' Definition of density matrices : '
    write(resultfile,*)' rho_K(a^+b) =   (Jf || (a^+ b)_K || Ji) / sqrt(2K+1) '
    write(resultfile,*)'  where reduced matrix element is convention of Edmonds '
    write(resultfile,*)' (note: if isospin is good symmetry, then '
    write(resultfile,*)'  doubly reduced/ divided by sqrt(2T+1) as well'
    write(resultfile,*)' '
    write(resultfile,*)' Note time-reversal symmetry relation: '
    write(resultfile,*)' rho_K(a^+b, Ji->Jf) = (-1)^(ja-jb + Ji-Jf) rho_K(b^+a, Jf->Ji) '
    write(resultfile,*)' For isospin, add in factor (-1)^(Ti-Tf)'
    write(resultfile,*)' '
    write(resultfile,*)'+++++++++++++++++++++++++++++++++++++++++++++'
    write(resultfile,*)' '
	
	return

end subroutine write_density_results
!===================================================================


end module lileigenpackage
