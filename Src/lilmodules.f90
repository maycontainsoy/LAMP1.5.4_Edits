!
!   PROGRAM lilLAMP
!   reads in matrices to compute various expectation values
!
!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      module little
		  implicit none

		  integer ::  nsps
		  integer :: numsd
		  integer :: numneut,numprot
		  
!		  integer :: numOfJ
		  logical :: parity_test
      type jaggedArray
          complex(kind=8), allocatable, dimension(:,:) :: MK
		  integer(kind=4) :: dimMK
      end type jaggedArray

	  
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
integer :: npair  ! NEED TO SET
		  
	      logical :: ParityTest  ! if two parities exist
		  logical :: doHam       ! ifcomputing the Hamiltonian
	      real (kind=4) :: Jmax
	      integer (kind=4) :: bigJmax,Jmax2 !2*Jmax+1
	      logical :: isOdd
	      integer (kind=4) :: numOfJ  !default is int(Jmax+1.)
		  integer (kind=4) :: numOfJObs  ! for scalar observables
		  integer (kind=4) :: numOfM
		  
!.............. IMPORTANT ADDITION TO 1.3.0: CHOOSING THE POSSIBLE Jvalues
          real, allocatable :: xJlist(:)  !list of chosen 2xJ values
		  integer, allocatable :: mapJ2indx(:)   ! maps 2xJ to the index
		  integer, allocatable :: startJ2list(:)   ! list of starting indices

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
!
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
		  logical, allocatable :: chooseMval(:) ! for M-values
	  	  integer :: nproblems,npossible
		  real(4) :: eigcrit      ! criterion number for improving eigenvalues
		  
	      logical :: diagp
	      character(LEN = 100) :: path
		  
		  
	  contains 
		  
		  !==============================================================================

		  !
		  ! user enters values for the Norm matrix tolerance and the number of points for inversion.
		  !
		  ! OUTPUT:
		  !   tolerance = cuttoff point for eliminating zeroes in the Norm matrix (during zhegvdiag)
		  !
		  subroutine inputTolerance(tolerance)

		      implicit none
		  !............OUTPUT......................
		      real (kind=4), intent(out) :: tolerance

		  !---------------------- USER ENTERS VALUES ------------------------
		      print *, ' Enter tolerance for norm (typical = 0.0001 to 0.001) '
		      read(*,*) tolerance

		  return; end subroutine inputTolerance
		  
		  
		  
	  end module little
      module sporbit
!C
!C  single-particle ORBIT information
!C
      implicit none

      integer numorb          ! # of s.p. orbits
      integer numruns

!------------ CREATE A DEFINED TYPE----------------
      type orb
        integer :: nr       ! radial quantum number
        integer :: j        ! 2 x j
        integer :: l        ! L
        integer :: w        ! excitation 
        integer :: spstart  ! where these orbits correspond to 
                            ! start in spqn
		integer :: par

      end type orb
	  logical :: allsameparity

      type (orb),allocatable,target :: orbqn(:)

      end module sporbit
