!
! ADDED IN 1.4.6  Sept 2019 by CWJ @ SDSU
!
! data and routines for extracting one-body densities from LAMP
!
!  KEY INTERMEDIATE STEP
!  original density matrices are in M-scheme rho(i,j) where i = ja, ma 
!  need to convert to a,b, Kab, Mab, which are easier to interpret
!  for each block a x b  there are (2ja+1)*(2jb+1) entries
!  these get converted to Kab = Kmin, Kmax, with 2Kab+1 Mab values each
!
!
!
!  Xrho1b(isd,fsd,iomega)%x1brhoJM(i,j)


module dens1body
	use sporbit
	use spstate
	use lamplight
	use psis
	implicit none
	
!	logical :: get1bodydens = .false.
	
	
	type rhobase
		complex(kind=8), allocatable :: prho(:,:),nrho(:,:)
		complex(kind=8), allocatable :: Pprho(:,:),Pnrho(:,:)   ! parity-inverted on initial
		
	end type rhobase	
	type rhobaseJt
		integer :: Jtmin,Jtmax
		type (rhobase), allocatable :: rhoJt(:)
	end type rhobaseJt
	
	type rhobaseMK
		integer :: dimMK
		type (rhobase), allocatable :: MK(:,:)
	end type rhobaseMK
	
	type rhobaseJiJf
		type(rhobase), allocatable :: mat(:,:)
		
	end type rhobaseJiJf
	
	type (rhobase), allocatable, target :: Xrho1b(:,:,:,:,:)   ! isd,fsd, alpha, beta,gamma
		
	type(rhobase), allocatable :: rho1b_jMK(:,:,:,:,:)
! rho_1b_J_MK = < alpha | [ a^+ x b ]_Jt Mt P^Jbeta_MK | beta >		
	type(rhobase),allocatable :: rhotild_Jprime_MK(:,:,:,:,:)
	type (rhobaseMK), allocatable :: rho1b_J_MK(:,:,:)	
		

!  rhotr(sdf,sdi,Jf,Ji,Kf,Ki)	
	type (rhobaseJt), allocatable,target  :: rhotr(:,:,:,:,:,:)	
	type (rhobaseJiJf), allocatable,target ::  rhofinal(:,:,:,:) 	
		
! NEED TO STORE INFORMATION ON SOLUTIONS

   type solnstore
	   integer :: nkeep  ! # of non-null norm eigenvalues
	   complex(kind=8), allocatable :: ninvrt(:,:)
	   complex(kind=8), allocatable :: gvec(:,:)
	   real(kind=8), allocatable :: eval(:)
	   
   end type solnstore		
   
   integer :: myarrayindex  ! used to communicate which set
   integer :: myparity
   
   type (solnstore), allocatable,target :: solutions(:,:)
	   
   logical :: finddense	   
	
contains

!......... SUBROUTINE TO	
!
! CALLED BY: SAMPLE4DENSITY
!
	subroutine setup_den1b_orig_array(nSDs,nEuler)
!		use system_parameters
        use sporbit
		implicit none
		integer, intent(in) :: nSDs  ! # of SDs
		integer, intent(in) :: nEuler ! number of angles
			
		integer :: a, b
		integer :: Jabmin,Jabmax,Jjj,Mmm
		integer :: isd,fsd,iEuler
		integer :: alp,bet,gam
		
		
		allocate(Xrho1b(nsds,nsds,J2max1,numOfJ,J2max1))

		do isd = 1,nsds
			do fsd = 1,nsds
				do alp = 1,J2max1
					do bet = 1,numOfJ
						do gam = 1,J2max1
					if(numprot > 0)allocate( Xrho1b(isd,fsd,alp,bet,gam)%prho(nsps,nsps))
					if(numneut > 0)allocate( Xrho1b(isd,fsd,alp,bet,gam)%nrho(nsps,nsps))
					if(.not.allsameparity)then
						if(numprot > 0)allocate( Xrho1b(isd,fsd,alp,bet,gam)%Pprho(nsps,nsps))
						if(numneut > 0)allocate( Xrho1b(isd,fsd,alp,bet,gam)%Pnrho(nsps,nsps))						

					end if
				       end do 
			        end do
					
				end do
			end do  ! fsd
		end do      ! isd
		
		print*,' ARRAYS SET UP FOR DENSITIES'

		return
	end subroutine setup_den1b_orig_array

!===========================================


!----------------------------------------------------------------------------------------------	
end module dens1body
	
	