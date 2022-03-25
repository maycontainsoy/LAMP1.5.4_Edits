!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  package monty_applyh
!C
!C  routines to compute matrix elements
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
MODULE applyh
	USE nodeinfo ! Added by SML
	IMPLICIT NONE
CONTAINS

!.....................................................................!
! Subroutine: TBMEmaster
! SML QUESTIONS: do rhopij and rhonij already exist at this point in 
! the code? Sherpa has a seperate subroutine to broadcast rhos that I 
! suspect will be needed in LAMP but am not sure yet where to place 
! or call it. 02/21/2022
!.....................................................................!
	SUBROUTINE TBMEmaster(nsps,rhopij,rhonij,ovlpp,ovlpn,vme)
		USE hamlib,ONLY:usepnform
		USE psis,ONLY:numprot,numneut
		IMPLICIT NONE

		INTEGER :: nsps
		COMPLEX(KIND = 8) :: rhopij(nsps,nsps),rhonij(nsps,nsps)
		COMPLEX(KIND = 8) :: ovlpp,ovlpn
		COMPLEX(KIND = 8) :: vme,vtbmePP,vtbmeNN,vtbmePN

		vtbmePP = (0.d0,0.d0)
		vtbmeNN = (0.d0,0.d0)
		vtbmePN = (0.d0,0.d0)

		IF (numprot > 1) THEN 
			CALL TBMExx(1,nsps,rhopij,vtbmePP)			! Located in this file
		END IF ! numprot 

		IF (numneut > 1 ) THEN 
			IF (usepnform) THEN
				CALL TBMExx(2,nsps,rhonij,vtbmeNN)		! Located in this file
			ELSE ! usepnform 
				CALL TBMExx(1,nsps,rhonij,vtbmeNN)		! Located in this file 
			END IF ! usepnform 
		END IF ! numneut

		IF (numprot > 0) THEN 
			CALL applySPE(1,nsps,rhopij,vtbmePP)		! Located in this file
		END IF ! numprot
			
		IF (numneut > 0) THEN 
			CALL applySPE(2,nsps,rhonij,vtbmeNN)		! Located in this file
		END IF ! numneut 

		IF (numprot > 0 .and. numneut > 0) THEN 
			CALL TBMEpn(nsps,rhopij,rhonij,vtbmePN)	! Located in this file
		END IF  

		! PRINT*, ' Node = ', myMPIrank, ' end of TBMEmaster'

		vme = vtbmePP+vtbmeNN + vtbmePN
		vme = vme*ovlpp*ovlpn

		RETURN 
	END SUBROUTINE TBMEmaster

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	SUBROUTINE TBMExx(it,nsps,rhoij,vtbmeXX)
!C
!C  compute PP/NN matrix elements
!C
		USE nodeinfo
		USE hamlib
		IMPLICIT NONE 

		INTEGER :: it    ! species 1 = protons 2 = neutrons
		INTEGER :: nsps
		COMPLEX(KIND = 8):: rhoij(nsps,nsps)
		COMPLEX(KIND = 8) :: vtbmeXX
		COMPLEX(KIND = 8) :: zsum
		REAL(KIND = 4) :: vtmp
		INTEGER :: itbme
		INTEGER :: a,b,c,d
		INTEGER :: ierr ! For MPI communications
!	  print*,rhoij
!C--------------- loop over matrix elements

!$OMP PARALLEL shared(hmatorbXX,hmatXX), private(a,b,c,d,vtmp,zsum)
!$OMP  do schedule(static) reduction(+:vtbmeXX)
      do itbme = 1,nmatXX
          a = hmatorbXX(it,itbme,1)  
          b = hmatorbXX(it,itbme,2)
          c = hmatorbXX(it,itbme,3)
          d = hmatorbXX(it,itbme,4)

          vtmp = hmatXX(it,itbme)
!C  ------------- FIND ALL PERMUTATIONS 

          zsum = rhoij(a,c)*rhoij(b,d) - rhoij(a,d)*rhoij(b,c) &
               +rhoij(c,a)*rhoij(d,b) - rhoij(d,a)*rhoij(c,b)
	 
!	   if(it==1)print*,' testing ',vtmp,zsum
          vtbmeXX = vtbmeXX + vtmp*zsum
      enddo  ! loop over itbme
!$OMP END  DO 
!$OMP END PARALLEL

      return
      end subroutine TBMExx

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine TBMEpn(nsps,rhopij,rhonij,vtbmePN)
!C
!C  compute PP/NN matrix elements
!C

      use hamlib
      implicit none

      integer nsps
      complex(kind = 8):: rhopij(nsps,nsps),rhonij(nsps,nsps)
      complex(kind = 8) :: vtbmePN
      complex(kind = 8) :: zsum
      real(kind = 4) :: vtmp
      integer itbme
      integer a,b,c,d

!C--------------- loop over matrix elements
!$OMP PARALLEL shared(hmatorbPN,hmatPN), private(a,b,c,d,vtmp,zsum)
!$OMP  do schedule(static) reduction(+:vtbmePN)
      do itbme = 1,nmatPN
          a = hmatorbPN(1,itbme)  
          b = hmatorbPN(2,itbme)
          c = hmatorbPN(3,itbme)
          d = hmatorbPN(4,itbme)

          vtmp = hmatPN(itbme)
!  ------------- FIND ALL PERMUTATIONS 

          zsum = rhopij(a,c)*rhonij(b,d) + rhopij(c,a)*rhonij(d,b) ! including time-reversal
          vtbmePN = vtbmePN + vtmp*zsum
      enddo  ! loop over itbme
!$OMP END  DO 
!$OMP END PARALLEL	  
      return
      end subroutine TBMEpn
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine applySPE(it,nsps,rhoij,vmeX)
!
!  adds single-particle energies
!
      use hamlib
      implicit none

      integer it    ! species 1 = protons 2 = neutrons
      integer nsps
      complex(kind = 8):: rhoij(nsps,nsps)
      complex(kind = 8) :: vmeX
      complex(kind = 8) :: zsum
      real(kind = 4) :: vtmp
      integer itbme
      integer a,b

      do a = 1,nsps
         do b = 1,nsps
           vmeX = vmeX + speunX(b,a,it)*rhoij(b,a)  ! ASSUME SPE symmetric
         enddo  ! loop over b
      enddo  ! loop over a

      return
      end subroutine applySPE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

end module applyh
