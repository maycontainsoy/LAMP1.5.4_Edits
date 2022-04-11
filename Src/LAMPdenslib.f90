!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  package MONTY_DENSLIB.f  
!
module densitylib

	contains

  subroutine makerhoij(it,np,sdf,sdi,ovlp,rhoij,errflag)
!
!  INPUT:
!    it : species = 1 proton, 2 neutron
!    np : = num of particles for this
!    sdf, sdi : final (left) and initial (right) SDs
!
!  OUTPUT:
!    ovlp: = det  sdf^+ sdi
!    rhoij : = density matrix
!           rho = sdi *( sdf^+ sdi )^-1 *sdf^+
!

  	use spstate
!      use system_parameters

  	implicit none
    integer it
    integer np
    complex(kind = 8) :: sdf(nsps,np),sdi(nsps,np)
    complex(kind = 8) :: rhoij(nsps,nsps)
    complex(kind = 8) :: ovlp
	  logical :: errflag

!---------------- INTERMEDIATES ----------------
!complex(kind = 8), allocatable :: ovrtmp(:,:),X(:,:)

		complex(kind = 8) :: ovrtmp(nsps,nsps),X(nsps,nsps)
		complex(kind = 8) :: zzero,tmp
		integer :: i,j,k,l

!------------ FOR DECOMPOSITION
		integer :: ipiv(nsps)
!      integer,allocatable :: ipiv(:)
		integer :: info

!----------- COMPUTE  sdf^+ * sdi = ovrtmp
		errflag = .false.
		zzero = (0.d0,0.d0)
		if(np == 0)then
			ovlp = (1.d0,0.d0)
			rhoij=(0.d0,0.d0)
			return
		endif
	  
!      if(allocated(ovrtmp))deallocate(ovrtmp) !this might be inefficient; fix with pointers
!      allocate(ovrtmp(np,np))

		do i = 1,np
			do j = 1,np
				tmp = zzero
				do k = 1,nsps
					tmp = tmp + dconjg(sdf(k,i))*sdi(k,j)
				enddo ! loop over k
				ovrtmp(i,j) = tmp
!			print*,i,j,real(tmp)
			enddo  ! loop over j
    enddo  ! loop over i 

!------------- INVERT ovrtmp and find determinant = overlap
!              best to use LAPACK routines
!============ First LU decomposition
!      if(allocated(ipiv))deallocate(ipiv)
!      allocate(ipiv(np))

      call zgetrf(np,np,ovrtmp,nsps,ipiv,info)

      if(info /= 0)then
        print*,' problem with LU decomposition ',info
		errflag = .true.
		return
!        stop
      endif

!---------- COMPUTE OVERLAP --------
      ovlp = (1.d0,0.d0)
      do i = 1,np
        ovlp = ovlp*ovrtmp(i,i)
        if(ipiv(i) /= i)ovlp = -ovlp
      enddo

!----------- INVERT. The most efficient way is to solve 
!    OVRTEMP*X = SDF^+ so that X = (OVRTMP)^-1 * SDF^+

      do i = 1,np
         do j = 1,nsps
            x(i,j) = dconjg(sdf(j,i))
         enddo
      enddo
	  
      call ZGETRS( 'N', np, nsps, ovrtmp, nsps,IPIV,X,   nsps, INFO )


!-------------- CONSTRUCT density matrix

      do i = 1,nsps
         do j = 1,nsps
            tmp = zzero
            do k = 1,np
                  tmp = tmp  + sdi(j,k)*X(k,i)

            enddo    ! loop over k
            rhoij(i,j) = tmp
         enddo  ! loop over j
		 101   format(6f8.4)
      enddo  ! loop over i
	       

      return
      end subroutine makerhoij
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc

end module densitylib