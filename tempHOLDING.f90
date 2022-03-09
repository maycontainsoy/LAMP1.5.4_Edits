module psis
    !
    !  wave functions--slater determinants
    !
    ! subsumes old module system_parameters 
                    USE nodeinfo
          implicit none
          
          integer numprot,numneut,numsd ! from old module system_parameters
    
          logical special_override_p, special_override_n
          integer :: p_override, n_override
          type sdmaster
              complex(kind=8),pointer :: psd(:,:),nsd(:,:)
          end type sdmaster
    
          type(sdmaster), pointer :: psir(:),psil(:)  ! of it = time
              
          complex (kind=8), allocatable :: psd(:,:),nsd(:,:)
          complex (kind=8), allocatable :: psdf(:,:,:),nsdf(:,:,:)
    
      contains
    
          subroutine readinsd(nsps,n,ifile,psi,errflag)
    !
    !  subroutine to read in SD from file written by SHERPA
    !  NB: probably will have other routines as well
    !
    !  INPUT:
    !  nsps   = # of single-particle states
    !  n      = # of particles
    !  ifile  = logical number of the file
    !  
    !  OUTPUT
    !   psi(i,j) = slater determinant (real)
    !   errflag  = true if there was a problem in reading
    !
    
    
          implicit none
          integer nsps,n
          real psi(nsps,n) ! slater determinant
    
          integer i,j
          integer ifile
          logical errflag
    
          errflag=.false.
    
          do i = 1,n
            read(ifile,err=103,end=103)(psi(j,i),j=1,nsps)
          enddo
          return
      103 continue
          errflag=.true.
          return
    
          end subroutine readinsd  
    !=====================================================================
        subroutine set_nuclide
            implicit none
                INTEGER :: ierr ! for MPI communication
    
            special_override_p =.false.
            special_override_n =.false.
            
            IF (myMPIrank == root) THEN
                PRINT*,' Enter Z, N '
                READ*,numprot,numneut			
                    END IF 
    
                    ! Update other ranks
                    CALL MPI_BARRIER(icomm,ierr)
                    CALL MPI_BCAST(numprot,1,MPI_INT,root,icomm,ierr)
                    CALL MPI_BCAST(numneut,1,MPI_INT,root,icomm,ierr)
                    CALL MPI_BARRIER(icomm,ierr)
            
            return
        end subroutine set_nuclide      	
    
    !=====================================================================
    
          subroutine GetSherpaSD(pSD,nSD)
    !
    !  master routine for retrieving SD written out by SHERPA
    !  must also convert them from real*4 to complex*8
    !
    !      use system_parameters
          use spstate
          implicit none
    
          integer ifile
          complex(kind = 8) ::pSD(nsps,numprot),nSD(nsps,numneut)
          real(kind=4),allocatable :: sdtmp(:,:)
    
          ifile = 73
    
          call OpenSherpaSD(ifile,nsps,numprot,numneut)
          
    
    
          if(allocated(sdtmp))deallocate(sdtmp)
          if(numprot > 0 .or. special_override_p)then
              
              if(special_override_p)then
                  allocate(sdtmp(nsps,p_override))
                  call ReadSherpa(ifile,nsps,p_override,sdtmp)  
                  ! DO NOT CONVERT
                  
              else
               allocate(sdtmp(nsps,numprot))
               call ReadSherpa(ifile,nsps,numprot,sdtmp)
               call ConvertSD(nsps,numprot,sdtmp,pSD)
              end if
    
    
            deallocate(sdtmp)
          endif
          if(numneut > 0 .or. special_override_n)then
              
              if(special_override_n)then
                  allocate(sdtmp(nsps,n_override))
                  call ReadSherpa(ifile,nsps,n_override,sdtmp)  
                  ! DO NOT CONVERT
                  
              else
                allocate(sdtmp(nsps,numneut))
                 call ReadSherpa(ifile,nsps,numneut,sdtmp)
                 call ConvertSD(nsps,numneut,sdtmp,nSD)
                end if
            deallocate(sdtmp)
          endif
          close(ifile)
    
          return
          end subroutine GetSherpaSD
    
    !=====================================================================      
    
          subroutine OpenSherpaSD(ifile,nsps,Z,N)
              
    !		  use system_parameters
    
          use sporbit
          implicit none
    
    !.......... 
          integer N,Z
          integer nsps
          integer ifile
    !..............FILE HANDLING..........................
    
    
          character ychar*1
          character filename*25  		! 
          integer ilast
          integer tempfile			! location of temporary file
          data tempfile/99/
          logical errflag
          character title*60
    
          logical openflag,failflag
    
    !-------------- DUMMIES 
          integer i,ii,j
          integer norb
          integer zz,nn
    
    !-------------- OPEN FILE ---------------
    
    1     continue
    
          openflag = .false.
          do while(.not.openflag)
                    IF (myMPIrank == root) THEN 
                write(6,*)' Enter input filename (.sd)'
                read(5,'(a)')filename
                ilast=index(filename,' ')
                if(ilast.ne.0)then
                  ilast=ilast-1
                else
                  ilast=15
                endif
    
                open(unit=ifile,file=filename(1:ilast)//'.sd',status='old', err=33, form='unformatted')
                openflag = .true.
    33      continue
    
                IF (.not.openflag) THEN
                  write(6,*)' That file does not exist; do you wish to try another file (y/n)?'
                   read(5,'(a)')ychar
                   if(ychar == 'N' .or. ychar == 'n')stop ! might need to add MPI abort here instead
                        END IF 
                    END IF ! myMPIrank == root
          END DO   ! while on openflag
    
    !..............READ IN HEADER INFO...............
    !              CHECK THAT MATCHES SYSTEMS
    !
    
          read(ifile)norb
          failflag = .false.
          if(norb /= numorb)then
        write(6,*)' # of orbits mismatch ',norb,numorb
            failflag = .true.
          endif
          do i = 1,numorb
        read(ifile)ii,norb,j
        if(norb /= orbqn(i)%nr .or. j /= orbqn(i)%j)then
          write(6,*)' mismatch n,j:',norb,orbqn(i)%nr,j,orbqn(i)%j
          failflag = .true.
            endif
          enddo
    
          if(failflag)then
    
            write(6,*)' The single particle space does not match ', &
         ' with that previously chosen.'
            write(6,*)' Exit (x) or choose another file (c)?'
            read(5,'(a)')ychar
            if(ychar.eq.'x' .or. ychar.eq.'X')stop
            close(ifile)
            goto 1
          endif
    
    !...............CHECK IF N,Z match........
    
          read(ifile)zz,nn
    !      if(n.ne.nn .or. z.ne.zz)then
              
    !.......... ADDED in 1.5.1..... special override
           if( n/= nn)then
               if(special_override_n)then
                   print*,' not reading in neutron SD '
                   n_override = nn
               elseif(n==0)then
                   print*,' Do you mean to ignore the neutron SD?'
                   read(5,'(a)')ychar
                   if(ychar=='y'.or. ychar=='Y')then
                       special_override_n=.true.
                       n_override = nn
                   end if
                   
               else
                   write(6,*)' mismatch N: Old: ',n,', new: ',nn
                   write(6,*)' Exit (x) or choose another file (c)?'
                   read(5,'(a)')ychar
                   if(ychar.eq.'x' .or. ychar.eq.'X')stop
                   close(ifile)
                   goto 1
                   
               end if
                   
             
           end if
           if( z/= zz )then
               if(special_override_p)then
                   print*,' not reading in proton SD '
                   p_override = zz
                   
               elseif(z==0)then
                   print*,' Do you mean to ignore the proton SD?'
                   read(5,'(a)')ychar
                   if(ychar=='y'.or. ychar=='Y')then
                       special_override_p=.true.
                       p_override = zz
                   end if
                   if(special_override_p .and. special_override_n)then
                       print*,'cannot have an empty nucleus!'
                       stop
                   end if
                   
               else
                   write(6,*)' mismatch Z: Old: ',z,', new: ',zz
                   write(6,*)' Exit (x) or choose another file (c)?'
                   read(5,'(a)')ychar
                   if(ychar.eq.'x' .or. ychar.eq.'X')stop
                   close(ifile)
                   goto 1
                   
               end if
                   
             
           end if		  
    
    
          read(ifile)title
    
            write(6,'(60a)')'Title card: "',title,'"'
          return
          end subroutine OpenSherpaSD
    
    !=====================================================================
          subroutine ReadSherpa(ifile,nsps,np,sd)
    
    !....... NOTE THIS ASSUMES real*4 
          implicit none
    
          integer nsps,np
          real(kind=4),intent(OUT) :: sd(nsps,np) ! slater determinant
    
          integer i,j
          integer ifile
          logical errflag
    
          errflag=.false.
    
          do i = 1,np
            read(ifile,err=103,end=103)(sd(j,i),j=1,nsps)
    
          enddo
          return
      103 continue
          errflag=.true.
          return
    
          return
          end subroutine ReadSherpa
    
    !=====================================================================
          
    
          subroutine ConvertSD(nsps,np,sd,zsd)
    
          implicit none
          integer nsps,np
    
          real(kind=4) :: sd(nsps,np)
          complex(kind = 8) :: zsd(nsps,np)
    
          integer i,j
    
          do i = 1,nsps
             do j = 1,np
                zsd(i,j) = dcmplx( sd(i,j),0.0)
             enddo     ! loop over j
          enddo  ! loop over i
    
          return
          end subroutine ConvertSD
    
    !=====================================================================
          subroutine read_sd_txt(psdtrial,nsdtrial)
    !=====================================================================
    ! Reads a (time-reversal) Slater determinant
    ! in text format
    !====================================================================
    
    !      use system_parameters
          use spstate
    
          implicit none
    
    
          integer  :: i,j,k,n,indx,ilast  
          complex(kind=8) :: psdtrial(nsps,numprot),nsdtrial(nsps,numneut)
          real vec(nsps)
    
          character :: ychar
          character*25 :: filename
    
      311 continue
          write(6,*)' Enter input filename (.tsd)'
          read(5,'(a)')filename
          ilast=index(filename,' ')
          if(ilast.ne.0)then
              ilast=ilast-1
          else
              ilast=15
          endif
    
          open(unit=10,file=filename(1:ilast)//'.tsd',status='old',err=33,form='formatted')
          goto 44
      33  continue
          write(6,*)' That file does not exist; do you wish to try another file (y/n)?'
          read(5,'(a)')ychar
          if(ychar.eq.'n' .or.ychar.eq.'N')then
             return
          else
            goto 311
          endif
      44  continue
    
          psdtrial(:,:)=dcmplx(0.d0,0.d0)
          if(numneut/=0)nsdtrial(:,:)=dcmplx(0.d0,0.d0)
    
          do n=1,numprot
              read(10,*)(vec(i),i=1,nsps)
              do i=1,nsps
                 psdtrial(i,n)=dcmplx(dble(vec(i)),0.d0)
              enddo
           enddo
           do n=1,numneut
              read(10,*)(vec(i),i=1,nsps)
              do i=1,nsps
                 nsdtrial(i,n)=dcmplx(dble(vec(i)),0.d0)
              enddo
    
          enddo
      122 continue
          close(10)
    
          return
          end subroutine read_sd_txt
    !
    !==============================================================================
    !
    ! allocate the slater determinants psd / nsd & psdf / nsdf
    !
    ! OUTPUT:
    !   psd, nsd = proton/neutron slater determinants
    !   psdf, nsdf = ^^ with a third dimension for multiple slater deterimants (i.e. numsd > 1)
    !   ParityTest = checks for proper parity [?]
    !
    ! SUBROUTINES CALLED:
    !   PairLog = checks for parity and flags PairLog appropriately
    !
    subroutine allocateSlaterDet 
    
        use spstate
    
    
        implicit none
    
    
    !............INTERNAL....................
        integer :: ii,jj,iii
        integer :: ifile
        logical :: errflag
        integer :: nSDs,isd,jsd,ksd
    !	real(kind=4),allocatable :: psdtmp(nsps,numprot),nsdtmp(nsps,numneut)
    
        real(kind=4),allocatable :: psdtmp(:,:),nsdtmp(:,:)
        complex(kind = 8) ::pSDx(nsps,numprot),nSDx(nsps,numneut)
        real :: eee
    
        INTEGER :: ierr ! for MPI communication
        
        ifile  = 73
    
    !----------------- ALLOCATE SLATER DETERMINANTS -------------------
      allocate(psd(nsps,numprot),nsd(nsps,numneut))
        allocate(psdtmp(nsps,numprot))
        allocate(nsdtmp(nsps,numneut))
        
      IF (myMPIrank == root) THEN 
            print*,' Enter number of Slater determinants '
            read*,numsd
        END IF ! myMPIrank 
    
        !CALL MPI_BARRIER(icomm,ierr)
        CALL MPI_BCAST(numsd,1,MPI_INT,root,icomm,ierr)
        CALL MPI_BARRIER(icomm,ierr)
    
        PRINT*, 'node # ', myMPIrank, ' here 1'
    
        ! PRINT*,' node # ', myMPIrank, ' num SDs ', numsd
    
    !	print*,' testing ',numsd,nsps,numprot,numneut, ' TESTING'
        allocate (psdf(numsd,nsps,numprot),nsdf(numsd,nsps,numneut))
        jj=0
        do ii = 1,numsd
            if(jj==numsd)return
    
               call OpenSherpaSD(ifile,nsps,numprot,numneut)
               if(special_override_p)then
                deallocate(psdtmp)
                   allocate(psdtmp(nsps,p_override))
    
               end if	
               if(special_override_n)then
                deallocate(nsdtmp)
                   allocate(nsdtmp(nsps,n_override))
               end if
    
    !... special override is for separating proton and neutron SDs....
           if(special_override_p)call count_SDs_in_file(ifile,nsps,p_override,numneut,nSDs)
           if(special_override_n)call count_SDs_in_file(ifile,nsps,numprot,n_override,nSDs)
    !........... 'NORMAL' CASE.....................	   
           if(.not.special_override_p .and. .not.special_override_n)then
               call count_SDs_in_file(ifile,nsps,numprot,numneut,nSDs)
           end if
           
           
           rewind(ifile)
           call read_SD_header(ifile,.false.,errflag)
           
           if(nSDs==1)then
               if(numprot > 0)then			   
                   call ReadSherpa(ifile,nsps,numprot,psdtmp)
                   call ConvertSD(nsps,numprot,psdtmp,pSDx)
               end if
               if(special_override_p)then			   
                   call ReadSherpa(ifile,nsps,p_override,psdtmp)
               end if
                   
               if(numneut > 0)then
                   call ReadSherpa(ifile,nsps,numneut,nsdtmp)
                   call ConvertSD(nsps,numneut,nsdtmp,nSDx)
               end if
               if(special_override_n)then
                   call ReadSherpa(ifile,nsps,n_override,nsdtmp)
               end if
               jj = jj+1
               psdf(jj,:,:) = psdx
               nsdf(jj,:,:) = nsdx
               close(ifile)
               
               
           else
               print*,' There are ',nSDs,' Slater determinants in this file '
               print*,' Enter a Slater determinant index (0 or -1 to stop )'
               print*,' To read from an ordered list, enter -999 '
               
               read*,isd
               if(isd==-999)then
                   open(unit=33,file="orderedSDlist.dat",status='old')
                   do iii=1,nSDs
                       read(33,*,end=9909)ksd,eee
                       print*,ksd,eee
                       rewind(ifile)
                       call read_SD_header(ifile,.false.,errflag)
                       do jsd = 1,ksd
                           if(numprot > 0)then
                               call ReadSherpa(ifile,nsps,numprot,psdtmp)
                               call ConvertSD(nsps,numprot,psdtmp,pSDx)
                           end if
                           if(special_override_p)then
                               call ReadSherpa(ifile,nsps,p_override,psdtmp)
                           end if
                           if(numneut > 0)then
                               call ReadSherpa(ifile,nsps,numneut,nsdtmp)
                               call ConvertSD(nsps,numneut,nsdtmp,nSDx)
                           end if
                           if(special_override_n)then
                               call ReadSherpa(ifile,nsps,n_override,nsdtmp)
                           end if
                       
                       end do
                       jj = jj+1
                       if(numprot>0)psdf(jj,:,:) = psdx
                       if(numneut>0)nsdf(jj,:,:) = nsdx
                   
                       if(jj==numsd)then
                           close(ifile)		   
                           return
                       end if
                   
                       
                   end do
    9909           continue
                                  
               end if
               
               do while(isd > 0 .and. isd <= nSDs)
                   rewind(ifile)
                   call read_SD_header(ifile,.false.,errflag)
                   do jsd = 1,isd
                       if(numprot > 0)then
                           call ReadSherpa(ifile,nsps,numprot,psdtmp)
                           call ConvertSD(nsps,numprot,psdtmp,pSDx)
                           
                       end if
                       if(special_override_p)then
                           call ReadSherpa(ifile,nsps,p_override,psdtmp)
                       end if
                       if(numneut > 0)then
                           call ReadSherpa(ifile,nsps,numneut,nsdtmp)
                           call ConvertSD(nsps,numneut,nsdtmp,nSDx)
                           
                       end if
                       if(special_override_n)then
                           call ReadSherpa(ifile,nsps,n_override,nsdtmp)
                           
                       end if
                       
                   end do
                   jj = jj+1
                   if(numprot>0)psdf(jj,:,:) = psdx
                   if(numneut>0)nsdf(jj,:,:) = nsdx
                   
                   if(jj==numsd)then
                       close(ifile)		   
                       return
                   end if
                   
                   
                   print*,' Enter a Slater determinant index (0 or -1 to stop )'
                   read*,isd
               end do
               
               
           end if
    
        end do
    
    return
    
    end subroutine allocateSlaterDet	  
    !======================================================
    subroutine read_SD_header(ifile,verbose,failflag)
        
        use sporbit
    !	use system_parameters
        
        implicit none
        integer,intent(IN) :: ifile
        logical,intent(IN) :: verbose
        logical,intent(OUT) :: failflag
        
        integer :: norb,i,ii,j
        integer :: nn,zz
        character title*60
        
        
        failflag = .false.
        
        read(ifile)norb
     
        if(norb /= numorb)then
             write(6,*)' # of orbits mismatch ',norb,numorb
             failflag = .true.
        endif
        do i = 1,numorb
           read(ifile)ii,norb,j
           if(norb /= orbqn(i)%nr .or. j /= orbqn(i)%j)then
              write(6,*)' mismatch n,j:',norb,orbqn(i)%nr,j,orbqn(i)%j
              failflag = .true.
            endif
        enddo
    
        if(failflag)then
            close(ifile)
            return
        end if
    !
    !...............CHECK IF N,Z match........
    
          read(ifile)zz,nn
          if(((numneut.ne.nn) .and. .not.special_override_n) .or. & 
             ((numprot.ne.zz) .and. .not. special_override_p))then
              write(6,*)' mismatch Z,N. Old: ',numprot,numneut, ', new: ',zz,nn
              failflag = .true.
              close(ifile)
              return
          endif
    
          read(ifile)title
    
          if(verbose)write(6,'(60a)')'Title card: "',title,'"'
            
        return
        
        
    end subroutine read_SD_header
    
    !===================================================================
    
    subroutine count_SDs_in_file(ifile,nsps,z,n,nSDs)
        implicit none
        integer,intent(IN) :: ifile
        integer,intent(IN) :: nsps,z,n
        integer,intent(OUT) :: nSDs
        integer :: isd
        logical :: errflag
        
        real(kind=4) :: psdtmp(nsps,z),nsdtmp(nsps,n)
        
        nSDs = 0
        
        do isd = 1,100000
            errflag=.false.
               if(z > 0)call ReadSherpa2(ifile,nsps,z,psdtmp,errflag)
            if(errflag)return
               if(n > 0)call ReadSherpa2(ifile,nsps,n,nsdtmp,errflag)
            if(errflag)return
            nSDs = nSDs+1
        end do
        print*,' Wait should not have gotten here ',nSDs
        stop
        
    
    end subroutine count_SDs_in_file
    !===================================================================
          subroutine ReadSherpa2(ifile,nsps,np,sd,errflag)
    
    !....... NOTE THIS ASSUMES real*4 
    !use spstate, only:spsqn
          implicit none
    
          integer nsps,np
          real(kind=4),intent(OUT) :: sd(nsps,np) ! slater determinant
    
          integer i,j
          integer ifile
          logical errflag
    
          errflag=.false.
    
          do i = 1,np
            read(ifile,err=103,end=103)(sd(j,i),j=1,nsps)
    !		print*,' sd ',sd(:,1)
    
          enddo
    !	  print*,spsqn(:)%m
          return
      103 continue
          errflag=.true.
          return
    
          return
          end subroutine ReadSherpa2
        
    
      end module psis
    