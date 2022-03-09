!
!   PROGRAM lilLAMP
!   reads in matrices to compute various expectation values
!
!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  

	  use little
	use lileigenpackage   ! in file LAMPeigen.f90
    use dens1b
	implicit none
    real (kind=4) :: tolerance, elapsedTime, Jtolerance
	real(kind=4)  :: newJmax
	integer(kind=8) :: nlevelmax  ! # of levels possible
	integer(kind=8) :: nftotal    ! actual number of levels???
	integer :: ij
		
    real (kind=8), allocatable :: jall(:), pallPair(:,:),obsall(:,:)!!
    real (kind=8),  allocatable  :: jlist(:)!!
    real (kind = 8) :: normSum, hamSum
    real (kind=8), allocatable  :: problist(:,:), hamlist(:,:),obslist(:,:)!!
	character :: ychar,achar
	
	logical :: chooseJ = .true.
	
	integer :: num_threads,mythread
	integer :: omp_get_thread_num,omp_get_num_threads

    real (kind=8) :: clock_start, clock_stop,x
	real (kind =8) :: norm_start,norm_stop,ham_start,ham_stop
	integer(4) :: numsdused  ! allows user to reduce number of SDs used
	integer(4)  :: nprint  ! how many to print out
	
	character :: eval_char,menu_char
	

   

    print*,' '
	print*,' Welcome to lilLAMP! ( Oct 2021 version )'
	print*,' I read in matrices from LAMP and compute expectation values '
	print*,' You must first read in the Hamiltonian .mat file '
	print*,' '
	
	call read_matfromfile(.true.)
	
	
	menu_char = 'n'
	
	do while(menu_char == 'n')
	   print*,'  Choose from the following: '
	   print*,' * (e) Computing energy spectrum '
	   print*,' * (x) Compute expectation value of scalar operator '
	   print*,' * (d) Compute one-body density matrices '
	   read(5,'(a)')menu_char
	   
	   select case (menu_char)
	   
	   case ('e','E')
	   
	   menu_char='e'
	   print*,' Evaluating energy spectrum '
	   
	   case('x','X')
	   menu_char='x'
	   print*,' Evaluating expectation values '
	   
	   case('d','D')
	   menu_char='d'
	   print*,' Evaluating one-body transitition density matrices '
	   
	   case default 
	   
	   print*,' I do not recognize the choice of ',menu_char
	   menu_char='n'
	   
	   
       end select
	
     end do
	
savesolutions=.false.	
	
select case (menu_char)

! .................... COMPUTE ENERGY-SPECTRUM.............
case('e')	

nlevelmax = 0
do ij = 1,numOfJ
	nlevelmax = nlevelmax + nint(2* xJlist(ij)+1)
end do
nlevelmax = nlevelmax*numsd
allocate(jall(nlevelmax),pallPair(2,nlevelmax),obsall(2,nlevelmax))
allocate(problist(2,numOfJ),hamlist(2,numOfJ))

print*,' How many energies to print to screen? (all will be written to file)'
read*,nprint

!........ ADDED IN 1.4.2... ABILITY TO DO CHANGES IN TOLERANCE ........
ychar = 'y'
call inputTolerance(tolerance)	
jtarget = -1.
numsdused = numsd
stdeigen=.true.
compute_expect = .false.
do while (ychar=='y' .or. ychar=='Y')
    call EigenSolverPackage(tolerance, nlevelmax,numsdused,nftotal,jall,pallPair,obsall,normSum,hamSum,problist,hamlist)

    print *, '' 
    print *,' Sum of norms = ', dble(normSum)
    print *,' Sum of trace(H) = ', dble(hamSum)

    call J_WriteResults(tolerance,nlevelmax,nftotal,jall,pallPair,obsall,problist,hamlist,ParityTest,nprint)!!

    print*,' '
    print*,' Do you want to run with other parameters? (y/n)'
    read(5,'(a)')ychar
    if(ychar=='y'.or. ychar=='Y')then
        call inputTolerance(tolerance)
		if(numsd > 1)then
			print*,' Enter number of SDs to used (max ',numsd,')'
			read*,numsdused
			if(numsdused > numsd)numsdused=numsd
		end if

			stdeigen=.true.   ! FOR NOW RUN WITH STANDARD EIGENSOLVER
        print*,' Choose a target J to print out (y/n)?'
	    read(5,'(a)')achar
	     if(achar=='y'.or.achar=='Y')then
			 print*,' Enter target J '
			 read*,jtarget
		 else
			 jtarget = -1.
		 end if
    end if
	

 end do ! while

!......................COMPUTE EXPECTATION VALUES...........

case ('x')
	
	print*,' '
	print*,' You must now read in the observable .mat file '
	print*,' '
	call read_matfromfile(.false.)
	
	
    nlevelmax = 0
	do ij = 1,numOfJ
		nlevelmax = nlevelmax + nint(2* xJlist(ij)+1)
	end do
	nlevelmax = nlevelmax*numsd
	allocate(jall(nlevelmax),pallPair(2,nlevelmax),obsall(2,nlevelmax))
	allocate(problist(2,numOfJ),hamlist(2,numOfJ))
	
	print*,' How many energies to print to screen? (all will be written to file)'
	read*,nprint

	!........ ADDED IN 1.4.2... ABILITY TO DO CHANGES IN TOLERANCE ........
	ychar = 'y'
	call inputTolerance(tolerance)	
	jtarget = -1.
	numsdused = numsd
	stdeigen=.true.
	compute_expect = .true.
	do while (ychar=='y' .or. ychar=='Y')
        call EigenSolverPackage(tolerance, nlevelmax,numsdused,nftotal,jall,pallPair,obsall,normSum,hamSum,problist,hamlist)

        print *, '' 
        print *,' Sum of norms = ', dble(normSum)
        print *,' Sum of trace(H) = ', dble(hamSum)

        call J_WriteResults(tolerance,nlevelmax,nftotal,jall,pallPair,obsall,problist,hamlist,ParityTest,nprint)!!

        print*,' '
        print*,' Do you want to run with other parameters? (y/n)'
        read(5,'(a)')ychar
        if(ychar=='y'.or. ychar=='Y')then
	        call inputTolerance(tolerance)
			if(numsd > 1)then
				print*,' Enter number of SDs to used (max ',numsd,')'
				read*,numsdused
				if(numsdused > numsd)numsdused=numsd
			end if

				stdeigen=.true.   ! FOR NOW RUN WITH STANDARD EIGENSOLVER
            print*,' Choose a target J to print out (y/n)?'
		    read(5,'(a)')achar
		     if(achar=='y'.or.achar=='Y')then
				 print*,' Enter target J '
				 read*,jtarget
			 else
				 jtarget = -1.
			 end if
        end if
		

     end do ! while

!.............COMPUTE DENSITY MATRICES....................

case ('d')

savesolutions = .true.
compute_expect = .false.
stdeigen =.true.
!......... READ IN UNTRANSFORMED DENSITY MATRICES..........

call readin_densitymat

!......... solve generalized eigenvalue problem and store ......
nlevelmax = 0
do ij = 1,numOfJ
	nlevelmax = nlevelmax + nint(2* xJlist(ij)+1)
end do
nlevelmax = nlevelmax*numsd
allocate(jall(nlevelmax),pallPair(2,nlevelmax),obsall(2,nlevelmax))
allocate(problist(2,numOfJ),hamlist(2,numOfJ))

print*,' How many energies to print to screen? '
read*,nprint

!........ ADDED IN 1.4.2... ABILITY TO DO CHANGES IN TOLERANCE ........
ychar = 'y'
call inputTolerance(tolerance)	
jtarget = -1.
numsdused = numsd

call EigenSolverPackage(tolerance, nlevelmax,numsdused,nftotal,jall,pallPair,obsall,normSum,hamSum,problist,hamlist)

print *, '' 
print *,' Sum of norms = ', dble(normSum)
if(dble(normSum)< 0.99)then
	print*,' (Note: if sum of norms < 1, it may be due to restriction on J)'
	print*,' (Check your original LAMP run if you have concerns )'
end if

print *,' Sum of trace(H) = ', dble(hamSum)


call J_WriteResults(tolerance,nlevelmax,nftotal,jall,pallPair,obsall,problist,hamlist,ParityTest,nprint)

call sortresults
call mapresults

call open_output(33,menu_char)




call write_density_results(33)


end select

end

!
!  ROUTINE TO WRITE GENERALIZED MATRICES TO FILE
!  added 1.4.4  July 2019 by CWJ @ SDSU
!
!  NOTE this version assumes some files have already been created
!
subroutine read_matfromfile(hamflag)
	use little
	use sporbit,only:numorb
	implicit none
	logical :: hamflag
	
	type (jaggedArray), pointer :: Nx(:,:,:),PNx(:,:,:),Hx(:,:,:),PHx(:,:,:)
	character*80 ::filename
	integer :: ilast
	
	integer :: intJ,matdim
	real :: floatJ
	integer :: ii,jj,sdX, sdY
	
	integer :: dum1,dum2,dum3
	logical :: ptest
	real :: xJ
	integer :: ij
	complex(kind=8) :: ctest
	
	logical :: binary_file
	
	binary_file=.false.
	print*,' '
	
	if(hamflag)then
	    print*,' Enter name of .mat file for Hamiltonian (do not enter extension)'
	else
	    print*,' Enter name of .mat file for observable (do not enter extension)'
	end if
	read(5,'(a)')filename
	ilast = index(filename,' ')-1
	open(unit=3,file=filename(1:ilast)//".mat",status='old',err=101)
	
	read(3,*,err=101)dum1,dum2
	
	goto 111
101 continue  

    print*,' File may be unformatted (binary); trying that '
	binary_file=.true.
	
	close(3)
	open(unit=3,file=filename(1:ilast)//".mat",status='old',form='unformatted')
		
	read(3)dum1,dum2
	
	
111 continue	
	
	
	if(hamflag)then
		numorb = dum1
		nsps   = dum2
	else
	
	    if(dum1 /= numorb)then
		   print*,' Error, mismatch in num of orbits ',numorb,dum1
		   stop
	    end if
	    if(dum2 /= nsps)then
		    print*,' Error, mismatch in num of s.p. states ',nsps,dum2
		    stop
	   end if
   end if
!	write(2,*)numorb,nsps

    if (binary_file)then
	    read(3)dum1,dum2,dum3
		
	else
	    read(3,*)dum1,dum2,dum3
		
	end if
	
	
	if(hamflag)then
		numprot=dum1
		numneut =dum2
		numsd  = dum3
	else
	   if(dum1 /= numprot)then
	   	   print*,' Error, mismatch in num prot',numprot,dum1
		   stop
	    end if
	    if(dum2 /= numneut)then
		   print*,' Error, mismatch in num neutrons',numneut,dum2
		   stop
	    end if
	    if(dum3 /= numsd)then
		   print*,' Error, mismatch in num SDs',numsd,dum3
		   stop
	    end if
	end if
	if(binary_file)then
	    read(3)dum1,numOfJObs,ptest
	else
	    read(3,*)dum1,numOfJObs,ptest
	end if


	if(hamflag)then
		ParityTest = .not.ptest
		if(ParityTest)then
			npair = 2
		else
			npair = 1
		end if
		numOfJ=numOfJObs
	else
	
        if(ptest .eqv. ParityTest)then
		   print*,' Error, mismatch in parity ',ParityTest,ptest
		   stop
	    end if
	    if(numOfJObs > numOfJ)then
		   print*,' Need to increase maximum J ',numOfJ,numOfJObs
		   stop
	    end if
	end if
	
!............ SET UP xJlist
    if(hamflag)then
       allocate(xJlist(numOfJ),mapJ2indx(0:numOfJ*2),startJ2list(0:numOfJ*2))

       if(mod(numprot+numneut,2) == 0)then
            isOdd = .false.
       else
           isOdd = .true.
       endif

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
       end do
   end if
!................... NOW ALLOCATE.....................

    if(hamflag)then
       allocate(N_J_MK(numsd,numsd,numOfJ))	
       allocate(H_J_MK(numsd,numsd,numOfJ))	
 	   if(Paritytest) allocate(PN_J_MK(numsd,numsd,numOfJ),PH_J_MK(numsd,numsd,numOfJ))	
	   
	   do intJ = 1,numOfJ
		
          floatJ = xJlist(intJ)
		  matdim = nint(2*floatJ)+1
		  N_J_MK(1,1,intJ)%dimMK=matdim
		  H_J_MK(1,1,intJ)%dimMK=matdim
		  if(ParityTest)then
			  PN_J_MK(1,1,intJ)%dimMK=matdim
			  PH_J_MK(1,1,intJ)%dimMK=matdim			  
			  
		  end if
		  
		  do sdX =1,numsd
			  do sdy = 1,numsd
				  allocate(N_J_MK(sdx,sdy,intJ)%MK(matdim,matdim))
				  allocate(H_J_MK(sdx,sdy,intJ)%MK(matdim,matdim))
				  if(ParityTest)allocate(PN_J_MK(sdx,sdy,intJ)%MK(matdim,matdim))
				  if(ParityTest)allocate(PH_J_MK(sdx,sdy,intJ)%MK(matdim,matdim))	
			   end do  !sdy
		  end do  !sdx
		
	   end do  !intJ
	   
	   Nx => N_J_MK
	   Hx => H_J_MK
	   if(ParityTest)then
		   PNx => PN_J_MK
		   PHx => PH_J_MK
	   end if
	   
		
	else
		
		
       allocate(Nobs_J_MK(numsd,numsd,numOfJ))	
       allocate(Hobs_J_MK(numsd,numsd,numOfJ))	
	   if(Paritytest) allocate(PNobs_J_MK(numsd,numsd,numOfJ),PHobs_J_MK(numsd,numsd,numOfJ))	

	   do intJ = 1,numOfJ
		
          floatJ = xJlist(intJ)
		  matdim=nint(2*floatJ)+1
		  Nobs_J_MK(1,1,intJ)%dimMK=matdim
		  Hobs_J_MK(1,1,intJ)%dimMK=matdim
		  if(ParityTest)then
			  PNobs_J_MK(1,1,intJ)%dimMK=matdim
			  PHobs_J_MK(1,1,intJ)%dimMK=matdim
			  
		  end if
		  do sdX =1,numsd
			  do sdy = 1,numsd
				  allocate(Nobs_J_MK(sdx,sdy,intJ)%MK(matdim,matdim))
				  allocate(Hobs_J_MK(sdx,sdy,intJ)%MK(matdim,matdim))
				  if(ParityTest)allocate(PNobs_J_MK(sdx,sdy,intJ)%MK(matdim,matdim))
				  if(ParityTest)allocate(PHobs_J_MK(sdx,sdy,intJ)%MK(matdim,matdim))	
			   end do  !sdy
		  end do  !sdx
		
	   end do  !intJ
	   
	   Nx => Nobs_J_MK
	   Hx => Hobs_J_MK
	   if(ParityTest)then
		   PNx => PNobs_J_MK
		   PHx => PHobs_J_MK
	   end if
	   
   end if
    do intJ = 1, numOfJ
        floatJ = xJlist(intJ)
		matdim = N_J_MK(1,1,intJ)%dimMK		
		do sdX = 1,numsd
			do sdY = 1,numsd

				do ii = 1,matdim
					if(binary_file)then
						read(3)(Nx(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
						
					else
						read(3,*)(Nx(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
						
					end if
				end do !ii
				
			end do !sdY
		end do ! sdX
		do sdX = 1,numsd
			do sdY = 1,numsd
				do ii = 1,matdim
					if(binary_file)then
						read(3)(Hx(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
						
					else
						read(3,*)(Hx(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
						
					end if
				end do !ii
			end do !sdY
		end do ! sdX
		if(.not.ParityTest)cycle
		do sdX = 1,numsd
			do sdY = 1,numsd
				do ii = 1,matdim
					if(binary_file)then
						read(3)(PNx(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
						
					else
						read(3,*)(PNx(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
						
					end if
				end do !ii
			end do !sdY
		end do ! sdX
		do sdX = 1,numsd
			do sdY = 1,numsd
				do ii = 1,matdim
					if(binary_file)then
						read(3)(PHx(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
						
					else
						read(3,*)(PHx(sdX,sdY,intJ)%MK(ii,jj),jj=1,matdim)
						
					end if
				end do !ii
			end do !sdY
		end do ! sdX		
		
	end do  ! int J

	close(3)
	return
end subroutine read_matfromfile


!=============================================



subroutine open_output(outfile,menu_char)
	use little
	implicit none
	integer :: outfile
	character :: menu_char
	
	character(40) :: filename
	integer :: ilast
	
	
	if(menu_char=='d')then
	   print*,' Enter name of .dres outputfile for densities (do not include extension) '
    else
 	   print*,' Enter name of .res outputfile for spectra (do not include extension) '
		
	end if   
	read(5,'(a)')filename
	ilast = index(filename,' ')-1
	
	if(menu_char=='d')then
   	   open(unit=outfile,file=filename(1:ilast)//".dres",status='unknown')
	   print*,' Results written to ', filename(1:ilast),'.dres '
	   write(outfile,*)' LAMP/lilLAMP version 1.5.1 May 2021 '
	   write(outfile,*)' '
	   write(outfile,'(2i6)')numprot,numneut
	   if( mod(numprot+numneut,2)==0)then
		   write(outfile,*)'0 +'
		   
	   else
		   write(outfile,*)'1 +'
		   
		   
	   endif
   else
   	   open(unit=outfile,file=filename(1:ilast)//".res",status='unknown')
	   print*,' Results written to ', filename(1:ilast),'.res '
	   write(outfile,'("# Results from LAMP version 1.4.9 ")')
	   
   end if
	
   return
	
end subroutine open_output

