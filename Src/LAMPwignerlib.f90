!
!  library of functions for computing wigner_d matrices
!

module wignersfriend
	implicit none
	integer :: maxn4fact
	real(kind=8), allocatable :: factarray(:)
contains

      function wigner_d(xjj,xmp,xm,theta)

!------------computes Wigner (little) d-function d^J_m'm(beta)
!------------cf Edmonds 
!     in double precision

      implicit none

      real(kind = 8) wigner_d
!---------INPUT--------------------------

      real(8) :: xjj,xmp,xm		! coefficients of d^j_m'm
      real(8) theta		! angle of rotation 

!---------INTERMEDIARIES-----------------
      
      integer i1,i2,i3,i4	! various combination of j-m etc
      real(8) alpha,beta		! alpha = m'-m,  beta = m'+m
      real phase		! from symmetries of d-function 
      				! must have alpha, beta > 0 
      real(8) :: xmm,xmmp 		! from reordering of m',m 
      real(8) :: uu			! argument = cos(theta)
      integer n			! order of jacobi polynomial = j-m'
!      real xjacobi		! function = value jacobi polynomial
      real(8) :: xnorm		! overall factor
!      real xnorm2
!      integer size		! dimension for factorial array 
!      parameter(size=80)
!      real fac_ar(0:size)
      integer imax

!-------------------------------------------------------------

      if(abs(xmp).gt.xjj)stop 'xmp > xjj in wigner_d '
      if(abs(xm).gt.xjj)stop 'xm > xjj in wigner_d '
      
!------------------------    put M values into order for
!                            Jacobi polynomials

      call jac_param(xjj,xmp,xm,xmm,xmmp,n,alpha,beta,phase)

      i1=nint(xjj+xmm)
      i2=nint(xjj-xmm)
      i3=nint(xjj+xmmp)
      i4=nint(xjj-xmmp)
     
      imax = max(i1,i2)
      imax = max(imax,i3)
      imax = max(imax,i4)
!      call lnfact(imax,size,fac_ar)

!      xnorm=exp((fac_ar(i1)+fac_ar(i2)-fac_ar(i3)-fac_ar(i4))/2.)
      xnorm=exp( (factarray(i1)+factarray(i2)  -factarray(i3)-factarray(i4) )/2.0)
!	  if(abs(xnorm-xnorm2)> 0.1)then
!		  print*,xnorm,xnorm2,i1,i2,i3,i4
!	  end if
		 
		 
      
      uu=cos(theta)
      i1=nint(beta)
      i2=nint(alpha)
	  
      wigner_d=phase*xnorm*cos(theta/2.)**i1*sin(theta/2.)**i2*xjacobi(n,alpha,beta,uu)

      return
      end function wigner_d
!===========================================================

      subroutine jac_param(xjj,xm,xmp,xmm,xmmp, n,alpha,beta,phase)

!------------- returns parameters and, if necessary, a phase for 
!              input into jacobi polynomials

      implicit none 

!------------INPUT------------------------------------------

      real(8) :: xjj,xm,xmp		! from wigner d^j_m'm
      
!-----------OUTPUT------------------------------------------

      real(8) :: xmm,xmmp		! reordering/symmetry of m',m
      real(8) :: alpha,beta		! alpha, beta >= 0
      				! alpha = m'-m,  beta = m'+m
      real phase		! from symmetries/reordering
      integer n			! order of polynomial = j - m'				

!-----------INTERMEDIARIES

      integer iph

!----------------------------------------------------------      
      n=0
      alpha=0.d0
      beta=0.d0
      phase=1.0

      xmm = xm
      xmmp = xmp
      if(xjj.eq.0.00)return
      if(xm.ge.0.00)then
         if(xmp.ge.0.00)then
             if(xm.ge.xmp)then
                xmm=xm
                xmmp=xmp
                phase=1.0
             else
                xmm=xmp
                xmmp=xm
                iph=nint(xm-xmp)
                phase=(-1.0)**iph
             end if
         elseif(xmp.lt.0.0)then
             if(xm.ge.abs(xmp))then
                xmm=xm
                xmmp=xmp
                phase=1.00
             else
                xmm=-xmp
                xmmp=-xm
                phase=1.0
             end if
         end if
      elseif(xm.lt.0.0)then
         if(xmp.ge.0.0)then
             if(abs(xm).ge.xmp)then
                xmm=-xm
                xmmp=-xmp
                iph=nint(xm-xmp)
                phase=(-1.0)**iph
             else
                xmm=xmp
                xmmp=xm
                iph=nint(xm-xmp)
                phase=(-1.00)**iph
             end if
         elseif(xmp.lt.0.0)then
             if(abs(xm).ge.abs(xmp))then
                xmm=-xm
                xmmp=-xmp
                iph=nint(xm-xmp)
                phase=(-1.0)**iph
             else
                xmm=-xmp
                xmmp=-xm
                phase=1.0
             end if
         end if
      end if

      alpha=xmm-xmmp
      beta=xmm+xmmp
      if(alpha.lt.0.0)then
          write(6,*)xm,xmp,xmm,xmmp
          stop 'error in wigner_d alpha < 0'
      end if
      if(beta.lt.0.0)then
          write(6,*)xm,xmp,xmm,xmmp
         stop 'error in wigner_d alpha < 0'
      end if
      n=nint(xjj-xmm)
      return
      end subroutine jac_param
!===========================================================

      function xjacobi(n,alpha,beta,x) 
!
!  computes jacobi polynomial P^(alpha, beta)_n(x)
!

      implicit none
      real(kind =8) :: xjacobi
!******************   use recurrence relation to compute jacobi polys
      integer,parameter :: maxval=60
      real(8) :: xjac(0:maxval)
      integer n		! order of polynomial
      real(8) :: x		! argument of polynomial
      real(8) :: alpha,beta

!----------- INTERMEDIATES -----------------------      

      integer i
      real(8) :: xn
      real(8) :: one_nab,two_nab
      real(8) :: a1,a2,a3,a4

!-------------------------------------------------
      if(n > maxval)then
          print*,' must increase maxval in rotlib.f from ',maxval
          print*,' to ',n
          stop
      end if
      xjac(0)=1.d0
      xjac(1)=0.5d0*(2.d00*(alpha+1.0d0)+(alpha+beta+2.0d0)* (x-1.0d0))
      do i=2,n
         xn=float(i-1)
         one_nab=xn+alpha+beta
         two_nab=2.0*xn+alpha+beta
         a1=2.d0*(xn+1.0)*(one_nab+1.d0)*two_nab
         a2=(two_nab+1.0)*(alpha**2-beta**2)
         a3=two_nab*(two_nab+1.0)*(two_nab+2.0)
         a4=2.0*(xn+alpha)*(xn+beta)*(two_nab+2.0)
         xjac(i)=((a2+a3*x)*xjac(i-1)-a4*xjac(i-2))/a1
      end do 
      xjacobi=xjac(n)
      return 
      end function xjacobi
!===========================================================
! obtains log of factorial
! to initialize array of ln factorials, set n < 0, n = -max n desired
! 
SUBROUTINE Factorial(n,fact)

IMPLICIT NONE

REAL (KIND = 8), INTENT(IN) :: n        !Integer to be factorialized
REAL (KIND = 8), INTENT(OUT) :: fact        !Variable to be exported -- factorial = n!
REAL (KIND = 8) :: temp, d_i            ! i -- index for DO loop, temp -- dummy variable for accumulation of n!
INTEGER :: i,i_n

i_n = INT(n)

if(i_n > maxn4fact)then
	print*,' Need to increase maxn4fact ',maxn4fact,i_n
	stop
end if
	
if(i_n >= 0)then
	fact = factarray(i_n)
	return
end if
!................. OTHERWISE SET UP
maxn4fact = - i_n
allocate(factarray(0:-i_n))

temp = DLOG(DBLE(1.0))

factarray(0)=0.d0
factarray(1)=0.d0

DO i = 2,-i_n
    d_i = DBLE(i)
    temp = temp+DLOG(d_i)
	factarray(i)=temp
END DO

!fact = temp
return

END SUBROUTINE Factorial
!===========================================================
subroutine make_rotmat(alpha,beta,gamma,RotMat)
	USE spstate
	!use wignersfriend

	IMPLICIT NONE

	REAL (KIND = 8), INTENT(IN) :: alpha,beta,gamma
	COMPLEX (KIND = 8), DIMENSION(nsps,nsps), INTENT(OUT) :: RotMat
	integer :: a,b
	integer :: ja,jb,ma,mb,na,nb,la,lb
	real*8 :: dtmp
	complex(kind=8) :: eye = (0.d0, 1.d0)
	
	rotmat(:,:)=dcmplx(0.d0,0.d0)
	do a = 1,nsps
		ja = spsqn(a)%j
		ma = spsqn(a)%m
		na = spsqn(a)%nr
		la = spsqn(a)%l
		do b = 1,nsps
			jb = spsqn(b)%j
			mb = spsqn(b)%m
			nb = spsqn(b)%nr
			lb = spsqn(b)%l
			if(ja/=jb)cycle
			if(na/=nb)cycle
			if(la/=lb)cycle
			dtmp = wigner_d(ja*0.5d0,ma*0.5d0,mb*0.5d0,beta)
!			print*,' wigner ',a,b,beta,dtmp
			rotmat(a,b)=exp(eye*0.5d0*ma*gamma)*dtmp*exp(eye*0.50*mb*alpha) ! Edmonds 4.1.12
			
		end do
		
	end do
	
	
	return
end subroutine make_rotmat


!===========================================================
end module wignersfriend



