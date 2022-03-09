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

!      use Projection_Lib,only:factarray
      implicit none

      real(kind = 4) wigner_d
!---------INPUT--------------------------

      real xjj,xmp,xm		! coefficients of d^j_m'm
      real theta		! angle of rotation 

!---------INTERMEDIARIES-----------------
      
      integer i1,i2,i3,i4	! various combination of j-m etc
      real alpha,beta		! alpha = m'-m,  beta = m'+m
      real phase		! from symmetries of d-function 
      				! must have alpha, beta > 0 
      real xmm,xmmp 		! from reordering of m',m 
      real uu			! argument = cos(theta)
      integer n			! order of jacobi polynomial = j-m'
!      real xjacobi		! function = value jacobi polynomial
      real xnorm		! overall factor
      real xnorm2
      integer size		! dimension for factorial array 
      parameter(size=80)
      real fac_ar(0:size)
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


      subroutine lnfact(imax,size,fac_ar)
!
!  returns array fac_ar of ln factorials
!      
      integer imax		! max factorial needed
      integer size
      real fac_ar(0:size)	! fac_ar(i) = ln(i!)
      integer i
      
      if(imax.gt.size)stop ' error in lnfact '
      fac_ar(0) = 0.0
      fac_ar(1) = 0.0
      do i = 2,imax
        fac_ar(i)=fac_ar(i-1)+log(float(i))
      enddo
      return
      end subroutine lnfact
      

!===========================================================


      subroutine jac_param(xjj,xm,xmp,xmm,xmmp, n,alpha,beta,phase)

!------------- returns parameters and, if necessary, a phase for 
!              input into jacobi polynomials

      implicit none 

!------------INPUT------------------------------------------

      real xjj,xm,xmp		! from wigner d^j_m'm
      
!-----------OUTPUT------------------------------------------

      real xmm,xmmp		! reordering/symmetry of m',m
      real alpha,beta		! alpha, beta >= 0
      				! alpha = m'-m,  beta = m'+m
      real phase		! from symmetries/reordering
      integer n			! order of polynomial = j - m'				

!-----------INTERMEDIARIES

      integer iph

!----------------------------------------------------------      
      n=0
      alpha=0.0
      beta=0.0
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
      real(kind =4) :: xjacobi
!******************   use recurrence relation to compute jacobi polys
      integer,parameter :: maxval=60
      real xjac(0:maxval)
      integer n		! order of polynomial
      real x		! argument of polynomial
      real alpha,beta

!----------- INTERMEDIATES -----------------------      

      integer i
      real xn
      real one_nab,two_nab
      real a1,a2,a3,a4

!-------------------------------------------------
      if(n > maxval)then
          print*,' must increase maxval in rotlib.f from ',maxval
          print*,' to ',n
          stop
      end if
      xjac(0)=1.0
      xjac(1)=0.50*(2.00*(alpha+1.00)+(alpha+beta+2.00)* (x-1.00))
      do i=2,n
         xn=float(i-1)
         one_nab=xn+alpha+beta
         two_nab=2.0*xn+alpha+beta
         a1=2.0*(xn+1.0)*(one_nab+1.0)*two_nab
         a2=(two_nab+1.0)*(alpha**2-beta**2)
         a3=two_nab*(two_nab+1.0)*(two_nab+2.0)
         a4=2.0*(xn+alpha)*(xn+beta)*(two_nab+2.0)
         xjac(i)=((a2+a3*x)*xjac(i-1)-a4*xjac(i-2))/a1
      end do 
      xjacobi=xjac(n)
      return 
      end function xjacobi
!===========================================================
end module wignersfriend



