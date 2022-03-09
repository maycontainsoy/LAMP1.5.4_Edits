module tracy
  implicit none
  
  real(kind=8), allocatable :: normtrace(:), pnormtrace(:),hamtrace(:),phamtrace(:)
  real(kind=8),allocatable :: mnormtrace(:),pmnormtrace(:)
	  real(kind=8),allocatable :: mknorm(:),pmknorm(:)
end module tracy  