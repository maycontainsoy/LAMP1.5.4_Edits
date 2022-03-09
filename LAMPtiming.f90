MODULE LAMPtiming

  CONTAINS 

!.....................................................................!
! Subroutine: clocker 
! Pulled from Bigstick file btimelib.f90 
!.....................................................................!
SUBROUTINE clocker(timer, starter)
	IMPLICIT NONE 

	CHARACTER(LEN=3) :: timer, starter 
	CHARACTER(LEN=12) :: real_clock(3)

	INTEGER(KIND=4) :: time_values(8)
	
	REAL(KIND=8) :: timenow 
	REAL(KIND=8) :: wall_time 
	REAL(KIND=8) :: smin, shour, sday 

	smin = 60.0d0
	shour = smin * 60.0d0
	sday = shour * 24.0d0

	

END SUBROUTINE clocker

END MODULE LAMPtiming 