! Added by SML
! Copied from Sherpa 
MODULE nodeinfo 
  IMPLICIT NONE 
	INCLUDE 'mpif.h'

	INTEGER :: nMPIranks, myMPIrank, icomm 
	INTEGER :: num_global_threads
	INTEGER :: root 

END MODULE nodeinfo 