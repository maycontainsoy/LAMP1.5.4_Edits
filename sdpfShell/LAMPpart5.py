import os 

os.system('mpiexec -n 1 lamp-mpi.x < Si44_LAMPinTDA > Si44_TDAout')
os.system('mpiexec -n 1 lamp-mpi.x < Si40_LAMPinTDA > Si40_TDAout')
