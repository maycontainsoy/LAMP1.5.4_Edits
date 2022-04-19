import os 

os.system('mpiexec -n 1 lamp-mpi.x < Al40_LAMPinTDA > Al40_TDAout')
os.system('mpiexec -n 1 lamp-mpi.x < Al39_LAMPinTDA > Al39_TDAout')
