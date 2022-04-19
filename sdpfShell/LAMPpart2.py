import os 

os.system('mpiexec -n 1 lamp-mpi.x < Al43_LAMPinTDA > Al43_TDAout')
os.system('mpiexec -n 1 lamp-mpi.x < Al42_LAMPinTDA > Al42_TDAout')
