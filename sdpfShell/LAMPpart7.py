import os 

os.system('mpiexec -n 1 lamp-mpi.x < P43_LAMPinTDA > P43_TDAout')
os.system('mpiexec -n 1 lamp-mpi.x < P42_LAMPinTDA > P42_TDAout')
