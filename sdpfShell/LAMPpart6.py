import os 

os.system('mpiexec -n 1 lamp-mpi.x < P46_LAMPinTDA > P46_TDAout')
os.system('mpiexec -n 1 lamp-mpi.x < P45_LAMPinTDA > P45_TDAout')
