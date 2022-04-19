import os 

os.system('mpiexec -n 1 lamp-mpi.x < Al38_LAMPinTDA > Al38_TDAout')
os.system('mpiexec -n 1 lamp-mpi.x < Al37_LAMPinTDA > Al37_TDAout')
