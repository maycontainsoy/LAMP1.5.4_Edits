import os 

os.system('mpiexec -n 1 lamp-mpi.x < S46_LAMPinTDA > S46_TDAout')
os.system('mpiexec -n 1 lamp-mpi.x < S44_LAMPinTDA > S44_TDAout')
os.system('mpiexec -n 1 lamp-mpi.x < Cl45_LAMPinTDA > Cl45_TDAout')
