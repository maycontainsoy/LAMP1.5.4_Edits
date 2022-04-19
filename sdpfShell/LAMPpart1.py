import os 

os.system('mpiexec -n 1 lamp-mpi.x < Na37_LAMPinTDA > Na37_TDAout')
os.system('mpiexec -n 1 lamp-mpi.x < Na35_LAMPinTDA > Na35_TDAout')
os.system('mpiexec -n 1 lamp-mpi.x < Mg40_LAMPinTDA > Mg40_TDAout')
os.system('mpiexec -n 1 lamp-mpi.x < Mg38_LAMPinTDA > Mg48_TDAout')
os.system('mpiexec -n 1 lamp-mpi.x < Mg37_LAMPinTDA > Mg37_TDAout')
os.system('mpiexec -n 1 lamp-mpi.x < Mg35_LAMPinTDA > Mg35_TDAout')
os.system('mpiexec -n 1 lamp-mpi.x < Mg33_LAMPinTDA > Mg33_TDAout')
os.system('mpiexec -n 1 lamp-mpi.x < Mg32_LAMPinTDA > Mg32_TDAout')
