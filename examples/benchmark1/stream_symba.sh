!#/bin/bash

root=../../exec/gnu
nframes=5000 
natoms=2
ndim=1
qSymbolic=1
Debug=1
nmc=1000
parallel=0
NP=2
frmt='xyz'

if [ parallel == 1 ]; then
     mpirun -np $NP ${root}/symbrun.exe ${nframes} ${natoms} ${ndim} ${qSymbolic} ${Debug} ${nmc} ${frmt}
else
     ${root}/symbrun.exe ${nframes} ${natoms} ${ndim} ${qSymbolic} ${Debug} ${nmc} ${frmt}
fi
