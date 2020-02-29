!#/bin/bash

root=../../exec/gnu
nframes=10000 
natoms=2
ndim=1
qMIMethod=1
qMIShuffle=0
Nshuffles=10
Rcut=0.02
statP=0.95
Debug=1
parallel=0
NP=2
frmt='csv'

if [ parallel == 1 ]; then
     mpirun -np $NP ${root}/mirun.exe ${nframes} ${natoms} ${ndim} ${qMIMethod} ${qMIShuffle} ${Nshuffles} ${Rcut} ${statP} ${Debug} ${frmt}
else
     ${root}/mirun.exe ${nframes} ${natoms} ${ndim} ${qMIMethod} ${qMIShuffle} ${Nshuffles} ${Rcut} ${statP} ${Debug} ${frmt}
fi
