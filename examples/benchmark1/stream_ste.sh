!#/bin/bash

root=../../exec/gnu
nframes=5000 
natoms=2
ndim=1
qTEMethod=1
qTENorm=0
qTEShuffle=1
Nshuffles=10
Rcut=0.02
statP=0.95
Debug=1
parallel=0
NP=2
frmt='xyz'

if [ parallel == 1 ]; then
     mpirun -np $NP ${root}/terun.exe ${nframes} ${natoms} ${ndim} ${qTEMethod} ${qTENorm} ${qTEShuffle} ${Nshuffles} ${Rcut} ${statP} ${Debug} ${frmt}
else
     ${root}/terun.exe ${nframes} ${natoms} ${ndim} ${qTEMethod} ${qTENorm} ${qTEShuffle} ${Nshuffles} ${Rcut} ${statP} ${Debug} ${frmt}
fi
