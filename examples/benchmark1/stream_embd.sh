!#/bin/bash

root=../../exec/gnu
nframes=5000 
natoms=2
ndim=1
T1=1
T2=1
M1=1
M2=1
Debug=1
parallel=0
NP=2
frmt='xyz'

if [ parallel == 1 ]; then
     mpirun -np $NP ${root}/embdrun.exe ${nframes} ${natoms} ${ndim} ${T1} ${T2} ${M1} ${M2} ${Debug} ${frmt}
else
     ${root}/embdrun.exe ${nframes} ${natoms} ${ndim} ${T1} ${T2} ${M1} ${M2} ${Debug} ${frmt}
fi
