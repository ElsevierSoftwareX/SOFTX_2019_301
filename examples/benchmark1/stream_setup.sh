!#/bin/bash

root=../../exec/gnu
natoms=2
ndim=1
Benchmark=1
Debug=1
iStart=1
iSkip=1
iStop=5000 
Offset=0
Ax=0.5 
Cxy=1.0 
Sigma_x=0.5 
Sigma_y=1.0
frmt='xyz'

if [ $Benchmark == 1 ]; then 
     ${root}/sysrun.exe ${iStart} ${iSkip} ${iStop} ${Offset} ${natoms} ${ndim} ${Debug} ${Benchmark} ${Ax} ${Cxy} ${Sigma_x} ${Sigma_y} ${frmt}
else
     ${root}/sysrun.exe ${iStart} ${iSkip} ${iStop} ${Offset} ${Natoms} ${Ndim} ${Debug} ${Benchmark} ${frmt}

fi  

