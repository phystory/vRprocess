#!/bin/bash

for f in `ls -d ../trajectories/trajectories*/`; do
  sed 's/Ye/Ye\n0.900000E+01    Start Time\n0.280000E+08    Stop Time\n1.000000E-12    Init Del t/' $f/trajectory_00001.dat >| trajectory_00001.dat 
  echo "  0.3141E+08 -0.3323E+03 -0.2678E+03 -0.8414E+02  0.1278E+00 -7.2000E+05  0.3867E-01" >> trajectory_00001.dat
  ./source/xnetp 
  OLDIFS=$IFS
  IFS=/
  set $f
  echo $3
  ./movedata.sh $3
  IFS=$OLDIFS
done
