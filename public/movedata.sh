#!/bin/bash

mkdir ../../P2data/$1
mv tso1 ../../P2data/$1/tso1  
mv ev1 ../../P2data/$1/ev1 
mv net_diag1 ../../P2data/$1/net_diag1
cp zna.dat ../../P2data/$1/
cp Stable_Nuclides.txt ../../P2data/$1/
cp read_ts_nc.py ../../P2data/$1/
cp read_ts_mf.py ../../P2data/$1/
