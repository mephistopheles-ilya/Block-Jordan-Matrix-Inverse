#!/bin/bash

for(( m=1; m <=6; m++ ))
do
    for(( p=1; p<=6; p++))
    do
        mpirun -np $p -x ASAN_OPTIONS=detect_leaks=0 ./a.out 6 $m 6 0 matrices/c.txt
    done
done
