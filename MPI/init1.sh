#!/bin/bash

flag=0
for(( n=1; n <= 10; n++))
do
    for(( m=1; m <= n; m++))
    do
        for(( r=0; r < 10; r++))
        do
            for(( s=1; s <= 4; s++))
            do
                for(( p=1; p <= 6; p++))
                do
                    mpirun -np $p -x ASAN_OPTIONS=detect_leaks=0 ./a.out $n $m $r $s > mpi_res.txt
                    ../Linear/a.out $n $m $r $s > lin_res.txt
                    diff mpi_res.txt lin_res.txt
                    if [ $? -ne 0 ]; then
                    echo "Test failed"
                    echo "p = $p n = $n m = $m r = $r s = $s"
                    flag=1
                    break
                    fi
                done
                if [ $flag -ne 0 ]; then
                    break;
                fi
            done
            if [ $flag -ne 0 ]; then
                break;
            fi
        done
        if [ $flag -ne 0 ]; then
            break;
        fi
    done
    if [ $flag -ne 0 ]; then
        break;
    fi
done
