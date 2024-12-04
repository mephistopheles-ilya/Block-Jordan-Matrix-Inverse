#! /bin/bash
for ((n=1; n<=12;n++))
do 
    for ((k=1;k<=4;k+=1))
    do 
        echo "n=$n k=$k --------------" ; 
        ./a.out $n 12 1e-15 $k
        ./a.out $n 12 1e-14 $k
        ./a.out $n 12 1e-13 $k 
    done
done
