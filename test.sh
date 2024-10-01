#!/bin/bash

cp matrices/matrix_10.txt matrix.txt

for j in {1..10}
do
    for i in {1..10}
    do
      ./a.out 10 $j $i 0 matrix.txt > mymatrices/mymatrix_$i.txt
      diff -w mymatrices/mymatrix_$i.txt matrices/matrix_$i.txt
    done
done
