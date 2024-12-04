#!/bin/bash
# 4 runs
echo "a.txt --------------"
./a.out 4 4 1e-13 0 matrices/a.txt;
./a.out 4 4 1e-14 0 matrices/a.txt;
./a.out 4 4 1e-15 0 matrices/a.txt;
echo "a20.txt --------------"
./a.out 4 4 1e-13 0 matrices/a20.txt
./a.out 4 4 1e-14 0 matrices/a20.txt
./a.out 4 4 1e-15 0 matrices/a20.txt
echo "a40.txt --------------"
./a.out 4 4 1e-13 0 matrices/a40.txt
./a.out 4 4 1e-14 0 matrices/a40.txt
./a.out 4 4 1e-15 0 matrices/a40.txt
echo "b.txt --------------"
./a.out 4 4 1e-13 0 matrices/b.txt
./a.out 4 4 1e-14 0 matrices/b.txt
./a.out 4 4 1e-15 0 matrices/b.txt
echo "c.txt --------------"
./a.out 6 6 1e-13 0 matrices/c.txt
./a.out 6 6 1e-14 0 matrices/c.txt
./a.out 6 6 1e-15 0 matrices/c.txt
echo "d.txt --------------"
./a.out 6 6 1e-13 0 matrices/d.txt
./a.out 6 6 1e-14 0 matrices/d.txt
./a.out 6 6 1e-15 0 matrices/d.txt
echo "e.txt --------------"
./a.out 6 6 1e-13 0 matrices/e.txt
./a.out 6 6 1e-14 0 matrices/e.txt
./a.out 6 6 1e-15 0 matrices/e.txt
echo "f.txt --------------"
./a.out 4 6 1e-14 0 matrices/f.txt
echo "g.txt --------------"
./a.out 4 6 1e-14 0 matrices/g.txt
